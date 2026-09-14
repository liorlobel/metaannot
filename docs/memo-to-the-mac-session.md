# Memo for the Mac session (console/GUI work)

Written 2026-09-14, against `main` at `1d5489c` (PR #59 merged, CI green).
Nothing here is urgent for correctness on macOS or Linux. Two of the three are
things you are better placed to fold in than I am, because they touch
`console/console.py` and you are the one editing it.

---

## 1. `console/console.py` cannot be imported on Windows at all

`console/console.py:39` is a bare `import fcntl` at module scope, with **no
platform guard anywhere in the file** — so on Windows the module dies at
import, before any of it runs.

Measured on this machine just now:

```
tests/test_console.py + tests/test_console_contract.py
  1 failed, 16 passed, 281 errors
```

The 281 errors are all the same import. `test_the_console_imports_only_the_standard_library`
is the one real failure, for the same reason. Linux CI is green and cannot see
any of this, because it never runs on Windows.

**It is a small fix.** There is only one genuine `fcntl` use site in 4,114
lines:

| line | thing | note |
|---|---|---|
| 39 | `import fcntl` | module scope, unguarded |
| 3688 | `fcntl.flock(fd, LOCK_EX \| LOCK_NB)` | the only real call |
| 3673 | `os.O_NOFOLLOW`, `os.O_NONBLOCK` | neither attribute exists on Windows |
| 353 | `os.O_NONBLOCK` | same |

`metaannot.py:151-155` already has the pattern to copy, with the reasoning
written out beside it:

```python
try:
    import fcntl
    import select
except ImportError:                                # pragma: no cover - Windows
    fcntl = select = None
```

and `getattr(os, "O_NOFOLLOW", 0)` for the two open flags.

**What it does NOT fix**, and should not pretend to: the advisory lock at 3688
is what stops two consoles serving one directory, and Windows has no `flock`.
Guarding the import makes the module importable and the other 281 tests
meaningful; it leaves that one guarantee weaker on Windows, which is worth
saying out loud in the code rather than papering over.

**I did not make this change**, deliberately — you are editing this file and I
did not want to hand you a conflict in it. If you would rather I did, say so.

---

## 2. The console's version-mismatch warning is now checkable, and wasn't before

`console/console.py:2943` does:

```python
if run.get("version") and run["version"] != v.get("engine_version"):
```

and the comment above it is right about why it exists: the engine behind
`--metaannot` need not be the engine that wrote the directory.

**On real data that check silently passes when it should fire.** The machine
that produced the 8-plex cohort holds five `metaannot.py` files that all report
`0.2.0` and are **four different builds**. The cohort's own log opens
`metaannot 0.2.0 starting`. Point the console at any of those four and
`run["version"] == engine_version` compares equal, so no warning appears — in
exactly the case the warning is for. (`uc_run/run.log` is worse: it says
`1.0.0`, a build no file on that machine still holds.)

PR #59 fixes the missing half. Every run now records the sha256 of the file
that ran:

- `_run.source_sha256` in `.metaannot_state.json` — full 64-char digest
- `metaannot_source` in **`describe --json`** and **`doctor --json`** — same digest
- the banner: `metaannot 0.7.2 starting (source 902824f154b3)` — the short
  form is the first 12 hex of the same digest, `SOURCE_DIGEST_SHORT`

So the check above can compare digests instead of version strings and become
sound. Roughly:

```python
ran = run.get("source_sha256")
have = v.get("engine_source")          # from describe --json metaannot_source
if ran and have and ran != have:
    ...different BUILD, whatever the version strings say...
```

Two notes if you take it:

- **Both keys are additions**, so `DESCRIBE_VERSION` and `DOCTOR_VERSION` did
  not move (their own documented rule: a bump is for a key REMOVED or a meaning
  CHANGED). A console that has never heard of `metaannot_source` is unaffected.
  `CONSOLE_VERSION` is yours and I have not touched it.
- **Older directories have no `source_sha256`.** Absent is the normal state for
  anything written before #59, and it is not a mismatch — fall back to the
  version comparison rather than warning on a missing key.

The header at 2934 would also be a reasonable place to show the short digest
beside `metaannot %s`, since "which build wrote this directory" is precisely
the question that header is answering.

---

## 3. What else changed under you in #59

Only one thing that a console reads, and it is the addition above. For
completeness:

- `MEASURED_455K` in `tests/test_scheduler.py` was corrected — it had SignalP's
  208,396 s filed in `interpro`'s cell. `order_s` on `STAGES` was **not**
  affected, so the dispatch order the console displays is unchanged.
- README / TUTORIAL / CLAUDE.md corrected on `unipept` and `taxonomy` (both have
  now run on real data — seventeen of the twenty-one stages have, not fifteen)
  and on the 8-plex run, which finished in 84.2 h.
- **Upgrading recomputes nothing**: all 21 stage signatures compare equal
  against the previous commit, and no stage's `keys` list changed.

---

## 4. One caveat if you diff re-runs

Re-deriving both cohorts under the current DIAMOND identity floors moved MEROPS
by +2 and TADB by −1 on UC — and **neither database is floored**. That is
DIAMOND being non-deterministic across runs, not a code change. If the console
ever diffs two runs of one directory, small deltas on unfloored databases are
noise and should not be presented as drift.
