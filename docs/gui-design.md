# A console for metaannot — design

Status: **Phase 0 and M1 shipped in v0.5.0**; M2 onwards is still design. This
document is the design record it was written to be, not a manual — `console/`
is documented for users in README.md. It is kept because the reasoning is what
a later milestone has to argue against, and it is annotated **Built:** where
what shipped settled a question differently from what is written below. Where
a decision here was made and honoured, it stands as written.

The goal is what FragPipe is to MSFragger: a front end that removes the setup
friction without becoming a second place where the science is decided.

---

## What it is, in one sentence

**A watcher and a config author, running on the host that writes the results
directory it is watching, serving one page to `127.0.0.1`.**

It never runs a stage. It never parents a run. It never holds a path belonging
to a different operating system than the pipeline.

> **Built:** the first half, and not to `127.0.0.1`. What shipped is a watcher
> only — M1 — and it binds a mode-0700 UNIX socket reached over `ssh -L`, never
> a TCP port. The mitigation floated at the end of "Front end" below became the
> design: a loopback listener on a shared lab server is reachable by every other
> account on that box, and the token alternative leaks into `ps`, shell history
> and the URL bar, so file permissions are the whole access control. One clause
> got stronger in the building: the shipped console writes **nothing at all**
> into a results directory, which is a rule the socket has to obey too — a
> `--socket` inside or under a watched tree is refused rather than accepted,
> because binding there creates both the socket and a lock file.

## What it is not

- **Not a scheduler.** metaannot already has one, with a dependency graph, a
  CPU/RAM budget and a GPU lease.
- **Not a second source of truth for config semantics.** It generates its form
  from `DEFAULT_CONFIG` and its stage list from `STAGES`. If it ever hard-codes
  what a key means, it will drift from the engine and start lying. That is the
  failure mode this whole design is arranged to prevent.
- **Not a run owner.** It attaches to nothing, so there is no reattach feature
  to write and nothing to lose when it is closed.

> **Built:** both honoured, the first by a stricter route than "reads
> `DEFAULT_CONFIG`". The shipped console does not import metaannot at all —
> stage order, dependency edges, the `_run` key and the paths it polls come out
> of one `describe --json`, shelled out at startup and cached. It also acquired
> a third refusal that belongs in this list: it writes nothing into a results
> directory, ever, which is what makes it safe to point at a running job.

---

## Target platforms

| platform | pipeline | console |
| --- | --- | --- |
| **Linux x86-64** | all 21 stages | runs here |
| **Windows + WSL2** | all 21 stages, inside WSL | runs inside WSL, viewed in a Windows browser |
| **macOS, Apple Silicon** | ~19 of 21 | runs here |
| Windows, native | **not supported** | read-only fallback view only |

Windows native is not a pipeline host and will not become one. HMMER requires
POSIX and has no native Windows build, and `pfam`, `dbcan`, `ncbifam`,
`jackhmmer` and `kofam` all invoke `hmmsearch` directly. WSL2 is the Windows
answer and it is a first-class one — 21 of 21 stages, identical to Linux.

### What a Mac cannot do, stated plainly

Native `osx-arm64` builds exist for HMMER, DIAMOND, MMseqs2, Foldseek and
hhsuite; KOfamScan and eggNOG-mapper are noarch; Bioconductor ships arm64
binaries for limma and QFeatures. So **metaannot's entire default-on stage set
runs natively on a Mac with no emulation**, and the R report is native too.

Two stages do not, and no amount of front-end design changes that:

- **`interpro`** — excluded by upstream policy, not by accident. InterProScan 6
  is a Nextflow pipeline whose images are all `linux/amd64`, so on Apple
  Silicon it would mean emulating the single most expensive stage.
- **`structure`** (ESMFold, and Foldseek with it) — gated by metaannot's own
  `if not torch.cuda.is_available(): die(...)`. Fixable in principle via
  Torch's MPS backend, but that is a rewrite of the most delicate stage, not a
  flag.

`topology` is available but impractical: TMbed has no MPS path and loads
ProtT5 in fp16 only inside its CUDA branch, so a Mac silently takes the CPU
fallback and runs for a very long time. It does not error — which is worse.

**So a Mac is a real development and small-dataset host, and a real console
host for watching a run on a Linux box. It is not where a 455k-protein run
happens.** The console should say this on its setup screen rather than letting
someone discover it after starting.

---

## Why a watcher, and why that is the whole design

Three facts about the engine, each verified in the source, force this shape.

### 1. The job model already exists on disk

`results/.metaannot_state.json` carries per-stage `{signature, status,
started, finished, seconds}`, written temp-file-then-`os.replace` under a lock.
`mark_running()` is stamped **before** a stage starts, and its docstring says
why: *"Recorded before the stage starts, so that a run killed mid-write leaves
a trace."* Beside it sit an append-only `metaannot.log` and a `.metaannot.lock`
holding `{pid, host, started}`.

Tool stdout goes to `/dev/null`, so there is no pipe worth owning.

A console that polls those files is a **watcher**. "Survives being closed for
three days" then stops being a feature to engineer and becomes a property of
never having attached. Closing the tab, rebooting the desktop, and opening cold
on day three are the same code path.

### 2. Fourteen of twenty-one stages hash a database path *string*

This is the subtle one, and it is worth getting exactly right because the
obvious version of it is wrong.

File **inputs** are hashed by content, not by path. `_stat()` says so directly:
*"Identity of one input: its content, not the text of its path"* — respelling
`results` as an absolute path used to recompute everything, and that was fixed
deliberately. So a GUI that rewrites the *protein FASTA* path is safe.

But `signature()` also hashes `{k: _dig(cfg, k) for k in stage["keys"]}`, by
value — and **14 of the 21 stages list a path-valued `db.*` key**:

```
emapper   db.eggnog_data          interpro   db.interproscan_sh
pfam      db.pfam_hmm             jackhmmer  db.jackhmmer_db
dbcan     db.dbcan_hmm            hhblits    db.hhblits_db
diamond   db.diamond              foldseek   db.foldseek_target, extra_targets
ncbifam   db.ncbifam_hmm          taxonomy   db.ncbi_taxonomy
kofam     db.kofam_profiles,      context    gff
          db.kofam_ko_list        unipept/join  quant_table, manifest
```

So a console that helpfully rewrites `D:\db\Pfam-A.hmm` into
`/mnt/d/db/Pfam-A.hmm` — same file, different string — invalidates those
stages and restarts InterProScan. Thirty-four hours, silently, for a cosmetic
edit.

**The fix is structural, not a rule someone has to remember.** If the file
picker is a server-side endpoint listing directories on the pipeline host, it
is *physically incapable* of returning a foreign path. The database paths are
exactly what the preflight screen edits, so this is not a hypothetical.

### 3. metaannot is importable and self-describing

`tests/conftest.py` already loads it by path and hundreds of tests drive it
in-process. `STAGES` is 21 entries carrying `name / enabled / deps / keys /
gpu`, and `requirements()` already returns the preflight screen's exact data
model. The UI generates itself from the engine's own definitions.

> **Built:** self-describing, yes; importable, deliberately unused. `describe
> --json` landed in Phase 0 and the console reads that instead, which is what
> lets it be one stdlib-only file that can be `scp`'d to a machine where
> metaannot lives at a path it was never told about. The contract is versioned
> (`DESCRIBE_VERSION`) so a console can say "I do not understand this shape"
> rather than guess, and it carries more than this section anticipated: the
> config vocabulary (`path_keys`, `db_path_keys`, `replace_blocks`,
> `freeform_keys`, `retired_keys`), the two keys the scheduler sorts each
> round's ready set by (`cost`, then `order_s` — the stage's seconds on the
> release's reference run, which breaks the eleven-way tie the rank alone
> leaves), and `paths` — the
> files a watcher polls, named rather than reconstructed, so the day one of
> them moves the watcher moves with it.

---

## Front end: a local page, with an optional native window

**v1 serves one page to `127.0.0.1` and is opened in a browser.** A `--window`
flag later renders the identical HTML in a native OS window via `pywebview`
(525 kB, BSD-3), giving a dock entry and an app icon, and degrading to the
browser when its system WebKit packages are absent.

> **Built:** a browser, yes; `127.0.0.1`, no, and no `--window` either. The
> reasoning that changed it is three paragraphs down, in this section, and the
> **Built:** note under it is where the decision ended up. Nothing that shipped
> opens a TCP port.

A native-only app was seriously considered and rejected for one reason: it
excludes the headless Linux server, which is metaannot's own documented
deployment — `examples/server-run-plan/` describes eight datasets over SSH and
tmux, assuming no display anywhere. On such a host every native toolkit needs
X11 forwarding, and `sshd`'s `X11Forwarding` defaults to `no`; `ssh -L` needs
nothing on either end and is the access the runbook already uses. It also
excludes Mac-at-desk-against-Linux-host: macOS has shipped no X server since
10.8, and an X-forwarded window dies with the connection and cannot be
reattached — the exact failure "run it in tmux" exists to prevent.

Most arguments for native did not survive examination. `http://localhost` is a
secure context, so a served page can raise real OS notifications — and
"run finished" belongs in the *engine* anyway, where it also serves the tmux
user with no UI open. Native file pickers are neutral, because the axis that
matters is which **host** runs the picker, not which toolkit draws it.

One argument did survive: a `127.0.0.1` listener on a shared lab server is
reachable by every other account on that box. Mitigation: bind a `0700` UNIX
socket (OpenSSH `-L` accepts socket forms), or port 0 plus a random token and
an `Origin`/`Host` check.

> **Built:** the first mitigation, unconditionally — there is no TCP mode to
> fall back to. A token was rejected outright rather than kept as an option: it
> leaks into `ps`, into shell history and into the URL bar, and it would be a
> second access-control mechanism to keep correct. File permissions are the
> whole of it, which is why the socket's directory has to be private and is
> checked for that, and why the console refuses to create one for you. The
> price is a documented one: a UNIX-socket forward target needs OpenSSH 6.7 or
> newer on the workstation side.
>
> `--window` and `pywebview` did not ship and are not scheduled. The socket
> plus `ssh -L` covers the headless server, the Mac-at-desk case and the local
> one with no second rendering path to keep working.

`tkinter` was considered because it is in the standard library and would fit
the project's instincts. It is present in this lab's WSL Python, but it is
routinely absent from Homebrew and pyenv Pythons on macOS, so "one file, no
package to install" would not survive contact with a Mac.

---

## Engine changes first

Each of these is worth making even if the console is never built. They land
individually, with a test, before any UI exists.

| change | why | size |
| --- | --- | --- |
| Handle `SIGTERM` as `SIGINT` is handled | Lock release is only `atexit`. Python's default SIGINT unwinds and runs it; SIGTERM does not. So `kill`, `systemctl stop` and `wsl --terminate` all leave a stale lock and the next run refuses to start. | ~4 lines |
| A top-level `_run` record in the state file — `{run_id, version, config_path, argv, host, pid, started, finished, final_status}` | The state file has no run-level metadata, so a results directory opened cold cannot say what produced it. With eight datasets on one machine that is the difference between a job list and a pile of directories. | ~15 lines |
| A `_run.last_seen` heartbeat on a ~30 s timer | `_holder_is_alive()` returns True unconditionally for another host and always on Windows, so nothing outside the process can tell "InterProScan, hour 19" from "the box died". | ~10 lines |
| `config.effective.yaml` written beside the results | Records the merged config actually used, defaults included. Reproducibility, not UI. | ~10 lines |
| `describe --json` | Emits `DEFAULT_CONFIG`'s shape, `STAGES` and `requirements()` as JSON, so the console reads a contract rather than importing private functions. | ~40 lines |

Until `describe --json` exists, v1 **imports** metaannot for description only and
**shells out** for every action, with `tests/test_console_contract.py` pinning
the names and shapes the console reads. That test is what stops the coupling
being invisible.

> **Built:** all five landed in Phase 0, and two of them differently from the
> row above.
>
> The `SIGTERM` row asked for "handled as `SIGINT` is handled", and that is
> what shipped first and was then changed. Unwinding means waiting for the
> stage pool's workers, and a tmbed chunk or an InterProScan stage is an hour;
> `systemd` reaches `TimeoutStopSec` long before that and sends `SIGKILL`, so
> the lock release the row exists for is exactly what unwinding loses. `run`
> and `all` now install a handler that removes the lock, writes one
> pre-encoded line to fd 2 and `os._exit(128 + N)`s. Ctrl-C still unwinds and
> still stamps `_run`; a killed run leaves the lock gone and the record still
> reading `running`, which is the state the heartbeat exists to let a watcher
> describe — and the reason the console must never call it death.
>
> The `_run` row's field list grew: `run_id`, `version`, `config_path`, `argv`,
> `host`, `pid`, `started`, `last_seen`, `last_seen_epoch`, `heartbeat_s`,
> `finished`, `final_status`. Two timestamps for the same instant, because two
> hosts sharing one filesystem cannot subtract each other's local clocks.
>
> `tests/test_console_contract.py` exists and does more than pin names: it
> proves by AST scan and by snapshotting a results directory around every route
> that the console writes nothing into one.

---

## Milestones

**Phase 0 — engine changes.** ✅ **shipped, v0.5.0.** No UI. Proves the engine
half stands alone: no orphaned locks after `kill`, results directories that say
what produced them, a config record that is reproducible.

**M1 — watch only.** ✅ **shipped, v0.5.0** as `console/console.py`
(`CONSOLE_VERSION 0.1.0`, its own number and not `__version__`). Project list
over N results directories, stage table from the state file, log tail by byte
offset, heartbeat with "last output N s ago". No config editing, no launching.
It can be pointed at a running three-day job on day one without touching it.
*If this pane is not useful, stop here and the loss is one file.*

> **Built:** the "without touching it" clause turned out to be the whole
> design and was hardened accordingly — see the note under "What it is not".
> `--root` scanning for results directories (three levels, re-scanned every 30
> s) was not in this milestone as written and is in it now, because eight
> datasets on one machine is the case the console exists for. What the
> heartbeat pane says is bounded on purpose: how long since the last stamp,
> what it is reading, and the `ps` line — never a verdict.

**M2 — preflight.** *Not in v0.5.0.* `requirements()` as a checklist with a server-side
directory picker and metaannot's own validators run in place. Three row types:
OK, MISS, and **MANUAL** — MANUAL rows carry a link-out and structurally no
button, because a GUI cannot accept a licence on your behalf.

**M3 — the ID-overlap gate.** One button runs `run --only emapper` — minutes,
not hours — and shows FASTA vs eggNOG vs quant identifier coverage side by
side, offering the `emapper_id_transform` or `emapper_strip_id_prefix` that
metaannot's own diagnostic names. `CLAUDE.md` calls identifier mismatch the most
common failure, and an unnoticed one produces a 95% dark bin that reads as
biology.

**M4 — config authoring.** ~14 project fields plus the 21-row stage table, not
189 inputs. Overlay-and-generate: hold only what differs from `DEFAULT_CONFIG`
and re-attach the rationale comments that `yaml.safe_dump` currently discards.
A config the console did not author is **read-only until adopted**. Diff shown
before any write, timestamped backup, refuse to write if the file changed on
disk since load — because the author also edits in vim.

**M5 — launch and supervise.** Launch into tmux on the pipeline host. Never as
a child of the console.

**M6 — the degraded reader.** The same file under Windows Python over
`\\wsl.localhost`, read-only. Exists for one real failure on this hardware:
after Modern Standby resume, new `wsl.exe` sessions die while the distro and
the job keep running. A console living inside WSL is unreachable in exactly
that state, so it cannot be the thing that tells you what is happening — and it
must refuse to offer `wsl --shutdown`, which would destroy the run.

---

## The hard parts, honestly

**The Windows→WSL bridge is not solved, it is concentrated.** Running inside
WSL genuinely deletes path translation, drvfs, 9p and cross-boundary process
launching *after startup*. What it cannot delete is bootstrap: "can the browser
reach the page" has several independent failure sources that all present
identically as a blank tab, to the user least equipped to debug one.

**Surviving three days is architecturally free and operationally not.** The
console holds no state, so closing it costs nothing. The cost is elsewhere: the
run must be parented inside Linux, because `nohup … &` launched via `wsl.exe`
dies when the call returns.

**The degraded reader is a second Windows code path exercised only during a
failure that appears unpredictably** — which is the definition of code that is
wrong when you finally need it. Accept that trade with open eyes.

---

## Decide before the first commit

1. **Local-only, or remote too?** If three-day runs belong on the Fedora server
   and WSL is the development path — which this machine's crash history argues
   — then M5 should move down and a server story should move up. A far-side
   agent is a *second product*; name it as one now or not at all.
2. **Who are the users beyond this lab?** A console assuming a working WSL
   distro and 200 GB of databases is right for the lab and useless for a
   stranger evaluating metaannot on a bare laptop. Those are different
   products. Pick one.
3. **Separate repo, or `console/` inside metaannot?** Inside is assumed, because
   one repo is what makes the contract test possible — but a console bug will
   then look like an engine bug in a user's report.
4. **Does the single-file philosophy extend to the console?** Assumed yes: one
   stdlib-only file, deployable by `scp`.

> **Built:** 3 and 4 went the way they were assumed to. `console/console.py` is
> one stdlib-only file inside this repository, and `tests/test_console_contract.py`
> is what that bought. The mitigation for the risk named in 3 is the banner and
> the error page, both of which say the console is the console: a render failure
> prints "Nothing was written to any results directory" and points at the
> console's own stderr, so a bug here does not read as an engine bug.
>
> 1 and 2 are still open, and M1 was deliberately scoped so that they could
> stay open. A watcher is right for the lab server and for a stranger's laptop
> alike, and it commits to neither answer; M5 (launch and supervise) is where
> the local-or-remote question has to be settled, and M2/M4 are where "which
> users" starts to decide what the screens say.
