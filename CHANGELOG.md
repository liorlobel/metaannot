# Changelog

## v0.5.0 — 2026-09-10

One question, asked in two halves: what can be known about a run that is still
going, from somewhere other than the tmux window that started it. Until now the
answer was `tail -f` and a state file parsed by hand, on a machine with the
results mounted — and eight datasets on one server is where that stops scaling.

The engine half (#24) puts the run's identity, its liveness and the
configuration it actually used on disk where a reader can find them, and makes
`describe --json` a contract rather than a convenience. The other half (#25) is
`console/console.py`: one stdlib-only file that reads exactly those files and
nothing else, and writes nothing at all. (A line count stood here and was stale
before the release was cut, because the fixes below landed in that same file. A
number nobody can keep true is worse than none.)

Two things were measured on the way. A signal landing between the lock's
`os.open(O_EXCL)` and the `atexit` registration that removes it stranded a lock
in **8 of 30** attempts — each one printing `interrupted.` and exiting cleanly,
so nothing said anything had gone wrong. And `describe --json` turned out to
have been a release behind without anyone noticing: `cost` has been on every
stage since v0.4.0 and never reached the contract, so every stage published
`cost: null` to whatever was reading it.

### Added

**A console: `console/console.py`, milestone M1 — a watcher.** It renders the
three files a results directory already writes — `.metaannot_state.json`, the
log and `.metaannot.lock` — as one page, for every project on the machine at
once, with a stage table, a log tail read by byte offset, and what the run
record says about itself. Python 3.9 and the standard library, so there is
nothing to install and `scp` is a deployment.

**It never writes a byte into a results directory**, and that is the property
everything else is arranged around rather than a side effect. Not a lockfile,
not a cache, not a temp file, not a log line: every reader opens `O_RDONLY`,
the console never `chdir()`s into a watched directory so not even a core dump
can land there, there is no `do_POST`, and the flock that stops two consoles
fighting over one socket lives in the console's own runtime directory — a
`--socket` inside or under a watched tree is refused, because binding there
creates both the socket and a lock file that outlives the process. That is what
makes pointing it at a job already three days in a non-decision, and it is
proved rather than intended: `tests/test_console_contract.py` checks it by AST
scan and by snapshotting a results directory around every route.

**It never imports metaannot either.** Stage order, dependency edges, the
`_run` key and the names of the files it polls all come from one
`describe --json`, shelled out at startup and cached. A front end that
hard-codes what a stage is called or where a state file lives drifts the day
either changes, and the drift is invisible until someone reads a stale answer
off a page that looks authoritative. Not importing is also what lets one file
be copied to a machine where the engine lives at a path it was never told
about; `--metaannot` and `--python` are how it is told.

**It binds a mode-0700 UNIX socket, never a TCP port.** A `127.0.0.1` listener
on a shared lab server is reachable by every other account on that box, and the
alternative — a random token — leaks into `ps`, into shell history and into the
URL bar, and would be a second access-control mechanism to keep correct. File
permissions are the whole of it, which is why the socket's directory must be
private and is checked for that, and why the console will not create one for
you. Reached with `ssh -N -L 8080:<socket> <host>`; the banner prints the exact
line for the socket it bound. The cost is documented: a UNIX-socket forward
target needs OpenSSH 6.7 or newer on the workstation side.

**And it will not tell you a run is dead.** The heartbeat is advisory — one
failed write ends the heartbeat thread while the run carries on — and the
engine's rule is that unprovable means alive, so a console that guessed would
be inviting `--force-unlock` on a live run, which is the one thing that
corrupts a results directory. It reports how long since the last sign of work,
names what it read, hands over the `ps` line for the pid in the lock, and
stops. Where the evidence is genuinely informative it says exactly that much: a
`FATAL` as the newest log line while the record still says `running` is called
the killed-mid-write shape, not death. There is no button that acts — no
`--force-unlock`, no launch, no config edit, no `do_POST`; the handler serves
`GET` and `HEAD`.

`CONSOLE_VERSION` is `0.1.0` and is deliberately not `__version__`. One
repository, two programs, two audiences: tying the numbers together would mean
either bumping one nobody asked about or lying about the other.

**The preflight checklist, the server-side directory picker and `doctor --json`
are M2 and are not in this release.** Neither is config authoring (M4) or
launching (M5). `docs/gui-design.md` describes all of them and is annotated
with what M1 actually settled; nothing there should be looked for in this tag.

**A `_run` record in the state file, with an advisory heartbeat.** The state
file recorded what each stage did and nothing about the run, so a results
directory opened cold — one of eight on a shared machine — could not say what
produced it. `_run` carries `run_id`, `version`, `config_path`, `argv`, `host`,
`pid`, `started`, `finished`, `final_status`, `heartbeat_s` and a `last_seen`
stamped every `heartbeat_s` seconds (30 by default; `0` turns it off).
`last_seen` is written twice, as a local-time string and as an epoch float,
because two hosts sharing one filesystem cannot subtract each other's local
clocks.

It is advisory and nothing reclaims a lock on the strength of it. A heartbeat
that stopped is not a process that stopped — one failed write ends the timer,
not the run — and from another host the two are indistinguishable, so acting on
it would trade a stale lock, which is a message and `--force-unlock`, for two
runs writing one directory, which is silent corruption. The heartbeat exists to
give a person the number; the decision stays theirs.

**`results/config.effective.yaml`, written under the lock.** The results kept
the stage records but not the settings behind them, and the config file beside
them is whatever it says today rather than what it said in March — nor does it
carry the defaults nobody wrote down, or the `--threads`/`--ram`/`--faa` the
command line added. It is written next to the `_run` record and under the same
lock, because a run must not be able to die between them: written apart, an
early fatal error left a `config.effective.yaml` describing a run the state
file had never heard of, and on a resume one that flatly contradicted the
previous `_run` — "finished ok, threads 7, against a FASTA that does not
exist", none of which had happened.

It is deliberately not a stage output and appears in no stage's inputs. A file
in the results root rewritten on every run would, the moment any stage listed
it, invalidate that stage on every run — thirty-four hours of InterProScan
spent recording what the run already knew.

**`describe --json` became a contract, with a `DESCRIBE_VERSION`.** Anything
reading that JSON is not in `metaannot.py` and cannot be fixed in the same
commit as the engine, so it needs to be able to say "I do not understand this
shape" instead of guessing. The document now carries the config vocabulary a
form generator would otherwise have to infer — `path_keys`, `db_path_keys`,
`replace_blocks`, `freeform_keys`, `retired_keys` — plus `stage_names`, `bins`,
`quant_formats`, `run_key`, and `paths`: the files a watcher polls, named
rather than reconstructed, so the day one of them moves the watcher moves with
it. Those paths are absolute even with no `--config`, since a watcher that
stored a relative one would poll whatever directory the describing process
happened to be in.

`DESCRIBE_VERSION` is bumped when a key is removed or its meaning changes, and
not when one is added, because a reader that breaks on an unknown key was going
to break anyway. It is its own number, like `SIGNATURE_VERSION`: a tool-version
bump must not invalidate a stage cache, and a schema bump must not either.

**`describe --json` now emits each stage's `cost`.** It was added to `STAGES`
in v0.4.0 — the 1/2/3 rank the scheduler sorts each round's ready set by — and
`describe` went on emitting its hardcoded seven fields, so every stage came
back `cost: null` and a front end could not have known the field existed.
Nothing failed, which is exactly the problem: the projection is a list of
names. It is pinned now from both directions by a test that compares the two
sets and fails if `STAGES` carries a field `describe` does not emit.
`DESCRIBE_VERSION` is deliberately not bumped for it — this is an added field.

### Fixed

**A signal between taking the lock and arming its release stranded it, 8 times
in 30.** `os.open(O_EXCL)` and the `atexit.register` that removes the lock used
to have the whole of `__enter__` between them. A signal landing in that window
left a fully written lock file behind with no hook to remove it — and the run
exited *cleanly*, printing `interrupted.`, so nothing anywhere said a lock had
been stranded. The next run then refused to start. The hook is armed **before**
`__enter__` is called now, and `__exit__` on a lock that was never taken is a
no-op, so arming it first costs nothing. The lock is also held from the instant
the file exists rather than from the instant it is written, which is the other
half of the same window: `__exit__` was a no-op while `self.held` was still
False, and the file is written within microseconds so the zero-byte grace never
applied either.

**`ResultsLock.is_still_ours()` has three answers, because two was one too
few.** A killed run can still be inside a stage when the operator decides it
has hung, `--force-unlock`s the directory and starts a replacement; whatever
the old run writes next must land on neither the new run's `_run` record nor
its lock. One gate answers that for both writers — but they agree only on two
of the three cases. "Ours" and "somebody else's" are the same answer to both.
"No lock at all" is not: there is nothing for `__exit__` to remove, while for
`RunRecord` a vacant path means nobody was superseded and its own final verdict
is still worth writing. Collapsing vacant into somebody-else's cost an
unsuperseded run its own verdict — `_run` stranded at `running` under a WARN
announcing a handover that never happened. `True` is ours, `False` is somebody
else's, `None` is vacant, and each caller reads it for itself. It is keyed on
the lock file's *content*, because a lock changes hands by remove-then-create
and a filesystem may hand the same inode straight back.

**The rest of this section was found by auditing the merge itself.** Two
branches that each passed their own tests met in one tree, and what follows is
what reading the result turned up. None of it was reported by a user; all of it
would have been.

**The tool told an operator that a DIAMOND identity floor was not in effect
while it was applying it.** `diamond_min_pidents` arrived in v0.4.0 mirroring
`diamond_evalues` key for key, and the hand-written list of free-form config
blocks was not extended to cover it — so `unknown_keys()` called
`diamond_min_pidents.mydb` a typo. Nothing was dropped, and that is the harm
rather than the mitigation: `report_unknown_keys()` only **logs**, and
`load_config()` then `deep_merge`s the user's config whatever it found, so the
floor was in force the whole time — `diamond_min_pident_for(cfg, "mydb")`
returns the user's number — while the log said of that very key "It is being
ignored, so this setting is NOT in effect". Told that, the one setting to
suspect when a database's hits look wrong is the one setting the operator has
been assured cannot be responsible. `doctor` then counted the same key as a
failure (`ok = False`) and exited **non-zero over a config that was valid and
working**, which stops any wrapper or CI job that gates on `doctor`. The list is
derived from the defaults now: anything named `diamond_<something>` whose
default is a dict has database tags for keys by construction —
`diamond_weights`, `diamond_evalues` and `diamond_min_pidents` today, with
`diamond_workers` excluded because it is an int — which makes the next one
free-form on the day it is added rather than on the day somebody notices.
`vfdb_category_weights` is deliberately still named by hand: its keys are
VFDB's own category codes, not database tags, so it is free-form for a
different reason.

**TMbed's ETA under-promised the wait on every resume, by the whole ratio of
resumed to new work.** The rate was the elapsed time divided by every residue
accounted for, and a resumed chunk contributes its residues without costing any
time — so a run resuming 90 chunks of 100 reported a tenth of the real time
remaining, and reported it before it had predicted a single residue. That is
the one number in this stage an operator acts on: it decides whether they wait,
go home, or kill the run. The rate is now divided by what *this process*
predicted, and there is no estimate at all until this process has finished a
chunk of its own — one line without a number, rather than a number that is
wrong.

**A resumed TMbed chunk was adopted by record count, so retuning the chunk size
duplicated and dropped predictions at once.** `tmbed_chunk_residues` is outside
every stage signature on purpose — it changes the order of the records and not
one prediction in them — which makes tuning it on a resume a sanctioned
operation. It also **re-plans** which protein goes in which chunk. Chunk 3 of
the new plan then covers a different set of proteins from the chunk 3 whose
file is on disk, and a count-only test adopted that file whenever it happened
to be long enough: the concatenation wrote two records for every protein in
both plans and none for the proteins in neither, which is a topology set that
is silently short *and* silently duplicated — precisely what `tmbed_failed.tsv`
and `tmbed_allow_partial` exist to make loud. Adoption is by identifier now,
and by equality rather than containment, and a chunk that does not match is
named in a WARN and discarded rather than left where the salvage path would
find it later. The question is also asked once, for every chunk, before
anything is written or run: it now has side effects — a warning and a removal —
and asking it again inside the loop would repeat the warning and re-examine a
file it had already taken away.

**`tier_coverage.tsv is absent for that reason` was said without looking at the
file.** Results directories get re-run: last month's config tiered and this
month's does not — `exclude_id_prefixes` grew to cover every namespace, or the
proteome was swapped for one whose ids carry no source prefix — and the
previous run's `tier_coverage.tsv` is then still sitting beside results that
are new, with the log asserting it is not there. Someone opening the directory
cold reads a table describing a proteome this run never had. The message now
looks first, and when the file is there it says so, with its size and its date,
at WARN rather than INFO. It is not deleted: this is the one file that is
deliberately outside `finalise`'s declared outputs, so nothing else would
notice it going, and metaannot does not delete what it did not just write.

**An unreadable lock file stranded the lock it was being read to release.**
`is_still_ours()` already treated bytes it could not *read* as "still ours",
because an I/O error is not evidence that a lock changed hands. Bytes it could
not *decode* took a different exit: `UnicodeDecodeError` is a `ValueError` and
not an `OSError`, so the exception left through `__exit__`, which had already
cleared `held` — stranding exactly the lock the release path exists to remove.
One truncated multi-byte character from a torn NFS write, or a page of nulls
where a crashed writer's payload should be, was enough. Undecodable content is
damage rather than a rival's payload, since everything this file writes is
ASCII JSON, so it gets the same answer for the same reason. A garbled but
decodable file is a different question and still reads as not ours.

**A superseded run could write its stale state over a replacement that had
already finished.** The ownership gate was asked afresh on every write, so the
verdict depended on what the replacement happened to be doing at that instant —
and the losing sequence needs no race at all. Run A is `--force-unlock`ed and
logs `this run no longer holds ...`; run B finishes and removes its own lock on
the way out; A's final stamp then finds the path **vacant**, which is `None`
rather than `False`, takes the branch that exists for a replacement which has
not started yet, and writes A's whole stale in-memory state dict over B's
completed results. The three-state gate is right and the branch's reasoning is
right; what was missing is that a run told once that it lost the directory has
lost it for good. The flag is latched now, and read before the lock is asked. A
run that was never superseded is untouched, which is the case that branch is
there for.

**The supersession of the SIGTERM handler was undocumented in the one place it
had to be.** `handle_sigterm_like_sigint()` still described a `SIGTERM` that
unwinds and stamps `_run`, which is what it does for `init`, `doctor`,
`describe`, `subset`, `report`, `object` and `run --dry-run` — and is not what
happens to a `run` or an `all`, because `cmd_run` registers a handler over it
the instant it has the lock, and the later registration wins. Both docstrings
now say which subcommands each handler covers and why the later one takes
precedence, and a test pins the order, so reversing it fails loudly instead of
quietly changing what `kill` does to a three-day run.

**The console's log pane grew without bound on the page it is meant to be left
open on.** Every line `/api/log` delivered became a `<div>` that stayed, and the
whole design is a tab left open for the length of a three-day run — so a stage
emitting a few lines a second put a quarter of a million nodes in one scroller
and the tab's memory climbed until the browser killed it. The pane keeps the
last 2,000 lines now, ten times what the first render shows, and it says so:
a pane holding a tail looks exactly like a pane holding the whole log, and a
reader who believes the second concludes a stage never logged something it
logged four hours ago. The cap is the console's, served to the page as
`data-max` so that how much of a three-day log this program will hold is
decided and read in one place; a page cached from a console that served no
`data-max` bounds itself at the same 2,000, since an unbounded pane is the one
outcome that is not acceptable. Trimming takes lines off the top, above a reader
who has scrolled up, so the scroll position is given back exactly what the
shortening took — and the "N earlier lines dropped" marker is cleared, not
carried, when the pane is emptied because the log rotated.

**The console ranked its NEXT rows in table order over a queue the engine takes
longest-first.** v0.4.0 made the scheduler sort each round's ready set by
`cost`, and `describe --json` emits that field "so a front end can order or
annotate the table the same way" — which this console then dropped. On the
455,571-protein run the page listed dbcan NEXT above ncbifam NEXT, the engine
dispatched ncbifam, and dbcan did not get a worker for 27 hours. The same run
that produced the scheduling fix produced this reading of it. The rows are
annotated rather than reordered — the table is in the engine's stage order and
a dependency graph read out of order is harder, not easier — and the claim is
narrowed to what a watcher can actually establish: `nothing is blocking it`
became `no dependency is blocking it`, with the queue said separately as a
**ranking** and nothing more. A NEXT row names up to four of the ready stages
ranked ahead of it and counts the rest, because on a fresh run every
dependency-free stage is ready at once and a row that names nine of them is a
directory listing where a sentence was wanted. What the console cannot see goes
under the table once rather than onto every row: how many stages start together
depends on `stage_workers` and on what is still running, and a stage that needs
the GPU waits for `gpu_workers` however it is ranked — `gpu_lease` defers it and
leaves it in `remaining` to be reconsidered next round. The ranking is refused
whole rather than done partially: if any ready row lacks a usable numeric
`cost` — an engine older than the field, a stage added to a newer one without
it, or a `cost` arriving as a bool, a string, a null or an infinity — no row is
ranked and the note under the table is not printed either. A partial sort over a
missing key would rank a stage by a number this file made up, and the whole
point of taking the order from `describe --json` is that it is the engine's
order.

**One poll of an `esmfold` run was a 455,000-entry directory walk, per open
tab.** The part-file cache is keyed on the output directory's mtime, which is
exactly right for a directory holding a handful of declared outputs and exactly
wrong for the one that costs the most: `esmfold`'s declared output is
`results/structures/.done` and the stage writes one `.pdb` per dark protein into
that same directory, so the mtime moved between every poll and the key never
hit. A part file is evidence and not a control, so the mtime key now has a floor
under it: a directory is rescanned only if its mtime has changed **and** 30
seconds have passed since the last scan, and in between the last answer stands,
re-`stat()`ed so the size it reports is still live. `time.monotonic()` rather
than `time()`, because an NTP step backwards over a lab server's first hour
would otherwise park the floor for as long as the step. Noticing a part file
half a minute late costs a sentence on one poll; the walk cost the box the run
was on.

**A dependency turned off in the config was reported as blocking, if it carried
an old record.** `decide()` returns `disabled (run.X)` for any stage whose flag
is falsy whatever the state file holds, `finish()` counts it as skipped, and
`unmet_deps()` passes it — so the engine walks straight past. The console asked
the wrong set: a dependency disabled in this config but carrying last week's
`failed`, or a killed run's `running`, made its dependents read WAIT, "waiting
on foldseek (failed)", over a queue that was not waiting for anything. Whether
the engine waits is decided by the config alone; a record only decides whether
that dependency's own row may read OFF, and those are two questions now.

**A `running` record from an earlier run was rendered as this run's live
stage.** `mark_running()` writes `{signature, status, started}` and no
`finished` at all — which is the point of it — so the era check compared `None`
against the run's start time, fell through to "this run", and drew a stage with
a duration counting up, for ever, from a clock that stopped days ago. That is
the exact shape a `kill` leaves behind, so the case is common rather than
exotic. A record is dated by its `finished` as before, and a record still marked
`running` — the one kind that has no `finished` at all — is dated by its
`started` instead. One that predates this run gets its own STALE state: not the running colour at full
strength, no duration, and a sentence saying the record was written by an
earlier run and is not evidence that anything is running now. It stops there, as
every verdict on this page does; what it adds is the engine's own reading, which
is that a run reaching a stage in this state recomputes it.

**"The engine carries on with every stage that does not depend on a casualty"
was wrong, in two places on the page.** The dispatch loop is
`while (remaining or futures) and not failure`: nothing new starts after the
first failure. What it then does is wait for the stages already running, which
it cannot interrupt — and an InterProScan that started an hour before the
failure has hours left in it. The observable is the same, a run that looks busy
long after a stage was lost, and the reason a reader would infer from the old
sentence is the opposite of the real one. Both the failure block and the index's
grouping rule now say what actually happens.

### Documentation

**Both documents described a `SIGTERM` that does not happen.** README and
TUTORIAL said a killed run takes Ctrl-C's path — "the same `interrupted`
message, the same `_run` stamp" — and that "the stage already running has to
finish first, so the process does not exit immediately". That was true of the
first version of the handler and of nothing since. Both passages are rewritten
against the code: `run` and `all` release the lock, write one line to stderr and
`os._exit`, and the consequence the README had backwards is now the right way
round — a `_run` record still reading `running` is the **ordinary** trace of a
`kill`, not evidence of a `SIGKILL` or a lost machine. Both now say what to do
with a directory in that state, including the one thing worth checking first:
the handler stops metaannot, not the InterProScan or DIAMOND it launched, and
those keep writing with no lock left to keep a second writer out.

**The test counts were off by roughly 490.** They are the only yardstick a
reader has for judging whether a checkout is sound, and a fifth of the suite
missing from them reads as a failed checkout rather than as a stale document.
Recounted, with the R tests explained as what they are — in the default run,
skipping rather than failing where R is absent — and pinned by a test that
compares each figure against a real collection of this suite. The Windows
paragraph was stale in both directions and is rewritten from the code: `doctor
--fix` is a refusal there rather than a recipe that emits `mkdir -p`, the signal
tests are skips rather than failures, and stale-lock reclamation is no longer a
Windows exception at all. The console's tests have not been run on Windows and
the paragraph says only that.

**`docs/gui-design.md` said the console did not exist.** It opened "Status:
**design, not built.** Nothing here exists yet" in a tree that ships it, and
described a front end serving to `127.0.0.1` — the one decision the build
reversed. It stays a design record rather than becoming a manual: the reasoning
is what a later milestone has to argue against. Where what shipped settled a
question differently it now carries a **Built:** note saying so, and where a
decision was made and honoured it stands as written.

**The console is documented where a reader would look for it.** A new README
section covers what it is, the three things it refuses to do and why each makes
it safe, the socket, the SSH forward and every flag; TUTORIAL phase 5b offers it
beside `tail -f` in the place an operator actually meets the question; and
CLAUDE.md gains a standing rule and an entry in the deliberate list, because the
refusals are the design and an agent reading that file first should not try to
help by adding a write path.

### Tests

The suite is 1154 passed, 1 skipped, 6 xfailed and 37 deselected on a default
run, in three to five minutes. `-m slow` is those 37; `-m R` selects 35 that the
default run already includes and that skip cleanly where `Rscript` is absent.

Six xfails are open, over five tests, and **not one of them is new.** All five
markers stand exactly as they did at the v0.4.0 release commit (`b846753`): one
in `test_config.py`, one in `test_parsers.py`, one in `test_scheduler.py` and
two in `test_taxonomy.py`, the `test_parsers.py` one parametrized `empty` /
`header_only` and so worth two of the six items. The one that reads newest — an
`emapper.annotations` file with no data rows raising a bare `KeyError` instead
of the "no `#query` header" message every other malformed file gets, reachable
by adopting a zero-row eggNOG table — has carried that marker, with that
parametrization, since `7e665e5` (#1); the file it lives in was last touched in
`2e26ca5` (#11). Both are ancestors of the v0.4.0 release. What changed at this release is only the count in the
README, which said **8** against a tree that has had 6 since v0.4.0.

The console arrives with two test files of its own. `test_console_contract.py`
is the one that matters most: it proves by AST scan and by snapshotting a
results directory around every route that nothing is written into one, and it
pins `describe --json` from both directions — including the check that failed on
`cost` and would fail again for the next field added to `STAGES` and not
emitted.

## v0.4.0 — 2026-09-10

Scheduling, resumability and thresholds, all of it forced by one 455,571-protein
run rather than by review. Two measurements set the agenda: InterProScan started
**27.1 hours in**, because it sits tenth in a table whose first wave included a
ten-minute stage; and TMbed wrote **nothing at all for 39 hours** — no progress
bar, no partial file — after two earlier attempts died at 2 h 36 min with
nothing recoverable.

Six defects were fixed along the way, two of which only appeared under test: a
signal handler that could hang a run instead of releasing its lock, and a
structure-resume guard placed where it could never be reached. The suite is
green on Windows for the first time (735 tests), and on every CI job from
Python 3.9 to 3.13.

### Added

**The scheduler dispatches longest-first, not in table order.** Every stage
with no dependencies is ready in the first round, so the first `stage_workers`
of them IN TABLE ORDER started and the rest queued. The 455,571-protein run
shows what that costs: with `stage_workers: 3` the first wave was emapper,
pfam and dbcan, and InterProScan — the longest stage in the pipeline — did not
start until **27.1 hours in**, when kofam finally freed a worker. dbcan, which
takes ten minutes, held one of the three from the first second.

Each stage now carries a coarse cost rank (3 = hours, 2 = minutes, 1 =
seconds) taken off two real runs rather than intuition, and each round's ready
set is sorted by it. Longest-processing-time-first is the standard greedy
answer and it costs one sort of a list that is never longer than 21. The rank
reaches no cache signature and no `keys` list — scheduling order cannot change
a stage's output — and `stage_priority` indexes rather than defaults, so a
stage added without a cost raises instead of silently ranking as trivial.

**TMbed runs in resumable, length-sorted chunks.** TMbed writes nothing until
it finishes — not "buffers a bit", nothing — so a single invocation over a
whole proteome was an all-or-nothing bet measured in days. On the
455,571-protein run it produced no progress bar and no partial file for 31
hours; the two attempts before it died at 2 h 36 min with nothing recoverable.

The input is now grouped into chunks of about `tmbed_chunk_residues`
(5,000,000 by default, roughly 17k average proteins) and each is committed as
it lands, so an interrupted run resumes from the last finished chunk. Under 5M
residues there is exactly one chunk and the behaviour is unchanged; `0`
disables the split.

Chunks are length-sorted longest-first for two reasons: ProtT5 pads every
sequence in a batch out to the longest one in it, so a chunk of similar
lengths wastes less work; and whatever is going to exhaust the device is then
in the FIRST chunk, where it costs one chunk to discover instead of the whole
stage. The count is bounded at both ends — the budget is a floor, 256 parts a
ceiling — since ProtT5 is loaded once per chunk. The floor adds the longest
sequence to `total / max_parts`, because greedy packing closes a chunk when
the next sequence would overflow it and can fall short by that much; without
the term, 5,000 × 100 residues planned 264 parts against a ceiling of 256.

A chunk that dies keeps whatever TMbed wrote. What reaches the committed file
is reconciled against what was handed over — from the files, because a chunk
can exit 0 and still come back short — and anything missing is named in
`results/topology/tmbed_failed.tsv`. `tmbed_allow_partial` (default false)
decides whether that is fatal, and `tmbed_max_consecutive_failures` (2) stops
a wedged card from failing every remaining chunk the same way, slowly. The
stage is now `empty_ok`: a proteome whose every sequence is over
`tmbed_max_len` leaves an empty prediction file on purpose and runs nothing.

`tmbed_chunk_residues` is deliberately not part of the signature — it changes
the order of the records and not one prediction in them, and listing it would
discard a 30-hour stage because someone tuned a checkpoint size. The two keys
that change what is in the file are listed.

**A per-database identity floor for DIAMOND: 50% for CARD, VFDB and BAGEL.**
Hits from those three are read as claims about a particular protein, not as a
family assignment, and the pipeline reported them down to the 30% global
default. On the 455,571-protein run that was 39,138 CARD hits and 107,219 VFDB
hits; above 50% they are 4,661 and 21,623. A 32%-identity match to a
beta-lactamase over half a query is a hit against the fold, not evidence that
the protein confers resistance.

`diamond_min_pidents` mirrors `diamond_evalues` and is applied twice on
purpose: as DIAMOND's `--id` during the search, so a floored database writes
thousands of rows rather than hundreds of thousands, and again in the reader,
because a `<tag>.tsv` adopted from elsewhere never saw `--id`. It also kills
the half-weight rule for a floored database — nothing below
`diamond_strong_pident` survives the filter — which the run now says once,
rather than leaving the config implying a grading that cannot happen.

**The identifier key, not just the source tag.** A tier tag says which SOURCE
a protein came from; the identifier key says what the id actually is, and it
lives behind the tag:

    uhgpL_MGYG000004906_01237   tier uhgpL_     key MGYG#_#
    uhgpSM_MGYG000009567_01280  tier uhgpSM_    key MGYG#_#
    OIDECCNN_00158              tier OIDECCNN_  key #
    ampS_AMP10.000_478          tier ampS_      key AMP#.#_#

Four tiers over three key spaces. `uhgpL_`, `uhgpSM_` and the `ent_`
entrapment set all wrap the same MGnify `MGYG` namespace — 31,782,562 of the
search database's 36,637,426 records — while `OIDECCNN_` is the cohort's own
Prokka run. Reporting only the tag hid that two tiers are one namespace under
two labels, which is precisely the case `emapper_strip_id_prefix` exists for:
one row of a precomputed table annotates a protein under every tag it carries,
and a tag left out of that list does not error — its whole tier reports as
unannotated, which reads as biology.

`tier_coverage.tsv` gains `key_shape` and `key_shape_pct`; the run names any
key shared by more than one tier; and `prepare_emapper` warns, while it can
still be acted on, when one sharing tier is listed and another is not. Digit
RUNS are masked rather than digits, because the Prokka tier runs
`OIDECCNN_00001` to `OIDECCNN_1712297` — 5-, 6- and 7-digit accessions — and
masking per digit would split one key space into three.

**A tier that matches nothing is a zero, not a low number.** The AMPSphere
tier of the 455,571-protein run matched 0 of its 2,168 proteins — no eggNOG
table on that machine is keyed on AMP/SPHERE ids — and it is 30% of the whole
dark fraction. The headline coverage was 98.4% (`uhgpL_` 98.94%, `OIDECCNN_`
98.74%, `uhgpSM_` 97.23%), so nothing said so, and every one of those proteins
bins `4_dark` for want of a join rather than for want of biology. A tier at
0% now gets its own line naming the two ways out.

**Coverage per identifier tier.** A merged search database is the normal case
and its tiers do not annotate alike: the run this was built on is `uhgpL_`
(395,467), `OIDECCNN_` (43,139), `uhgpSM_` (14,797) and `ampS_` (2,168), one
of which arrives with precomputed annotations and one of which is ORFs nobody
has ever seen. `finalise` writes `results/tier_coverage.tsv` — count, share,
percentage carrying each kind of evidence, percentage dark, median export
score — and `emapper` logs and records the same split for its own coverage.

`exclude_id_prefixes` is applied before the split. `proteins_faa` is not
always the identified subset, and run the whole FragPipe search database
through and its two LARGEST namespaces are the controls: 18,318,713 `rev_` and
2,280,823 `ent_` against 11,379,230 `uhgpL_`. A tier table whose top row is
the decoy set is not a description of the biology, and those namespaces would
also spend the twelve-tier budget on entries that are there to be ignored.
What is dropped is counted out loud rather than silently.

Tiers are detected from the identifier prefix, and declined in the two cases
where a prefix is not a source label: one prefix over everything, and more
than twelve. When declined, no file is written **and the run says why**, so an
absent table is never ambiguous. A column absent from the frame is skipped
rather than reported as 0%, because "the stage did not run" and "the stage
found nothing" must not share a cell.

**`doctor` recognises a motif seed set.** The BAGEL database this pipeline
built and searched — 262 entries, mean 15 residues, exactly 0 hits against
38,204 proteins — was not a bacteriocin sequence database. BAGEL4 ships two
things that look alike on disk, and what got built was the motif SEED set its
HMM step uses. The existing check answered that with "lower your `--evalue`",
which is wrong twice over: it sends the reader to tune a threshold on the
wrong kind of file, and the tuned search then produces meaningless hits
instead of meaningless silence.

Two signals, either sufficient: a typical sequence under 25 residues, and
seed-set markers in the headers of a source FASTA beside the `.dmnd`. 25 and
not 40, because mature nisin is 34 residues and a genuinely short bacteriocin
database must not be accused of this. The message replaces the e-value advice
rather than joining it; that check is not gone, only narrower.

### Fixed

**A structure interrupted mid-write counted as folded, for ever.** The resume
rule for `esmfold` is "the `.pdb` is there, so it is folded", and each
structure was written in place. `open()` truncates immediately, so a run
interrupted between the open and the flush — a wedged card, a bugcheck, a
Ctrl-C, all of which this stage has seen — left a 0-byte file that every later
run counted. The protein was then permanently absent from Foldseek with
nothing anywhere to say why.

Structures go through `atomic_out` now, whose temp is dot-prefixed and keeps
its extension, so `glob("*.pdb")` never sees a file something is still
writing. For the files already on disk, `folded_already()` checks for a
coordinate line — and not only in the fold loop: the early exit that decides
whether the GPU is touched at all asks the same question first, so a guard in
the loop alone would have been unreachable.

**The results lock survived SIGTERM and SIGHUP.** `atexit` does not run on
either — Python's default handler terminates the process outright — so a run
stopped by `kill`, by a scheduler hitting its time limit, or by a closing ssh
session left a lock naming a pid that no longer exists. On the same host the
next run can prove it is dead; from another node it cannot, and the resume
became a stale-lock refusal needing `--force-unlock`.

The handler releases the lock and exits 128+N. Nothing in it may take a
lock: a Python signal handler runs IN THE MAIN THREAD, between two bytecodes
of whatever that thread was doing, so any lock the interrupted frame holds is
still held and is not reentrant. The first version called `log()`, which goes
through `sys.stderr`, whose buffer lock is exactly that -- on a run emitting a
progress line per stage per minute the signal eventually lands mid-write and
the process HANGS instead of releasing the lock, which is strictly worse than
the stale lock this exists to prevent. CI caught it on one job of seven. The
message is pre-formatted and pre-encoded at registration and goes out through
`os.write(2, ...)`; the cost is that the final line reaches stderr but not the
log file, whose buffer cannot be safely touched from a handler. It does not
stop the tools already running — they are separate processes that outlive us —
and it uses `os._exit`, because `SystemExit` would unwind through the stage
pool's `with`, which waits for its workers, and a tmbed chunk can be an hour.
SIGINT stays unhandled on purpose: Python raises `KeyboardInterrupt` for it,
which unwinds the lock's `with` and runs `atexit`.

Writing the test found a live bug in the handler itself. `cmd_run` bound
`lock` twice — the `ResultsLock` at the top and a `threading.Lock` 230 lines
below — and a closure captures the name, so the handler called
`threading.Lock.__exit__`, died with "release unlocked lock", and took the run
out with an uncaught `RuntimeError` and exit 1 instead of releasing anything.
The mutex is `state_lock` now, and both the handler and `atexit` hold the
bound method rather than the name.

**`doctor --fix` ran POSIX shell through cmd.exe.** Every command it generates
is `mkdir -p` / `curl` / `tar` / `gunzip` / `hmmpress`, and
`subprocess(shell=True)` on Windows hands those to cmd.exe, where
`mkdir -p C:\db` creates a directory called `-p` and each step can exit 0
having done nothing — precisely the failure the post-install verification
exists to catch, happening for every item at once. It refuses there now and
names the alternative: `--install-plan` here, `wsl bash install.sh` where the
tools live. `--install-plan` itself still works on Windows.

**The results lock could never be reclaimed on Windows.** `_holder_is_alive`
returned True unconditionally there, because `os.kill(pid, 0)` calls
`TerminateProcess` and asking whether a process is alive that way would kill
it — but that made a crashed run permanently unresumable without
`--force-unlock`, and a crash is exactly when reclaiming has to work.
`OpenProcess` asks without touching the process. Unprovable still means alive:
only `ERROR_INVALID_PARAMETER`, the answer for a pid that does not exist at
all, is taken as proof of death.

### Documentation

The resource guide now carries the second run's measured numbers and says
plainly which three stages had not finished when it was written, rather than
rounding an unfinished stage into the table. The linear-scaling claim is
replaced with the measured spread: 11.9× the proteins gave 14–30× the time on
the stages that finished, with MMseqs2 the exception at 7.6×. A new table
gives the start time of every stage on that run, which is where the 27-hour
InterProScan delay is visible.

Corrected throughout: the Pittsburgh proteome is **455,571 proteins /
214.9M residues**, counted off the FASTA. Several comments written during this
work said 1.3M, which was a recollection rather than a measurement. Also
corrected in three places: BAGEL's 262 entries of median length 15 are its
motif seed set, not "262 bacteriocin sequences", and the documented
`stage_workers` default is 4.

### Tests

The suite is green on Windows for the first time: 723 passed, 0 failed. Three
failures there were platform assumptions rather than defects — a path test
comparing POSIX strings, and two tests using signals Windows cannot deliver to
a child — and two were the real Windows defects fixed above. The signal tests
are skipped on Windows and verified under Linux, where they pass.

### Added

**`doctor` reports the GPU.** It had no GPU check at all, so a machine with no
usable CUDA device turned on `run.structure`, was told "all checks passed",
waited hours, and learned the truth when `esmfold` finally ran and exited. With
`run.structure` or `run.topology` on, a `== gpu ==` block now says what the
stages will find.

`cuda_probe()` separates the two cases that look identical and are not: no card
at all, and a perfectly good card with a CPU-only `torch` wheel. The second is
the more confusing failure because the hardware is right there, so it is named
with its fix. It does not import torch when torch is absent — `doctor` has to
stay fast and has to run on machines that have none.

`structure` without CUDA is a MISS, because `stage_esmfold` exits rather than
fold on CPU; the message names the way out, which is to fold on a GPU host and
copy `results/structures/` back, since `foldseek` is CPU-only and searches
whatever models are there. `topology` is a WARN and not a MISS: SignalP 6 is
CPU-only and unaffected, and tmbed does run without a GPU — just one to two
orders of magnitude slower.

`stage_tmbed` now names the protein count before starting a CPU fallback.
"Fell back to CPU" is a footnote at 5,000 proteins and a two-day decision at
455,000, and tmbed writes nothing until it finishes, so that path cannot be
told apart from a hang while it is running.

### Fixed

**RAM detection answered 0 on everything without `/proc`.** `detect_ram_gb()`
tried `os.sysconf("SC_PHYS_PAGES")`, then `/proc/meminfo`, then gave up. macOS
*defines* `_SC_PHYS_PAGES` but `sysconf` returns EINVAL for it, and Darwin has
no `/proc`, so both probes fell through — and the caller reads 0 as "could not
detect", so the memory budget silently became zero and every stage ran with no
allocation: no DIAMOND `-b`, no `-Xmx` for InterProScan, no
`--split-memory-limit` for MMseqs2 or Foldseek. Nothing errored; the run was
just quietly unbudgeted. Windows had neither probe and returned 0 too, which
matters because `doctor`, `report` and `object` run there directly.

Added `hw.memsize` for Darwin and `GlobalMemoryStatusEx` for Windows, after the
existing probes rather than instead of them. Verified: 127 GB on Windows, 94 GB
in WSL2. The first probe also now checks `sysconf` returned a positive number
before trusting it. Invisible until now because CI is `ubuntu-latest` only.


**Foldseek's legacy-column fallback could never run.**
`stage_foldseek` asks for `qtmscore` and `qlen`
and falls back to the ten legacy columns when a build does not have them. It
caught `StageError`. `run_cmd` raises a plain `RuntimeError` on a non-zero
exit, and `StageError` is a **subclass** of `RuntimeError` — so the handler
could not catch the one failure it exists for. A Foldseek 5 build lost its
structure evidence outright instead of degrading. The test that covered it
injected `ma.StageError("Invalid selection: qtmscore")`, a message no tool
emits, so it passed over dead code.

`stage_foldseek` is the only site with this shape: the other eight
`except StageError` handlers are correct, because nothing on their paths
reaches `run_cmd`.

**A rejected format code costs seconds, not the search.** `easy-search`
validates `--format-output` in `getOutputFormat` (EasyStructureSearch.cpp
line 42) and does not create its temporary directory until line 59, so
foldseek exits before it prefilters anything and leaves no tree behind. The
retry is therefore cheap, and the gate is permissive in the direction that
matters: a rejection whose wording it cannot parse still falls back rather
than losing the structural evidence to a changed message.

**Only a rejected format code is retried, and only one the fallback can
drop.** `FOLDSEEK_COLS_LEGACY` is a strict subset of `FOLDSEEK_COLS`, so
falling back can only ever remove `qlen`, `tlen`, `qtmscore` and `ttmscore`.
The other ten fields are in **both** lists — a build that rejects `lddt` or
`theader` fails the retry identically, and the second error is then the one
the operator has to explain. Those re-raise with a line saying why. Both
halves of the phrase must appear on **one line** of the stderr tail:
`<path> does not exist` is stock MMseqs2 wording for a missing database, and
tested across a multi-line tail it would splice an unrelated line onto the
words "format code" and read a path as a rejected column.

**The scratch tree is now removed on every exit, not only success.** The
cleanup sat after the loop body, so every re-raise — a full disk, an OOM
kill, the undroppable-column path — left it behind. Against AFDB50 that is
tens to hundreds of GB, and `CLAUDE.md` already records
`results/foldseek/tmp*` as something that is never cleaned up. It is a
`try/finally` now.

**The warning was wrong about `qlen`.** Foldseek 5 and earlier accept `qlen`;
only `qtmscore` and `ttmscore` are missing. The message said "no
qtmscore/qlen" and now names the right two columns and the version boundary.

**Re-weighting VFDB was a silent no-op on a re-run.** `vfdb_category_weights`
was in `finalise`'s signature keys but not in `integrate`'s. `build_annotation`
is what applies the weighting, and `integrate` is what runs `build_annotation`
— so editing the map invalidated only the stage that could not act on it.
`integrate` stayed cached with the old scores; `finalise` re-ran, found no
structure or profile evidence, took its `reusing the first pass` branch, and
copied the stale `annotation_pass1.tsv` through verbatim. `annotation_final.tsv`
and `bin_summary.tsv` came out byte-identical and the run reported success.

That branch is the common case, not an edge: it is taken whenever
`run.structure`, `run.hhblits` and `run.jackhmmer` are all off, which is the
default and is what all eight configs under `examples/server-run-plan/` set.
The v0.3.0 upgrade itself was unaffected only by luck — `integrate`'s signature
changed anyway that release, because its input list lost `p.effectors`.

`foldseek_target_priority` was missing from `integrate` for the same reason and
is added with it. It reaches `build_annotation` on any re-run where a previous
Foldseek result is already on disk.

The rule is now stated where it can be checked: **every config key
`build_annotation` reads unconditionally must appear in the signature keys of
both stages that run it.** The three keys read only inside its `emit_dark`
branch — `exclude_id_prefixes`, `max_dark_structures`, `max_len_structure` —
belong to `integrate` alone, because `finalise` never writes `dark.faa`. A new
test derives that set from the source rather than listing it, so a key added to
`build_annotation` later cannot quietly skip the signature.

**Cost on an existing results directory.** Adding keys changes the signature
payload, so `integrate` re-runs once for everybody, and `finalise` with it.
Both work in process. No search stage is upstream-invalidated, and `esmfold`
takes `dark.faa` as its only input and hashes it by content — so an unchanged
work list keeps every fold. No InterProScan, KOfam, ESMFold or Foldseek compute
is discarded.

**`Ntox` could not fire, in the pattern list added to fix exactly that.**
v0.3.0 added `CdiA`, `LXG`, `Ntox`, `nuclease toxin` and `zeta toxin` to
`toxin_fold_patterns` because the shipped list had matched 0 of 38,204 real
proteins — the whole-word rule meaning `Tc toxin` could not match inside
`holotoxin`. `Ntox` was added with the same defect it was added to repair.
Every family in that set is `Ntox` followed by a number — Ntox15, Ntox28,
Ntox47 — and a digit is a word character, so the trailing `\b` meant a bare
`Ntox` matched none of them. It matched only the string `Ntox` standing alone,
which is not how the family is ever written. The pattern is now `Ntox\d*`.

The test that was meant to prove it worked passed anyway, because its one case
— `"Ntox47 nuclease toxin domain"` — is also matched by the neighbouring
`nuclease toxin` pattern. Each case in that test now names the single pattern
it exercises, and the test re-compiles the list **without** that pattern and
asserts the description stops matching, so a case carried by a neighbour fails
instead of reading as a pass.

Patterns in `toxin_fold_patterns` are joined into one alternation and are
therefore regexes rather than literals. That was always true and is now said
in the config comment, since the default list contains a metacharacter for the
first time.

`toxin_fold_patterns` is in the `finalise` signature keys, so `finalise`
re-runs on the next run and nothing else recomputes.

## v0.3.0 — 2026-09-08

Isobaric quantification, a cross-source agreement check, and the removal of a
stage that asked the user to do the work somewhere else. Validated on two real
datasets: the UC antibiotics label-free run end to end (38,204 proteins, all
fifteen applicable stages, report and R object built from it), and Pittsburgh
ICB-melanoma FragPipe TMT (8 plexes, 455,571 proteins).

`SIGNATURE_VERSION` stays at **1**. Nothing here changes what a *search* stage's
output means, so a resumed run keeps its InterProScan, KOfam and ESMFold
compute — which is the whole reason that number is not tied to `__version__`.

**Upgrading from 0.2.0.** Three things to know:

- **Re-run `foldseek` if you have results from 0.2.0.** The stage now asks for
  `qtmscore` and `qlen`, and without them the TM gate silently falls back to
  `alntmscore` — normalised by the alignment rather than the query, so a short
  local match inside a long protein can pass. `finalise` detects the old
  10-column table and says so, and unlike in 0.2.0 the advice now works. On the
  UC run the stricter gate moved `3s_structure_only` from 246 to 141: **43% of
  those calls were artefacts of the weaker gate.**
- `finalise` will re-run by itself, because `toxin_fold_patterns` changed.
  That is cheap and is what you want — the VFDB re-weighting and the repaired
  `toxin_fold` both land there.
- A config carrying `run.effectors`, `effector_predictions` or
  `effector_prediction_weight` still loads. Those keys are named as removed,
  with the reason, and `doctor` exits zero on them.

### Added

**FragPipe TMT (isobaric) input: `quant_format: fragpipe_tmt`.** `quant_table`
becomes the FragPipe run directory holding the per-plex `TMTn/` folders, and
the reader takes the reporter channels from each plex's `ion.tsv` or
`peptide.tsv` — never from `tmt-report/`, whose matrices are already log2 and
median-centred, are protein level (which deletes the shared-peptide rule,
`peptide_assignment`, `peptide_evidence.tsv` and the peptide assay of the R
object) and carry TMT-Integrator's own protein inference. Per-plex reporter
intensities are linear, so the existing log2 path, roll-up and
median-of-ratios size factor apply unchanged, and the reader emits exactly the
shape the label-free readers emit, plus a design carrying `sample`, `plex` and
`channel`.

The new `tmt:` block holds `plex_glob`, `level`, `annotation`,
`reference_name`, `reference_channel`, `use_reference_ratios`,
`condition_from_name`, `within_plex_normalise`, `min_plexes`,
`drop_empty_channels` and `min_purity`. Two things the real files forced:
the reference channel is resolved **per plex** (a pool sits at `131C` in six
plexes of a real run and at `131N` in the other two, so a single
`reference_channel` cannot describe it — `reference_name: "Pool*"` globs the
sample name instead), and the annotation file is `<PLEX>_annotation.txt`, not
`annotation.txt`.

Plexes are joined on a feature id comparable across them, and a feature not
identified in a plex is `NA` for every sample of that plex, never `0`.
Refusals: reporter columns that cannot be mapped to samples, two plexes
claiming one sample name, a plex with no annotation file, a missing per-plex
table (never a fallback to `tmt-report/`), and — under every `quant_format` —
a `tmt-report/` matrix or the TMT flavour of `msstats.csv`, each named for what
it is.

**A TMT design where the plex is a batch and not the hypothesis.**
`design_from_input.tsv` gains a `plex` column, and for `fragpipe_tmt` the
default `analysis.design_formula` becomes `~ 0 + group + plex` with
`factor_cols: "group,plex"`, so the contrasts come out over the condition
rather than over the plex. A formula the user wrote is left alone, and a
single-plex run keeps `~ 0 + group`. The report says, where it prints the
coefficients, that `plex` is a batch term; a plex perfectly confounded with the
condition now stops the run in `plex_confounding()` — before the rank check,
whose own message only ever named a coefficient — with the cross-tabulation
and the plex named, and the rank-deficiency message names `plex` too when it is
in the model.

The condition comes from the annotated sample names when they carry it
unambiguously (`tmt.condition_from_name`, default `auto`: every name must split
on the same separator into levels of at least two samples, and no other
separator may group them differently) and otherwise from `analysis.metadata`.
When it can be had neither way the run stops and prints the metadata file to
write, row by row, starting from `design_from_input.tsv`. It is never taken
from the plex.

**The reference channel is never a modelled sample.** Under the default
covariate treatment it is dropped from the sample columns and the plex carries
the batch; under `use_reference_ratios: true` every channel is divided by its
plex's reference. The ratio treatment propagates the reference's missingness —
a feature with no reference in a plex is `NA` for that whole plex — so the
reader counts and reports the values that costs (0.53% on the real 8-plex run)
and warns when it is large. Which treatment was used, where the condition came
from, and the reference channel of each plex are written to
`results/quant/design_notes.txt` and copied into `design_record.txt`.

**Within-plex normalisation, on by default.** `tmt.within_plex_normalise:
median` divides each channel by its own median and multiplies by the plex's
median channel, before the roll-up. The channels of one plex are the same
LC-MS run, so what differs between them is loading and labelling efficiency —
a per-channel constant the roll-up would otherwise sum into the protein — and
because the median commutes with log2 this is exactly a per-channel median
centring, done one plex at a time. The **between**-plex difference is left
alone on purpose: that is the batch, and the plex term (or the report's own
`normalise: median`) removes it. `none` keeps FragPipe's numbers and warns.
The choice, with the log2 factors applied, is written into `design_notes.txt`
and from there into `design_record.txt`.

**`min_purity` is implemented rather than refused.** FragPipe writes precursor
purity only into `psm.tsv`, so `tmt.min_purity` now reads that file per plex,
keys it on the columns the feature id was built from, and judges a feature by
the **median** purity of its PSMs — its reporter intensities are a sum over
them, so no single PSM's purity describes it. `0` (the default) does not read
psm.tsv at all, and psm.tsv is a stage input only when the filter is on. A
feature matching no PSM row is **kept** and counted, because an unmatched key
is a join failure (~1.5% of ion keys in a real plex), not a dirty spectrum.
This is an approximation of the per-PSM filter TMT-Integrator applies before
summarising, and is reported as one.

**A protein-level `min_plexes`, applied beside `min_valid_per_group`.**
`analysis.min_plexes` counts plexes where `min_valid_per_group` counts
samples: an isobaric run's missingness is plex-shaped, so a per-group sample
count can be satisfied entirely inside one batch. The report now applies both,
counts each against the same starting set and **prints them separately**, with
the histogram of proteins by number of plexes; with the filter at its default
of 1 it still reports how many retained proteins live in a single plex, which
is the number to set it on. It is inert without a plex column, so label-free
behaviour is unchanged, and asking for it where no per-sample plex exists
stops the report instead of passing everything. `min_plexes` joins
`design_record.txt`.

**The report prints the isobaric design it used.** A new section states, from
`design_notes.txt` verbatim, the input shape, where the condition came from,
which reference treatment was applied, each plex's reference channel and the
within-plex normalisation, then prints samples per plex and condition-by-plex
and names the plex a covariate — a TMT batch and never a condition.

**Ratio compression is stated, not corrected.** The report's limitations
section gains an isobaric-only entry: reporter ratios are pulled toward 1 by
co-isolated precursors, so log2 fold changes are lower bounds, the shrinkage
is worst exactly where a strain-redundant database puts near-identical
paralogues in the same isolation window, direction and ranking are more
trustworthy than magnitude, and `min_lfc` is stricter here than on label-free
data. metaannot applies no interference model and no purity-weighted
rescaling: every such correction divides by an estimate of the contamination
and turns a known bias into an unknown variance.

**The taxon size factor was checked for a plex leak, and reports its
exposure.** On a fixture with no biology, run twice with the same random draw
and a 7.5x loading spread between plexes, the plex effect lands entirely in
the per-sample part every taxon shares (1.575 log2 recovered against a true
1.585) and the taxon-by-taxon part is identical to 9.4e-05 log2 — so a plex
effect does not reach a taxon's ratio model. Plex-shaped **missingness** does:
a taxon with too few members observed in every plex loses the complete-case
reference, falls back to poscounts and mixes references computed inside one
plex with references spanning all of them (1.40 log2 of displacement in the
fixture). The join stage now counts those taxa and warns, naming both
`min_plexes` keys.

**A per-plex `protein.tsv` is refused instead of being read as one plex.** It
is the only isobaric file that carries neither of the markers the other
refusals key on — no `ReferenceIntensity`, no `Channel <mass>` — so
`quant_format: fragpipe` took every numeric column it had and quantified a
SINGLE plex as the whole experiment, with `Length`, `Protein Qvalue` and
`Razor Intensity` in the matrix beside the eleven channels, silently. It is
now refused naming the reporter columns it found. The discriminator is the
prefix form (`Intensity <sample>` for a reporter channel against
`<sample> Intensity` for a label-free run), so a label-free
`combined_protein.tsv` cannot match it and label-free behaviour is unchanged.

**The isobaric path is checked against planted numbers, not against "it
runs".** A two-plex fixture with a 2x condition effect (half up, half down), a
3x plex loading and a deliberately UNBALANCED condition split — 3 'a' + 1 'b'
in one plex and the reverse in the other, so modelling the plex and ignoring
it give different answers — is read, rolled up and fitted. With the plex in
the model the planted effect comes back at ±1.00 log2 and the unregulated
proteins at 0.00; without it every protein gains half the batch (+0.79 log2)
and the proteins that truly went **down** are reported as barely moving. The
same is asserted through real `limma` on the report's own formula, and the
ratio and covariate routes are shown to agree on the planted effect to within
0.10 log2.

**A running tool now says it is running: `progress_interval_s` (default 60,
`0` disables).** Every interval, the newest line the running tool wrote to
stderr is echoed under its stage tag with the elapsed time
(`tmbed running 2h36m | 61%|###### | 23310/38204 …`), and a tool that writes
nothing gets `no output yet on stderr` on the same schedule. On the run this
comes from, tmbed ran 2 h 36 min and then died — twice — and InterProScan for
2.9 h, with the log saying nothing at all between the command and the failure,
because `run_cmd` captured stderr only to quote it on failure. The only way to
tell either apart from a hang was to watch its CPU ticks accumulate in
`/proc`, and the tutorial tells people to expect 1–3 **day** runs.

Nothing is hoarded: stdout still goes to `/dev/null` — InterProScan emits tens
of MB of chatter and buffering it bought nothing — and stderr is now read as
it arrives into a 200-line ring, which is all the failure tail (last 15 lines)
ever quoted. `run_cmd`'s return value and its `RuntimeError` are unchanged;
every stage depends on both. Carriage-return redraws, which are what `tqdm`
emits, come out as one readable line rather than a wall of partial redraws:
CSI escapes and control characters are stripped and the line is capped at 160
characters. Stderr is also decoded with `errors="replace"` now, so a tool that
writes latin-1 there cannot turn its own failure into a `UnicodeDecodeError`.

**The GPU is scheduled: stages that need it take an exclusive lease
(`gpu_workers`, default 1).** Independent stages run concurrently and the CPU
and RAM budgets are split between them; the GPU was not modelled at all.
`tmbed` held 15.5 GB of a 16 GB card and ESMFold peaked at 13.3 GB on a single
short sequence, so on the run this comes from they only avoided each other by
accident — they happened to be in separate invocations. Enable `run.topology`
and `run.structure` together on a one-GPU machine and they contend, and the
failure looks like an unexplained CUDA OOM in whichever loses.

`tmbed` and `esmfold` are now marked `gpu=True` in the stage table, and the
scheduler dispatches at most `gpu_workers` of them at a time. Nothing else is
serialised: CPU-only stages are dispatched in the same round as before, and a
deferred GPU stage waits in the ready list rather than occupying a worker. A
stage that is waiting logs it once — naming the stage holding the card and the
`gpu_workers` value — so an enabled stage that has not started is explained
instead of looking stuck.

`gpu_workers` is a lease count, **not** a device map: every GPU stage is pinned
to the single `gpu_device`, so raising it on a two-card machine runs two stages
on the same card rather than one per card. One stage per device is not
implemented, and the config comment, the README and the startup log all say so
rather than presenting 1 as a physical limit.

**A per-database DIAMOND e-value: `diamond_evalues` (empty by default).**
Overrides `thresholds.diamond_evalue` for one tag, in the search and in the
filter `integrate` applies to the hit table — filtering the table back down to
the global threshold would have left the setting doing nothing. Each database
searched at a threshold other than the headline one says so in the log.

**`results/source_agreement.tsv`, and an `ncbifam_accs` column.** Coverage
overlap and concordance are different questions, and the tool answered only the
first. `finalise` now compares every pair of columns that shares an identifier
namespace — Pfam by hmmsearch against Pfam inside InterProScan, the same for
NCBIfam, eggNOG KOs against KOfamScan's, and Pfam names against eggNOG's —
counting identical, overlapping and disjoint calls and recording which side is
the superset when they differ. The report gains a section for it.

The NCBIfam comparison was impossible before, because `ncbifam_hits` kept only
family NAMES while InterProScan reports ACCESSIONS: comparing them scored 97.2%
"conflict" between two searches of one library. The accession was in the
`hmmsearch` tblout all along and was simply discarded; keeping it as
`ncbifam_accs` turns that into 100% agreement over the same 12,677 proteins.
A near-total disjoint rate is now reported as the namespace mismatch it almost
always is, rather than as a finding.

**`esmfold_vram_cap`, `esmfold_bytes_per_residue_pair`, `esmfold_vram_reserve_gb`.**
The fold work-list is now capped by the VRAM actually free once the weights are
resident, not only by the fixed `max_len_structure`. ESMFold does not raise an
out-of-memory error when it stops fitting — the driver pages device memory to
host RAM and the fold just becomes one to two orders of magnitude slower, so
the failure is silent and the run stops being finishable rather than stopping.
Measured over 1,819 folds on a 16 GB card: 22.1 s median at 470-478 aa, then
140 s at 481 aa and 2,053 s by 491 aa, for a 0.6% increase in length. Sequences
over the cap are reported as never attempted, with the log naming which of the
two limits bound and what to do with them. See the README for the formula and
the calibration.

**`esmfold_allow_partial` and `esmfold_max_consecutive_failures`.** What a
long fold run does when the card fails rather than the protein. See Fixed
below; both default to the old behaviour of stopping rather than going on with
a reduced structure set.

**`doctor` now has a `== tmt ==` block.** Under `quant_format: fragpipe_tmt`
it lists the plexes it discovered, checks each has the `ion.tsv`/`peptide.tsv`
its `tmt.level` names and a readable `<PLEX>_annotation.txt`, reports the
channel count per plex, warns when a sample name appears in more than one, and
checks that `tmt.reference_name` or `tmt.reference_channel` resolves in every
plex. Before this the only TMT line `doctor` printed lived inside
`== manifest ==`, and a TMT config normally sets no manifest — so the one
layout this format is most particular about went unchecked until the run died
on it.

### Removed

**The `effectors` stage, and with it `run.effectors`, `effector_predictions`
and `effector_prediction_weight`.** It predicted nothing. It read result files
the user had to obtain first from Bastion3/4/6, EffectiveDB, T4SEpp or
SecretomeP — web services this tool has no way to invoke — so the stage's real
requirement was that you go and do the work somewhere else and come back.

Every one of those services is now unreachable. Probed from two networks:
Bastion3/4/6 never complete a TCP handshake, BastionHub returns a 503,
`effectors.csb.univie.ac.at` is NXDOMAIN, EffectiveDB's own redirect target
404s, and SecretomeP is retired by DTU with only the 2004 *mammalian* standalone
ever distributed. Scale would have killed the online route independently:
455,000 proteins against a 100-sequence form is 4,550 scripted submissions to an
unpaid academic server that grants no permission for it.

Replacing it with an offline predictor was considered and rejected on the
evidence rather than the effort. Every published model in this field is trained
on 138–509 proteobacterial *pathogen* proteins — T4SEpp's 509 positives are 303
*Legionella* and 135 *Coxiella*, and sixteen common gut genera appear zero times
in its positives **and** its negatives, so a *Bacteroides* protein is outside
both classes at once. The world's entire supply of validated T6SS effectors is
331 sequences, four of them *Bacteroides*. On a gut metaproteome those numbers
would have entered the score as a confident float with no way to see it going
wrong.

**Removing it is a provable no-op for every run performed to date.**
`run.effectors` defaulted to `False` and `effector_predictions` to `{}`, and the
ingest was gated on the table having more than one column, so on both paths the
prediction columns were empty and the scoring loop added nothing. An older
config carrying any of the three keys now gets a line saying the setting was
removed and why — not a spelling suggestion, which would send you hunting for a
typo you did not make. `doctor` reports them as a WARN and still exits zero.

A `pred_*` column you join in yourself still reaches the candidate shortlist as
a column. It deliberately no longer feeds the score, which is the right status
for somebody else's model.

### Changed

**`effector_score` is now `export_score`.** Every term in it — signal peptide
class, β-barrel, LPXTG or SLH anchor, CAZy, small size, toxin-like fold,
mobile-element or secretion neighbourhood, no KO — asks whether a protein
*leaves the cell*, not whether it is an effector of a secretion system. Nothing
in the tool asks the second question any more, so the old name promised a claim
the evidence never supported. `effector_score` is still written, with identical
values, for one release. `bin_summary.tsv` gains `median_export_score`.

### Fixed

**VFDB hits are weighted by VFDB's own category.** A flat `+4` threw away the
one thing the hit already told you. On a real gut metaproteome, of 3,308 VFDB
hits the two largest categories were *Immune modulation* (965) and
*Nutritional/Metabolic factor* (903 — GroEL, ClpP, GuaA, LPS biosynthesis),
each collecting the largest DIAMOND weight in the config, while the categories
that actually name an exported effector — VFC0086 effector delivery (232) and
VFC0235 exotoxin (132) — were 11% of the signal. `vfdb_category_weights` is
keyed on the stable numeric code, because VFDB can reword a category name and
will not renumber it; an unlisted code falls back to the flat weight, and the
parse rate is logged once so a format change surfaces as a line rather than as
silence. All 3,308 hits parsed. Score mass from VFDB falls to 52%.

**`toxin_fold` could not fire.** It carries the joint-largest weight in the
score and matched **0 of 38,204** proteins on a real run. The whole-word rule
means the shipped `Tc toxin` pattern cannot match inside `holotoxin`, so PDB
2vse — a genuine Tc-family holotoxin, and the one real hit in that dataset — was
missed. `holotoxin` is now listed explicitly rather than loosening the existing
pattern, and the contact-dependent and T6SS families gut commensals actually
carry (`CdiA`, `LXG`, `Ntox`, `nuclease toxin`, `zeta toxin`) are added. Every
pattern stays whole-word: bare `deaminase` and `hemolysin` already burned this
once. The regex now also reads `hh_desc`, which is free and lifts the ceiling
above the proteins that got a Foldseek hit at all.


**A fold that fails is no longer allowed to destroy a stage that is nearly
done.** On a 1,912-protein dark set, ESMFold reached 1,785 and then died with
`CUDA driver error: device not ready`, taking the whole run with it. The cause
was not memory: a long sequence at a large `esmfold_chunk_size` runs a single
attention kernel for long enough that the display driver resets the device
under it, and that arrives as a plain `RuntimeError`, not an
`OutOfMemoryError`. The retry path tested for OOM specifically and re-raised
everything else, so the one fault it could actually have recovered from was
the one it did not catch.

Both faults now get the same treatment, because both answer to the same
remedy: retry once at half the chunk size, and skip a sequence that fails
twice rather than abandon the stage. A skipped protein is written to
`results/structures/esmfold_failed.tsv` with its length and the error, so the
shortfall message in `build_annotation` can now say which missing proteins
were attempted and failed and which were added to `dark.faa` after the last
fold — it used to have to name both causes in one hedged sentence, because
nothing recorded which had happened.

A card that has stopped responding altogether fails every sequence, and
walking the rest of the list takes hours to produce nothing, so after
`esmfold_max_consecutive_failures` (default 5) in a row the stage stops and
says how many structures it has. Every `.pdb` is written as it is folded, so
rerunning resumes. `esmfold_allow_partial: true` instead finishes with what
folded and lets Foldseek search it, which is a real reduction in evidence and
says so in the log.

**"Supply it in `analysis.metadata`" is no longer said to people who
already have.** When a TMT design carries no condition — sample names like
`MF0030` hold none, and the plex is a batch, so there is nothing to fall back
on — the warning advised setting `analysis.metadata`. It gave that advice
whether or not `analysis.metadata` was set and already named every sample, and
on a real run where it was, the line read as a defect rather than as the
division of labour it is: `design_from_input.tsv` records what the INPUT knew,
and the metadata is merged at report time. The warning now looks: it says the
metadata names every sample and this affects only that one file, or names the
samples it leaves out, or says the path does not exist, or says the
`analysis.sample_col` column is missing.

**Fold progress reports a rate and an ETA, not just a count.** Folding time
rises steeply with length and the work-list is sorted shortest first, so a
count says nothing about how long the rest will take. When a real run slowed
by roughly a factor of ten partway through its tail, the log gave no sign of
it and the slowdown had to be read off file timestamps. The line now carries
the current length, the chunk size in use, seconds per protein and hours
remaining, and comes every 10 rather than every 25 — at the long end of a real
dark set, 25 proteins is hours between lines.

**A finished fold stage no longer loads the weights.** It used to upload ~11 GB
to the card and only then discover that every sequence was already folded or
excluded by `max_len_structure`. The static limit and the already-on-disk check
now run first, so a resumed or complete stage returns without touching the GPU.
On a machine where the GPU is reached through a virtualisation layer that
upload is itself a risk, not merely waste.

**Tools are launched by the absolute path `PATH` resolves them to.** On
Windows, `subprocess` without a shell calls `CreateProcess`, which does its
own `PATH` search and only ever appends `.exe` — it ignores `PATHEXT`. So a
`.cmd` or `.bat` earlier on `PATH` lost to an `.exe` later on it, and the tool
that ran was not the tool `shutil.which` reported. `run_cmd` and the two
direct `subprocess.run` callers now resolve the name first, so the
availability check, the logged command line and the process that starts all
agree. On POSIX this changes nothing.
*(Corrected after v0.3.0. This paragraph said "the two direct
`subprocess.run` callers", but `stage_tmbed` was a third and was not resolving
its name — see the Unreleased entry, which routes it through `run_cmd`
instead. `_run_rscript` remains outside this, deliberately: it needs the
stderr `run_cmd` discards on success. Nothing else in the v0.3.0 entry has
been changed.)*

This surfaced as a test-suite failure, which is the more useful half: the
`stub_bin` fixture writes extension-less Python scripts, so on Windows it was
entirely inert and 18 tests silently exercised whatever real tools happened to
be installed. The fixture now writes a `.cmd` shim beside each stub.

**A DIAMOND database that cannot possibly hit is no longer searched in
silence.** Two real failures, one after the other, and both produce 0 hits —
which in `annotation_final.tsv` is indistinguishable from a real absence of
virulence factors.

A `diamond makedb` that had failed left a **zero-byte** `.dmnd` on disk. It
passed every existence check, and the stage would have searched it. `doctor`
now reports it as `MISS` and `stage_diamond` refuses before any search starts,
naming the file, its size and the `diamond makedb` line that rebuilds it.
Anything under 128 bytes is smaller than DIAMOND's own header (a one-sequence,
11-residue database is 143 bytes), and a file that `diamond dbinfo` reports as
holding no sequences is refused on the same grounds.

The BAGEL bacteriocin database built correctly — 262 sequences, **median length
15 residues** — and returned exactly 0 hits against 38,204 proteins at the
pipeline default of `--evalue 1e-10`. That is not a finding about the biology:
using DIAMOND's own BLOSUM62 constants (Lambda 0.267, K 0.041), the best
e-value a *perfect* 15-residue alignment can reach is about 1e-5, so the search
was incapable of a hit before it started. The run now warns, with that number,
and says what to change — a per-database `diamond_evalues` entry, or the
knowledge that this database needs different settings than the rest. It also
names the scoring weight the database is holding while it cannot hit
(`diamond_weights` `bagel: 3`, equal to TADB3), because a weight claims the
database can contribute to the effector ranking.

The typical sequence length is read from `diamond dbinfo` (its Sequences and
Letters, so a mean); when diamond is not installed the check falls back to a
source FASTA beside the database or named in `sources.diamond` (a median), and
when neither is available it says the length is unknown rather than guessing
one and acting on the guess.

**A log line can no longer kill a multi-hour run on Windows.** Tool output is
decoded with `errors="replace"`, which puts U+FFFD into descriptions, and
printing U+FFFD to a console in the Windows code page (cp1252 here) raises
`UnicodeEncodeError` — so a stage that logged a protein description containing
one replacement character could abort a run that had been going for hours, on
the log line rather than on the work. `stdout` and `stderr` are reconfigured to
UTF-8 with `errors="replace"` at startup where the platform allows it, and
every log write falls back to replacing what the stream cannot encode, for the
streams that are not ours to reconfigure. The log **file** was already UTF-8;
it was the console that was the hazard.

**"the rest were skipped (OOM)" is no longer said about structures that were
never folded.** The first structure run printed

> `WARN 0/1913 requested structures exist in …/results/structures; the rest
> were skipped (OOM) or never folded`

which is a correct sentence in the wrong situation: the directory was empty
because nothing had been folded **yet**, and a reader reasonably concluded the
run had already lost 1,913 models to the OOM killer. The three causes are now
separate messages. No `.done` marker, no `plddt.tsv` and no models is INFO —
the esmfold stage has not run, nothing has been lost, and the pass simply
carries no structural evidence. Models but no `.done` marker says the rest are
pending rather than lost. Only a finished esmfold blames OOM, and even then
only for the models that length does not already explain: proteins over
`max_len_structure` are excluded before ESMFold ever sees them and are counted
separately in the same line. The count is also now taken over the requested
proteins present in the table rather than over every `.pdb` on disk, so models
left from an earlier, larger `dark.faa` can no longer mask a real shortfall.

**`peptide_evidence.tsv` no longer crashes when nothing is assigned.** With
every feature shared across taxa — which is what a strain-redundant metagenome
database really produces, and the case `peptide_assignment: taxon_unique`
exists for — the evidence table's index came from the dropped side of an outer
join and was named `razor_protein`, so the rename to `protein_id` raised a bare
`KeyError('protein_id')` instead of writing an empty-but-well-formed table.
The index is now named before it is reset. Label-free and isobaric alike.

### Documentation

**README's "Limitation: FragPipe TMT output is not supported" is gone, and
what replaced it says which designs work.** The section now opens with what is
read (`TMTn/ion.tsv` or `TMTn/peptide.tsv`, `TMTn/<PLEX>_annotation.txt`, and
`TMTn/psm.tsv` only when `tmt.min_purity` is set) against what is refused by
name under every `quant_format` (the eight `tmt-report/` matrices, the TMT
flavour of `MSstats.csv`, the per-plex `protein.tsv`) — the three files a user
reaches for first — and states that nothing falls back to them.

Read: **reference-free** runs with no pool anywhere; **multi-plex without a
bridge**, where the plexes are linked by the plex coefficient and the report's
median normalisation rather than by a shared channel; **mixed plex sizes**, a
16-, an 11- and a 6-channel plex read as one experiment; **TMTpro 16/18 and
iTRAQ**, since channel labels come from the annotation and no channel set is
hard-coded (real data so far is TMT-11 only); and a **single plex**, which
keeps `~ 0 + group`.

Not supported-with-a-caveat but named as unfixable: a condition that does not
cross plexes (stops the run twice over), a formula that is rank deficient once
`plex` is in it, a condition crossing only some plexes, a plex holding one
sample, a protein quantified in one plex only, and TMT effect sizes compared
with label-free ones.

**Every new config key is documented with its default**, in one table:
`tmt.plex_glob` `"TMT*"`, `tmt.level` `"ion"`, `tmt.annotation`
`"{plex}_annotation.txt"`, `tmt.reference_name` `""`, `tmt.reference_channel`
`""`, `tmt.use_reference_ratios` `false`, `tmt.condition_from_name` `"auto"`,
`tmt.within_plex_normalise` `"median"`, `tmt.min_plexes` `1`,
`tmt.drop_empty_channels` `true`, `tmt.min_purity` `0`, and
`analysis.min_plexes` `1` — plus the two whose defaults change for
`fragpipe_tmt` with more than one plex, `analysis.design_formula`
(`"~ 0 + group"` -> `"~ 0 + group + plex"`) and `analysis.factor_cols`
(`"group"` -> `"group,plex"`).

The two schema traps a user meets on their own data are stated where they bite
rather than left to a log line: **the reference channel is not at a fixed
channel position across plexes** (a real 8-plex run puts the pool at `131C` in
six and `131N` in two, which is why `reference_name` globs the sample name and
why `reference_channel: "131C"` dies naming the plex it missed), and **purity
exists only at PSM level** (`psm.tsv` alone has a `Purity` column, so
`tmt.min_purity` is a join and judges a feature by the median purity of its
PSMs).

**The three new keys, with their defaults.** `progress_interval_s` `60`
(seconds between progress lines from a running tool; `0` disables),
`gpu_workers` `1` (GPU stages in flight at once — a lease count, not a device
map) and `diamond_evalues` `{}` (per-database `--evalue`, overriding
`thresholds.diamond_evalue` for the named tag in the search *and* in the filter
`integrate` applies to the hit table). All three default to the behaviour that
was there before, so an existing `config.yaml` needs no edit and
`SIGNATURE_VERSION` is unchanged; `diamond_evalues` is in the signature keys of
`diamond`, `integrate` and `finalise`, so setting one does recompute those.
Each is refused with a message rather than a traceback when its value is not
a number.

**What is watched, and what is still not.** The progress line is the tool's
newest line of **stderr** — stdout goes to `/dev/null` as before, so a tool
that reports there gets only the heartbeat — and it is not a checkpoint:
`esmfold` and `hhblits` resume protein by protein, every stage that shells out
to one long command does not, and a stage interrupted mid-write is recomputed
rather than adopted. Stages that work in-process (`integrate`, `finalise`,
`join`) emit no progress lines at all. `"skipped (OOM)"` likewise remains an
inference from absence: a protein whose fold attempts both fail leaves no
`.pdb` and no `plddt.tsv` row, and only the run log's `skipping <id>` records
it. All of this is in README rather than left implied.

## v0.2.0 — 2026-09-07

A test suite, and the first three defects it and the first real run turned up.
`SIGNATURE_VERSION` is unchanged at 1: no stage's output means anything
different, so an existing results directory stays valid and nothing recomputes.

### Added

**`tests/` and CI, where there were none.** 419 tests, one named for each
defect found by hand during development or by the first real run, every one
carrying a one-line comment stating the symptom it protects against. The
suite runs offline with no external tool: where a stage shells out to
hmmsearch, DIAMOND or MMseqs2 the binary is a stub on PATH, so the stage's own
plumbing — atomic writes, adoption, signatures, the results lock — is
exercised without the tool. Fixtures are generated from fixed seeds; nothing
binary is committed. The default suite is about a hundred seconds and excludes
a `slow` mark carrying the resume sweep and parallel-versus-serial across
every `quant_format` x `peptide_assignment`. R tests skip cleanly where R is
absent and include a real knit where it is present. CI covers Python 3.9-3.13.

This also settles two claims the README made and nothing checked: parallel and
serial execution now provably produce identical output, and five repeated runs
produce one hash.

Eight tests are `xfail(strict)` rather than passing. Each names a guard that
is still missing — among them `parse_emapper` raising a bare `KeyError` on a
zero-row table, and a gap in a Unipept lineage truncating a protein's
consensus to the rank before it. They are marked so a fix turns them green
rather than being quietly forgotten.

**`quant/sample_columns.txt`.** `build_object.R` and the report both read it
as the authoritative statement of which columns carry intensities, and nothing
wrote it, so every run fell through to `design_from_input.tsv` — or, with no
manifest, to guessing from column types. The `join` stage now records what it
detected, renamed and filtered. Not one of the stage's declared outputs, so an
existing results directory keeps its old fallback until `join` reruns;
`--force --only join` produces the file without recomputing anything else.

### Fixed

**A non-UTF-8 byte in a tool's output is no longer fatal.** One `0xa0` — a
latin-1 non-breaking space — in a search result killed the `integrate` stage
of the first full run on real data, after InterProScan had already spent three
hours, with a message that named neither the file nor the stage. Every parser
reads through `opener()`, and all six died on it: `parse_diamond`,
`parse_interproscan`, `parse_kofam`, `parse_hmm_tblout`, `parse_hmm_lib_desc`
and `read_fasta`. VFDB subject titles, InterPro descriptions, HMM `DESC` lines
and FASTA headers all carry latin-1 in the wild. Now `errors="replace"`, the
same as `read_delim_table` beside it, so a bad byte costs one character of a
description rather than the stage.

The gzip branch of the same function passed no encoding at all, so it used the
locale's codec and was never UTF-8 by construction. `emapper_precomputed` is
routinely a `.gz`, which made it the input most likely to be read differently
on the server than on the laptop.

**The effector shortlist is not "empty by construction" without SignalP.**
`docs/signalp-6.md` and README said `surface_or_secreted` is False for
everything when the topology stage is off, so an empty shortlist was a
missing-tool artefact. Two of the gate's four terms need no topology tool: the
LPxTG sortase motif, computed from the sequence, and an anchor domain from the
`pfam` stage. On the first real run 604 proteins passed the gate with topology
off, 155 of them KO-less; the shortlist was empty because nothing was
significant — none of the 212 groups that reached the model passed FDR 0.05.
The docs now say what SignalP genuinely adds, which is every secreted protein
carrying neither motif, and that a shortlist without it is not empty but
biased toward cell-wall-anchored surface proteins.

### Known limitations

Unchanged from v0.1.0, and still the ones that decide how a result is read:
isobaric labelling is refused rather than quantified, only the first contrast
gets the full analysis, and the InterProScan memory budget is inherited by
child JVMs in distributed mode. See the v0.1.0 entry below.

Issue #5 records a further one found by the same run: the report prints
retention by bin, including a bin that contributes zero proteins to the model,
but nothing escalates it — and on that run neither `4_dark` nor `3d_duf_only`
reached the statistics at all.

## v0.1.0 — 2026-09-07

First public release. This is the version that has actually been run end to end
on a real label-free metaproteomics dataset, which no earlier state of the code
had been.

### Verified on real data

A 36-sample label-free FragPipe dataset against a matched-metagenome Prokka
database with a precomputed eggNOG table: all identified proteins resolved to
the database, 122,278 peptide features rolled up to 16,070 protein groups
across 832 taxa, the R Markdown report knitted to HTML with figures and tables,
and the QFeatures object built with both assays linked. Bin percentages fell
inside the ranges TUTORIAL.md gives as sane.

Not yet exercised on real data: the external search tools (hmmsearch, DIAMOND,
InterProScan, KOfamScan, HHblits, jackhmmer, ESMFold, Foldseek) and the Unipept
API. Their parsers and integration logic are tested, and the tools install
cleanly, but a full run with them enabled has not been completed.

### Fixed

An audit of the previous state found 211 defects. 250 fixes landed across two
passes; the counts differ because several findings needed changes in more than
one place. The ones that changed results rather than tidying code:

- **Razor identifiers never matched the FASTA.** The quantification reader
  preferred FragPipe's `Protein ID` column, which for a non-UniProt database
  holds `<id> <description>`. Identifiers are now taken from `Protein` and
  reduced to their first whitespace token, the same key the FASTA uses.
- **eggNOG's own domain calls were ignored when binning.** Only hmmsearch Pfam
  hits counted, so with the Pfam stage off, proteins eggNOG had already
  assigned a domain to were reported as having no evidence at all.
- **Fractionated manifests were rejected.** Repeated sample names are how
  FragPipe denotes fractions; they are now collapsed into one sample per group
  instead of aborting the run.
- **Taxonomy identifiers were corrupted on reuse.** A cached table was re-read
  without string dtypes, so a numeric taxid round-tripped through float and was
  written back as `821.0`.
- **A search database that prefixes its identifiers could not be joined** to an
  annotation table that does not carry the prefix, and the diagnostic then
  reported the proteins as genuinely unannotated. New `emapper_strip_id_prefix`.
- **FragPipe's zero means "not quantified"** and was being summed as a real
  zero, turning missingness into fold change. Now missing by default, under
  `zero_intensity_is_missing`.
- **Text files were written with the platform codec**, so on Windows the report
  died on a Unicode character and left an empty file. Every text file is now
  UTF-8.
- **Interrupted writers left truncated outputs** that the next run adopted as
  complete. Every stage now writes atomically.
- **The dark-protein work-list was capped for structure prediction**, and the
  profile-search stages reused that capped file, so a GPU budget silently
  truncated an unrelated search. The two lists are now separate.
- **Two configured download URLs were dead**, and one returned an HTML landing
  page that `curl -fL` accepted with exit 0, writing a web page to disk as a
  database. Downloads are now rejected if they begin with HTML.
- **Report and object plumbing**: contrasts were derived from the wrong term of
  the design formula, level names were not sanitised the way R sanitises them,
  parameters reached the document unquoted, and a failing `Rscript` did not set
  a non-zero exit status.

Silent behaviour is now reported: proteins dropped by the feature-count filter,
proteins withheld from folding by the structure cap, decoy and contaminant rows
removed, taxa whose size factor falls back to a plain sum, and which taxonomic
rank is in force.

### Added, all defaulting to existing behaviour

`rollup_method` (a log-space alternative to the plain sum), a
`taxon_or_family_unique` peptide-assignment mode, `taxon_rank`,
`ncbifam_uninformative_test`, `foldseek_target_priority`, `full_content_digest`,
`peptide_only_reader`, and `exclude_id_prefixes`.

### Known limitations

- **Isobaric labelling is not supported.** Reporter-ion channels are not read.
  FragPipe TMT output is refused with an explanatory error rather than being
  quantified from the pooled MS1 column, which is what earlier code did
  silently. See the limitation section in README.md.
- **Only the first contrast gets the full analysis.** Differential abundance is
  fitted for every contrast, but the ratio model, the effector shortlist and
  the R object cover the first one only.
- **The InterProScan memory budget is passed through `_JAVA_OPTIONS`**, which
  the JVM applies at higher precedence than a command-line `-Xmx` and which
  every child JVM inherits. That is correct for the default standalone mode,
  where InterProScan's workers are threads inside one JVM, and it is how the
  stage budget actually reaches it. If you run InterProScan in distributed
  mode, its worker JVMs are configured for 9 GB each in
  `interproscan.properties` and would instead inherit the whole stage budget,
  so divide `--ram` by the worker count there.
- **SignalP 6.0 is licence-gated** and must be installed by hand; see
  `docs/signalp-6.md`. Without it, and without a GPU for tmbed, the topology
  stage is off, so the signal-peptide and beta-barrel terms of
  `surface_or_secreted` are dead. That narrows the effector shortlist's gate
  rather than closing it: an LPxTG sortase motif or an anchor domain still
  qualifies a protein, and neither needs a topology tool.
  *(Corrected after v0.2.0. This bullet originally said the shortlist is then
  "empty by construction rather than by result", which is wrong — see the
  v0.2.0 entry above. Nothing else in the v0.1.0 entry has been changed.)*
