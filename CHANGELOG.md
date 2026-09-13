# Changelog

## Unreleased

### Added

**The heartbeat now says how far through InterProScan is.** It is the longest
stage in the pipeline — 57 h of the 455,571-protein run — and it went all of
that blind. Its stderr says it is alive and never says how far through it is,
so the only way to get an ETA was to write a chunk-counting script against its
`-T` tree by hand while the run was going. The signal was there the whole
time: InterProScan leaves a `.fasta` per chunk it splits out and a `.raw`
beside that chunk once it has been analysed, both under the temp directory
this tool already gives it. `run_cmd` now takes an optional probe a stage can
supply when it knows something the tool's own output does not carry, and
`stage_interpro` supplies one that counts those files:

```
[31840.7s] INFO    interpro | interproscan.sh running 8h50m | chunk 312/380, ~2h08m left | ...
```

The scheduler entry above ends by saying that once the dispatch order is
right, what is left of a run's wall clock is `stage_workers` and InterProScan
itself. This does not make that stage shorter. It makes it legible, which on
two and a half days of waiting is the difference between watching a log and
watching `/proc`.

**What the line refuses to say is most of the change.** Each refusal is
there because the obvious stronger claim is one the files cannot support:

- The unit is CHUNKS. Not sequences, not a percentage. InterProScan chose the
  slices, they are not equal, and the merge and write after the last one are
  not counted at all — so the line reaches `380/380` with real work left, and
  a percentage computed off it would read 100% during that tail.
- No remaining time is offered while the denominator is still moving.
  InterProScan splits while it analyses, so a rate measured then is measured
  against a number about to grow; it is perfectly computable and perfectly
  wrong. The count must hold still for `_IPS_ETA_STABLE_TICKS` consecutive
  heartbeats before a rate is anchored. The test for that is a fake clock over
  a tree where chunks are finishing AND being created, because a test where
  only the denominator moves passes whether or not the gate exists — which the
  first version of it did.
- `.raw` beside `.fasta` is THE LAYOUT OF ONE OBSERVED RUN, read off a real
  `-T` tree rather than out of InterProScan's source, and no InterProScan was
  available here to check it against. So every branch that cannot be sure
  reports nothing rather than something wrong: an unrecognised tree, a tree
  too large to walk within `_IPS_SCAN_CAP` entries, a `.raw` that outlived its
  `.fasta` (which drops the denominator and keeps the count). A probe that
  never engaged is then reported at the END of the stage, where whether it
  ever saw a chunk file is settled, rather than guessed at from a timeout
  while the tool may simply still be splitting.

A probe that RAISES is dropped for the rest of that command and said once. It
reads a directory the tool is concurrently writing, where a file vanishing
between the readdir and the stat is ordinary, and an annotation is never worth
a stage — nor worth a WARN per minute for two days.

**The report now says what `analysis.min_plexes` bought, whichever way the run
set it.** The first full run on real data was filtered at 2 and the document
never said what that decided: it printed how many protein groups the filter
removed, in a clause subordinate to a sentence about removal, and the only
counterfactual a reader could build out of the page was to cumulate the
per-N-plex histogram. That reading is wrong, and the caption now says why —
`n_plex` is counted over every quantified row, *before* `min_valid_per_group`,
and the retained set is `keep_valid & keep_plex`, so the two filters interact
and a protein quantified in every plex can still be one `min_valid_per_group`
drops. With the setting above 1 the report states the conclusion instead:
`analysis.min_plexes: %d keeps %d of the %d protein group(s)
min_valid_per_group passed; the other %d reach no differential abundance,
enrichment or shortlist table below, and at 1 - the default, where this filter
does nothing - every one of them would.` At the default, the advice that was
already there carries its price rather than naming a setting and stopping:
`Set analysis.min_plexes: 2 to drop them, leaving %d of the %d retained here.`

The counts are left as their format specifiers here deliberately. A triple for
the first real run was written into this entry and taken back out: it was
transcribed from a reconstruction rather than measured, and it contradicted
this same section, which establishes from that run's own
`min_features_per_protein` line that `annotated_quant.tsv` held 1,282 rows -
so no figure larger than 1,282 can be what `min_valid_per_group` passed on it.
The same rule was already being applied one paragraph away, where a plex count
that had not been measured here was refused; applying it to one transcription
and not the other is how a number nobody checked ends up quoted back as
evidence. Both directions are priced with the same pair of counts, and the
first
of them is the retained count the next line prints, so the sentence can be
checked against the rest of the page.

One comparison and not a ladder, for the same reason the bin escalation prints
a conclusion beside its table rather than a second table: what 1, 2, 3 ... each
leave is the quiet table this change exists to replace, and every point past
the retained set is non-linear anyway — `drop_zero_variance` runs after the
subset and `eBayes` moderates across whatever survives — so a ladder implies a
smoothness the code cannot support. The comparison is against 1 because 1 is
not a point on a ladder: it is the shipped default, the run the user would have
had without touching the key, and the only counterfactual this chunk can
evaluate exactly, since `keep_valid` and `n_plex` are both in hand and neither
depends on the setting. Escalated on `COVERAGE_MIN_N`, the floor every other
coverage claim in the document uses: a filter that removed at least as much as
it kept has redefined the experiment rather than trimmed it, which is a GATE,
but the same judgement over a handful of protein groups is an anecdote and gets
the quiet tier. Silent where the setting took nothing `min_valid_per_group`
would have kept, on a label-free run, and on a one-plex run at the default,
where every protein is confined by definition and the advice would be to delete
the experiment. Not silent on a one-plex run with the filter set: that run is
about to hand limma an empty matrix, and `keeps 0 of the N` is the most useful
line on the page.

**And the report no longer lets `analysis.min_plexes` be read as covering
`tmt.min_plexes`.** The two share a name, one is protein level and applied in
the report and the other is feature level and applied in the reader before the
roll-up, and the report can be exact about the first and cannot say anything at
all about the second: changing it re-runs the join. Where the feature-level
filter is actually above 1 — read out of `quant/design_notes.txt`, the same way
the reference treatment already is — the sentence above is followed by one
saying so. At its default nothing was filtered before the roll-up, the protein
counts are complete, and there is nothing to disclaim, so nothing is said.

**The join stage now says how many proteins it quantified from nothing.**
`{n}/{total} protein(s) had every feature dropped under '<mode>': they have a
row in peptide_evidence.tsv, no number in annotated_quant.tsv, and appear in no
report table.` Nothing anywhere printed this before — the features line counts
FEATURES, not proteins, and the `min_features_per_protein` WARN's own
denominator has already excluded them — so a reader filtering
`peptide_evidence.tsv` had no way to learn that a majority of its rows held no
number. It is the count removed from the dominance rate's denominator, owed to
the reader on every run that has any, INCLUDING the runs where that WARN is
silent. Deliberately INFO and never WARN however large it gets: in a
strain-redundant database this is the rule the user chose doing exactly what it
says, and a WARN that fires on every real run is the line a reader learns to
skip — taking the real one with it.

**The report reads `taxon_unique_dominated`, which until now nothing did.** The
flag reached `annotated_quant.tsv` and the R object's `rowData` and had no
consumer, so the only place the finding existed was one line of stderr that a
reader of the HTML report never sees. The feature-support chunk now states the
rate over the proteins THIS DOCUMENT is about, after both
`min_features_per_protein` and `analysis.min_features` — on the first real run
`784/1,282`, 61.2%, where the join stage's line says 54.7%. The two
denominators differ on purpose and each names its own: a stage reports the
population the stage decided, and the report cannot borrow that number because
it filters again afterwards. Gated on `COVERAGE_MIN_N`, the same floor every
other coverage claim in the document uses, and the join stage's
`ASSESSABLE_MIN_N` is pinned equal to it so the two halves of the tool call the
same size "too few" — below it the counts are printed without the percentage,
because a fraction over a handful of proteins is what gets pasted into a
methods section.

**Longest-first dispatch now discriminates inside the hours class, which is the
only place it ever mattered.** `stage_priority(name)` returned
`STAGE_COSTS[name]` and nothing else — a three-level ordinal, 3 = hours, 2 =
minutes, 1 = seconds. Eleven stages declare cost 3 and the durations measured
inside that one rank run from seconds to more than two days, so the ordering
that v0.4.0 added to stop InterProScan queueing behind a ten-minute stage left
it queueing behind its own rank instead: Python's sort is stable, the eleven
tied, and the whole hours class was dispatched in TABLE ORDER with `interpro`
seventh of it. On the 455,571-protein run InterProScan STARTED at **27.1 h**,
the figure TUTORIAL's own start table publishes, and it finished last — so
every one of those 27.1 hours was an hour the longest stage in the pipeline
spent waiting for a worker, and an hour on the wall. (The run's total and
InterProScan's own duration are quoted elsewhere in this entry from the
replay; they are not repeated here, because 57.9 and 84.2 do not subtract to
27.1 and the start time is the one this paragraph is about.) The same tie stands in front of
all eight configs in `examples/server-run-plan`, each of which annotates a
cohort in its own results directory and therefore starts cold, and each of
which lists interpro LAST in its own hours class rather than seventh.

**Each stage now carries `order_s`, the seconds it took on the release's
reference run, and the sort key is the pair.** `stage_priority()` returns
`(STAGE_COSTS[name], STAGE_ORDER_S[name])`; `run_now.sort(key=stage_priority,
reverse=True)` is unchanged, character for character. The rank leads, and that
is the containment rather than a convenience: no figure in the table — and no
figure a later refinement might read from somewhere else — can promote a
seconds-class stage past an hours-class one, so
`test_a_short_stage_never_outranks_a_long_one_wherever_the_table_puts_it` and
the orphan rule beside it hold for every possible duration and not just for the
numbers in the table today. What the seconds can do is break the tie, and that
is the whole of the defect.

**What it is worth, replayed rather than asserted.** The dispatch loop is
replayed over the real stage graph on this repository's own published
durations, once under the key v0.4.0 shipped and once under this one. The
replay is `test_the_measured_durations_replay_to_a_shorter_run` in
`tests/test_scheduler.py`, which asserts the COMPARISON and never a remembered
figure; every number in the table below is then recomputed from that same
replay and checked against this text by
`test_the_changelog_ordering_figures_are_a_replay_not_a_recollection` in
`tests/test_docs.py`, so neither this entry nor those tests can drift from the
other:

| durations | selection | `stage_workers` | cost rank alone | with the tie broken |
|---|---|---|---|---|
| 455,571 proteins | the whole pipeline | 3 | 86.6 h | 72.9 h |
| 455,571 proteins | the whole pipeline | 4 (the default) | 66.7 h | 59.5 h |
| 455,571 proteins | the eight cohort configs' selection | 3 | 63.5 h | 57.9 h |
| 38,204 proteins | the whole pipeline | 4 (the default) | 4.9 h | 4.5 h |

Three of the rows above land exactly on the dependency graph's critical path,
which is the shortest any ordering of the same work can be: there is nothing
left for a cleverer priority rule to collect, and what remains is
`stage_workers` and InterProScan itself. The remaining row — the whole
pipeline at `stage_workers: 3` — is bounded by the work rather than by the
graph.
**Note which floor is quoted.** `sum(durations) / stage_workers` is the
independent-jobs bound and it sits BELOW the critical path here, because
`integrate` is a barrier all of its feeders must clear and `integrate → esmfold
→ foldseek → finalise → taxonomy → join` runs out on one worker at the end;
quoting it would have promised hours that cannot be collected.

**The row that carries the result is `interpro`'s, and the eight configs in
`examples/server-run-plan` are why the table ships in the source rather than
being learned from a results directory.** Those configs run with `topology:
false` and `structure: false`, so their hours class is a different set —
emapper, pfam, ncbifam, kofam, interpro, with interpro LAST in the stage table
— and each is a first run in its own results directory. A self-calibrating key
would have saved them nothing, because a stage whose record matches is `cached`
and never reaches the sort at all; a table that ships with the release is right
on the first run, on a fresh clone, on a read-only directory and under
`--force` alike.

**Only the ORDER of those numbers is claimed, which is what makes one run's
measurements a fair table for every run.** List scheduling reads the order of
its priority list and never the magnitudes, so any monotone rescaling of
`STAGE_ORDER_S` produces a byte-identical schedule — checked, under several
monotone rescalings, in
`test_the_order_survives_any_rescaling_of_the_reference_seconds`. A machine
four times faster and a dataset an order of magnitude larger therefore schedule
the same way, and the repository's own evidence agrees: its two documented runs
are 11.9× apart in protein count and 8×–66× apart per stage, and they order
every stage both of them timed identically. So the numbers are not normalised,
averaged or fitted, the comment at `STAGE_ORDER_S` says which rows are
measurements and which are placements for the stages that run never enabled,
and it says in full sentences why nobody should "fix" them into a constant per
rank: a constant is not a rescaling, it is the ordinal again.

**The placed rows are the one part of this that argues with itself, so their
sensitivity is measured too.** jackhmmer, hhblits, unipept, context, smorf and
taxonomy have never been enabled on a run this repository publishes durations
for, so they carry a placement below every measured stage of their own rank on
a stated argument — the hours-class ones query the dark set, a few percent of a
proteome, rather than all of it.
`test_the_saving_does_not_rest_on_the_stages_no_run_has_timed` scales those
placed rows in both directions and reports what happens: down, nothing changes
direction; up by ten, so that stages nothing has ever timed become the longest
things in their rank while still being dispatched near the bottom of it, the
new key loses a little on the smaller dataset and nothing at all on the larger
one. That is the real cost of a placement being wrong about ORDER, it is
bounded by the length of the one misplaced stage, and it is bounded at all only
because the rank leads.

**The cost rank survives untouched, and #30's division of the machine is
bit-for-bit what it was.** `stage_cost()` is split out for the two callers that
want a CLASS rather than an order — `share()`, which weights the CPU and RAM
cut, and the WARN about a stage that cannot grow into freed CPU — because a
rank is a RATIO of a machine where the seconds are only an ORDER of dispatch.
Weighting the box by durations would have given InterProScan 21 of 22 cores and
left each of its round-mates one, for the life of the run, since a tool's
thread count is fixed when it is launched. Within a rank the weights are equal,
so reordering two same-rank stages moves only which of them takes the remainder
of a division of equals; mixed rounds are unchanged. Handing either caller the
pair raises instead of dividing by a duration, and
`test_the_dispatch_key_is_not_a_quantity_the_machine_can_be_divided_by` pins
that. Neither number reaches a cache signature or any `keys` list — scheduling
order cannot change a stage's output — and both accessors still index rather
than default, so a stage added without a cost or without a duration raises
instead of being ranked as trivial or sorted last in silence.

**The two tests that were meant to prove the ordering worked had pinned the
half of it that did not.** `test_equal_cost_stages_keep_table_order` asserted
the tie order as though table order were a neutral tiebreak — it is not a
neutral one, it is a specifically bad one — and
`test_the_first_wave_no_longer_goes_to_the_shortest_stages_in_table_order`
asserted a first wave of `emapper, pfam, signalp, tmbed`: a first wave with the
longest stage in the pipeline ABSENT, pinned as correct. The mechanism was
tested and the discrimination never was, which is why a green suite could guard
a feature that did nothing. Both are rewritten to say what they now pin, the
first as `test_interpro_is_dispatched_before_every_measurably_shorter_stage`,
which fails against the ordinal-only key rather than passing beside it. Worth
recording: the old tie test had already become vacuous the moment the key
became a pair, because `stage_priority(n) == 3` is false for every stage and it
compared two empty lists. Stability is still required and is still pinned,
narrowed to the claim it can honestly make — two stages the reference run
cannot tell apart keep table order.

**The order is disclosed once per run, on the log, and published per stage in
`describe --json`.** Until now the within-rank order could be predicted by
reading `STAGES` top to bottom; it now comes from the seconds beside each
`cost`, so the run names the order it is taking before its first dispatch —
once, not per stage, because `stage_priority()` is a total order over stage
names and the order of any round's ready set is that order restricted to it. No
durations in that line: they are this release's figures for another dataset and
printed beside a stage about to start they would read as an estimate of the run
in hand, which is exactly what the resource guide tells an operator not to do.
`describe --json` grows `order_s` beside `cost` instead, so a front end can
reproduce the dispatch order from the release alone, and `--dry-run` prints
the order it would really take — which it can do exactly, and could not do at
all if the order depended on a results directory it is forbidden to write to
and may never have read. `DESCRIBE_VERSION` does
NOT move: the rule at the constant is that it moves when a key is removed or
its meaning changes, and `cost` still means the rank. What changed is that
`cost` is no longer SUFFICIENT to reproduce the order, and the honest answer to
that is the new key beside it rather than a bump no consumer gates on.

**The console had to move with it, and it is a companion change rather than a
consequence left lying.** `console.rank_ready()` ranked a NEXT row by `cost`
alone, which annotated the hours class in the console's own table order over a
queue the engine takes by duration — the page would have said the engine ranked
six of the ready stages ahead of InterProScan when it ranked none ahead of it,
inside exactly the eleven stages an operator watches, and `next_detail()`'s own
docstring records what that class of drift cost the last time (the page listed
dbcan NEXT above ncbifam NEXT, the engine dispatched ncbifam, and dbcan waited
27 hours). That is why the second key is in `describe --json` rather than only
in the engine: the console reads both fields and sorts by the pair, and
`ready_rows()` declines to rank at all when either is missing rather than
inventing one. No new xfail was taken anywhere for this — an engine-created
console defect recorded as a permanent expected failure would have been a worse
trade than fixing the page — and the contract test asserts the whole order
again, stage for stage, with the between-rank half kept beside it because that
half must hold even when a duration is wrong.

**The prose that advertised the feature had the defect written into it, and is
corrected.** README's scheduler section said of the ordering: "The order is a
starting order, not a schedule: it makes nothing faster, and InterProScan still
paces a large run." The first clause is true, and the second was the symptom
stated as if it were the design — a sentence written by people who had no
reason to believe the ordering mattered, in the paragraph that sells it.
TUTORIAL's sizing section had the same shape from the other end: it credited
v0.4.0 with the half that worked, that dbcan waits now, and said nothing about
the half that did not. Both now say which half v0.4.0 fixed and which half this
change fixes; README's Scale table says that its own figures ARE the order
table and which stages carry a placement instead; and CLAUDE.md's "things that
are deliberate, do not fix them" list gains a bullet on what not to "fix" here
— not the magnitudes, not the accessor split, and not into a reader of the
state file. The pinned suite counts in README's Tests section are re-derived
from a real run rather than adjusted by hand.

**Weighed and turned down.** Reading each stage's own `seconds` out of
`.metaannot_state.json` and sorting by that: it would make scheduling a
consumer of a document v0.6.0 and #39 have just finished making honest about
being missing, unparseable or written by another run, it pays nothing on
exactly the runs that need it most (a first run, a fresh cohort directory, the
GPU box whose other stages arrive as `adopt` and carry no duration at all), and
`seconds` is not a property of a stage but of the stage, the dataset, the box
and the CPU cut `share()` happened to give it — InterProScan's own figure was
measured at `-cpu 7` because signalp and tmbed were beside it — so feeding it
back closes a loop in which the scheduler's past decision is its present input.
It remains a defensible refinement LATER, and it needs this table anyway as the
per-stage fallback, which is strictly more information than the class median or
the "unknown means longest" sentinel a bare state-file key has to invent. Also
turned down: longest-remaining-path-first, which on this graph is the same
order — `integrate` is the single confluence, so every feeder's remaining path
is its own duration plus one shared constant, and
`test_remaining_path_length_would_order_this_graph_the_same_way` records that
precondition so the day it stops holding the question is asked again rather
than assumed away. And a WARN where a stage's measured duration contradicts its
declared rank: a good idea, and a separate change, because one of the
hours-class stages would fire it immediately and legitimately and a new WARN on
every run does not belong in the same commit as a scheduling change.

**A `cached` stage now says when the record it is being reused on is no longer
about the file that is there.** On the `cached` branch, and for an `ok` record
only, `decide()` compares each declared output's mtime against the `finished`
stamp of the record that describes it. Every output written more than
`OUTPUT_STAMP_SLACK_S` after its own record goes into ONE WARN for the whole
run, naming each stage, each file, when it was last written and the run whose
record it is — a report whose length grows with the number of late stages and
whose line count does not. The verdict does not move: every one of them is
still reused, no record is rewritten and nothing is deleted.

**The gap it reports is the one the concurrency fix left open and could not
close.** A run `SIGKILL`ed between `atomic_out`'s rename and `finish()`'s record
is *not* it — `mark_running` wrote a `"running"` record before the stage
started, and the next run recomputes on that. What is left is the case where the
record that survives belongs to a different run than the bytes do: A marks a
stage running, the operator `--force-unlock`s, B recomputes that stage and
records it `ok`, A's rename lands on top of B's file inside
`min(heartbeat_s, STATE_PROBE_S)`, and A is killed without writing anything at
all. `signature()` hashes inputs and config and never an output, so nothing else
in the program is capable of noticing. The same shape needs no kill whatever: a
table put there by hand, a parked `.superseded.*` moved back, an over-broad copy
back from the GPU box.

**Nothing in the state file moved for it.** The exact key-set assertion on an
`ok` record in `tests/test_console_contract.py` is unchanged, deliberately: the
check is built out of `finished`, which every `ok` record has always carried, so
a results directory written by an earlier build needs no migration and a record
this build cannot date is compared against nothing. The console is untouched
too — it already stats every declared output and shows both the file and the
`finished` stamp on the row, and the engine, which knows the signature it has
just matched, is what says the sentence.

**It is evidence and not protection**, and the limits are in the message itself
rather than only in the README. An output that is *not* newer proves nothing at
all: an overwrite that landed before the record, or one that kept the file's
timestamp, leaves no trace here. Only DECLARED outputs are compared, so for
`diamond`, `hhblits` and `esmfold`, whose declared output is a `.done` sentinel,
a sentinel that dates cleanly says nothing about the per-database tables, the
per-query `.hhr` files, `plddt.tsv`, `esmfold_failed.tsv` or the per-protein
PDBs beside it. (Those tables and those PDBs go through `atomic_out` and *are*
covered by the ownership gate; they are outside this check because they are not
declared.) And `emapper`'s live branch writes its annotations file with no temp
and no rename, which makes it the one declared output a kill can genuinely
truncate — a timestamp says nothing about truncation, and this does not claim
to.

**Weighed and turned down: recording each output's SIZE in its record instead.**
A size survives a copy where an mtime does not, but it is blind in the very
window this exists for — A and B ran the same stage over the same inputs, so
byte-identical output is the expected case there, not a coincidence — and it
would have moved a pinned record shape, put a key no reader acts on into a
published document, and left a permanent warning on a legitimately re-`rsync`ed
adopted output with no supported way to silence it. The size-shaped hole is
named in the README instead: a substitution that preserves a timestamp is
invisible here.

**One case does report itself on every run, and nothing silences it.** A stage
this box ran and recorded `ok`, over which the operator later `rsync`s a newer
output from the GPU box, is a record that really is no longer about the file
that is there. It is reported, correctly, and nothing re-records it — and the
remedy the message names, `--force --only <stage>`, would recompute a GPU stage
on the wrong machine. Re-recording a stage without recomputing it would answer
it and is deliberately not in this change.

**Where several stages are late at once the cause is not named, because it
cannot be had from a timestamp.** A directory copied with `cp -r`, unpacked from
an archive or restored from a backup has every mtime in it rewritten at one
moment, and so do that many separate replacements: measured on this project's
own results directory, a `cp -r` puts every late output inside twenty
milliseconds of every other, and so does overwriting the two outputs a
`--only pfam dbcan` run leaves. The report offers both readings, asserts
neither, and says that rebuilding stage by stage would spend hours of compute
for nothing if it is the first.

**The false warnings are pinned as hard as the true one.** `--force --only`
compares nothing, having already decided to recompute. An `adopted` record is
never dated against its file, because `finished` there is when this box noticed
it and `rsync -a` preserves the source mtime — dating it would report the
documented GPU hand-off on every run for ever. A `finished` stamp that names no
single instant is declined rather than dated, and the run says which record. An
ordinary resume of a whole directory reports nothing. And there is one
suppression, because it is the only one that is a MEASUREMENT: a filesystem
whose clock leads this machine's, taken on a file the run has just written,
makes every output look newer than its record, and two clocks compared against
each other say nothing at all — so the run says that once and compares nothing
further. A `--dry-run` writes nothing, so it cannot take that measurement, and
it says the question went unasked rather than leaving a reader to infer it was
asked and answered no.

**A run now names every record that never reached the state file, and says
what the next run will do with it.** A declined state write used to be silent
to everything but the return value nobody reads: the record was gone, the
outputs were on disk, and the operator found out at the start of the next run.
Each declined write is now ledgered by key and named twice — once at the moment
of the loss, and once in an account `stamp_run()` gives before the run exits,
which names the stage records that never landed, quotes the error verbatim, and
says which of the two residues the next run will act on: a stage this run had
already recorded `running` is RECOMPUTED, a stage with no record at all is
ADOPTED with the warning that says outright it cannot tell a finished file from
an interrupted one. The ledger is a ledger and never a work queue — nothing
reads it to decide what to write, which is the whole difference between it and
the deferral this replaced. A key name carried forward is a claim about this
run's own I/O and stays true; a RECORD carried forward would be written under a
later answer about who owns the directory than the one it was decided under,
which is the shape of the defect the missing/unparseable arm of
`update_state()` was written against.

**The account is said twice because one of the two exits cannot be covered.**
`stamp_run()` is the funnel every orderly ending goes through — `cmd_run`'s
"ok" and "failed", `main()`'s "failed" for a `die()`, `main()`'s "interrupted"
for Ctrl-C or SIGTERM — and the final stamp is itself a merged write, so it can
be the last thing to fail and has to be able to appear in its own account. A
`kill` reaches `_release_lock_on_signal()`, which may not format a string, take
a lock or touch buffered I/O, and `TerminateProcess` on Windows runs nothing at
all: there is no stamp and no account on those paths, which is why the identity
of each loss is also said when it happens rather than only at the end.

### Fixed

twenty-four entries, in six groups. Each heading carries its own count and a
test counts the entries under it.

#### Four defects in the check this change set added

**The report is one WARN for the run, because suppressing a per-stage one was
all-or-nothing and the remedy it recommends destroyed it.** The first build
logged a line per late stage and suppressed the lot whenever EVERY `ok` record
in the document was late, reading that as a copy rather than a replaced file.
Driven with ordinary commands: run a project to completion, wait past the
slack, `cp -r` it, and the copy says one line; then `--force --only pfam` — the
action every one of those lines names — re-records one stage, "every record is
late" stops holding, and every run after that says the per-stage sentence about
each of the stages left, for ever, with nothing that silences them. That is
exactly the failure `OUTPUT_STAMP_SLACK_S`'s own comment claims to prevent,
reached by following the advice. The count is now the property: one line before
the `--force --only` and one line after it, and
`tests/test_scheduler.py` drives all four steps rather than backdating every
record, which is what the test that stood there did and is why it could not see
this.

**Two `ok` records is two records.** `--only pfam dbcan` leaves the minimum any
`--only` run leaves, and overwriting both outputs is the shape this check
exists for — a superseded run renaming inside `min(heartbeat_s, STATE_PROBE_S)`.
The suppression answered it with the directory-wide sentence about a directory
"copied, extracted or restored": neither stage named and a cause asserted that
nothing could know. Both are named now. The cause is not claimed at all,
because three candidate ways to claim it were driven over real directories and
none of them works — the mtimes of the late outputs cluster just as tightly for
two overwrites as for a `cp -r`; the state file's own mtime is rewritten by
every run, including the `--force --only` above; and a broken dependency order
is what a single replaced input produces by definition, so suppressing on it
swallows the case the check is for. The three are written out at the check.

**A `--dry-run` cannot measure the filesystem's clock, and now says so.**
`_fs_clock_ahead()` is taken on a file the run has just written, a dry run
writes nothing on purpose, and the skew branch was therefore unreachable there:
with a throwaway engine reporting a lead of a thousand seconds over a directory
holding one late output, a real run read no mtime as evidence and a `--dry-run`
reported the stage. Neither the message nor the README said the question had
gone unasked; the README called a dry run an audit and left it there. The dry
run's report now names the gap, and the README says it plainly.

**A stage recorded in the hour a zone repeats is declined rather than reported
for ever.** `_stamp_epoch()` is `time.mktime(time.strptime(...))` on a naive
local time and cannot resolve a repeated local hour: with `TZ=America/New_York`
the second 01:30 of 2026-11-01 is epoch 1793514600, `strftime` writes
`2026-11-01T01:30:00`, and reading that back gives 1793511000 — an hour early,
twelve hundred times the slack. `_stamp_instant()` asks PEP 495's `fold`
instead, which answers for the skipped hour and for half-hour shifts too, and
returns nothing where the text names two moments or none. Round-tripping the
text through `localtime` was tried first and does not work: the earlier reading
of a repeated hour formats back to the same text, so the round trip succeeds on
exactly the case that needs catching. The record is named rather than passed
over in silence, because a stage that went unjudged must not look like one that
dated cleanly. `_stamp_epoch()` is left as it was, since it is what the adoption
paths have always read.

#### Seven defects in the state file's I/O error handling

**A state write that cannot be made no longer kills the run, and no longer
reports the stage that SUCCEEDED as the one that failed.** `update_state()` had
exactly one `try` in it and it was around the read-back; `payload =
save_state(path, merged)` was bare. A read needs no blocks, so ENOSPC — the
errno this is about — could never reach the `unreadable` arm that declines
politely: it arrives at `open(tmp,"w")`, at `fh.write`, at the implicit
`close`, or at `os.replace`. Driven without mocks, twice: a real ram-disk
filesystem with no free blocks, and a real read-only results directory. The
OSError left `update_state`, left `finish()`, and from the dispatch loop took
the run — at hour thirty, for an InterProScan whose output was already complete
on disk. From the DRAIN loop it was worse than a traceback: `finish()` sits
inside an `except BaseException` there, so the bookkeeping error of a stage
that had succeeded was logged as `stage 'X' failed` and the run exited 1 — the
standing rule that a succeeded stage is never failed by bookkeeping, inverted
in shipped code. Both halves of the I/O decline now, and the guard is
`except OSError` and never `except Exception`: a `TypeError` out of
`json.dumps` is a record this program should not have built and still comes out
loud at the call site.

**A log file that cannot take a line costs the line and not the run.**
`_write_safely()` caught only `UnicodeEncodeError` and the `_LOGFH.flush()`
beside it was bare, while `p.state` and `p.logfile` are both under the results
root and share a filesystem. So the ENOSPC that stopped a record being written
also stopped the line explaining it — and since `update_state()` says what it
declined through `log()`, the state guard above would have been decoration
without this one. Driven on the same full filesystem, and the measurement
decided the shape: the buffered `write()` returned normally and the `flush()`
that followed it raised, so a guard on the write alone — the obvious one to
write — would have caught nothing. Lines are still attempted afterwards rather
than the handle being dropped, so a log resumes by itself once space is freed,
and the notice is said once. stderr is deliberately left raising: a log file
nobody is watching live and the channel `tmux` keeps are different facts.

**The pre-write read is opened twice before it is called unreadable.** One
`STATE_READ_TRIES`-bounded retry, with no sleep, inside
`_read_state_for_merge()`. An NFS `ESTALE` is a handle the client revalidates
on the next `open()`, and a Windows sharing violation is an indexer or a backup
agent holding the file for a few milliseconds; costing a record for either was
the cheapest thing here to stop. It has its own budget rather than a share of
`STATE_WRITE_TRIES`, whose comment defines it as the anti-clobber budget, and
it sits in the reader so that the rename path's ownership probe gets the second
open too. No sleep, because there is nothing to wait for that a second
`open()` does not do and because a sleep there would run inside `_STATELOCK`;
an error that genuinely needs time is declined, ledgered and reported instead.

**The heartbeat's give-up warning is armed again, and only by I/O.** `_beat()`
counted exceptions, `_tick()` returned `not self.superseded`, and a `_save()`
that returned False therefore reset the counter to zero — so on a read-side
`EACCES` or `EIO` the warning had never been able to fire, and guarding the
write half would have made that true of every I/O failure there is.
`_run.last_seen` would then freeze in the file for the rest of a three-day run
with nothing said, which is the one conclusion a live run must not invite. The
counter reads the ledger now, so only an I/O decline arms it: a superseded run,
a lock somebody else holds, a VACANT lock and a tick that cannot get the
scheduler's lock all still decline without putting a disk warning in front of
an operator who caused none.

**A state write that fails leaves no temp behind.** `save_state()` unlinks the
`.part` it minted microseconds earlier, on the failure path, inside the lock it
already holds, and never touches the state file itself. The reason is the
BLOCKS and not the clutter: on a full disk a half-written temp holds the space
the next attempt needs. The clutter reading was measured and is wrong — the
name carries the pid and the thread, so five failed writes from one thread left
ONE file — and the docstring's `find results -name '.*.part.*'` sweep still
says what it always said, because it is about a writer KILLED mid-update and no
`except` clause runs for that.

**A `running` record is no longer reported as an interrupted writer without
qualification.** `mark_running()` writes `running` before the stage starts and
`finish()` writes the real record hours later, so a state file that goes bad in
between leaves exactly that record for a stage that finished perfectly well —
and `decide()` said flatly that "the previous run was interrupted while this
stage was writing". From the next run's position those are the same bytes, and
the sentence now names both readings. The VERDICT does not move: what is on
disk is as likely to be half a file as a whole one, and weakening that record
is the obvious next idea and the wrong one.

**The line about an unreadable state file no longer states a policy that the
next write falsifies.** It said "nothing is written to it until a read
succeeds", which stopped being true at the next successful write, and it was
said on the FIRST failed open — where `_said_once()` consumes its key for the
life of the run, so a single blip spent the one line a filesystem going
permanently bad an hour later would need. It is said after the read budget is
out, and it describes the write it declined.

#### Two places that measured the wrong thing

**`_unreleased_fixed_groups()` read a shipped release.** Cutting v0.6.0 left an
`## Unreleased` holding only `### Added`, and `txt.index("### Fixed", start)`
walked past that section and landed under `## v0.6.0` — every group of a shipped
release read out as though it were in flight. The three tests built on it went
on passing, because a shipped section agrees with its own prose for ever and
will keep doing so, which is the worst way for a scan to be wrong: they had
stopped checking the section in flight and said nothing about it. `_newest_section()`
and `_fixed_groups()` take the bound in their arguments now, a test pins the
boundary on a document built to have the trap in it, and each consumer asserts
the PAIRING — groups without the sentence that counts them, or that sentence
without the groups, both fail — so a section with nothing to correct asserts
that it counts nothing rather than asserting nothing at all.

**The README's test counts were a release behind.** The Tests section quoted a
passing count from before this change set added its tests. The band in `_near()`
is wide enough that nothing failed, which is the point of the band and not a
reason to leave a counted number wrong: the README's own claim about those
numbers is that they are counted rather than estimated. Re-counted, with the
no-R pair re-derived from them.

#### Six defects in what a stop does to the tools a run launched

**A `kill` of a run now stops the tools it launched, and it stops them before
it lets go of the lock.** The signal handler exits through `os._exit()`, which
runs no `finally` and no `except BaseException`, so `run_cmd`'s own cleanup was
never reached and *every* child of a killed run survived — documented
behaviour with an uncosted price. On the first full run on real data that price
was three concurrent `hmmsearch --cpu 7` against the same 455,571 proteins from
13:03 to 20:04, the box 1.6x oversubscribed during `pfam`, and two of the three
finished results discarded because a dot-prefixed `.part` file is something
no stage can adopt.

Each tool is now launched through one helper with `start_new_session=True`, so
it is a session and process-group leader (`pgid == pid`, no extra syscall) and
one `os.killpg` takes it and everything it spawned. That is also the only thing
that ever reached a launcher's children: `proc.kill()` left them alive with
`ppid 1` and still inside metaannot's OWN process group, where no group kill
could reach them without suicide — measured, and it is the shape of
`interproscan.sh` -> java, `emapper` -> its children, and torch -> its
dataloader workers holding VRAM. The kill happens on every stop path: the
raw signal handler (a bare `SIGKILL` to each group, since a handler may not
wait), an unwinding `run_cmd`, the dispatch loop, and `main()`'s
`KeyboardInterrupt`; the three that unwind send `SIGTERM`, wait a few seconds
so a tool can remove its own scratch tree, and then `SIGKILL`. The success path
sweeps the group too, after asking with signal 0 — a tool can exit 0 leaving a
live member behind, and "the tool said it was done" is not the same statement
as "nothing of the tool is left".

**Ctrl-C had to be fixed in the same change or it would have got worse.** What
stops a tool on a keyboard Ctrl-C today is not any line in this program — it is
the tty broadcasting `SIGINT` to the whole foreground process group, and
`run_cmd`'s `except BaseException` is unreachable during dispatch because it
runs in a worker thread while Python delivers signals to the main thread only.
A session of its own removes that broadcast, so the interrupt is now caught
*inside* the stage pool's `with` (not around it: `__exit__` joins, so a kill
placed after it waits for the very stage it is stopping) and a latch stops the
workers from starting the next tool while the run unwinds. Two pre-existing
holes close with it: a `kill -INT` from a script, with no terminal, never
reached a tool at all and left the main thread in a join measured in hours; and
a Ctrl-C mid-`tmbed` failed one chunk and then ran every remaining chunk to
completion, because that stage catches `RuntimeError` per chunk — which is why
the refusal is a `StageError`, the one exception it re-raises.

**And a leftover census at the head of every run, because the measured incident
ran no handler at all.** The 12:58 log line `removing a stale lock from pid 336`
is written only on the branch that finds a lock FILE still on disk, and the
handler removes that file — so pid 336 never reached its handler, and nothing
in a handler could have prevented those 98 core-hours. `SIGKILL`, the OOM
reaper and a host reset all leave the same state, and the group kill makes it
marginally worse, since a tool in its own session no longer dies with a closing
tmux pane either.

So after the lock is taken and before anything is dispatched, every run scans
its declared output directories for names matching `atomic_out`'s convention
(derived from `ATOMIC_SUFFIX`, so it cannot drift: a leading dot, a stem, two
all-digit fields, `.part`), aggregates them by directory and pid — mandatory,
not cosmetic, since an interrupted `esmfold` legitimately leaves one per dark
protein — says one WARN naming the sizes, the newest mtime and which stage
writes those names, and then stats them again two seconds later. A file
whose **size** changes is positive proof of a second writer, and the run is
refused with the `lsof`/`fuser` that finds it; a file that does not change is a corpse
and gets one line. Size is the primary evidence because it is an integer,
monotone for an append-only tblout, needs no clock, and is immune to the
two-second SMB and FAT mtime granularity `OUTPUT_STAMP_SLACK_S` already
documents; the mtime is compared against `.metaannot_state.json`'s own mtime on
the same filesystem, stat'ed before this run rewrites it, and suppressed above
the skew threshold `report_output_stamps` already uses. The two seconds are
spent only when a candidate exists.

**Nothing is deleted, renamed, moved, truncated or read past `stat`, and
nothing is ever signalled on the strength of a filename.** The pid in
`.pfam.336.140234.part.tblout` is metaannot's — `atomic_out` mints it with
`os.getpid()` — so it names the run that opened the file and says nothing about
the `hmmsearch` that was writing it; it carries no host, and it can have been
recycled. The census therefore asks the process table exactly one question, and
it is `== os.getpid()`. The 237 MB tblout from the incident is still in
`results/hmm/` and stays there: rule 5, and because it is the only artefact
showing three hmmsearches ran, because an `unlink` would not free the space while
the writer holds the inode — it would only hide the bytes from `du` and destroy
`lsof <path>`, the one handle from the file back to the process — and because
a complete tblout may be real work.

**What this does not cover, stated so nobody reads silence as health.** No
handler runs for `SIGKILL`, an OOM kill or a host reset, and nothing in a
process that has been killed uncooperatively can act; the census is the whole
of the answer there. Under `systemd` the question does not arise —
`KillMode=control-group` is the default, so a supervised run's cgroup is torn
down whatever we do — and the leak was always specific to a bare `kill` from a
shell or a tmux pane. A tool that calls `setsid` itself escapes the group kill;
none in this tool set does, and that cannot be proved for a tool it has not
met. An orphan that leaves no pid-bearing name is invisible to the census, and
the list is the one `_park_superseded` already gives: `hhblits`' per-query
`.hhr` and its own `.hhr.part`, `esmfold`'s `plddt.tsv` and
`esmfold_failed.tsv` appended in place, the `.done` sentinels. Windows is
uncovered: `TerminateProcess` runs no handler and there are no POSIX process
groups, so no slot is claimed and the kill loop cannot run. On Windows the real
answer is a Job Object with `JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE`, which is not
testable from here and is not attempted.

**Two costs, both real.** `emapper`'s live branch hands `emapper.py` the output
path and lets it write `eggnog/emapper.emapper.annotations` in place, with no
temp and no rename — the one declared output a kill can leave genuinely short.
A group `SIGKILL` now truncates it in cases where an orphan would have finished
it, and what stands between that and adoption is only the `running` record
`mark_running` wrote first. That is the documented, correct outcome
(recompute), but it is reached more often; routing `emapper` through
`atomic_out` is the real fix and is a separate change. And Ctrl-Z now suspends
metaannot while its tools keep running, because a tool in its own session is
not in the foreground group — forwarding `SIGSTOP`/`SIGCONT` would fix it and
is a third and fourth signal in the function whose history includes a measured
hang, so it is deliberately left out and written down instead.

#### Four defects an audit of the two changes above found

**The recycle window the one dangerous syscall rests on was a CLOCK, and the
clock was five seconds.** The comment above the handler's kill loop stated the
residual hole "as a count and not as 'small'": between the wait that reaps and
the store that clears the registry slot "there are a few bytecodes", so for a
signal landing there to reach a stranger the pid space would have to wrap —
99998 creations on macOS — inside them. That was false, and false by the whole
of a `reader.join(timeout=5)`: the reap is `run_cmd`'s own `proc.wait()`, and
the slot was cleared by the spawn helper's `finally`, which is reached only
after that join. Worse, it was the entire timeout in exactly the case the
success-path sweep exists for — a surviving child inherits the tool's stderr,
so `iter(stream.readline, "")` never sees EOF. Measured, with a tool that exits
0 leaving one child: reap to sweep 5.02s, reap to slot cleared 5.07s, the
leader already out of the process table and the slot still advertising its pid,
at which `_end_tool_group` aimed `killpg(pg, 0)` and then `SIGTERM` and
`SIGKILL`. The sweep and the store now both happen AT the reap — the helper
hands the caller that does its own reap an idempotent `finish_group()`, and
`run_cmd` calls it between the wait and the join — so a group already empty
there is forgotten a few bytecodes later having been sent nothing at all, and a
group with a survivor in it still holds its own pgid, which makes the number
unrecyclable until the sweep's poll notices it go. The residual window is one
of those polls, 0.05 s, and the comment now says clock where it said count.
A consequence rather than the reason: the join returns instead of expiring,
because the survivor holding the pipe open has just been killed.

**The census made `--force-unlock-live` inert, and it now tells the two writers
apart instead.** That flag's help text is "take the directory even when this
host can see the holder is still running", and a live holder's tool is
precisely what leaves a `.part` file growing — so the growth check refused the
flag's only case, every time. Measured: run A live inside `pfam` with a growing
`.part`, then B with `--only pfam --force-unlock-live` — B exited 1, refused by
the census, and dispatched nothing, with no escape anywhere. It took the
documented handover with it, since `another run holds this results directory`
and the whole `_park_superseded` apparatus the TUTORIAL describes then named a
state nothing could reach. The suite missed it because the handover tests use
stubs that have already exited.

What tells the cases apart is the only thing a leftover's NAME carries: the pid
that MINTED it, which is metaannot's own and therefore names a RUN. Where that
is the pid of the lock this run has just taken from a holder **this host proved
alive**, the file is the displaced holder's tool and the run proceeds behind a
loud `WARN` naming the pid, the growth and the fact that there are now two
writers in the directory on purpose. Anything else still refuses — including a
file lying beside it, because the exemption is per FILE and the flag is an
assertion about one process rather than about the directory. It is keyed on
what the lock really DISPLACED and not on the flag being typed, so
`--force-unlock-live` over a vacant, garbled, remote or provably dead holder
excuses nothing and the `kill -9` orphan is refused exactly as before. Nothing
here asks the process table anything: `ResultsLock` did the proving before the
lock changed hands and the census only consumes its answer, which keeps the
rule that a pid out of a filename is never read as evidence about a process.
The one case it cannot separate is a file minted by a dead run whose pid was
later recycled onto the live holder; the message prints the number and points
at `lsof` rather than calling it proof.

**`CLAUDE.md` is now one of the surfaces the count scanner reads.** It was
outside it entirely, and a standing-instruction file is the worst place for an
unpinned count: its whole job is to stop the next reader reintroducing a
defect, so a number that is wrong there is a rule that is wrong where it is
most read. The rule this change adds to it names the spawn-helper allowlist,
and that number is DERIVED from the allowlist itself — counted as its
`subprocess.run(` entries, so it cannot be "length minus one". Two counts were
caught by the scanner in the drafting of this very change, which is the
argument for the surface; no released `CLAUDE.md` ever carried a wrong one,
and an earlier draft of this entry said otherwise. Making the file scannable
cost one further classification — "three parts", a trio the bullets under it
spell out.

**The test named for `_run_rscript`'s spawn site did not pin that spawn site.**
It was cited as covering the Mac's whole exposure to the group kill, on the
grounds that its spy had moved from `subprocess.run` to `subprocess.Popen`.
Measured: with `metaannot.py` reverted to the commit before this change set and
that test file unchanged, it PASSES — `subprocess.run` calls `Popen` itself, so
a spy there sees the launch either way, and moving it changed what the test
watches and not what it proves. Its name is about the resolved PATH, which it
does pin, so it keeps that job and says so. What actually changed at that call
site is the session and the registration, and a second test now asserts those:
`start_new_session=True` in the launch, a tool slot claimed for the launched
pid — only the helper claims one — and `stdin` on `/dev/null`. It fails on the
old spawn site, which is what the one beside it never did.

#### One defect in the join stage's most quoted number

**The taxon-unique dominance rate was denominated on rows that could never be
counted in it, and reported a majority finding as a minor caveat.** `ev` is
built by an OUTER join of the per-protein counts of KEPT features with the
per-protein counts of DROPPED ones, so a razor protein whose every feature was
dropped arrives as a row whose `n_features_used`, `n_unique`,
`n_taxon_unique` and `n_family_unique` are all 0. The flag is
`(n_taxon_unique + n_family_unique) > n_unique`, and `(0 + 0) > 0` is False,
so such a row can never reach the numerator while still swelling the
denominator — and those rows hold no measurement at all: they have a row in
`peptide_evidence.tsv`, which is where this change's new INFO line counts them,
and no number in the quant matrix and no place in any report table. On the first full run on real data 5,039 of the
8,238 rows were such rows, and the line read **1,749/8,238 — 21.2%**. Over the
proteins the assignment rule actually decided something about it
reads **1,749/3,199 — 54.7%**, which is the same numerator and the same
`peptide_evidence.tsv`, byte for byte: any archived run can be recomputed both
ways from the file it already wrote. The denominator is not a new invention —
`before` in the `min_features_per_protein` WARN two lines below is `len(quant)`,
i.e. exactly this count, so that line already reported `1,282/3,199` beside the
old line's 8,238 with nothing saying they were different populations.

The per-protein flag is unchanged, and that is the point: it was never wrong.
`taxon_unique_dominated` is true of exactly the rows whose taxon- plus
family-unique features outnumber their own unique ones, a tie is False because
"rests more on" is strict, and under `protein_unique` it is structurally
always False rather than broken. No per-protein claim, and no published
column, moves. What moves is one summary line, in one direction: the loudest
number in the join log more than doubles on that run, with no change to the
classifier, so two logs from the same data are no longer comparable line for
line.

Not fixed here, and not folded in: under `peptide_assignment: razor` nothing is
dropped, so a protein quantified entirely from cross-taxon `shared` features
sits in the denominator with the flag False — worse than dominated and counted
as clean. The flag is not misclassifying it; its coverage is incomplete, which
wants either a second flag or a decision to widen a predicate whose current
meaning is exactly right. Under `razor` the rate is therefore a floor.


## v0.6.0 — 2026-09-12

Everything here came out of one question — what can a program, rather than a
person, find out about a run before and while it happens — and out of answering
it honestly enough that the answers could be checked.

`doctor --json` is the headline: the other half of the contract `describe
--json` started, and the piece the console's M2 preflight was blocked on.
Holding it to its own stated scope turned up **twenty-four live false verdicts**
in the shipped `doctor`, four states of a path on which it printed no document
at all, and one on which it never returned.

The defects fixed alongside it were each found by driving rather than by
review. A **FIFO** at any operator-supplied path hung a run for ever with no
output and no exit status, and the obvious fix was wrong twice over before the
right one appeared — the rule in the end was that a pipe may be read once per
run, derived from the code rather than assumed, because `mkfifo p; zcat
big.faa.gz > p &` never worked at a path a run opens more than once. **Two
concurrent runs** could erase each other's stage records, and that took as long
as it did because each guard in turn proved unreachable on the very path it was
written for. And the **report** now says out loud when a bin contributes
nothing to the model, instead of printing two zeros in a tibble and leaving the
reader to work out that the statistics cover none of the fraction this tool
exists to study.

`SIGNATURE_VERSION` does not move, and that was checked rather than assumed: no
stage's `keys` list changed, so upgrading recomputes nothing.

`doctor --json`, the other half of the contract `describe --json` started. The
console's milestone M2 is a preflight checklist — can this config actually run,
here, now — and it was blocked on this: the console never imports metaannot, so
everything it knows it gets by shelling out and parsing stdout, and `doctor`
had no stdout a program could read.

The scope was fixed before anything was written, and both halves of it were
arrived at the hard way:

> Doctor checks that a required input EXISTS, is the RIGHT KIND of thing, and
> is NOT EMPTY. It does not parse. A file that passes doctor can still be wrong
> inside, and doctor says so.
>
> Doctor fails an input exactly when an enabled stage will die on it.

Holding the second sentence turned up **twenty-four** live false verdicts, in
four waves, and every one of them has an entry under Fixed. The fourth wave
also turned up something worse than a wrong verdict, and it has its own group:
four states of a path on which `doctor` printed NO DOCUMENT AT ALL, one of
them by never returning.

Seven were in the shipped `doctor`, all of them the same bug: the verdict was
one boolean — `ok &= r["ok"]`, `ok &= good or label == "gff"` — folded over
every check, with no idea which stage cares about which. Correcting those
changes `doctor`'s exit status for real configs, which is the point — the first
of them failed the ordinary annotate-without-MS-quant run that `run` exits 0
on.

Seven more were in the first draft of this document, found by reading every
row against the stage function it describes a second time, and they are the
more interesting half: a row that names the stages it kills can be
confidently, specifically wrong in a way one boolean never could. Two claimed
stages that return or fall back before they ever open the file — `stage_unipept`
with `unipept.result` set, and the peptide readers, which re-read with the
peptide-only reader when the full one refuses. Two branched on
`found.present`, or on wording written for "absent", while asserting they had
checked the KIND. One guarded a check on the wrong config key, so a stage that
dies went unnamed. One probed for a backend the stage stopped requiring. And
one produced no row at all — the most expensive false pass there is, because
there is nothing in the document to be wrong and so nothing can notice.

Five more survived both readings and were found by RUNNING every row instead:
one of them printed no document at all, because a directory at `manifest`
reached `read_manifest`, which opens the path with `open()` and raised an
`IsADirectoryError` that the `except StageError` beside it does not catch. A
command that prints nothing has no contract; every other row in the document
went down with the one that could not be read.

**The counts in these paragraphs are the FALSE VERDICTS, and they are not the
number of entries under Fixed.** That section also carries one attribution
that was never a verdict, twenty-eight corrections to what the document says
about ITSELF, one ENGINE defect a verdict uncovered, four states of a path
that produced no document, ten defects the sixth reading corrected, twelve
defects the seventh reading corrected, seventeen defects the eighth reading
corrected, and seven defects the ninth reading corrected —
eighty entries that change no config's
`doctor` verdict. The engine defect was missing from this sentence while it accounted
for thirty-two of thirty-three; a test now sums every clause of it against the
groups that are not false verdicts, so a new group cannot be left out of it. The sentence here used to
claim the two numbers were the same, which they never were — so the section is
now grouped under headings that carry their own counts, and a test in
`tests/test_docs.py` counts the entries under each heading against the number
written in it. The tests behind the verdicts themselves are one per verdict,
each driving the stage — or, for `proteins_faa`, the `run` command — into the
state the row claims.

### Added

**The report now says out loud when a bin reaches the statistics with nothing
in it (#5).** Every number was already there. `retention by bin` printed the
zeros, the join stage already warned that the proteins `min_features_per_protein`
dropped are absent from `annotated_quant.tsv` and from every report table, and
`DEFAULT_CONFIG` already said what a default of two costs the population this
tool exists to study. Nothing escalated any of it, while the same document
raised a GATE over findings orders of magnitude smaller, so the reader had to
notice a cell in a printed tibble and work out what it implied. The retention
chunk now draws the conclusion instead. A bin that was quantified and then
filtered out entirely gets a GATE naming it and the count behind it, and saying
that no table, no enrichment and no shortlist below is about it. The KO-less
bins get a second GATE of their own, stating how much of that whole population
the statistics cover — the sentence the issue asked for, and the one this tool
exists to be able to make. A NOTE then names the knobs that decided it:
`analysis.min_valid_per_group` with its value and the design it was applied
over, `analysis.min_plexes` where an isobaric run applies that as well — the
retained set is `keep_valid & keep_plex`, so on the suite's own TMT fixture the
plex filter removes ten of fifty against `min_valid_per_group`'s none, and a
line calling `min_valid_per_group` "what the counts above measure" would be
untrue there — and `min_features_per_protein`, which chose the denominators
before the report saw anything. That key is TOP-LEVEL and is named without a
prefix: `join` is a stage name and `run.join` is a boolean, so a reader who
wrote the `join:` block an earlier draft of this NOTE implied was told
`unrecognised key 'join'` and had the setting silently ignored. Every dotted
config path the document prints is now checked against `DEFAULT_CONFIG` by the
suite, in both directions. On the run this was found on it reads
`0/170 3d_duf_only and 0/810 4_dark quantified group(s) reach the model`,
followed by the KO-less coverage, which was eleven of 3,894 quantified groups.

**The model the fraction is read against is the WHOLE model.** The per-bin
analysis drops the quantified groups with no annotation row — they are not a
bin, and they are reported by their own NOTE further up — but they ARE fitted,
and being the well-covered ones they pass `min_valid_per_group` where the
sparse dark ones fail. So the KO-less GATE counts them, `sum(cov_tab$n_tested)`
rather than the binned subtotal, and the retention table's percentage column is
`pct_tested` rather than `pct`, since with `drop_zero_variance` set `n_kept` and
`n_tested` differ and a name after neither leaves the reader to guess. Read
against the binned subtotal the sentence said `1 of the 34 group(s) in the
model` on a page that had already printed both the unbinned groups it left
out and a model of 48 — an error the flattering way, in the loudest sentence
the report has.

**Neither default moved**, and that is deliberate rather than unfinished.
`min_features_per_protein` is already 1, already carries the warning, and the
run that produced this issue overrode it knowingly. Making the consequence
impossible to miss is a different thing from preventing the choice.

**The escalation is by DENOMINATOR**, which is the only part of this that
needed deciding. A bin that retained nothing because it HELD nothing is not a
finding: `3p_profile_only` is fed only by `hhblits` and `jackhmmer`, so on a
default run it is quantified-nothing, tested-nothing, for ever, and a gate that
fires on it every time teaches the reader to skip the one line the report most
needs them to read. So a bin with no quantified group is silent — its zero is
already in the bin composition table — a bin under `COVERAGE_MIN_N` gets a NOTE
rather than a GATE, because "none of them" over a handful of proteins is an
anecdote and not a statement about a population, and a bin at or above it gets
the GATE. The percentage floor, `COVERAGE_MIN_PCT`, is applied to the KO-less
bins TOGETHER and to no single bin: whether one small bin kept a tenth or a
fiftieth of itself is noise, while what fraction of the KO-less proteome the
statistics cover is the claim the document is for. It is deliberately not a
test of whether a bin is depleted RELATIVE to the others — with that run's
missingness the whole proteome kept about one group in eighty, and against that
baseline a bin of several hundred coming out at zero is not even surprising.
Surprise is the wrong question; coverage is the question, and its answer is a
count rather than a p-value.

**The ratio model and the effector shortlist say it too**, rather than
silently having nothing to rank. A ratio model whose usable set contains no
KO-less protein now says so beside the count it already printed, because that
model exists on the argument that a KO-less protein is likelier than a mapped
enzyme to be riding its organism's abundance — run without one, it is answering
that question for the proteins it was never the argument for. It says WHICH of
two different things happened, because they have different causes and
different fixes: either no KO-less group reached the model at all, which is the
coverage failure gated above, upstream of this model and nothing to do with
taxonomy, or the KO-less groups that did reach it have no usable taxon, which
is a limit of the taxonomy and is where `taxon_min_proteins_for_factor` or a
better assignment would help. And it is under the same floor as everything
else here: the second case is a GATE only over `COVERAGE_MIN_N` KO-less groups
in the model and a NOTE under it. Without that floor it fired on a knit with
every bin fully retained and nothing wrong with it, which is the exact noise
the denominator rule exists to keep out.

An empty shortlist now prints the population it was drawn from: `0 candidates`
reads as a negative result, which is exactly what it was on the first real run,
and is indistinguishable on the page from a shortlist that had nobody to rank.
Only the second is a GATE. Where the KO-less population was tested, the NOTE
names how many groups it was drawn from; where nothing KO-less reached the
model but too few were quantified for the GATE, the NOTE says that nothing was
rankable. It never says the list is "a statement about those 0", which is the
missing population dressed as a weak negative result and the reading this whole
change exists to prevent. What the denominator changes is the TIER, never the
claim — which is the same defect, one report section over, as the one the
v0.2.0 entry corrected about what an empty shortlist means.

The rule is pinned by a real knit over a run reproducing the issue's shape and
by a knit of a healthy run that must stay quiet, not by a string match on the
template, and by the report's own coverage block, its ratio-model escalation
and its empty-shortlist branch, each lifted out of the document and run
against the counts that reach it — the issue's per-bin ones included. The
healthy-knit test asserts what the report DOES say as well as what it does
not, because a report that cannot raise a line satisfies every assertion that
it was not raised.

**`--force-unlock-live`, and a `--force-unlock` that refuses a holder it can
see running.** Registered on `run` and on `all`, because `all` is the command
the tutorial leads with and a flag the refusal names has to exist wherever the
refusal can fire. The process-table reading is now three answers rather than a
boolean — `PROVEN_ALIVE`, `PROVEN_DEAD`, `UNPROVABLE` — and the two callers of
it want opposite defaults: `_holder_is_alive()` stays `!= PROVEN_DEAD`, branch
for branch, so unprovable still means alive and nothing about when a lock is
reclaimed has moved; the refusal reads `== PROVEN_ALIVE` only. That asymmetry
is the whole care of it. A lock written on another node of the array cannot be
disproved from here, and a refusal keyed on "not provably dead" would refuse
there — removing the documented escape hatch on exactly the machine it exists
for, and on every lock on Windows besides. The message prints what an operator
would otherwise go and assemble: pid, host, when the run started, `last_seen`
and how long ago that was, which stages are recorded `running`, the command
line, and the `ps -p` to run. `_run` stays advisory throughout: it is read as
text for a person and nothing branches on it.

**`metaannot.py doctor --json`.** One object on stdout and nothing else,
exactly as `describe --json` does it. The header is describe's header key for
key, and `requirements` is describe's `requirements` array element for element
— the same `requirements(cfg, p)` call, not a re-shaping of it, so `cmds`,
`size_gb`, `disk_gb` and `note` have one home and a check reaches them through
`requirement_id`. What doctor adds is `checks`: one entry per line the printed
report prints, in print order, under the same `== ... ==` headings, with
`detail` being **the printed sentence itself** rather than a second wording of
it. The two renderings cannot drift, and a test asserts every sentence in the
document appears in the terminal output.

**Every check names which enabled stages it kills.** `blocks` is a list of
stage names that joins straight to `describe --json`'s `stage_names`;
`blocks_commands` carries `run`, `report` and `object` for the failures where a
whole command refuses and no stage dies; `degrades` carries the stages that run
anyway and produce less. The exit status is a fold over those:
`status == "fail"` if and only if `fails_reason` is set, and
`fails_reason == "stage_or_command_dies"` if and only if
`blocks`/`blocks_commands` is non-empty. `verdict.rule` states it as a sentence
in the document, so a consumer that understands nothing else can still compute
the status a newer doctor produced. `fails_reason` has **two** values, not the
three the design sketch proposed: every "config contradiction" this vocabulary
can express turned out to kill something outright — `run.taxonomy` without
`run.unipept` is `stage_taxonomy` dying on a `unipept_lca` nothing wrote — so
it is reported as `stage_or_command_dies` with that stage, or that command,
named, and a value nothing emits is a branch nobody exercises. The value is
named for both halves because it was first named for one; see the entry under
Fixed.

**`requires` and `degraded_by` on every stage**, and on `describe --json`'s
per-stage object. This is the attribution `blocks` is built from, and it could
not be read off `requirements()`, which gates on `run:` flags with inline
boolean logic (`if R.get("pfam") or R.get("dbcan") or R.get("ncbifam") or
R.get("jackhmmer")` for hmmer) and so knows that *something* wants hmmsearch
but not which stage dies without it. Adding it to `describe --json` does **not**
move `DESCRIBE_VERSION`: the rule at that constant is that a version moves when
a key is removed or its meaning changes, and nothing was. It is a deliberate
decision rather than an accident, which is what
`test_describe_emits_every_field_a_stage_dict_carries` and the exact key-set
pin in `tests/test_config.py` exist to force.

**`scope.statement`** — the first scope sentence above, as a literal string the
engine authors and a front end renders verbatim. "and doctor says so" is half
the rule, and a console that has to write that sentence itself becomes its
second home.

**A per-check `depth`** (`config`, `existence`, `kind`, `header`, `parsed`,
`probe`) and a `caveat` on every row that went past a file's existence.
"doctor does not parse" was already false as an unqualified claim — the DIAMOND
usability check, each plex's TMT annotation file, `read_manifest` on the
manifest itself and `pd.read_csv(nrows=0)` on the quant header are four checks
that read into files, all of them shipped and all of them legitimate. Now each
says so on its own row. The enum is documented as a **ratchet rather than a
menu**: a fifth needs a reason written down at the constant, and the ratchet
bites in both directions — see Fixed.

**`found.kind` names the state instead of leaving it to arithmetic**: `file`,
`dir`, `empty_file`, `empty_dir`, `symlink_broken`, `absent`, `unset`, `other`
and `unreadable` — plus `null`, which is a value and means "this build did not
compute it"; see Fixed for both of those. A dangling symlink has always been
rejected "for the right reason by accident" — `os.path.exists` follows the link
and answers False — and a volume that is not mounted is how that happens here.
It is now reported as what it is, in the printed report as well.

**`expect` says what the right kind of thing is, and which key decided.**
`expect.derived_from: ["quant_format"]` on the quant table, because `file`
versus `dir` is config-dependent: only `fragpipe_tmt` reads a directory. A
front end renders the reason without carrying its own copy of the eight-format
table.

**`remedy`: `auto` / `manual` / `config` / `input` / `none`.**
`requirements()[].manual` is overloaded — `dbentry()` sets it to
"no path configured" for a key nobody filled in, while SignalP 6 and
InterProScan set it to a licence and a version-specific distribution. The first
wants an edit box, the second a link-out and structurally no button. The split
is in `remedy`; `requirements()[].manual` is left exactly as it was, and so is
the printed `MANUAL` mark, which an operator reads as "doctor will not fetch
this" and which is true of both.

**`totals` says what it leaves out.** `counted` / `not_counted` (MANUAL items
are in neither figure, so both totals are a floor rather than the whole job)
and `unsized` — the ids where the `0.0` in `requirements[]` means "not
estimated" rather than "free", which is how a large Java distribution came to
look like no download at all. `accuracy: "indicative"` and
`basis: "hand-maintained"` are there to branch on, and the sentence CLAUDE.md
tells an operator to repeat about those numbers lives in `totals.detail` rather
than in whatever is rendering them.

**`DOCTOR_VERSION`, with the clause `DESCRIBE_VERSION`'s rule does not have.**
Same rule — bumped when a key is removed or its meaning changes, never when one
is added — plus: **adding a value to a closed enum is a meaning change and is a
bump.** The closed sets are `DOCTOR_STATUSES`, `DOCTOR_REMEDIES`,
`DOCTOR_FAIL_REASONS`, `DOCTOR_DEPTHS`, `DOCTOR_COMMANDS`,
`DOCTOR_FOUND_KINDS` and `DOCTOR_EXPECT_KINDS` — seven, the last two added
under Fixed below after they were found outside the set; `finding`, the stage
names and every `detail` are open and may grow freely. Pinned by an exact
key-set test in the style of the one that protects `describe --json`, which
exists because a version constant that pins nothing let nine documented fields
be deleted while the suite stayed green.

**`doctor` now checks `contigs_fna`.** It never did, though `stage_smorf` dies
with "run.smorf needs contigs_fna pointing at the assembly" — so any
`run.smorf: true` config passed doctor and died at run time.

### Fixed

104 entries, in fifteen groups. Each heading carries its own count and
a test counts the entries under it, because the prose above this section has
already claimed a number twice that the list below it did not have — and the
reason the two never matched is that "false verdict" and "entry" are not the
same thing: a correction to what the document says about itself changes no
config's exit status and is grouped separately here.

#### Seven false verdicts in the shipped `doctor`

**A missing `quant_table` no longer fails an annotate-only run.** `stage_join`
logs "quant table not found, skipping join" and *returns*, and `run.join` is
true in `DEFAULT_CONFIG`, so an ordinary annotate-without-MS-quant config
completed and exited 0 while `doctor` exited 1 on it. It is now a warning that
names `join` in `degrades` — and it becomes a failure the moment `run.unipept`
or `run.taxonomy` is on, because `peptide_features()` reads the table and has
no skip branch. Same file, same check id, opposite verdict: severity is a
property of the config, not of the input. A **zero-byte** quant table is fatal
either way, because `os.path.exists` is true so `stage_join` never takes its
skip branch and dies in the reader instead. Missing is survivable; empty is
not, and only a per-state check can say that.

**A missing DIAMOND database is a warning, not a failure**, and the fatal case
is the opposite one. `stage_diamond` logs "diamond database missing, skipping"
for each absent database and searches the rest, writing an empty `.done` marker
when none of them exist; what it dies on is a database that is PRESENT and
cannot answer — `diamond_db_check`'s refusal, which is a zero-byte `.dmnd` from
a failed `makedb`. Those are now two checks about one file, with different
depths and different verdicts. On the default config this alone was four
spurious failures.

**A set-but-missing `gff` now fails.** `cmd_doctor` did `ok &= good or label ==
"gff"`, so a `gff` pointing at nothing printed `MISS` and still exited 0 —
while `stage_context` calls `die()` on exactly that. An *empty* `gff` is
correctly a warning: the stage logs, writes an empty context table and returns.

**`emapper_precomputed` is no longer checked when `run.eggnog` is off.** It was
folded into the verdict with no guard, so a stale path failed doctor even when
the only stage that reads it was not running.

**A `manifest` that is set and absent is contingent, and says so.**
`read_manifest` is reached from `read_feature_table` and directly from
`stage_join` — there is no `read_protein_table` in this tool; the
protein-level path is inline in `stage_join` — and both are downstream of
`stage_join`'s early return — so with the quant table gone, nothing opens the manifest. It is
a warning with `depends_on: ["input:quant_table"]` and a sentence saying it
turns fatal the moment the quant table is restored, which a test drives both
ways.

**The `== tmt ==` block no longer fails a config nothing reads.** It ran on
`quant_format == "fragpipe_tmt"` alone and every `MISS` inside it set the
verdict, so a TMT config with `run.join`, `run.unipept` and `run.taxonomy` all
off exited 1 — on a config where no enabled stage opens the quant tree at all.
It is now skipped, with one line saying why.

**`taxonomy_source` with `run.taxonomy` off is a warning.** `resolve_taxonomy()` logs
"falling back to eggnog" and carries on, so nothing dies; doctor was failing it.

#### One attribution, which corrected no verdict at all

The `== R ==` block is new in this release, so there was no shipped verdict to
correct. What it fixes is the ATTRIBUTION — where a failure that kills no stage
is recorded — and it sat above the wave marker for a release with nothing
saying it was not one of the seven.

**A missing required R package blocks `report` and `object`, which are
commands and not stages**, and each of the fourteen packages is its own row
with its own install line — `BiocManager::install` or `install.packages`,
decided by Bioconductor membership — because a preflight checklist offers its
affordance per row.

#### Seven false verdicts in the first draft of this document

Verdicts this document itself got wrong, each found by reading the row against
the stage function it describes and fixed with a test that drives that
function.

**A `taxon_rank` with `run.taxonomy` off asked for no taxdump at all.**
`requirements()` gated the `ncbi_taxonomy` entry on `run.taxonomy`, so on
`run.join: true` + `taxon_rank: genus` the row did not EXIST — and
`requirement_effect()`'s taxon_rank condition, which is written for exactly
this case, was reading an entry nothing had produced. `stage_join` →
`resolve_taxonomy` → `collapse_taxon_rank` dies there with `taxon_rank='genus'
needs db.ncbi_taxonomy`, and `collapse_taxon_rank`'s own log is what walks an
operator into that config: on every run with no rank it explains what raw seed
taxids cost and tells them to set one. The gate is now `run.taxonomy or
(run.join and taxon_rank)`, fixed where the row is PRODUCED rather than where
its effect is read, and tested in both directions.

**A zero-byte or directory `unipept.result`, and the same for the pept2lca
cache, passed both Unipept rows.** They asserted `depth: "kind"` — "also asked
what kind of thing it is and whether it is empty" — computed that answer, and
then branched on `found.present`, which a zero-byte file and a directory both
satisfy. `stage_unipept`'s own test is `os.path.exists()` too, so it takes the
ingest branch and `read_unipept_result` dies on "Could not determine
delimiter"; a zero-byte CACHE is worse, because `getsize > 0` means it is not
adopted, every peptide is still to do, and with `allow_http` false the stage
dies on "the cache is incomplete" — the very message that row exists to
pre-empt. This was the first scope sentence broken by doctor's own data.

**The `== manifest ==` block claimed stages that never open a manifest.** It
asked `quant_consumers()`, which is "who names the quant table", where the
question is "who calls `read_manifest`" — so with `unipept.result` set, the
manifest gone and join and taxonomy off, `doctor` exited 1 naming `unipept`
while `run` exited 0. `read_manifest` is reached from two places only:
`read_feature_table`'s fragpipe_peptide/fragpipe_ion branch and `stage_join`'s
own diann/fragpipe branch. On the MSstats formats the design is in the table,
and `read_fragpipe_tmt` logs "`manifest` is ignored for quant_format
'fragpipe_tmt'" and reads the annotation files instead — so on four of the
eight formats a missing manifest kills nothing at all, and it now says so on
one row instead of failing the config.

**The `== tmt ==` refusal rows claimed the two stages that fall back.**
`peptide_features()` catches `StageError` from the full reader and re-reads
with the peptide-only reader — "Quantification is still unavailable, only the
taxonomy work continues" — and `read_fragpipe_tmt_peptides` never opens an
annotation file, resolves a reference or checks a sample name for collisions.
So a reference conflict, an unreadable annotation, a duplicated sample name and
an unmatched reference are join's death and nobody else's, unless
`peptide_only_reader: never` says refuse rather than fall back. Three lists
come out of that where one was being used: who opens the tree at all (a
missing plex directory or level file is fatal to every reader), who dies when
the FULL reader refuses, and whose OUTPUT a tree that reads perfectly but
describes a ragged design degrades — which is join alone, because the taxonomy
stages take peptides out of it and never touch an intensity. On a config with
join off those rows are now warnings, which is also what stops them tripping
`_check()`'s own "a failure names what it costs" invariant and taking the
command down with an AssertionError.

**`req:esmfold` failed a host where ESMFold runs.** `requirements()` probed
`_pyhas("esm")` while `stage_esmfold` tries fair-esm and THEN transformers and
dies only if both fail — and on a current card the transformers route is the
only one that builds, because `fair-esm[esmfold]` needs an openfold pinned to a
2022 commit whose CUDA kernels do not compile against a modern toolkit at all.
The requirement now accepts either backend, its install line names
transformers, and its note still names fair-esm as what the stage prefers when
it is importable.

**A directory at `proteins_faa` claimed that `run` refuses.** Absent, dangling
symlink and directory were one branch whose wording was written for "absent":
`detail` said "`run` refuses before any stage starts" and the row carried
`blocks_commands: ["run"]` and `depth: "existence"`. `cmd_run`'s guard is
`os.path.exists()`, which is TRUE for a directory — `run` does not refuse, it
schedules everything and each stage dies as it opens the path — so the one
machine-readable claim a preflight screen would gate the `run` button on was
false, over a row that had stat'd the path and knew better. The same wording
bug (with the verdict right) was on `gff` and `contigs_fna`, whose stages guard
with `os.path.exists` too. And `_finding()` answered `empty` for an EMPTY
DIRECTORY where a file was expected, indistinguishable in the document from a
zero-byte file — which on `gff` is the difference between a warning and a
failure. It now takes the expected kind and answers `wrong_kind`.

**`taxonomy:format` was guarded on `unipept.result` for both stages.** Right
for unipept, which returns after ingesting the export and never reaches its own
format check; wrong for taxonomy, whose `die("the taxonomy comparison needs
peptide-level input")` is guarded on nothing at all — it reads `p.unipept_lca`
whoever wrote it and then calls `peptide_features()` itself. A protein-level
`quant_format` with `unipept.result` set and `run.taxonomy` on was a false pass
on a stage that dies in its first seconds.

#### Six claims the first draft made about ITSELF

Not verdicts, but the half a consumer cannot check:

**A row's `caveat` is now derived from what the check opened, not declared
beside it.** Both DIAMOND usability rows hard-coded `depth: "header"` and "read
from the DIAMOND header and a sequence-length profile; the sequences themselves
were not parsed", whatever `diamond_db_profile()` had actually done. With
`diamond` not on PATH — the ordinary state the first time anyone runs doctor —
it falls through to `_fasta_lengths()`, which reads up to 200,000 FASTA
records, sequence lines and all, and `motif_seed_evidence()` then reads 200
deflines; the caveat asserted the opposite of what had just happened, on the
one field whose whole job is to qualify "doctor does not parse". A zero-byte
`.dmnd`, meanwhile, is refused on `os.path.getsize` alone and opens nothing.
`diamond_db_check()` now returns what it read as tokens and doctor turns those
into the depth and the sentence, so the three cases say three different things.

**`depth` is a ratchet in BOTH directions.** It was enforced upward only —
`header` and `parsed` must carry a caveat — and nothing stopped a row claiming
it had looked deeper than it had, which every "no path configured" requirement
row was doing with `depth: "kind"` over a path nobody had named. `_check()` now
refuses an `existence`/`kind` row with a null `found`, and those rows report
`depth: "config"` and a null `found`, which is what "this build did not compute
it" is supposed to look like.

**`expect.because` said "found on PATH by `have()`" for two requirements that
are not on PATH.** InterProScan is tested with `_exists(db["interproscan_sh"])`
and ESMFold with `_pyhas()`, so a preflight screen told an operator to put
InterProScan on PATH when the fix is a config key. Both now say what was
actually probed, carry the config key in `config_keys`, and — for InterProScan
— go through `dbentry()`'s "no path configured" marker, so an unfilled
`db.interproscan_sh` gets `remedy: "config"` (an edit box) instead of landing
on the licence side of the split `_remedy_for()` exists to make. A path that IS
set and not there is still `manual`: that one really is a version-specific Java
distribution.

**`expect.members` is populated, and says why where it is empty.** It was
never passed a value anywhere in the file, though the cases it was written for
were all right there: the four `hmmpress` siblings `_pressed()` requires,
`foldseek_target`'s `.dbtype` and `.index`, the `nodes.dmp` inside a taxdump,
the `eggnog.db` inside an eggNOG data directory. It carries the engine's OWN
test, never the files a tool would merely like to find — which is why
`hhblits_db`'s is empty and its `because` says so: `_prefix_exists()` accepts
any non-empty sibling of the stem.

**`found.kind` and `expect.kind` are inside the versioned set.** `_found()`'s
docstring and the README both tell a consumer to switch on `found.kind` and
never to compute `bytes > 0` for itself, and `expect.kind`'s values were
enumerated in no constant, no document and no test — so adding `fifo` or
`unreadable` to one, or a sixth kind to the other, would not have been a bump
and a console's switch would have fallen through in silence. Both are now
closed sets beside the version constant under the same rule, enforced in
`_check()` where the rows are made. Closing a vocabulary is not itself a bump:
it adds a promise rather than changing a meaning.

**`expect.kind: "probe"` on a config row is now `"setting"`.** It collided with
`depth: "probe"` and meant the opposite of it: those rows carry `depth:
"config"` and `found: null` — they consult nothing outside the config — while a
`depth: "probe"` row asks the HOST. `probe` now means the same thing in both
fields, and a test asserts that a row whose `expect.kind` is `probe` has
`depth: "probe"` too.


#### Five more false verdicts, found by running every row

The third wave, and the reason it exists: the first two were readings. These
were found by driving every row and comparing what `doctor` said with what
`run` did on the same config.

**A directory at `manifest` produced NO DOCUMENT AT ALL.** `_manifest_checks`
branched on `found.present`, which a directory satisfies, and fell through to
`read_manifest(path)` guarded only by `except StageError` — while
`read_manifest` opens through `opener()`, which is `open()`, so a directory
raised `IsADirectoryError`: `rc=1`, a traceback on stderr and nothing on
stdout. Every other row went with it — the tools, the databases, the resources
and the whole R block — and a consumer got no verdict about any of them. The
manifest is the one input the earlier directory fix did not reach. It now has
the same split the others have, made BEFORE `read_manifest`, with two
different costs behind it: an EMPTY manifest dies inside `read_manifest` with a
`StageError` that `peptide_features()` catches and re-reads around, so it costs
`join` alone, while any other wrong kind raises an `OSError` that nothing
catches, so it costs every reader. The `except` is widened to `OSError` as
well, as a backstop the branch above should now make unreachable — an
unreachable backstop costs one row in a document that survives, and a missing
one costs the document. A test walks every configured input path through all
nine states — absent, empty, directory, dangling symlink, FIFO, socket, device
node, unreadable directory and unreadable file — and asserts stdout parses as
JSON in all fifty-four cells.

**`req:ncbi_taxonomy` failed a `join` that skips.** `stage_join`'s first act is
`if not os.path.exists(quant_table): log("quant table not found, skipping
join"); return`, before `resolve_taxonomy` and before anything else it
requires — so with the quant table absent, nothing `join` declares can kill it.
With `taxon_rank: genus`, no taxdump and no quant table, `doctor` exited 1
naming `req:ncbi_taxonomy` while `run` exited 0 on "done: 4 run, 0 adopted, 17
skipped", and `run.taxonomy` is off by default, so that is the MORE COMMON half
of the configs. doctor already held the contingency twice —
`input:quant_table` says `degrades: ["join"]` on the same document, and
`_manifest_checks` gates its whole verdict on the same `os.path.exists` — and
applied it nowhere near the requirement rows. It is one predicate now, and the
rows it decides carry `depends_on: ["input:quant_table"]`.

**`input:contigs_fna` claimed `blocks: ["smorf"]` whatever is installed.**
`stage_smorf` dies at the top only for an absent or unset assembly; past that
guard the BINARY decides. With an ORF finder on PATH a directory is handed to
it and the stage dies; with none installed the path is never opened at all —
the stage logs "neither smorf(inder) nor macrel is installed", writes an empty
`smorf.faa` and returns, and `run` exits 0. The row's comment about the
inversion against its own tools is true of the absent branch and had been
copied onto this one. Measured on a machine with no `smorf`, no `smorfinder`
and no `macrel`: `doctor` 1, `run` 0 — and the test that pinned it ran on that
same machine, so the suite was green-lighting a verdict false in its own
environment. The verdict is read off the probe now, and the test drives both
sides of it with a stub on PATH rather than trusting the laptop.

**Two tool probes disagreed with their own stages.** `req:kofamscan` probed
`have("exec_annotation")` while `stage_kofam` accepts `exec_annotation` OR
`kofamscan`, so a host with only the second was told `blocks: ["kofam"]` about
a stage that runs. `req:smorf` probed `have("smorf") or have("macrel")` while
`stage_smorf` also accepts `smorfinder`. Both probes now accept what their
stage accepts, `smorf_tools()` is the single probe the stage, `requirements()`
and `doctor` all share, and a test asserts each documented name is enough on
its own.

**A FIFO at `proteins_faa` claimed that `run` refuses.** `_found()` had no
branch for anything that is neither a regular file nor a directory, so a FIFO,
a socket or a device node came out as `kind: "absent", present: false` — and
the row then said `blocks_commands: ["run"]` and "run refuses before any stage
starts", while `cmd_run`'s `os.path.exists()` is TRUE for one and schedules
everything. The identical false claim a directory used to produce, one path
state over. `present` is `os.path.exists()`'s own answer now, and the state has
a name.

#### One engine defect a verdict uncovered

**`stage_context` wrote a context table with NO COLUMNS.** The `input:gff` row
says an empty GFF "writes a table with no rows, exactly as the unset case
does". It did not: the unset branch writes
`pd.DataFrame(columns=["protein_id"])`, while a GFF that parsed to nothing fell
through to `pd.DataFrame([]).to_csv()`, which writes one newline and no header
at all — so `context.tsv` was one byte, `nonempty()` was true, and
`parse_context`'s `pd.read_csv` raised `EmptyDataError: No columns to parse
from file`. `integrate` ALWAYS runs, so that killed every run with
`run.context` on whose GFF produced no rows: `doctor` 0, `run` 1, "FATAL stage
'integrate' failed: No columns to parse from file". A HEADER-ONLY GFF with no
CDS records is the same defect through a different door, and that one is
outside `doctor`'s scope entirely — the file is neither missing nor empty, so
`doctor` reports `ok` for it and must, which is why the fix had to be the
engine's. Both branches write the column now, and a test asserts the two files
are byte-identical.

#### Six claims the third reading corrected

**`verdict.rule` dropped the whole `blocks_commands` class**, and it is the
statement of the rule a consumer actually parses, because it is IN the
document. It said a check fails exactly when an enabled STAGE dies on it, or
the declared `setting_ignored` exception — neither of which covers
`r:package:limma`, which is `status: "fail"`, `blocks: []`, `blocks_commands:
["report", "object"]` and no stage anywhere. `TUTORIAL.md` said the same thing;
`doctor()`'s docstring and the README's closing sentence had it right, so three
renderings of one rule disagreed and the loosest was the machine-readable one.
The reason value is renamed from `stage_dies` to `stage_or_command_dies` in the
same pass: a name that has to be explained away is a name that gets mis-copied,
and it had been, three times. That is a renamed value in a closed enum and
would be a `doctor_version` bump under this file's own rule — except that
`doctor_version` 1 is in this same unreleased change set and has never shipped,
so there is no consumer to break. The config sweep could not have caught
this — and it is `CONFIGS` in `tests/test_doctor_json.py`, a named list the
tests parametrise over, not the 768-way product this sentence used to claim:
`no_r_env()` strips `Rscript` from PATH for every subprocess `doctor` call, so
no R failure is ever in a swept document and `blocks_commands` never appeared
in one. The sweep carries a `no_faa` config now, and the R rows get
the sweep's own invariants in process.

**`expect.because` was still assumed for every `tool` row.** The blanket branch
returned "found on PATH by `have()`; every stage that uses it dies first" for
every tool — and on `req:smorf` it returned that beside `blocks: []` and
`degrades: ["smorf"]`, so one object contradicted itself. This is the defect the
`interproscan` and `esmfold` fixes above set out to kill, left standing on the
branch they did not take. Both halves are derived now: the probe from a table
driven against `requirements()` one executable name at a time, and the cost
from the same `requirement_effect()` call that fills `blocks` and `degrades`.

**`found.kind`'s null is documented, and the two missing states have names.**
The README enumerated seven values "so a consumer never has to do arithmetic on
bytes" while `_check()` explicitly admits one more — its guard is `not in
(None,) + DOCTOR_FOUND_KINDS` — and every requirement row, the CUDA probe and
the whole R block emit a non-null `found` whose `kind` is `null`, which in the
default document is the MAJORITY of the rows. The null is
deliberate and means "this build did not compute it"; the defect was that
nothing said so while both the constant and the README told a consumer to
switch on the field. `other` and `unreadable` join the set at the same time —
the comment beside the constant had named both as values that "would not have
been a bump", while `_found()` was answering `absent` for both.

**The depth ratchet did not bite downward for `skip` rows.** `config` is
defined as "consulted nothing outside the config", and every skipped input row
was stat'ing its path to fill `found` and then reporting `config` — while a TMT
root row was claiming `kind` over a `quant_table` nobody had configured. The
guard is symmetric now in both senses: an unset `found` is not evidence of a
stat, and a real one forbids the `config` claim. `_stat_depth()` derives it, so
no row writes its own.

**The count of checks that read INTO a file is derived.** The document said
"three checks" in PUBLISHED text — the `caveat` on every
`tmt:<plex>:annotation` row — while its own `scope.statement` said four, and
"go past a header" was not even the same predicate, since the quant-table check
reads the header line and nothing else. Four places were counting by hand.
`DOCTOR_DEEP_CHECKS` is the tuple now; nothing counts it, and a test asserts no
hand-written count of "checks" is left in the file.

**`read_protein_table` and `taxon_map` do not exist.** Four sentences added
by this change set named a `read_protein_table` that has never been in this
tool — `read_manifest` is reached from `read_feature_table` and from
`stage_join`, whose protein-level path is inline `pd.read_csv` plus column
detection — and two more named a `taxon_map`, where the function that logs
"falling back to eggnog" is `resolve_taxonomy()`. A test now walks every
backticked call this section writes and refuses one the tool and the suite do
not define, which is why the two invented names are written here without their
parentheses.

#### Five more false verdicts, from a differential sweep against `run`

The fourth wave, and the reason it exists: the three before it read the rows
against the engine. This one RAN both — every configured input through every
state a path can be in, `doctor` against the real `run`, and compared the exit
statuses. Four of the five below came out of a single cell of that table.

**A mode-000 input passed every check that branches on `kind == "file"`.**
`_found()` reached `unreadable` only when `os.listdir()` or
`os.path.getsize()` raised, and `getsize()` is a STAT: it needs search
permission on the parent directory and nothing at all on the file, so a
regular file nobody may open answered its own size cheerfully and came back
`kind: "file"`. Measured, four rows said `ok` while `run` exited 1 on
"Permission denied" — `proteins_faa` (`emapper`), `quant_table` (`join`),
`gff` (`context`) and `emapper_precomputed` (`emapper`). Shared FragPipe or
eggNOG output on a cluster is the ordinary way an input gets into this state,
and the remedy is a permission, which is exactly what `unreadable` exists to
say and could not: the state had a name in `DOCTOR_FOUND_KINDS`, a sentence in
`_present_kind_phrase()` and a `finding`, and was reachable for a DIRECTORY
only. `os.access()` now answers alongside `getsize()`, `PATH_STATES` in the
suite gained the `unreadable_file` it never had, and the four rows are tested
by driving `run` and asserting the stage it dies in is the stage the row
NAMES.

**`gpu:topology` was identical for all three values of the key it names.**
`stage_tmbed` maps `tmbed_use_gpu` onto three different command lines — `auto`
to `--use-gpu --cpu-fallback`, `true` to `--use-gpu --no-cpu-fallback`, `false`
to `--no-use-gpu` — and the row was `warn` / `blocks: []` /
`degrades: ["tmbed"]` for all three. Two of those are wrong in opposite
directions. TMbed tolerates a missing or failing GPU only under
`--cpu-fallback`, so under `true` on a CPU-only host every chunk fails, every
protein lands in `tmbed_failed.tsv`, and the stage dies on "tmbed finished 0
of N prediction(s)" — unless `tmbed_allow_partial` is on, in which case it
returns with an EMPTY prediction file, which is a degradation and not a death.
So the verdict is fatal, and which of the two it is depends on a SECOND key
the row never mentioned. Under `false` nothing is fatal, but the row still
advised "Set `tmbed_use_gpu: false`" to an operator who had already set it —
and the slowness there is not a fallback, it is the configuration. The three
command lines are `TMBED_GPU_MODES` now, read by the stage and by the row, and
a stub `tmbed` that refuses `--no-cpu-fallback` drives the fatal half.

**`input:contigs_fna` EMPTY asserted an outcome only the host decides.** It
said `stage_smorf` "hands the empty assembly to whichever ORF finder is
installed and writes an empty candidate list" — two claims welded together,
and which one happens depends on whether an ORF finder is on PATH: with none
the assembly is never opened at all. That is the exact fault the DIRECTORY
branch one state over had already been corrected for, left standing on the
neighbouring branch. It reads `smorf_tools()` now and carries
`depends_on: ["req:smorf"]` like its neighbour.

**A requirement row promised that restoring the quant table would make it
fatal, where it would not.** TWO conditions take `join` out of a requirement's
`blocks` — an absent quant table, and, for `ncbi_taxonomy`, an empty
`taxon_rank` — and the contingency sentence is about only the first.
`depends_on` was derived beside `requirement_effect()` rather than from it
("join is on, the quant table is gone, and join's `requires` tuple names this
id"), so on `taxon_rank: ""` it sent an operator to restore a file that would
change nothing. `requirement_effect()` takes the contingency as a parameter
now and the row asks it the actual question: would restoring the quant table
put `join` in `blocks`.

**The TMT level file's only test was `os.path.exists()`.** `tmt:<plex>:
level_file` fired when the file was absent and said nothing about any other
state — while `read_fragpipe_tmt`'s guard is the same `os.path.exists()`, so
it does not die there either. Driven one state at a time on that reader: a
DIRECTORY raises `IsADirectoryError` out of `header_columns()`, a chmod-000
file raises `PermissionError` there, an EMPTY one gets past it and dies in
`read_delim_table` with pandas' `EmptyDataError`, and a FIFO does not raise at
all — it blocks. None is a `StageError`, none is caught downstream, and all of
them cost every consumer of the tree. The row is a three-way split now, like
every other input in the file.

#### Four states of a path that produced no document at all

Worse than a wrong verdict, and the reason the fix below is structural rather
than a fifth branch. A command that prints nothing has no contract: every
other row — the tools, the databases, the resources, the R block — goes down
with the one that could not be read, and there is nothing in the document to
be wrong, so nothing can notice. Three previous rounds each fixed the site
that had just been caught.

**A FIFO at a plex's annotation file made `doctor` HANG INDEFINITELY.** Not
fail — hang. `tmt_annotation_path()` guards only `os.path.exists()`,
`read_tmt_annotation()` goes straight to `opener()`, which at the time was a
bare `open()`, and a read-only open of a FIFO with no writer used to block
with nothing to interrupt it. No document AND no
exit, which is worse than any traceback, because a caller waiting on the
process has nothing to time out against either.

**A chmod-000 annotation file, and a directory at the same path, printed a
traceback and zero bytes of stdout.** `PermissionError` and
`IsADirectoryError` are `OSError`s, and the `except` beside the call was
`StageError` alone — the identical shape as the `manifest` directory two waves
ago, at the one site that fix had not reached.

**A chmod-000 plex directory took the document down from an `os.listdir()`
inside an f-string.** The level-file row was building a helpful "`{pdir}`
holds …" sentence, and a directory nobody may list raises `PermissionError`
out of the `sorted(os.listdir(pdir))[:8]` in the middle of it. A sentence's
garnish may not be able to fail; it goes through a `_dir_sample()` that cannot
raise.

**A FIFO at `--config` hung every subcommand, not just `doctor`.**
`load_config()` guards `os.path.exists()` and `os.path.isdir()`, and neither is
false for a FIFO, so `run`, `describe` and `doctor` all blocked on the open.
It is a named refusal now.

And the fix is one gate, not a fifth branch. `regular_readable()` proves a path
with an `os.open(O_RDONLY | O_NONBLOCK)` and an `fstat()` on the descriptor:
`O_NONBLOCK` is what makes a FIFO return instead of wait, `fstat()` is the only
"is this a regular file" test that cannot be raced by the path changing under
a stat, and opening at all is the only test a permission cannot lie to.
`_deep_readable()` wraps it for `doctor` and hands back the sentence as well
as the verdict; every check that reads into a file goes through it, and each
keeps a wide `except` as the backstop for the state the gate did not
anticipate. A test walks EVERY path `doctor` can be pointed at — the six
configured input keys, `emapper_precomputed` among them, and the TMT tree's
plex directories, level files and annotation files — through absent, empty,
directory, dangling symlink, FIFO, socket, device node, unreadable directory
and unreadable file, asserting in all 81 cases that stdout parses as JSON and
that the process EXITS, with a timeout so a future hang fails the test instead
of wedging the suite. `emapper_precomputed` was outside that sweep while the
sentence said EVERY, which is the same shape of claim as the counts below: it
is in it now, and its nine cells were clean. A second test parses
`metaannot.py` and asserts the RULE: a function that emits rows and reads into
a file must call `_deep_readable()` and must carry the backstop.

#### Ten claims the fourth reading corrected

**A FIFO at an input does not make the stage DIE — it makes `run` HANG, and
the rows said "dies".** The differential sweep that found the mode-000 defect
above also ran `run` on a FIFO at `proteins_faa`, `quant_table`, `manifest`
and `gff`, one at a time, and it never returned on any of them: `read_fasta`,
`read_manifest` and `parse_gff` all reach `open()` eventually, and a read-only
open of a FIFO with no writer blocks. `blocks` is still the honest
machine-readable answer — the stage does not complete either way — but "not
one of them can get a FASTA out of this" sends an operator looking through a
log for a failure that is never coming, and a run left overnight on it looks
exactly like a long stage. The rows say what really happens now. `opener()` is
deliberately NOT changed: a FIFO with a writer on the other end reads
perfectly well, and refusing one would remove something that works. `doctor`
is the command that may never block, which is why only `doctor` guards the
open.

**"NULL IS AN EIGHTH VALUE", beside a tuple of nine.** The third consecutive
round to ship a wrong hand-written count, and the second of them inside the
comment on the very constant that enforces the set being counted. The comment
states no position at all now; `DOCTOR_NULL_KIND_NOTE` derives the ordinal
from `DOCTOR_FOUND_KINDS` and a test recomputes the README's.

**"the DIAMOND usability rows" were named in the null-`kind` list, in three
places, and never belonged there.** `db:diamond:<tag>:usable` is built with
`found=_found(path)` and carries a real kind, which makes it the document's
own counterexample to the sentence it was named in. The families that really
carry a null are `DOCTOR_NULL_KIND_ROWS` — `(prose, id glob)` pairs, not a
sentence — and a test drives a document containing every one of them and holds
the list against the rows that really carry one.

**The count of non-verdict entries accounted for thirty-two of thirty-three.**
"One attribution that was never a verdict and twelve corrections" left out the
ENGINE defect, which is also an entry that changes no config's verdict. The
sentence names every bucket now and a test sums its clauses against the groups
that are not false verdicts.

**`input:contigs_fna` wrote its own `depth`.** The present-and-wrong-kind
branch hard-coded `depth="kind"`, which made the previous round's "`_stat_depth()`
derives it, so no row writes its own" false — invisibly, because for a
DIRECTORY the two agree. They stop agreeing the moment `unreadable` is
reachable for a file, which is the fix above. `_check()` now enforces the
ratchet in that direction as well: a row carrying a real `found` may not claim
a depth DEEPER than its `found.kind` justifies without a `caveat` saying what
it read. That is the structural version — the branch is fixed, and the next
one cannot be written wrong without failing at the point the row is made.

**"The 768-config sweep could not have caught this."** There is no 768-config
sweep. `CONFIGS` in `tests/test_doctor_json.py` is a named list the doctor
tests parametrise over, and the number was the size of a cartesian product
nobody ever ran. The sentence names the list now, and a test checks the list
is still there.

**`smorf.faa` is `smorf_proteins.faa`.** The `contigs_fna` row told an
operator to look for a file at a path this tool has never written;
`Paths.smorf_faa` is `results/smorf/smorf_proteins.faa`.

**`expect.derived_from` named a `tmt.annotation_glob` that has never been a
config key.** The key is `tmt.annotation` — a pattern or a `{plex: path}` map —
and `derived_from` exists precisely so a front end can render "this wants a
file because <key> decided", which a key that does not exist cannot do.

**`TUTORIAL.md` said the `== R ==` block "probes four packages ... and checks
neither knitr nor pandoc".** It probes fourteen packages and checks `pandoc`, and has
since the release that widened it; the sentence that described the old probe
outlived it. This is what the count sweep below was for.

**Every count in prose is now DERIVED or PINNED.** Three consecutive rounds
shipped a wrong hand-written count, so the fix is the class: a test scans all
six surfaces a count can live in — metaannot.py's comments and docstrings, the
document `doctor --json` emits, the test suite's own comments and docstrings,
the README, the TUTORIAL and this file's Unreleased section — for
every "<number> <countable noun>", and requires each to be classified as
DERIVED (recomputed from the source, right there, and compared), MEASURED (a
fact about a dataset or a host: not recomputable, so the pin is that the
SENTENCE STILL EXISTS and an edit cannot strand it) or PROSE (not a
cardinality of anything this codebase has, with the reason written down). A
new count that is none of those fails the test and is named with its context.
`ASSIGNMENT_CLASSES` was added so the README's "five values" has something to
be derived FROM, and the join log line is built from it.

#### Six claims the fifth reading corrected

**Five of the seven closed vocabularies were enforced nowhere but in the eleven
configs the tests walk.** `_check()` guarded `found.kind` and `expect.kind` and
wrote down the reason — a typo here should be a loud failure and not a value a
consumer's switch falls through on — and that reason never applied to only two
of them. Fed a near-miss straight into `_check()`, `status` accepted `"faill"`,
`remedy` accepted `"confgi"`, `depth` accepted `"kindd"`, `fails_reason`
accepted `"setting_ignorred"`, `blocks_commands` accepted `["repport"]`, and
`blocks` — whose names join to `describe --json`'s `stage_names` — accepted
`["run"]`, which is a COMMAND. `status` is the expensive one: `exit_status` is
`1 if fails else 0` counted over `status == "fail"`, so a row that MEANS to
fail and misspells it is a silent zero, in a branch no config in the suite
reaches, with nothing anywhere objecting. Every one of them is guarded now, at
the place the rows are made and before anything downstream reasons about the
value; `blocks` and `blocks_commands` are checked against EACH OTHER, so a
command in the stage list is refused with the sentence that ends the mistake
rather than with "unknown stage"; and `degrades`, which is a stage list too,
got the same guard. A table-driven test drives a near-miss into every
vocabulary, and a second test asserts that every closed `DOCTOR_*` tuple has a
row in that table — because the defect was an omission, and an omission is
caught by an enumeration or not at all.

**Three published sentences said a stage "dies" on a path where `run` HANGS.**
The differential sweep's only `doctor`-against-`run` disagreements were seven
FIFO cells, and in every one `doctor` was right and fast while `run` never
returned — and three sentences still described that as a death.
`precomputed_emapper:<i>` said "it dies reading this"; driven, `prepare_emapper()`
on a FIFO was still scanning after six seconds and returned only when a writer
arrived. `tmt:<plex>:level_file` said the reader "opens the path for the header
and cannot read a table out of it", which is a sentence about a call that
RETURNS; driven, `header_columns()` on a FIFO was still alive after five, while
its SIBLING row twenty lines below in the same loop already said "or, on a FIFO
with no writer, never returns at all". `TUTORIAL.md` told an operator that a
directory "dies as it opens the path" and that "a FIFO or a socket behaves the
same way". The verb is DERIVED now — `_open_outcome()` reads it off
`found.kind`, exactly as `_stat_depth()` derives the depth and
`_present_kind_phrase()` derives what is there — because a verb is a claim
about a path state, and a hand-written claim about a path state is true of the
state its author had in mind. Two more were found with it: the
`emapper_precomputed` failure row was the one fail row in the document with no
`_present_kind_phrase()` on it, so a FIFO there was not even named as one; and
the `quant_table` row appended the FIFO paragraph to its `fragpipe_tmt` arm,
where `tmt_plex_dirs()` only globs and nothing is ever opened — a hang claimed
where there is none, which is the same defect pointing the other way.

**`run` hangs on a FIFO, that is an engine defect, and it stays.** The decision
is recorded where the primitive is rather than only here, because a CHANGELOG
entry ages out of the file a maintainer is reading. `regular_readable()` would
turn all seven hangs into refusals with a message, and it is deliberately not
called from a stage, for two MEASURED reasons. It cannot tell a FIFO that will
never be written from one being written right now — `O_NONBLOCK` returns at
once either way and `fstat()` says `S_ISFIFO` either way — so a gate would
refuse `mkfifo p; zcat big.faa.gz > p &`, which is how a disk-constrained
cluster feeds this tool and which works today. And the probe is not a
read-only observation: it opens the read end and closes it, and the producer
gets `EPIPE` for it — driven, the writer died with `BrokenPipeError` while
`regular_readable()` was deciding. A gate in front of a stage's open would
destroy the stream it was meant to protect. `doctor` may never block, because a
command whose job is to answer before the run is worthless if it hangs; `run`
may block, because blocking on a pipe is what reading a pipe is. Two tests pin
the capability, so a future round that decides otherwise has to retire them on
purpose.

**"a consumer reading the eight fields it knew about is unaffected by a NINTH
appearing", over a per-stage dict that emits ten.** The fourth consecutive
round to ship a wrong hand-written count, and this one sat above BOTH `requires`
and `degraded_by`, which arrived together in this change set.
`tests/test_config.py` used the identical construction for `cost` correctly and
`DESCRIBE_PER_STAGE` in that same file lists all ten, so two files disagreed
and the one in the shipped source was the wrong one. The comment states the
rule without a hand-count now, and the one number it does carry is recomputed
from the dict literal itself.

**The count scanner could not see a number written in digits, did not read
`tests/`, and read four constants where the emitted document is wider.** Three
structural reasons the sweep that was built to end this class did not catch the
fourth instance, all fixed here rather than in the sentence. `_COUNT_RE`
alternated over number WORDS only, so every count written 200,000 / 200 / 9 /
12 / 14 / 768 was invisible — 768 had been caught by a bespoke test written for
that one string. `tests/` was not one of the surfaces, and two of the three
stale counts this round found live there, one of them in the docstring of the
test that counts this file's own entries. And `DMND_READS` publishes "up to
200,000 records" and "up to 200 deflines" as caveat text on real
`db:diamond:<tag>:usable` rows — IN the emitted document, and the sentences
that qualify its own "doctor does not parse" — hand-written, one function away
from the literals that decided them; both are constants now, the caveats are
built from them, and the scan reads them. The noun alternation was missing
`rows`, `keys`, `paths`, `fields`, `files`, `surfaces`, `records`, `deflines`
and `reasons`, two of which were the noun of a live wrong count while they were
missing. Reading digits without drowning in them took one rule: a digit that is
part of a larger token is not a count, so `2.4 MB` is not "4 MB", `mode-000` is
not "000" and `cost-3 stages` is not "3 stages", and a unit between the number
and the noun disqualifies the pair, so "a 4 MB state file" is not four states.

**Tests that passed without the fix they name.** A test that cannot fail is
worse than no test, because the next reader believes the thing is pinned. The
null-kind test asserted only that every row carrying a null `found.kind` is
named by `DOCTOR_NULL_KIND_ROWS` — so putting "the DIAMOND usability rows" back
into that constant, the exact false entry that stood for three rounds, passed
it unchanged. It now asserts the other direction too: a glob may only name rows
that really carry a null, and it has to name some. `_fifo_clause()`
reaching "every row a FIFO can reach" was a claim a report made and nothing
checked; the rows are a list now and a parametrised test drives a FIFO into
each and holds the row against what the reader really does — with the mirror
test that a DIRECTORY at the same path is described differently, so the fix
cannot be one hard-coded verb. A duplicate
FIFO-streaming test added this round was deleted in favour of the one in
`tests/test_stages.py` that already drove it. And `NUMBER_WORDS` was defined
twice in `tests/test_docs.py`, the second shadowing the first, so three
assertions written against a table containing `"no": 0` were running against
one that did not have it.

#### Ten defects the sixth reading corrected, one of them the hang itself

**A FIFO with no writer no longer hangs `run`: it waits, says so, and then
dies naming the path.** Seven operator-supplied paths used to hang forever
when a FIFO sat at them — `proteins_faa`, `quant_table`, `manifest`, `gff`,
`emapper_precomputed`, and a TMT plex's `ion.tsv` and its annotation — all
seven measured against the real `run`, which neither finished nor failed,
printed nothing, and left no exit status to time out against. A previous round
declined to fix it and wrote the reasoning into `regular_readable()`'s
docstring, and BOTH of its reasons were correct and both were driven:
`regular_readable()` cannot tell a FIFO nobody will write from one being
written right now, so gating a stage on it would refuse `mkfifo p; zcat
big.faa.gz > p &`, which is how a disk-constrained cluster feeds this tool and
which works today; and worse, the probe is not a read-only observation — it
opens the read end and closes it, and the writer on the other side gets
`EPIPE`, measured, with the writer thread dying while the probe was still
deciding.

The fix keeps both of those and pays neither, by removing the thing that
causes them: it OPENS ONCE AND KEEPS THE DESCRIPTOR. `opener()` — the single
choke point every parser in this file reads through — now opens with
`os.open(O_RDONLY | O_NONBLOCK)`, which on a FIFO returns at once instead of
waiting AND releases a writer that is blocked in its own `open()`, and then
`fstat()`s the DESCRIPTOR rather than the path, so nothing can be swapped
underneath the decision. A regular file has the flag cleared and is read, one
open, exactly as before. A FIFO gets a log line naming the path and saying the
read is waiting — that alone removes the silence, which was the worst part of
the old behaviour — and then a bounded wait for the descriptor to become
readable; if a writer appears, `O_NONBLOCK` comes off and the read goes
through the SAME descriptor, so the live-writer workflow is byte for byte what
it was. If the wait expires the stage dies naming the path, the kind, and the
setting. A socket, a device node or a directory is refused at once, naming
what is really there. Clearing `O_NONBLOCK` before the read is not
housekeeping: a non-blocking descriptor raises `EAGAIN` mid-stream the moment
a writer pauses, which would have turned a working pipe into an intermittent
failure — worse than the hang, because a hang is at least reproducible. The
gzip branch reads the descriptor too, rather than re-opening by name;
`emapper_precomputed` is routinely a `.gz`, so that branch is the one the
commonest piped input goes through.

**`header_columns()` opened its path with a bare `open()`, one function away
from the gate.** It is what `refuse_isobaric_matrix()` calls before any quant
reader runs, so a FIFO at a TMT plex's `ion.tsv` — and at `quant_table` — hung
THERE rather than in the reader the row named. One reader outside the choke
point is one path still hanging, which is the whole argument for there being a
choke point. It reads through `opener()` now; `newline=""` went with the
change and nothing moved, because the only use of the line is a
`rstrip("\r\n")` that strips a CRLF either way.

**`fifo_wait_s`, six hours, and the number is an argument rather than a
taste.** The two costs are not symmetric. Waiting too long costs only the tail
of a mistake already made, and it now costs it visibly, since the wait is
announced before it starts. Waiting too briefly costs a workflow that works,
and the writer that has to survive is not the shell one-liner — it is a
producer that is itself a queued job on a shared cluster, where waits are
measured in hours, so anything in minutes would refuse a correct setup while
claiming to protect it. The upper bound is what makes it six and not sixty: a
run holds an exclusive lock on its results directory while it waits, and a job
started at the end of a working day should have failed with a message by the
next morning rather than still be sitting on the open. `0` refuses a FIFO
outright, for an operator who pipes nothing in. `doctor` is unchanged and
still never waits at all: it reads every path through `_deep_readable()`,
which is the asymmetry the two commands are supposed to have.

**`section` was the eighth closed vocabulary and `_check()` did not guard
it.** A typo took the TEXT command down with a bare `KeyError` out of
`print_doctor()`'s `titles[seen]` — no document at all, which is the most
expensive failure this command has — while under `--json` the row was dropped
from `sections` in silence and left in `checks`. The enumeration tripwire in
the suite could not see the omission either, because it collects the `DOCTOR_*`
constants that are TUPLES and `DOCTOR_SECTIONS` was a list. It is a tuple now,
`_check()` has the guard arm, the near-miss table has its row, and the
tripwire works for it.

**A malformed `found` or `expect` was accepted or crashed, rather than
refused.** `found` was read through `.get("kind")`, so a dict built by hand
without one answered `None` — which is a LEGAL value meaning "this build did
not compute it" — and the row shipped a `found` that every consumer switches
on and that says nothing. `expect` was indexed directly, so the same mistake
was a bare `KeyError` raised from inside the function whose entire job is to
make a malformed row impossible, and it cost the whole document. Both are
named now, with the sentence saying which builder to use.

**README.md said the suite ships "thirteen test modules"; it ships sixteen.**
Published, and made worse by this change set rather than by drift alone.

**The state sweep's own counts were one state and twelve cells stale.** The
entry above described it as six states and thirty cells while the sweep had
already outgrown both numbers, and it now covers `emapper_precomputed` as
well — which is the second half of the same correction: another entry claimed
the sweep walked EVERY path `doctor` can be pointed at while that one key was
outside it.

**`metaannot.py` pointed a reader at a `_diamond_source_fasta` this file does
not have.** The function is `diamond_source_fasta()`; the underscored spelling
appeared exactly once in the repository, in the docstring that named it, and
it is quoted here without its parentheses because a test refuses a changelog
that names a function nothing defines — which is the same rule, one surface
over.

**Two docstrings disagreed about which paths had been measured against
`run`.** One said flatly that seven paths hang `run` indefinitely while
`_open_outcome()`'s hedged about which of them had been measured against the
command and which only against the reader. All seven have now been measured
against `run`; the two say the same, stronger thing.

**The count scanner's noun list was a WHITELIST, which is the defect class
rather than the nouns missing from it.** It could not see modules,
vocabularies, tests, callers, assertions, branches, sentences, helpers,
functions, guards or arms — and the two wrong counts above are exactly what a
whitelist cannot catch. The rule is negative space now: a number followed by
something that reads as a countable noun must be classified — derived,
measured or exempt — or the scan fails and names the sentence. It surfaced a
pile of counts that had never been visible, each of which is now classified.
`_emitted_prose()` is built from a real `doctor --json` document as well,
rather than from a hand-listed set of constants, so a number in a sentence the
document publishes is in the scan by construction rather than by somebody
remembering to add the constant it came from.

#### Twelve defects the seventh reading corrected, and the premise under them all

**`mkfifo p; zcat big.faa.gz > p &` does not work at `proteins_faa` and never
did — the premise two rounds of this change set reasoned from was false.** A
verifier drove that workflow end to end through the real CLI at each of the
operator-supplied inputs, with live writers. It works at `emapper_precomputed`,
plain and gzipped. At `proteins_faa` the run FAILS at `integrate`; at
`quant_table` it HANGS FOREVER holding the results lock while the writer takes
`BrokenPipeError`; at `manifest` it fails after burning the whole wait. The
reason is one sentence and nobody checked it: an ordinary run opens each of
those inputs MORE THAN ONCE, and a FIFO can be drained exactly once. Traced
with a `sitecustomize` shim over `builtins.open` and `os.open`, an ordinary run
opened `proteins_faa` four times, `quant_table` four, `manifest` three and
`emapper_precomputed` twice. So the previous round waited six hours in order to
then fail, and quoted as advice, in the failure, a workflow that cannot work at
that path — which is worse than the hang it replaced, because it is
confidently wrong.

**The read count is DERIVED and a test recomputes it from a driven run.**
`INPUT_READ_SITES` names, per configured input, every site a run opens it at
and the condition under which that read happens; `set_read_plan()` evaluates it
against the config in hand, because the count is a property of the RUN — a
project with `run.unipept` on reads its quant table once more than one without,
so a pipe that works in the first config does not work in the second. It is not
trusted: a test drives the real CLI under an open-tracing shim and fails
naming the new reader when one appears. (It drove two configs whose answers
differ, which was not enough; see the eighth reading below, where it became
`test_the_read_plan_and_the_choke_point_hold_across_the_config_space`.) A hand-written list of readers is
exactly how this went wrong the first time.

**A FIFO at a path the run opens more than once is refused at the FIRST open,
naming the count and the reads.** Not after `fifo_wait_s`, because waiting
cannot change the arithmetic: the first reader drains the pipe and the next
finds an empty one. The message lists the reads, says that waiting longer,
starting the writer earlier and making it faster all change nothing, and names
the inputs that DO take a pipe in this config — read off the plan rather than
asserted, so it cannot go stale against a config it was not written for.

**The bare opens the previous round's report said did not exist.** That report
stated `opener()` was the only reader of the quant table. `read_delim_table()`
opened it again for its delimiter sniff and pandas opened it a third time, and
that bare open at `metaannot.py:7092` is where `quant_table` ACTUALLY hung —
the `SIGABRT` stack is `read_delim_table -> read_feature_table -> stage_join`,
and the diff never touched it. `_fasta_lengths()`, `_count_fasta()` and a bare
`pd.read_csv(qpath)` in the protein-level branch of `stage_join` were three
more. All four go through `opener()` now, and a test drives a run and fails
on any open of a configured input that does not — the CLASS, rather than the
instances that were found by hand. (It drove one config, in which the reader
it was written for is never reached; it is part of
`test_the_read_plan_and_the_choke_point_hold_across_the_config_space` below.)

**`quant_table` and `manifest` are single-read now, so a pipe works at them.**
`read_named_table()` reads the header and the body from ONE open, so the
recogniser that must run before pandas does no longer costs a second open;
`read_delim_table()` puts its sniffed line back with `_HeadRestored()` instead
of letting pandas re-open the path, which keeps every byte pandas would have
seen in the order it would have seen it — asserted against a table with a
duplicate column name and an embedded quoted delimiter, because `names=` would
not have; `stage_join` reads one header for both of its refusals; and
`read_feature_table` reads the manifest once and uses it for both the column
mapping and the design. `proteins_faa` stays multi-read and is refused: the
emapper stage and the integrate stage read it in different stages at different
times, and one of them can be adopted from cache without the other running, so
there is nothing clean to collapse.

**`header_columns()` probed and closed in front of `read_delim_table()`'s own
open, in the same function.** That is precisely the probe-then-close pattern
the single-open design exists to avoid, reintroduced by the fix that was meant
to remove it: driven, a re-opening writer took three `EPIPE`s. The
NO-EPIPE claim was falsified by driving. There is one open on that path now.

**The wait's message was FALSE in a reachable case.** A writer attached at 0s
whose first byte lands at 5s of a 3s wait got "no writer appeared ... nothing
has opened the other end" — both halves false — and then died of `EPIPE` at
5.01s because the reader had given up. A non-blocking one-byte read tells the
two apart for nothing: a pipe with no writer answers end-of-file, one with a
silent writer answers `EAGAIN`. The refusal now says which it is, and offers
the remedy that goes with it — start a writer, or raise `fifo_wait_s`.

**A stream that ended early was not an error at all, which is the worst
outcome here.** A writer that wrote half a FASTA and died yielded
`[("p1", "MKV"), ("p2", "MK")]` with nothing to say the input was a fragment: a
short read looks exactly like a complete file to every reader above the stream,
so the number at the end of the run is wrong and looks right. The gzip branch
was already correct — gzip is FRAMED, so a missing end-of-stream marker raises
`EOFError`, which is now translated into a `StageError` that names the path. A
plain stream has no frame, so what it can see is that the writer closed without
writing a byte, or that the stream stopped in the middle of a line. A writer
killed exactly on a line boundary is a case no amount of looking at the bytes
can catch, and that limit is written down rather than left to be found: pipe
`.gz` when you have the choice.

**`fifo_wait_s` bounded only the first byte.** An idle writer holding the write
end open — an `O_RDWR` keeper, a producer blocked on its own input — put the
run straight back into the hang the setting exists to prevent, one read further
in, and back onto the results lock it holds while it waits. The bound is an
IDLE timeout now: every read waits at most `fifo_wait_s` for the NEXT byte. A
writer that is streaming never comes near it, which is pinned from the writer's
side as well as the reader's.

**`metaannot.py` published a Windows fallback that does not exist.** "On
Windows there is neither ... so `opener()` falls back there to exactly the
plain `open()` it has always used" — there is no such branch; `opener()` calls
`_open_for_read()` unconditionally, demonstrated with `fcntl = select = None`
and a count: one `os.open`, zero builtin `open`. It matters past the false
sentence, because CPython's builtin `open()` adds `O_BINARY` on Windows and a
bare `os.open()` does not, so the gzip branch — `os.fdopen(fd, "rb")`, and
`emapper_precomputed` is routinely a `.gz` — would have read compressed bytes
through a CRT text-mode descriptor: CRLF translation and `0x1A` treated as end
of file. The author reached for `getattr(os, "O_NONBLOCK", 0)` and stopped one
flag short. `_OPEN_FLAGS` carries both, and a test reaches it by giving the
platform an `os.O_BINARY` to find.

**The derived verb was published for sockets and device nodes, where it is
false.** `_raises_promptly()` was `f["kind"] != "other"`, and `other` is one
bucket for a FIFO, a UNIX socket AND a device node — so a socket at
`proteins_faa` published "dies, but not at once ... A FIFO in particular is not
refused on sight ... waits `fifo_wait_s`", while a socket cannot be opened as
a file at all and `os.open()` itself fails on it, which is what the stage-side
test asserts. The two halves
of one change set contradicted each other in the emitted document, inside the
derivation that was built to remove exactly that defect class. `found` carries
an `other_kind` now — a new KEY, not a new `kind`, so no `DOCTOR_VERSION` bump
— every sentence is keyed on it, and `PATH_STATES` has a `socket` and a
`device_node` state, whose absence is why nothing caught this.

**A test that could wedge the whole suite, and a DERIVED count that derived
nothing.** `test_the_wait_is_the_configured_one_and_zero_refuses_at_once`
called `list(read_fasta(fifo))` on the main thread with no thread, no join and
no timeout, twenty lines below a sibling whose docstring carries exactly that
discipline; against a reverted `opener()` it hung until `SIGKILL`. Every FIFO
test goes through one helper with a join now, so a regression FAILS rather than
wedges. Separately, `COUNT_PROSE["three states"]` was a DERIVED row whose
recompute was the literal `3` — a rule with nothing under it, one level up from
the staleness check that exists to prevent them. It derives the states from the
sentence that enumerates them, and `test_no_derived_count_recomputes_a_literal`
now refuses any DERIVED entry whose recompute reads nothing at all.

#### Seventeen defects the eighth reading corrected, and the class under all of them

The class is this: **the read plan was a hand-written table, and a
hand-written table states a fact once per place it happens to be true.** Every
defect in this group is that shape or a test that could not see it.

**The plan under-counted `manifest`, and `doctor --json` published the wrong
promise before the run.** `INPUT_READ_SITES["manifest"]` named one site,
`stage_join -> read_manifest`, and missed that `read_feature_table()` opens
the manifest itself — so it is reached from `stage_join` AND from
`peptide_features()` in each peptide stage. Driven with `run.unipept` and
`run.taxonomy` on and a manifest on a FIFO with a live writer, the run logged
"this run reads it exactly once, so it is read as a live stream" and then died
at the second open with "this run has already read it once". The `_FIFO_TAKEN`
backstop worked exactly as designed — a message, not a hang — and everything
above it was wrong. Worse, `doctor` had already published "A FIFO here is NOT
refused on sight, because this run reads this input exactly once … a supported
way to feed this tool at a single-read input" for the same config: the command
an operator reads before committing days of compute blessed a workflow the run
refuses. In the toy project the taxonomy stage won the race and `join` died at
0.3s; in a real run `join` is late, so the cost is the whole pipeline.

**The plan over-counted, and the refusal quoted reads the run will not
make.** With `unipept.result` set, `stage_unipept` ingests the export
and *returns* before `peptide_features()`, so `quant_table` was planned 2 and
read 1 (3 and 2 with `run.taxonomy` on). With `run.eggnog` off, nothing opens
`emapper_precomputed` and the plan said 1. And on the MSstats formats nothing
opens the manifest at all while the plan said 1. Over-counting refuses
a pipe that would have worked, with a message naming reads that never happen.

**The sites are derived from the same predicates `doctor` reasons with.**
`quant_consumers()`, `_peptide_readers()` and `_manifest_readers()` already
knew every one of those facts, because `doctor`'s rows are built on them —
which is why one surface was right and the other wrong. The plan calls them
now, so a row and a refusal cannot answer "who opens this input" differently.

**`db.ncbi_taxonomy`'s `.dmp` files hung `run` with the results lock held.**
Four bare `open()` calls inside `NCBITaxonomy`, on paths an operator types
into a config exactly as they type `quant_table` — 45 seconds, no output, no
exit status, lock not released. This is the original failure mode of the whole
change set, intact, in the one class of path nobody had swept: an `input:` row
made a path feel operator-supplied and a `db:` row made it feel like
infrastructure, and an `open()` cannot tell the difference. They go through
`opener()` now and they are in the plan, which also makes the count real: a
config with `run.taxonomy` AND `taxon_rank` builds two `NCBITaxonomy` objects
and reads every `.dmp` twice, so a pipe cannot work there and the refusal says
so; with one of them, it can.

**`unipept.result` was read through a bare `pd.read_csv(sep=None,
engine="python")`.** A sixth reader of an operator-supplied path outside the
choke point, found by tracing every open a real run makes of every path its
config names rather than by reading the list of readers somebody had written.
It goes through `read_delim_table()` now, which is one open, one sniff and the
C parser instead of the slow one.

**`subset --quant` hung on a FIFO, and `subset` had never been swept at all.**
`cmd_subset()`'s protein-level branch was `pd.read_csv(args.quant, …)` — a
path a user types on a command line, which is as operator-supplied as a path
gets — and it blocked with no exit status. Every sweep in this suite had
driven `run` and nothing else.

**`doctor` re-opened the quant table by name after proving it.**
`_manifest_checks()` probed with `regular_readable()` and then let pandas open
the path again: a stat and an open with instructions in between, which is the
race every other check in that file was rewritten to remove. It also split the
header with pandas' own sniffer while the run splits it with `_sniff_sep()`,
so `doctor` could pre-check a mapping against a different set of columns from
the one the run would see. `_open_regular_text()` is `doctor`'s opener now —
`_open_for_read()` with the waiting arm removed — and the header goes through
the shared `_header_cells()`.

**The plan test drove two configs and the choke-point test drove one.** A
verifier reverted `read_delim_table()`'s `opener(path)` back to a bare
`open()` — the very reader the last report named as "where `quant_table`
actually hung" — and the whole suite came back byte-identical. It is not dead
code: `read_delim_table` is the second reader of `quant_table` under
`quant_format: msstats_protein` and of `analysis.metadata`. Both tests are one
sweep over the CONFIG SPACE now — every quant format crossed with the run
flags that move the answer — and the paths it watches are taken from the
CONFIG rather than from the table under test, so an input nobody listed cannot
hide. The same revert now fails seven of its points.

**An input with no plan that is read twice is a failure, not a gap.** A path
with no entry is treated as single-read, which is the hang with the refusal
switched off — so the sweep asserts it, and that is the assertion
`db.ncbi_taxonomy` would have failed for as long as it has existed.

**The `_OPEN_FLAGS` test could only agree with itself.** It monkeypatched the
constant with its own recomputation and then asserted that its own expression
had the `O_BINARY` bit: driven, deleting `| getattr(os, "O_BINARY", 0)` from
`metaannot.py` left both C1 tests passing. Same shape as the C4 defect this
change set had just fixed, one level up, in the test written for it. It stubs
`os.O_BINARY` BEFORE importing the module and reads `_OPEN_FLAGS` off it,
which separates the trees: `0x8004` with the term, `0x4` without.

**A character device published a sentence that contradicted its own row.**
`_open_for_read()` falls through to `die()` for a device node and `die()`
raises a `StageError`, which `peptide_features()` catches — measured, a
taxonomy run over a char-device manifest logged "falling back to the
peptide-only reader" and returned 30 rows. The `input:manifest` row keyed its
sentence on `_raises_promptly()`, true of a character device because a
character device answers at once, and so published "it dies with an OSError …
every stage that opens a quant table on this format goes with it" beside its
own `blocks: ["join"]`. `_dies_with_a_stage_error()` is the predicate now, the
sentence and `blocks` are both computed from it, and it fixes an under-claim
in the other direction too: a socket and a block device share the `other`
bucket with the FIFO and end in an `OSError`, which the old
`empty or kind == "other"` test read as survivable.

**Two FIFO tests still wedged rather than failing.**
`test_header_columns_on_a_fifo_waits_and_then_dies` and
`test_prepare_emapper_on_a_fifo_waits_and_then_dies` each called the reader on
the MAIN thread inside `pytest.raises` — the assertion that cannot fail, since
against an unbounded opener the call does not raise, it blocks. Found the only
way this class can be found: by reverting the bound and running the suite
under a watchdog that tells a wedge from a slow test by where the main thread
is standing.

**"A socket fails ENXIO at `os.open()`" is wrong on the platform it was
measured on.** It is `EOPNOTSUPP`. The errno is the OS's to choose and is
published nowhere now; what is asserted, in the code and in the prose, is that
an `OSError` comes out of the open, which is the half that decides who dies.

**`cmd_doctor()` computed the run's read plan and not the run's wait.**
`set_read_plan()` was called and `set_fifo_wait()` was not, so `doctor`
reasoned about a wait the run would not use. Latent — the FIFO clause reaches
for `DEFAULT_CONFIG['fifo_wait_s']` explicitly — and fixed anyway, because
"doctor reasons about the run's settings" has to be true of all of them or it
is not a rule.

**Two `doctor` rows quoted a pandas message this program can no longer
produce.** They promised "Could not determine delimiter" for a zero-byte
`unipept.result`, which was pandas' python-engine sniffer talking; with the
reader through `opener()` it is "No columns to parse from file". The reader
refuses in this file's own words now, the sentence lives in one constant both
rows quote, and a test drives the reader and holds the row against what came
out.

**The README had the two before-counts swapped.** `quant_table` was read three
times and `manifest` twice, not the other way round.

**`DESCRIBE_PER_STAGE` grew `requires` and `degraded_by` without an entry
here.** The change is right and was explained in a comment beside it — the two
fields are on `STAGES`, `describe --json` emits them, and
`test_describe_emits_every_field_a_stage_dict_carries` in the console contract
is what forces the decision — but a tripwire that asks for a deliberate
decision is not satisfied by a comment the changelog does not carry.
`DESCRIBE_VERSION` does not move: nothing was removed or redefined.

#### Seven defects the ninth reading corrected, and the class under all of them

The class is a **claim about a read that nothing measured**: a read the plan
did not know about, a read no sweep ever drove, or a refusal this file assumed
somebody else's program would make. The eighth reading replaced the
hand-written table with functions; this one drove the commands and the configs
the functions had still never been asked about.

**The read plan had no entry of any kind for a `fragpipe_tmt` tree, and that
is `db.ncbi_taxonomy`'s defect one input over.** `quant_table` names a
DIRECTORY on that format, and a directory is not something any reader opens —
the per-plex `ion.tsv`/`peptide.tsv`, each plex's annotation file and, when
`tmt.min_purity` reads one, its `psm.tsv` are. `_quant_table_sites()` returned
`()` for the format on the grounds that a pipe at the setting itself never
reaches an open, which is true and beside the point. Driven, a two-plex tree
with `run.taxonomy` and `run.join` on opened every file under it twice with no
plan entry to say so — an unplanned path is treated as single-read, so this
was the hang with the refusal that exists to prevent it switched off. The
paths now expand the way `_taxdump_paths()` expands its directory, and the
files only the FULL reader opens are their own input, because the peptide-only
reader it falls back to never opens an annotation and a site claiming it did
would be a false read in a refusal message.

**The plan under-counted the quant table wherever the peptide-only fallback
can run.** `peptide_features()` calls the full reader, CATCHES its refusal and
re-reads the same table with `read_feature_peptides()` — a second open of one
operator-supplied path — and the plan said one. Driven with a manifest whose
one run matches no column, `run.taxonomy` on and `run.join` off: two opens,
plan said one, and `doctor --json` published the single-read promise for it.
The site is planned now and MARKED, because it is the one read in the plan
that turns on the FILE rather than on the config: `planned_reads()` counts it,
which is the FIFO question — a pipe is drained by its first reader whether or
not the second open was conditional — and `certain_reads()` does not, which is
what a measurement of a completed run is held against.

**`auto_contrasts()` read `analysis.metadata` through a bare `pd.read_csv()`,
and `metaannot report` used to HANG on it with no exit status.** On a TMT project the
report takes its contrasts from the operator's metadata, and that path was
read twice: once by `header_columns()` in `_tmt_report_design()`, which is
gated and which printed "this is a FIFO, and this run reads it exactly once,
so it is read as a live stream", and then here. The second open found a
drained pipe with no writer on it and blocked there — no further output,
nothing to time out against, and `fifo_wait_s: 0`, which exists to refuse a
FIFO outright, made no difference at all, because it bounds `opener()` and a
bare pandas open has no wait to bound. The promise was printed one line before
the silence. The read goes through `opener()`; the second open is refused by
name instead, which is what `_FIFO_TAKEN` is for. A command that can hang is
exactly as bad as a stage that can hang, and `report` holds no results lock
only by accident.

**`report` and `object` had never been swept, which is why the reader above
survived.** The sweep that found a bare `pd.read_csv()` in `cmd_subset` drove
`run`, `doctor` and `subset`, and a command nobody drives is a command whose
readers nobody sees. It drives `describe`, `report` and `object` too now.
`describe` and `object` open no operator-supplied path at all — the measured
answer for them, asserted as one, so the day either grows a reader the sweep
says so instead of passing on an empty trace.

**`doctor`'s manifest promise was guarded by nothing in
`tests/test_doctor_json.py`.** The correction this whole round was called for
— no longer publishing "this run reads this input exactly once" for a manifest
a taxonomy run reads three times — was measured by putting `_manifest_sites()`
back the way it was: every other test in that file stayed green
while the document published the false sentence again, and the only thing that
fell was a plan comparison in `tests/test_stages.py`, which is a different
file catching it as a side effect. A promise is published in the document and
has to be pinned against the document. The new test takes its count from
`_manifest_readers()` and not from the plan under test, because a count
recomputed from the thing under test moves with it.

**`_manifest_checks()`'s missing-manifest branch asserted who it kills instead
of deriving it.** The wrong-kind branch beside it computes that from
`_dies_with_a_stage_error()` — a `StageError` is a refusal `peptide_features()`
re-reads around, an `OSError` is caught by nothing — and this one wrote the
answer out by hand, which is how the wrong-kind branch itself once came to
contradict its own `blocks` field. The answer does not change (`os.open()`
raises `FileNotFoundError` on an absent path, and a driven run with the
taxonomy stage on dies with "[Errno 2] No such file or directory" and no
fallback); it is now the same computation rather than a second statement of
it.

**`doctor`'s `contigs_fna` row and the read plan's own comment claimed what a
third-party ORF finder does.** The row said a path that is not a file is
handed to `smorf`/`macrel` "which cannot read an assembly out of it", and the
plan said a FIFO there "is refused by smorf or macrel, not by anything here".
Nothing in this codebase establishes either, and neither is this codebase's to
promise. The alternative — fstat the descriptor before `run_cmd` launches the
tool — was considered and rejected, for the two facts already measured in
`regular_readable()`'s docstring: every other check in this file opens ONCE
and keeps the descriptor, which is what makes it both safe and free, and a
path handed to a subprocess BY NAME cannot be checked that way, because the
tool does its own open. Probing first is open-close-reopen, which sends a live
writer `EPIPE` while we are still deciding, and it cannot tell a pipe with a
writer from one without — so it would break, and then refuse, `mkfifo p; zcat
contigs.fna.gz > p &` in order to catch a misconfiguration. Both sentences say
what is known instead: the path is handed over as a filename, this process
never opens it, and what happens next is that tool's answer and not this
one's.

### Changed

**The state file is no longer written from one process's snapshot.** `#23`.
`save_state(path, state)` took a whole in-memory dict and rewrote the document
with it, and `_STATELOCK` beside it is a `threading.Lock`, which serialises the
threads of one process and nothing between processes. Two runs sharing a
results directory therefore each rewrote the file from their own view and the
last writer won, silently. A reviewer reproduced the consequence: run B
completed `pfam`, `dbcan`, `diamond` and `cluster`, run A wrote its own older
one-stage dict wholesale, and run C then printed `adopting output this run did
not produce` for the stages that had gone — the warning that says outright it
cannot tell a finished file from an interrupted one.

Every write goes through `update_state(path, state, keys, claim, drop)` now and
names the keys it changed — which at every site that changes the document is
exactly one key, and always was. It
re-reads the document immediately before writing, merges, writes, and reads
back to check the file still says what it wrote, redoing the merge against the
current document when it does not. `save_state()` keeps its name, its
signature and its whole-document semantics as the write step, and is still what
a failure-injection test patches; a state file that is MISSING and one whose
bytes will not parse are handled like any other now — the write creates a
document holding the keys it names and nothing else, and rebuilds no part of
what was there. The pre-write read is
deliberately not `load_state()`: that one WARNs "every stage will be
recomputed" and returns an empty `_State`, so per write it would say that on
every write and merge into `{}` — deleting every other record because the
filesystem hiccuped once.

None of this is exclusion and no sentence in it says otherwise. `os.replace`
was already atomic, so torn files were never the defect; what the merge does is
shrink the lost-update window from the length of a run to the gap between one
read and one rename. What is left inside that gap is not noticed: the read-back
compares the file against the payload just renamed, so it catches a writer that
lands AFTER the rename and redoes the merge, while one that landed BEFORE it is
simply not in what was merged — this run's complete document goes over it, the
read-back matches, and nothing is logged. On NFS, attribute
caching caps even that, which the comment above `update_state()` says in so
many words.

**A run whose state file names a run it has never seen stops writing.** The
merge read produces the succession check for free: `_run`'s identity is
compared against this run's own and that of the run it took the directory over
from, captured before `_run` was overwritten so that an ordinary resume — where
a foreign `_run` is the normal thing to find — keeps writing for the whole of
its life. A THIRD identity is a replacement, and the run stands down, latching
the same `superseded` flag `ResultsLock.is_still_ours()` raises. Any doubt at
all is "still ours": no `_run`, no `run_id`, a record that is not a dict, a
read that failed. The identity is `(run_id, host, pid, started)` rather than
`run_id` alone, because `run_id` is `%Y%m%dT%H%M%S-<pid>` with no host in it
and two nodes of a shared array can mint the same one.

**A superseded run's stage output is parked instead of renamed over the live
run's.** The variant that never touched the state file at all: a stage still
inside `st["fn"]` when its run is superseded runs to completion, and
`atomic_out()` renamed its result into place at the end, under the replacement's
valid signature if the replacement had already recorded that stage. The rename
now asks `_directory_still_ours()`, which is an in-memory flag plus at most one
read of the state file per `STATE_PROBE_S` — `atomic_out()` is called once per
item by several stages, and the console solved the same shape in the same
directory with `PART_RESCAN_S`. There is exactly one path to not renaming and
it needs a read that SUCCEEDED and returned a foreign, parsed identity; every
failure renames, so the trade the issue rejects — a transient error aborting a
stage that is legitimately finishing — is not taken. Declared outputs are
parked as `.superseded.<stem>.<run_id><ext>`, with the marker in FRONT of the
stem so a console looking for `.<stem>.` does not offer a dead run's leftovers
as the live run's work in progress. Per-item outputs are not parked at all —
one PDB per dark protein would be one undeletable file per protein — and keep
the `.part` convention a crash already leaves.

**What that is a bound on, stated where the claim is made.** Two halves, and
they are different kinds of thing. A stage that STARTS after the handover
cannot rename, full stop: `mark_running()` is a merged write and the succession
check refuses it before the stage begins. A stage ALREADY RUNNING when the
directory changed hands is only NOTICED, at whichever comes first of the next
heartbeat tick and the fallback probe — `min(heartbeat_s, STATE_PROBE_S)`, 30 s
as shipped — and inside that window nothing detects the handover at all: a
takeover that completes in under a second is caught by nothing and that stage's
output lands on the live run's. It is a clock, not exclusion, and the earlier
wording of `CLAUDE.md` rule 4 and of the TUTORIAL — a superseded run "does not
rename its outputs over them" — was false for exactly the in-flight stage `#23`
is about. Both now say which half is which. Outputs written WITHOUT
`atomic_out` are not covered at all, and three stages write some: `hhblits`'
per-query `<id>.hhr`, `esmfold`'s `plddt.tsv` and `esmfold_failed.tsv`, and the
`.done` sentinel of `diamond`, `hhblits` and `esmfold`. That is in the
NOT-DONE list rather than fixed here, because routing them through
`atomic_out` changes the temp names the `.part` sweep and the console match on.

**Stage records carry the `run_id` that wrote them.** It makes a document two
runs have both written into auditable rather than inferred from timestamps,
which is the chief cost of merging. `tests/test_console_contract.py` moves with
it, since the `ok` record's key set is pinned there exactly.

**A state write whose pre-write read FAILED writes nothing.** Merging into the
empty dict a failed read returns would write a document holding this run's keys
and nothing else, deleting every other stage's record because the filesystem
hiccuped once — the worst failure this design could have. So the write is
declined and said once in the log. It costs at most the one record that call
carried, which `decide()` then reads as an output with no record: adopt or
recompute, never corruption. The next write of the same key — the stage's own
next record, or the next heartbeat tick — reaches the file normally.

**A write no longer rebuilds a whole document, and `_run` is written only on
proof of ownership.** `#23` again, after three revisions of the same fix.
Where the pre-write read came back MISSING or UNPARSEABLE, `update_state()`
rebuilt the whole document from `state` — the writing process's in-memory
snapshot — on the reasoning that a document which is not there has nothing in
it to lose. It has: the records of whatever run holds the directory now. Three
revisions kept that write and tried to gate it, and every gate was broken in
turn.

`_STATE_WATCH["lost"]` is raised only by `_judge_ownership()`, which runs only
on a read that SUCCEEDED and PARSED, so on the exact failure that reaches the
branch it is guaranteed silent. Asking the lock latch as well was not enough
either: nothing in `finish()`, `mark_running()` or `record()` reads the lock,
so in the real ordering neither half is up, and the test that passed did so
only because its fixture called `RunRecord._save()` first. And no lock read can
close it, which is the finding that settled it — A parks; B takes the directory
with `--force-unlock-live`, records `dbcan` and `diamond`, and EXITS, removing
its own lock; the document is removed; A's tail, or A's heartbeat with no stage
of its own finished at all, rebuilds `{_run: A, pfam: A}`. `is_still_ours()`
answers `None` on a vacant path and the lock latch can never come up.
Reproduced at `heartbeat_s` 0, 1 and 30. **A vacant lock is exactly what a
FINISHED replacement leaves behind, so vacancy can never license rebuilding a
document.**

So the write is gone rather than guarded. A run writes only the keys it is
updating: into a readable document it merges them, and where the document is
missing or unparseable it creates one holding just those keys. Every other
record stays gone. That is the trade this codebase already takes everywhere —
losing a stage RECORD costs a recomputation, because the OUTPUT is still on
disk and only its provenance is gone, and `signature()` prefers
over-invalidating to under-invalidating for the same reason. Losing a live
run's records to a dead run's snapshot costs a results table that looks fine
and is not. A run whose document goes missing under it says so, once, naming
the path, and carries on recording.

`_run` needed the same treatment one level up, because it is not a record of
work but a claim about who owns the directory. It is written only on
`is_still_ours() is True` — this run's own lock, or one it could not read,
which is no evidence that anything changed hands. `False` was already refused.
**`None` — a VACANT path — is now refused too, and the first version of this
change only pretended to.** That version declined the write that would CREATE
the document and allowed one into a document already there, on the reasoning
that such a document "is the one this run has been writing all along, and its
`_run` is already this run's". It is not: a run's own STAGE record creates the
document one call earlier, holding a stage key and no `_run` at all, and `_run`
then merged into it. So the gate never fired for any run that had recorded a
stage — which is every real run. Driven with real processes: A parks inside
`pfam`, B takes the directory with `--force-unlock-live`, records its stages
and exits, the document is removed, A's `pfam` returns — and A exits `0` leaving
`{_run: A, pfam: A}` having latched nothing and warned about nothing.
Reproduced in unit form on a removed document and on a zeroed one, at
`heartbeat_s` 0, 1 and 30.

**What the real gate costs is paid by the run that triggers it and by no other,
and that was measured before it was chosen.** A run whose lock is vacant writes
no `_run` again: no `last_seen`, no `final_status`, no `finished`. Its record
stops where the last heartbeat left it, `final_status: "running"`, which a
console reads as a run that never ended — so the refusal says so in the log
rather than going quiet. The same documentation had claimed this cost twice
already while the code was not paying it: a solo run that lost BOTH its lock and
its document still exited `0` with `final_status: "ok"` and `finished` set,
because the stage record recreated the file and `_run` rode in behind it. No
ordinary run pays it now either: the results lock is released by an `atexit`
hook that runs AFTER `cmd_run`'s `stamp_run("ok"/"failed")` and after `main()`'s
`stamp_run("interrupted")`, and `cmd_run`'s SIGTERM handler releases the lock
and `os._exit()`s without stamping at all. Driven end to end, `is_still_ours()`
answered `True` at the final stamp of a clean run and of a Ctrl-C'd one.

`cmd_run`'s own first `_run` write now goes through `RunRecord._save()` rather
than calling `update_state()` with `RUN_KEY` directly. It is three lines after
`lock.__enter__()`, so the gate cannot answer anything but "ours" there and
nothing about the run changes — the point is that "`_run` is written only on
positive proof of ownership" has no exception in it and `RUN_KEY` has exactly
one writer. Two ways to write that key, one of them gated, is how the gate this
replaced came to be inoperative in the first place.

STAGE keys are deliberately not gated that way — a stage record claims nothing
about the directory and the document it creates holds nothing to corrupt, so
gating it would take the solo case from "loses the records it had already
written" to "records nothing ever again". `_save()`'s "the final stamp is not
skippable" reasoning is unchanged for the case it was written for: a run that
still holds its own lock, or one whose lock is unreadable.

**`--force`'s discard no longer deletes a key the write does not name.** The
last place in the file where a write touched a record it had no knowledge of,
and it needed no race to do it. `--force` pops the records of the selected
stages from ONE read at the start of the run, and every merged write then
carried those NAMES and popped whatever was under them — so a superseded
`--force` run writing `pfam` was watched deleting a live run's `dbcan` record,
silently, out of a document holding no `_run` to stand it down. The discard now
carries the RECORD it was decided about and deletes only while the document
still holds that record; where it does not, the record is left alone and the
reason is logged once. A document with no `_run` in it is the ordinary shape a
run recreates after its document is removed, and on a replacement running with
`heartbeat_s: 0` that shape lasts until that run's final stamp, so this was not
a window to be inside.

`CLAUDE.md` rule 4, README and the TUTORIAL published the old claim as a
guarantee — "whether or not it has noticed that it was superseded", "That part
does not depend on timing", "That is a guarantee, and no part of it is a
clock" — and the records half is not one. What IS timing-free is that a write
names its keys and rebuilds nothing. What is a CLOCK, and the window is a
read-modify-rename gap rather than the thirty seconds the outputs half runs on,
is that a record another process writes between this run's pre-write read and
its `os.replace` is dropped rather than merged, and the read-back catches only
a writer that lands after the rename — `update_state()`'s own docstring says
this and the published sentences contradicted it. The mechanism is demonstrated
by injecting a write into that gap; the race was not won at shipped speeds, so
its width is unmeasured rather than small. All of them now say which parts are
guarantees, which are clocks and which are not covered, and the TUTORIAL's
troubleshooting row no longer tells an operator that a superseded run leaves
the live run's outputs untouched.

**What else the same reading turned up.** `_watch_state()` took an `interval`
and never read it, so it no longer takes one. The comment on `heartbeat_s` said
"nothing BRANCHES on this number", which is true of the value copied into
`_run` and false of the config key it sits on: `watch()` returns without
starting the thread when it is 0, `_beat()` waits on it, and it is half of
`min(heartbeat_s, STATE_PROBE_S)`. Correcting that sentence then OVERSHOT into
"nothing anywhere reads `_run.heartbeat_s` and decides something on it", which
is false in the other direction: `console/console.py`'s `beat_of()` reads that
exact field and its fresh/late/long bands are arithmetic on it — 3x and 20x —
and a record carrying `0`, the documented way to turn the heartbeat off, gets a
band and a sentence of its own. The old wording, "nothing in **this file**",
was the true one, and it is back, with the console named as the reader. What
goes into the document is a fact a reader interprets; only the config key is
private to `metaannot.py`. And `_directory_still_ours()` said in one
sentence that it "essentially never" opens the state file and in the next that
it "never opens the file at all" — with `heartbeat_s` equal to `STATE_PROBE_S`
the margin is tick jitter, so a late tick does let one probe through, which
costs a read and cannot change an answer. Both sentences now say the same
thing.

**The state file's temp follows the `.part` convention.** It was
`<path>.<pid>.<tid>.tmp` in the results root — the one leftover
`find results -name '.*.part.*'` could not see, in the directory an operator
opens cold.

**Two improvements were taken back OUT of this change, to be filed on their
own.** Neither is `#23` and each is a clean standalone change; both were
written while this one was being written, and bundling them made the
concurrency fix bigger than the defect and stretched the guarantee the docs
publish over mechanisms that do not support it.

- *Output sizes in the stage record.* Every `ok`/`adopted` record carried
  `outputs`: `[{path, size}]`, the size each declared output had when the
  record was made, and `decide()` WARNed when a cached output no longer
  matched. It is the only evidence available for the case no ownership check
  can reach — a run `SIGKILL`ed between `atomic_out`'s rename and its record,
  which writes nothing and can prove nothing. It is also warn-only, silent
  whenever the sizes happen to match, a new field in a record whose key set is
  pinned exactly by `tests/test_console_contract.py`, and about sixty lines.
  Worth having; not this issue.
- *Deferred state writes.* An `OSError` in the pre-write read held the keys in
  a module-level `pending` set and let the next merged write carry them,
  rather than the record simply not being written. That is error handling for
  an unrelated failure mode — an `ENOSPC`, an NFS `EIO` — bundled into a
  concurrency fix, and about thirty lines including the carry into every retry
  attempt. Declining the write is what ships; deferring it can be filed with
  the rest of the state file's I/O error handling.

**`doctor --fix` cannot be combined with `--json`**, and the refusal is an
argparse mutually-exclusive group, so it **exits 2**. `die()` exits 1, which is
also "problems found", so a caller could not have told a refused invocation
from a failed check. `--install-plan` composes fine: the file is still written
under `--json` and the document records it; only the line announcing it is
suppressed, since it would otherwise land in stdout beside the JSON.

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
