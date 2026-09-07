"""Generators for every input metaannot reads.

Nothing here is committed as data: each function writes a file from a fixed
seed, so a fixture is reproducible without a binary in git. Where a real tool
would have produced the file (hmmsearch, DIAMOND, Foldseek, ...), the writer
reproduces that tool's on-disk format exactly, because the parsers are pinned
to column positions.
"""
from __future__ import annotations

import os
import random
from dataclasses import dataclass, field

AA = "ACDEFGHIKLMNPQRSTVWY"

# eggnog-mapper 2.1.x column order, which is also metaannot's CANONICAL_V2.
EMAPPER_COLS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
]


# ----------------------------------------------------------------------
# proteins
# ----------------------------------------------------------------------
@dataclass
class Protein:
    """One synthetic protein plus the evidence it is meant to carry."""
    pid: str
    seq: str
    ko: str = ""                 # KEGG_ko cell, e.g. "ko:K01234"
    pathway: str = ""            # KEGG_Pathway cell, e.g. "ko00010,map00010"
    pfams_emapper: str = ""      # eggNOG's own PFAMs cell
    cazy: str = ""
    seed_taxid: str = ""
    description: str = ""
    expect_bin: str = "4_dark"   # what build_annotation should decide
    in_emapper: bool = True      # False = no row at all in the eggNOG table


def _seq(rng, n):
    return "".join(rng.choice(AA) for _ in range(n))


def protein_set(seed=1, n_extra=0):
    """A protein per bin, plus optional filler.

    The bins are decided by the evidence attached here, so a change to the
    classifier shows up as a named protein landing in the wrong bin rather
    than as a percentage moving.
    """
    rng = random.Random(seed)
    ps = [
        # KO with a specific (non-global) map
        Protein("P_ko_path", _seq(rng, 240), ko="ko:K01234",
                pathway="ko00010,map00010", seed_taxid="820",
                description="phosphoglucomutase", expect_bin="1_ko_pathway"),
        # KO whose ONLY map is a global one -> must not count as a pathway
        Protein("P_ko_global", _seq(rng, 230), ko="ko:K02345",
                pathway="ko01100,map01100", seed_taxid="820",
                description="global map only", expect_bin="2_ko_orphan"),
        # KO with no map at all
        Protein("P_ko_orphan", _seq(rng, 210), ko="ko:K03456", pathway="",
                seed_taxid="821", description="orphan KO",
                expect_bin="2_ko_orphan"),
        # no KO, an informative eggNOG Pfam only
        Protein("P_eggpfam", _seq(rng, 200), pfams_emapper="Peptidase_S8",
                seed_taxid="821", description="subtilase",
                expect_bin="3_annotated_no_ko"),
        # no KO, eggNOG Pfam that is a DUF -> the DUF bin, not the annotated one
        Protein("P_eggduf", _seq(rng, 190), pfams_emapper="DUF1234",
                seed_taxid="821", description="DUF1234 family",
                expect_bin="3d_duf_only"),
        # no KO, eggNOG Pfam that is a UPF -> same rule as DUF
        Protein("P_eggupf", _seq(rng, 185), pfams_emapper="UPF0102",
                seed_taxid="822", expect_bin="3d_duf_only"),
        # CAZy only
        Protein("P_cazy", _seq(rng, 300), cazy="GH13", seed_taxid="822",
                expect_bin="3_annotated_no_ko"),
        # nothing at all, but present in the eggNOG table
        Protein("P_dark1", _seq(rng, 95), seed_taxid="823", expect_bin="4_dark"),
        # nothing at all, absent from the eggNOG table
        Protein("P_dark2", _seq(rng, 88), in_emapper=False, expect_bin="4_dark"),
        # dark, small, and taxon-less: the shared-peptide veto case
        Protein("P_dark3", _seq(rng, 70), in_emapper=False, expect_bin="4_dark"),
    ]
    for i in range(n_extra):
        ps.append(Protein(f"P_fill{i:04d}", _seq(rng, rng.randint(60, 400)),
                          ko="ko:K0%04d" % (1000 + i),
                          pathway="ko00020,map00020",
                          seed_taxid=str(830 + (i % 7)),
                          expect_bin="1_ko_pathway"))
    return ps


def write_fasta(path, proteins, description=True):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for p in proteins:
            desc = f" {p.description}" if description and p.description else ""
            fh.write(f">{p.pid}{desc}\n{p.seq}\n")
    return path


def write_emapper(path, proteins, id_prefix="", header=True):
    """A precomputed .emapper.annotations table.

    `id_prefix` is REMOVED from the ids written here, so the table looks like
    one built before a search database prefixed its identifiers.
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("## emapper-2.1.12\n")
        if header:
            fh.write("#" + "\t".join(EMAPPER_COLS) + "\n")
        for p in proteins:
            if not p.in_emapper:
                continue
            pid = p.pid[len(id_prefix):] if id_prefix and p.pid.startswith(id_prefix) \
                else p.pid
            row = {c: "-" for c in EMAPPER_COLS}
            row["query"] = pid
            row["seed_ortholog"] = f"{p.seed_taxid or '820'}.SEED{pid}"
            row["evalue"] = "1e-50"
            row["score"] = "200.0"
            row["eggNOG_OGs"] = "COG0001@1|root,COG0001@2|Bacteria"
            row["max_annot_lvl"] = "2|Bacteria"
            row["COG_category"] = "S"
            row["Description"] = p.description or "-"
            row["KEGG_ko"] = p.ko or "-"
            row["KEGG_Pathway"] = p.pathway or "-"
            row["CAZy"] = p.cazy or "-"
            row["PFAMs"] = p.pfams_emapper or "-"
            fh.write("\t".join(row[c] for c in EMAPPER_COLS) + "\n")
    return path


# ----------------------------------------------------------------------
# quantification
# ----------------------------------------------------------------------
def peptides_for(proteins, per_protein=3, seed=2):
    """(peptide, razor, [candidates]) with a couple of deliberately shared."""
    rng = random.Random(seed)
    rows = []
    for p in proteins:
        for k in range(per_protein):
            pep = "".join(rng.choice(AA) for _ in range(9)) + "K"
            rows.append({"peptide": f"{pep}{k}", "razor": p.pid,
                         "candidates": [p.pid]})
    return rows


def write_peptide_table(path, proteins, samples, rows=None, seed=3,
                        zeros=(), shared=(), protein_id_with_desc=True,
                        suffix=" Intensity", extra_cols=None):
    """A FragPipe combined_peptide.tsv.

    `zeros` is a set of (row_index, sample) written as a literal 0, which is
    what FragPipe writes for "not quantified".
    `shared` is a list of (row_index, [extra protein ids]) put into
    'Mapped Proteins'.
    """
    rng = random.Random(seed)
    rows = rows if rows is not None else peptides_for(proteins, seed=seed)
    extra = dict(shared)
    cols = (["Peptide Sequence", "Protein", "Protein ID", "Mapped Proteins"]
            + list((extra_cols or {}).keys())
            + [f"{s}{suffix}" for s in samples])
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for i, r in enumerate(rows):
            desc = " hypothetical protein" if protein_id_with_desc else ""
            vals = []
            for s in samples:
                if (i, s) in set(zeros):
                    vals.append("0")
                else:
                    vals.append(f"{rng.randint(10000, 100000)}.0")
            rec = [r["peptide"], r["razor"], r["razor"] + desc,
                   ";".join(extra.get(i, []))]
            rec += [str(col[i]) for col in (extra_cols or {}).values()]
            fh.write("\t".join(rec + vals) + "\n")
    return path


def write_tmt_peptide_table(path, proteins, channels, seed=4):
    """The shape a FragPipe TMT per-plex table has: reporter channels named
    'Intensity <sample>' and one bare MS1 'Intensity'."""
    rng = random.Random(seed)
    cols = ["Peptide Sequence", "Protein", "Intensity"] + \
           [f"Intensity {c}" for c in channels]
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for i, p in enumerate(proteins):
            vals = [str(rng.randint(1000, 9000)) for _ in range(len(channels) + 1)]
            fh.write("\t".join([f"PEPTIDEK{i}", p.pid] + vals) + "\n")
    return path


def write_manifest(path, entries):
    """entries: list of (file, experiment, bioreplicate, data_type)."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for e in entries:
            fh.write("\t".join(str(x) for x in e) + "\n")
    return path


def split_sample(name):
    """'A_1' -> ('A', '1'). FragPipe names a quant column
    <experiment>_<bioreplicate>, so the fixture has to write the two halves
    that reassemble into the column name."""
    exp, _, rep = str(name).rpartition("_")
    return (exp, rep) if exp else (name, "")


def simple_manifest(path, samples, groups=None, fractions=1, data_type="DDA"):
    """One row per raw FILE. With fractions > 1 a sample name repeats across
    its fraction rows, which is exactly how FragPipe denotes a fractionated
    acquisition."""
    entries = []
    for s in samples:
        exp, rep = split_sample(s)
        if groups and s in groups:
            exp = groups[s]
        dt = data_type[s] if isinstance(data_type, dict) else data_type
        for f in range(fractions):
            entries.append((f"/raw/{s}_f{f}.mzML", exp, rep, dt))
    return write_manifest(path, entries)


# ----------------------------------------------------------------------
# tool outputs
# ----------------------------------------------------------------------
def tblout_line(target, query_name, query_acc, evalue=1e-30, score=100.0,
                desc="-"):
    """One hmmsearch --tblout row: 18 fixed columns then a free-text
    description of the TARGET."""
    f = [target, "-", query_name, query_acc, f"{evalue:g}", f"{score:.1f}",
         "0.0", f"{evalue:g}", f"{score:.1f}", "0.0",
         "1.0", "1", "0", "0", "1", "1", "1", "1", desc]
    return " ".join(f)


def write_tblout(path, hits, header=True):
    """hits: list of (target, query_name, query_acc[, desc])."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        if header:
            fh.write("#                                                     "
                     "--- full sequence ---\n")
            fh.write("# target name  accession  query name  accession"
                     "    E-value  score  bias\n")
        for h in hits:
            target, qn, qa = h[0], h[1], h[2]
            desc = h[3] if len(h) > 3 else "-"
            fh.write(tblout_line(target, qn, qa, desc=desc) + "\n")
        fh.write("#\n# Program:         hmmsearch\n")
    return path


def write_domtblout(path, hits):
    """hits: list of (target, query_name, qlen, hmm_from, hmm_to, i_evalue).
    23-column hmmsearch --domtblout."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# target name accession tlen query name accession qlen\n")
        for target, qn, qlen, hf, ht, ev in hits:
            f = [target, "-", "300", qn, "-", str(qlen),
                 "1e-20", "80.0", "0.0", "1", "1", f"{ev:g}", f"{ev:g}",
                 "80.0", "0.0", str(hf), str(ht), "10", "200", "10", "200",
                 "0.95", "-"]
            fh.write(" ".join(f) + "\n")
    return path


def write_hmm_library(path, models):
    """models: list of (name, acc, desc). An HMMER3 text library, enough for
    parse_hmm_lib_desc."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for name, acc, desc in models:
            fh.write("HMMER3/f [3.4 | Aug 2023]\n")
            fh.write(f"NAME  {name}\n")
            if acc:
                fh.write(f"ACC   {acc}\n")
            fh.write(f"DESC  {desc}\n")
            fh.write("LENG  100\n//\n")
    return path


def write_diamond(path, hits):
    """hits: list of (qseqid, sseqid, pident, evalue, bitscore, qcov, stitle).
    The 9-column custom format stage_diamond asks for."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for q, s, pid, ev, bs, qcov, title in hits:
            fh.write("\t".join([q, s, f"{pid}", "150", f"{ev:g}", f"{bs}",
                                f"{qcov}", "80", title]) + "\n")
    return path


def write_signalp(path, calls):
    """calls: {protein_id: class}. SignalP 6.0 prediction_results.txt."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# SignalP-6.0\tOrganism: other\n")
        fh.write("# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n")
        for pid, cls in calls.items():
            cs = "CS pos: 22-23. AQA-QT" if cls != "OTHER" else ""
            fh.write(f"{pid}\t{cls}\t0.01\t0.99\t{cs}\n")
    return path


def write_tmbed(path, preds):
    """preds: {protein_id: label_string}. tmbed --out-format 0 (3 lines)."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for pid, labels in preds.items():
            fh.write(f">{pid}\n{'A' * len(labels)}\n{labels}\n")
    return path


FOLDSEEK_FIELDS = ["query", "target", "fident", "alnlen", "evalue", "bits",
                   "prob", "alntmscore", "lddt", "theader"]


def write_foldseek(path, hits, with_db_col=True, header=True):
    """hits: list of (query, target, evalue, bits, prob, tm, theader, db).
    The 10-field legacy layout stage_foldseek actually writes, plus the
    target_db column it appends."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        if header:
            cols = FOLDSEEK_FIELDS + (["target_db"] if with_db_col else [])
            fh.write("#" + "\t".join(cols) + "\n")
        for q, t, ev, bits, prob, tm, theader, db in hits:
            row = [q, t, "0.35", "180", f"{ev:g}", f"{bits}", f"{prob}",
                   f"{tm}", "0.8", theader]
            if with_db_col:
                row.append(db)
            fh.write("\t".join(row) + "\n")
    return path


def write_kofam(path, calls):
    """calls: list of (gene, ko, score, evalue, definition). KOfamScan
    detail-tsv, where a leading '*' marks a call above the family threshold."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("#\tgene name\tKO\tthrshld\tscore\tE-value\tKO definition\n")
        for gene, ko, score, ev, defn in calls:
            fh.write("\t".join(["*", gene, ko, "100.0", f"{score}", f"{ev:g}",
                                f'"{defn}"']) + "\n")
    return path


def write_interpro(path, rows):
    """rows: list of (protein, analysis, signature, desc, ipr, ipr_desc, go).
    A 14-column InterProScan TSV (-iprlookup -goterms)."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for prot, ana, sig, desc, ipr, iprd, go in rows:
            fh.write("\t".join([prot, "d41d8cd9", "300", ana, sig, desc,
                                "10", "120", "1.0E-20", "T", "01-01-2026",
                                ipr, iprd, go]) + "\n")
    return path


def write_hhr(directory, query, hits):
    """hits: list of (hit, desc, prob, evalue). One .hhr per query."""
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, f"{query}.hhr")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"Query         {query}\nMatch_columns 200\n\n")
        fh.write(" No Hit                             Prob E-value P-value  "
                 "Score    SS Cols Query HMM  Template HMM\n")
        for i, (hit, desc, prob, ev) in enumerate(hits, 1):
            name = f"{hit} {desc}".strip()
            fh.write(f"{i:3d} {name:<30.30s} {prob:5.1f} {ev:7g} {ev:7g} "
                     f"{100.0:6.1f} {0.0:5.1f} {150:4d}  10-160    5-155 (200)\n")
        fh.write("\n")
    return path


def write_cluster(path, members):
    """members: {member: representative}. MMseqs easy-cluster tsv is
    rep <TAB> member."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for member, rep in members.items():
            fh.write(f"{rep}\t{member}\n")
    return path


def write_unipept(path, rows, modern=True, sep=","):
    """rows: list of (peptide, taxon_id, rank, {rank: taxid}).

    `modern` writes the current CLI's `domain_id`; otherwise the legacy gem's
    `superkingdom_id`.
    """
    top = "domain_id" if modern else "superkingdom_id"
    cols = ["peptide", "taxon_id", "taxon_name", "taxon_rank", top,
            "phylum_id", "class_id", "order_id", "family_id", "genus_id",
            "species_id"]
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(sep.join(cols) + "\n")
        for pep, tid, rank, lin in rows:
            rec = [pep, str(tid), "Name", rank,
                   str(lin.get("domain", "")), str(lin.get("phylum", "")),
                   str(lin.get("class", "")), str(lin.get("order", "")),
                   str(lin.get("family", "")), str(lin.get("genus", "")),
                   str(lin.get("species", ""))]
            fh.write(sep.join(rec) + "\n")
    return path


def write_taxdump(directory, nodes, names, merged=None, deleted=()):
    """nodes: {taxid: (parent, rank)}; names: {taxid: scientific name}."""
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, "nodes.dmp"), "w", encoding="utf-8") as fh:
        for tid, (parent, rank) in nodes.items():
            fh.write(f"{tid}\t|\t{parent}\t|\t{rank}\t|\n")
    with open(os.path.join(directory, "names.dmp"), "w", encoding="utf-8") as fh:
        for tid, nm in names.items():
            fh.write(f"{tid}\t|\t{nm}\t|\t\t|\tscientific name\t|\n")
    if merged:
        with open(os.path.join(directory, "merged.dmp"), "w",
                  encoding="utf-8") as fh:
            for old, new in merged.items():
                fh.write(f"{old}\t|\t{new}\t|\n")
    if deleted:
        with open(os.path.join(directory, "delnodes.dmp"), "w",
                  encoding="utf-8") as fh:
            for t in deleted:
                fh.write(f"{t}\t|\n")
    return directory


def toy_taxdump(directory):
    """A three-genus tree deep enough for the lineage comparison."""
    nodes = {
        "1": ("1", "no rank"),
        "2": ("1", "domain"),
        "1239": ("2", "phylum"),
        "91061": ("1239", "class"),
        "186826": ("91061", "order"),
        "81852": ("186826", "family"),
        "1350": ("81852", "genus"),        # Enterococcus
        "1351": ("1350", "species"),       # E. faecalis
        "1352": ("1350", "species"),       # E. faecium
        "816": ("2", "genus"),             # Bacteroides
        "820": ("816", "species"),         # B. uniformis
        "821": ("816", "species"),         # B. thetaiotaomicron
        "822": ("816", "species"),
        "823": ("816", "species"),
    }
    names = {t: f"taxon{t}" for t in nodes}
    names.update({"1350": "Enterococcus", "1351": "Enterococcus faecalis",
                  "816": "Bacteroides", "820": "Bacteroides uniformis"})
    return write_taxdump(directory, nodes, names,
                         merged={"9999": "820"}, deleted=["8888"])


# ----------------------------------------------------------------------
# ground-truth quantitative datasets
# ----------------------------------------------------------------------
@dataclass
class GroundTruth:
    """A protein-level matrix with a known organism shift and known
    regulation, for the taxon-reference tests."""
    samples: list
    groups: dict
    taxon_of: dict
    truly_regulated: set
    passengers: set
    shift_log2: float
    frame: object = field(default=None)


def ground_truth_taxon_shift(n_members=12, n_regulated=3, n_per_group=6,
                             shift_log2=1.0, effect_log2=2.0, seed=11,
                             other_taxa=3, noise=0.03):
    """One taxon doubles between groups; a named subset of it is genuinely
    regulated on top of that shift.

    Returns (DataFrame of intensities, GroundTruth). The passengers move by
    exactly `shift_log2`; the regulated proteins by `shift_log2 + effect_log2`.
    A summed taxon reference is dragged by the regulated members, so the
    passengers do NOT collapse to zero under it — which is the defect the
    median of ratios fixes.
    """
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    samples = [f"S{i+1}" for i in range(2 * n_per_group)]
    groups = {s: ("ctrl" if i < n_per_group else "case")
              for i, s in enumerate(samples)}
    is_case = np.array([groups[s] == "case" for s in samples], dtype=float)

    rows, taxon_of, regulated, passengers = {}, {}, set(), set()
    for i in range(n_members):
        pid = f"T1_p{i:02d}"
        taxon_of[pid] = "T1"
        base = 12.0 + rng.normal(0, 1.0)
        lfc = shift_log2
        if i < n_regulated:
            lfc += effect_log2
            regulated.add(pid)
        else:
            passengers.add(pid)
        vals = base + lfc * is_case + rng.normal(0, noise, len(samples))
        rows[pid] = 2.0 ** vals
    for t in range(other_taxa):
        for i in range(n_members):
            pid = f"T{t+2}_p{i:02d}"
            taxon_of[pid] = f"T{t+2}"
            base = 12.0 + rng.normal(0, 1.0)
            rows[pid] = 2.0 ** (base + rng.normal(0, noise, len(samples)))
    df = pd.DataFrame(rows, index=samples).T
    df.index.name = "group_id"
    return df, GroundTruth(samples=samples, groups=groups, taxon_of=taxon_of,
                           truly_regulated=regulated, passengers=passengers,
                           shift_log2=shift_log2, frame=df)


# ----------------------------------------------------------------------
# the other quantification formats
# ----------------------------------------------------------------------
def write_diann_matrix(path, proteins, samples, seed=21):
    """DIA-NN report.pg_matrix.tsv."""
    rng = random.Random(seed)
    cols = ["Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
            "First.Protein.Description"] + list(samples)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for p in proteins:
            vals = [str(rng.randint(10000, 100000)) for _ in samples]
            fh.write("\t".join([p.pid, p.pid, p.pid, "", p.description or ""]
                               + vals) + "\n")
    return path


def write_fragpipe_protein(path, proteins, samples, seed=22):
    """FragPipe combined_protein.tsv."""
    rng = random.Random(seed)
    cols = ["Protein", "Protein ID", "Gene", "Description",
            "Indistinguishable Proteins", "Protein Length"] + \
           [f"{s} Intensity" for s in samples]
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for p in proteins:
            vals = [str(rng.randint(10000, 100000)) for _ in samples]
            fh.write("\t".join([p.pid, f"{p.pid} {p.description}", "",
                                p.description or "", "", str(len(p.seq))]
                               + vals) + "\n")
    return path


def write_ion_table(path, proteins, samples, seed=23):
    """FragPipe combined_ion.tsv: one row per (modified peptide, charge)."""
    rng = random.Random(seed)
    cols = ["Modified Sequence", "Peptide Sequence", "Charge", "Protein",
            "Protein ID", "Mapped Proteins"] + \
           [f"{s} Intensity" for s in samples]
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for p in proteins:
            for k in range(3):
                pep = "".join(rng.choice(AA) for _ in range(9)) + "K"
                for z in (2, 3):
                    vals = [str(rng.randint(10000, 100000)) for _ in samples]
                    fh.write("\t".join([f"{pep}[{k}]", pep, str(z), p.pid,
                                        f"{p.pid} {p.description}", ""]
                                       + vals) + "\n")
    return path


def _long_rows(proteins, samples, groups, seed):
    rng = random.Random(seed)
    for p in proteins:
        for k in range(3):
            pep = "".join(rng.choice(AA) for _ in range(9)) + "K"
            for s in samples:
                yield p.pid, pep, s, groups[s], rng.randint(10000, 100000)


def write_msstats_csv(path, proteins, samples, groups, seed=24):
    """FragPipe label-free MSstats.csv (long, feature level)."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,"
                 "ProductCharge,IsotopeLabelType,Condition,BioReplicate,Run,"
                 "Intensity\n")
        for prot, pep, run, grp, val in _long_rows(proteins, samples, groups,
                                                   seed):
            fh.write(f"{prot},{pep},2,NA,NA,L,{grp},{run},{run},{val}\n")
    return path


def write_msstats_feature(path, proteins, samples, groups, seed=25):
    """dataProcess()$FeatureLevelData exported as TSV."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("PROTEIN\tPEPTIDE\tFEATURE\tRUN\tGROUP_ORIGINAL\t"
                 "SUBJECT_ORIGINAL\tINTENSITY\n")
        for prot, pep, run, grp, val in _long_rows(proteins, samples, groups,
                                                   seed):
            fh.write(f"{prot}\t{pep}\t{pep}_2\t{run}\t{grp}\t{run}\t{val}\n")
    return path


def write_msstats_protein(path, proteins, samples, groups, seed=26):
    """dataProcess()$ProteinLevelData exported as TSV (log2 intensities)."""
    import math
    rng = random.Random(seed)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("Protein\toriginalRUN\tGROUP_ORIGINAL\tSUBJECT_ORIGINAL\t"
                 "LogIntensities\n")
        for p in proteins:
            for s in samples:
                v = math.log2(rng.randint(10000, 100000))
                fh.write(f"{p.pid}\t{s}\t{groups[s]}\t{s}\t{v:.6f}\n")
    return path
