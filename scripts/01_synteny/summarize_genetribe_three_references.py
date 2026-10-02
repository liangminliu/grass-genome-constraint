#!/usr/bin/env python3
import csv
import gzip
import argparse
from pathlib import Path
from collections import defaultdict

BED_GENE_CACHE = {}
OUTDIR_INDEX_CACHE = {}


def norm_name(s: str) -> str:
    return s.strip().lower()


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def infer_default_sample_list():
    if Path("stat_94sample.list").exists():
        return "stat_94sample.list"
    if Path("draw.sample.list").exists():
        return "draw.sample.list"
    return "draw.sample.list"


def load_genome_sizes(genome_size_file):
    gs = {}
    with open_text(genome_size_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 2:
                continue
            sp = cols[0]
            size_str = cols[1].replace(",", "")
            try:
                gs[norm_name(sp)] = (sp, int(size_str))
            except ValueError:
                continue
    return gs


def build_dir_index(root):
    idx = {}
    for p in Path(root).iterdir():
        if p.is_dir():
            idx[norm_name(p.name)] = p
    return idx


def build_data_index(data_dir):
    idx = {}
    data_dir = Path(data_dir)
    if not data_dir.exists():
        return idx

    for p in data_dir.glob("*.bed"):
        idx[norm_name(p.stem)] = p

    for p in data_dir.glob("*.bed.gz"):
        stem = p.name[:-7]
        idx[norm_name(stem)] = p

    return idx


def canonical_species_name(name, genome_idx, dir_idx, data_bed_idx):
    key = norm_name(name)
    if key in genome_idx:
        return genome_idx[key][0]
    if key in data_bed_idx:
        return data_bed_idx[key].stem.replace(".bed", "")
    if key in dir_idx:
        return dir_idx[key].name
    return name


def resolve_dir(species_name, dir_idx):
    return dir_idx.get(norm_name(species_name))


def resolve_bed(species_name, data_bed_idx, dir_idx):
    key = norm_name(species_name)
    if key in data_bed_idx:
        return data_bed_idx[key]

    sp_dir = dir_idx.get(key)
    if sp_dir:
        candidates = {}
        for p in sp_dir.iterdir():
            if p.is_file():
                candidates[norm_name(p.name)] = p

        exact1 = f"{species_name}.bed"
        exact2 = f"{species_name}.bed.gz"
        if norm_name(exact1) in candidates:
            return candidates[norm_name(exact1)]
        if norm_name(exact2) in candidates:
            return candidates[norm_name(exact2)]

        for p in sp_dir.iterdir():
            if not p.is_file():
                continue
            if p.suffix == ".bed" and norm_name(p.stem) == key:
                return p
            if p.name.endswith(".bed.gz"):
                stem = p.name[:-7]
                if norm_name(stem) == key:
                    return p
    return None


def outdir_file_index(outdir):
    outdir = Path(outdir).resolve()
    key = str(outdir)
    if key in OUTDIR_INDEX_CACHE:
        return OUTDIR_INDEX_CACHE[key]

    idx = {}
    if outdir.exists():
        for p in outdir.iterdir():
            if p.is_file():
                idx[norm_name(p.name)] = p
    OUTDIR_INDEX_CACHE[key] = idx
    return idx


def pick_file_case_insensitive(outdir, filename):
    idx = outdir_file_index(outdir)
    return idx.get(norm_name(filename))


def locate_pair_files(query_name, ref_name, dir_idx):
    """
    Primary:
      {query_dir}/ref_{ref}_output/
    Fallback:
      {ref_dir}/ref_{query}_output/

    File naming inside the chosen output dir may exist in either orientation:
      {ref}_{query}.*
      {query}_{ref}.*
    """
    qkey = norm_name(query_name)
    rkey = norm_name(ref_name)

    query_dir = resolve_dir(query_name, dir_idx)
    ref_dir = resolve_dir(ref_name, dir_idx)

    candidates = []

    if query_dir is not None:
        candidates.append((
            query_dir,
            query_dir / f"ref_{ref_name}_output",
            f"{ref_name}_{query_name}",
            f"{query_name}_{ref_name}",
            "query_dir_primary"
        ))

    if ref_dir is not None and qkey != rkey:
        candidates.append((
            ref_dir,
            ref_dir / f"ref_{query_name}_output",
            f"{ref_name}_{query_name}",
            f"{query_name}_{ref_name}",
            "ref_dir_fallback"
        ))

    best = None
    for owner_dir, outdir, ref_query_base, query_ref_base, source_mode in candidates:
        if not outdir.exists():
            continue

        result = {
            "owner_dir": owner_dir,
            "outdir": outdir,
            "source_mode": source_mode,
            "block_pos_ref_query": pick_file_case_insensitive(outdir, f"{ref_query_base}.block_pos"),
            "block_pos_query_ref": pick_file_case_insensitive(outdir, f"{query_ref_base}.block_pos"),
                        # RBH should normally be ref_query; query_ref is only a fallback.
            "RBH_ref_query": pick_file_case_insensitive(outdir, f"{ref_query_base}.RBH"),
            "RBH_query_ref": pick_file_case_insensitive(outdir, f"{query_ref_base}.RBH"),
            "ref_one2many": pick_file_case_insensitive(outdir, f"{ref_query_base}.one2many"),
            "query_one2many": pick_file_case_insensitive(outdir, f"{query_ref_base}.one2many"),
            "ref_singleton": pick_file_case_insensitive(outdir, f"{ref_query_base}.singleton"),
            "query_singleton": pick_file_case_insensitive(outdir, f"{query_ref_base}.singleton"),
        }

        score = sum(x is not None for x in [
            result["block_pos_ref_query"],
            result["block_pos_query_ref"],
            result["RBH_ref_query"],
            result["RBH_query_ref"],
            result["ref_one2many"],
            result["query_one2many"],
            result["ref_singleton"],
            result["query_singleton"],
        ])

        if best is None or score > best["score"]:
            best = {"score": score, "result": result}

    return None if best is None else best["result"]


def bed_gene_total(bed_path):
    bed_path = Path(bed_path).resolve()
    key = str(bed_path)
    if key in BED_GENE_CACHE:
        return BED_GENE_CACHE[key]

    genes = set()
    with open_text(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) >= 4:
                genes.add(cols[3])

    BED_GENE_CACHE[key] = len(genes)
    return len(genes)


def merge_intervals(intervals):
    if not intervals:
        return 0

    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    total = 0
    cur_s, cur_e = intervals[0]

    for s, e in intervals[1:]:
        if s <= cur_e + 1:
            if e > cur_e:
                cur_e = e
        else:
            total += (cur_e - cur_s + 1)
            cur_s, cur_e = s, e

    total += (cur_e - cur_s + 1)
    return total


def parse_block_pos(block_pos_file):
    n_blocks = 0
    raw_ref_len = 0
    raw_query_len = 0
    ref_intervals = defaultdict(list)
    query_intervals = defaultdict(list)

    with open_text(block_pos_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 6:
                continue

            ref_chr = cols[0]
            qry_chr = cols[3]

            try:
                ref_start = int(cols[1])
                ref_end = int(cols[2])
                qry_start = int(cols[4])
                qry_end = int(cols[5])
            except ValueError:
                continue

            n_blocks += 1
            raw_ref_len += (ref_end - ref_start + 1)
            raw_query_len += (qry_end - qry_start + 1)
            ref_intervals[ref_chr].append((ref_start, ref_end))
            query_intervals[qry_chr].append((qry_start, qry_end))

    merged_ref_len = sum(merge_intervals(v) for v in ref_intervals.values())
    merged_query_len = sum(merge_intervals(v) for v in query_intervals.values())

    return n_blocks, raw_ref_len, raw_query_len, merged_ref_len, merged_query_len


def count_unique_first_col(pair_file):
    genes = set()
    with open_text(pair_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) >= 1:
                genes.add(cols[0])
    return len(genes)


def parse_rbh(rbh_file):
    ref_genes = set()
    query_genes = set()

    with open_text(rbh_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 2:
                continue
            ref_genes.add(cols[0])
            query_genes.add(cols[1])

    return len(ref_genes), len(query_genes)


def count_singleton(singleton_file):
    n = 0
    with open_text(singleton_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            n += 1
    return n


def choose_block_pos_file(located):
    if located.get("block_pos_ref_query") is not None:
        return located["block_pos_ref_query"], "ref_query"
    if located.get("block_pos_query_ref") is not None:
        return located["block_pos_query_ref"], "query_ref"
    return None, None


def parse_block_pos_oriented(block_pos_file, orientation):
    n_blocks, raw_a_len, raw_b_len, merged_a_len, merged_b_len = parse_block_pos(block_pos_file)
    if orientation == "ref_query":
        return n_blocks, raw_a_len, raw_b_len, merged_a_len, merged_b_len
    if orientation == "query_ref":
        return n_blocks, raw_b_len, raw_a_len, merged_b_len, merged_a_len
    raise ValueError(f"Unknown block_pos orientation: {orientation}")


def choose_rbh_file(located):
    # RBH should normally come from ref_query; use query_ref only as fallback.
    if located.get("RBH_ref_query") is not None:
        return located["RBH_ref_query"], "ref_query"
    if located.get("RBH_query_ref") is not None:
        return located["RBH_query_ref"], "query_ref"
    return None, None


def parse_rbh_oriented(rbh_file, orientation):
    n_first, n_second = parse_rbh(rbh_file)
    if orientation == "ref_query":
        return n_first, n_second
    if orientation == "query_ref":
        return n_second, n_first
    raise ValueError(f"Unknown RBH orientation: {orientation}")


def safe_prop(n, d):
    if n in (None, "NA") or d in (None, "NA", 0):
        return "NA"
    return f"{n / d:.6f}"


def main():
    parser = argparse.ArgumentParser(
        description="Summarize GeneTribe main table using directional one2many/RBH logic and NA self-comparisons"
    )
    parser.add_argument(
        "-l", "--sample-list",
        default=infer_default_sample_list(),
        help="sample list file"
    )
    parser.add_argument(
        "-g", "--genome-size",
        default="genome.size",
        help="genome size file"
    )
    parser.add_argument(
        "-d", "--data-dir",
        default="../data",
        help="directory containing *.bed"
    )
    parser.add_argument(
        "-o", "--output",
        default="genetribe_main_summary.tsv",
        help="output TSV"
    )
    parser.add_argument(
        "--refs",
        nargs="+",
        default=["Aegilops_tauschii", "Phyllostachys_edulis", "Zea_mays"],
        help="reference species list"
    )
    args = parser.parse_args()

    root = Path(".").resolve()
    genome_idx = load_genome_sizes(args.genome_size)
    dir_idx = build_dir_index(root)
    data_bed_idx = build_data_index(args.data_dir)

    with open(args.sample_list) as f:
        raw_queries = [x.strip() for x in f if x.strip()]

    out_header = [
        "query_species",
        "ref_species",
        "ref_genome_size_bp",
        "ref_gene_total",
        "query_genome_size_bp",
        "query_gene_total",
        "n_syntenic_blocks",
        "total_syntenic_block_length_ref_bp",
        "total_syntenic_block_length_query_bp",
        "raw_syntenic_block_length_ref_bp",
        "raw_syntenic_block_length_query_bp",
        "n_block_genes_ref",
        "n_block_genes_query",
        "prop_syntenic_genes_ref",
        "prop_syntenic_genes_query",
        "n_RBH_ref",
        "n_RBH_query",
        "prop_RBH_ref",
        "prop_RBH_query",
        "n_singleton_ref",
        "n_singleton_query",
        "prop_singleton_ref",
        "prop_singleton_query",
        "result_source_mode",
        "result_owner_dir",
        "status",
    ]

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(out_header)

        for raw_query in raw_queries:
            query = canonical_species_name(raw_query, genome_idx, dir_idx, data_bed_idx)

            for raw_ref in args.refs:
                ref = canonical_species_name(raw_ref, genome_idx, dir_idx, data_bed_idx)
                status = []

                ref_genome_size = genome_idx.get(norm_name(ref), (ref, "NA"))[1]
                qry_genome_size = genome_idx.get(norm_name(query), (query, "NA"))[1]

                ref_bed = resolve_bed(ref, data_bed_idx, dir_idx)
                qry_bed = resolve_bed(query, data_bed_idx, dir_idx)

                if ref_genome_size == "NA":
                    status.append("missing_ref_genome_size")
                if qry_genome_size == "NA":
                    status.append("missing_query_genome_size")
                if ref_bed is None:
                    status.append("missing_ref_bed")
                if qry_bed is None:
                    status.append("missing_query_bed")

                # -------- self comparison: keep result fields as NA --------
                if norm_name(query) == norm_name(ref):
                    ref_gene_total = "NA" if ref_bed is None else bed_gene_total(ref_bed)
                    qry_gene_total = "NA" if qry_bed is None else bed_gene_total(qry_bed)

                    writer.writerow([
                        query, ref,
                        ref_genome_size, ref_gene_total,
                        qry_genome_size, qry_gene_total,
                        "NA", "NA", "NA", "NA", "NA",
                        "NA", "NA", "NA", "NA",
                        "NA", "NA", "NA", "NA",
                        "NA", "NA", "NA", "NA",
                        "self_comparison",
                        "self",
                        "self_comparison_NA" if not status else "self_comparison_NA;" + ";".join(status)
                    ])
                    continue

                located = locate_pair_files(query, ref, dir_idx)

                if located is None:
                    writer.writerow([
                        query, ref,
                        ref_genome_size, "NA" if ref_bed is None else bed_gene_total(ref_bed),
                        qry_genome_size, "NA" if qry_bed is None else bed_gene_total(qry_bed),
                        "NA", "NA", "NA", "NA", "NA",
                        "NA", "NA", "NA", "NA",
                        "NA", "NA", "NA", "NA",
                        "NA", "NA", "NA", "NA",
                        "NA", "NA",
                        "pair_output_not_found" + (";" + ";".join(status) if status else "")
                    ])
                    continue

                owner_dir = located["owner_dir"].name
                source_mode = located["source_mode"]

                block_pos_file, block_pos_orientation = choose_block_pos_file(located)
                rbh_file, rbh_orientation = choose_rbh_file(located)
                ref_one2many_file = located["ref_one2many"]
                qry_one2many_file = located["query_one2many"]
                ref_singleton_file = located["ref_singleton"]
                qry_singleton_file = located["query_singleton"]

                if block_pos_file is None:
                    status.append("missing_block_pos")

                if rbh_file is None:
                    status.append("missing_RBH")
                if ref_one2many_file is None:
                    status.append("missing_ref_one2many")
                if qry_one2many_file is None:
                    status.append("missing_query_one2many")
                if ref_singleton_file is None:
                    status.append("missing_ref_singleton")
                if qry_singleton_file is None:
                    status.append("missing_query_singleton")

                hard_missing = (
                    ref_bed is None or
                    qry_bed is None or
                    block_pos_file is None or
                    rbh_file is None or
                    ref_one2many_file is None or
                    qry_one2many_file is None
                )

                if hard_missing:
                    writer.writerow([
                        query, ref,
                        ref_genome_size, "NA" if ref_bed is None else bed_gene_total(ref_bed),
                        qry_genome_size, "NA" if qry_bed is None else bed_gene_total(qry_bed),
                        "NA", "NA", "NA", "NA", "NA",
                        "NA", "NA", "NA", "NA",
                        "NA", "NA", "NA", "NA",
                        "NA" if ref_singleton_file is None else count_singleton(ref_singleton_file),
                        "NA" if qry_singleton_file is None else count_singleton(qry_singleton_file),
                        "NA", "NA",
                        source_mode,
                        owner_dir,
                        ";".join(status) if status else "incomplete_files"
                    ])
                    continue

                try:
                    ref_gene_total = bed_gene_total(ref_bed)
                except Exception as e:
                    ref_gene_total = "NA"
                    status.append(f"ref_bed_error={e}")

                try:
                    qry_gene_total = bed_gene_total(qry_bed)
                except Exception as e:
                    qry_gene_total = "NA"
                    status.append(f"qry_bed_error={e}")

                try:
                    n_blocks, raw_ref_len, raw_qry_len, merged_ref_len, merged_qry_len = parse_block_pos_oriented(
                        block_pos_file, block_pos_orientation
                    )
                except Exception as e:
                    n_blocks, raw_ref_len, raw_qry_len, merged_ref_len, merged_qry_len = ("NA",) * 5
                    status.append(f"block_pos_error={e}")

                try:
                    n_block_genes_ref = count_unique_first_col(ref_one2many_file)
                except Exception as e:
                    n_block_genes_ref = "NA"
                    status.append(f"ref_one2many_error={e}")

                try:
                    n_block_genes_qry = count_unique_first_col(qry_one2many_file)
                except Exception as e:
                    n_block_genes_qry = "NA"
                    status.append(f"query_one2many_error={e}")

                try:
                    n_rbh_ref, n_rbh_qry = parse_rbh_oriented(rbh_file, rbh_orientation)
                except Exception as e:
                    n_rbh_ref, n_rbh_qry = "NA", "NA"
                    status.append(f"rbh_error={e}")

                try:
                    n_singleton_ref = count_singleton(ref_singleton_file) if ref_singleton_file else "NA"
                except Exception as e:
                    n_singleton_ref = "NA"
                    status.append(f"ref_singleton_error={e}")

                try:
                    n_singleton_qry = count_singleton(qry_singleton_file) if qry_singleton_file else "NA"
                except Exception as e:
                    n_singleton_qry = "NA"
                    status.append(f"qry_singleton_error={e}")

                prop_syn_ref = safe_prop(n_block_genes_ref, ref_gene_total)
                prop_syn_qry = safe_prop(n_block_genes_qry, qry_gene_total)
                prop_rbh_ref = safe_prop(n_rbh_ref, ref_gene_total)
                prop_rbh_qry = safe_prop(n_rbh_qry, qry_gene_total)
                prop_singleton_ref = safe_prop(n_singleton_ref, ref_gene_total)
                prop_singleton_qry = safe_prop(n_singleton_qry, qry_gene_total)

                writer.writerow([
                    query,
                    ref,
                    ref_genome_size,
                    ref_gene_total,
                    qry_genome_size,
                    qry_gene_total,
                    n_blocks,
                    merged_ref_len,
                    merged_qry_len,
                    raw_ref_len,
                    raw_qry_len,
                    n_block_genes_ref,
                    n_block_genes_qry,
                    prop_syn_ref,
                    prop_syn_qry,
                    n_rbh_ref,
                    n_rbh_qry,
                    prop_rbh_ref,
                    prop_rbh_qry,
                    n_singleton_ref,
                    n_singleton_qry,
                    prop_singleton_ref,
                    prop_singleton_qry,
                    source_mode,
                    owner_dir,
                    "OK" if not status else ";".join(status),
                ])


if __name__ == "__main__":
    main()