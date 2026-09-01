#!/usr/bin/env python3
"""HBA1/HBA2/HBB transcript subset GTF by tag priority (§12.4)."""
import argparse, sys

GLOBIN_GENES = {"HBA1", "HBA2", "HBB"}
MANE_TAGS = {"MANE_Select", "MANE_Plus_Clinical"}
FALLBACK_TAG = "Ensembl_canonical"

def parse_attrs(s):
    out = {}
    for part in s.rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        k, _, v = part.partition(" ")
        out.setdefault(k.strip(), []).append(v.strip().strip('"'))
    return out

def get_one(a, key):
    v = a.get(key, [])
    return v[0] if v else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-g", "--gtf", required=True)
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    selected = {}    # gname -> list of (tid, mode)
    genes_by_name = {}  # gname -> set of gid
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "transcript":
                continue
            a = parse_attrs(c[8])
            gname = get_one(a, "gene_name")
            if gname not in GLOBIN_GENES:
                continue
            tid = get_one(a, "transcript_id")
            if not tid:
                continue
            gid = get_one(a, "gene_id") or ""
            genes_by_name.setdefault(gname, set()).add(gid)
            tags = set()
            for t in a.get("tag", []):
                tags.update(t.split(","))
            if tags & MANE_TAGS:            # set intersection
                mode = "MANE"
            elif FALLBACK_TAG in tags:
                mode = "Ensembl_canonical"
            else:
                mode = "bare_gene_name"
            selected.setdefault(gname, []).append((tid, mode))

    if not selected:
        print("WARNING: no HBA1/HBA2/HBB transcripts found in " + args.gtf, file=sys.stderr)

    chosen_by_gene = {}
    for gname, entries in selected.items():
        manes = [tid for tid, m in entries if m == "MANE"]
        canon = [tid for tid, m in entries if m == "Ensembl_canonical"]
        if manes:
            chosen_by_gene[gname] = (manes, "MANE")
        elif canon:
            chosen_by_gene[gname] = (canon, "Ensembl_canonical")
        else:
            chosen_by_gene[gname] = ([tid for tid, _ in entries], "bare_gene_name")
            print(f"WARNING[{gname}]: no MANE/Ensembl_canonical tag; "
                  f"selected {len(chosen_by_gene[gname][0])} transcript(s) by bare "
                  f"gene_name — multiplexes transcripts per gene (SPEC §12.4", file=sys.stderr)

    chosen_tids = set()
    for tids, _ in chosen_by_gene.values():
        chosen_tids.update(tids)
    chosen_gids = set()
    for gname, _ in chosen_by_gene.items():
        chosen_gids.update(genes_by_name.get(gname, set()))

    with open(args.gtf) as fh, open(args.out, "w")as out:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            a = parse_attrs(c[8])
            tidc = get_one(a, "transcript_id")
            gidc = get_one(a, "gene_id")
            if c[2] == "gene":
                keep = gidc in chosen_gids
            elif c[2] == "transcript":
                keep = tidc in chosen_tids
            else:
                keep = tidc in chosen_tids
            if keep:
                out.write(line)

    n = len(chosen_tids)
    print(f"Selected {n} globin transcript(s)")
    for gname in chosen_by_gene:
        tids_sel = chosen_by_gene[gname][0]
        mode_sel = chosen_by_gene[gname][1]
        print("  " + gname + ": " + str(len(tids_sel)) + " transcript(s), mode=" + mode_sel)

if __name__ == "__main__":
    main()
