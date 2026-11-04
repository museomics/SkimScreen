#!/usr/bin/env python3
"""
Plot Kraken2 report results across N samples.

Expects one report per sample sub-directory:
    <indir>/<sample>/kraken_report.txt

Outputs (in --outdir):
    kraken_reports_long.tsv            every report line, all samples (long format)
    kraken_<rank>_reads.tsv            samples x taxa clade read counts at --rank
    kraken_<rank>_pct.tsv              same, as % of total reads
    kraken_overview.tsv                unclassified / domain read counts and %
    stacked_bar_<rank>.png/.pdf        top --top taxa per sample at --rank
    heatmap_<rank>.png/.pdf            samples x top --heatmap-top taxa (log colour scale)
    classification_overview.png/.pdf   unclassified + Eukaryota/Bacteria/Archaea/Viruses

All percentages are of TOTAL reads (classified + unclassified), so bars sum to 100%.

Usage:
    plot_kraken_reports.py --indir kraken_out --outdir kraken_plots --rank G --top 20

Requires: pandas, numpy, matplotlib (scipy optional, for heatmap clustering).
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

plt.rcParams["pdf.fonttype"] = 42  # editable text in PDFs
plt.rcParams["ps.fonttype"] = 42

RANK_ALIASES = {
    "domain": "D", "superkingdom": "D", "kingdom": "K", "phylum": "P",
    "class": "C", "order": "O", "family": "F", "genus": "G", "species": "S",
}
RANK_NAMES = {"D": "domain", "K": "kingdom", "P": "phylum", "C": "class",
              "O": "order", "F": "family", "G": "genus", "S": "species"}
DOMAINS = {2759: "Eukaryota", 2: "Bacteria", 2157: "Archaea", 10239: "Viruses"}

COL_UNCLASSIFIED = "#d9d9d9"
COL_NO_RANK = "#969696"
COL_OTHER = "#525252"


# ── Parsing ───────────────────────────────────────────────────────────────────
def parse_report(path):
    """Parse a Kraken2 report (standard 6-col or --report-minimizer-data 8-col)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) == 6:
                pct, clade, direct, rank, taxid, name = parts
            elif len(parts) == 8:
                pct, clade, direct, _, _, rank, taxid, name = parts
            else:
                continue
            try:
                pct, clade, direct, taxid = float(pct), int(clade), int(direct), int(taxid)
            except ValueError:
                continue  # header line or malformed
            depth = (len(name) - len(name.lstrip(" "))) // 2
            rows.append((pct, clade, direct, rank.strip(), taxid, name.strip(), depth))
    return pd.DataFrame(rows, columns=["pct", "clade_reads", "direct_reads",
                                       "rank", "taxid", "name", "depth"])


def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def load_reports(indir, report_name):
    paths = sorted(Path(indir).glob(f"*/{report_name}"), key=lambda p: natural_key(p.parent.name))
    if not paths:
        sys.exit(f"ERROR: no '{report_name}' files found in {indir}/*/")
    frames = []
    for p in paths:
        df = parse_report(p)
        if df.empty:
            print(f"WARNING: no parseable lines in {p} — skipping", file=sys.stderr)
            continue
        df.insert(0, "sample", p.parent.name)
        frames.append(df)
    if not frames:
        sys.exit("ERROR: none of the reports could be parsed")
    return pd.concat(frames, ignore_index=True)


def sample_totals(long):
    """Total, unclassified and classified (root clade) reads per sample."""
    def get(taxid):
        return (long[long.taxid == taxid].groupby("sample")["clade_reads"].sum()
                .reindex(samples_in_order(long), fill_value=0))
    unclassified, classified = get(0), get(1)
    return pd.DataFrame({"total": unclassified + classified,
                         "unclassified": unclassified, "classified": classified})


def samples_in_order(long):
    return list(dict.fromkeys(long["sample"]))


# ── Tables ────────────────────────────────────────────────────────────────────
def rank_table(long, rank, totals):
    """samples x taxa clade reads at an exact rank code, plus unassigned categories."""
    sub = long[long["rank"] == rank].copy()
    # Same name at different taxids (e.g. Bacillus the bacterium vs the stick insect)
    dup = sub.groupby("name")["taxid"].nunique()
    dup_names = set(dup[dup > 1].index)
    sub["label"] = [f"{n} ({t})" if n in dup_names else n for n, t in zip(sub.name, sub.taxid)]
    reads = (sub.pivot_table(index="sample", columns="label", values="clade_reads",
                             aggfunc="sum", fill_value=0)
             .reindex(totals.index, fill_value=0))
    no_rank = (totals["classified"] - reads.sum(axis=1)).clip(lower=0)
    return reads, no_rank


def overview_table(long, totals):
    out = pd.DataFrame(index=totals.index)
    out["Unclassified"] = totals["unclassified"]
    for taxid, name in DOMAINS.items():
        out[name] = (long[long.taxid == taxid].groupby("sample")["clade_reads"].sum()
                     .reindex(totals.index, fill_value=0))
    out["Other classified"] = (totals["classified"] - out[list(DOMAINS.values())].sum(axis=1)).clip(lower=0)
    return out


def to_pct(df, totals):
    return df.div(totals["total"].replace(0, np.nan), axis=0).fillna(0) * 100


# ── Plot helpers ──────────────────────────────────────────────────────────────
def taxa_colours(n):
    # Greys are reserved for Other / no-rank / Unclassified, so drop them from the taxon palette
    grey_idx = {"tab20": {14, 15}, "tab20b": set(), "tab20c": {16, 17, 18, 19}}
    pool = [plt.get_cmap(c)(i) for c in ("tab20", "tab20b", "tab20c")
            for i in range(20) if i not in grey_idx[c]]
    return [pool[i % len(pool)] for i in range(n)]


def tick_fontsize(n):
    return 9 if n <= 40 else 7 if n <= 100 else 5 if n <= 250 else 3


def save(fig, outdir, stem):
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {stem}.png / .pdf")


def stacked_bar(pct_df, colours, title, ylabel, outdir, stem):
    """pct_df columns are plotted bottom-to-top; colours is a matching list."""
    n = len(pct_df)
    width = min(max(7, (0.22 if n <= 100 else 0.12) * n + 4), 60)
    fig, ax = plt.subplots(figsize=(width, 5.5))
    x = np.arange(n)
    bottom = np.zeros(n)
    for col, colour in zip(pct_df.columns, colours):
        vals = pct_df[col].to_numpy()
        ax.bar(x, vals, bottom=bottom, color=colour, width=0.85,
               edgecolor="white", linewidth=0.3 if n <= 100 else 0)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(pct_df.index, rotation=90, fontsize=tick_fontsize(n))
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(0, 100)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.spines[["top", "right"]].set_visible(False)
    handles = [Patch(facecolor=c, label=l) for l, c in zip(pct_df.columns, colours)][::-1]
    ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, frameon=False, ncol=1 if len(handles) <= 30 else 2)
    save(fig, outdir, stem)


def order_samples(pct_rank, how):
    if how == "name":
        return list(pct_rank.index)
    # "dominant": group samples by their top taxon, then by its abundance
    top = pct_rank.idxmax(axis=1)
    top_val = pct_rank.max(axis=1)
    taxon_rank = {t: i for i, t in enumerate(pct_rank.mean().sort_values(ascending=False).index)}
    key = pd.DataFrame({"t": top.map(taxon_rank), "v": -top_val})
    return list(key.sort_values(["t", "v"]).index)


def cluster_order(mat):
    """Hierarchical clustering order of rows (euclidean on log10 %, average linkage)."""
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
    except ImportError:
        return None
    if mat.shape[0] < 3:
        return None
    return leaves_list(linkage(np.log10(mat + 0.01), method="average", metric="euclidean"))


def heatmap(pct_rank, top_n, rank_label, cluster, outdir, stem):
    cols = pct_rank.mean().sort_values(ascending=False).index[:top_n]
    mat = pct_rank[cols]
    if cluster:
        r, c = cluster_order(mat.to_numpy()), cluster_order(mat.to_numpy().T)
        if r is not None:
            mat = mat.iloc[r]
        if c is not None:
            mat = mat.iloc[:, c]
        elif r is None:
            print("  (scipy not available or too few samples — heatmap not clustered)")
    n_s, n_t = mat.shape
    fig, ax = plt.subplots(figsize=(min(max(6, 0.32 * n_t + 3), 40),
                                    min(max(4, (0.2 if n_s <= 100 else 0.12) * n_s + 2), 100)))
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_under("white")
    vmax = max(100.0, float(mat.to_numpy().max()))
    im = ax.imshow(mat.to_numpy(), aspect="auto", interpolation="nearest", cmap=cmap,
                   norm=LogNorm(vmin=0.01, vmax=vmax))
    ax.set_xticks(range(n_t))
    ax.set_xticklabels(mat.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(n_s))
    ax.set_yticklabels(mat.index, fontsize=tick_fontsize(n_s))
    ax.set_title(f"Top {n_t} taxa at {rank_label} rank (% of total reads)")
    cb = fig.colorbar(im, ax=ax, shrink=min(1.0, 8 / max(n_s * 0.2, 1)) if n_s > 40 else 0.8,
                      extend="min", pad=0.02)
    cb.set_label("% of total reads (log scale; white < 0.01%)")
    save(fig, outdir, stem)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", required=True, help="Directory containing <sample>/kraken_report.txt")
    ap.add_argument("--outdir", required=True, help="Output directory for tables and plots")
    ap.add_argument("--rank", default="G",
                    help="Rank for stacked bar + heatmap: D,P,C,O,F,G,S or name, e.g. 'order' (default: G)")
    ap.add_argument("--top", type=int, default=20, help="Taxa shown in stacked bar; rest = Other (default: 20)")
    ap.add_argument("--heatmap-top", type=int, default=30, help="Taxa shown in heatmap (default: 30)")
    ap.add_argument("--sort", choices=["name", "dominant"], default="name",
                    help="Sample order in bar plots: by name, or grouped by dominant taxon (default: name)")
    ap.add_argument("--no-cluster", action="store_true", help="Do not cluster heatmap rows/columns")
    ap.add_argument("--report-name", default="kraken_report.txt", help="Report filename (default: kraken_report.txt)")
    args = ap.parse_args()

    rank = RANK_ALIASES.get(args.rank.lower(), args.rank.upper())
    if rank not in RANK_NAMES:
        sys.exit(f"ERROR: unsupported rank '{args.rank}'. Use one of: {', '.join(RANK_NAMES)}")
    rank_label = RANK_NAMES[rank]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    long = load_reports(args.indir, args.report_name)
    totals = sample_totals(long)
    empty = totals.index[totals["total"] == 0]
    if len(empty):
        print(f"WARNING: {len(empty)} sample(s) have 0 reads: {', '.join(empty)}", file=sys.stderr)
    print(f"Parsed {len(totals)} samples from {args.indir}")

    # Tables
    long.to_csv(outdir / "kraken_reports_long.tsv", sep="\t", index=False)
    reads, no_rank = rank_table(long, rank, totals)
    if reads.shape[1] == 0:
        sys.exit(f"ERROR: no taxa at rank '{rank}' in any report")
    pct_rank = to_pct(reads, totals)
    reads_out = reads.copy()
    reads_out[f"no_{rank_label}_assigned"] = no_rank
    reads_out["unclassified"] = totals["unclassified"]
    reads_out["total_reads"] = totals["total"]
    reads_out.to_csv(outdir / f"kraken_{rank_label}_reads.tsv", sep="\t")
    pct_out = to_pct(reads_out.drop(columns="total_reads"), totals).round(4)
    pct_out.to_csv(outdir / f"kraken_{rank_label}_pct.tsv", sep="\t")
    ov = overview_table(long, totals)
    ov_pct = to_pct(ov, totals)
    pd.concat([ov.add_suffix("_reads"), ov_pct.round(4).add_suffix("_pct"),
               totals["total"].rename("total_reads")], axis=1).to_csv(outdir / "kraken_overview.tsv", sep="\t")
    print(f"  wrote kraken_reports_long.tsv, kraken_{rank_label}_reads.tsv, "
          f"kraken_{rank_label}_pct.tsv, kraken_overview.tsv")

    order = order_samples(pct_rank, args.sort)

    # Stacked bar at rank
    top_taxa = list(pct_rank.mean().sort_values(ascending=False).index[:args.top])
    other = pct_rank.drop(columns=top_taxa).sum(axis=1)
    bar = pct_rank[top_taxa].copy()
    if (other > 0).any():
        bar["Other"] = other
    bar[f"Classified, no {rank_label} assigned"] = to_pct(no_rank.to_frame(), totals).iloc[:, 0]
    bar["Unclassified"] = to_pct(totals[["unclassified"]], totals)["unclassified"]
    colours = taxa_colours(len(top_taxa)) + ([COL_OTHER] if "Other" in bar else []) + [COL_NO_RANK, COL_UNCLASSIFIED]
    stacked_bar(bar.loc[order], colours,
                f"Kraken2 classification at {rank_label} rank (top {len(top_taxa)} taxa)",
                "% of total reads", outdir, f"stacked_bar_{rank_label}")

    # Heatmap
    heatmap(pct_rank, args.heatmap_top, rank_label, not args.no_cluster, outdir, f"heatmap_{rank_label}")

    # Classification overview (unclassified at top, like the stacked bar)
    ov_cols = ["Eukaryota", "Bacteria", "Archaea", "Viruses", "Other classified", "Unclassified"]
    ov_colours = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", COL_NO_RANK, COL_UNCLASSIFIED]
    stacked_bar(ov_pct.loc[order, ov_cols], ov_colours, "Kraken2 classification overview",
                "% of total reads", outdir, "classification_overview")


if __name__ == "__main__":
    main()
