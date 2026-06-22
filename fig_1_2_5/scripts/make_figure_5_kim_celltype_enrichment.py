from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib_cache")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/xdg_cache")
os.environ.setdefault("NUMBA_CACHE_DIR", "/private/tmp/numba_cache")
for key in ["MPLCONFIGDIR", "XDG_CACHE_HOME", "NUMBA_CACHE_DIR"]:
    Path(os.environ[key]).mkdir(parents=True, exist_ok=True)

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

import make_figure_5_kim_region_enrichment as base


THIS_DIR = Path(__file__).resolve().parent
ROOT = THIS_DIR.parent
OUT_DIR = ROOT / "outputs" / "figure_5_cross_species" / "kim_celltype_enrichment"
FIG_DIR = OUT_DIR / "figures"
CURRENT_H5AD = ROOT / "inputs" / "clean_objects" / "whole_habenula.h5ad"

MARKER_LOGFC_MIN = 0.25
MARKER_FDR_MAX = 0.05
MIN_CELLS_PER_CELLTYPE = 50
TOP_DOTPLOT_GENES_PER_CELLTYPE = 8

CELLTYPE_ORDER = [
    "Endothelial",
    "Pericytes",
    "Fibroblasts",
    "Astrocytes",
    "Microglia",
    "Macrophages",
    "Polydendrocytes",
    "Oligodendrocytes",
    "Diff. Oligodendrocytes",
    "LHb",
    "MHb",
]


def load_celltype_data() -> tuple[ad.AnnData, sparse.csr_matrix, list[str], list[str]]:
    a = ad.read_h5ad(CURRENT_H5AD)
    obs = a.obs.copy()
    obs["celltype"] = obs["gs"].astype(str)
    counts_by_ct = obs["celltype"].value_counts()
    categories = [ct for ct in CELLTYPE_ORDER if counts_by_ct.get(ct, 0) >= MIN_CELLS_PER_CELLTYPE]
    categories.extend(
        [
            ct
            for ct in sorted(counts_by_ct.index)
            if ct not in categories and counts_by_ct.get(ct, 0) >= MIN_CELLS_PER_CELLTYPE
        ]
    )
    keep = obs["celltype"].isin(categories).to_numpy()
    obs = obs.loc[keep].copy()
    obs["celltype"] = pd.Categorical(obs["celltype"], categories=categories, ordered=True)
    counts = sparse.csr_matrix(a.layers["counts"])[keep, :].tocsr()
    out = ad.AnnData(X=counts.copy(), obs=obs[["celltype"]].copy(), var=pd.DataFrame(index=a.var_names.astype(str)))
    out.obs_names = [f"cell_{i}" for i in range(out.n_obs)]
    out.layers["counts"] = counts.copy()
    sc.pp.normalize_total(out, target_sum=1e4)
    out.uns.pop("log1p", None)
    sc.pp.log1p(out)
    norm = sparse.csr_matrix(out.X)
    return out, norm, out.var_names.astype(str).tolist(), categories


def compute_celltype_markers(a: ad.AnnData, categories: list[str]) -> pd.DataFrame:
    sc.tl.rank_genes_groups(
        a,
        groupby="celltype",
        groups=categories,
        method="wilcoxon",
        reference="rest",
        pts=True,
        rankby_abs=False,
        key_added="celltype_markers",
    )
    pieces = []
    for celltype in categories:
        df = sc.get.rank_genes_groups_df(a, group=celltype, key="celltype_markers")
        df["celltype"] = celltype
        pieces.append(df)
    markers = pd.concat(pieces, ignore_index=True)
    markers = markers.rename(columns={"names": "mouse_gene", "logfoldchanges": "logFC_vs_other_celltypes"})
    markers["mouse_gene_upper"] = markers["mouse_gene"].astype(str).str.upper()
    markers["is_celltype_marker"] = markers["logFC_vs_other_celltypes"].gt(MARKER_LOGFC_MIN) & markers["pvals_adj"].lt(MARKER_FDR_MAX)
    markers.to_csv(OUT_DIR / "mouse_habenula_celltype_markers_wilcoxon.csv", index=False)
    markers.loc[markers["is_celltype_marker"]].to_csv(
        OUT_DIR / "mouse_habenula_celltype_markers_filtered_logfc0p25_fdr0p05.csv",
        index=False,
    )
    return markers


def average_expression_table(a: ad.AnnData, norm: sparse.csr_matrix, var_names: list[str], categories: list[str]) -> pd.DataFrame:
    rows = []
    counts = sparse.csr_matrix(a.layers["counts"])
    celltypes = a.obs["celltype"].astype(str).to_numpy()
    for celltype in categories:
        mask = celltypes == celltype
        mean_expr = np.asarray(norm[mask, :].mean(axis=0)).ravel()
        pct_expr = np.asarray((counts[mask, :] > 0).mean(axis=0)).ravel()
        rows.append(
            pd.DataFrame(
                {
                    "mouse_gene": var_names,
                    "celltype": celltype,
                    "mean_log1p_norm_expr": mean_expr,
                    "pct_expr": pct_expr,
                }
            )
        )
    expr = pd.concat(rows, ignore_index=True)
    expr.to_csv(OUT_DIR / "mouse_habenula_celltype_average_expression_long.csv", index=False)
    return expr


def map_kim_to_mouse(kim: pd.DataFrame, var_names: list[str]) -> pd.DataFrame:
    mapped = base.map_kim_to_mouse(kim, var_names)
    mapped.to_csv(OUT_DIR / "kim_human_degs_mapped_to_mouse_genes_celltype_universe.csv", index=False)
    return mapped


def categorize_kim_degs(mapped: pd.DataFrame, markers: pd.DataFrame, avg_expr_long: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    marker_sets = {
        celltype: set(markers.loc[markers["celltype"].eq(celltype) & markers["is_celltype_marker"], "mouse_gene_upper"])
        for celltype in markers["celltype"].drop_duplicates()
    }
    marker_info = markers.loc[markers["is_celltype_marker"]].copy()
    marker_info = marker_info[
        [
            "mouse_gene",
            "mouse_gene_upper",
            "celltype",
            "logFC_vs_other_celltypes",
            "pvals_adj",
            "pct_nz_group",
            "pct_nz_reference",
        ]
    ].rename(columns={"celltype": "marker_celltype", "pvals_adj": "celltype_marker_fdr"})

    rows = []
    for row in mapped.itertuples(index=False):
        hits = [celltype for celltype, genes in marker_sets.items() if row.mouse_gene_upper in genes]
        if hits:
            hit_info = marker_info[marker_info["mouse_gene_upper"].eq(row.mouse_gene_upper) & marker_info["marker_celltype"].isin(hits)]
            hit_info = hit_info.sort_values(["logFC_vs_other_celltypes", "celltype_marker_fdr"], ascending=[False, True])
            primary = str(hit_info.iloc[0]["marker_celltype"])
            all_hits = ";".join(hit_info["marker_celltype"].astype(str).tolist())
            primary_logfc = float(hit_info.iloc[0]["logFC_vs_other_celltypes"])
            primary_fdr = float(hit_info.iloc[0]["celltype_marker_fdr"])
        else:
            primary = "not_celltype_enriched"
            all_hits = ""
            primary_logfc = np.nan
            primary_fdr = np.nan
        rec = row._asdict()
        rec["primary_enriched_celltype"] = primary
        rec["all_enriched_celltypes"] = all_hits
        rec["primary_celltype_logFC"] = primary_logfc
        rec["primary_celltype_marker_fdr"] = primary_fdr
        rows.append(rec)

    cat = pd.DataFrame(rows)
    avg_wide = avg_expr_long.pivot(index="mouse_gene", columns="celltype", values="mean_log1p_norm_expr")
    avg_wide.columns = [f"avg_expr_{col}" for col in avg_wide.columns]
    pct_wide = avg_expr_long.pivot(index="mouse_gene", columns="celltype", values="pct_expr")
    pct_wide.columns = [f"pct_expr_{col}" for col in pct_wide.columns]
    expr_wide = pd.concat([avg_wide, pct_wide], axis=1).reset_index()
    cat = cat.merge(expr_wide, on="mouse_gene", how="left")
    cat = cat.sort_values(["primary_enriched_celltype", "p_value", "human_gene"])
    cat.to_csv(OUT_DIR / "kim_degs_categorized_by_mouse_habenula_celltype_enrichment.csv", index=False)

    unique = cat.drop_duplicates(["human_gene", "mouse_gene"])
    summary = (
        unique.groupby(["primary_enriched_celltype", "direction"], dropna=False)
        .size()
        .reset_index(name="n_mapped_deg_genes")
    )
    total = unique.shape[0]
    summary["percent_of_mapped_degs"] = summary["n_mapped_deg_genes"] / max(total, 1) * 100
    summary.to_csv(OUT_DIR / "kim_deg_celltype_enrichment_category_counts.csv", index=False)
    return cat, summary


def fisher_celltype_overlap(mapped: pd.DataFrame, markers: pd.DataFrame, var_names: list[str], categories: list[str]) -> pd.DataFrame:
    mapped_genes = set(mapped["mouse_gene_upper"])
    universe = set(pd.Series(var_names).str.upper())
    rows = []
    for celltype in categories:
        marker_set = set(markers.loc[markers["celltype"].eq(celltype) & markers["is_celltype_marker"], "mouse_gene_upper"]) & universe
        overlap = sorted(marker_set & mapped_genes)
        a_count = len(overlap)
        b_count = len(marker_set - mapped_genes)
        c_count = len(mapped_genes - marker_set)
        d_count = max(len(universe) - a_count - b_count - c_count, 0)
        odds_ratio, p_value = stats.fisher_exact([[a_count, b_count], [c_count, d_count]], alternative="greater")
        rows.append(
            {
                "celltype": celltype,
                "n_celltype_markers": len(marker_set),
                "n_kim_mapped_deg_genes": len(mapped_genes),
                "n_overlap": a_count,
                "odds_ratio": odds_ratio,
                "p_value": p_value,
                "overlap_mouse_genes": ";".join(overlap),
            }
        )
    out = pd.DataFrame(rows)
    out["p_adj"] = multipletests(out["p_value"], method="fdr_bh")[1]
    out["neg_log10_fdr"] = -np.log10(np.maximum(out["p_adj"], 1e-300))
    out = out.sort_values(["p_adj", "p_value", "n_overlap"], ascending=[True, True, False])
    out.to_csv(OUT_DIR / "kim_deg_x_mouse_habenula_celltype_marker_overlap_fisher.csv", index=False)
    return out


def plot_summary(summary: pd.DataFrame, overlap: pd.DataFrame, categories: list[str]) -> dict[str, Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")

    top_order = (
        summary.groupby("primary_enriched_celltype")["n_mapped_deg_genes"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    )
    top_order = [ct for ct in categories if ct in top_order] + [ct for ct in top_order if ct not in categories]
    if "not_celltype_enriched" in top_order:
        top_order = [ct for ct in top_order if ct != "not_celltype_enriched"] + ["not_celltype_enriched"]

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.4), gridspec_kw={"width_ratios": [1.45, 1.0]})
    sns.barplot(
        data=summary,
        x="primary_enriched_celltype",
        y="n_mapped_deg_genes",
        hue="direction",
        order=top_order,
        palette={"up_in_suicide": "#C44E52", "down_in_suicide": "#4C72B0"},
        ax=axes[0],
    )
    axes[0].set_xlabel("Primary mouse cell-type marker category")
    axes[0].set_ylabel("Kim mapped DEG gene count")
    axes[0].set_title("Kim suicide Hb DEGs categorized by mouse habenula cell-type enrichment")
    axes[0].tick_params(axis="x", labelrotation=35)

    heat = overlap.set_index("celltype")[["neg_log10_fdr"]].reindex(categories)
    annot = overlap.set_index("celltype")[["n_overlap"]].reindex(categories)
    vmax = max(1.0, float(np.nanmax(heat.to_numpy(dtype=float))))
    sns.heatmap(
        heat,
        cmap="mako",
        vmin=0,
        vmax=vmax,
        linewidths=0.35,
        linecolor="white",
        annot=annot,
        fmt=".0f",
        cbar_kws={"label": "-log10 FDR"},
        ax=axes[1],
    )
    axes[1].set_xlabel("Kim DEG marker-set overlap")
    axes[1].set_ylabel("Mouse habenula cell-type marker set")
    axes[1].set_title("Fisher enrichment of Kim DEGs in cell-type markers")

    fig.tight_layout()
    png = FIG_DIR / "kim_celltype_enrichment_summary.png"
    pdf = FIG_DIR / "kim_celltype_enrichment_summary.pdf"
    fig.savefig(png, dpi=240)
    fig.savefig(pdf)
    plt.close(fig)
    return {"summary_png": png, "summary_pdf": pdf}


def plot_dotplot(cat: pd.DataFrame, categories: list[str]) -> dict[str, Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    enriched = cat[cat["primary_enriched_celltype"].ne("not_celltype_enriched")].drop_duplicates(["human_gene", "mouse_gene"]).copy()
    if enriched.empty:
        enriched = cat.drop_duplicates(["human_gene", "mouse_gene"]).sort_values("p_value").head(50).copy()
    else:
        pieces = []
        ranked_celltypes = (
            enriched.groupby("primary_enriched_celltype")
            .size()
            .sort_values(ascending=False)
            .index.tolist()
        )
        for celltype in ranked_celltypes:
            pieces.append(enriched[enriched["primary_enriched_celltype"].eq(celltype)].sort_values("p_value").head(TOP_DOTPLOT_GENES_PER_CELLTYPE))
        enriched = pd.concat(pieces, ignore_index=True)

    enriched["label"] = enriched["human_gene"] + " / " + enriched["mouse_gene"]
    labels = enriched["label"].tolist()
    records = []
    for _, row in enriched.iterrows():
        for celltype in categories:
            records.append(
                {
                    "label": row["label"],
                    "human_gene": row["human_gene"],
                    "mouse_gene": row["mouse_gene"],
                    "primary_enriched_celltype": row["primary_enriched_celltype"],
                    "direction": row["direction"],
                    "celltype": celltype,
                    "mean_log1p_norm_expr": row[f"avg_expr_{celltype}"],
                    "pct_expr": row[f"pct_expr_{celltype}"] * 100,
                }
            )
    dot = pd.DataFrame(records)
    mean_by_gene = dot.groupby("label", observed=True)["mean_log1p_norm_expr"].transform("mean")
    sd_by_gene = dot.groupby("label", observed=True)["mean_log1p_norm_expr"].transform("std").replace(0, np.nan)
    dot["expr_z_by_gene"] = ((dot["mean_log1p_norm_expr"] - mean_by_gene) / sd_by_gene).fillna(0).clip(-2.5, 2.5)
    dot.to_csv(OUT_DIR / "kim_celltype_enriched_deg_dotplot_values.csv", index=False)

    x_pos = {celltype: i for i, celltype in enumerate(categories)}
    y_pos = {label: i for i, label in enumerate(labels)}
    fig, ax = plt.subplots(figsize=(10.8, max(7.5, 0.28 * len(labels) + 1.6)))
    scatter = ax.scatter(
        dot["celltype"].map(x_pos),
        dot["label"].map(y_pos),
        s=np.maximum(dot["pct_expr"], 1) * 2.5,
        c=dot["expr_z_by_gene"],
        cmap="vlag",
        vmin=-2.5,
        vmax=2.5,
        edgecolor="0.35",
        linewidth=0.2,
    )
    ax.set_xticks(range(len(categories)), categories, rotation=35, ha="right")
    ax.set_yticks(range(len(labels)), labels)
    ax.invert_yaxis()
    ax.set_xlabel("Mouse habenula cell type")
    ax.set_ylabel("Kim DEG / mapped mouse gene")
    ax.set_title("Kim DEGs with mouse habenula cell-type-enriched expression")
    cbar = fig.colorbar(scatter, ax=ax, pad=0.015)
    cbar.set_label("Mean expression z-score across cell types")
    for size in [10, 30, 60]:
        ax.scatter([], [], s=size * 2.5, c="lightgray", edgecolor="0.35", label=f"{size}%")
    ax.legend(title="Pct expressed", bbox_to_anchor=(1.08, 1.0), loc="upper left", frameon=False)
    fig.tight_layout()
    png = FIG_DIR / "kim_celltype_enriched_deg_dotplot.png"
    pdf = FIG_DIR / "kim_celltype_enriched_deg_dotplot.pdf"
    fig.savefig(png, dpi=240)
    fig.savefig(pdf)
    plt.close(fig)
    return {"dotplot_png": png, "dotplot_pdf": pdf}


def write_summary(
    kim: pd.DataFrame,
    mapped: pd.DataFrame,
    markers: pd.DataFrame,
    cat: pd.DataFrame,
    summary: pd.DataFrame,
    overlap: pd.DataFrame,
    categories: list[str],
) -> Path:
    unique = cat.drop_duplicates(["human_gene", "mouse_gene"])
    n_mapped = mapped["mouse_gene"].nunique()
    n_enriched = unique[unique["primary_enriched_celltype"].ne("not_celltype_enriched")].shape[0]
    marker_counts = markers[markers["is_celltype_marker"]].groupby("celltype").size().reindex(categories).fillna(0).astype(int)
    top_overlap = overlap.sort_values(["p_adj", "p_value", "n_overlap"], ascending=[True, True, False]).head(5)

    lines = [
        "# Kim 2022 Cell-Type Enrichment Replication",
        "",
        "This analysis repeats the Kim et al. cell-type enrichment concept using this project's broad `gs` cell-type annotations from `inputs/clean_objects/whole_habenula.h5ad`.",
        "",
        "## Method",
        "",
        "- Kim et al. Table S2 human suicide-vs-control Hb DEGs were parsed from Additional File 2.",
        "- Human symbols were mapped to mouse genes using the local mouse-human symbol link table, with case-insensitive fallback.",
        "- Mouse cell-type markers were computed from raw counts in `inputs/clean_objects/whole_habenula.h5ad` using Wilcoxon rank-sum tests via Scanpy.",
        f"- Marker threshold: `logFC > {MARKER_LOGFC_MIN}` and `adjusted p < {MARKER_FDR_MAX}`.",
        "- Each mapped Kim DEG was categorized by its strongest positive cell-type marker assignment, if any.",
        "- Fisher's exact tests evaluated overlap between the Kim mapped DEG set and each mouse cell-type marker set.",
        "",
        "## Key Counts",
        "",
        f"- Unique expanded Kim human DEG symbols: {kim['human_gene'].nunique()}",
        f"- Unique mapped mouse Kim DEG genes present in the whole h5ad object: {n_mapped}",
        f"- Kim mapped DEG genes with at least one cell-type-enriched marker assignment: {n_enriched}",
        "",
        "## Cell-Type Marker Counts",
        "",
        "| celltype | n_marker_genes |",
        "| --- | ---: |",
    ]
    for celltype, count in marker_counts.items():
        lines.append(f"| {celltype} | {count} |")

    lines.extend(
        [
            "",
            "## Top Fisher Overlap Results",
            "",
            "| celltype | n_markers | n_overlap | odds_ratio | FDR | overlap_mouse_genes |",
            "| --- | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in top_overlap.itertuples(index=False):
        genes = row.overlap_mouse_genes if isinstance(row.overlap_mouse_genes, str) else ""
        if len(genes) > 180:
            genes = genes[:177] + "..."
        lines.append(
            f"| {row.celltype} | {row.n_celltype_markers} | {row.n_overlap} | {row.odds_ratio:.3g} | {row.p_adj:.3g} | {genes} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is the faithful companion to the Kim-style analysis because it tests cell-type marker enrichment rather than only broad LHb/MHb region enrichment. If endothelial/vascular markers dominate this table, it supports the paper's interpretation in this local mouse single-cell dataset. If not, it suggests the Kim endothelial signal is not recovered under these annotations and thresholds.",
        ]
    )
    return OUT_DIR


def main() -> dict[str, Path | int]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    base.OUT_DIR.mkdir(parents=True, exist_ok=True)

    kim = base.load_kim_degs()
    adata, norm, var_names, categories = load_celltype_data()
    markers = compute_celltype_markers(adata, categories)
    avg_expr = average_expression_table(adata, norm, var_names, categories)
    mapped = map_kim_to_mouse(kim, var_names)
    cat, summary = categorize_kim_degs(mapped, markers, avg_expr)
    overlap = fisher_celltype_overlap(mapped, markers, var_names, categories)
    fig_paths = plot_summary(summary, overlap, categories)
    dot_paths = plot_dotplot(cat, categories)
    summary_path = write_summary(kim, mapped, markers, cat, summary, overlap, categories)

    print(f"Wrote {OUT_DIR}")
    print(f"Mapped Kim DEG mouse genes: {mapped['mouse_gene'].nunique()}")
    print(f"Filtered cell-type markers: {int(markers['is_celltype_marker'].sum())}")
    print(overlap.head(12).to_string(index=False))
    return {
        "output_dir": OUT_DIR,
        "summary": summary_path,
        **fig_paths,
        **dot_paths,
    }


if __name__ == "__main__":
    main()
