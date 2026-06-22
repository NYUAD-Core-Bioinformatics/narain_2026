from __future__ import annotations

import os
from itertools import product
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib_cache")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/xdg_cache")
for key in ["MPLCONFIGDIR", "XDG_CACHE_HOME"]:
    Path(os.environ[key]).mkdir(parents=True, exist_ok=True)

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


HUMAN_MARKER_TABLE = Path("inputs/figure_5_cross_species/yalcinbas_single_nucleus_marker_genes.csv")
LHB_H5AD = Path("inputs/clean_objects/lateral_habenula_regions.h5ad")
MHB_H5AD = Path("inputs/clean_objects/medial_habenula_regions.h5ad")
OUT_DIR = Path("outputs/figure_5_cross_species/yalcinbas_sensitivity_analysis")

HUMAN_HB_POPULATIONS = ["LHb.1", "LHb.2", "LHb.3", "LHb.4", "LHb.5", "LHb.6", "LHb.7", "MHb.1", "MHb.2", "MHb.3"]
MODULE_DEPTHS = [10, 15, 25, 50, 100]
HUMAN_OVERLAP_DEPTHS = [25, 50, 100, 200]
MOUSE_OVERLAP_DEPTHS = [50, 100, 200, 500]

HUMAN_CLUSTER_PALETTE = {
    "LHb.1": "#3B7EA1",
    "LHb.2": "#D95F02",
    "LHb.3": "#5AA469",
    "LHb.4": "#8E5EA2",
    "LHb.5": "#C44E52",
    "LHb.6": "#E6AB02",
    "LHb.7": "#4C6A9C",
    "MHb.1": "#2A9D8F",
    "MHb.2": "#E76F51",
    "MHb.3": "#7F7F7F",
}

LHB_REGION_ORDER = ["Hbx", "Lateral", "Marginal", "Oval/Medial"]
MHB_REGION_MAP = {
    "0": "Lateral",
    "1": "Ventral",
    "2": "Dorsal",
    "3": "Ventral/Dorsal",
    "4": "Superior/Dorsal",
    "5": "Ventral/Lateral/Dorsal",
}
MHB_REGION_ORDER = ["Ventral", "Ventral/Dorsal", "Lateral", "Dorsal", "Superior/Dorsal", "Ventral/Lateral/Dorsal"]

BAD_GENE_PREFIXES = ("Rpl", "Rps", "mt-", "Hbb", "Hba")


def natural_key(value: object) -> list[object]:
    text = str(value)
    out: list[object] = []
    current = ""
    is_digit = False
    for char in text:
        if char.isdigit() == is_digit:
            current += char
        else:
            if current:
                out.append(int(current) if is_digit else current)
            current = char
            is_digit = char.isdigit()
    if current:
        out.append(int(current) if is_digit else current)
    return out


def as_csr(x) -> sparse.csr_matrix:
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(x)


def normalize_log1p_counts(x: sparse.csr_matrix, target_sum: float = 1e4) -> sparse.csr_matrix:
    x = x.astype(np.float32, copy=True).tocsr()
    lib = np.asarray(x.sum(axis=1)).ravel()
    scale = np.zeros_like(lib, dtype=np.float32)
    keep = lib > 0
    scale[keep] = target_sum / lib[keep]
    x = x.multiply(scale[:, None]).tocsr()
    x.data = np.log1p(x.data)
    return x


def aggregate_by_group(x: sparse.csr_matrix, groups: pd.Series, categories: list[str], binary: bool = False) -> np.ndarray:
    if binary:
        x = x.copy()
        x.data = np.ones_like(x.data)
    values = groups.astype(str).to_numpy()
    rows = []
    for category in categories:
        mask = values == category
        rows.append(np.asarray(x[mask, :].mean(axis=0)).ravel())
    return np.vstack(rows)


def zscore_dense(x: np.ndarray) -> np.ndarray:
    return (x - x.mean(axis=0, keepdims=True)) / (x.std(axis=0, keepdims=True) + 1e-6)


def gene_indices(var_names: list[str], genes: list[str]) -> dict[str, int]:
    lookup = {gene: i for i, gene in enumerate(var_names)}
    return {gene: lookup[gene] for gene in genes if gene in lookup}


def is_excluded_gene(gene: str) -> bool:
    return gene.startswith(BAD_GENE_PREFIXES) or gene in {"Malat1", "Meg3"}


def load_human_markers() -> pd.DataFrame:
    markers = pd.read_csv(HUMAN_MARKER_TABLE)
    symbol = pd.Series(index=markers.index, dtype=object)
    for candidate in ["Symbol", "ID.x", "ID", "gene_id"]:
        if candidate in markers.columns:
            values = markers[candidate].astype(object)
            values = values.where(values.notna() & values.astype(str).ne("nan") & values.astype(str).ne(""))
            symbol = symbol.fillna(values)
    markers = markers.copy()
    markers["human_gene"] = symbol.astype(str)
    markers = markers[markers["human_gene"].notna() & markers["human_gene"].ne("nan")].copy()
    markers["gene_upper"] = markers["human_gene"].str.upper()
    markers = markers.sort_values(["human_cluster", "adj.P.Val", "P.Value", "logFC"], ascending=[True, True, True, False])
    return markers


def map_human_to_mouse(markers: pd.DataFrame, var_names: list[str]) -> pd.DataFrame:
    mouse_lookup = pd.DataFrame({"mouse_gene": var_names})
    mouse_lookup["gene_upper"] = mouse_lookup["mouse_gene"].str.upper()
    mouse_lookup = mouse_lookup.drop_duplicates("gene_upper")
    mapped = markers.merge(mouse_lookup, on="gene_upper", how="inner")
    mapped = mapped.sort_values(["human_cluster", "adj.P.Val", "P.Value", "logFC"], ascending=[True, True, True, False])
    return mapped


def load_compartment(compartment: str) -> tuple[pd.DataFrame, sparse.csr_matrix, sparse.csr_matrix, list[str], list[str]]:
    if compartment == "LHb":
        a = ad.read_h5ad(LHB_H5AD)
        obs = a.obs.copy()
        obs["mouse_region"] = obs["cl"].astype(str).replace({"HbX": "Hbx", "Oval-Medial": "Oval/Medial"})
        categories = [x for x in LHB_REGION_ORDER if x in set(obs["mouse_region"])]
    elif compartment == "MHb":
        a = ad.read_h5ad(MHB_H5AD)
        obs = a.obs.copy()
        obs["mouse_region"] = obs["leiden"].astype(str).map(MHB_REGION_MAP).fillna("Unknown")
        categories = [x for x in MHB_REGION_ORDER if x in set(obs["mouse_region"])]
    else:
        raise ValueError(f"Unknown compartment: {compartment}")

    obs = obs[obs["mouse_region"].isin(categories)].copy()
    keep = np.asarray(a.obs_names.isin(obs.index))
    counts = as_csr(a.layers["counts"])[keep, :].tocsr()
    norm = normalize_log1p_counts(counts)
    var_names = a.var_names.astype(str).tolist()
    obs["mouse_region"] = pd.Categorical(obs["mouse_region"], categories=categories, ordered=True)
    return obs, counts, norm, var_names, categories


def compute_mouse_gene_scores(
    compartment: str,
    obs: pd.DataFrame,
    counts: sparse.csr_matrix,
    norm: sparse.csr_matrix,
    var_names: list[str],
    categories: list[str],
) -> pd.DataFrame:
    group_sizes = obs["mouse_region"].value_counts().reindex(categories).to_numpy()
    mean_by_group = aggregate_by_group(norm, obs["mouse_region"], categories)
    pct_by_group = aggregate_by_group(counts, obs["mouse_region"], categories, binary=True)
    total_sum = np.asarray(norm.sum(axis=0)).ravel()
    total_n = norm.shape[0]

    records = []
    for i, region in enumerate(categories):
        rest_n = total_n - group_sizes[i]
        rest_mean = (total_sum - mean_by_group[i] * group_sizes[i]) / max(rest_n, 1)
        df = pd.DataFrame(
            {
                "compartment": compartment,
                "mouse_region": region,
                "mouse_gene": var_names,
                "gene_upper": pd.Series(var_names).str.upper().to_numpy(),
                "mean_log_expr": mean_by_group[i],
                "pct_expr": pct_by_group[i],
                "mouse_logFC_vs_other_regions": mean_by_group[i] - rest_mean,
            }
        )
        df["excluded_gene"] = df["mouse_gene"].map(is_excluded_gene)
        records.append(df)
    scores = pd.concat(records, ignore_index=True)
    scores.to_csv(OUT_DIR / f"{compartment.lower()}_mouse_region_gene_scores_all_genes.csv", index=False)
    return scores


def top_human_set(mapped: pd.DataFrame, human_cluster: str, depth: int, universe: set[str]) -> set[str]:
    sub = mapped[mapped["human_cluster"].eq(human_cluster) & mapped["logFC"].gt(0)].drop_duplicates("gene_upper")
    return set(sub.head(depth)["gene_upper"]) & universe


def top_mouse_set(scores: pd.DataFrame, region: str, depth: int, universe: set[str]) -> set[str]:
    sub = scores[
        scores["mouse_region"].astype(str).eq(region)
        & scores["mouse_logFC_vs_other_regions"].gt(0)
        & scores["pct_expr"].ge(0.05)
        & ~scores["excluded_gene"]
    ].sort_values("mouse_logFC_vs_other_regions", ascending=False)
    return set(sub.head(depth)["gene_upper"]) & universe


def run_overlap_sensitivity(
    compartment: str,
    mapped: pd.DataFrame,
    scores: pd.DataFrame,
    categories: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_universe = set(scores.loc[~scores["excluded_gene"], "gene_upper"])
    mapped_universe = set(mapped["gene_upper"])
    universe = score_universe & mapped_universe

    rows = []
    for human_depth, mouse_depth in product(HUMAN_OVERLAP_DEPTHS, MOUSE_OVERLAP_DEPTHS):
        grid_rows = []
        for region in categories:
            mouse_set = top_mouse_set(scores, region, mouse_depth, universe)
            for human_cluster in HUMAN_HB_POPULATIONS:
                human_set = top_human_set(mapped, human_cluster, human_depth, universe)
                overlap = sorted(mouse_set & human_set)
                a_count = len(overlap)
                b_count = len(mouse_set - human_set)
                c_count = len(human_set - mouse_set)
                d_count = max(len(universe) - a_count - b_count - c_count, 0)
                odds, p_value = stats.fisher_exact([[a_count, b_count], [c_count, d_count]], alternative="greater")
                grid_rows.append(
                    {
                        "compartment": compartment,
                        "mouse_region": region,
                        "human_cluster": human_cluster,
                        "human_marker_depth": human_depth,
                        "mouse_marker_depth": mouse_depth,
                        "n_mouse_markers": len(mouse_set),
                        "n_human_markers": len(human_set),
                        "n_overlap": a_count,
                        "odds_ratio": odds,
                        "p_value": p_value,
                        "overlap_genes": ";".join(overlap),
                    }
                )
        grid_df = pd.DataFrame(grid_rows)
        grid_df["p_adj_grid"] = multipletests(grid_df["p_value"], method="fdr_bh")[1]
        rows.append(grid_df)

    sensitivity = pd.concat(rows, ignore_index=True)
    sensitivity["p_adj_all_sensitivity_tests"] = multipletests(sensitivity["p_value"], method="fdr_bh")[1]
    sensitivity["neg_log10_p_adj_grid"] = -np.log10(np.maximum(sensitivity["p_adj_grid"], 1e-300))
    sensitivity.to_csv(OUT_DIR / f"{compartment.lower()}_overlap_threshold_sensitivity.csv", index=False)

    best = (
        sensitivity.sort_values(["compartment", "mouse_region", "human_marker_depth", "mouse_marker_depth", "p_adj_grid", "p_value", "n_overlap"], ascending=[True, True, True, True, True, True, False])
        .groupby(["compartment", "mouse_region", "human_marker_depth", "mouse_marker_depth"], as_index=False)
        .head(1)
    )
    n_grid = len(HUMAN_OVERLAP_DEPTHS) * len(MOUSE_OVERLAP_DEPTHS)
    best_freq = (
        best.groupby(["compartment", "mouse_region", "human_cluster"], observed=True)
        .size()
        .reset_index(name="n_grid_best_match")
    )
    best_freq["best_match_fraction"] = best_freq["n_grid_best_match"] / n_grid

    sig_freq = (
        sensitivity.assign(significant=sensitivity["p_adj_grid"].lt(0.05))
        .groupby(["compartment", "mouse_region", "human_cluster"], observed=True)
        .agg(
            n_grid_significant=("significant", "sum"),
            median_neg_log10_fdr=("neg_log10_p_adj_grid", "median"),
            min_fdr=("p_adj_grid", "min"),
            max_overlap=("n_overlap", "max"),
        )
        .reset_index()
    )
    sig_freq["significant_fraction"] = sig_freq["n_grid_significant"] / n_grid
    summary = sig_freq.merge(best_freq, on=["compartment", "mouse_region", "human_cluster"], how="left")
    summary["n_grid_best_match"] = summary["n_grid_best_match"].fillna(0).astype(int)
    summary["best_match_fraction"] = summary["best_match_fraction"].fillna(0.0)
    summary = summary.sort_values(["compartment", "mouse_region", "best_match_fraction", "significant_fraction", "median_neg_log10_fdr"], ascending=[True, True, False, False, False])
    summary.to_csv(OUT_DIR / f"{compartment.lower()}_overlap_threshold_sensitivity_summary.csv", index=False)
    return sensitivity, summary


def run_module_depth_sensitivity(
    compartment: str,
    mapped: pd.DataFrame,
    obs: pd.DataFrame,
    norm: sparse.csr_matrix,
    var_names: list[str],
    categories: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for depth in MODULE_DEPTHS:
        human_sets: dict[str, list[str]] = {}
        union_genes: list[str] = []
        for human_cluster in HUMAN_HB_POPULATIONS:
            sub = mapped[mapped["human_cluster"].eq(human_cluster) & mapped["logFC"].gt(0)].drop_duplicates("gene_upper")
            genes = sub.head(depth)["mouse_gene"].tolist()
            human_sets[human_cluster] = genes
            union_genes.extend(genes)

        union_genes = sorted(set(union_genes))
        idx = gene_indices(var_names, union_genes)
        union_genes = [gene for gene in union_genes if gene in idx]
        gene_position = {gene: i for i, gene in enumerate(union_genes)}
        if not union_genes:
            continue

        dense = norm[:, [idx[gene] for gene in union_genes]].toarray().astype(np.float32)
        dense = zscore_dense(dense)

        for human_cluster in HUMAN_HB_POPULATIONS:
            genes = [gene for gene in human_sets[human_cluster] if gene in gene_position]
            if not genes:
                continue
            score = dense[:, [gene_position[gene] for gene in genes]].mean(axis=1)
            for region in categories:
                mask = obs["mouse_region"].astype(str).eq(region).to_numpy()
                rows.append(
                    {
                        "compartment": compartment,
                        "mouse_region": region,
                        "human_cluster": human_cluster,
                        "module_marker_depth": depth,
                        "n_marker_genes_mapped": len(genes),
                        "mean_module_score": float(np.mean(score[mask])),
                    }
                )

    module = pd.DataFrame(rows)
    module.to_csv(OUT_DIR / f"{compartment.lower()}_module_depth_sensitivity.csv", index=False)

    best = (
        module.sort_values(["compartment", "mouse_region", "module_marker_depth", "mean_module_score"], ascending=[True, True, True, False])
        .groupby(["compartment", "mouse_region", "module_marker_depth"], as_index=False)
        .head(1)
    )
    n_depth = len(MODULE_DEPTHS)
    summary = (
        best.groupby(["compartment", "mouse_region", "human_cluster"], observed=True)
        .size()
        .reset_index(name="n_depth_best_module")
    )
    summary["best_module_fraction"] = summary["n_depth_best_module"] / n_depth
    mean_scores = (
        module.groupby(["compartment", "mouse_region", "human_cluster"], observed=True)["mean_module_score"]
        .mean()
        .reset_index(name="mean_score_across_depths")
    )
    summary = summary.merge(mean_scores, on=["compartment", "mouse_region", "human_cluster"], how="right")
    summary["n_depth_best_module"] = summary["n_depth_best_module"].fillna(0).astype(int)
    summary["best_module_fraction"] = summary["best_module_fraction"].fillna(0.0)
    summary = summary.sort_values(["compartment", "mouse_region", "best_module_fraction", "mean_score_across_depths"], ascending=[True, True, False, False])
    summary.to_csv(OUT_DIR / f"{compartment.lower()}_module_depth_sensitivity_summary.csv", index=False)
    return module, summary


def run_rank_correlation(
    compartment: str,
    mapped: pd.DataFrame,
    scores: pd.DataFrame,
    categories: list[str],
) -> pd.DataFrame:
    rows = []
    score_base = scores.loc[~scores["excluded_gene"], ["mouse_region", "gene_upper", "mouse_gene", "mouse_logFC_vs_other_regions"]].copy()
    for human_cluster in HUMAN_HB_POPULATIONS:
        human = mapped[mapped["human_cluster"].eq(human_cluster)].copy()
        human = human.sort_values(["adj.P.Val", "P.Value", "logFC"], ascending=[True, True, False]).drop_duplicates("gene_upper")
        human = human[["gene_upper", "human_gene", "logFC", "P.Value", "adj.P.Val"]].rename(columns={"logFC": "human_logFC"})
        for region in categories:
            mouse = score_base[score_base["mouse_region"].astype(str).eq(region)].copy()
            merged = mouse.merge(human, on="gene_upper", how="inner")
            if merged.shape[0] < 20 or merged["mouse_logFC_vs_other_regions"].nunique() < 2 or merged["human_logFC"].nunique() < 2:
                rho, p_value = np.nan, 1.0
            else:
                rho, p_value = stats.spearmanr(merged["mouse_logFC_vs_other_regions"], merged["human_logFC"])
                if not np.isfinite(rho) or not np.isfinite(p_value):
                    rho, p_value = np.nan, 1.0
            rows.append(
                {
                    "compartment": compartment,
                    "mouse_region": region,
                    "human_cluster": human_cluster,
                    "n_ranked_genes": int(merged.shape[0]),
                    "spearman_r": rho,
                    "p_value": p_value,
                }
            )
    corr = pd.DataFrame(rows)
    corr["p_adj"] = multipletests(corr["p_value"], method="fdr_bh")[1]
    corr["neg_log10_p_adj"] = -np.log10(np.maximum(corr["p_adj"], 1e-300))
    corr = corr.sort_values(["compartment", "mouse_region", "p_adj", "spearman_r"], ascending=[True, True, True, False])
    corr.to_csv(OUT_DIR / f"{compartment.lower()}_threshold_free_rank_correlation.csv", index=False)
    return corr


def plot_pair_heatmap(
    df: pd.DataFrame,
    compartment: str,
    value_col: str,
    output_name: str,
    title: str,
    center: float | None = None,
    cmap: str = "mako",
    annot_col: str | None = None,
    human_clusters: list[str] | None = None,
) -> None:
    heat = df.pivot(index="mouse_region", columns="human_cluster", values=value_col)
    region_order = LHB_REGION_ORDER if compartment == "LHb" else MHB_REGION_ORDER
    human_clusters = human_clusters or HUMAN_HB_POPULATIONS
    heat = heat.reindex([x for x in region_order if x in heat.index], columns=human_clusters)
    annot = None
    if annot_col:
        annot = df.pivot(index="mouse_region", columns="human_cluster", values=annot_col)
        annot = annot.reindex(index=heat.index, columns=heat.columns)
    fig, ax = plt.subplots(figsize=(11, max(3.2, 0.5 * heat.shape[0] + 1.5)))
    sns.heatmap(
        heat,
        cmap=cmap,
        center=center,
        linewidths=0.3,
        linecolor="white",
        annot=annot,
        fmt=".2g" if annot_col else "",
        cbar_kws={"label": value_col},
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xlabel("Human Yalcinbas Hb population")
    ax.set_ylabel(f"Mouse {compartment} recovered region")
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"{output_name}.pdf")
    fig.savefig(OUT_DIR / f"{output_name}.png", dpi=220)
    plt.close(fig)


def ordered_summary_rows(summary: pd.DataFrame) -> pd.DataFrame:
    order = [("LHb", region) for region in LHB_REGION_ORDER]
    order.extend(("MHb", region) for region in MHB_REGION_ORDER)
    keyed = summary.set_index(["compartment", "mouse_region"], drop=False)
    present = [key for key in order if key in keyed.index]
    out = keyed.loc[present].reset_index(drop=True)
    out["row_label"] = out["compartment"] + " " + out["mouse_region"]
    return out


def method_agreement(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    method_cols = [
        ("overlap", "best_overlap_human_cluster"),
        ("module", "best_module_human_cluster"),
        ("rankcorr", "best_positive_rankcorr_human_cluster"),
    ]
    for row in summary.itertuples(index=False):
        for method, col in method_cols:
            rows.append(
                {
                    "compartment": row.compartment,
                    "mouse_region": row.mouse_region,
                    "row_label": f"{row.compartment} {row.mouse_region}",
                    "human_cluster": getattr(row, col),
                    "method": method,
                    "agreement_count": 1,
                }
            )
    agreement = pd.DataFrame(rows)
    if agreement.empty:
        return agreement
    agreement = (
        agreement.groupby(["compartment", "mouse_region", "row_label", "human_cluster"], observed=True)["agreement_count"]
        .sum()
        .reset_index()
    )
    return agreement


def plot_calibrated_best_match_summary(summary: pd.DataFrame, output_name: str, title: str) -> None:
    summary = ordered_summary_rows(summary)
    n_rows = summary.shape[0]
    y_positions = np.arange(n_rows)
    row_labels = summary["row_label"].tolist()
    height = max(5.2, 0.55 * n_rows + 1.7)
    fig, axes = plt.subplots(1, 3, figsize=(15.5, height), sharey=True)

    panels = [
        (
            "Overlap grid",
            "overlap_best_fraction",
            "best_overlap_human_cluster",
            "fraction of threshold grid",
            0.0,
            1.05,
        ),
        (
            "Module score depths",
            "module_best_fraction",
            "best_module_human_cluster",
            "fraction of module depths",
            0.0,
            1.05,
        ),
        (
            "Rank correlation",
            "positive_rankcorr_spearman_r",
            "best_positive_rankcorr_human_cluster",
            "Spearman r",
            0.0,
            max(0.35, float(summary["positive_rankcorr_spearman_r"].max()) * 1.15),
        ),
    ]

    for ax, (panel_title, value_col, cluster_col, x_label, x_min, x_max) in zip(axes, panels):
        ax.axvline(0, color="#333333", linewidth=0.8)
        for y, row in zip(y_positions, summary.itertuples(index=False)):
            value = float(getattr(row, value_col))
            cluster = getattr(row, cluster_col)
            color = HUMAN_CLUSTER_PALETTE.get(cluster, "#666666")
            ax.scatter(value, y, s=135, color=color, edgecolor="white", linewidth=1.1, zorder=3)
            if value > (x_max * 0.78):
                text_x = value - (x_max * 0.035)
                ha = "right"
            else:
                text_x = value + (x_max * 0.035)
                ha = "left"
            ax.text(text_x, y, cluster, va="center", ha=ha, fontsize=9, color="#222222")

        ax.set_title(panel_title, fontsize=12, pad=10)
        ax.set_xlabel(x_label)
        ax.set_xlim(x_min, x_max)
        ax.grid(axis="x", color="#E0E0E0", linewidth=0.8)
        ax.set_axisbelow(True)
        ax.invert_yaxis()
        ax.tick_params(axis="y", length=0)

    axes[0].set_yticks(y_positions)
    axes[0].set_yticklabels(row_labels)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)

    compartments = summary["compartment"].tolist()
    for i in range(1, n_rows):
        if compartments[i] != compartments[i - 1]:
            for ax in axes:
                ax.axhline(i - 0.5, color="#444444", linewidth=1.2)

    used_clusters = sorted(
        {
            cluster
            for col in [
                "best_overlap_human_cluster",
                "best_module_human_cluster",
                "best_positive_rankcorr_human_cluster",
            ]
            for cluster in summary[col].dropna().astype(str)
        },
        key=natural_key,
    )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=7,
            markerfacecolor=HUMAN_CLUSTER_PALETTE.get(cluster, "#666666"),
            markeredgecolor="white",
            label=cluster,
        )
        for cluster in used_clusters
    ]
    fig.legend(handles=handles, loc="lower center", ncol=min(8, len(handles)), frameon=False, bbox_to_anchor=(0.5, 0.01))
    fig.suptitle(title, fontsize=14, y=0.98)
    fig.tight_layout(rect=[0, 0.06, 1, 0.95])
    fig.savefig(OUT_DIR / f"{output_name}.pdf")
    fig.savefig(OUT_DIR / f"{output_name}.png", dpi=240)
    plt.close(fig)


def plot_method_agreement_heatmap(summary: pd.DataFrame, output_name: str, title: str) -> None:
    summary = ordered_summary_rows(summary)
    agreement = method_agreement(summary)
    agreement.to_csv(OUT_DIR / f"{output_name}_values.csv", index=False)
    row_order = summary["row_label"].tolist()
    heat = agreement.pivot(index="row_label", columns="human_cluster", values="agreement_count")
    columns = [cluster for cluster in HUMAN_HB_POPULATIONS if cluster in heat.columns]
    heat = heat.reindex(index=row_order, columns=columns).fillna(0)
    annot = heat.astype(int).astype(str).replace("0", "")

    fig, ax = plt.subplots(figsize=(10.8, max(5.0, 0.5 * heat.shape[0] + 1.8)))
    sns.heatmap(
        heat,
        cmap=sns.color_palette(["#F4F4F2", "#BFD7EA", "#67A9CF", "#1B6CA8"], as_cmap=True),
        vmin=0,
        vmax=3,
        linewidths=0.4,
        linecolor="white",
        annot=annot,
        fmt="",
        cbar_kws={"label": "number of methods selecting population", "ticks": [0, 1, 2, 3]},
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xlabel("Human Yalcinbas Hb population")
    ax.set_ylabel("Mouse recovered region")

    labels = summary["compartment"].tolist()
    for i in range(1, len(labels)):
        if labels[i] != labels[i - 1]:
            ax.axhline(i, color="#444444", linewidth=1.2)

    fig.tight_layout()
    fig.savefig(OUT_DIR / f"{output_name}.pdf")
    fig.savefig(OUT_DIR / f"{output_name}.png", dpi=240)
    plt.close(fig)


def write_compartment_summary(
    compartment: str,
    overlap_summary: pd.DataFrame,
    module_summary: pd.DataFrame,
    rank_corr: pd.DataFrame,
) -> pd.DataFrame:
    best_overlap = (
        overlap_summary.sort_values(["mouse_region", "best_match_fraction", "significant_fraction", "median_neg_log10_fdr"], ascending=[True, False, False, False])
        .groupby(["compartment", "mouse_region"], as_index=False)
        .head(1)
        .rename(
            columns={
                "human_cluster": "best_overlap_human_cluster",
                "best_match_fraction": "overlap_best_fraction",
                "significant_fraction": "overlap_significant_fraction",
                "min_fdr": "overlap_min_fdr",
                "max_overlap": "overlap_max_genes",
            }
        )
    )
    best_module = (
        module_summary.sort_values(["mouse_region", "best_module_fraction", "mean_score_across_depths"], ascending=[True, False, False])
        .groupby(["compartment", "mouse_region"], as_index=False)
        .head(1)
        .rename(
            columns={
                "human_cluster": "best_module_human_cluster",
                "best_module_fraction": "module_best_fraction",
                "mean_score_across_depths": "module_mean_score_across_depths",
            }
        )
    )
    positive_rank = rank_corr[rank_corr["spearman_r"].gt(0)].copy()
    if positive_rank.empty:
        positive_rank = rank_corr.copy()
    best_rank = (
        positive_rank.sort_values(["mouse_region", "spearman_r", "p_adj"], ascending=[True, False, True])
        .groupby(["compartment", "mouse_region"], as_index=False)
        .head(1)
        .rename(
            columns={
                "human_cluster": "best_positive_rankcorr_human_cluster",
                "spearman_r": "positive_rankcorr_spearman_r",
                "p_adj": "positive_rankcorr_fdr",
            }
        )
    )

    out = best_overlap[
        [
            "compartment",
            "mouse_region",
            "best_overlap_human_cluster",
            "overlap_best_fraction",
            "overlap_significant_fraction",
            "overlap_min_fdr",
            "overlap_max_genes",
        ]
    ].merge(
        best_module[
            [
                "compartment",
                "mouse_region",
                "best_module_human_cluster",
                "module_best_fraction",
                "module_mean_score_across_depths",
            ]
        ],
        on=["compartment", "mouse_region"],
        how="outer",
    ).merge(
        best_rank[
            [
                "compartment",
                "mouse_region",
                "best_positive_rankcorr_human_cluster",
                "positive_rankcorr_spearman_r",
                "positive_rankcorr_fdr",
                "n_ranked_genes",
            ]
        ],
        on=["compartment", "mouse_region"],
        how="outer",
    )
    out.to_csv(OUT_DIR / f"{compartment.lower()}_calibrated_best_match_summary.csv", index=False)
    return out


def allowed_human_clusters_for_compartment(compartment: str) -> list[str]:
    if compartment == "LHb":
        return [cluster for cluster in HUMAN_HB_POPULATIONS if cluster.startswith("LHb.")]
    if compartment == "MHb":
        return [cluster for cluster in HUMAN_HB_POPULATIONS if cluster.startswith("MHb.")]
    return HUMAN_HB_POPULATIONS


def write_same_compartment_summary(
    compartment: str,
    overlap_summary: pd.DataFrame,
    module_summary: pd.DataFrame,
    rank_corr: pd.DataFrame,
) -> pd.DataFrame:
    allowed = allowed_human_clusters_for_compartment(compartment)
    overlap_summary = overlap_summary[overlap_summary["human_cluster"].isin(allowed)].copy()
    module_summary = module_summary[module_summary["human_cluster"].isin(allowed)].copy()
    rank_corr = rank_corr[rank_corr["human_cluster"].isin(allowed)].copy()
    return write_compartment_summary(
        f"{compartment}_same_compartment",
        overlap_summary,
        module_summary,
        rank_corr,
    )


def run_compartment(compartment: str, human_markers: pd.DataFrame) -> pd.DataFrame:
    print(f"Running calibration for {compartment}")
    obs, counts, norm, var_names, categories = load_compartment(compartment)
    mapped = map_human_to_mouse(human_markers, var_names)
    mapped.to_csv(OUT_DIR / f"{compartment.lower()}_human_markers_mapped_to_mouse_symbols.csv", index=False)
    scores = compute_mouse_gene_scores(compartment, obs, counts, norm, var_names, categories)

    overlap_sensitivity, overlap_summary = run_overlap_sensitivity(compartment, mapped, scores, categories)
    module_sensitivity, module_summary = run_module_depth_sensitivity(compartment, mapped, obs, norm, var_names, categories)
    rank_corr = run_rank_correlation(compartment, mapped, scores, categories)

    plot_pair_heatmap(
        overlap_summary,
        compartment,
        "best_match_fraction",
        f"{compartment.lower()}_overlap_best_match_stability_heatmap",
        f"{compartment}: fraction of threshold grid where pair is best overlap match",
        cmap="mako",
        annot_col="significant_fraction",
    )
    plot_pair_heatmap(
        module_summary,
        compartment,
        "best_module_fraction",
        f"{compartment.lower()}_module_best_match_stability_heatmap",
        f"{compartment}: fraction of module depths where pair is top module score",
        cmap="crest",
        annot_col="mean_score_across_depths",
    )
    plot_pair_heatmap(
        rank_corr,
        compartment,
        "spearman_r",
        f"{compartment.lower()}_threshold_free_rank_correlation_heatmap",
        f"{compartment}: threshold-free human-vs-mouse marker score Spearman correlation",
        center=0,
        cmap="vlag",
        annot_col="p_adj",
    )

    same_clusters = allowed_human_clusters_for_compartment(compartment)
    same_overlap = overlap_summary[overlap_summary["human_cluster"].isin(same_clusters)].copy()
    same_module = module_summary[module_summary["human_cluster"].isin(same_clusters)].copy()
    same_rank = rank_corr[rank_corr["human_cluster"].isin(same_clusters)].copy()
    plot_pair_heatmap(
        same_overlap,
        compartment,
        "best_match_fraction",
        f"{compartment.lower()}_same_compartment_overlap_best_match_stability_heatmap",
        f"{compartment}: same-compartment overlap best-match stability",
        cmap="mako",
        annot_col="significant_fraction",
        human_clusters=same_clusters,
    )
    plot_pair_heatmap(
        same_module,
        compartment,
        "best_module_fraction",
        f"{compartment.lower()}_same_compartment_module_best_match_stability_heatmap",
        f"{compartment}: same-compartment module-score best-match stability",
        cmap="crest",
        annot_col="mean_score_across_depths",
        human_clusters=same_clusters,
    )
    plot_pair_heatmap(
        same_rank,
        compartment,
        "spearman_r",
        f"{compartment.lower()}_same_compartment_threshold_free_rank_correlation_heatmap",
        f"{compartment}: same-compartment threshold-free rank correlation",
        center=0,
        cmap="vlag",
        annot_col="p_adj",
        human_clusters=same_clusters,
    )

    unconstrained = write_compartment_summary(compartment, overlap_summary, module_summary, rank_corr)
    constrained = write_same_compartment_summary(compartment, overlap_summary, module_summary, rank_corr)
    constrained.to_csv(OUT_DIR / f"{compartment.lower()}_same_compartment_calibrated_best_match_summary.csv", index=False)
    return unconstrained


def append_summary_table(lines: list[str], summary: pd.DataFrame) -> None:
    lines.append("| compartment | mouse_region | overlap_best | overlap_best_fraction | overlap_significant_fraction | module_best | module_best_fraction | rankcorr_best | rankcorr_r | rankcorr_fdr |")
    lines.append("| --- | --- | --- | ---: | ---: | --- | ---: | --- | ---: | ---: |")
    for row in summary.itertuples():
        lines.append(
            f"| {row.compartment} | {row.mouse_region} | {row.best_overlap_human_cluster} | "
            f"{row.overlap_best_fraction:.2f} | {row.overlap_significant_fraction:.2f} | "
            f"{row.best_module_human_cluster} | {row.module_best_fraction:.2f} | "
            f"{row.best_positive_rankcorr_human_cluster} | {row.positive_rankcorr_spearman_r:.3f} | {row.positive_rankcorr_fdr:.3g} |"
        )


def write_overall_summary(all_summary: pd.DataFrame, same_compartment_summary: pd.DataFrame) -> None:
    all_summary.to_csv(OUT_DIR / "lhb_mhb_calibrated_best_match_summary.csv", index=False)
    same_compartment_summary.to_csv(OUT_DIR / "lhb_mhb_same_compartment_calibrated_best_match_summary.csv", index=False)
    plot_calibrated_best_match_summary(
        same_compartment_summary,
        "lhb_mhb_same_compartment_calibrated_best_match_summary_figure",
        "Calibrated Yalcinbas Human Hb Matches For Recovered Mouse LHb And MHb Regions",
    )
    plot_method_agreement_heatmap(
        same_compartment_summary,
        "lhb_mhb_same_compartment_method_agreement_heatmap",
        "Agreement Across Overlap, Module Score, And Rank-Correlation Matches",
    )
    lines = [
        "# Yalcinbas Threshold Calibration And MHb Extension",
        "",
        "This report adds threshold-sensitivity and threshold-free rank-correlation analyses to the fixed-threshold Yalcinbas comparison.",
        "",
        "Threshold grid:",
        f"- human overlap marker depths: {HUMAN_OVERLAP_DEPTHS}",
        f"- mouse overlap marker depths: {MOUSE_OVERLAP_DEPTHS}",
        f"- module marker depths: {MODULE_DEPTHS}",
        "",
        "Interpretation guide:",
        "- `overlap_best_fraction`: fraction of overlap-threshold grid settings where a human population was the best match for a mouse region.",
        "- `overlap_significant_fraction`: fraction of overlap-threshold grid settings where the pair was significant at grid-level FDR < 0.05.",
        "- `module_best_fraction`: fraction of module depths where a human population had the highest mean module score in a mouse region.",
        "- `rankcorr_spearman_r`: threshold-free correlation between human one-vs-rest logFC and mouse region-vs-rest logFC across all mapped genes.",
        "",
        "## Best calibrated matches across all human Hb populations",
        "",
    ]
    append_summary_table(lines, all_summary)
    lines.extend(
        [
            "",
            "## Same-compartment calibrated matches",
            "",
            "This table restricts mouse LHb regions to human `LHb.*` populations and mouse MHb regions to human `MHb.*` populations. It is the preferred table when asking specifically for LHb-to-LHb or MHb-to-MHb conservation.",
            "",
        ]
    )
    append_summary_table(lines, same_compartment_summary)
    lines.extend(
        [
            "",
            "Notes:",
            "- These calibration analyses do not replace the fixed marker-depth results; they test whether those results are stable to marker-depth choices.",
            "- The MHb extension maps cleaned medial habenula Leiden clusters to broad MHb regional marker programs.",
            "- Human `MHb.3` has very few nuclei in the Yalcinbas single-nucleus object, so MHb.3-related matches should be interpreted cautiously.",
        ]
    )
    # Summary Markdown output intentionally omitted from the publication package.


def main() -> pd.DataFrame:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    human_markers = load_human_markers()
    summaries = [run_compartment("LHb", human_markers), run_compartment("MHb", human_markers)]
    all_summary = pd.concat(summaries, ignore_index=True)
    same_compartment_summary = pd.concat(
        [
            pd.read_csv(OUT_DIR / "lhb_same_compartment_calibrated_best_match_summary.csv"),
            pd.read_csv(OUT_DIR / "mhb_same_compartment_calibrated_best_match_summary.csv"),
        ],
        ignore_index=True,
    )
    write_overall_summary(all_summary, same_compartment_summary)
    print(all_summary.to_string(index=False))
    print("\nSame-compartment summary:")
    print(same_compartment_summary.to_string(index=False))
    return all_summary


if __name__ == "__main__":
    main()
