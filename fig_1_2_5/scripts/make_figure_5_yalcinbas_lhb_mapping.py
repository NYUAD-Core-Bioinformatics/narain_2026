from __future__ import annotations

import os
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
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


CURRENT_H5AD = Path("inputs/clean_objects/whole_habenula.h5ad")
LHB_H5AD = Path("inputs/clean_objects/lateral_habenula_regions.h5ad")
HUMAN_MARKERS = Path("inputs/figure_5_cross_species/yalcinbas_single_nucleus_marker_genes.csv")
HUMAN_BULK_DE = Path("inputs/figure_5_cross_species/yalcinbas_bulk_schizophrenia_ranked_genes.csv")
CELLTYPE_PHENOTYPE_DE = Path("inputs/figure_2_mouse_differential_expression/celltype_phenotype_bootstrap_de.csv")
LHB_REGION_PHENOTYPE_DE = Path("inputs/figure_2_mouse_differential_expression/lhb_region_phenotype_de.csv")
OUT_DIR = Path("outputs/figure_5_cross_species/yalcinbas_marker_mapping")

HUMAN_HB_POPULATIONS = ["LHb.1", "LHb.2", "LHb.3", "LHb.4", "LHb.5", "LHb.6", "LHb.7", "MHb.1", "MHb.2", "MHb.3"]
LHB_REGION_ORDER = ["Hbx", "Lateral", "Marginal", "Oval/Medial"]
TOP_HUMAN_MARKERS_FOR_DOTPLOT = 6
TOP_HUMAN_MARKERS_FOR_MODULE = 25
TOP_HUMAN_MARKERS_FOR_OVERLAP = 50
TOP_MOUSE_MARKERS_FOR_OVERLAP = 200
HUMAN_DISEASE_TOP_N = 173


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


def aggregate_by_group(x: sparse.csr_matrix, groups: pd.Series, categories: list[str], binary: bool = False) -> np.ndarray:
    if binary:
        x = x.copy()
        x.data = np.ones_like(x.data)
    rows = []
    group_values = groups.astype(str).to_numpy()
    for category in categories:
        mask = group_values == category
        rows.append(np.asarray(x[mask, :].mean(axis=0)).ravel())
    return np.vstack(rows)


def gene_indices(var_names: list[str], genes: list[str]) -> dict[str, int]:
    lookup = {gene: i for i, gene in enumerate(var_names)}
    return {gene: lookup[gene] for gene in genes if gene in lookup}


def load_lhb() -> tuple[pd.DataFrame, sparse.csr_matrix, sparse.csr_matrix, list[str]]:
    a = ad.read_h5ad(LHB_H5AD)
    obs = a.obs.copy()
    obs["mouse_lhb_region"] = obs["cl"].astype(str)
    obs["mouse_lhb_region"] = obs["mouse_lhb_region"].replace({"Oval-Medial": "Oval/Medial", "HbX": "Hbx"})
    categories = [x for x in LHB_REGION_ORDER if x in set(obs["mouse_lhb_region"])]
    obs["mouse_lhb_region"] = pd.Categorical(obs["mouse_lhb_region"], categories=categories, ordered=True)
    counts = as_csr(a.layers["counts"])
    norm = as_csr(a.layers["l1p"])
    var_names = a.var_names.astype(str).tolist()
    return obs, counts, norm, var_names


def write_current_crosswalk(lhb_obs: pd.DataFrame) -> None:
    current = ad.read_h5ad(CURRENT_H5AD, backed="r")
    cur_obs = current.obs[["l04", "lcl", "gs", "pheno", "dataset"]].copy()
    current.file.close()

    joined = cur_obs.join(lhb_obs[["mouse_lhb_region"]], how="inner")
    joined.to_csv(OUT_DIR / "mouse_lhb_cells_with_region_labels.csv")

    counts = pd.crosstab(joined["l04"].astype(str), joined["mouse_lhb_region"].astype(str))
    counts = counts.reindex(sorted(counts.index, key=natural_key))
    counts.to_csv(OUT_DIR / "mouse_l04_by_lhb_region_counts.csv")
    row_pct = counts.div(counts.sum(axis=1), axis=0).fillna(0)
    row_pct.to_csv(OUT_DIR / "mouse_l04_by_lhb_region_fraction.csv")


def map_human_markers(var_names: list[str]) -> pd.DataFrame:
    markers = pd.read_csv(HUMAN_MARKERS)
    mouse_lookup = pd.DataFrame({"mouse_gene": var_names})
    mouse_lookup["gene_upper"] = mouse_lookup["mouse_gene"].str.upper()
    mouse_lookup = mouse_lookup.drop_duplicates("gene_upper")

    symbol = pd.Series(index=markers.index, dtype=object)
    for candidate in ["Symbol", "ID.x", "ID", "gene_id"]:
        if candidate in markers.columns:
            values = markers[candidate].astype(object)
            values = values.where(values.notna() & values.astype(str).ne("nan") & values.astype(str).ne(""))
            symbol = symbol.fillna(values)
    mapped = markers.copy()
    mapped["human_gene"] = symbol.astype(str)
    mapped = mapped[mapped["human_gene"].notna() & mapped["human_gene"].ne("nan")].copy()
    mapped["gene_upper"] = mapped["human_gene"].str.upper()
    mapped = mapped.merge(mouse_lookup, on="gene_upper", how="inner")
    mapped = mapped.sort_values(["human_cluster", "adj.P.Val", "P.Value", "logFC"], ascending=[True, True, True, False])
    mapped.to_csv(OUT_DIR / "yalcinbas_lhb_markers_mapped_to_mouse.csv", index=False)
    return mapped


def top_human_marker_sets(mapped: pd.DataFrame, n: int) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for cluster in HUMAN_HB_POPULATIONS:
        sub = mapped[mapped["human_cluster"].eq(cluster) & mapped["logFC"].gt(0)].drop_duplicates("gene_upper")
        out[cluster] = sub.head(n)["mouse_gene"].tolist()
    return out


def compute_mouse_region_markers(
    lhb_obs: pd.DataFrame,
    counts: sparse.csr_matrix,
    norm: sparse.csr_matrix,
    var_names: list[str],
) -> dict[str, set[str]]:
    categories = list(lhb_obs["mouse_lhb_region"].cat.categories)
    cluster_sizes = lhb_obs["mouse_lhb_region"].value_counts().reindex(categories).to_numpy()
    mean_by_cluster = aggregate_by_group(norm, lhb_obs["mouse_lhb_region"], categories)
    pct_by_cluster = aggregate_by_group(counts, lhb_obs["mouse_lhb_region"], categories, binary=True)
    total_sum = np.asarray(norm.sum(axis=0)).ravel()
    total_n = norm.shape[0]

    records = []
    marker_sets: dict[str, set[str]] = {}
    for i, category in enumerate(categories):
        rest_n = total_n - cluster_sizes[i]
        rest_mean = (total_sum - mean_by_cluster[i] * cluster_sizes[i]) / max(rest_n, 1)
        df = pd.DataFrame(
            {
                "mouse_lhb_region": category,
                "mouse_gene": var_names,
                "avg_log_expr": mean_by_cluster[i],
                "pct_expr": pct_by_cluster[i],
                "logFC_vs_other_lhb_regions": mean_by_cluster[i] - rest_mean,
            }
        )
        df = df[
            df["logFC_vs_other_lhb_regions"].gt(0)
            & df["pct_expr"].ge(0.05)
            & ~df["mouse_gene"].str.startswith(("Rp", "mt-"))
        ].copy()
        df = df.sort_values("logFC_vs_other_lhb_regions", ascending=False)
        marker_sets[category] = set(df.head(TOP_MOUSE_MARKERS_FOR_OVERLAP)["mouse_gene"].str.upper())
        records.append(df.head(1000))

    pd.concat(records, ignore_index=True).to_csv(OUT_DIR / "mouse_lhb_region_marker_genes.csv", index=False)
    return marker_sets


def plot_dotplot(
    lhb_obs: pd.DataFrame,
    counts: sparse.csr_matrix,
    norm: sparse.csr_matrix,
    var_names: list[str],
    human_sets: dict[str, list[str]],
) -> None:
    categories = list(lhb_obs["mouse_lhb_region"].cat.categories)
    selected_rows = []
    selected_genes = []
    for human_cluster in HUMAN_HB_POPULATIONS:
        for gene in human_sets.get(human_cluster, [])[:TOP_HUMAN_MARKERS_FOR_DOTPLOT]:
            if gene not in selected_genes:
                selected_genes.append(gene)
                selected_rows.append((human_cluster, gene))

    idx = gene_indices(var_names, selected_genes)
    selected_genes = [gene for gene in selected_genes if gene in idx]
    selected_rows = [(cluster, gene) for cluster, gene in selected_rows if gene in idx]
    positions = [idx[gene] for gene in selected_genes]

    mean_expr = aggregate_by_group(norm[:, positions], lhb_obs["mouse_lhb_region"], categories)
    pct_expr = aggregate_by_group(counts[:, positions], lhb_obs["mouse_lhb_region"], categories, binary=True) * 100
    z_expr = (mean_expr - mean_expr.mean(axis=0, keepdims=True)) / (mean_expr.std(axis=0, keepdims=True) + 1e-6)
    z_expr = np.clip(z_expr, -2.5, 2.5)

    gene_to_label = {gene: f"{cluster}: {gene}" for cluster, gene in selected_rows}
    records = []
    for gene_i, gene in enumerate(selected_genes):
        for cluster_i, region in enumerate(categories):
            records.append(
                {
                    "human_marker_gene": gene_to_label[gene],
                    "mouse_lhb_region": region,
                    "mean_expr_z": z_expr[cluster_i, gene_i],
                    "pct_expr": pct_expr[cluster_i, gene_i],
                }
            )
    dot = pd.DataFrame(records)
    dot.to_csv(OUT_DIR / "yalcinbas_lhb_marker_dotplot_values.csv", index=False)

    y_labels = [gene_to_label[gene] for gene in selected_genes]
    y_pos = {label: i for i, label in enumerate(y_labels)}
    x_pos = {region: i for i, region in enumerate(categories)}

    fig, ax = plt.subplots(figsize=(8, max(8, len(y_labels) * 0.24)))
    scatter = ax.scatter(
        dot["mouse_lhb_region"].map(x_pos),
        dot["human_marker_gene"].map(y_pos),
        s=np.maximum(dot["pct_expr"], 1) * 3.0,
        c=dot["mean_expr_z"],
        cmap="vlag",
        vmin=-2.5,
        vmax=2.5,
        linewidth=0.2,
        edgecolor="0.35",
    )
    ax.set_xticks(range(len(categories)), categories, rotation=30, ha="right")
    ax.set_yticks(range(len(y_labels)), y_labels)
    ax.invert_yaxis()
    ax.set_xlabel("Recovered mouse LHb region")
    ax.set_ylabel("Mapped human Hb marker genes")
    ax.set_title("Analysis A: human Hb markers across recovered mouse LHb regions")
    cbar = fig.colorbar(scatter, ax=ax, pad=0.01)
    cbar.set_label("Mean expression z-score across mouse regions")
    for size in [10, 30, 60]:
        ax.scatter([], [], s=size * 3.0, c="lightgray", edgecolor="0.35", label=f"{size}%")
    ax.legend(title="Pct expressed", bbox_to_anchor=(1.18, 1.0), loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "yalcinbas_lhb_marker_dotplot.png", dpi=220)
    fig.savefig(OUT_DIR / "yalcinbas_lhb_marker_dotplot.pdf")
    plt.close(fig)


def plot_module_scores(
    lhb_obs: pd.DataFrame,
    norm: sparse.csr_matrix,
    var_names: list[str],
    human_sets: dict[str, list[str]],
) -> pd.DataFrame:
    categories = list(lhb_obs["mouse_lhb_region"].cat.categories)
    union_genes = sorted({gene for genes in human_sets.values() for gene in genes})
    idx = gene_indices(var_names, union_genes)
    union_genes = [gene for gene in union_genes if gene in idx]
    positions = [idx[gene] for gene in union_genes]
    gene_position = {gene: i for i, gene in enumerate(union_genes)}

    dense = norm[:, positions].toarray().astype(np.float32)
    dense = (dense - dense.mean(axis=0, keepdims=True)) / (dense.std(axis=0, keepdims=True) + 1e-6)

    score_records = []
    score_matrix = pd.DataFrame(index=HUMAN_HB_POPULATIONS, columns=categories, dtype=float)
    for human_cluster in HUMAN_HB_POPULATIONS:
        genes = [gene for gene in human_sets.get(human_cluster, []) if gene in gene_position]
        if not genes:
            continue
        score = dense[:, [gene_position[gene] for gene in genes]].mean(axis=1)
        for region in categories:
            mask = lhb_obs["mouse_lhb_region"].astype(str).eq(region).to_numpy()
            value = float(np.mean(score[mask]))
            score_matrix.loc[human_cluster, region] = value
            score_records.append(
                {
                    "human_cluster": human_cluster,
                    "mouse_lhb_region": region,
                    "mean_module_score": value,
                    "n_marker_genes": len(genes),
                }
            )

    score_table = pd.DataFrame(score_records)
    score_table.to_csv(OUT_DIR / "yalcinbas_lhb_module_scores.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 6.2))
    sns.heatmap(
        score_matrix.astype(float),
        cmap="vlag",
        center=0,
        linewidths=0.3,
        linecolor="white",
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Mean marker module score"},
        ax=ax,
    )
    ax.set_xlabel("Recovered mouse LHb region")
    ax.set_ylabel("Human Hb marker module")
    ax.set_title("Analysis A: human marker module scores by recovered LHb region")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "yalcinbas_lhb_module_score_heatmap.png", dpi=220)
    fig.savefig(OUT_DIR / "yalcinbas_lhb_module_score_heatmap.pdf")
    plt.close(fig)
    return score_table


def plot_marker_overlap(mapped: pd.DataFrame, mouse_marker_sets: dict[str, set[str]]) -> pd.DataFrame:
    human_sets: dict[str, set[str]] = {}
    for human_cluster in HUMAN_HB_POPULATIONS:
        sub = mapped[mapped["human_cluster"].eq(human_cluster) & mapped["logFC"].gt(0)].drop_duplicates("gene_upper")
        human_sets[human_cluster] = set(sub.head(TOP_HUMAN_MARKERS_FOR_OVERLAP)["mouse_gene"].str.upper())

    universe = set(mapped["mouse_gene"].str.upper())
    rows = []
    for region, mouse_set in mouse_marker_sets.items():
        mouse_set = mouse_set & universe
        for human_cluster, human_set in human_sets.items():
            human_set = human_set & universe
            overlap = sorted(mouse_set & human_set)
            a = len(overlap)
            b = len(mouse_set - human_set)
            c = len(human_set - mouse_set)
            d = max(len(universe) - a - b - c, 0)
            odds, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append(
                {
                    "mouse_lhb_region": region,
                    "human_cluster": human_cluster,
                    "n_mouse_markers": len(mouse_set),
                    "n_human_markers": len(human_set),
                    "n_overlap": a,
                    "odds_ratio": odds,
                    "p_value": p_value,
                    "overlap_genes": ";".join(overlap),
                }
            )
    overlap = pd.DataFrame(rows)
    overlap["p_adj"] = multipletests(overlap["p_value"], method="fdr_bh")[1]
    overlap["neg_log10_p_adj"] = -np.log10(np.maximum(overlap["p_adj"], 1e-300))
    overlap.to_csv(OUT_DIR / "yalcinbas_lhb_region_marker_overlap.csv", index=False)

    heat = overlap.pivot(index="human_cluster", columns="mouse_lhb_region", values="neg_log10_p_adj").reindex(HUMAN_HB_POPULATIONS)
    annot = overlap.pivot(index="human_cluster", columns="mouse_lhb_region", values="n_overlap").reindex(HUMAN_HB_POPULATIONS)
    fig, ax = plt.subplots(figsize=(8, 6.2))
    sns.heatmap(
        heat,
        cmap="mako",
        linewidths=0.3,
        linecolor="white",
        annot=annot,
        fmt=".0f",
        cbar_kws={"label": "-log10 FDR, Fisher overlap"},
        ax=ax,
    )
    ax.set_xlabel("Recovered mouse LHb region marker set")
    ax.set_ylabel("Human Hb marker set")
    ax.set_title("Analysis A: marker overlap using recovered LHb regions")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "yalcinbas_lhb_region_marker_overlap_heatmap.png", dpi=220)
    fig.savefig(OUT_DIR / "yalcinbas_lhb_region_marker_overlap_heatmap.pdf")
    plt.close(fig)

    best = overlap.sort_values(["mouse_lhb_region", "p_adj", "n_overlap"], ascending=[True, True, False])
    best = best.groupby("mouse_lhb_region", as_index=False).head(1)
    best.to_csv(OUT_DIR / "yalcinbas_top_human_cluster_per_lhb_region.csv", index=False)
    return best


def human_scz_gene_sets() -> tuple[set[str], set[str], set[str], pd.DataFrame]:
    human = pd.read_csv(HUMAN_BULK_DE)
    human["human_gene"] = human["Symbol"].fillna(human["MGI_Symbol"]).astype(str)
    human["gene_upper"] = human["human_gene"].str.upper()
    human = human[human["human_gene"].notna() & human["human_gene"].ne("nan")].drop_duplicates("gene_upper")
    ranked = human.sort_values("P.Value").head(HUMAN_DISEASE_TOP_N).copy()
    ranked.to_csv(OUT_DIR / "yalcinbas_bulk_schizophrenia_top173_genes.csv", index=False)
    all_set = set(ranked["gene_upper"])
    up_set = set(ranked.loc[ranked["logFC"].gt(0), "gene_upper"])
    down_set = set(ranked.loc[ranked["logFC"].lt(0), "gene_upper"])
    return all_set, up_set, down_set, human


def overlap_mouse_de(
    mouse_de: pd.DataFrame,
    group_cols: list[str],
    gene_col: str,
    logfc_col: str,
    out_prefix: str,
) -> pd.DataFrame:
    human_all, human_up, human_down, human = human_scz_gene_sets()
    mouse_de = mouse_de.copy()
    mouse_de["gene_upper"] = mouse_de[gene_col].astype(str).str.upper()
    mouse_de = mouse_de[mouse_de["gene_upper"].ne("NAN")].copy()
    universe = set(human["gene_upper"]) | set(mouse_de["gene_upper"])

    rows = []
    for key, group in mouse_de.groupby(group_cols, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        mouse_all = set(group["gene_upper"])
        mouse_up = set(group.loc[group[logfc_col].gt(0), "gene_upper"])
        mouse_down = set(group.loc[group[logfc_col].lt(0), "gene_upper"])
        overlap = sorted(mouse_all & human_all)
        concordant = sorted((mouse_up & human_up) | (mouse_down & human_down))
        discordant = sorted((mouse_up & human_down) | (mouse_down & human_up))
        a = len(overlap)
        b = len(mouse_all - human_all)
        c = len(human_all - mouse_all)
        d = max(len(universe) - a - b - c, 0)
        odds, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        row = {col: value for col, value in zip(group_cols, key)}
        row.update(
            {
                "n_mouse_de_genes": len(mouse_all),
                "n_human_scz_genes": len(human_all),
                "n_overlap": a,
                "n_concordant_direction": len(concordant),
                "n_discordant_direction": len(discordant),
                "odds_ratio": odds,
                "p_value": p_value,
                "overlap_genes": ";".join(overlap),
                "concordant_genes": ";".join(concordant),
                "discordant_genes": ";".join(discordant),
            }
        )
        rows.append(row)
    out = pd.DataFrame(rows)
    out["p_adj"] = multipletests(out["p_value"], method="fdr_bh")[1]
    out["neg_log10_p_adj"] = -np.log10(np.maximum(out["p_adj"], 1e-300))
    out = out.sort_values(["p_adj", "p_value", "n_overlap"], ascending=[True, True, False])
    out.to_csv(OUT_DIR / f"{out_prefix}.csv", index=False)

    label_col = group_cols[0]
    heat = out.pivot(index=label_col, columns=group_cols[-1], values="neg_log10_p_adj")
    annot = out.pivot(index=label_col, columns=group_cols[-1], values="n_overlap")
    fig, ax = plt.subplots(figsize=(max(7, heat.shape[1] * 1.2), max(4, heat.shape[0] * 0.45)))
    sns.heatmap(
        heat,
        cmap="rocket_r",
        linewidths=0.3,
        linecolor="white",
        annot=annot,
        fmt=".0f",
        cbar_kws={"label": "-log10 FDR, Fisher overlap"},
        ax=ax,
    )
    ax.set_title(out_prefix.replace("_", " "))
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"{out_prefix}_heatmap.png", dpi=220)
    fig.savefig(OUT_DIR / f"{out_prefix}_heatmap.pdf")
    plt.close(fig)
    return out


def run_analysis_b() -> tuple[pd.DataFrame, pd.DataFrame]:
    boot = pd.read_csv(CELLTYPE_PHENOTYPE_DE)
    boot = boot.rename(columns={"ct": "celltype", "names": "mouse_gene", "avg_logFC": "logFC"})
    boot_overlap = overlap_mouse_de(
        boot,
        group_cols=["celltype", "pheno"],
        gene_col="mouse_gene",
        logfc_col="logFC",
        out_prefix="yalcinbas_bulk_schizophrenia_celltype_overlap",
    )

    zone = pd.read_csv(LHB_REGION_PHENOTYPE_DE)
    zone = zone.rename(columns={"names": "mouse_gene", "logfoldchanges": "logFC"})
    zone["region"] = zone["region"].replace({"Oval-Medial": "Oval/Medial", "HbX": "Hbx"})
    zone_overlap = overlap_mouse_de(
        zone,
        group_cols=["region", "pheno"],
        gene_col="mouse_gene",
        logfc_col="logFC",
        out_prefix="yalcinbas_bulk_schizophrenia_lhb_region_overlap",
    )
    return boot_overlap, zone_overlap


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    lhb_obs, counts, norm, var_names = load_lhb()
    write_current_crosswalk(lhb_obs)

    mapped = map_human_markers(var_names)
    human_dot_sets = top_human_marker_sets(mapped, TOP_HUMAN_MARKERS_FOR_DOTPLOT)
    human_module_sets = top_human_marker_sets(mapped, TOP_HUMAN_MARKERS_FOR_MODULE)

    mouse_marker_sets = compute_mouse_region_markers(lhb_obs, counts, norm, var_names)
    plot_dotplot(lhb_obs, counts, norm, var_names, human_dot_sets)
    module_scores = plot_module_scores(lhb_obs, norm, var_names, human_module_sets)
    best_matches = plot_marker_overlap(mapped, mouse_marker_sets)
    boot_overlap, zone_overlap = run_analysis_b()

    summary = {
        "lhb_cells": int(lhb_obs.shape[0]),
        "lhb_regions": ";".join([f"{k}={v}" for k, v in lhb_obs["mouse_lhb_region"].value_counts().reindex(lhb_obs["mouse_lhb_region"].cat.categories).items()]),
        "best_analysis_a_matches": ";".join(
            f"{r.mouse_lhb_region}->{r.human_cluster} FDR={r.p_adj:.3g} overlap={int(r.n_overlap)}"
            for r in best_matches.itertuples()
        ),
        "analysis_b_bootstrap_min_fdr": float(boot_overlap["p_adj"].min()),
        "analysis_b_bootstrap_fdr_lt_0_1": int(boot_overlap["p_adj"].lt(0.1).sum()),
        "analysis_b_lhb_region_min_fdr": float(zone_overlap["p_adj"].min()),
        "analysis_b_lhb_region_fdr_lt_0_1": int(zone_overlap["p_adj"].lt(0.1).sum()),
        "strongest_human_module_by_region": ";".join(
            module_scores.sort_values(["mouse_lhb_region", "mean_module_score"], ascending=[True, False])
            .groupby("mouse_lhb_region", as_index=False)
            .head(1)
            .apply(lambda r: f"{r['mouse_lhb_region']}->{r['human_cluster']} score={r['mean_module_score']:.3f}", axis=1)
        ),
    }
    pd.Series(summary).to_csv(OUT_DIR / "yalcinbas_lhb_mapping_summary.csv", header=["value"])
    print(pd.Series(summary).to_string())


if __name__ == "__main__":
    main()
