from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/private/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib_cache")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/xdg_cache")
for key in ["NUMBA_CACHE_DIR", "MPLCONFIGDIR", "XDG_CACHE_HOME"]:
    Path(os.environ[key]).mkdir(parents=True, exist_ok=True)

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse


H5AD_PATH = Path("inputs/clean_objects/whole_habenula.h5ad")
OUT_DIR = Path("outputs/figure_1_mouse_habenula_atlas")

BROAD_ORDER = [
    "LHb",
    "MHb",
    "Oligodendrocytes",
    "Diff. Oligodendrocytes",
    "Polydendrocytes",
    "Fibroblasts",
    "Astrocytes",
    "Endothelial",
    "Pericytes",
    "Macrophages",
    "Microglia",
]

BROAD_MARKERS = [
    "Snap25",
    "Celf4",
    "Slc17a6",
    "Tac2",
    "Gap43",
    "Mog",
    "Olig1",
    "Gpr17",
    "Pdgfra",
    "Col3a1",
    "Slc6a11",
    "Slc6a1",
    "Slc4a4",
    "Cldn5",
    "Abcc9",
    "Pdgfrb",
    "Cx3cr1",
    "C1qa",
    "C1qb",
    "C1qc",
    "Mrc1",
    "Apoe",
]

MHB_ASSIGNMENT_MARKERS = {
    "Ventral": ["Lmo3"],
    "Ventral 2/3": ["Slc18a3"],
    "Lateral": ["Sema3d"],
    "Dorsal": ["Wif1"],
    "Superior": ["Avil"],
}

MHB_HEATMAP_MARKERS = ["Tac2", "Lmo3", "Slc18a3", "Sema3d", "Wif1", "Avil"]
MHB_HEATMAP_ORDER = ["Dorsal", "Ventral 2/3", "Unknown", "Lateral", "Superior", "Ventral"]

LHB_ASSIGNMENT_MARKERS = {
    "HbX": ["Cartpt", "Chrnb3"],
    "Lateral": ["Plch1", "Peg10", "Pbx3"],
    "Marginal": ["Vgf"],
    "Oval/Medial": ["Chrm3"],
}

LHB_HEATMAP_MARKERS = [
    "Gpd2",
    "Paqr8",
    "Parm1",
    "Gpr151",
    "Sst",
    "Bcl11b",
    "Th",
    "Chrm3",
    "Vgf",
    "Plch1",
    "Peg10",
    "Pbx3",
    "Cartpt",
    "Chrnb3",
]
LHB_HEATMAP_ORDER = ["HbX", "Lateral", "Marginal", "Oval/Medial"]

PALETTE = {
    "LHb": "#2a6fbb",
    "MHb": "#2ca25f",
    "Oligodendrocytes": "#8c510a",
    "Diff. Oligodendrocytes": "#80cdc1",
    "Polydendrocytes": "#b2df8a",
    "Fibroblasts": "#f1b6da",
    "Astrocytes": "#f6e600",
    "Endothelial": "#d01c8b",
    "Pericytes": "#253494",
    "Macrophages": "#b2182b",
    "Microglia": "#fbb4ae",
}

MHB_PALETTE = {
    "Ventral": "#f46d43",
    "Ventral 2/3": "#2c7fb8",
    "Unknown": "#41ab5d",
    "Lateral": "#d95fbc",
    "Dorsal": "#8c564b",
    "Superior": "#e31a1c",
}

LHB_PALETTE = {
    "HbX": "#2c7fb8",
    "Lateral": "#f46d43",
    "Marginal": "#1f9e89",
    "Oval/Medial": "#e6c800",
}


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


def load_h5ad_metadata(h5ad_path: Path = H5AD_PATH) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    adata = sc.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs.copy()
    umap = pd.DataFrame(adata.obsm["X_umap"], index=obs.index, columns=["UMAP1", "UMAP2"])
    var_names = adata.var_names.astype(str).tolist()
    adata.file.close()
    return obs, umap, var_names


def read_counts_layer(h5ad_path: Path = H5AD_PATH) -> sparse.csr_matrix:
    with h5py.File(h5ad_path, "r") as h5:
        group = h5["layers"]["counts"]
        shape = tuple(int(x) for x in group.attrs["shape"])
        counts = sparse.csr_matrix(
            (group["data"][:], group["indices"][:], group["indptr"][:]),
            shape=shape,
        )
    counts.sum_duplicates()
    return counts


def gene_positions(var_names: list[str], genes: list[str]) -> tuple[list[str], list[int]]:
    lookup = {gene.upper(): (i, gene) for i, gene in enumerate(var_names)}
    mapped = []
    positions = []
    for gene in genes:
        item = lookup.get(gene.upper())
        if item is not None:
            pos, actual = item
            positions.append(pos)
            mapped.append(actual)
    return mapped, positions


def normalized_marker_matrix(counts: sparse.csr_matrix, positions: list[int], target_sum: float = 1e4) -> sparse.csr_matrix:
    x = counts[:, positions].astype(np.float32, copy=True).tocsr()
    lib = np.asarray(counts.sum(axis=1)).ravel()
    scale = np.zeros_like(lib, dtype=np.float32)
    keep = lib > 0
    scale[keep] = target_sum / lib[keep]
    x = x.multiply(scale[:, None]).tocsr()
    x.data = np.log1p(x.data)
    return x


def zscore_dense(x: np.ndarray) -> np.ndarray:
    return (x - x.mean(axis=0, keepdims=True)) / (x.std(axis=0, keepdims=True) + 1e-6)


def plot_panel_label(ax, label: str) -> None:
    ax.text(-0.08, 1.04, label, transform=ax.transAxes, fontsize=24, fontweight="bold", va="bottom")


def add_cluster_labels(
    ax,
    df: pd.DataFrame,
    label_col: str,
    min_cells: int = 20,
    fontsize: int = 12,
    offsets: dict[str, tuple[float, float]] | None = None,
) -> None:
    offsets = offsets or {}
    for label, group in df.groupby(label_col, observed=True):
        if group.shape[0] < min_cells:
            continue
        xy = group[["UMAP1", "UMAP2"]].median().to_numpy()
        dx, dy = offsets.get(str(label), (0.0, 0.0))
        xy[0] += dx
        xy[1] += dy
        ax.text(
            xy[0],
            xy[1],
            str(label),
            fontsize=fontsize,
            fontweight="bold",
            ha="center",
            va="center",
            color="black",
            path_effects=[pe.withStroke(linewidth=3, foreground="white")],
        )


def plot_umap(
    ax,
    df: pd.DataFrame,
    label_col: str,
    order: list[str],
    palette: dict[str, str],
    title: str = "",
    label_fontsize: int = 12,
    label_offsets: dict[str, tuple[float, float]] | None = None,
) -> None:
    for label in order:
        mask = df[label_col].astype(str).eq(label)
        if not mask.any():
            continue
        ax.scatter(
            df.loc[mask, "UMAP1"],
            df.loc[mask, "UMAP2"],
            s=0.9,
            c=palette.get(label, "lightgray"),
            alpha=0.85,
            linewidths=0,
            rasterized=True,
        )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title, fontsize=12)
    for spine in ax.spines.values():
        spine.set_visible(False)
    add_cluster_labels(ax, df, label_col, fontsize=label_fontsize, offsets=label_offsets)


def broad_dotplot_values(obs: pd.DataFrame, counts: sparse.csr_matrix, var_names: list[str]) -> pd.DataFrame:
    genes, positions = gene_positions(var_names, BROAD_MARKERS)
    x_norm = normalized_marker_matrix(counts, positions)
    x_bin = counts[:, positions].copy()
    x_bin.data = np.ones_like(x_bin.data)

    records = []
    for group in BROAD_ORDER:
        mask = obs["gs"].astype(str).eq(group).to_numpy()
        if not mask.any():
            continue
        means = np.asarray(x_norm[mask, :].mean(axis=0)).ravel()
        pcts = np.asarray(x_bin[mask, :].mean(axis=0)).ravel() * 100
        for gene, mean, pct in zip(genes, means, pcts):
            records.append({"group": group, "gene": gene, "mean_expr": mean, "pct_expr": pct})

    df = pd.DataFrame(records)
    df["mean_expr_scaled"] = df.groupby("gene")["mean_expr"].transform(
        lambda s: (s - s.min()) / (s.max() - s.min() + 1e-6)
    )
    return df


def plot_dotplot(ax, dot: pd.DataFrame) -> None:
    genes = [g for g in BROAD_MARKERS if g in set(dot["gene"])]
    groups = [g for g in BROAD_ORDER if g in set(dot["group"])]
    x_pos = {gene: i for i, gene in enumerate(genes)}
    y_pos = {group: i for i, group in enumerate(groups)}

    sca = ax.scatter(
        dot["gene"].map(x_pos),
        dot["group"].map(y_pos),
        s=np.maximum(dot["pct_expr"], 1) * 2.2,
        c=dot["mean_expr_scaled"],
        cmap="Reds",
        vmin=0,
        vmax=1,
        edgecolor="0.35",
        linewidth=0.25,
    )
    ax.set_xticks(range(len(genes)), genes, rotation=90, fontsize=7)
    ax.set_yticks(range(len(groups)), groups, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("")
    ax.set_ylabel("")
    for size in [10, 30, 50, 70, 90]:
        ax.scatter([], [], s=size * 2.2, c="lightgray", edgecolor="0.35", label=f"{size}")
    leg = ax.legend(title="Fraction of cells\nin group (%)", frameon=False, bbox_to_anchor=(1.05, 0.72), loc="upper left", fontsize=7)
    leg.get_title().set_fontsize(7)
    cbar = plt.colorbar(sca, ax=ax, fraction=0.045, pad=0.22)
    cbar.set_label("Mean expression\nin group", fontsize=7)
    cbar.ax.tick_params(labelsize=7)


def module_score_table(
    obs_subset: pd.DataFrame,
    counts_subset: sparse.csr_matrix,
    var_names: list[str],
    marker_sets: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_cols = {}
    mapped_records = []
    for label, genes in marker_sets.items():
        mapped, positions = gene_positions(var_names, genes)
        mapped_records.append({"label": label, "requested_genes": ";".join(genes), "mapped_genes": ";".join(mapped)})
        if not positions:
            score_cols[label] = np.zeros(obs_subset.shape[0], dtype=np.float32)
            continue
        x = normalized_marker_matrix(counts_subset, positions).toarray().astype(np.float32)
        score_cols[label] = zscore_dense(x).mean(axis=1)

    scores = pd.DataFrame(score_cols, index=obs_subset.index)
    mapped = pd.DataFrame(mapped_records)
    return scores, mapped


def assign_subregions(
    obs_subset: pd.DataFrame,
    counts_subset: sparse.csr_matrix,
    var_names: list[str],
    marker_sets: dict[str, list[str]],
    output_col: str,
    min_cluster_cells: int = 100,
    unknown_threshold: float | None = None,
    unknown_margin: float = 0.02,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    obs_subset = obs_subset.copy()
    cluster_counts = obs_subset["l04"].astype(str).value_counts()
    keep_clusters = sorted(cluster_counts[cluster_counts.ge(min_cluster_cells)].index.tolist(), key=natural_key)
    keep_mask = obs_subset["l04"].astype(str).isin(keep_clusters).to_numpy()
    obs_subset = obs_subset.loc[keep_mask].copy()
    counts_subset = counts_subset[keep_mask, :].tocsr()

    scores, mapped = module_score_table(obs_subset, counts_subset, var_names, marker_sets)
    cluster_scores = (
        scores.assign(l04=obs_subset["l04"].astype(str).to_numpy())
        .groupby("l04", observed=True)
        .mean()
        .reindex(keep_clusters)
    )

    cluster_assignments = []
    for cluster, row in cluster_scores.iterrows():
        ranked = row.sort_values(ascending=False)
        best_label = str(ranked.index[0])
        best_score = float(ranked.iloc[0])
        second_score = float(ranked.iloc[1]) if ranked.shape[0] > 1 else -np.inf
        assigned = best_label
        if unknown_threshold is not None:
            if best_score < unknown_threshold or (best_score - second_score) < unknown_margin:
                assigned = "Unknown"
        cluster_assignments.append(
            {
                "l04": cluster,
                output_col: assigned,
                "best_marker_label": best_label,
                "best_marker_score": best_score,
                "second_marker_score": second_score,
                "n_cells": int(cluster_counts.loc[cluster]),
            }
        )
    assignment = pd.DataFrame(cluster_assignments)
    obs_subset[output_col] = obs_subset["l04"].astype(str).map(dict(zip(assignment["l04"], assignment[output_col])))
    return obs_subset, assignment, mapped


def expression_heatmap_values(
    obs_subset: pd.DataFrame,
    counts_subset: sparse.csr_matrix,
    var_names: list[str],
    group_col: str,
    genes: list[str],
    group_order: list[str],
) -> pd.DataFrame:
    mapped, positions = gene_positions(var_names, genes)
    x = normalized_marker_matrix(counts_subset, positions).toarray().astype(np.float32)
    x_z = zscore_dense(x)
    rows = []
    for group in group_order:
        mask = obs_subset[group_col].astype(str).eq(group).to_numpy()
        if not mask.any():
            rows.append(pd.Series(np.nan, index=mapped, name=group))
        else:
            rows.append(pd.Series(x_z[mask, :].mean(axis=0), index=mapped, name=group))
    return pd.DataFrame(rows)


def plot_heatmap(ax, values: pd.DataFrame, title: str, cbar_label: str = "mean z-score") -> None:
    sns.heatmap(
        values,
        cmap="bwr",
        center=0,
        vmin=-1,
        vmax=1,
        linewidths=0.35,
        linecolor="white",
        cbar_kws={"label": cbar_label, "shrink": 0.55},
        ax=ax,
    )
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=90, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)


def recreate_figure1(h5ad_path: Path = H5AD_PATH, out_dir: Path = OUT_DIR) -> dict[str, Path]:
    out_dir.mkdir(exist_ok=True)
    obs, umap, var_names = load_h5ad_metadata(h5ad_path)
    counts = read_counts_layer(h5ad_path)
    plot_df = pd.concat([obs[["gs", "l04"]], umap], axis=1)

    broad_dot = broad_dotplot_values(obs, counts, var_names)
    broad_dot.to_csv(out_dir / "panel_b_broad_dotplot_values.csv", index=False)

    mhb_mask = obs["gs"].astype(str).eq("MHb").to_numpy()
    lhb_mask = obs["gs"].astype(str).eq("LHb").to_numpy()

    mhb_obs = obs.loc[mhb_mask].copy()
    mhb_umap = umap.loc[mhb_obs.index].copy()
    mhb_counts = counts[mhb_mask, :].tocsr()
    mhb_obs, mhb_assignment, mhb_mapped = assign_subregions(
        mhb_obs,
        mhb_counts,
        var_names,
        MHB_ASSIGNMENT_MARKERS,
        "figure1_mhb_subregion",
        unknown_threshold=0.03,
        unknown_margin=0.015,
    )
    mhb_umap = mhb_umap.loc[mhb_obs.index].copy()
    mhb_plot = pd.concat([mhb_obs[["figure1_mhb_subregion", "l04"]], mhb_umap], axis=1)
    mhb_counts_used = counts[obs.index.isin(mhb_obs.index), :].tocsr()
    mhb_heat = expression_heatmap_values(
        mhb_obs,
        mhb_counts_used,
        var_names,
        "figure1_mhb_subregion",
        MHB_HEATMAP_MARKERS,
        [g for g in MHB_HEATMAP_ORDER if g in set(mhb_obs["figure1_mhb_subregion"])],
    )

    lhb_obs = obs.loc[lhb_mask].copy()
    lhb_umap = umap.loc[lhb_obs.index].copy()
    lhb_counts = counts[lhb_mask, :].tocsr()
    lhb_obs, lhb_assignment, lhb_mapped = assign_subregions(
        lhb_obs,
        lhb_counts,
        var_names,
        LHB_ASSIGNMENT_MARKERS,
        "figure1_lhb_subregion",
        unknown_threshold=None,
    )
    lhb_umap = lhb_umap.loc[lhb_obs.index].copy()
    lhb_plot = pd.concat([lhb_obs[["figure1_lhb_subregion", "l04"]], lhb_umap], axis=1)
    lhb_counts_used = counts[obs.index.isin(lhb_obs.index), :].tocsr()
    lhb_heat = expression_heatmap_values(
        lhb_obs,
        lhb_counts_used,
        var_names,
        "figure1_lhb_subregion",
        LHB_HEATMAP_MARKERS,
        LHB_HEATMAP_ORDER,
    )

    mhb_assignment.to_csv(out_dir / "panel_c_mhb_l04_to_figure1_subregion.csv", index=False)
    lhb_assignment.to_csv(out_dir / "panel_e_lhb_l04_to_figure1_subregion.csv", index=False)
    pd.concat([mhb_mapped.assign(panel="C/D"), lhb_mapped.assign(panel="E/F")], ignore_index=True).to_csv(
        out_dir / "figure1_marker_genes_mapped_to_h5ad.csv", index=False
    )
    mhb_heat.to_csv(out_dir / "panel_d_mhb_marker_heatmap_values.csv")
    lhb_heat.to_csv(out_dir / "panel_f_lhb_marker_heatmap_values.csv")

    fig = plt.figure(figsize=(16.5, 22.5), constrained_layout=False)
    gs = fig.add_gridspec(3, 2, width_ratios=[1.05, 1.0], height_ratios=[1.0, 1.0, 1.0], wspace=0.30, hspace=0.42)

    ax_a = fig.add_subplot(gs[0, 0])
    broad_label_offsets = {
        "Polydendrocytes": (-1.2, 0.8),
        "Diff. Oligodendrocytes": (1.2, -0.3),
        "Oligodendrocytes": (1.7, 0.7),
        "Fibroblasts": (-0.4, -0.7),
        "Microglia": (0.5, -0.4),
        "Endothelial": (0.8, -0.4),
        "Macrophages": (-0.3, -0.5),
        "Pericytes": (-0.4, -0.5),
    }
    plot_umap(ax_a, plot_df, "gs", BROAD_ORDER, PALETTE, label_fontsize=10, label_offsets=broad_label_offsets)
    plot_panel_label(ax_a, "A")

    ax_b = fig.add_subplot(gs[0, 1])
    plot_dotplot(ax_b, broad_dot)
    plot_panel_label(ax_b, "B")

    ax_c = fig.add_subplot(gs[1, 0])
    mhb_order = [g for g in MHB_HEATMAP_ORDER if g in set(mhb_plot["figure1_mhb_subregion"])]
    plot_umap(ax_c, mhb_plot, "figure1_mhb_subregion", mhb_order, MHB_PALETTE)
    plot_panel_label(ax_c, "C")

    ax_d = fig.add_subplot(gs[1, 1])
    plot_heatmap(ax_d, mhb_heat, "MHb reconstructed subregion markers")
    plot_panel_label(ax_d, "D")

    ax_e = fig.add_subplot(gs[2, 0])
    plot_umap(ax_e, lhb_plot, "figure1_lhb_subregion", LHB_HEATMAP_ORDER, LHB_PALETTE)
    plot_panel_label(ax_e, "E")

    ax_f = fig.add_subplot(gs[2, 1])
    plot_heatmap(ax_f, lhb_heat, "LHb reconstructed subregion markers")
    plot_panel_label(ax_f, "F")

    fig.suptitle("Mouse habenula cell classes and reconstructed subregions", fontsize=18, y=0.995)
    png = out_dir / "figure_1_mouse_habenula_atlas.png"
    pdf = out_dir / "figure_1_mouse_habenula_atlas.pdf"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    return {
        "png": png,
        "pdf": pdf,
        "broad_dotplot_values": out_dir / "panel_b_broad_dotplot_values.csv",
        "mhb_assignment": out_dir / "panel_c_mhb_l04_to_figure1_subregion.csv",
        "lhb_assignment": out_dir / "panel_e_lhb_l04_to_figure1_subregion.csv",
    }


if __name__ == "__main__":
    paths = recreate_figure1()
    for key, value in paths.items():
        print(f"{key}: {value}")
