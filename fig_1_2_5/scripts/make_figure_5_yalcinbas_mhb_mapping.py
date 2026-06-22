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


LHB_H5AD = Path("inputs/clean_objects/lateral_habenula_regions.h5ad")
MHB_H5AD = Path("inputs/clean_objects/medial_habenula_regions.h5ad")
HUMAN_MARKERS = Path("inputs/figure_5_cross_species/yalcinbas_single_nucleus_marker_genes.csv")
HUMAN_BULK_DE = Path("inputs/figure_5_cross_species/yalcinbas_bulk_schizophrenia_ranked_genes.csv")

MHB_REGION_PHENOTYPE_DE = Path("inputs/figure_2_mouse_differential_expression/mhb_region_phenotype_de.csv")
OUT_DIR = Path("outputs/figure_5_cross_species")
SUPP_DIR = OUT_DIR / "yalcinbas_marker_mapping"
HUMAN_OUT_DIR = OUT_DIR / "yalcinbas_marker_mapping"
MHB_DE_OUT_DIR = OUT_DIR.parent / "figure_2_mouse_differential_expression"

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
PHENO_ORDER = ["si", "susceptible", "resilient", "spt", "tst"]
HUMAN_LHB_POPULATIONS = ["LHb.1", "LHb.2", "LHb.3", "LHb.4", "LHb.5", "LHb.6", "LHb.7"]
HUMAN_MHB_POPULATIONS = ["MHb.1", "MHb.2", "MHb.3"]
TOP_MARKERS_FOR_DOTPLOT = 6
HUMAN_DISEASE_TOP_N = 173


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


def gene_indices(var_names: list[str], genes: list[str]) -> dict[str, int]:
    lookup = {gene: i for i, gene in enumerate(var_names)}
    return {gene: lookup[gene] for gene in genes if gene in lookup}


def load_human_markers() -> pd.DataFrame:
    markers = pd.read_csv(HUMAN_MARKERS)
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


def map_human_to_mouse(markers: pd.DataFrame, var_names: list[str], out_path: Path) -> pd.DataFrame:
    mouse_lookup = pd.DataFrame({"mouse_gene": var_names})
    mouse_lookup["gene_upper"] = mouse_lookup["mouse_gene"].str.upper()
    mouse_lookup = mouse_lookup.drop_duplicates("gene_upper")
    mapped = markers.merge(mouse_lookup, on="gene_upper", how="inner")
    mapped = mapped.sort_values(["human_cluster", "adj.P.Val", "P.Value", "logFC"], ascending=[True, True, True, False])
    mapped.to_csv(out_path, index=False)
    return mapped


def top_human_marker_sets(mapped: pd.DataFrame, clusters: list[str], n: int) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for cluster in clusters:
        sub = mapped[mapped["human_cluster"].eq(cluster) & mapped["logFC"].gt(0)].drop_duplicates("gene_upper")
        out[cluster] = sub.head(n)["mouse_gene"].tolist()
    return out


def load_lhb_for_dotplot() -> tuple[pd.DataFrame, sparse.csr_matrix, sparse.csr_matrix, list[str], list[str]]:
    a = ad.read_h5ad(LHB_H5AD)
    obs = a.obs.copy()
    obs["mouse_region"] = obs["cl"].astype(str).replace({"HbX": "Hbx", "Oval-Medial": "Oval/Medial"})
    categories = [region for region in LHB_REGION_ORDER if region in set(obs["mouse_region"])]
    obs["mouse_region"] = pd.Categorical(obs["mouse_region"], categories=categories, ordered=True)
    counts = as_csr(a.layers["counts"])
    norm = as_csr(a.layers["l1p"]) if "l1p" in a.layers else normalize_log1p_counts(counts)
    return obs, counts, norm, a.var_names.astype(str).tolist(), categories


def load_mhb_for_dotplot() -> tuple[pd.DataFrame, sparse.csr_matrix, sparse.csr_matrix, list[str], list[str]]:
    a = ad.read_h5ad(MHB_H5AD)
    obs = a.obs.copy()
    obs["mouse_region"] = obs["leiden"].astype(str).map(MHB_REGION_MAP).fillna("Unknown")
    categories = [region for region in MHB_REGION_ORDER if region in set(obs["mouse_region"])]
    obs = obs[obs["mouse_region"].isin(categories)].copy()
    keep = np.asarray(a.obs_names.isin(obs.index))
    counts = as_csr(a.layers["counts"])[keep, :].tocsr()
    norm = normalize_log1p_counts(counts)
    obs["mouse_region"] = pd.Categorical(obs["mouse_region"], categories=categories, ordered=True)
    return obs, counts, norm, a.var_names.astype(str).tolist(), categories


def dotplot_values(
    obs: pd.DataFrame,
    counts: sparse.csr_matrix,
    norm: sparse.csr_matrix,
    var_names: list[str],
    categories: list[str],
    human_sets: dict[str, list[str]],
    human_clusters: list[str],
) -> tuple[pd.DataFrame, list[str]]:
    selected_rows = []
    selected_genes = []
    for human_cluster in human_clusters:
        for gene in human_sets.get(human_cluster, [])[:TOP_MARKERS_FOR_DOTPLOT]:
            if gene not in selected_genes:
                selected_genes.append(gene)
                selected_rows.append((human_cluster, gene))

    idx = gene_indices(var_names, selected_genes)
    selected_genes = [gene for gene in selected_genes if gene in idx]
    selected_rows = [(cluster, gene) for cluster, gene in selected_rows if gene in idx]
    positions = [idx[gene] for gene in selected_genes]

    mean_expr = aggregate_by_group(norm[:, positions], obs["mouse_region"], categories)
    pct_expr = aggregate_by_group(counts[:, positions], obs["mouse_region"], categories, binary=True) * 100
    z_expr = (mean_expr - mean_expr.mean(axis=0, keepdims=True)) / (mean_expr.std(axis=0, keepdims=True) + 1e-6)
    z_expr = np.clip(z_expr, -2.5, 2.5)

    gene_to_label = {gene: f"{cluster}: {gene}" for cluster, gene in selected_rows}
    records = []
    for gene_i, gene in enumerate(selected_genes):
        for cluster_i, region in enumerate(categories):
            records.append(
                {
                    "human_marker_gene": gene_to_label[gene],
                    "human_cluster": selected_rows[gene_i][0],
                    "mouse_region": region,
                    "mean_expr_z": z_expr[cluster_i, gene_i],
                    "pct_expr": pct_expr[cluster_i, gene_i],
                }
            )
    dot = pd.DataFrame(records)
    y_labels = [gene_to_label[gene] for gene in selected_genes]
    return dot, y_labels


def draw_dotplot(
    ax: plt.Axes,
    dot: pd.DataFrame,
    y_labels: list[str],
    categories: list[str],
    title: str,
) -> plt.Collection:
    y_pos = {label: i for i, label in enumerate(y_labels)}
    x_pos = {region: i for i, region in enumerate(categories)}
    scatter = ax.scatter(
        dot["mouse_region"].map(x_pos),
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
    ax.set_xlabel("Recovered mouse region")
    ax.set_ylabel("Mapped human marker genes")
    ax.set_title(title)
    ax.grid(axis="x", color="#E7E7E7", linewidth=0.5)
    ax.set_axisbelow(True)
    return scatter


def save_dotplot(
    dot: pd.DataFrame,
    y_labels: list[str],
    categories: list[str],
    title: str,
    out_prefix: str,
) -> dict[str, Path]:
    fig, ax = plt.subplots(figsize=(8, max(5.5, len(y_labels) * 0.24)))
    scatter = draw_dotplot(ax, dot, y_labels, categories, title)
    cbar = fig.colorbar(scatter, ax=ax, pad=0.01)
    cbar.set_label("Mean expression z-score across mouse regions")
    for size in [10, 30, 60]:
        ax.scatter([], [], s=size * 3.0, c="lightgray", edgecolor="0.35", label=f"{size}%")
    ax.legend(title="Pct expressed", bbox_to_anchor=(1.18, 1.0), loc="upper left", frameon=False)
    fig.tight_layout()
    png = SUPP_DIR / f"{out_prefix}.png"
    pdf = SUPP_DIR / f"{out_prefix}.pdf"
    fig.savefig(png, dpi=240)
    fig.savefig(pdf)
    plt.close(fig)
    return {"png": png, "pdf": pdf}


def make_supplementary_dotplots() -> dict[str, Path]:
    SUPP_DIR.mkdir(parents=True, exist_ok=True)
    markers = load_human_markers()

    lhb_obs, lhb_counts, lhb_norm, lhb_var, lhb_categories = load_lhb_for_dotplot()
    lhb_mapped = map_human_to_mouse(markers, lhb_var, SUPP_DIR / "yalcinbas_lhb_supplement_markers_mapped_to_mouse.csv")
    lhb_sets = top_human_marker_sets(lhb_mapped, HUMAN_LHB_POPULATIONS, TOP_MARKERS_FOR_DOTPLOT)
    lhb_dot, lhb_y = dotplot_values(lhb_obs, lhb_counts, lhb_norm, lhb_var, lhb_categories, lhb_sets, HUMAN_LHB_POPULATIONS)
    lhb_dot.to_csv(SUPP_DIR / "yalcinbas_lhb_supplement_marker_dotplot_values.csv", index=False)
    lhb_paths = save_dotplot(
        lhb_dot,
        lhb_y,
        lhb_categories,
        "Supplementary Analysis A: human LHb marker genes across mouse LHb regions",
        "yalcinbas_lhb_marker_dotplot_supplement",
    )

    mhb_obs, mhb_counts, mhb_norm, mhb_var, mhb_categories = load_mhb_for_dotplot()
    mhb_mapped = map_human_to_mouse(markers, mhb_var, SUPP_DIR / "yalcinbas_mhb_markers_mapped_to_mouse.csv")
    mhb_sets = top_human_marker_sets(mhb_mapped, HUMAN_MHB_POPULATIONS, TOP_MARKERS_FOR_DOTPLOT)
    mhb_dot, mhb_y = dotplot_values(mhb_obs, mhb_counts, mhb_norm, mhb_var, mhb_categories, mhb_sets, HUMAN_MHB_POPULATIONS)
    mhb_dot.to_csv(SUPP_DIR / "yalcinbas_mhb_marker_dotplot_values.csv", index=False)
    mhb_paths = save_dotplot(
        mhb_dot,
        mhb_y,
        mhb_categories,
        "Supplementary Analysis A: human MHb marker genes across mouse MHb regions",
        "yalcinbas_mhb_marker_dotplot",
    )

    fig, axes = plt.subplots(1, 2, figsize=(17, max(9.5, len(lhb_y) * 0.24)), gridspec_kw={"width_ratios": [1.1, 0.95]})
    scatter = draw_dotplot(
        axes[0],
        lhb_dot,
        lhb_y,
        lhb_categories,
        "Human LHb marker genes across recovered mouse LHb regions",
    )
    draw_dotplot(
        axes[1],
        mhb_dot,
        mhb_y,
        mhb_categories,
        "Human MHb marker genes across recovered mouse MHb regions",
    )
    cbar = fig.colorbar(scatter, ax=axes, pad=0.015, fraction=0.025)
    cbar.set_label("Mean expression z-score across mouse regions")
    for size in [10, 30, 60]:
        axes[1].scatter([], [], s=size * 3.0, c="lightgray", edgecolor="0.35", label=f"{size}%")
    axes[1].legend(title="Pct expressed", bbox_to_anchor=(1.02, 1.0), loc="upper left", frameon=False)
    fig.suptitle("Supplementary Human Marker Gene Dot Plots For Analysis A", fontsize=14, y=0.995)
    fig.tight_layout(rect=[0, 0, 0.96, 0.97])
    combined_png = SUPP_DIR / "yalcinbas_lhb_mhb_marker_dotplot.png"
    combined_pdf = SUPP_DIR / "yalcinbas_lhb_mhb_marker_dotplot.pdf"
    fig.savefig(combined_png, dpi=240)
    fig.savefig(combined_pdf)
    plt.close(fig)

    return {
        "lhb_dotplot_png": lhb_paths["png"],
        "lhb_dotplot_pdf": lhb_paths["pdf"],
        "mhb_dotplot_png": mhb_paths["png"],
        "mhb_dotplot_pdf": mhb_paths["pdf"],
        "combined_dotplot_png": combined_png,
        "combined_dotplot_pdf": combined_pdf,
    }


def compute_mhb_region_de(force: bool = False) -> pd.DataFrame:
    MHB_DE_OUT_DIR.mkdir(parents=True, exist_ok=True)
    if MHB_REGION_PHENOTYPE_DE.exists() and not force:
        return pd.read_csv(MHB_REGION_PHENOTYPE_DE)

    a = ad.read_h5ad(MHB_H5AD)
    a.obs["region"] = a.obs["leiden"].astype(str).map(MHB_REGION_MAP).fillna("Unknown")
    keep_regions = [region for region in MHB_REGION_ORDER if region in set(a.obs["region"])]
    a = a[a.obs["region"].isin(keep_regions)].copy()

    records = []
    for region in keep_regions:
        sub = a[a.obs["region"].astype(str).eq(region)].copy()
        phenotypes = [pheno for pheno in PHENO_ORDER if pheno in set(sub.obs["pheno"])]
        for pheno in phenotypes:
            if "control" not in set(sub.obs["pheno"]):
                continue
            adx = sub[sub.obs["pheno"].isin(["control", pheno])].copy()
            n_control = int(np.sum(adx.obs["pheno"].astype(str).eq("control")))
            n_pheno = int(np.sum(adx.obs["pheno"].astype(str).eq(pheno)))
            if n_control == 0 or n_pheno == 0:
                continue
            adx.X = adx.layers["counts"].copy()
            sc.pp.normalize_total(adx, target_sum=1e4)
            adx.uns.pop("log1p", None)
            sc.pp.log1p(adx)
            key = f"mhb_{region.replace('/', '_')}_{pheno}"
            sc.tl.rank_genes_groups(
                adx,
                groupby="pheno",
                method="wilcoxon",
                reference="control",
                groups=[pheno],
                rankby_abs=True,
                pts=True,
                key_added=key,
            )
            df = sc.get.rank_genes_groups_df(adx, group=pheno, key=key)
            df["pheno"] = pheno
            df["lfc"] = np.abs(df["logfoldchanges"])
            df = df[df["pvals_adj"].le(0.05)].copy()
            df["log_adj_pval"] = -1 * np.log(np.maximum(df["pvals_adj"], 1e-300))
            df["region"] = region
            df["phenotype"] = pheno
            df["n_control_cells"] = n_control
            df["n_pheno_cells"] = n_pheno
            records.append(df)

    if records:
        de = pd.concat(records, ignore_index=True)
    else:
        de = pd.DataFrame(
            columns=[
                "names",
                "scores",
                "logfoldchanges",
                "pvals",
                "pvals_adj",
                "pct_nz_group",
                "pct_nz_reference",
                "pheno",
                "lfc",
                "log_adj_pval",
                "region",
                "phenotype",
                "n_control_cells",
                "n_pheno_cells",
            ]
        )

    preferred = [
        "names",
        "scores",
        "logfoldchanges",
        "pvals",
        "pvals_adj",
        "pct_nz_group",
        "pct_nz_reference",
        "pheno",
        "lfc",
        "log_adj_pval",
        "region",
        "phenotype",
        "n_control_cells",
        "n_pheno_cells",
    ]
    cols = [col for col in preferred if col in de.columns] + [col for col in de.columns if col not in preferred]
    de = de[cols].sort_values(["region", "pheno", "pvals_adj", "pvals"], ascending=[True, True, True, True])

    counts = (
        de.groupby(["region", "pheno"], observed=True)
        .agg(n_region_de_genes=("names", "nunique"), n_rows=("names", "size"))
        .reset_index()
    )
    counts.to_csv(MHB_DE_OUT_DIR / "mhb_region_phenotype_de_gene_counts.csv", index=False)
    return de


def human_scz_gene_sets() -> tuple[set[str], set[str], set[str], pd.DataFrame]:
    human = pd.read_csv(HUMAN_BULK_DE)
    human["human_gene"] = human["Symbol"].fillna(human["MGI_Symbol"]).astype(str)
    human["gene_upper"] = human["human_gene"].str.upper()
    human = human[human["human_gene"].notna() & human["human_gene"].ne("nan")].drop_duplicates("gene_upper")
    ranked = human.sort_values("P.Value").head(HUMAN_DISEASE_TOP_N).copy()
    all_set = set(ranked["gene_upper"])
    up_set = set(ranked.loc[ranked["logFC"].gt(0), "gene_upper"])
    down_set = set(ranked.loc[ranked["logFC"].lt(0), "gene_upper"])
    return all_set, up_set, down_set, human


def overlap_mouse_de(
    mouse_de: pd.DataFrame,
    group_cols: list[str],
    gene_col: str,
    logfc_col: str,
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
        a_count = len(overlap)
        b_count = len(mouse_all - human_all)
        c_count = len(human_all - mouse_all)
        d_count = max(len(universe) - a_count - b_count - c_count, 0)
        odds_ratio, p_value = stats.fisher_exact([[a_count, b_count], [c_count, d_count]], alternative="greater")
        row = {col: value for col, value in zip(group_cols, key)}
        row.update(
            {
                "n_mouse_de_genes": len(mouse_all),
                "n_human_scz_genes": len(human_all),
                "n_overlap": a_count,
                "n_concordant_direction": len(concordant),
                "n_discordant_direction": len(discordant),
                "odds_ratio": odds_ratio,
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
    return out.sort_values(["p_adj", "p_value", "n_overlap"], ascending=[True, True, False])


def plot_mhb_region_disease_overlap(overlap: pd.DataFrame) -> dict[str, Path]:
    heat = overlap.pivot(index="region", columns="pheno", values="neg_log10_p_adj")
    heat = heat.reindex([region for region in MHB_REGION_ORDER if region in heat.index], columns=PHENO_ORDER).fillna(0)
    annot = overlap.copy()
    annot["label"] = annot.apply(lambda row: f"{int(row.n_overlap)}*" if row.p_adj < 0.05 else str(int(row.n_overlap)), axis=1)
    labels = annot.pivot(index="region", columns="pheno", values="label").reindex(index=heat.index, columns=heat.columns).fillna("")

    fig, ax = plt.subplots(figsize=(8, 4.8))
    vmax = max(1.0, float(np.nanmax(heat.to_numpy(dtype=float))))
    sns.heatmap(
        heat,
        cmap="rocket_r",
        vmin=0,
        vmax=vmax,
        linewidths=0.35,
        linecolor="white",
        annot=labels,
        fmt="",
        cbar_kws={"label": "-log10 FDR, Fisher overlap"},
        ax=ax,
    )
    ax.set_title("Analysis B: recovered MHb subregion phenotype DE overlap with Yalcinbas SCZ Hb genes")
    ax.set_xlabel("Mouse phenotype DE list")
    ax.set_ylabel("Recovered mouse MHb region")
    ax.tick_params(axis="x", labelrotation=35)
    fig.tight_layout()
    png = HUMAN_OUT_DIR / "yalcinbas_bulk_schizophrenia_mhb_region_overlap_heatmap.png"
    pdf = HUMAN_OUT_DIR / "yalcinbas_bulk_schizophrenia_mhb_region_overlap_heatmap.pdf"
    fig.savefig(png, dpi=240)
    fig.savefig(pdf)
    plt.close(fig)
    return {"png": png, "pdf": pdf}


def run_mhb_region_disease_overlap(mhb_de: pd.DataFrame) -> pd.DataFrame:
    HUMAN_OUT_DIR.mkdir(parents=True, exist_ok=True)
    overlap = overlap_mouse_de(
        mhb_de.rename(columns={"names": "mouse_gene", "logfoldchanges": "logFC"}),
        ["region", "pheno"],
        "mouse_gene",
        "logFC",
    )
    overlap.to_csv(HUMAN_OUT_DIR / "yalcinbas_bulk_schizophrenia_mhb_region_overlap.csv", index=False)
    plot_mhb_region_disease_overlap(overlap)
    return overlap


def main(force_de: bool = False) -> dict[str, object]:
    SUPP_DIR.mkdir(parents=True, exist_ok=True)
    HUMAN_OUT_DIR.mkdir(parents=True, exist_ok=True)
    MHB_DE_OUT_DIR.mkdir(parents=True, exist_ok=True)

    mhb_de = compute_mhb_region_de(force=force_de)
    mhb_overlap = run_mhb_region_disease_overlap(mhb_de)
    dotplot_paths = make_supplementary_dotplots()

    print(f"Loaded MHb region phenotype DE table ({mhb_de.shape[0]} rows)")
    print(f"Wrote {HUMAN_OUT_DIR / 'yalcinbas_bulk_schizophrenia_mhb_region_overlap.csv'}")
    print(f"Wrote {dotplot_paths['combined_dotplot_png']}")
    return {
        "mhb_de_rows": int(mhb_de.shape[0]),
        "mhb_de_rows": int(mhb_de.shape[0]),
        "mhb_overlap": HUMAN_OUT_DIR / "yalcinbas_bulk_schizophrenia_mhb_region_overlap.csv",
        **dotplot_paths,
    }


if __name__ == "__main__":
    main(force_de=False)
