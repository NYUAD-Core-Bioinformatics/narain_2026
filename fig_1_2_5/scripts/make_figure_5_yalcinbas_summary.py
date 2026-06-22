from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib_cache")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/xdg_cache")
for key in ["MPLCONFIGDIR", "XDG_CACHE_HOME"]:
    Path(os.environ[key]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(".")
HUMAN_COMPARISON_DIR = ROOT / "outputs" / "figure_5_cross_species" / "yalcinbas_marker_mapping"
CALIBRATION_DIR = ROOT / "outputs" / "figure_5_cross_species" / "yalcinbas_sensitivity_analysis"
OUT_DIR = ROOT / "outputs" / "figure_5_cross_species" / "yalcinbas_summary"
HUMAN_BULK_DE = ROOT / "inputs" / "figure_5_cross_species" / "yalcinbas_bulk_schizophrenia_ranked_genes.csv"
LHB_OVERALL_DE = ROOT / "inputs" / "figure_2_mouse_differential_expression" / "lhb_region_markers.csv"
MHB_OVERALL_DE = ROOT / "inputs" / "figure_2_mouse_differential_expression" / "mhb_region_markers.csv"
MHB_REGION_OVERLAP = HUMAN_COMPARISON_DIR / "yalcinbas_bulk_schizophrenia_mhb_region_overlap.csv"

HUMAN_DISEASE_TOP_N = 173

LHB_REGION_ORDER = ["Hbx", "Lateral", "Marginal", "Oval/Medial"]
MHB_REGION_ORDER = ["Ventral", "Ventral/Dorsal", "Lateral", "Dorsal", "Superior/Dorsal"]
ROW_ORDER = [f"LHb {region}" for region in LHB_REGION_ORDER] + [f"MHb {region}" for region in MHB_REGION_ORDER]
HUMAN_CLUSTER_ORDER = ["LHb.1", "LHb.2", "LHb.3", "LHb.4", "LHb.5", "LHb.6", "LHb.7", "MHb.1", "MHb.2", "MHb.3"]
PHENO_ORDER = ["si", "susceptible", "resilient", "spt", "tst"]
BROAD_CELLTYPE_ORDER = [
    "Endothelial",
    "Pericytes",
    "Astrocytes",
    "Fibroblasts",
    "Macrophages",
    "Microglia",
    "Oligodendrocytes",
    "LHb",
    "MHb",
]


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
        a = len(overlap)
        b = len(mouse_all - human_all)
        c = len(human_all - mouse_all)
        d = max(len(universe) - a - b - c, 0)
        odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        row = {col: value for col, value in zip(group_cols, key)}
        row.update(
            {
                "n_mouse_de_genes": len(mouse_all),
                "n_human_scz_genes": len(human_all),
                "n_overlap": a,
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


def build_hb_compartment_overlap() -> pd.DataFrame:
    lhb = pd.read_csv(LHB_OVERALL_DE).rename(columns={"names": "mouse_gene", "logfoldchanges": "logFC"})
    lhb["compartment"] = "LHb overall"
    mhb = pd.read_csv(MHB_OVERALL_DE).rename(columns={"names": "mouse_gene", "logfoldchanges": "logFC"})
    mhb["compartment"] = "MHb overall"
    de = pd.concat([lhb, mhb], ignore_index=True)
    out = overlap_mouse_de(de, ["compartment", "pheno"], "mouse_gene", "logFC")
    out.to_csv(OUT_DIR / "yalcinbas_bulk_schizophrenia_habenula_compartment_overlap.csv", index=False)
    return out


def add_row_label(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["row_label"] = df["compartment"] + " " + df["mouse_region"].astype(str)
    return df


def load_analysis_a_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    overlap_parts = []
    module_parts = []
    for compartment in ["LHb", "MHb"]:
        prefix = compartment.lower()
        allowed_prefix = f"{compartment}."
        overlap = pd.read_csv(CALIBRATION_DIR / f"{prefix}_overlap_threshold_sensitivity_summary.csv")
        overlap = overlap[overlap["human_cluster"].astype(str).str.startswith(allowed_prefix)].copy()
        overlap_parts.append(add_row_label(overlap))

        module = pd.read_csv(CALIBRATION_DIR / f"{prefix}_module_depth_sensitivity_summary.csv")
        module = module[module["human_cluster"].astype(str).str.startswith(allowed_prefix)].copy()
        module_parts.append(add_row_label(module))

    overlap = pd.concat(overlap_parts, ignore_index=True)
    module = pd.concat(module_parts, ignore_index=True)
    summary = pd.read_csv(CALIBRATION_DIR / "lhb_mhb_same_compartment_calibrated_best_match_summary.csv")
    summary["row_label"] = summary["compartment"] + " " + summary["mouse_region"].astype(str)
    summary.to_csv(OUT_DIR / "yalcinbas_lhb_mhb_top_region_matches.csv", index=False)
    return overlap, module, summary


def pivot_metric(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    heat = df.pivot(index="row_label", columns="human_cluster", values=value_col)
    return heat.reindex(index=ROW_ORDER, columns=HUMAN_CLUSTER_ORDER)


def format_disease_annotations(df: pd.DataFrame, row_col: str, row_order: list[str]) -> pd.DataFrame:
    ann = df.copy()
    ann["label"] = ann.apply(
        lambda row: f"{int(row.n_overlap)}*" if row.p_adj < 0.05 else str(int(row.n_overlap)),
        axis=1,
    )
    heat = ann.pivot(index=row_col, columns="pheno", values="label")
    return heat.reindex(index=row_order, columns=PHENO_ORDER).fillna("")


def plot_analysis_a_heatmap(
    ax: plt.Axes,
    heat: pd.DataFrame,
    title: str,
    cmap: str,
    cbar_label: str,
    center: float | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
) -> None:
    annot = heat.apply(lambda col: col.map(lambda value: "" if pd.isna(value) else f"{value:.2f}"))
    sns.heatmap(
        heat,
        mask=heat.isna(),
        cmap=cmap,
        center=center,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.35,
        linecolor="white",
        annot=annot,
        fmt="",
        cbar_kws={"label": cbar_label},
        ax=ax,
    )
    ax.axhline(len(LHB_REGION_ORDER), color="#333333", linewidth=1.2)
    ax.axvline(7, color="#333333", linewidth=1.2)
    ax.set_title(title, fontsize=12, pad=8)
    ax.set_xlabel("Human Yalcinbas Hb population")
    ax.set_ylabel("Mouse recovered region")
    ax.tick_params(axis="x", labelrotation=45)
    ax.tick_params(axis="y", labelrotation=0)


def plot_top_match_table(ax: plt.Axes, summary: pd.DataFrame) -> None:
    ax.axis("off")
    summary = summary.set_index("row_label").reindex(ROW_ORDER).reset_index()
    rows = []
    for row in summary.itertuples(index=False):
        rows.append(
            [
                row.row_label,
                f"{row.best_overlap_human_cluster} ({row.overlap_best_fraction:.2f})",
                f"{row.best_module_human_cluster} ({row.module_best_fraction:.2f})",
                f"{row.best_positive_rankcorr_human_cluster} (r={row.positive_rankcorr_spearman_r:.2f})",
            ]
        )
    table = ax.table(
        cellText=rows,
        colLabels=["Mouse region", "Overlap best", "Module best", "Rank-corr best"],
        cellLoc="left",
        colLoc="left",
        loc="center",
        colWidths=[0.28, 0.22, 0.22, 0.28],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1, 1.35)
    for (row_idx, _col_idx), cell in table.get_celld().items():
        cell.set_edgecolor("#DADADA")
        if row_idx == 0:
            cell.set_facecolor("#F0F0F0")
            cell.set_text_props(weight="bold")
        elif row_idx <= len(LHB_REGION_ORDER):
            cell.set_facecolor("#F8FBFD")
        else:
            cell.set_facecolor("#FDF8F5")
    ax.set_title("Top calibrated human match per mouse subtype", fontsize=12, pad=8)


def plot_disease_heatmap(
    ax: plt.Axes,
    df: pd.DataFrame,
    row_col: str,
    row_order: list[str],
    title: str,
    vmax: float,
) -> None:
    heat = df.pivot(index=row_col, columns="pheno", values="neg_log10_p_adj")
    heat = heat.reindex(index=row_order, columns=PHENO_ORDER).fillna(0)
    annot = format_disease_annotations(df, row_col, row_order)
    sns.heatmap(
        heat,
        cmap="rocket_r",
        vmin=0,
        vmax=vmax,
        linewidths=0.35,
        linecolor="white",
        annot=annot,
        fmt="",
        cbar_kws={"label": "-log10 FDR"},
        ax=ax,
    )
    ax.set_title(title, fontsize=12, pad=8)
    ax.set_xlabel("Mouse phenotype DE list")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=35)
    ax.tick_params(axis="y", labelrotation=0)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=15,
        fontweight="bold",
        va="top",
        ha="right",
    )


def make_comprehensive_figure() -> dict[str, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="white", font_scale=0.9)

    overlap, module, summary = load_analysis_a_tables()
    broad = pd.read_csv(HUMAN_COMPARISON_DIR / "yalcinbas_bulk_schizophrenia_celltype_overlap.csv")
    lhb_region = pd.read_csv(HUMAN_COMPARISON_DIR / "yalcinbas_bulk_schizophrenia_lhb_region_overlap.csv")
    lhb_region["region"] = lhb_region["region"].replace({"Oval-Medial": "Oval/Medial", "HbX": "Hbx"})
    hb_compartment = build_hb_compartment_overlap()
    if MHB_REGION_OVERLAP.exists():
        mhb_panel = pd.read_csv(MHB_REGION_OVERLAP)
        mhb_panel_row_col = "region"
        mhb_panel_order = [region for region in MHB_REGION_ORDER if region in set(mhb_panel["region"])]
        mhb_panel_title = "Analysis B: recovered MHb subregion phenotype DE overlap"
    else:
        mhb_panel = hb_compartment
        mhb_panel_row_col = "compartment"
        mhb_panel_order = ["LHb overall", "MHb overall"]
        mhb_panel_title = "Analysis B: LHb/MHb phenotype DE overlap"

    broad_order = [cell for cell in BROAD_CELLTYPE_ORDER if cell in set(broad["celltype"])]
    lhb_order = [region for region in LHB_REGION_ORDER if region in set(lhb_region["region"])]
    disease_vmax = max(
        1.0,
        float(broad["neg_log10_p_adj"].max()),
        float(lhb_region["neg_log10_p_adj"].max()),
        float(hb_compartment["neg_log10_p_adj"].max()),
        float(mhb_panel["neg_log10_p_adj"].max()),
    )

    overlap_heat = pivot_metric(overlap, "best_match_fraction")
    module_heat = pivot_metric(module, "mean_score_across_depths")
    module_abs = np.nanmax(np.abs(module_heat.to_numpy(dtype=float)))
    module_abs = max(0.05, float(module_abs))

    fig = plt.figure(figsize=(22, 19), constrained_layout=False)
    grid = fig.add_gridspec(
        nrows=3,
        ncols=2,
        height_ratios=[1.12, 0.95, 0.9],
        width_ratios=[1.05, 1.0],
        hspace=0.44,
        wspace=0.24,
    )

    ax_a = fig.add_subplot(grid[0, 0])
    plot_analysis_a_heatmap(
        ax_a,
        overlap_heat,
        "Analysis A: marker-overlap best-match stability",
        cmap="mako",
        cbar_label="fraction of threshold grid",
        vmin=0,
        vmax=1,
    )

    ax_b = fig.add_subplot(grid[0, 1])
    plot_analysis_a_heatmap(
        ax_b,
        module_heat,
        "Analysis A: human marker module-score conservation",
        cmap="vlag",
        cbar_label="mean z-scored module score",
        center=0,
        vmin=-module_abs,
        vmax=module_abs,
    )

    ax_c = fig.add_subplot(grid[1, 0])
    plot_top_match_table(ax_c, summary)

    ax_d = fig.add_subplot(grid[1, 1])
    plot_disease_heatmap(
        ax_d,
        broad,
        "celltype",
        broad_order,
        "Analysis B: broad mouse cell-type DE overlap with Yalcinbas SCZ Hb genes",
        vmax=disease_vmax,
    )

    ax_e = fig.add_subplot(grid[2, 0])
    plot_disease_heatmap(
        ax_e,
        lhb_region,
        "region",
        lhb_order,
        "Analysis B: recovered LHb subregion phenotype DE overlap",
        vmax=disease_vmax,
    )

    ax_f = fig.add_subplot(grid[2, 1])
    plot_disease_heatmap(
        ax_f,
        mhb_panel,
        mhb_panel_row_col,
        mhb_panel_order,
        mhb_panel_title,
        vmax=disease_vmax,
    )

    for ax, label in zip([ax_a, ax_b, ax_c, ax_d, ax_e, ax_f], list("ABCDEF")):
        add_panel_label(ax, label)

    fig.suptitle(
        "Yalcinbas Human Habenula Comparison: Mouse LHb/MHb Conservation And Disease-Overlap Summary",
        fontsize=17,
        y=0.985,
    )
    fig.text(
        0.5,
        0.018,
        "Analysis A uses mapped human snRNA-seq Hb markers after mouse-human symbol matching. "
        "Analysis B tests overlap with the top 173 Yalcinbas bulk schizophrenia-ranked Hb genes; heatmap numbers are overlapping genes and * marks FDR < 0.05. "
        "MHb subregion phenotype DE was rebuilt with the same per-cell Wilcoxon phenotype-vs-control workflow used for the LHb region DE table.",
        ha="center",
        va="bottom",
        fontsize=9.5,
    )

    png = OUT_DIR / "figure_5_yalcinbas_summary.png"
    pdf = OUT_DIR / "figure_5_yalcinbas_summary.pdf"
    fig.savefig(png, dpi=240, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    return {"png": png, "pdf": pdf}


def main() -> dict[str, Path]:
    paths = make_comprehensive_figure()
    print(f"Wrote {paths['png']}")
    print(f"Wrote {paths['pdf']}")
    return paths


if __name__ == "__main__":
    main()
