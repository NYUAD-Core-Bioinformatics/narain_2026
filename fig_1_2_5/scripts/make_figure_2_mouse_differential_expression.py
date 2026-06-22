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
import pandas as pd
import seaborn as sns

INPUT_DIR = Path("inputs/figure_2_mouse_differential_expression")
OUT_DIR = Path("outputs/figure_2_mouse_differential_expression")

PHENO_ORDER = ["si", "susceptible", "resilient", "spt", "tst"]
CELLTYPE_ORDER = [
    "Astrocytes",
    "Endothelial",
    "Fibroblasts",
    "LHb",
    "MHb",
    "Macrophages",
    "Microglia",
    "Oligodendrocytes",
    "Pericytes",
]
LHB_REGION_ORDER = ["Hbx", "Lateral", "Marginal", "Oval/Medial"]
MHB_REGION_ORDER = ["Ventral", "Ventral/Dorsal", "Lateral", "Dorsal", "Superior/Dorsal", "Ventral/Lateral/Dorsal"]


def count_genes(path: Path, group_cols: list[str], out_col: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    gene_col = "names" if "names" in df.columns else "gene"
    out = df.groupby(group_cols, observed=True)[gene_col].nunique().reset_index(name=out_col)
    return out


def plot_heatmap(df: pd.DataFrame, row_col: str, value_col: str, row_order: list[str], out_prefix: str, title: str) -> None:
    table = df.pivot_table(index=row_col, columns="pheno", values=value_col, fill_value=0, aggfunc="sum")
    table = table.reindex(index=[x for x in row_order if x in table.index], columns=[x for x in PHENO_ORDER if x in table.columns])
    fig_w = max(5.5, 0.55 * len(table.columns) + 2.5)
    fig_h = max(3.2, 0.34 * len(table.index) + 1.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(table, ax=ax, cmap="viridis", linewidths=0.4, linecolor="white", cbar_kws={"label": "DE genes"})
    ax.set_title(title)
    ax.set_xlabel("Phenotype")
    ax.set_ylabel("")
    fig.tight_layout()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{out_prefix}.png", dpi=240)
    fig.savefig(OUT_DIR / f"{out_prefix}.pdf")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    celltype_counts = count_genes(
        INPUT_DIR / "celltype_phenotype_bootstrap_de.csv",
        ["ct", "pheno"],
        "n_stable_de_genes",
    )
    lhb_counts = count_genes(
        INPUT_DIR / "lhb_region_phenotype_de.csv",
        ["region", "pheno"],
        "n_region_de_genes",
    )
    mhb_counts = count_genes(
        INPUT_DIR / "mhb_region_phenotype_de.csv",
        ["region", "pheno"],
        "n_region_de_genes",
    )

    celltype_counts.to_csv(OUT_DIR / "celltype_phenotype_de_gene_counts.csv", index=False)
    lhb_counts.to_csv(OUT_DIR / "lhb_region_phenotype_de_gene_counts.csv", index=False)
    mhb_counts.to_csv(OUT_DIR / "mhb_region_phenotype_de_gene_counts.csv", index=False)

    plot_heatmap(celltype_counts, "ct", "n_stable_de_genes", CELLTYPE_ORDER, "celltype_phenotype_de_gene_counts", "Cell-type phenotype DE genes")
    plot_heatmap(lhb_counts, "region", "n_region_de_genes", LHB_REGION_ORDER, "lhb_region_phenotype_de_gene_counts", "LHb region phenotype DE genes")
    plot_heatmap(mhb_counts, "region", "n_region_de_genes", MHB_REGION_ORDER, "mhb_region_phenotype_de_gene_counts", "MHb region phenotype DE genes")


if __name__ == "__main__":
    main()
