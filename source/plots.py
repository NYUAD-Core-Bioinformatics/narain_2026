import os

from helpers import (
    bootstrap_de_by_celltype,
    cell_counts_and_proportions,
    compute_stable_gene_proportions,
    differential_expression_by_phenotype,
    differential_expression_by_region,
    plot_venn_by_phenotype,
    plot_venn_by_region,
)

adx = sc.read_h5ad("adx_for_boostraps.h5ad")
lhb = sc.read_h5ad("lhb_apr18.h5ad")
combo_l = sc.read_h5ad("combo_l_may8.h5ad")
mhb = sc.read_h5ad("lhb_apr18.h5ad")
# ============================================================
# MAIN ANALYSIS PIPELINE
# ============================================================

# ----------------------------
# Inputs / assumptions
# ----------------------------
# lhb : AnnData object containing LHb cells only
# adx : AnnData object containing integrated dataset with cell-type labels
# Required obs fields:
#   lhb.obs["pheno"]       – behavioural phenotype
#   lhb.obs["cl"]          – LHb anatomical region
#   adx.obs["pheno"]       – behavioural phenotype
#   adx.obs["subclass_gs"] – cell-type annotation
# All objects must contain raw counts in .layers["counts"]

OUTDIR = "publication_figures"
os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# 1. Differential expression: phenotype vs control (whole LHb)
# ============================================================

dedf_by_pheno = differential_expression_by_phenotype(
    adata=lhb,
    pheno_key="pheno",
    control_label="control",
)

dedf_by_pheno.to_csv(f"{OUTDIR}/lhb_DE_by_phenotype.csv", index=False)


# ============================================================
# 2. Differential expression: phenotype vs control by region
# ============================================================

dedf_pheno_byregion = differential_expression_by_region(
    adata=lhb,
    region_key="cl",
    pheno_key="pheno",
    control_label="control",
)

dedf_pheno_byregion.to_csv(f"{OUTDIR}/lhb_DE_by_phenotype_and_region.csv", index=False)


# ============================================================
# 3. Bootstrap DE (cell-number–controlled, per cell type)
# ============================================================

dedf_boot = bootstrap_de_by_celltype(
    adata=adx,
    celltype_key="subclass_gs",
    pheno_key="pheno",
    control_label="control",
    n_bootstrap=10000,
    n_jobs=14,
)

# Save raw bootstrap results (optional but recommended)
for celltype, df in dedf_boot.items():
    df.to_csv(
        f"{OUTDIR}/bootstrap_DE_{celltype}.csv",
        index=False,
    )


# ============================================================
# 4. Stable DE genes and proportions
# ============================================================

df_stable_props = compute_stable_gene_proportions(dedf_boot)

df_stable_props.to_csv(
    f"{OUTDIR}/stable_DE_gene_proportions.csv",
    index=False,
)


# ============================================================
# 5. Venn diagrams
# ============================================================

# 5a. Across phenotypes (whole LHb)
plot_venn_by_phenotype(
    dedf_by_pheno,
    out_file=f"{OUTDIR}/lhb_phenotype_venn.pdf",
)

# 5b. Across regions, per phenotype
for ph in dedf_pheno_byregion["phenotype"].unique():
    plot_venn_by_region(
        dedf_region=dedf_pheno_byregion,
        phenotype=ph,
        out_file=f"{OUTDIR}/lhb_region_venn_{ph}.pdf",
    )


# ============================================================
# 6. Cell counts and proportions
# ============================================================

cell_comp = cell_counts_and_proportions(
    adata=lhb,
    group_keys=("pheno", "cl"),
)

cell_comp.to_csv(
    f"{OUTDIR}/lhb_cell_counts_by_pheno_and_region.csv",
    index=False,
)


# ============================================================
# 7. (Optional) Whole-dataset cell-type composition
# ============================================================

celltype_counts = (
    adx.obs["subclass_gs"]
    .value_counts()
    .reset_index()
    .rename(columns={"index": "celltype", "subclass_gs": "n_cells"})
)

celltype_counts["prop"] = celltype_counts["n_cells"] / celltype_counts["n_cells"].sum()

celltype_counts.to_csv(
    f"{OUTDIR}/whole_dataset_celltype_composition.csv",
    index=False,
)
