# ============================================================
# Imports
# ============================================================

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import tqdm
from joblib import Parallel, delayed
from venn import venn

# ============================================================
# Differential expression: phenotype vs control (whole LHb)
# ============================================================


def differential_expression_by_phenotype(
    adata: sc.AnnData,
    pheno_key: str = "pheno",
    control_label: str = "control",
    layer: str = "counts",
    target_sum: float = 1e4,
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    phenotypes = [p for p in adata.obs[pheno_key].unique() if p != control_label]
    de_tables = []

    for ph in phenotypes:
        sub = adata[adata.obs[pheno_key].isin([control_label, ph])].copy()

        if layer is not None:
            sub.X = sub.layers[layer].copy()

        sc.pp.normalize_total(sub, target_sum=target_sum)
        sc.pp.log1p(sub)

        key = f"de_{ph}"
        sc.tl.rank_genes_groups(
            sub,
            groupby=pheno_key,
            method="wilcoxon",
            reference=control_label,
            rankby_abs=True,
            pts=True,
            key_added=key,
        )

        df = sc.get.rank_genes_groups_df(sub, group=ph, key=key)
        df["phenotype"] = ph
        de_tables.append(df)

    dedf = pd.concat(de_tables, ignore_index=True)
    dedf["lfc"] = dedf["logfoldchanges"].abs()
    dedf = dedf[dedf["pvals_adj"] <= fdr_threshold].copy()
    dedf["log_adj_pval"] = -np.log(dedf["pvals_adj"])

    return dedf


# ============================================================
# Differential expression: phenotype vs control within regions
# ============================================================


def differential_expression_by_region(
    adata: sc.AnnData,
    region_key: str = "cl",
    pheno_key: str = "pheno",
    control_label: str = "control",
    layer: str = "counts",
    target_sum: float = 1e4,
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    results = []

    for region in adata.obs[region_key].unique():
        sub_region = adata[adata.obs[region_key] == region].copy()
        phenotypes = [
            p for p in sub_region.obs[pheno_key].unique() if p != control_label
        ]

        for ph in phenotypes:
            sub = sub_region[sub_region.obs[pheno_key].isin([control_label, ph])].copy()

            if layer is not None:
                sub.X = sub.layers[layer].copy()

            sc.pp.normalize_total(sub, target_sum=target_sum)
            sc.pp.log1p(sub)

            key = f"de_{region}_{ph}"
            sc.tl.rank_genes_groups(
                sub,
                groupby=pheno_key,
                method="wilcoxon",
                reference=control_label,
                rankby_abs=True,
                pts=True,
                key_added=key,
            )

            df = sc.get.rank_genes_groups_df(sub, group=ph, key=key)
            df["phenotype"] = ph
            df["region"] = region
            results.append(df)

    dedf = pd.concat(results, ignore_index=True)
    dedf["lfc"] = dedf["logfoldchanges"].abs()
    dedf = dedf[dedf["pvals_adj"] <= fdr_threshold].copy()
    dedf["log_adj_pval"] = -np.log(dedf["pvals_adj"])

    return dedf


# ============================================================
# Bootstrap DE controlling for cell number differences
# ============================================================


def _run_one_bootstrap(
    ct: sc.AnnData,
    pheno_key: str,
    control_label: str,
    conditions: np.ndarray,
    n_sample: int,
    iteration: int,
    layer: str = "counts",
    target_sum: float = 1e4,
) -> pd.DataFrame:
    rng = np.random.RandomState(iteration)
    sampled_indices = []

    for cond in conditions:
        idx = ct.obs.index[ct.obs[pheno_key] == cond].values
        sampled = rng.choice(idx, size=n_sample, replace=False)
        sampled_indices.append(sampled)

    sampled_indices = np.concatenate(sampled_indices)
    sub = ct[sampled_indices].copy()

    if layer is not None:
        sub.X = sub.layers[layer].copy()

    sc.pp.normalize_total(sub, target_sum=target_sum)
    sc.pp.log1p(sub)

    key = f"boot_{iteration}"
    sc.tl.rank_genes_groups(
        sub,
        groupby=pheno_key,
        method="wilcoxon",
        reference=control_label,
        rankby_abs=True,
        pts=True,
        key_added=key,
    )

    dfs = []
    for ph in conditions:
        if ph == control_label:
            continue
        df = sc.get.rank_genes_groups_df(sub, group=ph, key=key)
        df["phenotype"] = ph
        df["iteration"] = iteration
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


def bootstrap_de_by_celltype(
    adata: sc.AnnData,
    celltype_key: str = "subclass_gs",
    pheno_key: str = "pheno",
    control_label: str = "control",
    layer: str = "counts",
    n_bootstrap: int = 10000,
    n_jobs: int = 8,
    fdr_threshold: float = 0.05,
) -> dict:
    results = {}

    for celltype in adata.obs[celltype_key].unique():
        ct = adata[adata.obs[celltype_key] == celltype].copy()
        conditions = ct.obs[pheno_key].unique()

        n_sample = int(ct.obs[pheno_key].value_counts().loc[conditions].min())
        print(f"{celltype}: {n_sample} cells per condition")

        all_iters = Parallel(n_jobs=n_jobs)(
            delayed(_run_one_bootstrap)(
                ct=ct,
                pheno_key=pheno_key,
                control_label=control_label,
                conditions=conditions,
                n_sample=n_sample,
                iteration=i,
                layer=layer,
            )
            for i in tqdm.tqdm(range(n_bootstrap))
        )

        big_df = pd.concat(all_iters, ignore_index=True)
        big_df = big_df[big_df["pvals_adj"] <= fdr_threshold].copy()
        big_df["lfc"] = big_df["logfoldchanges"].abs()
        big_df[celltype_key] = celltype

        results[celltype] = big_df

    return results


# ============================================================
# Bootstrap summarisation and proportions
# ============================================================


def summarize_bootstrap_results(
    df: pd.DataFrame,
    freq_cutoff: float = 0.5,
) -> pd.DataFrame:
    n_iterations = df["iteration"].nunique()

    summary = (
        df.groupby("names")
        .agg(
            freq=("iteration", "nunique"),
            avg_logFC=("logfoldchanges", "mean"),
            std_logFC=("logfoldchanges", "std"),
            avg_pval=("pvals_adj", "mean"),
        )
        .reset_index()
    )

    summary["freq_ratio"] = summary["freq"] / n_iterations
    stable = summary[summary["freq_ratio"] >= freq_cutoff].copy()

    return stable.sort_values("freq_ratio", ascending=False)


def compute_stable_gene_proportions(dedf_boot: dict) -> pd.DataFrame:
    rows = []

    for celltype, df in dedf_boot.items():
        for ph, ph_df in df.groupby("phenotype"):
            stable = summarize_bootstrap_results(ph_df)
            stable = stable[
                ~(
                    stable["names"].str.startswith("Rp")
                    | stable["names"].str.startswith("mt-")
                )
            ]
            rows.append(
                {
                    "celltype": celltype,
                    "phenotype": ph,
                    "count": stable.shape[0],
                }
            )

    df_counts = pd.DataFrame(rows)
    df_counts["total"] = df_counts.groupby("celltype")["count"].transform("sum")
    df_counts["percentage"] = (df_counts["count"] / df_counts["total"]) * 100

    return df_counts


# ============================================================
# Venn diagrams of shared DE genes
# ============================================================


def plot_venn_by_phenotype(dedf_by_pheno: pd.DataFrame, out_file: str):
    sets = {
        ph: set(df["names"].tolist()) for ph, df in dedf_by_pheno.groupby("phenotype")
    }
    plt.figure(figsize=(8, 8))
    venn(sets, fmt="{size}")
    plt.title("DE genes: phenotype vs control (LHb)")
    plt.tight_layout()
    plt.savefig(out_file)
    plt.show()


def plot_venn_by_region(
    dedf_region: pd.DataFrame,
    phenotype: str,
    out_file: str,
):
    df = dedf_region[dedf_region["phenotype"] == phenotype]
    sets = {region: set(sub["names"].tolist()) for region, sub in df.groupby("region")}
    plt.figure(figsize=(8, 8))
    venn(sets, fmt="{size}")
    plt.title(f"{phenotype}: DE genes across LHb regions")
    plt.tight_layout()
    plt.savefig(out_file)
    plt.show()


# ============================================================
# Cell counts and proportions
# ============================================================


def cell_counts_and_proportions(
    adata: sc.AnnData,
    group_keys=("pheno", "cl"),
) -> pd.DataFrame:
    counts = adata.obs.groupby(list(group_keys)).size().reset_index(name="n_cells")

    if len(group_keys) > 1:
        counts["prop"] = counts["n_cells"] / counts.groupby(group_keys[0])[
            "n_cells"
        ].transform("sum")
    else:
        counts["prop"] = counts["n_cells"] / counts["n_cells"].sum()

    return counts
