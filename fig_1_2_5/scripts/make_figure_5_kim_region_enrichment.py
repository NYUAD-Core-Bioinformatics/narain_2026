from __future__ import annotations

import os
import re
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

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


THIS_DIR = Path(__file__).resolve().parent
ROOT = THIS_DIR.parent
OUT_DIR = ROOT / "outputs" / "figure_5_cross_species" / "kim_region_enrichment"
FIG_DIR = OUT_DIR / "figures"

KIM_SUPPLEMENT = ROOT / "inputs" / "figure_5_cross_species" / "kim_suicide_habenula_degs.xlsx"
ORTHOLOG_LINK = ROOT / "inputs" / "figure_5_cross_species" / "human_mouse_gene_symbols.csv"
LHB_H5AD = ROOT / "inputs" / "clean_objects" / "lateral_habenula_regions.h5ad"
MHB_H5AD = ROOT / "inputs" / "clean_objects" / "medial_habenula_regions.h5ad"

REGION_ORDER = ["LHb", "MHb"]
MARKER_LOGFC_MIN = 0.25
MARKER_FDR_MAX = 0.05
TOP_DOTPLOT_GENES_PER_REGION = 50


def col_to_index(cell_ref: str) -> int:
    match = re.match(r"([A-Z]+)", cell_ref)
    if match is None:
        return 0
    out = 0
    for char in match.group(1):
        out = out * 26 + ord(char) - 64
    return out - 1


def read_xlsx_sheet(path: Path, sheet_name: str) -> pd.DataFrame:
    ns = {
        "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    }
    with zipfile.ZipFile(path) as zf:
        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        rels_root = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rels = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels_root}

        target = None
        for sheet in workbook.findall(".//a:sheet", ns):
            if sheet.attrib.get("name") == sheet_name:
                rel_id = sheet.attrib[f"{{{ns['r']}}}id"]
                target = rels[rel_id]
                break
        if target is None:
            raise ValueError(f"Sheet not found: {sheet_name}")
        if not target.startswith("xl/"):
            target = f"xl/{target}"

        shared_strings = []
        if "xl/sharedStrings.xml" in zf.namelist():
            ss_root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in ss_root.findall("a:si", ns):
                shared_strings.append("".join(t.text or "" for t in si.findall(".//a:t", ns)))

        sheet_root = ET.fromstring(zf.read(target))
        rows: list[list[str]] = []
        for row in sheet_root.findall(".//a:row", ns):
            values: dict[int, str] = {}
            max_col = -1
            for cell in row.findall("a:c", ns):
                col_idx = col_to_index(cell.attrib["r"])
                max_col = max(max_col, col_idx)
                value_el = cell.find("a:v", ns)
                value = "" if value_el is None or value_el.text is None else value_el.text
                if cell.attrib.get("t") == "s" and value:
                    value = shared_strings[int(value)]
                values[col_idx] = value
            rows.append([values.get(i, "") for i in range(max_col + 1)])

    if len(rows) < 2:
        return pd.DataFrame()
    title = rows[0][0]
    header = rows[1]
    data = rows[2:]
    width = len(header)
    data = [row + [""] * (width - len(row)) for row in data if any(str(x).strip() for x in row)]
    df = pd.DataFrame([row[:width] for row in data], columns=header)
    df.attrs["title"] = title
    return df


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.astype(str).str.strip().replace({"": np.nan}), errors="coerce")


def clean_gene_symbol(value: object) -> str:
    text = str(value).strip()
    text = re.sub(r"\s+", "", text)
    return text


def split_human_symbols(symbol: str) -> list[str]:
    symbol = clean_gene_symbol(symbol)
    if not symbol or symbol.lower() == "nan":
        return []
    if symbol == "HSPA1A/B":
        return ["HSPA1A", "HSPA1B"]
    parts = re.split(r"[;,]", symbol)
    out = []
    for part in parts:
        part = part.strip()
        if part:
            out.append(part)
    return out


def load_kim_degs() -> pd.DataFrame:
    deg = read_xlsx_sheet(KIM_SUPPLEMENT, "Table S2")
    deg = deg.rename(
        columns={
            "Probe set ID": "probe_set_id",
            "Gene symbol": "human_gene",
            "Gene Accession": "gene_accession",
            "Control (mean)a": "control_mean",
            "Suicide (mean)b": "suicide_mean",
            "Ratio (S-to-C)": "ratio_s_to_c",
            "P-value": "p_value",
        }
    )
    deg["human_gene"] = deg["human_gene"].map(clean_gene_symbol)
    for col in ["control_mean", "suicide_mean", "ratio_s_to_c", "p_value"]:
        deg[col] = numeric(deg[col])
    deg["log2_ratio_s_to_c"] = np.log2(deg["ratio_s_to_c"])
    deg["direction"] = np.where(deg["ratio_s_to_c"].ge(1), "up_in_suicide", "down_in_suicide")

    expanded = []
    for row in deg.itertuples(index=False):
        for human_gene in split_human_symbols(row.human_gene):
            record = row._asdict()
            record["human_gene_original"] = row.human_gene
            record["human_gene"] = human_gene
            record["human_gene_upper"] = human_gene.upper()
            expanded.append(record)
    out = pd.DataFrame(expanded)
    out.to_csv(OUT_DIR / "kim_table_s2_human_degs_expanded.csv", index=False)
    return out


def load_kim_original_table_s3() -> pd.DataFrame:
    s3 = read_xlsx_sheet(KIM_SUPPLEMENT, "Table S3")
    rows = []
    current_columns = []
    for i, row in s3.iterrows():
        if i == 0:
            current_columns = [
                "gene_symbol",
                "oligodendrocyte",
                "mhb_neuron",
                "lhb_neuron",
                "microglia",
                "ependymal",
                "astrocyte",
                "opc",
                "mural",
                "endothelial",
                "kim_enriched_cell_type",
            ]
            continue
        values = row.tolist()
        if not values or not str(values[0]).strip():
            continue
        if str(values[0]).strip().startswith("a"):
            continue
        rows.append(values[: len(current_columns)])
    out = pd.DataFrame(rows, columns=current_columns)
    out["gene_symbol"] = out["gene_symbol"].map(clean_gene_symbol)
    for col in current_columns[1:-1]:
        out[col] = numeric(out[col])
    out.to_csv(OUT_DIR / "kim_original_table_s3_mouse_hb_celltype_enrichment.csv", index=False)
    return out


def load_region_data() -> tuple[ad.AnnData, sparse.csr_matrix, list[str]]:
    lhb = ad.read_h5ad(LHB_H5AD)
    mhb = ad.read_h5ad(MHB_H5AD)
    common = sorted(set(lhb.var_names.astype(str)) & set(mhb.var_names.astype(str)))
    lhb_idx = pd.Index(lhb.var_names.astype(str)).get_indexer(common)
    mhb_idx = pd.Index(mhb.var_names.astype(str)).get_indexer(common)

    lhb_counts = sparse.csr_matrix(lhb.layers["counts"])[:, lhb_idx]
    mhb_counts = sparse.csr_matrix(mhb.layers["counts"])[:, mhb_idx]
    counts = sparse.vstack([lhb_counts, mhb_counts], format="csr")
    obs = pd.DataFrame(
        {
            "region": pd.Categorical(["LHb"] * lhb.n_obs + ["MHb"] * mhb.n_obs, categories=REGION_ORDER, ordered=True)
        },
        index=[f"LHb_{idx}" for idx in lhb.obs_names.astype(str)] + [f"MHb_{idx}" for idx in mhb.obs_names.astype(str)],
    )
    a = ad.AnnData(X=counts.copy(), obs=obs, var=pd.DataFrame(index=common))
    a.layers["counts"] = counts.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    a.uns.pop("log1p", None)
    sc.pp.log1p(a)
    norm = sparse.csr_matrix(a.X)
    return a, norm, common


def compute_region_markers(a: ad.AnnData) -> pd.DataFrame:
    sc.tl.rank_genes_groups(
        a,
        groupby="region",
        groups=REGION_ORDER,
        method="wilcoxon",
        reference="rest",
        pts=True,
        rankby_abs=False,
        key_added="region_markers",
    )
    pieces = []
    for region in REGION_ORDER:
        df = sc.get.rank_genes_groups_df(a, group=region, key="region_markers")
        df["region"] = region
        pieces.append(df)
    markers = pd.concat(pieces, ignore_index=True)
    markers = markers.rename(columns={"names": "mouse_gene", "logfoldchanges": "logFC_vs_other_region"})
    markers["mouse_gene_upper"] = markers["mouse_gene"].astype(str).str.upper()
    markers["is_region_marker"] = markers["logFC_vs_other_region"].gt(MARKER_LOGFC_MIN) & markers["pvals_adj"].lt(MARKER_FDR_MAX)
    markers.to_csv(OUT_DIR / "mouse_lhb_mhb_region_markers_wilcoxon.csv", index=False)
    markers.loc[markers["is_region_marker"]].to_csv(OUT_DIR / "mouse_lhb_mhb_region_markers_filtered_logfc0p25_fdr0p05.csv", index=False)
    return markers


def average_expression_table(a: ad.AnnData, norm: sparse.csr_matrix, var_names: list[str]) -> pd.DataFrame:
    rows = []
    for region in REGION_ORDER:
        mask = a.obs["region"].astype(str).to_numpy() == region
        mean_expr = np.asarray(norm[mask, :].mean(axis=0)).ravel()
        pct_expr = np.asarray((a.layers["counts"][mask, :] > 0).mean(axis=0)).ravel()
        rows.append(
            pd.DataFrame(
                {
                    "mouse_gene": var_names,
                    "region": region,
                    "mean_log1p_norm_expr": mean_expr,
                    "pct_expr": pct_expr,
                }
            )
        )
    expr = pd.concat(rows, ignore_index=True)
    wide = expr.pivot(index="mouse_gene", columns="region", values="mean_log1p_norm_expr").reset_index()
    wide.columns.name = None
    wide = wide.rename(columns={"LHb": "avg_expr_LHb", "MHb": "avg_expr_MHb"})
    pct = expr.pivot(index="mouse_gene", columns="region", values="pct_expr").reset_index()
    pct.columns.name = None
    pct = pct.rename(columns={"LHb": "pct_expr_LHb", "MHb": "pct_expr_MHb"})
    out = wide.merge(pct, on="mouse_gene", how="left")
    out["mouse_gene_upper"] = out["mouse_gene"].str.upper()
    out.to_csv(OUT_DIR / "mouse_lhb_mhb_average_expression_all_genes.csv", index=False)
    return out


def map_kim_to_mouse(kim: pd.DataFrame, var_names: list[str]) -> pd.DataFrame:
    link = pd.read_csv(ORTHOLOG_LINK)
    link["human_gene_upper"] = link["human_gene_upper"].astype(str).str.upper()
    link["mouse_gene"] = link["mouse_gene"].astype(str)
    link = link.drop_duplicates(["human_gene_upper", "mouse_gene"])

    var_lookup = pd.DataFrame({"mouse_gene": var_names})
    var_lookup["mouse_gene_upper"] = var_lookup["mouse_gene"].str.upper()
    var_lookup = var_lookup.drop_duplicates("mouse_gene_upper")

    mapped = kim.merge(link[["human_gene_upper", "mouse_gene"]], on="human_gene_upper", how="left")
    needs_fallback = mapped["mouse_gene"].isna()
    fallback = mapped.loc[needs_fallback].drop(columns=["mouse_gene"]).merge(
        var_lookup.rename(columns={"mouse_gene_upper": "human_gene_upper"}),
        on="human_gene_upper",
        how="left",
    )
    mapped = pd.concat([mapped.loc[~needs_fallback], fallback], ignore_index=True)
    mapped = mapped[mapped["mouse_gene"].notna()].copy()
    mapped = mapped.merge(var_lookup[["mouse_gene"]], on="mouse_gene", how="inner")
    mapped = mapped.drop_duplicates(["probe_set_id", "human_gene", "mouse_gene"])
    mapped["mouse_gene_upper"] = mapped["mouse_gene"].str.upper()
    mapped.to_csv(OUT_DIR / "kim_human_degs_mapped_to_mouse_genes.csv", index=False)
    return mapped


def categorize_kim_degs(mapped: pd.DataFrame, markers: pd.DataFrame, avg_expr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    marker_sets = {
        region: set(markers.loc[markers["region"].eq(region) & markers["is_region_marker"], "mouse_gene_upper"])
        for region in REGION_ORDER
    }
    marker_info = markers[markers["is_region_marker"]][
        ["mouse_gene", "mouse_gene_upper", "region", "logFC_vs_other_region", "pvals_adj", "pct_nz_group", "pct_nz_reference"]
    ].rename(columns={"region": "enriched_region", "pvals_adj": "region_marker_fdr"})

    records = []
    for row in mapped.itertuples(index=False):
        regions = [region for region in REGION_ORDER if row.mouse_gene_upper in marker_sets[region]]
        if len(regions) == 1:
            category = f"{regions[0]}-enriched"
            enriched_region = regions[0]
        elif len(regions) > 1:
            category = "both"
            enriched_region = "both"
        else:
            category = "not_region_enriched"
            enriched_region = ""
        rec = row._asdict()
        rec["kim_region_enrichment_category"] = category
        rec["enriched_region"] = enriched_region
        records.append(rec)
    cat = pd.DataFrame(records)
    cat = cat.merge(avg_expr, on=["mouse_gene", "mouse_gene_upper"], how="left")
    cat = cat.merge(marker_info, on=["mouse_gene", "mouse_gene_upper", "enriched_region"], how="left")
    cat = cat.sort_values(["kim_region_enrichment_category", "p_value", "human_gene"])
    cat.to_csv(OUT_DIR / "kim_degs_categorized_by_mouse_lhb_mhb_region_enrichment.csv", index=False)

    summary = (
        cat.drop_duplicates(["human_gene", "mouse_gene"])
        .groupby(["kim_region_enrichment_category", "direction"], dropna=False)
        .size()
        .reset_index(name="n_mapped_deg_genes")
    )
    total = cat.drop_duplicates(["human_gene", "mouse_gene"]).shape[0]
    restricted = cat[cat["kim_region_enrichment_category"].isin(["LHb-enriched", "MHb-enriched"])].drop_duplicates(["human_gene", "mouse_gene"]).shape[0]
    summary["percent_of_mapped_degs"] = summary["n_mapped_deg_genes"] / max(total, 1) * 100
    summary.attrs["n_total_mapped_gene_pairs"] = total
    summary.attrs["n_region_restricted_mapped_gene_pairs"] = restricted
    summary.to_csv(OUT_DIR / "kim_deg_region_enrichment_category_counts.csv", index=False)
    return cat, summary


def fisher_region_overlap(mapped: pd.DataFrame, markers: pd.DataFrame, var_names: list[str]) -> pd.DataFrame:
    mapped_genes = set(mapped["mouse_gene_upper"])
    universe = set(pd.Series(var_names).str.upper())
    rows = []
    for region in REGION_ORDER:
        marker_set = set(markers.loc[markers["region"].eq(region) & markers["is_region_marker"], "mouse_gene_upper"]) & universe
        overlap = sorted(marker_set & mapped_genes)
        a_count = len(overlap)
        b_count = len(marker_set - mapped_genes)
        c_count = len(mapped_genes - marker_set)
        d_count = max(len(universe) - a_count - b_count - c_count, 0)
        odds_ratio, p_value = stats.fisher_exact([[a_count, b_count], [c_count, d_count]], alternative="greater")
        rows.append(
            {
                "region": region,
                "n_region_markers": len(marker_set),
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
    out.to_csv(OUT_DIR / "kim_deg_x_mouse_lhb_mhb_marker_overlap_fisher.csv", index=False)
    return out


def plot_summary(cat: pd.DataFrame, summary: pd.DataFrame, overlap: pd.DataFrame) -> dict[str, Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")

    plot_counts = summary.copy()
    order = ["LHb-enriched", "MHb-enriched", "not_region_enriched", "both"]
    plot_counts["kim_region_enrichment_category"] = pd.Categorical(plot_counts["kim_region_enrichment_category"], categories=order, ordered=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1.25, 0.9]})
    sns.barplot(
        data=plot_counts,
        x="kim_region_enrichment_category",
        y="n_mapped_deg_genes",
        hue="direction",
        order=[x for x in order if x in set(plot_counts["kim_region_enrichment_category"].astype(str))],
        palette={"up_in_suicide": "#C44E52", "down_in_suicide": "#4C72B0"},
        ax=axes[0],
    )
    axes[0].set_xlabel("Mouse LHb/MHb region enrichment category")
    axes[0].set_ylabel("Kim mapped DEG gene count")
    axes[0].set_title("Kim suicide Hb DEGs categorized by broad mouse region enrichment")
    axes[0].tick_params(axis="x", labelrotation=25)

    heat = overlap.set_index("region")[["neg_log10_fdr"]].reindex(REGION_ORDER)
    annot = overlap.set_index("region")[["n_overlap"]].reindex(REGION_ORDER)
    sns.heatmap(
        heat,
        cmap="mako",
        vmin=0,
        linewidths=0.4,
        linecolor="white",
        annot=annot,
        fmt=".0f",
        cbar_kws={"label": "-log10 FDR"},
        ax=axes[1],
    )
    axes[1].set_xlabel("Kim DEG marker-set overlap")
    axes[1].set_ylabel("Mouse broad region marker set")
    axes[1].set_title("Fisher enrichment of Kim DEGs in region markers")
    fig.tight_layout()
    out = FIG_DIR / "kim_region_enrichment_summary.png"
    pdf = FIG_DIR / "kim_region_enrichment_summary.pdf"
    fig.savefig(out, dpi=240)
    fig.savefig(pdf)
    plt.close(fig)

    return {"summary_png": out, "summary_pdf": pdf}


def plot_dotplot(cat: pd.DataFrame) -> dict[str, Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    region_hits = cat[cat["kim_region_enrichment_category"].isin(["LHb-enriched", "MHb-enriched"])].copy()
    region_hits = region_hits.drop_duplicates(["mouse_gene", "human_gene"])
    if region_hits.empty:
        region_hits = cat.drop_duplicates(["mouse_gene", "human_gene"]).sort_values("p_value").head(40).copy()
    else:
        pieces = []
        for region in REGION_ORDER:
            sub = region_hits[region_hits["enriched_region"].eq(region)].sort_values("p_value").head(TOP_DOTPLOT_GENES_PER_REGION)
            pieces.append(sub)
        region_hits = pd.concat(pieces, ignore_index=True)

    region_hits["label"] = region_hits["human_gene"] + " / " + region_hits["mouse_gene"]
    region_hits["label"] = pd.Categorical(region_hits["label"], categories=region_hits["label"].tolist(), ordered=True)

    records = []
    for row in region_hits.itertuples(index=False):
        for region in REGION_ORDER:
            records.append(
                {
                    "label": row.label,
                    "human_gene": row.human_gene,
                    "mouse_gene": row.mouse_gene,
                    "enriched_region": row.enriched_region,
                    "direction": row.direction,
                    "region": region,
                    "mean_log1p_norm_expr": getattr(row, f"avg_expr_{region}"),
                    "pct_expr": getattr(row, f"pct_expr_{region}") * 100,
                }
            )
    dot = pd.DataFrame(records)
    mean_by_gene = dot.groupby("label", observed=True)["mean_log1p_norm_expr"].transform("mean")
    sd_by_gene = dot.groupby("label", observed=True)["mean_log1p_norm_expr"].transform("std").replace(0, np.nan)
    dot["expr_z_by_gene"] = ((dot["mean_log1p_norm_expr"] - mean_by_gene) / sd_by_gene).fillna(0).clip(-2.5, 2.5)
    dot.to_csv(OUT_DIR / "kim_region_enriched_deg_dotplot_values.csv", index=False)

    x_pos = {region: i for i, region in enumerate(REGION_ORDER)}
    y_labels = region_hits["label"].tolist()
    y_pos = {label: i for i, label in enumerate(y_labels)}
    fig, ax = plt.subplots(figsize=(5.8, max(5.2, 0.22 * len(y_labels) + 1.5)))
    scatter = ax.scatter(
        dot["region"].map(x_pos),
        dot["label"].map(y_pos),
        s=np.maximum(dot["pct_expr"], 1) * 3,
        c=dot["expr_z_by_gene"],
        cmap="vlag",
        vmin=-2.5,
        vmax=2.5,
        edgecolor="0.35",
        linewidth=0.2,
    )
    ax.set_xticks(range(len(REGION_ORDER)), REGION_ORDER)
    ax.set_yticks(range(len(y_labels)), y_labels)
    ax.invert_yaxis()
    ax.set_xlabel("Mouse broad habenula region")
    ax.set_ylabel("Kim DEG / mapped mouse gene")
    ax.set_title("Kim DEGs with broad LHb/MHb-enriched expression")
    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label("Mean expression z-score across regions")
    for size in [10, 30, 60]:
        ax.scatter([], [], s=size * 3, c="lightgray", edgecolor="0.35", label=f"{size}%")
    ax.legend(title="Pct expressed", bbox_to_anchor=(1.25, 1.0), loc="upper left", frameon=False)
    fig.tight_layout()
    png = FIG_DIR / "kim_region_enriched_deg_lhb_mhb_dotplot.png"
    pdf = FIG_DIR / "kim_region_enriched_deg_lhb_mhb_dotplot.pdf"
    fig.savefig(png, dpi=240)
    fig.savefig(pdf)
    plt.close(fig)
    return {"dotplot_png": png, "dotplot_pdf": pdf}


def write_methods_and_summary(
    kim: pd.DataFrame,
    mapped: pd.DataFrame,
    markers: pd.DataFrame,
    cat: pd.DataFrame,
    summary: pd.DataFrame,
    overlap: pd.DataFrame,
) -> Path:
    n_human = kim["human_gene"].nunique()
    n_mapped = mapped["mouse_gene"].nunique()
    n_lhb_markers = int(markers["region"].eq("LHb").mul(markers["is_region_marker"]).sum())
    n_mhb_markers = int(markers["region"].eq("MHb").mul(markers["is_region_marker"]).sum())
    restricted = cat[cat["kim_region_enrichment_category"].isin(["LHb-enriched", "MHb-enriched"])].drop_duplicates(["human_gene", "mouse_gene"])
    n_restricted = restricted.shape[0]
    lhb_hits = restricted[restricted["enriched_region"].eq("LHb")].shape[0]
    mhb_hits = restricted[restricted["enriched_region"].eq("MHb")].shape[0]

    lines = [
        "# Kim 2022 Region-Level Enrichment Replication",
        "",
        "This analysis replicates the Kim et al. cell-type enrichment concept using this project's broad mouse habenula regions rather than the original GSE137478 cell-type labels.",
        "",
        "## Inputs",
        "",
        "- Kim et al. Additional File 2, Table S2: human suicide-vs-control habenula DEGs.",
        "- Local mouse single-cell objects: `inputs/clean_objects/lateral_habenula_regions.h5ad` and `inputs/clean_objects/medial_habenula_regions.h5ad`.",
        "- Local mouse-human gene symbol link table from `inputs/figure_5_cross_species/human_mouse_gene_symbols.csv` with case-insensitive fallback to mouse symbols present in the h5ad objects.",
        "",
        "## Method",
        "",
        "- Human Kim DEG symbols were expanded where necessary, mapped to mouse genes, and restricted to genes present in both LHb and MHb h5ad objects.",
        "- LHb and MHb cells were combined at broad region level only.",
        "- Raw counts were normalized to 10,000 counts per cell and log1p transformed.",
        "- LHb and MHb region markers were computed using Wilcoxon rank-sum tests via Scanpy, analogous to the Seurat FindMarkers Wilcoxon step in the paper.",
        f"- Region markers were defined as `logFC > {MARKER_LOGFC_MIN}` and `adjusted p < {MARKER_FDR_MAX}`.",
        "- Kim DEGs were categorized by whether their mapped mouse gene was an LHb marker, an MHb marker, both, or not region-enriched.",
        "- Fisher's exact tests tested whether Kim DEGs were enriched among LHb or MHb region marker sets.",
        "",
        "## Key Counts",
        "",
        f"- Unique expanded Kim human DEG symbols: {n_human}",
        f"- Unique mapped mouse Kim DEG genes present in this LHb/MHb object pair: {n_mapped}",
        f"- LHb marker genes: {n_lhb_markers}",
        f"- MHb marker genes: {n_mhb_markers}",
        f"- Kim mapped DEG genes with broad LHb/MHb-restricted enrichment: {n_restricted}",
        f"- LHb-enriched Kim mapped DEG genes: {lhb_hits}",
        f"- MHb-enriched Kim mapped DEG genes: {mhb_hits}",
        "",
        "## Fisher Overlap",
        "",
        "| region | n_region_markers | n_overlap | odds_ratio | FDR | overlap_mouse_genes |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in overlap.itertuples(index=False):
        genes = row.overlap_mouse_genes if isinstance(row.overlap_mouse_genes, str) else ""
        if len(genes) > 160:
            genes = genes[:157] + "..."
        lines.append(
            f"| {row.region} | {row.n_region_markers} | {row.n_overlap} | {row.odds_ratio:.3g} | {row.p_adj:.3g} | {genes} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This version asks a narrower question than Kim et al.: whether the paper's human suicide Hb DEGs preferentially fall into genes enriched in broad mouse LHb or broad mouse MHb neurons/regions. It does not reproduce the original endothelial/mural/glial/neuron cell-type classification because the requested grouping is broad LHb vs MHb region level.",
        ]
    )
    return OUT_DIR


def main() -> dict[str, Path | int]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    kim = load_kim_degs()
    load_kim_original_table_s3()
    adata, norm, var_names = load_region_data()
    markers = compute_region_markers(adata)
    avg_expr = average_expression_table(adata, norm, var_names)
    mapped = map_kim_to_mouse(kim, var_names)
    cat, summary = categorize_kim_degs(mapped, markers, avg_expr)
    overlap = fisher_region_overlap(mapped, markers, var_names)
    fig_paths = plot_summary(cat, summary, overlap)
    dot_paths = plot_dotplot(cat)
    summary_path = write_methods_and_summary(kim, mapped, markers, cat, summary, overlap)

    print(f"Wrote {OUT_DIR}")
    print(f"Mapped Kim DEG mouse genes: {mapped['mouse_gene'].nunique()}")
    print(f"Filtered region markers: {int(markers['is_region_marker'].sum())}")
    print(overlap.to_string(index=False))
    return {
        "output_dir": OUT_DIR,
        "summary": summary_path,
        **fig_paths,
        **dot_paths,
    }


if __name__ == "__main__":
    main()
