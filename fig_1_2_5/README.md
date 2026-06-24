# Habenula publication analysis package

This repository contains the cleaned inputs, analysis scripts, and figure outputs for the mouse habenula stress-phenotype analyses and human habenula comparison analyses. The package starts from cleaned AnnData objects and curated figure-specific input tables.

## Layout

- `inputs/clean_objects/`: cleaned mouse AnnData objects used by the figure scripts.
- `inputs/figure_1_mouse_habenula_atlas/`: curated Figure 1 MHb panel C and LHb panel E region labels and UMAP coordinates.
- `inputs/figure_2_mouse_differential_expression/`: mouse differential-expression tables used for phenotype and region DE summaries.
- `inputs/figure_5_cross_species/`: Yalcinbas and Kim human reference inputs plus the human-mouse gene-symbol map.
- `scripts/`: one script per figure-generation step.
- `outputs/figure_1_mouse_habenula_atlas/`: mouse habenula atlas figure and source tables.
- `outputs/figure_2_mouse_differential_expression/`: phenotype and region DE summary figures and tables.
- `outputs/figure_5_cross_species/`: Yalcinbas and Kim comparison figures and source tables.
- `docs/methods.md`: manuscript-style methods for generating the figures from the cleaned inputs.

## Environment

```bash
conda env create -f environment.yml
conda activate habenula-publication
```

The `.h5ad` inputs are configured for Git LFS in `.gitattributes`.

## Recreate outputs

Run all figure-generation scripts from the repository root:

```bash
./run_all.sh
```

Run individual steps with:

```bash
python scripts/make_figure_1_mouse_habenula_atlas.py
python scripts/make_figure_2_mouse_differential_expression.py
python scripts/make_figure_5_yalcinbas_lhb_mapping.py
python scripts/make_figure_5_yalcinbas_mhb_mapping.py
python scripts/make_figure_5_yalcinbas_sensitivity.py
python scripts/make_figure_5_yalcinbas_summary.py
python scripts/make_figure_5_kim_region_enrichment.py
python scripts/make_figure_5_kim_celltype_enrichment.py
```
