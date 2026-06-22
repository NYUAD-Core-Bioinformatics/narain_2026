# Reproducibility

Create the conda environment and run the full figure workflow from the repository root:

```bash
conda env create -f environment.yml
conda activate habenula-publication
./run_all.sh
```

The workflow starts from the cleaned AnnData objects and curated figure-specific input tables under `inputs/` and writes all regenerated outputs under `outputs/`.
