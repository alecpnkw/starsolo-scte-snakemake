# starsolo-scte-snakemake
Alec Pankow
2022-10-04

A [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for running scTE (single-cell transposable element) 
quantitation from scRNAseq data with a single command. The scTE pipeline is as described in [He et al. 2021 (Nature Communications)](http://dx.doi.org/10.1038/s41467-021-21808-x).

### Quick start

Clone and install snakemake (with conda environment)

```bash
git clone https://github.com/alecpnkw/starsolo-scte-snakemake.git
cd starsolo-scte-snakemake
conda env create --file environment.yaml
```

Modify the configuration file (`config/config.yaml`) to suit your run:

```yaml
genome: 
  name: "hg38"
  fasta: "resources/genome.fa"
  gtf: "resources/annot.gtf"
fastqs: "<path-to-fastqs-comma-separated>"
umi_whitelist: "<path-to-umi-whitelist>"
```

Preview and run snakemake (see [documentation](https://snakemake.readthedocs.io/en/stable/) for full list of options)

```bash
# preview
snakemake --dry-run

# run on cluster using --profile with conda envs
snakemake \
  --jobs <n> \
  --use-conda \
  --profile <cluster-profile> \
  --keep-going \
  --conda-prefix <path-to-conda-envs-dir>
```

See [this](https://github.com/Snakemake-Profiles/doc) page for further documentation on Snakemake profiles. 

### Acknowledgements

Based on previous work by Roosheel Patel (@roosheelpatel)

### References

He, Jiangping, Isaac A. Babarinde, Li Sun, Shuyang Xu, Ruhai Chen, Junjie Shi, Yuanjie Wei, et al. 2021. “Identifying Transposable Element Expression Dynamics and Heterogeneity during Development at the Single-Cell Level with a Processing Pipeline scTE.” Nature Communications 12 (1): 1456. https://doi.org/10.1038/s41467-021-21808-x.
