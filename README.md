# starsolo-scte-snakemake
Alec Pankow
2025-09-11

A [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for running TE
quantitation from scRNAseq data with a single command. Previously this used the scTE pipeline is as described in 
[He et al. 2021 (Nature Communications)](http://dx.doi.org/10.1038/s41467-021-21808-x), but now runs [IRescue](https://github.com/bodegalab/irescue).

### Quick start

Clone and install snakemake (with conda environment)

```bash
git clone https://github.com/alecpnkw/starsolo-scte-snakemake.git
cd starsolo-scte-snakemake
conda env create --file environment.yaml
```

Modify the configuration file (`config/config.yaml`) to suit your run:

```yaml
samples: "config/samples.csv"

# paths to genome to use for mapping
genome: 
  name: "hg38"
  fasta: "resources/genome.fa"
  gtf: "resources/annot.gtf"

# starsolo cell barcode / UMI configuration
soloCBstart: 1
soloCBlen: 16
soloUMIstart: 17
soloUMIlen: 12
soloBarcodeReadLength: 29
umi_whitelist: "<path-to-umi-whitelist>"
```

Preview and run snakemake (see [documentation](https://snakemake.readthedocs.io/en/stable/) for full list of options)

```bash
# preview
snakemake --dry-run

# currently configured to be run on lsf with the snakemake lsf executor by defailt
snakemake \
  --jobs <n> \
  --use-conda \
  --keep-going
```

See [this page](https://github.com/Snakemake-Profiles/doc) for further documentation on Snakemake profiles. See [this page](https://github.com/BEFH/snakemake-executor-plugin-lsf) for information on the recommended the snakemake lsf executor plugin. 

### Acknowledgements

Based on previous work by Roosheel Patel (@roosheelpatel)

### References

He, Jiangping, Isaac A. Babarinde, Li Sun, Shuyang Xu, Ruhai Chen, Junjie Shi, Yuanjie Wei, et al. 2021. “Identifying Transposable Element Expression Dynamics and Heterogeneity during Development at the Single-Cell Level with a Processing Pipeline scTE.” Nature Communications 12 (1): 1456. https://doi.org/10.1038/s41467-021-21808-x.

Polimeni, Benedetto, Federica Marasca, Valeria Ranzani, and Beatrice Bodega. 2024. “IRescue: Uncertainty-Aware Quantification of Transposable Elements Expression at Single Cell Level.” Nucleic Acids Research 52 (19): e93. https://doi.org/10.1093/nar/gkae793.


