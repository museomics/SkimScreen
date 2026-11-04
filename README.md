# SkimScreen
Tools for screening genome skims for non-endogenous/contaminant reads.


## Kraken2 classification and plotting
**Stages the database once per job.** The `.k2d` files are copied to node-local SSD (`/tmp`). The job checks for enough free space first and always deletes the copy on exit, including after failure, `scancel` or a time-limit kill.
**Classifies with `--memory-mapping`.** The database is read from local disk instead of being loaded into RAM, so `--mem` doesn't need to exceed the size of `hash.k2d`. Every sample reuses the same staged copy. This is slightly slower, but more memory efficient for running Kraken2 agains large DBs.

Two scripts to classify paired-end reads from many samples against a large Kraken2 database on a SLURM cluster, then summarise the results across samples.

| Script | What it does |
|---|---|
| `kraken2_run.sh` | SLURM job: runs `k2 classify` on every sample in a manifest, with the database staged to node-local disk |
| `plot_kraken_reports.py` | Merges the per-sample Kraken2 reports into tables and plots |

---

### `kraken2_run.sh`
#### Usage
```bash
sbatch kraken2_run.sh <db_dir> <out_dir> <manifest.csv>
```

- `db_dir`: extracted Kraken2 database containing `hash.k2d`, `opts.k2d` and `taxo.k2d`
- `out_dir`: output root; one subdirectory is created per sample
- `manifest.csv`: comma-separated sample sheet (see below)

#### Manifest
Columns are matched by **header name**, so column order doesn't matter and extra columns are ignored. Required columns:

| Column | Accepted names |
|---|---|
| Sample ID | `ID` |
| Forward reads | `forward` or `R1` |
| Reverse reads | `reverse` or `R2` |

```csv
ID,forward,reverse
UK001-A01,/path/to/UK001-A01_R1.fastq.gz,/path/to/UK001-A01_R2.fastq.gz
```
Fields must not contain commas.


#### Output
```
<out_dir>/<ID>/
├── kraken_output.txt        # per-read classifications
├── kraken_report.txt        # per-taxon summary report
├── classified_{1,2}.fastq
└── unclassified_{1,2}.fastq
```

#### Notes
- Defaults: `--partition=himem`, `--mem=250G`, `--cpus-per-task=8`. Edit the `#SBATCH` lines to suit your cluster.
- Staging is a one-off cost per job. For example, the 455 GB `k2_NCBI_reference_20251007` database took about 45 minutes to copy from BeeGFS. Run all samples in one job (i.e. one manifest.csv) rather than one job per sample.
- The conda path and environment name (`read_qc`) are set in the *Environment* section of the script.
- Avoid memory-mapping the database directly from a network filesystem (BeeGFS/Lustre/NFS). Kraken2's random lookups perform very poorly there.

---

### `plot_kraken_reports.py`
#### Usage
```bash
python plot_kraken_reports.py --indir <out_dir> --outdir kraken_plots --rank genus --top 20
```

Reads every `<indir>/<sample>/kraken_report.txt` and uses the subdirectory name as the sample name. Both the standard 6-column report and the 8-column `--report-minimizer-data` report are supported.

| Option | Default | Description |
|---|---|---|
| `--indir` | *required* | Directory of per-sample subdirectories |
| `--outdir` | *required* | Output directory |
| `--rank` | `G` | Rank for the stacked bar and heatmap: `D P C O F G S`, or a name such as `order` |
| `--top` | `20` | Taxa coloured in the stacked bar; the rest are grouped as *Other* |
| `--heatmap-top` | `30` | Taxa shown in the heatmap |
| `--sort` | `name` | Sample order: `name` (natural sort) or `dominant` (grouped by each sample's main taxon) |
| `--no-cluster` | off | Disable hierarchical clustering of the heatmap |
| `--report-name` | `kraken_report.txt` | Report filename to look for |

#### Output
| File | Contents |
|---|---|
| `stacked_bar_<rank>.png/.pdf` | Top taxa per sample, plus *Other*, *Classified, no `<rank>` assigned* and *Unclassified* |
| `heatmap_<rank>.png/.pdf` | Samples × top taxa, log-scale colours, clustered |
| `classification_overview.png/.pdf` | Unclassified vs Eukaryota / Bacteria / Archaea / Viruses / other classified |
| `kraken_reports_long.tsv` | Every report line from every sample, in long format |
| `kraken_<rank>_reads.tsv` / `_pct.tsv` | Samples × taxa matrix (read counts / %) |
| `kraken_overview.tsv` | Data behind the overview plot |

#### Interpreting the plots
- All percentages are of **total reads**, including unclassified, so each bar sums to 100%.
- **Other**: reads assigned at the chosen rank to taxa outside the top `--top`. Top taxa are chosen by mean % across samples, so a taxon that dominates only a few samples can end up in *Other*. Check the heatmap or the `_pct.tsv` table.
- **Classified, no `<rank>` assigned**: reads classified only to a higher rank, or to a lineage with no taxon at that rank.
- Taxa that share a name are labelled with their taxid, e.g. `Bacillus (1386)` the bacterium and `Bacillus (7010)` the stick insect.

---

### Requirements
- Kraken2, including the `k2` wrapper command (`k2 classify`)
- SLURM, with node-local disk larger than the database
- Python ≥ 3.9 with `pandas`, `numpy` and `matplotlib`; `scipy` is optional (heatmap clustering)

```bash
conda install -c bioconda -c conda-forge kraken2 pandas numpy matplotlib scipy
```
