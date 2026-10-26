# IBS7048 2026 Fall Midterm (Advanced Metagenomics)

Code and outputs for Problems 1-6. The dummy FASTQ files are not included; `scripts/00_make_dummy.sh` regenerates them.

## Layout
- `scripts/` : code for each problem
- `results/` : outputs and figures produced by the scripts

| Problem | Script | Output in `results/` |
|---|---|---|
| 1 dummy files + Bash sorting | `00_make_dummy.sh`, `01_sort_bash.sh` | `output_p1.txt` |
| 2 R sorting + comparison | `02_sort.R` | `output_p2.txt` |
| 3 budget calculation | `03_budget.R` | `output_p3.txt` |
| 4a DADA2 pipeline | `04a_dada2.R` | `quality_forward.png`, `quality_reverse.png`, `read_tracking.csv` |
| 4b truncLen = 120/120 | `04b_trunc120.R` | `output_p4b_merged_rows_truncLen120.txt` |
| 4c top-5 ASV NJ tree | `04c_tree.R` | `asv_tree.png`, `top5_ASV_NJ.nwk` |
| 5c Taxon_X fold change | `05c_taxonX.R` | `output_p5c.txt` |
| 6b OR / CASI rescaling | `06b_score.R` | `output_p6b.txt` |

Problems 5a, 5b, 5d, 5e and 6a, 6c, 6d are written answers only (no code needed).

## How to run

### Problems 1-3 (Bash, R)
Use Linux, macOS or WSL.
```bash
bash scripts/00_make_dummy.sh      # creates IBS7048_2026/raw (40 dummy FASTQ files)
cd IBS7048_2026
bash ../scripts/01_sort_bash.sh    # Problem 1 -> sorted_bash/
Rscript ../scripts/02_sort.R       # Problem 2 -> sorted_R/ (needs sorted_bash from Problem 1)
Rscript ../scripts/03_budget.R     # Problem 3 (independent of the files)
```

### Problem 4 (DADA2)
R 4.x with Bioconductor packages:
```r
install.packages(c("BiocManager", "ggplot2", "phangorn", "ape"))
BiocManager::install(c("dada2", "DECIPHER"))
```
Run from an empty working folder. Each script writes its files into the current folder.
```bash
Rscript scripts/04a_dada2.R     # -> quality plots, read_tracking.csv, seqtab_nochim.rds
Rscript scripts/04b_trunc120.R  # prints merged read counts for truncLen = c(120,120)
Rscript scripts/04c_tree.R      # needs seqtab_nochim.rds from 04a -> asv_tree.png, top5_ASV_NJ.nwk
```
The scripts use the example reads shipped with DADA2, so no data download is needed.

### Problems 5c and 6b
Base R only, no extra packages:
```bash
Rscript scripts/05c_taxonX.R
Rscript scripts/06b_score.R
```
