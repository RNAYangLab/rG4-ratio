# rG4 ratio
- rG4 ratio includes Python scripts for quantifying transcriptome-level SHALiEP-seq reverse transcript ribosomal stalling counts and characterizing signal distributions within RNA G-quadruplex (rG4) regions by G4 raio. See Lan et al, 2026.
- The workflow contains two main steps:
  1. Generate transcript-level stalling counts and coverage from aligned sequencing reads as Bam file with corresponding species transcript fasta file.
  2. Calculate rG4 region siganal metrics, including average rG4 region counts, GINI index and G4 ratio.
## Requirement
- Python 3 is recommend.
- Required Python packages: pip install pysam biopython
## Step 1: Generate stalling counts and coverage
The input BAM file should have the associated index files.  
For each transcript in input transcript FASTA file, this script records:
- the transcript sequence.
- the number of ribosomal stalling counts at each nucleotide position.
- the sequencing coverage at each nucleotide position.
```bash
python get_stalling_counts.py --threads threads --fasta transcript.fa --bam sample.bam --output sample.counts
```
  
## Step 2：Calculate rG4 region siganal metrics
For each annotated rG4 region, it extracts the corresponding sequence and stalling counts, then calculates average rG4 region counts, GINI index and G4 ratio:  
### average rG4 region counts：
The average stalling count across all guanine (G) positions within the annotated rG4 region.
### GINI index (`gini`)
Gini index was calculated from the SHALiPE-seq libraries with reads continuous number of G residues in G-tract as described.

$$
\frac{\sum_{i=1}^{n}\sum_{j=1}^{n}\left|r_i-r_j\right|}
{2n^2\bar{r}}
$$

where $n$ denotes the number of G residues in the G-tracts (continuous runs of guanine in the G-rich region), $r_i$ denotes the reads number in SHALiPE profiling at position $i$.
### G4 ratio
G4 ratio was calculated from the SHALiPE-seq libraries as the ratio of reads counts on the 3’-guanine (max G) to the summed read counts across the corresponding G-cluster (max G + former G).

```math
\mathrm{G4\ Ratio} =
\begin{cases}
\displaystyle
\frac{1}{n}\sum_{i=1}^{n}
\frac{G_{\mathrm{max},i}}
{G_{\mathrm{former},i}+G_{\mathrm{max},i}},
& n \ge 4 \\[6pt]
0,
& n < 4
\end{cases}
```
where $n$ denotes the number of G-tracts, and and $G_{\mathrm{max},i}$ and $G_{\mathrm{former},i}$ denote the read count at 3’-guanine and preceding guanine. The G4 ratio was set to 0 when fewer than four structural G-tracts were identified.
```bash
python Caculate_GINI_G4ratio.py -RG4list predicted_RG4_list.tsv -C sample.counts -O sample_rG4.tsv -t threads
```
