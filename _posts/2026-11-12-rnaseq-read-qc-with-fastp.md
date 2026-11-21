---
title: Quality control of RNA-seq reads with FastQC and fastp
description: What Phred scores mean, how to trim adapters and low-quality bases with fastp, and which FastQC warnings you can safely ignore for RNA-seq.
tags: [rna-seq, qc, tutorial]
math: true
---

Before aligning a single read, it pays to check what came off the sequencer.
This post walks through a minimal QC pass for paired-end RNA-seq: look at the
raw reads, trim them, and summarise everything in one report.

## What a quality score means

Each base in a FASTQ file carries a Phred quality score $Q$, which encodes the
probability $P$ that the base call is wrong:

$$
Q = -10 \log_{10} P \qquad\Longleftrightarrow\qquad P = 10^{-Q/10}
$$

| Phred score | Error probability | Base call accuracy |
|------------:|------------------:|-------------------:|
| 10          | 1 in 10           | 90%                |
| 20          | 1 in 100          | 99%                |
| 30          | 1 in 1,000        | 99.9%              |
| 40          | 1 in 10,000       | 99.99%             |

Scores are stored as single ASCII characters with an offset of 33
("Phred+33"), so `I` (ASCII 73) means $Q = 40$. Decoding a line in Python:

```python
qual_line = "IIIIHHGG#"
scores = [ord(ch) - 33 for ch in qual_line]
print(scores)  # [40, 40, 40, 40, 39, 39, 38, 38, 2]
```

## Step 1: look at the raw reads

```bash
mkdir -p qc/raw
fastqc --threads 4 --outdir qc/raw reads/*.fastq.gz
```

## Step 2: trim with fastp

fastp detects and removes adapters, trims low-quality ends, and drops reads
that end up too short, all in one pass.

```bash
mkdir -p trimmed qc/fastp

for r1 in reads/*_R1.fastq.gz; do
  sample=$(basename "$r1" _R1.fastq.gz)
  fastp \
    --in1 "$r1" --in2 "reads/${sample}_R2.fastq.gz" \
    --out1 "trimmed/${sample}_R1.fastq.gz" --out2 "trimmed/${sample}_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --cut_right --cut_window_size 4 --cut_mean_quality 20 \
    --length_required 36 \
    --thread 4 \
    --html "qc/fastp/${sample}.html" --json "qc/fastp/${sample}.json"
done
```

`--cut_right` slides a 4-base window from the 5′ end and cuts the read as soon
as the window's mean quality falls below 20. Reads shorter than 36 bases after
trimming are discarded, together with their mate, so the pairs stay in sync.

## Step 3: one report for everything

```bash
multiqc qc/ --outdir qc/multiqc
```

MultiQC picks up both the FastQC output and the fastp JSON files and puts every
sample side by side, which makes outliers easy to spot.

## FastQC warnings that are normal for RNA-seq

FastQC was designed with genomic DNA in mind, so a few of its red flags are
expected for RNA-seq:

- **Per base sequence content** often fails over the first 10 to 12 bases.
  Random hexamer priming isn't truly random, and this bias is a known artefact.
- **Sequence duplication levels** are high because highly expressed genes
  produce many identical fragments. That is biology, not PCR trouble.
- **Overrepresented sequences** may be mitochondrial or ribosomal transcripts.
  Worth noting, but not a reason to throw the sample away.

What *is* worth worrying about: a sharp quality drop across all positions,
adapter content above a few percent after trimming, or one sample that looks
nothing like its replicates.

Write down the versions you used (`fastp --version`, `fastqc --version`).
Default settings change between releases.
{: .note}
