---
title: "kallisto vs STAR: pseudoalignment or full alignment for RNA-seq?"
description: How the two most common RNA-seq quantification tools differ, how each handles spliced reads, how their output feeds featureCounts and tximport, and complete commands to run both on the same paired-end samples.
tags: [rna-seq, alignment, quantification, tutorial]
---

## Two different questions

**STAR** (Dobin et al., 2013) is a *spliced aligner*. For every read it
finds where in the genome the read came from, base by base, including reads
that span exon–exon junctions. The output is a BAM file with exact
coordinates, mismatches and splice sites for each read. Counting reads per
gene is a separate step that happens afterwards.

**kallisto** (Bray et al., 2016) does *pseudoalignment*. It never works out
where in a transcript a read sits. It breaks the read into k-mers (31 bases
by default), looks them up in an index built from the transcriptome, and
records only *which transcripts the read is compatible with*. An
expectation–maximisation (EM) algorithm then shares out reads that fit
several transcripts, giving an estimated count and TPM for every transcript.

So the core trade-off is:

- STAR tells you **where** each read is. That's slow and memory-hungry, but
  you get rich information.
- kallisto tells you **how much** of each transcript is there. It's fast and
  light, but that's all it tells you.

## Splice awareness

RNA-seq reads come from mature mRNA, where the introns have been cut out and
the exons joined. A 100 bp read that starts near the end of exon 2 can finish
in exon 3, but in the genome those two pieces may sit thousands of bases
apart. Any tool that quantifies RNA-seq has to deal with these
*junction-spanning* reads, which make up a large share of reads in a typical
library. STAR and kallisto deal with them in completely different ways.

### STAR: splice-aware alignment

A DNA aligner such as BWA or Bowtie2 expects each read to match one
continuous stretch of the genome. Given a junction read, it either clips off
the part in the next exon or fails to place the read at all. STAR is built to
split reads across introns:

1. **Seed search.** STAR finds the longest prefix of the read that matches
   the genome exactly, its *maximal mappable prefix*, then repeats the
   search from where that match stopped. A junction read naturally breaks
   into two seeds, one per exon.
2. **Stitching.** The seeds are joined into one alignment, allowing a large
   gap between them: the intron. Gaps that start and end with the canonical
   GT…AG splice motif score better than non-canonical ones.
3. **Known junctions.** The junctions in the GTF are built into the index at
   `genomeGenerate`, which is why it needs `--sjdbGTFfile` and
   `--sjdbOverhang`. Reads that only overhang an exon by a few bases can still
   be placed across an annotated junction.
4. **New junctions.** STAR also finds junctions that aren't in the
   annotation. With `--twopassMode Basic` it collects the junctions found in a
   first pass, adds them to the index, and aligns again, which improves
   alignment of reads across those new junctions.

In the BAM, the intron shows up as an `N` in the CIGAR string. `51M2370N50M`
means 51 bases match, then 2,370 bases of intron are skipped, then 50 more
bases match. To see some spliced reads and count them:

```bash
BAM=star/SAMPLE1/Aligned.sortedByCoord.out.bam

# First few spliced alignments: read name, chromosome, position, CIGAR
samtools view $BAM | awk '$6 ~ /N/ {print $1, $3, $4, $6}' | head

# How many alignments are spliced?
samtools view -c $BAM
samtools view $BAM | awk '$6 ~ /N/' | wc -l
```

STAR also writes every junction it found to `SJ.out.tab`, one per line. The
columns are chromosome, intron start, intron end, strand (0 undefined,
1 plus, 2 minus), intron motif, whether it's annotated (1) or new (0),
uniquely mapping reads, multi-mapping reads, and maximum overhang. Listing the
best-supported unannotated junctions:

```bash
awk '$6 == 0 && $7 >= 10' star/SAMPLE1/SJ.out.tab | sort -k7,7nr | head
```

### kallisto: splicing is already in the index

kallisto never sees the genome. Its index is built from *transcript*
sequences, which are already spliced, so exon 2 and exon 3 sit side by side in
the reference. A junction read is just a continuous run of k-mers along the
transcript, with no gap to bridge. In effect, kallisto gets splicing for free.

The catch is that it only works for splicing **someone has already
annotated**:

- A read across an unannotated junction contains k-mers that don't exist in
  the index. It may fail to pseudoalign, or be assigned only from the k-mers
  that do match.
- Reads from introns, such as pre-mRNA, retained introns or genomic DNA
  contamination, have no transcript to match. Usually they simply fail to
  pseudoalign. Occasionally they share enough k-mers with a transcript to be
  wrongly assigned to it. Recent kallisto versions can reduce this by masking
  k-mers that also occur elsewhere in the genome (the `--d-list` option of
  `kallisto index`).
- kallisto can't tell you anything about the junctions themselves. It reports
  how much of each *annotated* isoform it thinks is present, never which
  junctions the reads crossed.

So if splicing is your question, for example differential splicing, intron
retention or new isoforms, you need a splice-aware aligner like STAR. If
splicing is only something that must be handled correctly on the way to
expression estimates, kallisto's approach is enough for a well-annotated
transcriptome.

## Side by side

| | STAR | kallisto |
|:--|:--|:--|
| Method | Spliced alignment to the genome | Pseudoalignment to the transcriptome |
| Reference | Genome FASTA + GTF annotation | Transcript FASTA |
| Main output | BAM file (+ optional gene counts) | Transcript counts and TPM |
| Natural level | Gene | Transcript (sum to genes afterwards) |
| Memory for human | About 32 GB to build and load the index | A few GB; runs on a laptop |
| Speed | Fast for an aligner, minutes per sample on a server | Very fast, often a few minutes on a laptop |
| Finds unannotated junctions, novel isoforms | Yes | No; only quantifies transcripts in the index |
| Variants, fusions, IGV visualisation | Yes, from the BAM | No |
| Uncertainty estimates | No | Yes, via bootstraps (`-b`) |
| Multi-mapping reads | Excluded from gene counts by default | Shared out between transcripts by the EM |

## When to use which

Pick **kallisto** when:

- you want differential expression on a well-annotated organism (human,
  mouse, yeast…), especially at the transcript level;
- you have many samples and limited compute;
- you'd like bootstrap uncertainty, e.g. for sleuth.

Pick **STAR** when:

- you need a BAM file: to view reads in IGV, call variants, check coverage,
  or run QC tools such as RSeQC or Qualimap;
- you care about splicing (rMATS, LeafCutter), fusions (STAR-Fusion) or
  unannotated junctions;
- your organism's annotation is incomplete, so a transcriptome-only index
  would miss reads.

For a standard gene-level differential expression analysis on a
well-annotated genome, the two approaches generally lead to very similar
conclusions. The differences tend to concentrate in genes with many similar
paralogs, in short genes, and in lowly expressed genes. Many pipelines get the
best of both: nf-core/rnaseq, for example, aligns with STAR by default and
then quantifies with salmon, a pseudoaligner similar to kallisto.

## Setup

Both tools are on Bioconda. The commands below were written for STAR 2.7.11b
and kallisto 0.51.1. The `subread` package provides featureCounts, used
later.

```bash
mamba create -n rnaseq -c conda-forge -c bioconda star=2.7.11b kallisto=0.51.1 subread samtools multiqc
mamba activate rnaseq
```

The examples assume paired-end reads that have already been trimmed, named
`trimmed/SAMPLE_R1.fastq.gz` and `trimmed/SAMPLE_R2.fastq.gz`, and a file
`samples.txt` with one sample name per line.

### Reference files

Use one annotation release for everything, so gene and transcript IDs match
across tools. Here that's GENCODE release 50 on GRCh38. Check
[gencodegenes.org](https://www.gencodegenes.org/human/) for newer releases.

```bash
mkdir -p ref && cd ref
BASE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50

# For STAR: genome + matching annotation (primary assembly versions)
wget $BASE/GRCh38.primary_assembly.genome.fa.gz
wget $BASE/gencode.v50.primary_assembly.annotation.gtf.gz

# For kallisto: transcript sequences
wget $BASE/gencode.v50.transcripts.fa.gz

# STAR can't read gzipped references
gunzip GRCh38.primary_assembly.genome.fa.gz gencode.v50.primary_assembly.annotation.gtf.gz
cd ..
```

## Running STAR

### 1. Build the index (once per genome and annotation)

`--sjdbOverhang` should be read length minus 1, so 100 for 101 bp reads.
This step needs roughly 32 GB of RAM for the human genome and takes about
an hour.

```bash
STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir ref/star_index \
     --genomeFastaFiles ref/GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile ref/gencode.v50.primary_assembly.annotation.gtf \
     --sjdbOverhang 100
```

### 2. Align each sample

```bash
mkdir -p star
while read s; do
  STAR --runThreadN 8 \
       --genomeDir ref/star_index \
       --readFilesIn trimmed/${s}_R1.fastq.gz trimmed/${s}_R2.fastq.gz \
       --readFilesCommand zcat \
       --outSAMtype BAM SortedByCoordinate \
       --quantMode GeneCounts \
       --outFileNamePrefix star/${s}/
  samtools index star/${s}/Aligned.sortedByCoord.out.bam
done < samples.txt
```

On macOS, use `--readFilesCommand gzcat` instead of `zcat`.

Each sample folder now contains:

- `Aligned.sortedByCoord.out.bam`: the alignments, ready for IGV;
- `ReadsPerGene.out.tab`: read counts per gene;
- `Log.final.out`: mapping statistics. Check "Uniquely mapped reads %" first.

### 3. Pick the right count column

`ReadsPerGene.out.tab` has four columns: gene ID, then counts for an
unstranded library, a forward-stranded library and a reverse-stranded
library. The first four rows are summary lines (`N_unmapped` and so on), not
genes.

Most modern kits, such as Illumina TruSeq Stranded and other dUTP-based
protocols, are **reverse-stranded, so use column 4**. If you're not sure,
compare the column totals: in a stranded library one of columns 3 and 4 will
be far larger than the other.

```bash
# Sum columns 3 and 4, skipping the 4 summary rows
tail -n +5 star/SAMPLE1/ReadsPerGene.out.tab | awk '{f+=$3; r+=$4} END {print "forward:", f, "reverse:", r}'
```

## Running kallisto

### 1. Build the index (once per transcriptome)

This takes a few minutes and modest memory.

```bash
kallisto index -i ref/kallisto_gencode50.idx -t 8 ref/gencode.v50.transcripts.fa.gz
```

### 2. Quantify each sample

```bash
mkdir -p kallisto
while read s; do
  kallisto quant \
       -i ref/kallisto_gencode50.idx \
       -o kallisto/${s} \
       -b 100 \
       -t 8 \
       --rf-stranded \
       trimmed/${s}_R1.fastq.gz trimmed/${s}_R2.fastq.gz \
       2> kallisto/${s}.log
done < samples.txt
```

The `2>` at the end saves kallisto's run summary, which it prints to the
terminal, so MultiQC can read it later.

- `--rf-stranded` is for reverse-stranded (dUTP) libraries, the same ones that
  use column 4 in STAR. Use `--fr-stranded` for forward-stranded kits, or drop
  the flag for unstranded libraries. If you choose wrong, the share of
  pseudoaligned reads in `run_info.json` collapses. That makes it a quick check
  worth doing on one sample first.
- `-b 100` runs 100 bootstrap rounds. They're needed for sleuth, but you can
  drop them if you only use DESeq2 or edgeR.
- For single-end reads, add `--single -l 200 -s 20`, with your library's
  estimated fragment length mean and SD. kallisto can't estimate these from
  single-end data.

Each sample folder contains `abundance.tsv` (estimated counts and TPM per
transcript), `abundance.h5` (the same plus bootstraps) and `run_info.json`
(including the pseudoalignment rate).

Newer kallisto documentation also describes a `kallisto bus` workflow,
wrapped by the `kb` tool, for bulk data. `kallisto quant` still works and is
the form most downstream tools expect.

## From reads to a count matrix

DESeq2, edgeR and limma all start from a table of **counts**: genes in rows,
samples in columns. TPM or FPKM values won't work as input, because the
statistics need to know how many fragments each number is based on. STAR and
kallisto reach that table by different routes, and the routes don't cross:

| Route | What the tool gives you | Counting / import step | Into DESeq2 with |
|:--|:--|:--|:--|
| STAR + featureCounts | Genome BAM | featureCounts | `DESeqDataSetFromMatrix` |
| STAR alone | `ReadsPerGene.out.tab` | Built in (`--quantMode GeneCounts`) | `DESeqDataSetFromMatrix` |
| STAR + salmon | Transcriptome BAM | salmon, then tximport | `DESeqDataSetFromTximport` |
| kallisto | `abundance.tsv` per sample | tximport | `DESeqDataSetFromTximport` |

featureCounts needs reads placed on genome coordinates, which only STAR
produces; current kallisto writes no BAM at all. tximport reads *estimated
abundances*, which STAR's plain counts aren't. So the pairings are fixed:
**STAR goes with featureCounts, kallisto goes with tximport.**

### STAR → featureCounts

featureCounts (from the Subread package) takes BAM files and a GTF, and
assigns each fragment to the gene whose exons it overlaps. Counting is
splice-aware too. For a read with an `N` in its CIGAR, only the aligned blocks
are compared with the annotation, so the skipped intron is never counted as
overlap. A read spanning exon 2 and exon 3 of the same gene adds one count to
that gene.

```bash
mkdir -p counts
featureCounts -T 8 \
    -p --countReadPairs \
    -s 2 \
    -t exon -g gene_id \
    -a ref/gencode.v50.primary_assembly.annotation.gtf \
    -o counts/featureCounts.txt \
    star/*/Aligned.sortedByCoord.out.bam
```

- `-p --countReadPairs` counts **fragments** rather than individual reads.
  In recent Subread versions `-p` only declares the data paired-end; without
  `--countReadPairs` each mate is counted separately and your counts roughly
  double.
- `-s 2` means reverse-stranded: the same library type as STAR's column 4
  and kallisto's `--rf-stranded`. Use `0` for unstranded and `1` for
  forward-stranded.
- `-t exon -g gene_id` counts overlaps with exons and adds them up per
  `gene_id`.
- By default, multi-mapping reads (`-M` to include them) and reads
  overlapping more than one gene (`-O` to allow) are **not** counted.

Always look at `counts/featureCounts.txt.summary`. It shows how many fragments
were `Assigned` and why the rest weren't. A large `Unassigned_NoFeatures`
count usually means intronic or intergenic reads, or a wrong `-s` setting.

featureCounts can also count reads per junction. Add `-J` and it writes
`counts/featureCounts.txt.jcounts`, a junction-level table that's a handy
starting point for splicing questions.

**featureCounts or STAR's `GeneCounts`?** For unique reads with default
settings they agree closely: the STAR manual states its counts match those of
htseq-count run with default settings. Use featureCounts when you need more
control (multi-mappers, exon-level counts with `-f`, junction counts, other
feature types) or when your BAMs come from a different aligner.

### kallisto → tximport

kallisto's `est_counts` are not whole numbers. The EM shares ambiguous reads
between transcripts, so a transcript can end up with 41.7 reads. kallisto
also reports each transcript's *effective length*, which tximport uses.
tximport does three things:

1. reads every sample's `abundance.tsv`;
2. adds transcript counts up to gene level, using the transcript-to-gene table;
3. calculates, for every gene and sample, the **average transcript length**,
   weighted by how much of each isoform is present.

The third step is the important one. Imagine a gene that makes the same number
of mRNA molecules in both conditions, but switches from a short isoform to a
long one. Longer molecules produce more fragments, so a plain count would
call the gene "upregulated". `DESeqDataSetFromTximport` turns the average
lengths into a per-sample correction (an *offset*), so DESeq2 corrects for
this automatically. It also rounds the estimated counts to integers for you.

tximport's `countsFromAbundance` argument controls what goes in the count
matrix:

- `"no"` (the default): kallisto's estimated counts plus the length matrix.
  Use this with `DESeqDataSetFromTximport`, which applies the offset itself.
- `"lengthScaledTPM"`: counts rebuilt from TPM so that the length correction
  is built in. Use this when the next tool won't take an offset, e.g.
  limma-voom or a plain edgeR count matrix.
- `"scaledTPM"` and `"dtuScaledTPM"` are for transcript-level work such as
  differential transcript usage (with `txOut = TRUE`).

GENCODE transcript names are long `|`-separated headers.
`ignoreAfterBar = TRUE` trims them to the first field (`ENST…`), so they match
the transcript-to-gene table. Build that table once from the same FASTA you
indexed:

```bash
zcat ref/gencode.v50.transcripts.fa.gz | grep '^>' | sed 's/^>//' \
  | awk -F'|' 'BEGIN {OFS="\t"} {print $1, $2}' > ref/tx2gene.tsv
```

Because everything comes from the same GENCODE release, the gene IDs
(`ENSG…` with version numbers) are identical in the featureCounts table and
in the tximport output. You can put the two matrices side by side without any
ID conversion.

### Getting both: STAR + salmon

If you want STAR's BAM *and* tximport-style estimates with a length offset,
ask STAR to also write alignments in transcript coordinates. Add
`TranscriptomeSAM` to the align command:

```bash
--quantMode TranscriptomeSAM GeneCounts
```

Each sample then gets an `Aligned.toTranscriptome.out.bam` as well. salmon,
in alignment mode, quantifies it like kallisto would quantify reads, and
tximport reads the result with `type = "salmon"`. This is how nf-core/rnaseq
works by default. The transcript FASTA has to match the GTF used for the STAR
index exactly, so build it from the GTF instead of downloading it separately:

```bash
mamba install -c conda-forge -c bioconda salmon gffread

gffread -w ref/transcripts_from_gtf.fa \
        -g ref/GRCh38.primary_assembly.genome.fa \
        ref/gencode.v50.primary_assembly.annotation.gtf

salmon quant -t ref/transcripts_from_gtf.fa -l A \
             -a star/SAMPLE1/Aligned.toTranscriptome.out.bam \
             -p 8 -o salmon/SAMPLE1
```

### Into DESeq2

`samples.csv` holds a `sample` column and a `condition` column.

```r
library(DESeq2)
library(tximport)

samples <- read.csv("samples.csv")
samples$condition <- factor(samples$condition)

## Route 1: kallisto -> tximport (counts + length offset)
files <- file.path("kallisto", samples$sample, "abundance.tsv")
names(files) <- samples$sample
tx2gene <- read.delim("ref/tx2gene.tsv", header = FALSE)

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene,
                ignoreAfterBar = TRUE)
dds_kallisto <- DESeqDataSetFromTximport(txi, colData = samples,
                                         design = ~ condition)

## Route 2: STAR -> featureCounts
fc <- read.delim("counts/featureCounts.txt", comment.char = "#",
                 check.names = FALSE)
counts_fc <- as.matrix(fc[, 7:ncol(fc)])        # first 6 columns are annotation
rownames(counts_fc) <- fc$Geneid
colnames(counts_fc) <- basename(dirname(colnames(counts_fc)))  # star/S1/... -> S1
counts_fc <- counts_fc[, samples$sample]        # same order as samples.csv
dds_star <- DESeqDataSetFromMatrix(counts_fc, colData = samples,
                                   design = ~ condition)

## Alternative to route 2: STAR's own GeneCounts (column 4 = reverse-stranded)
read_star <- function(f, col = 4) {
  x <- read.delim(f, header = FALSE, skip = 4)   # skip the 4 summary rows
  setNames(x[[col]], x[[1]])
}
counts_gc <- sapply(file.path("star", samples$sample, "ReadsPerGene.out.tab"),
                    read_star)
colnames(counts_gc) <- samples$sample

## Same analysis for both
res_kallisto <- results(DESeq(dds_kallisto), alpha = 0.05)
res_star     <- results(DESeq(dds_star),     alpha = 0.05)

## How well do the two routes agree?
common <- intersect(rownames(res_kallisto), rownames(res_star))
cor(res_kallisto[common, "log2FoldChange"], res_star[common, "log2FoldChange"],
    use = "complete.obs")
table(kallisto = res_kallisto[common, "padj"] < 0.05,
      star     = res_star[common, "padj"] < 0.05, useNA = "ifany")
```

Most significant genes should agree. The ones that don't are worth a look:
they're often genes with close paralogs, where kallisto's EM shares reads that
featureCounts throws away as multi-mapping or ambiguous.

### Mistakes to avoid

- **TPM or FPKM as DESeq2 input.** Always give it counts.
- **Rounding kallisto's counts and using `DESeqDataSetFromMatrix`.** It runs,
  but you lose the length offset. Use `DESeqDataSetFromTximport`.
- **Mismatched strandedness.** STAR column 4, featureCounts `-s 2` and
  kallisto `--rf-stranded` all describe the same library. If one step is
  wrong, that route's counts quietly drop.
- **Paired-end featureCounts without `--countReadPairs`.** You count reads, not
  fragments.
- **Mixing annotation releases.** Use one GENCODE release for the genome
  index, the GTF, the transcript FASTA and the transcript-to-gene table.

## Quick QC for both

MultiQC reads both STAR's `Log.final.out` and kallisto's logs. Run it over the
whole project to see mapping rates for every sample in one report:

```bash
multiqc star/ kallisto/ -o qc/multiqc
```

## Summary

- **kallisto**: fast, light and transcript-level. A great default for
  expression and differential expression analysis on well-annotated
  organisms.
- **STAR**: slower and heavier, but you get real alignments, which you need
  for anything to do with sequence, structure or visualisation.
- **Splicing:** STAR finds junctions in the genome; kallisto relies on
  junctions already present in annotated transcripts.
- **Counting:** STAR → featureCounts → `DESeqDataSetFromMatrix`;
  kallisto → tximport → `DESeqDataSetFromTximport`.
- If you're unsure, run both. The kallisto step adds only minutes.

## References

- Dobin A. *et al.* (2013). STAR: ultrafast universal RNA-seq aligner.
  *Bioinformatics* 29(1):15–21.
  [doi:10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)
- Bray N.L. *et al.* (2016). Near-optimal probabilistic RNA-seq quantification.
  *Nature Biotechnology* 34:525–527.
  [doi:10.1038/nbt.3519](https://doi.org/10.1038/nbt.3519)
- Liao Y., Smyth G.K., Shi W. (2014). featureCounts: an efficient general
  purpose program for assigning sequence reads to genomic features.
  *Bioinformatics* 30(7):923–930.
  [doi:10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656)
- Soneson C., Love M.I., Robinson M.D. (2015). Differential analyses for
  RNA-seq: transcript-level estimates improve gene-level inferences.
  *F1000Research* 4:1521.
  [doi:10.12688/f1000research.7563.1](https://doi.org/10.12688/f1000research.7563.1)
- [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
  and [kallisto documentation](https://pachterlab.github.io/kallisto/manual)
