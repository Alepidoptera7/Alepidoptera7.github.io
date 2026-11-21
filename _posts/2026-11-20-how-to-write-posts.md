---
title: How to write posts on this site
description: A quick reference for code blocks, equations, tables, figures and callouts. Delete this post once you're comfortable.
tags: [meta]
math: true
---

Every post is a Markdown file in the `_posts/` folder, named
`YYYY-MM-DD-short-title.md`. The date in the name sets the publish date, and
the rest becomes the URL. The three posts on this site are examples, so
replace them with your own.

* TOC
{:toc}

## Front matter

Each post starts with a block like this:

```yaml
---
title: Aligning long reads with minimap2
description: One-sentence summary shown in post lists and link previews.
tags: [long-reads, alignment]
math: true        # only if the post has equations
---
```

## Code

Fence code with three backticks and name the language for highlighting.
Readers get a Copy button on every block.

```bash
samtools view -b -q 30 -F 0x904 sample.bam > sample.filtered.bam
samtools index sample.filtered.bam
```

```r
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
```

Inline code works too: `bcftools view -i 'QUAL>30'`.

## Equations

Set `math: true` in the front matter, then write LaTeX. Inline:
$p = 10^{-Q/10}$. Display:

$$
\mathrm{TPM}_i = \frac{c_i / \ell_i}{\sum_j c_j / \ell_j} \times 10^6
$$

## Tables

| Tool    | Task               | Input   |
|:--------|:-------------------|:--------|
| STAR    | Spliced alignment  | FASTQ   |
| salmon  | Quantification     | FASTQ   |
| DESeq2  | Differential expr. | Counts  |

## Figures

Put images in `assets/img/` and reference them with a caption:

```markdown
{% raw %}![Volcano plot of DE genes]({{ '/assets/img/volcano.png' | relative_url }}){% endraw %}
*Figure 1. 412 genes pass padj < 0.05.*
```

## Callouts

Add `{: .note}` or `{: .warn}` on the line after a paragraph:

Reference genome versions matter: never mix GRCh37 and GRCh38 coordinates in one analysis.
{: .warn}

Saving your environment with `conda env export > environment.yml` makes a post reproducible.
{: .note}

## Drafts

Put unfinished posts in a `_drafts/` folder (no date in the file name). They
show up locally with `bundle exec jekyll serve --drafts` and are never
published.
