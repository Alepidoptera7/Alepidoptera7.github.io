---
title: Reading FASTA files and computing N50 in plain Python
description: A twenty-line FASTA reader, GC content that ignores Ns, and what N50 actually measures, with no dependencies.
tags: [python, genome-assembly, tutorial]
---

Biopython is great, but sometimes you're on a cluster node without it, or you
just want to understand what's happening. Here's everything you need to get
basic statistics from an assembly using only the standard library.

## A streaming FASTA reader

FASTA records can wrap sequence over many lines, so the reader collects lines
until it sees the next `>` header. Using a generator means even a 3 GB genome
never sits in memory all at once.

```python
import gzip

def read_fasta(path):
    """Yield (name, sequence) for each record. Handles .gz files."""
    opener = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)
```

## GC content

Ambiguous bases (`N`) shouldn't count toward the denominator, or a scaffold
full of gaps will look AT-rich.

```python
def gc_content(seq):
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    acgt = gc + seq.count("A") + seq.count("T")
    return gc / acgt if acgt else 0.0
```

## N50

N50 is the length of the shortest contig in the smallest set of longest contigs
that together cover half the assembly. In code: sort lengths from longest to
shortest and walk down until the running total reaches half.

```python
def n50(lengths):
    lengths = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    running = 0
    for length in lengths:
        running += length
        if running >= half:
            return length
    return 0
```

A worked example with five contigs of 100, 80, 60, 40 and 20 bp:

| Contig | Length | Running total | Reached 150 bp? |
|:-------|-------:|--------------:|:----------------|
| 1      | 100    | 100           | no              |
| 2      | 80     | 180           | yes, so N50 = 80 |
| 3      | 60     | 240           |                 |
| 4      | 40     | 280           |                 |
| 5      | 20     | 300           |                 |

N50 rewards contiguity but says nothing about correctness. An assembler that
wrongly joins two chromosomes gets a *better* N50. Pair it with a
completeness check such as BUSCO, or with alignment to a reference.
{: .warn}

## Putting it together

```python
import sys

lengths, gcs = [], []
for name, seq in read_fasta(sys.argv[1]):
    lengths.append(len(seq))
    gcs.append(gc_content(seq))

total = sum(lengths)
print(f"contigs:  {len(lengths):,}")
print(f"total bp: {total:,}")
print(f"longest:  {max(lengths):,}")
print(f"N50:      {n50(lengths):,}")
print(f"mean GC:  {sum(g * l for g, l in zip(gcs, lengths)) / total:.1%}")
```

```bash
python asm_stats.py assembly.fasta.gz
```

## The same thing with existing tools

If Biopython is available, `SeqIO.parse(path, "fasta")` replaces the reader and
`Bio.SeqUtils.gc_fraction` replaces the GC function. On the command line,
`seqkit stats -a assembly.fasta` reports N50 alongside other statistics.
