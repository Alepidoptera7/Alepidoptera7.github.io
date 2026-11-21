---
title: "Nextflow vs Snakemake: choosing a workflow manager for bioinformatics"
date: 2026-11-20 21:00:00 -0800
description: How Nextflow and Snakemake think about pipelines, where each one shines, and the same small RNA-seq pipeline (fastp, kallisto, MultiQC) written in both.
tags: [workflows, nextflow, snakemake, reproducibility, tutorial]
---

## Why use a workflow manager at all?

A shell script with a `for` loop works until it doesn't. Halfway through 48
samples a job runs out of memory, and you have to work out which outputs are
complete, which are half-written, and how to restart without redoing
everything. Then a reviewer asks you to rerun with a newer annotation, and you
no longer remember which tool versions you used.

Workflow managers handle that for you:

- **Restarting:** finished steps are skipped and only missing or outdated work
  is rerun.
- **Parallelism:** independent samples run at the same time, on a laptop, a
  cluster or the cloud, without changing the pipeline.
- **Software:** each step declares its own conda environment or container, so
  the versions are written down and reproducible.
- **Provenance:** logs, resource usage and the exact commands run are kept for
  every task.

**Nextflow** and **Snakemake** are the two most widely used in
bioinformatics. Both do all of the above, but they're built on different ideas,
and that changes how you write and think about a pipeline.

## Two different mental models

### Snakemake: start from the files you want

Snakemake works like GNU Make. You write **rules** that say how to make an
output file from input files, using wildcards like `{sample}` in file names.
Then you ask for the final files, and Snakemake works **backwards**: to make
`results/kallisto/A/abundance.tsv` it needs trimmed reads for sample A, which
need the raw FASTQ for sample A, and so on. It builds that dependency graph
(the DAG) before running anything.

File names are the interface between steps. If you can describe each step as
"this file pattern becomes that file pattern", Snakemake will feel natural.
It's written in Python, and any Python code can go in the Snakefile.

### Nextflow: data flows through channels

Nextflow works **forwards**. Each step is a **process**: a script with
declared inputs and outputs. Processes are connected by **channels**, queues
that carry values or files from one process to the next. When a sample pair
arrives on a channel, Nextflow starts a task for it. When the task finishes,
its outputs flow on to the next process.

File names matter much less. Every task runs in its own isolated folder under
`work/`, with its inputs linked in, so two tasks can both write `out.bam`
without clashing. You choose which results to copy into your results folder
at the end. The language is a Groovy-based DSL. Since Nextflow 26.04 the
stricter, better-checked "strict syntax" is the default.

## Side by side

| | Snakemake | Nextflow |
|:--|:--|:--|
| Model | Rules, pulled backwards from target files (like Make) | Processes, connected by channels, pushed forwards (dataflow) |
| Language | Python plus rule syntax | Groovy-based DSL |
| Where tasks run | In your project folder; outputs go where you name them | Each task in its own `work/` subfolder; results published at the end |
| Restarting | Reruns whatever's missing or outdated (changed inputs, code, parameters, software environment) | `-resume` reuses cached tasks whose inputs and script haven't changed |
| Software | Conda per rule; containers through Apptainer/Singularity | Conda, Docker, Apptainer/Singularity, Podman per process |
| Clusters and cloud | Executor plugins: SLURM, LSF, Kubernetes, Google Batch and others | Built-in executors: SLURM, SGE, PBS, LSF, AWS Batch, Google Batch, Azure Batch, Kubernetes |
| Ready-made pipelines | Snakemake Workflow Catalog, snakemake-wrappers | nf-core: a large curated set of community pipelines with shared standards |
| Easiest when | You know Python and think in files | You want to run nf-core pipelines, or scale to clusters or the cloud |

## The same pipeline in both

To make the comparison concrete, here's one small RNA-seq pipeline written
twice:

1. **fastp** trims each paired-end sample;
2. **kallisto** builds an index once, then quantifies each trimmed sample;
3. **MultiQC** collects every fastp and kallisto log into one report.

This covers the patterns that matter: one step per sample (fan-out), a shared
reference, and one step that waits for all samples (fan-in). Both versions
expect the same layout:

```
data/SAMPLE_R1.fastq.gz, data/SAMPLE_R2.fastq.gz   # raw reads
ref/gencode.v50.transcripts.fa.gz                  # transcriptome
```

Install both tools from Bioconda. The code below was written for Nextflow
26.10 and Snakemake 9.

```bash
mamba create -n workflows -c conda-forge -c bioconda nextflow snakemake
mamba activate workflows
```

### Snakemake version

Project layout:

```
Snakefile
config.yaml
envs/fastp.yaml
envs/kallisto.yaml
envs/multiqc.yaml
```

`config.yaml`:

```yaml
transcriptome: ref/gencode.v50.transcripts.fa.gz
```

`envs/kallisto.yaml` (the other two are the same, with `fastp` or `multiqc`):

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - kallisto
```

`Snakefile`:

```python
configfile: "config.yaml"

# Find samples from the raw file names: data/A_R1.fastq.gz -> "A"
SAMPLES = glob_wildcards("data/{sample}_R1.fastq.gz").sample


# The final files we want. Snakemake works backwards from these.
rule all:
    input:
        expand("results/kallisto/{sample}/abundance.tsv", sample=SAMPLES),
        "results/multiqc/multiqc_report.html",


rule fastp:
    input:
        r1="data/{sample}_R1.fastq.gz",
        r2="data/{sample}_R2.fastq.gz",
    output:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        json="results/qc/{sample}.fastp.json",
        html="results/qc/{sample}.fastp.html",
    threads: 4
    conda:
        "envs/fastp.yaml"
    shell:
        "fastp --in1 {input.r1} --in2 {input.r2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "--detect_adapter_for_pe --thread {threads} "
        "--json {output.json} --html {output.html}"


rule kallisto_index:
    input:
        config["transcriptome"],
    output:
        "results/ref/kallisto.idx",
    threads: 8
    resources:
        mem_mb=16000,
    conda:
        "envs/kallisto.yaml"
    shell:
        "kallisto index -t {threads} -i {output} {input}"


rule kallisto_quant:
    input:
        idx="results/ref/kallisto.idx",
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
    output:
        "results/kallisto/{sample}/abundance.tsv",
    log:
        "results/logs/kallisto/{sample}.log",
    threads: 4
    conda:
        "envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.idx} -o results/kallisto/{wildcards.sample} "
        "-t {threads} --rf-stranded {input.r1} {input.r2} 2> {log}"


# Fan-in: needs the logs from every sample
rule multiqc:
    input:
        expand("results/qc/{sample}.fastp.json", sample=SAMPLES),
        expand("results/logs/kallisto/{sample}.log", sample=SAMPLES),
    output:
        "results/multiqc/multiqc_report.html",
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc --force -o results/multiqc {input}"
```

Notice that there's no explicit wiring between steps. `kallisto_quant` asks
for `results/trimmed/{sample}_R1.fastq.gz`, and Snakemake works out that
`fastp` makes that file. The `expand()` calls in `multiqc` create the fan-in:
that rule can't start until every sample's logs exist.

Running it:

```bash
snakemake -n                                   # dry run: show what would run
snakemake --cores 8 --sdm conda                # run locally, conda envs per rule
snakemake --dag | dot -Tsvg > dag.svg          # draw the dependency graph
```

`--sdm` is short for `--software-deployment-method`. Run the same command
again after it finishes and Snakemake reports nothing to do. Delete one
sample's trimmed reads, or edit the fastp command, and only the affected jobs
rerun.

### Nextflow version

Project layout:

```
main.nf
nextflow.config
```

`nextflow.config`:

```groovy
params {
    reads         = 'data/*_R{1,2}.fastq.gz'
    transcriptome = 'ref/gencode.v50.transcripts.fa.gz'
}

profiles {
    conda {
        conda.enabled = true
    }
    docker {
        docker.enabled = true
    }
}
```

`main.nf`:

```groovy
process FASTP {
    tag "${sample}"
    cpus 4
    conda 'bioconda::fastp'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_R{1,2}.trimmed.fastq.gz"), emit: reads
    path "${sample}.fastp.json", emit: json

    script:
    """
    fastp --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 ${sample}_R1.trimmed.fastq.gz --out2 ${sample}_R2.trimmed.fastq.gz \\
        --detect_adapter_for_pe --thread ${task.cpus} \\
        --json ${sample}.fastp.json --html ${sample}.fastp.html
    """
}

process KALLISTO_INDEX {
    cpus 8
    memory '16 GB'
    conda 'bioconda::kallisto'

    input:
    path transcriptome

    output:
    path 'kallisto.idx'

    script:
    """
    kallisto index -t ${task.cpus} -i kallisto.idx ${transcriptome}
    """
}

process KALLISTO_QUANT {
    tag "${sample}"
    cpus 4
    conda 'bioconda::kallisto'

    input:
    tuple val(sample), path(reads)
    path index

    output:
    path "${sample}", emit: quant
    path "${sample}.kallisto.log", emit: log

    script:
    """
    kallisto quant -i ${index} -o ${sample} -t ${task.cpus} --rf-stranded \\
        ${reads[0]} ${reads[1]} 2> ${sample}.kallisto.log
    """
}

process MULTIQC {
    conda 'bioconda::multiqc'

    input:
    path logs

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    main:
    // Emits one item per sample: [ "A", [A_R1.fastq.gz, A_R2.fastq.gz] ]
    read_pairs = channel.fromFilePairs(params.reads, checkIfExists: true)

    index = KALLISTO_INDEX(file(params.transcriptome))
    FASTP(read_pairs)
    KALLISTO_QUANT(FASTP.out.reads, index)

    // Fan-in: gather every log into one list, then run MultiQC once
    all_logs = FASTP.out.json.mix(KALLISTO_QUANT.out.log).collect()
    MULTIQC(all_logs)

    publish:
    quant   = KALLISTO_QUANT.out.quant
    multiqc = MULTIQC.out
}

output {
    quant {
        path 'kallisto'
    }
    multiqc {
        path 'multiqc'
    }
}
```

Here the wiring is explicit. `FASTP.out.reads` is passed straight into
`KALLISTO_QUANT`, and file names only matter inside each task. A few details
worth knowing:

- `channel.fromFilePairs` groups `A_R1`/`A_R2` into one item with the sample
  name attached. That item is the "sample" that flows through the pipeline.
- `KALLISTO_INDEX` is called with a single file, so its output is a *value*
  that every `KALLISTO_QUANT` task can reuse, rather than being used up by the
  first sample.
- `.mix()` merges two channels, and `.collect()` waits for everything and
  delivers it as one list. That's how fan-in works in Nextflow.
- The `publish:` section and `output` block copy final results into
  `results/kallisto` and `results/multiqc`. This is the newer replacement for
  the `publishDir` directive you'll see in older pipelines.

Running it:

```bash
nextflow run main.nf -profile conda            # run locally with conda
nextflow run main.nf -profile conda -resume    # rerun, reusing finished tasks
nextflow run main.nf -profile conda -with-report report.html -with-dag dag.html
```

`-resume` checks each task's inputs and script against its cache in `work/`
and only reruns what changed. `-with-report` writes an HTML report of CPU,
memory and time per task, which is very useful for setting resources before
moving to a cluster.

## Moving to a cluster

Both tools run the same pipeline on an HPC cluster without changing the
pipeline code; only the configuration changes.

**Snakemake** uses an executor plugin. For SLURM, install
`snakemake-executor-plugin-slurm` and put the settings in a profile,
`profiles/slurm/config.yaml`:

```yaml
executor: slurm
jobs: 50
software-deployment-method: conda
default-resources:
  slurm_partition: "standard"
  mem_mb: 4000
  runtime: 60
```

```bash
snakemake --profile profiles/slurm
```

Each rule becomes a SLURM job, using its `threads` and `resources`.

**Nextflow** has SLURM support built in. Add a profile to
`nextflow.config`:

```groovy
profiles {
    slurm {
        process.executor = 'slurm'
        process.queue    = 'standard'
    }
}
```

```bash
nextflow run main.nf -profile conda,slurm -resume
```

Each task is submitted as a SLURM job with the `cpus` and `memory` set in its
process. On a cluster, run Nextflow itself inside a small, long-running job or
a `tmux` session, because the main process has to stay alive to submit and
track tasks.

## Don't write it if you don't have to

For standard analyses, someone has probably written and tested the pipeline
already. The best-known collection is **nf-core**, with peer-reviewed
Nextflow pipelines for RNA-seq, variant calling, ATAC-seq, single-cell and
more. They share conventions, test datasets and documentation.
`nf-core/rnaseq` alone covers trimming, STAR alignment, salmon quantification
and extensive QC:

```bash
nextflow run nf-core/rnaseq -r <release> -profile docker \
    --input samplesheet.csv --outdir results \
    --fasta GRCh38.primary_assembly.genome.fa.gz \
    --gtf gencode.v50.primary_assembly.annotation.gtf.gz \
    --gencode
```

Pick `<release>` from [nf-co.re/rnaseq](https://nf-co.re/rnaseq) and always
pin it, so a rerun next year uses the same pipeline. `--gencode` tells the
pipeline to expect GENCODE-style annotation.

On the Snakemake side, the
[Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/)
lists community workflows, and the
[snakemake-wrappers](https://snakemake-wrappers.readthedocs.io/) project
provides tested, reusable rules for hundreds of common tools.

## Which should you choose?

Pick **Snakemake** if:

- you're comfortable in Python and think about analyses as files turning
  into other files;
- you're building a custom analysis, or one that will keep changing, where
  being quick to write matters more than scale;
- you mostly run on a workstation or an HPC cluster.

Pick **Nextflow** if:

- you want to run nf-core pipelines, or contribute to them;
- you need to run on the cloud (AWS, Google, Azure) or Kubernetes, or at very
  large scale;
- you want every task isolated in its own folder, with containers as the
  default.

Honestly, either is far better than a pile of shell scripts. The concepts
carry over: once you understand rules and wildcards, or processes and
channels, learning the other tool mostly means learning new syntax.

## References

- Köster J., Rahmann S. (2012). Snakemake: a scalable bioinformatics workflow
  engine. *Bioinformatics* 28(19):2520–2522.
  [doi:10.1093/bioinformatics/bts480](https://doi.org/10.1093/bioinformatics/bts480)
- Mölder F. *et al.* (2021). Sustainable data analysis with Snakemake.
  *F1000Research* 10:33.
  [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2)
- Di Tommaso P. *et al.* (2017). Nextflow enables reproducible computational
  workflows. *Nature Biotechnology* 35:316–319.
  [doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)
- Ewels P.A. *et al.* (2020). The nf-core framework for community-curated
  bioinformatics pipelines. *Nature Biotechnology* 38:276–278.
  [doi:10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)
- [Snakemake documentation](https://snakemake.readthedocs.io/) and
  [Nextflow documentation](https://www.nextflow.io/docs/latest/)
