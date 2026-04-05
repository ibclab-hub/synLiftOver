# synLiftOver

**synLiftOver** is a synteny-aware coordinate conversion tool for transferring genomic annotations across genome assemblies. It improves annotation transfer robustness by combining synteny-guided alignment preprocessing with consistency-aware coordinate projection across single or multiple compatible alignment blocks.

## Overview

synLiftOver processes annotations through two steps:

1. **Synteny-guided alignment preprocessing (Step 1)**  
   Raw AXT alignment blocks are optionally filtered and/or merged using synteny block information.

2. **Coordinate projection (Step 2)**  
   Each query interval is projected onto the target assembly using one of two cases:
   - **Case 1**: single-block projection
   - **Case 2**: multi-block projection across compatible alignment blocks

This design enables robust annotation transfer in regions where conventional liftover approaches may produce fragmented or inconsistent results due to locally inconsistent alignments.

## Requirements

- Python 3
- Perl
- C++ compiler
- UCSC utilities: `chainNet`,`chainSplit`, `netSplit`, `netToAxt`, `faSize`, `faToTwoBit`
- Pairwise genome alignment resources:
  - AXT files (for projection)
  - chain/net files (for synteny generation and AXT conversion)

## Input formats

### Query annotation BED

Input annotations must be provided in tab-delimited BED format. The first three columns must be chromosome, start, and end. Additional columns are preserved in the output.

```text
chr1    10000    10200    gene1|exon1
```

### AXT alignment files

AXT files must be organized by chromosome:

```text
example/AXT/chr1.axt
```

### Target assembly size file

A chromosome size file for the target assembly is required. This can be generated using:

```bash
faSize --detailed target.fa > target.fa.size
```

### Synteny block BED file

Synteny blocks must be provided in BED-like format:

```text
1   chr1   246414   297604   +   chr7   0   51210   +
```

Fields: synteny ID, reference chromosome, reference start, reference end, reference strand, target chromosome, target start, target end, target strand.

> **Note**  
> In the manuscript workflow, synteny blocks were generated using MakeBlocks (embedded in RACA), followed by BED conversion using an in-house script. Precomputed example synteny files and wrapper scripts for MakeBlocks-based synteny generation are provided in the `example/` directory for reference.

## AXT generation

If chain/net files are available, AXT files can be generated using `netToAxt` on a per-chromosome basis before running synLiftOver.

## Step 1: Preprocess AXT alignments

The following commands assume execution from within the `example/` directory.

```bash
python3 ../src/step1/preprocess_axtFiltMerge.py \
    ./synteny/Conserved.Segments.bed \
    ./AXT/ \
    ./processed_AXT/ \
    1000 \
    0 \
    1 &> preprocess.log
```


#### Arguments:
| Position | Description |
|---|---|
| 1 | Synteny BED file |
| 2 | Input AXT directory |
| 3 | Output processed AXT directory |
| 4 | Minimum gap size for merging |
| 5 | Filter flag (0: disabled, 1: enabled) |
| 6 | Merge flag (0: disabled, 1: enabled) |

#### Preprocessing modes

| filter_flag | merge_flag | Mode |
|---|---|---|
| 0 | 0 | Raw |
| 1 | 0 | Filter |
| 0 | 1 | Merge (default) |
| 1 | 1 | Filter + Merge |

## Step 2: Run coordinate projection

```bash
perl ../src/step2/synLiftOver.pl \
    ./input.bed \
    ./fa/hs1.fa.size \
    ./outdir \
    ./AXT/ \
    0.95 \
    out \
    0 \
    0 &> log
```

#### Arguments:
| Position | Description |
|---|---|
| 1 | Input BED file |
| 2 | Target assembly size file |
| 3 | Output directory |
| 4 | AXT directory |
| 5 | Minimum-match threshold |
| 6 | Output prefix |
| 7 | Duplication flag (0: no duplicates, 1: allow duplicates) |
| 8 | Case 2 mode flag (0: strict, 1: relaxed) |

**Case 2 mode guidance:**
- Use **strict mode** (0) for same-species assembly pairs
- Use **relaxed mode** (1) for cross-species transfers, where greater alignment fragmentation is expected


## Synteny generation (Optional)

### Step 1. Prepare split chain/net files

MakeBlocks expects the ChainNet directory organized as:

```text
query_species/
└── target_species/
    ├── chain/
    │   └── chr*.chain
    └── net/
        └── chr*.net
```

### Step 2. Run MakeBlocks

```bash
../../scripts/run.Makeblocks Makeblocks_config >& log
```

This produces a `Conserved.Segments` file.

### Step 3. Convert to BED format

```bash
perl ../../scripts/csTobed.pl ./Conserved.Segments hg38 hs1 > ./Conserved.Segments.bed
```

## Example

Example input files used in this study are available on Figshare:
https://doi.org/10.6084/m9.figshare.31884754

You may download the files manually from the Figshare record and place them into the following directory structure:

```bash
example/
├── AXT/
│   └── chr1.axt
├── ChainNet/
│   └── hg38/hs1/
│       ├── chain/
│       │   └── chr1.chain
│       └── net/
│           └── chr1.net
└── processed_AXT/
    └── chr1.axt   # rename from processed_chr1.axt
```

Notes
- processed_chr1.axt should be renamed to chr1.axt after download.
- The provided files correspond to a representative chromosome-level example used for demonstration and testing.
- These files are sufficient to reproduce the example workflow described in this repository.

## Contact

For questions or issues, please open an issue on GitHub or contact the authors (qkrskdud0805@gmail.com).
