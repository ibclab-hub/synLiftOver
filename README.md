# synLiftOver

synLiftOver transfers genomic intervals across genome assemblies and species using synteny-guided alignment preprocessing and single-block or multi-block coordinate projection.

## Requirements

- Python 3, Perl, and a C++11 compiler
- Pairwise alignments in AXT format and target chromosome sizes
- Paired synteny intervals for optional alignment filtering and merging

The core runtime uses Python's standard library. Preparing alignment and synteny inputs may additionally require UCSC utilities and MakeBlocks, as described below.

## Build

From the repository root:

```bash
bash compile.sh
```

This creates `src/step2/synLiftOver.out`. Recompile on your machine before running.

## Repository layout

- `src/step1/`: synteny-guided AXT preprocessing
- `src/step2/`: coordinate projection
- `scripts/`: existing MakeBlocks executable and Conserved.Segments conversion script
- `example/`: example commands, configuration, and input layout

## Inputs

### Query BED

Use headerless, tab-delimited BED with 0-based, half-open, positive-length intervals:

```text
chr1    10000    10200    feature1    0    +
```

Additional metadata is preserved. A recognized `+` or `-` in BED column 6 is converted to target strand. For exon annotations, provide each exon as a separate interval. Coordinate-bearing extra fields, including BED12 block coordinates, are passed through rather than recalculated.

### Chromosome AXT files

Place uncompressed AXT alignments in `AXT_DIR/chr1.axt`, with filenames matching source BED chromosome names. Keep each source chromosome in its own file and preserve unique numeric block IDs and block order. Use chromosome filenames without extra periods before `.axt` for preprocessing.

AXT coordinates are interpreted as 1-based inclusive. Negative-strand target coordinates use the oriented target coordinate system. Paired alignment strings must have equal lengths.

### Target chromosome sizes

Provide a two-column file of target chromosome names and lengths. For example:

```bash
faSize --detailed target.fa > target.fa.size
```

### Paired synteny intervals

Use the nine-column format shown in `example/synteny/Conserved.Segments.bed`:

```text
1   chr1   246414   297604   +   chr7   0   51210   +
```

Columns are synteny ID, source chromosome/start/end/strand, and target chromosome/start/end/strand. Coordinates are forward 0-based, half-open. IDs must be unique. The relative alignment direction is positive when the two strand fields agree and negative otherwise.

## Step 1: Preprocess AXT alignments

Run from `example/`:

```bash
python3 ../src/step1/preprocess_axtFiltMerge.py \
    ./synteny/Conserved.Segments.bed \
    ./AXT/ ./processed_AXT/ \
    1000 0 1 ./fa/hs1.fa.size \
    > preprocess.out 2> preprocess.log
```

Use a new output directory for each preprocessing run.

| Position | Argument |
|---|---|
| 1 | Paired synteny file |
| 2 | Input AXT directory |
| 3 | Output AXT directory |
| 4 | Exclusive upper limit for the target gap (`GAP`) |
| 5 | Filter flag |
| 6 | Merge flag |
| 7 | Target chromosome sizes |

| Filter | Merge | Mode |
|---:|---:|---|
| 0 | 0 | Retain blocks without filtering or merging |
| 1 | 0 | Filter |
| 0 | 1 | Merge |
| 1 | 1 | Filter+Merge |

For the original-alignment baseline, use `AXT/` directly in Step 2.

**Filter** retains blocks contained within the same paired synteny interval on both genomes, with matching target chromosome and relative direction.

**Merge** joins consecutive original block IDs assigned to the same synteny interval. Source blocks must be directly adjacent, target chromosome and direction must match, and the oriented target gap must satisfy `0 <= gap < GAP`. With `GAP=1000`, the maximum permitted target gap is 999 bp. Gaps are represented by target `N` and source `-` characters.

**Filter+Merge** combines these operations. Filtering restricts the alignment context and may reduce both incorrect mappings and recovery. Merging can recover additional intervals or simplify an existing multi-block projection without changing its output.

## Step 2: Project coordinates

From `example/`, using the processed alignments:

```bash
perl ../src/step2/synLiftOver.pl \
    ./input.bed ./fa/hs1.fa.size ./outdir ./processed_AXT/ \
    0.95 out 0 0 > lift.out 2> lift.log
```

| Position | Argument |
|---|---|
| 1 | Source BED |
| 2 | Target chromosome sizes |
| 3 | Output directory |
| 4 | AXT directory |
| 5 | Minimum-match threshold |
| 6 | Output filename |
| 7 | Allow multiple candidate outputs: `0` (default) or `1` |
| 8 | Multi-block mode: `0` strict (default), `1` relaxed |

The minimum-match parameter applies to source overlap for partial single-block candidates, using a strict `>` comparison. It is not an overall sequence-identity or target-coverage threshold.

Strict and relaxed modes require compatible endpoint blocks on the same target chromosome and direction. Strict mode also excludes pairs with another block overlapping the intervening source or compatible target gap. Relaxed mode omits that exclusion. Strict mode is a starting point for same-species assembly transfers; relaxed mode can improve recovery in fragmented cross-species alignments. Choose according to the input alignments and desired recovery/consistency tradeoff.

Multi-block candidates are attempted first. With multiple outputs disabled, an accepted multi-block result suppresses single-block output. Otherwise, fully contained single-block candidates precede partial candidates. Allowing multiple outputs reports candidate projections, not independently established biological orthologs.

The mapped BED is written to `outdir/out`. Per-source-chromosome inputs, outputs, mapping logs, and warnings are kept under `outdir/tmp/`. Unmapped intervals are absent from the mapped BED. The wrapper runs one process per source chromosome and combines outputs in source-chromosome order. Use separate output directories for separate runs.

## Preparing alignments and synteny

### AXT from chain/net

Existing example commands are in `example/cmds_all`. With split chain/net files and source/target two-bit genomes, run `netToAxt` for each source chromosome. UCSC utilities such as `chainNet`, `chainSplit`, `netSplit`, `faSize`, and `faToTwoBit` can be used when preparing these inputs.

### MakeBlocks synteny workflow

The manuscript workflow uses MakeBlocks from RACA. The existing `scripts/run.Makeblocks` and `scripts/csTobed.pl` support the following example.

Organize split chain/net inputs as:

```text
example/ChainNet/hg38/hs1/chain/chr1.chain
example/ChainNet/hg38/hs1/net/chr1.net
```

From `example/synteny/`, edit `Makeblocks_config` for the input directories, species names, and resolution. The example uses `150000` bp.

```bash
../../scripts/run.Makeblocks Makeblocks_config > log 2>&1
perl ../../scripts/csTobed.pl Conserved.Segments hg38 hs1 > Conserved.Segments.bed
```

Use the resulting paired synteny file in Step 1. Precomputed example files are also available. Other synteny sources can be used after converting them to the paired forward-coordinate format above.

## Example data

Download the representative chromosome-level data from Figshare:
https://doi.org/10.6084/m9.figshare.31884754

The AXT and chain/net files in this repository are empty placeholders. Populate `example/AXT/` and `example/ChainNet/` with the downloaded data. To use the provided processed alignment directly, rename `processed_chr1.axt` to `chr1.axt` and place it under `example/processed_AXT/`. To generate processed alignments yourself, choose a fresh output directory for Step 1.

## Contact

Open a GitHub issue or contact qkrskdud0805@gmail.com.
