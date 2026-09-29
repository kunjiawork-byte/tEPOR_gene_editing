# Targeted amplicon sequencing analysis

Recovered from the supplied NGS code notes. Original sequences, deepseq demultiplexing parameters, quantification windows, conversion options, and process counts are retained. Repeated variable assignments and mutually exclusive table-writing blocks have been separated into explicit choices. Input paths are configurable, and file enumeration handles spaces.

## Workflow

1. `01_demultiplex_r2.sh`: provide `BASE_DIR`, containing `KJ28-D4-*` sample folders (or set `SAMPLE_PATTERN`). Each folder needs exactly one `*_R1_001.fastq.gz`, one `*_R2_001.fastq.gz`, and the original sample-specific `barcodes.fasta`. Retains cutadapt `-e 0.1 -q 20,20` and anchored R2 barcode matching. R2 is passed first and output as `.R2`; its paired R1 is output as `.R1`.
2. `02_fastqc_multiqc.sh`: set `QC_DIR` to the directory tree containing demultiplexed FASTQs. Retains FastQC 8 threads and MultiQC output `multiqc_report`.
3. `03_make_batch_input.py`: choose one target and either `r1` or `paired`. Run separately for each target-specific FASTQ directory; no filename-based target assignment is inferred. Paired mode skips missing R2 files as in the original notes.
4. `04_run_crispresso_batch.sh`: choose one recorded parameter set. Run only the setting actually used for that experiment, in a separate analysis directory if comparing settings.

Example: BCL11A CBE, R1-only (replace paths):

```bash
BASE_DIR=/path/to/00_fastq bash 01_demultiplex_r2.sh
QC_DIR=/path/to/00_fastq bash 02_fastqc_multiqc.sh
python3 03_make_batch_input.py --input-dir /path/to/BCL11A1620 --target BCL11A_CBE --reads r1 --output CRISPResso_batch_input.txt
bash 04_run_crispresso_batch.sh cbe
```

For paired reads use `--reads paired`. Batch table output contains absolute FASTQ paths resolved at runtime, not internal server paths.

## Target definitions

| Target option | Amplicon | Guide |
| --- | --- | --- |
| HBG_CBE | HBG2 | sgHBG-1 |
| EPOR_sg1 | EPOR | EPOR-sg1 |
| EPOR_sg2 | EPOR | EPOR-sg2 |
| BCL11A_CBE | BCL11A | BCL11A-1620 |
| BCL11A_Cas9 | BCL11A | BCL11A-1617 |

The HBG CBE batch block supplied HBG2 alone; HBG1 is retained as a reference sequence but is not silently added to that analysis. The BCL11A Cas9 pairing is taken from the BCL11A amplicon and its explicitly labeled Cas9 guide in the notes; no separate Cas9 batch-generation block was supplied. EPOR guide 1157 and HBG-175 are also preserved as reference-only guides.

`amplicons.fasta` and `guides.fasta` are reference sequences, NOT sample barcode files. The sample-specific barcode sequences were not supplied and must be provided separately.

## Recorded CRISPResso modes

| Mode | Center | Size | Additional options |
| --- | --- | --- | --- |
| cbe | -10 | 20 | --base_editor_output |
| cbe_ga | -10 | 20 | --base_editor_output --conversion_nuc_from G --conversion_nuc_to A |
| cas9_cut | -3 | 1 | none |
| cas9_wide | -10 | 20 | none |

All modes retain `--n_processes 4`. Both Cas9 windows were present as alternatives; the supplied notes do not establish which was used for each final figure. The source CBE command left conversion nucleotides at installed-version defaults; the cleaned `cbe` mode does the same.

## Excluded exploratory snippets and missing steps

The opening generic examples use cutadapt `-e 1` or dual-barcode `-e 0.15 --no-indels`. They are distinct from the later deepseq block and were not substituted into it. The generic R1 demultiplex example also labels its first-read output `.2`; it was not adopted. Incomplete CRISPResso/CRISPRessoPooled examples were excluded.

The notes mention `calculate_indel.py` and `calculate_CT_conversion.py` without their contents. Their postprocessing and any final plotting/statistical analysis cannot be reconstructed from this file.

Requires Bash, cutadapt, FastQC, MultiQC, Python 3, and CRISPResso2/CRISPRessoBatch. Record the exact installed versions used in the study; they were not included in the supplied notes. Activate the appropriate environment for each stage.

## Validation

Prepared scripts were syntax-checked; generated batch rows and recorded CRISPResso arguments were checked with small fixtures. No original sequencing data were rerun.
