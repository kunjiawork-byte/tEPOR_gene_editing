#!/usr/bin/env bash
# tEPOR gene editing study: fusion analysis
# Analysis parameters are retained from the original study script.
# Set the path variables below before running; requires Bash on Linux.
set -euo pipefail

############################
# Configuration (edit paths or set environment variables)
############################
THREADS=${THREADS:-16}

# K1-K24 paired-end input and output
DATADIR="${DATADIR:-/path/to/RNAseqK1-24/trimmed}"
OUTROOT="${OUTROOT:-/path/to/RNAseqK1-24/fusion_out}"

# Reference genome and Arriba resources
STAR_INDEX="${STAR_INDEX:-/path/to/star_index}"
GENOME_FA="${GENOME_FA:-/path/to/GRCh38.primary_assembly.genome.fa}"
GTF="${GTF:-/path/to/gencode.v40.annotation.gtf}"
BLACKLIST="${BLACKLIST:-/path/to/arriba_blacklist.tsv}"
KNOWN="${KNOWN:-/path/to/arriba_known_fusions.tsv}"
DOMAINS="${DOMAINS:-/path/to/arriba_protein_domains.gff3.gz}"
# Cytobands are used only for plotting
CYTO="${CYTO:-/path/to/arriba_cytobands.tsv.gz}"

SAMPLE_GLOB="K*_1_paired.fq.gz"
SKIP_SAMPLES=${SKIP_SAMPLES:-""}
SKIP_IF_DONE=${SKIP_IF_DONE:-1}

# Export the top N candidates per sample (IGV loci and optional PDF)
TOP_N=${TOP_N:-3}
############################

log(){ echo "[$(date '+%F %T')] $*" >&2; }

need=(STAR arriba samtools awk)
for b in "${need[@]}"; do
  command -v "$b" >/dev/null 2>&1 || { echo "ERROR: $b is not on PATH"; exit 1; }
done

mkdir -p "$OUTROOT"
cd "$DATADIR"

# Batch alignment and fusion calling
mapfile -t R1_LIST < <(ls ${SAMPLE_GLOB} 2>/dev/null | sort -V || true)
((${#R1_LIST[@]})) || { echo "No files found matching ${SAMPLE_GLOB}"; exit 1; }

for R1 in "${R1_LIST[@]}"; do
  SAMPLE="${R1%_1_paired.fq.gz}"
  R2="${SAMPLE}_2_paired.fq.gz"

  if [[ -n "$SKIP_SAMPLES" ]] && printf '%s\n' $SKIP_SAMPLES | grep -qx "$SAMPLE"; then
    log "Skipping $SAMPLE"; continue
  fi
  [[ -s "$R1" && -s "$R2" ]] || { log "WARNING: Missing paired FASTQ files：$R1 / $R2"; continue; }

  OUTDIR="${OUTROOT}/${SAMPLE}"
  mkdir -p "$OUTDIR"

  if (( SKIP_IF_DONE )) && [[ -s "$OUTDIR/${SAMPLE}_fusions.tsv" ]]; then
    log "Results already exist; skipping: $SAMPLE"; continue
  fi

  log ">>> Processing $SAMPLE"
  STAR --runThreadN "$THREADS" \
       --genomeDir "$STAR_INDEX" \
       --readFilesIn "$R1" "$R2" \
       --readFilesCommand zcat \
       --outSAMtype BAM SortedByCoordinate \
       --twopassMode Basic \
       --chimSegmentMin 12 \
       --chimJunctionOverhangMin 12 \
       --chimOutType Junctions SeparateSAMold \
       --alignSJDBoverhangMin 10 \
       --outFileNamePrefix "$OUTDIR/${SAMPLE}_"

  BAM="$OUTDIR/${SAMPLE}_Aligned.sortedByCoord.out.bam"
  CHIM_SAM="$OUTDIR/${SAMPLE}_Chimeric.out.sam"
  samtools index -@ "$THREADS" "$BAM"

  arriba \
    -x "$BAM" \
    -c "$CHIM_SAM" \
    -o "$OUTDIR/${SAMPLE}_fusions.tsv" \
    -O "$OUTDIR/${SAMPLE}_fusions.discarded.tsv" \
    -a "$GENOME_FA" \
    -g "$GTF" \
    -b "$BLACKLIST" \
    -k "$KNOWN" \
    -p "$DOMAINS"

  log "Completed $SAMPLE"
done

# Summarize and filter fusion calls
log ">>> Generating summaries and filtered tables"
FIXED="$OUTROOT/fusion_summary.fixed.tsv"
printf "sample\tgene1\tgene2\tjunction_reads\tspanning_pairs\tconfidence\treading_frame\ttype\n" > "$FIXED"

find "$OUTROOT" -maxdepth 2 -type f -name '*_fusions.tsv' | sort -V | while read -r f; do
  s=$(basename "$f" _fusions.tsv)
  awk -v sample="$s" 'BEGIN{FS=OFS="\t"}
    NR==1{
      for(i=1;i<=NF;i++){ h[tolower($i)]=i }
      gi=h["gene1"]; gj=h["gene2"];
      sr1=h["split_reads1"]; if(!sr1) sr1=h["junctionreads1"]; if(!sr1) sr1=h["jr1"];
      sr2=h["split_reads2"]; if(!sr2) sr2=h["junctionreads2"]; if(!sr2) sr2=h["jr2"];
      dm =h["discordant_mates"]; if(!dm) dm=h["spanningfrags"]; if(!dm) dm=h["spanning_pairs"]; if(!dm) dm=h["spanning_frags"];
      cf =h["confidence"];
      rf =h["reading_frame"]; if(!rf) rf=h["frame"];
      ty =h["type"];
      next
    }
    {
      g1=(gi?$gi:""); g2=(gj?$gj:"");
      v1=(sr1?($sr1+0):0); v2=(sr2?($sr2+0):0); jr=v1+v2;
      sp=(dm?($dm+0):0);
      c =(cf?$cf:""); rfv=(rf?$rf:""); t=(ty?$ty:"");
      print sample, g1, g2, jr, sp, c, rfv, t
    }' "$f" >> "$FIXED"
done

awk 'BEGIN{FS=OFS="\t"} NR==1 || ( ($4+0)>=2 && ($5+0)>=2 && tolower($6)!~/^(low|poor)$/ )' \
  "$FIXED" > "$OUTROOT/fusion_summary.strict.tsv"

awk 'BEGIN{FS=OFS="\t"} NR==1 || ( ( ($4+0)>=2 || ($5+0)>=2 ) && tolower($6)!~/^(low|poor)$/ )' \
  "$FIXED" > "$OUTROOT/fusion_summary.sensitive.tsv"

awk 'BEGIN{FS=OFS="\t"} NR==1 || ( ($4+0)>=2 )' \
  "$FIXED" > "$OUTROOT/fusion_summary.relaxed.tsv"

awk 'BEGIN{FS=OFS="\t"}
  NR>1{ s=$1; jr=$4+0; sp=$5+0; conf=tolower($6);
        tot[s]++; if(jr>=2&&sp>=2&&conf!~/(^low$|^poor$)/) str[s]++;
        if((jr>=2||sp>=2)&&conf!~/(^low$|^poor$)/) sen[s]++;
        if(jr>=2) rel[s]++ }
  END{ print "sample","total","strict","sensitive","relaxed";
       for(s in tot) print s, (tot[s]+0), (str[s]+0), (sen[s]+0), (rel[s]+0) }' \
  "$FIXED" | sort -k2,2nr > "$OUTROOT/fusion_per_sample_counts.tsv"

# Export top N candidates per sample
log ">>> Exporting top ${TOP_N} candidates per sample: IGV loci and optional PDFs"

if command -v draw_fusions.R >/dev/null 2>&1; then
  HAVE_DRAW=1
else
  HAVE_DRAW=0
  log "WARN: draw_fusions.R was not found; exporting loci only. Add the Arriba plotting script to PATH to generate PDFs."
fi

# Iterate over sample directories
find "$OUTROOT" -maxdepth 1 -type d -name 'K*' | sort -V | while read -r SD; do
  SAMPLE=$(basename "$SD")
  FUS="$SD/${SAMPLE}_fusions.tsv"
  BAM="$SD/${SAMPLE}_Aligned.sortedByCoord.out.bam"
  [[ -s "$FUS" ]] || { log "WARN: No fusion table; skipping $SAMPLE"; continue; }

  # Read the header
  HEADER=$(head -n1 "$FUS")

  # Rank by split-read support, spanning pairs, then confidence
  awk 'BEGIN{FS=OFS="\t"}
    NR==1{
      for(i=1;i<=NF;i++){ hl[tolower($i)]=i; real[i]=$i }
      sr1=hl["split_reads1"]; if(!sr1) sr1=hl["junctionreads1"]; if(!sr1) sr1=hl["jr1"];
      sr2=hl["split_reads2"]; if(!sr2) sr2=hl["junctionreads2"]; if(!sr2) sr2=hl["jr2"];
      dm =hl["discordant_mates"]; if(!dm) dm=hl["spanningfrags"]; if(!dm) dm=hl["spanning_pairs"]; if(!dm) dm=hl["spanning_frags"];
      cf =hl["confidence"];
      next
    }
    {
      jr=((sr1?$(sr1):0)+0)+((sr2?$(sr2):0)+0)
      sp=(dm?$(dm)+0:0)
      conf=tolower((cf?$(cf):""))
      w=(conf=="high"?10:(conf=="medium"?5:0))
      score = jr*1000 + sp*100 + w
      print score"\t"$0
    }' "$FUS" \
  | sort -k1,1nr \
  | awk -v N="$TOP_N" -v header="$HEADER" 'BEGIN{FS=OFS="\t"; print header} NR<=N { $1=""; sub(/^\t/,""); print }' \
  > "$SD/${SAMPLE}_fusions.top${TOP_N}.tsv"

  # Export the original breakpoint strings for IGV review
  if awk 'BEGIN{FS=OFS="\t"} NR==1{
            for(i=1;i<=NF;i++){ lc=tolower($i); if(lc=="breakpoint1") b1=i; if(lc=="breakpoint2") b2=i }
            if(!b1||!b2){ exit 1 }
            next
          }
          { print $b1"-150"; print $b2"-150" }' \
      "$SD/${SAMPLE}_fusions.top${TOP_N}.tsv" > "$SD/${SAMPLE}_top${TOP_N}_loci.txt"; then
    :
  else
    log "WARN: $SAMPLE is missing breakpoint1/2 columns; loci were not generated"
  fi

  # Draw PDFs when draw_fusions.R and an indexed BAM are available
  if (( HAVE_DRAW )) && [[ -s "$BAM" && -s "$BAM.bai" ]]; then
    CMD=(draw_fusions.R
         --fusions="$SD/${SAMPLE}_fusions.top${TOP_N}.tsv"
         --alignments="$BAM"
         --output="$SD/${SAMPLE}_fusions.top${TOP_N}.pdf"
         --annotation="$GTF"
         --proteinDomains="$DOMAINS")
    [[ -s "$CYTO" ]] && CMD+=(--cytobands="$CYTO")
    "${CMD[@]}" || log "WARN: draw_fusions.R failed: $SAMPLE"
  fi
done

log "Completed. Output directory: $OUTROOT"
log "Main summary: $OUTROOT/fusion_summary.fixed.tsv"
log "Per-sample counts: $OUTROOT/fusion_per_sample_counts.tsv"
log "Top N exports in each sample directory: *_fusions.top${TOP_N}.tsv / *_top${TOP_N}_loci.txt / (optional) *_fusions.top${TOP_N}.pdf"
