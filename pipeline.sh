#!/usr/bin/env bash

# wrapper script for running hisat2, stringtie, gffreads, and salmon
# along with CONFIG.sh parameters
# this script is based on the rnaseq_pipeline.sh script from the Perteas in the
# Salzberg lab, but with many changes


# Author: Michael Kesling

######
# check that I'm in the correct location relative to the $BASEDIR defined in CONFIG.sh

##############
# set the output directory to be pwd unless given as 1st argument
# of this script
OUTDIR="."
if [[ "$1" ]]; then
 if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
  exit 1
 fi
 OUTDIR=$1
fi

###############
# check for the existence of the configuration file
if [[ ! -f ./CONFIG.sh ]]; then
  usage
  echo "Error: configuration file (CONFIG.sh) is missing!"
  exit 1
fi

###############
# source the configuration file, which defined many important variables
source ./CONFIG.sh

###############
# set various output directory variables
SCRIPTARGS="$@"
ALIGNLOC=./hisat2
TCLOC=./transcriptome
LOGFILE=./run.log
WRKDIR=$(pwd -P)
errprog=""
STIELOC=./stringtie2

###############
# check that the needed software can be found by this script
# and that they have executable permission
if [[ ! -x $SAMTOOLS ]]; then
    errprog="samtools"
fi
if [[ ! -x $HISAT2 ]]; then
    errprog="hisat2"
fi
if [[ ! -x $GFFREAD ]]; then
    errprog="gffread"
fi
if [[ ! -x $SALMON ]]; then
    errprog="salmon"
fi
if [[ ! -x $STRINGTIE ]]; then
    errprog="stringtie"
fi
if [[ "$errprog" ]]; then
  echo "ERROR: $errprog program not found, please edit the configuration script."
  exit 1
fi

###############
# ensure that the needed directories have been created
if [[ $OUTDIR != "." ]]; then
  mkdir -p $OUTDIR
  cd $OUTDIR
fi

for d in "$TEMPLOC" "$ALIGNLOC" "$TCLOC" "$STIELOC" ; do
 if [ ! -d $d ]; then
    mkdir -p $d
 fi
done

##############
#determine samtools version
newsamtools=$( ($SAMTOOLS 2>&1) | grep 'Version: 1\.')


################
## list of samples
## (only paired reads, must follow _1.*/_2.* file naming convention)
reads1=(${FASTQLOC}/*_1.*)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_1./_2.}")


# perform a check that all samples are matched


# main script block
pipeline() {

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $SCRIPTARGS

for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    sample="${reads1[$i]%%.*}"
    sample="${sample%_*}"
    stime=`date +"%Y-%m-%d %H:%M:%S"`
    echo "[$stime] Processing sample: $sample"
    echo [$stime] "   * Alignment of reads to genome (HISAT2)"


    $HISAT2 -p $NUMCPUS --dta -x ${GENOMEIDX} \
	    -1 ${FASTQLOC}/${reads1[$i]} \
	    -2 ${FASTQLOC}/${reads2[$i]} \
	    -S ${TEMPLOC}/${sample}.sam 2>${ALIGNLOC}/${sample}.alnstats

    echo [`date +"%Y-%m-%d %H:%M:%S"`] "Sorting and binary conversion of SAM files"


    if [[ "$newsamtools" ]]; then
	$SAMTOOLS view -S -b ${TEMPLOC}/${sample}.sam | \
	    $SAMTOOLS sort -@ $NUMCPUS -o ${ALIGNLOC}/${sample}.bam - # -O BAM
    else
	$SAMTOOLS view -S -b ${TEMPLOC}/${sample}.sam | \
	    $SAMTOOLS sort -@ $NUMCPUS - ${ALIGNLOC}/${sample}
    fi
	    
    $SAMTOOLS index ${ALIGNLOC}/${sample}.bam
    $SAMTOOLS flagstat ${ALIGNLOC}/${sample}.bam

    echo "..removing intermediary files"
    \rm ${TEMPLOC}/${sample}.sam


    # Run StringTie to ID relevant Tcs.
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Assemble transcripts (StringTie)"
    $STRINGTIE -p $NUMCPUS -G ${GTFFILE} -o ${ALIGNLOC}/${sample}.gtf \
	       -l ${sample} ${ALIGNLOC}/${sample}.bam

done

###########
### Stringtie-merge all .gtf files into a single .gtf file
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Merge all transcripts (StringTie)"
ls ${ALIGNLOC}/*.gtf > ${ALIGNLOC}/mergelist.txt     # (rm 'ls -l' from original)

$STRINGTIE --merge -p $NUMCPUS -G ${GTFFILE} \
	   -o stringtie_merged.gtf ${ALIGNLOC}/mergelist.txt

###########
# OPTIONAL: re-run StringTie on all transcripts--including newly found ones
# I do this in order to compare results with those of Salmon

echo [`date + "%Y-%m-%d %H:%M:%S"`] "#> Estimate abundance for each sample (StringTie)"

for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    sample="${reads1[$i]%%.*}"
    dsample="${sample%%_*}"
    sample="${sample%_*}"
    if [ ! -d ${STIELOC}/${dsample} ]; then
       mkdir -p ${STIELOC}/${dsample}
    fi
    $STRINGTIE -e -B -p $NUMCPUS -G stringtie_merged.gtf \
    -o ${STIELOC}/${dsample}/${dsample}.gtf ${ALIGNLOC}/${sample}.bam
done






###########
# create fasta sequences of all identified and reference transcript from merge
# GENE NAMES NOT INFORMATIVE                                                    #######################################
$GFFREAD stringtie_merged.gtf -g ${GENOMELOC}/chrX.fa -w ${TCLOC}/${STUDY}.fa #set genome.fa in CONFIG?



###########
# index salmon
# want to capture error in a file
# may want to add the --gencode option to keep track of version
$SALMON index -t ${TCLOC}/${STUDY}.fa -i ${TCLOC}/${STUDY}.idx


# quantitate using salmon
# want to capture error in file
# may want to add --seqBias option
# in this option, we are aligning and quantitating via a single command


################
## Re-create list of reads1 and reads2                        ### MAY NOT NEED THESE STATEMENTS
reads1=(${FASTQLOC}/*_1.*)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_1./_2.}")


echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Salmon mapping and quantification: "


for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    sample="${reads1[$i]%%.*}"
    sample="${sample%_*}"
    $SALMON quant -i ${TCLOC}/${STUDY}.idx -l A \
						    -1 ${FASTQLOC}/${reads1[$i]} \
						    -2 ${FASTQLOC}/${reads2[$i]} \
						    -p $NUMCPUS --validateMappings -o salmonQuants/${sample}.quant
done


###############
# Sum transcript-level TPM for every quant.sf file to create gene-level TPMs and put everything
# into a single matrix.
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Merging Salmon TPM values at gene levels for all samples "

cat hisat2/mergelist.txt | sed 's/.*.\//salmonQuants\//g' | sed 's/gtf$/quant\/quant.sf/' > salmonQfiles.txt
Rscript post-process-salmon.R TRUE


}

# salmon outputs TPM, but edgeR and DESeq2 want Count data.

pipeline 2>&1 | tee $LOGFILE

# Using the tximport package, you can import salmon’s transcript-level quantifications and optionally aggregate them to the gene level for gene-level differential expression analysis. You can read more about how to import salmon’s results into DESeq2 by reading the tximport section of the excellent DESeq2 vignette. For instructions on importing for use with edgeR or limma, see the tximport vignette. For preparing salmon output for use with sleuth, see the wasabi package.

# will want to use sample-to-sample normalization via TMM
