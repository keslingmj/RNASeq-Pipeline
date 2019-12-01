## This is the configuration file for running hisat2 and other scripts

NUMCPUS=3
HISAT2=$(which hisat2)
STRINGTIE=$(which stringtie)
GFFREAD=$(which gffread)
SALMON=$(which salmon)
SAMTOOLS=$(which samtools)
BASEDIR="chrX_data"
FASTQLOC="$BASEDIR/samples"
GENOMEIDX="$BASEDIR/indexes/chrX_tran"
GTFFILE="$BASEDIR/genes/chrX.gtf"
GENOMELOC="$BASEDIR/genome"
TEMPLOC="./tmp"
STUDY="geuvadis"
