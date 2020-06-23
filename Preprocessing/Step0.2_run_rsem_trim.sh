# Complete pipeline submit script. Run on HiPerGator

## This was done for every dataset!

module load seqtk
module load bowtie
module load rsem


REF_NAME=REF/Homo_sapiens-GRCh38/ref_ensembl_coding/human
CODE_DIR=EQUALIZATION_PAPER/PIPELINE
COLLAB_DIR=EQUALIZATION_PAPER/


RAWDATA_DIR=EQUALIZATION_PAPER/FASTQ/EQ
RSEMDATA_DIR=EQUALIZATION_PAPER/RSEMDATA/EQ

## Do for read 1 and 2:

cd $RAWDATA_DIR
folders=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p) #list of folders
cd $folders

gunzip *_R1_001.fastq.gz
fastfile=$(find *R1*.fastq)
seqtk trimfq -e 27 $fastfile > trimmed_$fastfile
trimfast=trimmed_$fastfile
cd ..
rsem-calculate-expression -p 5 $RAWDATA_DIR/$folders/$trimfast \
          $REF_NAME \
          $RSEMDATA_DIR/$trimfast --append-names --no-bam-output --seed-length 28

cd $RAWDATA_DIR
folders=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p) #list of folders
cd $folders
gunzip *_R2_001.fastq.gz
fastfile=$(find *R2*.fastq)
seqtk trimfq -e 27 $fastfile > trimmed_$fastfile
trimfast=trimmed_$fastfile
cd ..
rsem-calculate-expression -p 5 $RAWDATA_DIR/$folders/$trimfast \
          $REF_NAME \
          $RSEMDATA_DIR/$trimfast --append-names --no-bam-output --seed-length 28


