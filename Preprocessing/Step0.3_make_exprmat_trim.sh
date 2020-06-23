## Run on HiPerGator
module load rsem

# This was for every dataset:
RUN_DIR=EQUALIZATION_PAPER/RSEMDATA/EQ
cd $RUN_DIR
rsem-generate-data-matrix trimmed_*.genes.results > exprmat_EQ_trimmed_codingRef.txt
