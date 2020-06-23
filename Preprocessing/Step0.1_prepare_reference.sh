module load rsem

cd REF/Homo_sapiens-GRCh38

rsem-prepare-reference --gtf Homo_sapiens.GRCh38.83.gtf \
Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--gff3-RNA-patterns 'protein_coding' \
--trusted-sources 'ensembl_havana' \
--bowtie \
ref_ensembl_coding/human

