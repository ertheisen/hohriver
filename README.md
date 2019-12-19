# appalachian.hg19v1 container for early NGS pipeline applications
How to run pipeline:

 singularity run -H $HOME:/home/$USER --app [appname] --bind [directory_info] container_name.simg

Path to genome in container needs to be:
 '/genomes/test/Sequence/Bowtie2Index/genome'
-or-
 '/genomes/spike/dm6/Sequence/Bowtie2Index/genome'
-or-
 '/genomes/spike/sacCer3/Sequence/Bowtie2Index/genome'
-or-
 '/genomes/STAR/[STAR index files]'
-or-
 '/genomes/anno/[gtf or gff annotation file]'

For Baker Cluster Users:
 directory to mount for fly spike-in is: /gpfs0/home/gdlessnicklab/share/trails/genomes/Drosophila_melanogaster/UCSC
 directory to mount for yeast spike-in is: /gpfs0/home/gdlessnicklab/share/trails/genomes/Saccharomyces_cerevisiae/UCSC
 directory to mount for E. coli spike-in is: /gpfs0/home/gdlessnicklab/share/trails/genomes/Escherichia_coli_K_12_DH10B/NCBI
 directory to mount for hg19 is: /reference/homo_sapiens/hg19/ucsc_assembly/illumina_download/

The correct bind commands for test and spike are as follows:
 --bind /reference/homo_sapiens/hg19/ucsc_assembly/illumina_download/:/genomes/test
 --bind /gpfs0/home/gdlessnicklab/share/trails/genomes/Escherichia_coli_K_12_DH10B/NCBI:/genomes/spike


To run interactive terminal in container

 singularity shell -H $HOME:/home/$USER --bind [directory_info] appalachian_hg19.simg

To run command from a program in the container

 singularity exec -H $HOME:/home/$USER --bind [directory_info] appalachian_hg19.simg [command information]

 ##Be sure to bind your data, your genome, and any script files you may want to run

Apps:
 cutnrun_full : full pipeline from fastqc to mapped bam, bed, normalized bedgraph, and normalized bigwig with mapping qc using plotFingerprint from deeptools -  bowtie2 alignment assumes 150 bp PE reads and allows dovetail and overlap in alignment
 cutntag : full pipeline from fastqc to mapped bam, bed, normalized bedgraph, and normalized bigwig with mapping qc using plotFingerprint from deeptools -  bowtie2 alignment assumes 150 bp PE reads and allows dovetail and overlap in alignment
 cutnrun_hardtrim : full pipeline from fastqc to mapped bam, bed, normalized bedgraph, and normalized bigwig with mapping qc using plotFingerprint from deeptools - bowtie2 alignment. First step trims 150 bp reads from HiSeq to 60 bp.
 cutntag_hardtrim : full pipeline from fastqc to mapped bam, bed, normalized bedgraph, and normalized bigwig with mapping qc using plotFingerprint from deeptools - bowtie2 alignment. First step trims 150 bp reads from HiSeq to 60 bp.
 fastqc_pre : FastQC pre-trim
 trim_qc : Fastq trimming with Trim Galore and FastQC post trim
 bowtie2_alignment : Alignment for CUT&RUN samples; insert 10 bp-700 bp, assumes paired end 150 bp reads, makes bams/beds/bg/bw with no spike in.
 bams : Generate downstream files for starting at bams using same pipeline - assumes hg19 alignment
 mapqc_fingerprint : Use deeptools to assess enrichment quality

Data should be stored in a "project" directory with a sub-directory called "fastq". The "fastq" subdirectory should have all of your paired end fastq or fastq.gz files. Then the correct bind for data is:
 --bind /path_to/project/:/data
 


