# GHG-cycling-genes
Pipeline to analyze functional gene amplicons from Illumina MiSeq runs with QIIME 2 and other commandline tools - from n00bs for n00bs

The raw sequence reads for this analysis can be found in the [Raw_Illumina_Seq_Reads](https://github.com/alissacox/GHG-cycling-genes/tree/master/Raw_Illumina_Seq_Reads) folder.

The [QIIME2](https://github.com/alissacox/GHG-cycling-genes/tree/master/QIIME2) folder contains separate pages for the [script required to run the analysis of Illumina MiSeq amplicon data in QIIME2](https://github.com/alissacox/GHG-cycling-genes/blob/master/QIIME2/Analysis_pipeline.md), directions & script to [create a custom database from NCBI fasta files](https://github.com/alissacox/GHG-cycling-genes/blob/master/QIIME2/Custom_Database_Creation) (since only 16S/18S databases are available 'premade' for QIIME2 users and there was little helpful detail available on the web on how to actually do this), and directions & script to [use the Entrez Direct E-Utils commandline tool to make a QIIME2 compatible taxonomy file](https://github.com/alissacox/GHG-cycling-genes/blob/master/QIIME2/Custom_Database_Taxonomy) for import into QIIME2.

Note: This analysis pipeline requires a working Linux terminal. This pipeline ran successfully on several Windows 10 machines with the [Windows Subsystem for Linux enabled](https://www.windowscentral.com/install-windows-subsystem-linux-windows-10) and using Ubuntu (18.04). These commands also ran without problems on a native Ubuntu (18.04) installation.
