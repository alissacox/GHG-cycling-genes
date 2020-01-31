Here is the pipeline for analyzing raw Illumina MiSeq reads of functional gene amplicons! These samples came from the surface soils (0-~25cm) of three drainfield types (see [Metadata file](https://github.com/alissacox/GHG-cycling-genes/blob/master/R_code/raw_data/200113_AHC_sequencing_sample_GHG_metadata.txt)). We extracted the DNA from these soils and amplified extracted DNA for both *pmoA* (particulate methane monooxygenase) and *nosZ* (nitrous oxide reductase). The run that generated these sequences had the two different gene amplicons (*pmoA* & *nosZ*) pooled equimolarly per sample. 
* Samples AHC90-96 contain *nosZ* only.

# Start QIIME2
* Need to have installed [QIIME2 version 2019.10](https://docs.qiime2.org/2019.10/install/native/) natively installed via [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Open a Linux terminal and type:
```
conda activate qiime2-2019.10
# You can test to make sure the QIIME2 is running propery by running: 
qiime --help
# To shut down QIIME2 @ the end of your analysis...
conda deactivate qiime2-2019.10
```
# Read in Miseq Fastq files (with embedded barcodes)
## A little 'homework'
* based on [QIIME2 import tutorial](https://docs.qiime2.org/2019.10/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq) and [processing Illumina data tutorial](http://qiime.org/1.3.0/tutorials/processing_illumina_data.html)
* Need to figure out format of our fastq files…. 
  * Indexes appear to be the last several BPs (8+8) of the 1st line in the FastQ file (AHC01):
> @M00763:347:000000000-CP6CN:1:1101:13406:1878 1:N:0:ACTCGCTA+TCGACTAG■	GGTGACTGGGACTTCTGGGTTGACTGGAAGGATCGCCGTATGTGGCCGACGGTTGTGCCGATTCTGGGCGTGACCTTCTGCGCGGCGACGCAGGCGTTTTTCTGGGTGAACTTCCGTCTGCCGTTTGGCGCGGTGTTCGCGGCGCTGGGCCTGCTGATCGGCGAGTGGATCAACCGCTACGTGAACTTCTGGGGTTGGACCTATTTCCCGATCTCGCTGGTGTTCCCGTCGGCTCTGATGGTTCCGGCGATCTGGCTTGACGTGATCCTTCTGCTTTCGGGCTCCTATGTGATCACGGCGA
> +
> CCCCCFGFFGG8CFGGFCECFE@CFFFFGGDG9FCCCCFEGGGGGG,@CGCFFGDGGGGGGGGGGGGGGGFGGGGGGDAFFGGGGGEGGGCGGGGB=>=F@FDFCGGGGGGGGFGGFGFFGFG>FDGG@FFGEFGECCEEGG>F7*CFGGFDDGFCC?FGGGGGGEFGF?EFG8FC8CGEGFF;;<E?FGGGGGGCED??FFGGGGGFGGGCEFGF>:E>5CFGFGFGDDDGG8=6?77*69@7EBD:>>>C54E465770/**77>*92>7C)))/6>EFE)48CC7+824::<F7(49(
 * it appears we have Casava 1.8 demultiplexed dual-end read formatted files based on file name & contents of the QIIME casava example -- only difference is that our 1st line contains the indexers, but QIIME just ignores that.
## Start Analysis in QIIME2: read in files
Change directory to QIIME2 folder containing folder of sequences:
```
cd /mnt/c/Users/xlibb/Desktop/QIIME2
```
* If your files are in multiple subfolders, QIIME2 won't be able to find them. If this is the case, move all the (fastq) files from sub-directories into AHC1-43_90-96 folder:
```
find -type f -print0 | xargs -0 mv -t AHC1-43_90-96
```
Tell QIIME to import sequences & demultiplex them:
```
qiime tools import \
	  --type 'SampleData[PairedEndSequencesWithQuality]' \
	  --input-path AHC1-43_90-96 \
	  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
	  --output-path demux-paired-end.qza
qiime demux summarize \
	  --i-data demux-paired-end.qza \
	  --o-visualization demux-paired-end.qzv
```
# Remove Primers and separate the pooled genes from one another -- use ['cutadapt'](https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage) pre read-joining
* we are doing this step to separate the *pmoA* and *nosZ* sequences from one another, because they were pooled in the MiSeq run. DADA2 includes parameters to trim off primers, but then we'd still have the genes mixed together which causes problems for DADA2...
* the sequencing overhangs are already removed in “raw reads” from MiSeq run
* In our sequences, the adapters should NOT be linked in raw reads because raw reads of 300bps are shorter than any of our target gene sequences (>300 bps). So need to find each primer and its reverse complement separately. 
* Cutadapt can take ALL wildcard bases. [Useful tool](http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html) to auto-create reverse complements
## To separate *pmoA* sequences from "all" the sequences using the primers we used for initial PCR (still embedded in the raw fasta files we imported from the MiSeq run)
* Make a pmoA folder to keep the rest of the pmoA analyis in - this will make things less confusing later on since some of the outputs will have the same names as the outputs from the same nosZ analyses!
Because our PCR amplicons were sheared to 300bp length for the MiSeq run, our sequences will contain both the original primers and their reverse complements. So we need to trim the primer off beginning of each strand & the reverse complement of the other primer at the end of each strand
* Our *pmoA* primers:
	* F primer (A189gcF): GGNGACTGGGACTTCTGG
		* Reverse complement = CCAGAAGTCCCAGTCNCC
	* R primer (mb661R): CCGGMGCAACGTCYTTACC
		* Reverse complement = GGTAARGACGTTGCKCCGG
* Cutadapt terminology:
	* -...adapter-* = the sequence found at 3’ END of sequence (reverse complement to opposite primer)
	* -..front-* = sequence found at 5’ START of sequence (actual primer)
To cut the the *pmoA* primers from the samples and only keeping reads that were trimmed:
```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux-paired-end.qza \
  --p-front-f GGNGACTGGGACTTCTGG \
  --p-adapter-f GGTAARGACGTTGCKCCGG \
  --p-front-r CCGGNGCAACGTCNTTACC \
  --p-adapter-r CCAGAAGTCCCAGTCNCC \
  --p-cores 8 \
  --p-discard-untrimmed \
  --p-overlap 10 \
  --p-minimum-length 150 \
  --p-error-rate 0.3 \
  --o-trimmed-sequences pmoA/trimmed-pmoA-unjoined-seqs.qza
# This step takes ~10 mins to run
qiime demux summarize \
  --i-data pmoA/trimmed-pmoA-unjoined-seqs.qza \
  --o-visualization pmoA/trimmed-pmoA-unjoined-seqs.qzv
```
Based on the output, this seems to have worked:
* <100 sequences show up in AHC 90-96, which is good because these samples only contain *nosZ*, and <100 is an order of magnitude below the next sample which actually has pmoA (showing ~1K sequences). We’ll see if they get filtered out later…
* The import & demux step resulted in 1.3M sequences across all samples, after this primer filtering step only 309K *pmoA* sequences remain. Seems reasonable, since *pmoA* is probably somewhat rare! 
# nosZ
* Make nosZ folder to keep this output in to make things less confusing!
* Our *nosZ* primers:
	* F primer (nosZ1F): CGYTGTTCMTCGACAGCCAG
		* Reverse complement: CTGGCTGTCGAKGAACARCG
	* R primer (nosZ1662R): CGSACCTTSTTGCCSTYGCG
		* Reverse complement: CGCRASGGCAASAAGGTSCG

Code for cutting the nosZ primers from the samples and only keeping reads that were trimmed:
```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux-paired-end.qza \
  --p-front-f CGYTGTTCMTCGACAGCCAG \
  --p-adapter-f CGCRASGGCAASAAGGTSCG \
  --p-front-r CGSACCTTSTTGCCSTYGCG \
  --p-adapter-r CTGGCTGTCGAKGAACARCG \
  --p-cores 8 \
  --p-discard-untrimmed \
  --p-overlap 10 \
  --p-minimum-length 150 \
  --p-error-rate 0.3 \
  --o-trimmed-sequences nosZ/trimmed-nosZ-unjoined-seqs.qza
qiime demux summarize \
  --i-data nosZ/trimmed-nosZ-unjoined-seqs.qza \
  --o-visualization nosZ/trimmed-nosZ-unjoined-seqs.qzv
```
Based on the outputs, this seems to have worked nicely too:
* The import & demux step resulted in 1.3M sequences, after this primer filtering step we ended up with ~1M *nosZ* sequences. YES! That adds up nicely to 1.3M with the 309K *pmoA* sequences! #winning. Longest seq = 420, 20% @ 414. Trim @ 414

# Dada2: join reads, quality filter and denoise
* Modified from [Atacama soils tutorial](https://docs.qiime2.org/2019.10/tutorials/atacama-soils/#atacama-demux)

First, we need to assess the quality of the primer-trimmed demultiplexed sequences for each gene. Then we can denoise each gene's sequences with Dada2

## pmoA trimmed sequences (cutadapt outputs from step above)
* We already trimmed the primers with 'cutadabt', so we don't need to trim primers off in this step.
* The ends of reads need to be trimmed based on Q plots of trimmed sequences (trimmed-pmoA-unjoined-seqs.qzv):
	* Forward read quality drops after ~265 bases
	* Reverse read quality drops after ~195 bases
```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-pmoA-unjoined-seqs.qza \
  --p-trunc-len-f 265 \
  --p-trunc-len-r 195 \
  --p-chimera-method pooled \
  --p-n-threads 0 \
  --p-max-ee-f 20 \
  --p-max-ee-r 20 \
  --o-table dada2-trimmed-pmoA-table.qza \
  --o-representative-sequences dada2-trimmed-pmoA-rep-seqs.qza \
  --o-denoising-stats dada2-trimmed-pmoA-denoising-stats.qza
qiime feature-table summarize \
  --i-table dada2-trimmed-pmoA-table.qza \
  --o-visualization dada2-trimmed-pmoA-table.qzv \
  --m-sample-metadata-file 191221_AHC_sequencing_sample_GHG_metadata.txt
qiime feature-table tabulate-seqs \
  --i-data dada2-trimmed-pmoA-rep-seqs.qza \
  --o-visualization dada2-trimmed-pmoA-rep-seqs.qzv
qiime metadata tabulate \
  --m-input-file dada2-trimmed-pmoA-denoising-stats.qza \
  --o-visualization dada2-trimmed-pmoA-denoising-stats.qzv
```
This step results in 'only' 161 unique *pmoA* sequences, but that's on par with other studies. Samples AHC>90 don't have any sequences (makes sense because these don't have *pmoA*), and two other samples (AHC 36 & 42) don't appear to have any *pmoA* sequences either (must not have had a lot of pmoA?).

Need to filter out ‘empty’ samples so we can do some of the downstream analyses without errors:
```
qiime feature-table filter-samples \
  --i-table dada2-trimmed-pmoA-table.qza \
  --p-min-features 1 \
  --o-filtered-table filtered-dada2-trimmed-pmoA-table.qza
qiime metadata tabulate \
  --m-input-file filtered-dada2-trimmed-pmoA-table.qza \
  --o-visualization filtered-dada2-trimmed-pmoA-table.qzv
  ```
This resulted in removing samples AHC>90 and samples 36 & 42 ...

To export these data (feature/ASV/"OTU" table) for use in R with the phyloseq package, or to look at them in a spreadsheet program:
```
qiime tools export   \
  --input-path filtered-dada2-trimmed-pmoA-table.qza  \
  --output-path exp-filt-pmoA-table-dada2
# To convert this file to a .tsv:
biom convert -i exp-filt-pmoA-table-dada2/feature-table.biom -o exp-filt-pmoA-table-dada2/pmoA-filt-table-dada2.tsv --to-tsv
```

## nosZ trimmed sequences (outputs from cutadapt step)
* Again, no trimming of primers required, trimming ends based on Q plots of trimmed sequences:
```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-nosZ-unjoined-seqs.qza \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 184 \
  --p-chimera-method pooled \
  --p-n-threads 0 \
  --p-max-ee-f 15 \
  --p-max-ee-r 15 \
  --o-table dada2-trimmed-nosZ-table.qza \
  --o-representative-sequences dada2-trimmed-nosZ-rep-seqs.qza \
  --o-denoising-stats dada2-trimmed-nosZ-denoising-stats.qza
qiime feature-table summarize \
  --i-table dada2-trimmed-nosZ-table.qza \
  --o-visualization dada2-trimmed-nosZ-table.qzv \
  --m-sample-metadata-file 191221_AHC_sequencing_sample_GHG_metadata.txt
qiime feature-table tabulate-seqs \
  --i-data dada2-trimmed-nosZ-rep-seqs.qza \
  --o-visualization dada2-trimmed-nosZ-rep-seqs.qzv
qiime metadata tabulate \
  --m-input-file dada2-trimmed-nosZ-denoising-stats.qza \
  --o-visualization dada2-trimmed-nosZ-denoising-stats.qzv
```
After this step, we are left with the total # of unique sequences = 7539, with a total frequency = 673,960… on all 50 samples.  Most common sequence comes back as nosZ on BLAST.

To export this feature table...
```
qiime tools export   \
  --input-path dada2-trimmed-nosZ-table.qza  \
  --output-path exp-nosZ-table-dada2
# Convert biom feature table to .tsv
biom convert -i exp-nosZ-table-dada2/feature-table.biom -o nosZ-table-dada2.tsv --to-tsv

