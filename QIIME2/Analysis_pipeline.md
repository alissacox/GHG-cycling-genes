Here is the pipeline for analyzing raw Illumina MiSeq reads of functional gene amplicons! The run that generated these genes had two different gene amplicons (*pmoA* & *nosZ*) pooled equimolarly per sample. 
* Samples AHC90-96 contain *nosZ* only.

# Start QIIME2
* Need to have installed [QIIME2 version 2019.10](https://docs.qiime2.org/2019.10/install/native/) natively installed via [miniconda](https://docs.conda.io/en/latest/miniconda.html)
Open a Linux terminal and type:
```
conda activate qiime2-2019.10
conda deactivate qiime2-2019.10 (to shut down)
# You can test to make sure the install worked by running: 
qiime --help
```
# Read in Miseq Fastq files (with embedded barcodes)
## A little 'homework'
* based on [QIIME2 Tutorial](https://docs.qiime2.org/2019.10/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq http://qiime.org/1.3.0/tutorials/processing_illumina_data.html)
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
* Our primers:
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
# This step Takes ~10 mins to run
qiime demux summarize \
  --i-data pmoA/trimmed-pmoA-unjoined-seqs.qza \
  --o-visualization pmoA/trimmed-pmoA-unjoined-seqs.qzv
```
Based on the output, this seems to have worked. <100 sequences show up in AHC 90-96, which is good because these samples only contain *nosZ*, and <100 is an order of magnitude below the next sample which actually has pmoA (showing ~1K sequences). We’ll see if they get filtered out later…

