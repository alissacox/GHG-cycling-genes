Here is the pipeline for analyzing raw Illumina MiSeq reads of functional gene amplicons!

#Start QIIME2
* Need to have installed [QIIME2 version 2019.10](https://docs.qiime2.org/2019.10/install/native/) natively installed via [miniconda](https://docs.conda.io/en/latest/miniconda.html)
Open a Linux terminal and type:
```
conda activate qiime2-2019.10
conda deactivate qiime2-2019.10 (to shut down)
# You can test to make sure the install worked by running: 
qiime --help

```
#Read in Miseq Fastq files (with embedded barcodes)
* based on [QIIME Tutorial](https://docs.qiime2.org/2019.10/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq http://qiime.org/1.3.0/tutorials/processing_illumina_data.html )
* Need to figure out format of our fastq files…. 
  * Indexes appear to be the last several BPs (8+8) of the 1st line in the FastQ file (AHC01):
> @M00763:347:000000000-CP6CN:1:1101:13406:1878 1:N:0:ACTCGCTA+TCGACTAG■	GGTGACTGGGACTTCTGGGTTGACTGGAAGGATCGCCGTATGTGGCCGACGGTTGTGCCGATTCTGGGCGTGACCTTCTGCGCGGCGACGCAGGCGTTTTTCTGGGTGAACTTCCGTCTGCCGTTTGGCGCGGTGTTCGCGGCGCTGGGCCTGCTGATCGGCGAGTGGATCAACCGCTACGTGAACTTCTGGGGTTGGACCTATTTCCCGATCTCGCTGGTGTTCCCGTCGGCTCTGATGGTTCCGGCGATCTGGCTTGACGTGATCCTTCTGCTTTCGGGCTCCTATGTGATCACGGCGA
> +
> CCCCCFGFFGG8CFGGFCECFE@CFFFFGGDG9FCCCCFEGGGGGG,@CGCFFGDGGGGGGGGGGGGGGGFGGGGGGDAFFGGGGGEGGGCGGGGB=>=F@FDFCGGGGGGGGFGGFGFFGFG>FDGG@FFGEFGECCEEGG>F7*CFGGFDDGFCC?FGGGGGGEFGF?EFG8FC8CGEGFF;;<E?FGGGGGGCED??FFGGGGGFGGGCEFGF>:E>5CFGFGFGDDDGG8=6?77*69@7EBD:>>>C54E465770/**77>*92>7C)))/6>EFE)48CC7+824::<F7(49(


