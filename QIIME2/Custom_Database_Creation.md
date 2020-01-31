Here are the steps to create the "custom" reference database for *pmoA* and *nosZ* to assign taxonomy in the [QIIME2 pipeline](https://github.com/alissacox/GHG-cycling-genes/blob/master/QIIME2/Analysis_pipeline.md).

# Constructing reference database
* based on [this forum post](https://forum.qiime2.org/t/creating-a-custom-reference-database/3488) but with lots of detail that was missing
## a note on reference database selection...
* There are lots of options for functional gene reference databases. The [FunGene database]() is awesome, but downloading lots of their sequences at once is troublesome from the web interface, AND because their sequences are based on protein accession numbers, there are lots of duplicate gene accession numbers (at least for *pmoA* and *nosZ*) which make them impossible to import into QIIME2. There is probably some way around this (suggestions welcome).
* We ended up using the [NCBI nuleotide database](https://www.ncbi.nlm.nih.gov/nucleotide/) to download our "reference sequences" because it's easier to download all sequences matching a search term from the web interface at once. There are surely more elegant solutions to this (suggestions welcome).

# Download reference sequences (FASTA format) from [NCBI nucleotide database](https://www.ncbi.nlm.nih.gov/nucleotide/)
Can easily download fasta files from NCBI nucleotide search for all “hits”:
* Search nucleotide databade for gene you want (pmoA, nosZ)
  * For *nosZ* restrict hits to 'Bacteria' and the 'INSDC Genbank' database using filters on left or the database will misbehave later from giant Eukaryotic sequences matching *nosZ* hits ... no need to narrow down selections for *pmoA*, but you can if you want to...
* Once you have the search results you want:
  * Send To > Complete Record > File > FASTA -- don’t include GI numbers

To examine structure of file (which will be huge and a text editor may not allow you to open it) using the Linux commandline: 
* od prints non-printing characters (useful in identifying how lines are split/terminated)
* sed… command can print certain lines in text doc (to examine)
```
head -20 [filename]
od -N300 -c [filename]
sed -n '[startline#],[endline#]p;[endline#+1]q' [filename]
	Example:show lines 44-65: sed -n '44,65p;66q' file.txt
``` 
NCBI FASTA format = 1 line of descriptor starting with a “>”, then various lines of 70 characters of the sequence.
* Problem: .fasta file from NCBI  is not formatted correctly for QIIME2. QIIME2 requires that there is 1 line for each sequence with identifier INFO and 1 line **_only_** for gene sequence. So, we need to do some text editing to "clean up" the NCBI Fasta file so QIIME will be able to read it in... There are lots and lots and lots of ways to do this. The way that worked best for us is using perl to edit the FASTA text file. 

## *pmoA* NCBI Fasta file
Navigate to the NCBI FASTA file you downloaded - in this case we stuck our file ("191220_NCBI_pmoA_sequences.fasta") in a folder called "NCBI_pmoA" that is inside a "Reference_database" folder. To navigate there on the commandline (which may be different on your system):
``` 
cd mnt/c/Users/xlibb/Desktop/QIIME2/Reference_database/NCBI_pmoA
``` 
First, we need to remove linebreaks (\n) after the 70 characters of ACTG in each line of sequences. Then we need to make sure there are no blank lines in the file:
```
perl -pe 's/([ACGT]{70})\n/\1/g' < 191220_NCBI_pmoA_sequences.fasta > 191220_pmoA_NCBI_single_line_sequences.fasta
perl -pe 's/^\s$//g' < 191220_pmoA_NCBI_single_line_sequences.fasta > 191220_clean_pmoA_NCBI_sequences.fasta
```
The "191220_clean_pmoA_NCBI_sequences.fasta" file can be read into QIIME2 - see next step below.  You could use some of the file examination commands above to verify that the format is 1 line of sequence ID info followed by 1 line of Sequence (ACTG)....

## *nosZ* NCBI Fasta file -- downloaded ONLY bacterial nosZ sequences from INSDC:
Navigate to the NCBI FASTA file you downloaded - in this case we stuck our file ("200102_NCBI_nosZ_bact_INSDC_sequences.fasta") in a folder called "NCBI_nosZ" that is inside a "Reference_database" folder. To navigate there on the commandline (which may be different on your system):
``` 
cd mnt/c/Users/xlibb/Desktop/QIIME2/Reference_database/NCBI_nosZ
``` 
To remove linebreaks between sequences and remove blank lines:
``` 
perl -pe 's/([ACGT]{70})\n/\1/g' < 200102_NCBI_nosZ_bact_INSDC_sequences.fasta > 200102_nosZ_NCBI_single_line_sequences.fasta
perl -pe 's/^\s$//g' < 200102_nosZ_NCBI_single_line_sequences.fasta > 200102_clean_nosZ_NCBI_sequences.fasta
``` 
The "200102_clean_nosZ_NCBI_sequences.fasta" file is ready for import into QIIME2. You could use some of the file examination commands above to verify this....
# Import your “Clean” FASTA files (nicely formatted to consist of {(1 header row + 1 sequence line) times # sequences} into QIIME2
activate QIIME2 if you haven't already:
``` 
conda activate qiime2-2019.10
``` 
## pmoA
First navigate to your *pmoA* folder:
``` 
cd mnt/c/Users/xlibb/Desktop/QIIME2/pmoA
``` 
Now import your 'clean' (nicely formatted) Fasta file of *pmoA* reference sequences:
``` 
qiime tools import \
  --input-path Reference_database/NCBI_pmoA/191220_clean_pmoA_NCBI_sequences.fasta \
  --input-format DNAFASTAFormat \
  --output-path pmoA/ref_NCBI_pmoA_seqs.qza \
  --type 'FeatureData[Sequence]'
#To check if this is formatted correctly, run:
qiime tools validate ref_NCBI_pmoA_seqs.qza
#To examine file:
qiime feature-table tabulate-seqs \
  --i-data ref_NCBI_pmoA_seqs.qza \
  --o-visualization ref_NCBI_pmoA_seqs.qzv
``` 	
## nosZ
First navigate to your *nosZ* folder:
``` 
cd mnt/c/Users/xlibb/Desktop/QIIME2/nosZ
``` 
Import your fasta file of *nosZ* reference sequences:
``` 
qiime tools import \
  --input-path Reference_database/NCBI_nosZ/200102_clean_nosZ_NCBI_sequences.fasta\
  --input-format DNAFASTAFormat \
  --output-path nosZ/ref_NCBI_nosZ_seqs.qza \
  --type 'FeatureData[Sequence]'  
#To check if this is formatted correctly, run:
qiime tools validate nosZ/ref_NCBI_nosZ_seqs.qza
``` 
Now you have two reference databse files! However, these do not seem to be well-suited to the Deblur [denoise-other](https://docs.qiime2.org/2019.10/plugins/available/deblur/denoise-other/) alignment and filtering protocol. When we tried to use these databases that way, the script always failed (after running FOREVER). It wasn't obvious to us why that happened. But the databases do work to assign taxonomy (see [pipeline](https://github.com/alissacox/GHG-cycling-genes/blob/master/QIIME2/Analysis_pipeline.md)) as long as you have [fetched taxonomy strings](https://github.com/alissacox/GHG-cycling-genes/blob/master/QIIME2/Custom_Database_Taxonomy.md) to match the accession numbers at the beginning of each line of the FASTA files you just formatted and imported into QIIME2.
