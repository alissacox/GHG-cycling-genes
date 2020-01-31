Here are the steps to create the taxonomy file corresponding to a [custom database](https://github.com/alissacox/GHG-cycling-genes/blob/master/QIIME2/Custom_Database_Creation.md) you created for use in QIIME2.

# Combining Accession numbers matching custom reference databse with taxonomy strings from NCBI databases

# Need to install Entrez Direct commandline tools
* Install Entrez for linux commandline: [see NCBI instructions](https://www.ncbi.nlm.nih.gov/books/NBK179288/ https://dataguide.nlm.nih.gov/edirect/install.html) 
  * on Windows 10 running Ubuntu as the Linux subsystem for Windows we had some problems during/post-installation... Most of the error messages were related to perl @INC folder missing various SSL and HTML modules, so had to install the follwing linux packages via the commandine: 
```
sudo apt-get install libio-socket-ssl-perl 
sudo apt-get install libhtml-html5-entities-perl
sudo apt-get install libhtml-parser-perl
sudo apt-get install libwww-perl
sudo apt-get install libxml-simple-perl
```
Ultimately, we had to re-install E-direct after installing those packages to get everything to work... This was not a problem when installing the E-utils in a native Ubuntu distribution!

Now that you have E-Utils direct installed...
# Download list of accession numbers from NCBI nucleotide search matching your [custom database](https://github.com/alissacox/GHG-cycling-genes/blob/master/QIIME2/Custom_Database_Creation.md): 
* Search for whatever you need on the NCBI nucleotide database
  * Send To > Complete Record > File > Accession list --- this creates a file with each accession # on its own line
# Using Eutils command line - use Acc list to fetch taxonomy strings
●	According to the [EDirect Cookbook](https://ncbi-hackathons.github.io/EDirectCookbook/), you can’t query more than 50,000 hits with E-Utils. Luckily we have fewer than 50K *pmoA* and *nosZ* sequences in our custom taxonomy databases
  * So if your custom database has more, you'll need to split up #  of accession numbers into 2 files and then paste outputs back together (should be reasonably easy with commandline unix text editing?)

## *pmoA*
* searched nucleotide database for “pmoa” and download all accession numbers (34k records - can process all at once) 

*pmoA* Accession number list: 191219_NCBI_pmoA_sequence_Accs.seq -- one # per line

Make a file (in notepad or in some commandline function) called “Entrez_fetch_pmoA_Tax_from_Accs_script.sh” containing:
```
for next in $(cat 191219_NCBI_pmoA_sequence_Accs.seq); do LINEAGE=$(efetch -db nucleotide -id $next -format gbc | xtract -insd INSDSeq_taxonomy); echo -e "$next\t$LINEAGE"; done
```
Then run the follwing to have E-Utils use each accession number to search the taxonomy string from the database and make a file with the accession # and the taxonomy string: 
```
bash Entrez_fetch_pmoA_Tax_from_Accs_script.sh > pmoA_Taxonomy_output.txt
```
Note: you must have an internet connection for this to work.

IGNORE the errors that pop up if the cursor is still “thinking” - the errors refer to those accession #s just are “empty” with no organism/lineage string. No big deal. This step took ~8 hrs to run total.

## nosZ  - you must have an internet connection for this to work
* search nucleotide database for “nosZ” - select “bacteria” and the Genbank INSDC sequences - and download all accession numbers. (~29k records so can process all at once)

*nosZ*	Accession number list: 200103_NCBI_nosZ_bact_INSDC_acc.seq -- one # per line

Make a file (in notepad or in some commandline function) called “Entrez_fetch_nosZ_Tax_from_Accs_script.sh” containing:
```
for next in $(cat 1191221_NCBI_nosZ_sequence_accs_.seq); do LINEAGE=$(efetch -db nucleotide -id $next -format gbc | xtract -insd INSDSeq_taxonomy); echo -e "$next\t$LINEAGE"; done
```
Then run: 
```
bash Entrez_fetch_nosZ_Tax_from_Accs_script.sh > nosZ_Taxonomy_output.txt
```
Again, you need to have a functioning internet connection, and you can ignore errors until the command finishes (~7 hours).
