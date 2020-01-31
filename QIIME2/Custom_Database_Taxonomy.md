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

*pmoA* Accession number list: 191219_NCBI_pmoA_sequence_Accs.seq -- one # per line. File is saved in Reference_databse/NCBI_pmoA.

Navigate into your directory containing the accession list:
```
cd mnt/c/Users/xlibb/Desktop/QIIME2/Reference_database/NCBI_pmoA
```
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

## nosZ
* search nucleotide database for “nosZ” - select “bacteria” and the Genbank INSDC sequences - and download all accession numbers. (~29k records so can process all at once)

*nosZ*	Accession number list: 200103_NCBI_nosZ_bact_INSDC_acc.seq -- one # per line. Saved in Reference_database/NCBI_nosZ.

Navigate into your directory containing the accession list:
```
cd mnt/c/Users/xlibb/Desktop/QIIME2/Reference_database/NCBI_nosZ
```
Make a file (in notepad or in some commandline function) called “Entrez_fetch_nosZ_Tax_from_Accs_script.sh” containing:
```
for next in $(cat 200103_NCBI_nosZ_bact_INSDC_acc.seq); do LINEAGE=$(efetch -db nucleotide -id $next -format gbc | xtract -insd INSDSeq_taxonomy); echo -e "$next\t$LINEAGE"; done
```
Then run: 
```
bash Entrez_fetch_nosZ_Tax_from_Accs_script.sh > nosZ_Taxonomy_output.txt
```
Again, you need to have a functioning internet connection, and you can ignore errors until the command finishes (~7 hours).

# Formatting taxonomy file for import into QIIME2
Of course, the file created by this script is not in the right format for QIIME2. To import a taxonomy file into QIIME2, you need a tab-delimited file with accession number and semi-colon separated strings. 

The {gene}_Taxonomy_output.txt file we just made from the E-Utils download has 2 columns of accession numbers, so first we need to delete the second column. We also need to delete the spaces following each semi-colon, and make “environmental samples” be “environmental_samples”. 

## *pmoA*
First navigate to the directory containing the taxonomy output file:
```
cd mnt/c/Users/xlibb/Desktop/QIIME2/Reference_database/NCBI_pmoA
```
To remove the 2nd column of duplicate accession #s:
```
awk '{$2=""; print $0}' 191219_pmoA_Taxonomy_output.txt >191219_pmoA_Taxonomy_singleacc.txt
```
Next we need to: 
* Remove spaces after semicolons 
* turn “environmental samples” to “environmental_samples” & other space-separated lowercase word strings
* replace remaining spaces with just a single tab
* delete lines that have accession numbers followed by nothing
* delete lines that have Eukaryota
```
perl -pe 's/(;\s)/;/g;s/([a-z]) /\1_/g;s/ ([a-z|C-Z]{3,})/_\1/g;s/ {1,3}/\t/g;s/\t\n/\tUnknown\n/g;s/^.+Eukaryota.+\n//g' < 191219_pmoA_Taxonomy_singleacc.txt > 191219_pmoA_Taxonomy_clean.txt
```
At this point you can import into QIIME2 without errors, but if you want to be able to do anything by taxonomic level, you need to make sure there are the same number of taxonomic levels IN EACH LINE. 

There is surely a way to to count the number of ';' in each line and then add an appropriate # of ';' so that all lines have the same number of ';', but this eluded us (suggestions welcome). 
A method that works is to first find # of ';' in each line: 
```
sed 's/[^;]//g' 191219_nosZ_Taxonomy_clean.txt | awk '{ print length }' 
```
And scroll through the output and find max # (again, there's probably a better way to do this. Suggestions welcome). In this case, the maximum number of ';' appears to be **7**.

So, to add ';'s so there are 7 in each line following at tab (\t) & have “environmental samples” be the very last “taxonomic designation” (if present - otherwise it shows up as a genus or family or class or whatever): 
```
perl -pe 's/(\t(([^;]+;){7}[a-zA-Z0-9_]+))\n/\1\n/g;s/(\t(([^;]+;){6}[a-zA-Z0-9_]+))\n/\1;\n/g;s/(\t(([^;]+;){5}[a-zA-Z0-9_]+))\n/\1;;\n/g;s/(\t(([^;]+;){4}[a-zA-Z0-9_]+))\n/\1;;;\n/g;s/(\t(([^;]+;){3}[a-zA-Z0-9_]+))\n/\1;;;;\n/g;s/(\t(([^;]+;){2}[a-zA-Z0-9_]+))\n/\1;;;;;\n/g;s/(\t(([^;]+;){1}[a-zA-Z0-9_]+))\n/\1;;;;;;\n/g;s/(\t(([^;]+;){0}[a-zA-Z0-9_]+))\n/\1;;;;;;;\n/g;s/(environmental_samples)([;]+)/\2\1/g;s/ //g' 191219_pmoA_Taxonomy_clean.txt > 191219_pmoA_Taxonomy_fullstrings.txt
```
File is now ready for import into QIIME2!
## *nosZ*
First navigate to the directory containing the taxonomy output file:
```
cd mnt/c/Users/xlibb/Desktop/QIIME2/Reference_database/NCBI_nosZ
```
To remove the 2nd column of duplicate accession #s:
```
awk '{$2=""; print $0}' 200103_nosZ_Taxonomy_output.txt >200103_nosZ_Taxonomy_singleacc.txt
```
Remove spaces after semicolons, turn “environmental samples” to “environmental_samples” & other space-separated lowercase word strings, replace remaining spaces with just a single tab, delete lines that have accession numbers followed by nothing, delete lines that have Eukaryota (there shouldn't be any if you donwloaded only 'Bacteria' from NCBI, but just in case):
```
perl -pe 's/(;\s)/;/g;s/([a-z]) /\1_/g;s/ ([a-z|C-Z]{3,})/_\1/g;s/ {1,3}/\t/g;s/\t\n/\tUnknown\n/g;s/^.+Eukaryota.+\n//g' < 200103_nosZ_Taxonomy_singleacc.txt >200103_nosZ_Taxonomy_clean.txt
```
Find # of ';' in each line: 
```
sed 's/[^;]//g' 191219_nosZ_Taxonomy_clean.txt | awk '{ print length }' 
```
And scroll through the output and find max # , which appears to be **7**.

Add ;s so there are 7 in each line & have “environmental samples” be the very last “taxonomic designation” (if present): 
```
perl -pe 's/(\t(([^;]+;){7}[a-zA-Z0-9_]+))\n/\1\n/g;s/(\t(([^;]+;){6}[a-zA-Z0-9_]+))\n/\1;\n/g;s/(\t(([^;]+;){5}[a-zA-Z0-9_]+))\n/\1;;\n/g;s/(\t(([^;]+;){4}[a-zA-Z0-9_]+))\n/\1;;;\n/g;s/(\t(([^;]+;){3}[a-zA-Z0-9_]+))\n/\1;;;;\n/g;s/(\t(([^;]+;){2}[a-zA-Z0-9_]+))\n/\1;;;;;\n/g;s/(\t(([^;]+;){1}[a-zA-Z0-9_]+))\n/\1;;;;;;\n/g;s/(\t(([^;]+;){0}[a-zA-Z0-9_]+))\n/\1;;;;;;;\n/g;s/(environmental_samples)([;]+)/\2\1/g;s/ //g' 200103_nosZ_Taxonomy_clean.txt > 200103_nosZ_Taxonomy_fullstrings.txt
```
File is now ready for import into QIIME2!
# Finally import properly formatted taxonomy file into QIIME2:
Navigate to the "home" folder that contains your *nosZ* and *pmoA* folders
```
cd mnt/c/Users/xlibb/Desktop/QIIME2
```
And import the taxonomy files to the appropriate folders for each gene:
```
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path Reference_database/NCBI_pmoA/191219_pmoA_Taxonomy_fullstrings.txt \
  --output-path pmoA/ref-NCBI_pmoA-taxonomy.qza
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path Reference_database/NCBI_nosZ/200103_nosZ_Taxonomy_fullstrings.txt \
  --output-path nosZ/ref-NCBI_nosZ-taxonomy.qza
```
Voila! Taxonomy file created to complement your custom reference database. Proceed with your analyses!
