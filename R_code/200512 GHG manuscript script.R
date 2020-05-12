##Code for analyzing data from GHG sampling / soil cores - Acushnet in JUly 2019
##Includes code for calculating alpha & beta diversity from .tsv files from QIIME using phyloseq

#### Data Analysis setup ##############################################################################################################################
#Set Working Directory:
setwd("G:/My Drive/PhD/R_work/GHG manuscript")


#Load Packages:
library(psych)
library(tidyverse)
library(ggiraphExtra) #lets you do fancy graphs with interactive linear models
library(viridis)
library(ggpubr)
library(phyloseq) #to install: BiocManager::install("phyloseq")
library(vegan) # included in phyloseq
library(ggdendro) # for making dendrogram
library(ggrepel)
library(SpiecEasi) #to install: library(devtools) & install_github("zdk123/SpiecEasi")
library(extrafont) # needed to change fonts in figures
#   Need to install Ghostscript (https://www.ghostscript.com/download/gsdnld.html) for this to work on Windows Machine
#   font_install('fontcm') to install Computer Modern fonts
#   Also need to run font_import() after first installation of 'extrafont' package, 
#       or anytime new fonts installed on computer... Also need to run loadfonts(device = "win") after that


sessionInfo()
#R version 4.0.0 (2020-04-24)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows >= 8 x64 (build 9200)
#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggpubr_0.2.5       magrittr_1.5       ggrepel_0.8.2      ggdendro_0.1-20    phyloseq_1.31.0    viridis_0.5.1      viridisLite_0.3.0 
#[8] ggiraphExtra_0.2.9 forcats_0.5.0      stringr_1.4.0      dplyr_0.8.5        purrr_0.3.4        readr_1.3.1        tidyr_1.0.2       
#[15] tibble_3.0.0       ggplot2_3.3.0      tidyverse_1.3.0    psych_1.9.12.31    SpiecEasi_1.0.7   


#Read in files:
BD_Veg_flux_data <- read.csv("raw_data/190913 GHG Soil Summary.csv")

nosZ_taxonomy_table <- read.csv("raw_data/200107 nosZ taxonomy table.csv") #samples v taxonomy...cleaned up taxonomy from original QIIME2 files
pmoA_taxonomy_table <- read.csv("raw_data/200122 pmoA taxonomy table.csv") #samples v taxonomy...cleaned up taxonomy from original QIIME2 files -- including replacing "streptomyces" with "Bacteria"
nosZ_ASVs_in <- read.csv("raw_data/200107 nosZ ASV table.csv")
pmoA_ASVs_in <- read.csv("raw_data/200107 pmoA ASV table.csv") #filtered from from original dada2 outputs files to only keep samples with >=1 feature
nosZ_taxonomy_in <- read.csv("raw_data/200421_nosZ_taxonomy_strings_complete.csv") #exported from original QIIME2 files (dada2_nosZ_taxonomy.qzv - extract .tsv version), separated taxonomy strings into separate columns
pmoA_taxonomy_in <- read.csv("raw_data/200421_pmoA_taxonomy_strings_complete.csv") #exported from original QIIME2 files (dada2_nosZ_taxonomy.qzv - extract .tsv version), separated taxonomy strings into separate columns

nosZ_rtree <- read_tree("raw_data/nosZ_rooted_tree.nwk")
pmoA_rtree <- read_tree("raw_data/pmoA_rooted_tree.nwk")
gene_metadata <- read.delim("raw_data/200113_AHC_sequencing_sample_GHG_metadata.txt", sep = "\t", row.names = 1)


#change ASV & taxonomy tables to be phyloseq-compatible (numeric matix)
nosZ_ASVs <- nosZ_ASVs_in %>%
  select(-OTU_ID_nosZ) %>%
  as.matrix
rownames(nosZ_ASVs) <- nosZ_ASVs_in$OTU_ID_nosZ
nosZ_taxonomy <- nosZ_taxonomy_in %>%
  select(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Source")) %>%
  as.matrix
rownames(nosZ_taxonomy) <- nosZ_taxonomy_in$Feature.ID

pmoA_ASVs <- pmoA_ASVs_in %>%
  select(-OTU_ID_pmoA) %>%
  as.matrix
rownames(pmoA_ASVs) <- pmoA_ASVs_in$OTU_ID_pmoA
pmoA_taxonomy <- pmoA_taxonomy_in %>%
  select(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Source")) %>%
  as.matrix
rownames(pmoA_taxonomy) <- pmoA_taxonomy_in$Feature.ID

# for dendrogram / non-phyloseq analyses
nosZ_ASVs_t <- t(nosZ_ASVs)
pmoA_ASVs_t <- t(pmoA_ASVs)

#use phyloseq() to combine text files (otu tables, taxonomy, sample_data etc.) into phyloseq object... https://joey711.github.io/phyloseq/import-data.html
  #structure: # phyloseq(otu_table(GP), phy_tree(GP), tax_table(GP), sample_data(GP))
nosZ_phylo_raw <- phyloseq(otu_table(nosZ_ASVs, taxa_are_rows = T), phy_tree(nosZ_rtree), tax_table(nosZ_taxonomy), sample_data(gene_metadata))
pmoA_phylo_raw <- phyloseq(otu_table(pmoA_ASVs, taxa_are_rows = T), phy_tree(pmoA_rtree), tax_table(pmoA_taxonomy), sample_data(gene_metadata))

sample_variables(pmoA_phylo_raw) #list metadata column names

#Example to remove samples with low read counts <- prune_samples(sampleSums(GP.chl)>=20, GP.chl)
sample_sums(pmoA_phylo_raw) #-- 5 samples have single digits (2-5)... next highest have 12 counts. Max = 550
summary(sample_sums(pmoA_phylo_raw)) #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.00   18.00   44.00   76.49   77.00  550.00 
boxplot(sample_sums(pmoA_phylo_raw))
taxa_sums(pmoA_phylo_raw) # -- lots of ASVs (~half) represented in single digits reads. one taxon has 831 reads
summary(taxa_sums(pmoA_phylo_raw)) # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    4.00    7.00   19.48   15.00  831.00
boxplot(taxa_sums(pmoA_phylo_raw))

sample_sums(nosZ_phylo_raw) #-- AHC94&95 have 2956 & 3577 reads... next highest has 4591 and they increase "normally" from there...
summary(sample_sums(nosZ_phylo_raw))  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2956   10178   13886   13479   16647   31163 
boxplot(sample_sums(nosZ_phylo_raw))
taxa_sums(nosZ_phylo_raw) # -- lots of ASVs (~1700) represented in single digits reads. one taxon has 21k reads
summary(taxa_sums(nosZ_phylo_raw)) # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.0     9.0    18.0   117.4    37.0 21186.0 
boxplot(taxa_sums(nosZ_phylo_raw))

#Example to remove singlteon taxa <- prune_taxa(taxa_sums(notree_musk_surface_pruned) > 1, notree_musk_surface_pruned) 
pmoA_phylo_raw_rm3s <- prune_taxa(taxa_sums(pmoA_phylo_raw) > 3, pmoA_phylo_raw) #removed taxa with 3 or fewer occurrences => dropped from 161 to 130 taxa
#nosZ_phylo_raw_rm3s <- prune_taxa(taxa_sums(nosZ_phylo_raw) > 3, nosZ_phylo_raw) #removed taxa with 3 or fewer occurrences => dropped from 5739 to 5236 taxa

#pmoA_phylo_raw_rm10s <- prune_taxa(taxa_sums(pmoA_phylo_raw) > 10, pmoA_phylo_raw) #removed taxa with 10 or fewer occurrences => dropped from 161 to 53 taxa
nosZ_phylo_raw_rm10s <- prune_taxa(taxa_sums(nosZ_phylo_raw) > 10, nosZ_phylo_raw) #removed taxa with 5 or fewer occurrences => dropped from 5739 to 4051 taxa


#Quality Filtering Decisions:
#pmoA: KEEP all taxa (not just with >3 reads ("pmoA_phylo_raw_rm3s")); rarify @ 10 ASVs (drops ~5 samples)
#nosZ: KEEP taxa with >10 reads ("nosZ_phylo_raw_rm10s"); rarify @4510 ASVs (drops 2 samples)

#RARIFY data:
pmoA_phylo_rare = rarefy_even_depth(pmoA_phylo_raw, rngseed = 16, sample.size = 10) # removes: AHC13AHC24AHC27AHC31AHC37 & 69OTUs
taxa_sums(pmoA_phylo_rare) #use this as it more conservatively preserves ASVs...
#RENAME RARIFIED OBJECT as "pmoA_phylo" to use existing code below:
pmoA_phylo <- pmoA_phylo_rare 
nosZ_phylo_10s_rare = rarefy_even_depth(nosZ_phylo_raw_rm10s, rngseed = 16, sample.size = 4510) # removes: AHC94AHC95 & 93 OTUs
#RENAME RARIFIED OBJECT as "nosZ_phylo" to use existing code below:
nosZ_phylo <- nosZ_phylo_10s_rare 

#Relative abundance (by PCT) -- based on rarified data:
nosZ_phylo_pct <- transform_sample_counts(nosZ_phylo_10s_rare, function(x) 100 * x/sum(x))
pmoA_phylo_pct <- transform_sample_counts(pmoA_phylo_rare, function(x) 100 * x/sum(x))

#do some basic ggplot themesettings:
theme_set(theme_classic(base_size = 15) + theme(text = element_text(family = "Lato")))

##################################################### START ANALYSIS #################################################################################

#### Soil Core Bulk Density & Vegetation Data ########################################################################################################
#descriptive stats (psych packa)
describeBy(BD_Veg_flux_data, group = BD_Veg_flux_data$Location)

#no significant differences among groups for BD, OM or Biomass:
kruskal.test(BD_Veg_flux_data$Core_BD_g_cm3, g = BD_Veg_flux_data$Location) #p = 0.051
pairwise.wilcox.test(x = BD_Veg_flux_data$Core_BD_g_cm3, g = BD_Veg_flux_data$Location,
                     p.adjust.method = "BH")
kruskal.test(BD_Veg_flux_data$Perc_OM, g = BD_Veg_flux_data$Location) #p = 0.837
pairwise.wilcox.test(x = BD_Veg_flux_data$Perc_OM, g = BD_Veg_flux_data$Location,
                     p.adjust.method = "BH")
kruskal.test(BD_Veg_flux_data$Biomass_g_cm2, g = BD_Veg_flux_data$Location) #p = 0.058
pairwise.wilcox.test(x = BD_Veg_flux_data$Biomass_g_cm2, g = BD_Veg_flux_data$Location,
                     p.adjust.method = "BH")

#boxplots of BD, OM, Veg Biomass
BD_pt <- ggplot(data=BD_Veg_flux_data, aes(x=Location, y=Core_BD_g_cm3)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15) +
  labs(y=expression("Bulk density " (g/cm^3)))
BD_pt

OM_pt <- ggplot(data=BD_Veg_flux_data, aes(x=Location, y=Perc_OM)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15)+
  labs(y="Organic matter (%)")
OM_pt

Veg_pt <- ggplot(data=BD_Veg_flux_data, aes(x=Location, y=Biomass_g_cm2)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15)+
  labs(y=expression("Above-ground biomass "(g/cm^2)))
Veg_pt
BD_OM_Veg_pt <- ggarrange(BD_pt, OM_pt, Veg_pt, ncol = 1, align = "v")
BD_OM_Veg_pt

ggsave("200109 BD OM Veg Boxes.png", plot = BD_OM_Veg_box, device = "png", path = NULL,
       scale = 1, width = 4, height = 7.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200109 BD OM Veg Boxes.pdf", plot = BD_OM_Veg_box, device = cairo_pdf, 
       scale = 1, width = 6.5, height = 7.5, units = "in",
       dpi = 600)


#Correlations among soil properties:
BioM_v_BD_plot <- ggplot(data=BD_Veg_flux_data, aes(x=Core_BD_g_cm3, y=Biomass_g_cm2)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic(base_size = 15) +
  ylab(expression("Above-ground biomass "(g/cm^2))) + xlab(expression("Bulk Density "(g/cm^3)))
BioM_v_BD_plot # SIGNIFICANT p = 0.006

BD_v_OM_lm <- lm(Perc_OM ~ Core_BD_g_cm3, data = BD_Veg_flux_data)
summary(BD_v_OM_lm) # NOT significant 
BioM_v_BD_lm <- lm(Core_BD_g_cm3 ~ Biomass_g_cm2, data = BD_Veg_flux_data)
summary(BioM_v_BD_lm) ## SIGNIFICANT p = 0.006
OM_v_BioM_lm <- lm(Biomass_g_cm2 ~ Perc_OM, data = BD_Veg_flux_data)
summary(OM_v_BioM_lm) # NOT significant


#### Surface GHG flux data ##########################################################################################################################
#no significant differences among groups for CH4, CO2 or N2O: -- Not really a valid test sice this is pseudoreplication....
kruskal.test(BD_Veg_flux_data$CH4_Flux_umol_m.2_h, g = BD_Veg_flux_data$Location) #p = 0.329
kruskal.test(BD_Veg_flux_data$CO2_Flux_umol_m.2_h, g = BD_Veg_flux_data$Location) #p = 0.066
kruskal.test(BD_Veg_flux_data$N2O_Flux_umol_m.2_h, g = BD_Veg_flux_data$Location) #p = 0.051

#boxplots of GHGs -- only 3 values so doing points instead
CH4_pt <- ggplot(data=BD_Veg_flux_data, aes(x=Location, y=CH4_Flux_umol_m.2_h)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15) +
  theme(text = element_text(family = "Lato")) +
  labs(y=expression(paste(CH[4]," flux (", mu,mol/m^2/h,")"))) +
  geom_hline(yintercept = 0, color = 'grey', linetype = 2)
CH4_pt

CO2_pt <- ggplot(data=BD_Veg_flux_data, aes(x=Location, y=CO2_Flux_umol_m.2_h)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15)+
  theme(text = element_text(family = "Lato")) +
  labs(y=expression(paste(CO[2]," flux (", mu,mol/m^2/h,")"))) +
  ggtitle("Flux")
CO2_pt

N2O_pt <- ggplot(data=BD_Veg_flux_data, aes(x=Location, y=N2O_Flux_umol_m.2_h)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15) +
  theme(text = element_text(family = "Lato")) +
  labs(y=expression(paste(N[2]*O," flux (", mu,mol/m^2/h,")")))
N2O_pt
GHGs_pt <- ggarrange(CO2_pt, CH4_pt, N2O_pt, ncol = 1, align = "v")
GHGs_pt

ggsave("200121 GHG Pts.png", plot = GHGs_pt, device = "png", path = NULL,
       scale = 1, width = 4, height = 7.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200512 GHG Pts.pdf", plot = GHGs_pt, device = cairo_pdf, 
       scale = 1, width = 6.5, height = 7.5, units = "in",
       dpi = 600)

##Correlations of GHG flux with soil properties:
#Methane
CH4_flux_v_BD_plot <- ggplot(data = BD_Veg_flux_data, aes(x = Core_BD_g_cm3, y = CH4_Flux_umol_m.2_h)) +
  geom_point() +
  theme_classic(base_size = 15)+
  geom_smooth(method = 'lm')
CH4_flux_v_BD_plot

CH4_flux_v_OM_plot <- ggplot(data = BD_Veg_flux_data, aes(x = Perc_OM, y = CH4_Flux_umol_m.2_h)) +
  geom_point() +
  theme_classic(base_size = 15)+
  geom_smooth(method = 'lm')
CH4_flux_v_OM_plot #clear inverse relationship -- not quite significant: p=0.069

CH4_flux_v_BioM_plot <- ggplot(data = BD_Veg_flux_data, aes(x = Biomass_g_cm2, y = CH4_Flux_umol_m.2_h)) +
  geom_point() +
  theme_classic(base_size = 15)+
  geom_smooth(method = 'lm')
CH4_flux_v_BioM_plot #inverse relationship

CH4_flux_v_BD_lm <- lm(CH4_Flux_umol_m.2_h ~ Core_BD_g_cm3, data = BD_Veg_flux_data)
summary(CH4_flux_v_BD_lm) #not sig
CH4_flux_v_OM_lm <- lm(CH4_Flux_umol_m.2_h ~ Perc_OM, data = BD_Veg_flux_data)
summary(CH4_flux_v_OM_lm) #not sig
CH4_flux_v_BioM_lm <- lm(CH4_Flux_umol_m.2_h ~ Biomass_g_cm2, data = BD_Veg_flux_data)
summary(CH4_flux_v_BioM_lm) #not sig
CH4_flux_v_BioM_BD_lm <- lm(CH4_Flux_umol_m.2_h ~ (Biomass_g_cm2+Core_BD_g_cm3), data = BD_Veg_flux_data)
summary(CH4_flux_v_BioM_BD_lm) #SIGNIFICANT! p = 0.032 
CH4_flux_v_BioM_BD_lm2 <- lm(CH4_Flux_umol_m.2_h ~ (Biomass_g_cm2*Core_BD_g_cm3), data = BD_Veg_flux_data)
summary(CH4_flux_v_BioM_BD_lm2) #NOT SIGNIFICANT! p = 0.08

ggPredict(CH4_flux_v_BioM_BD_lm,interactive=TRUE)


#Nitrous Oxide
N2O_flux_v_BD_plot <- ggplot(data = BD_Veg_flux_data, aes(x = Core_BD_g_cm3, y = N2O_Flux_umol_m.2_h)) +
  geom_point() +
  theme_classic(base_size = 15)+
  geom_smooth(method = 'lm')
N2O_flux_v_BD_plot #clear inverse relationship

N2O_flux_v_OM_plot <- ggplot(data = BD_Veg_flux_data, aes(x = Perc_OM, y = N2O_Flux_umol_m.2_h)) +
  geom_point() +
  theme_classic(base_size = 15)+
  geom_smooth(method = 'lm')
N2O_flux_v_OM_plot #clear positive relationship -- not significant

N2O_flux_v_BioM_plot <- ggplot(data = BD_Veg_flux_data, aes(x = Biomass_g_cm2, y = N2O_Flux_umol_m.2_h)) +
  geom_point() +
  theme_classic(base_size = 15)+
  geom_smooth(method = 'lm')
N2O_flux_v_BioM_plot #v strong positive relationship -- significant (p=0.006) 

N2O_flux_v_BD_lm <- lm(N2O_Flux_umol_m.2_h ~ Core_BD_g_cm3, data = BD_Veg_flux_data)
summary(N2O_flux_v_BD_lm) #not sig p = 0.062
N2O_flux_v_OM_lm <- lm(N2O_Flux_umol_m.2_h ~ Perc_OM, data = BD_Veg_flux_data)
summary(N2O_flux_v_OM_lm) #not sig 
N2O_flux_v_BioM_lm <- lm(N2O_Flux_umol_m.2_h ~ Biomass_g_cm2, data = BD_Veg_flux_data)
summary(N2O_flux_v_BioM_lm) #SIGNIFICANT! p = 0.006 
N2O_flux_v_BioM_BD_lm <- lm(N2O_Flux_umol_m.2_h ~ (Biomass_g_cm2+Core_BD_g_cm3), data = BD_Veg_flux_data)
summary(N2O_flux_v_BioM_BD_lm) #SIGNIFICANT! p = 0.030
N2O_flux_v_BioM_BD_lm2 <- lm(N2O_Flux_umol_m.2_h ~ (Biomass_g_cm2*Core_BD_g_cm3), data = BD_Veg_flux_data)
summary(N2O_flux_v_BioM_BD_lm2) #SIGNIFICANT! p = 0.02

ggPredict(N2O_flux_v_BioM_BD_lm,interactive=TRUE)


#Carbon Dioxide
CO2_flux_v_BD_plot <- ggplot(data = BD_Veg_flux_data, aes(x = Core_BD_g_cm3, y = CO2_Flux_umol_m.2_h)) +
  geom_point() +
  theme_classic(base_size = 15)+
  geom_smooth(method = 'lm')
CO2_flux_v_BD_plot #clear inverse relationship - significant! p=0.03

CO2_flux_v_OM_plot <- ggplot(data = BD_Veg_flux_data, aes(x = Perc_OM, y = CO2_Flux_umol_m.2_h)) +
  geom_point() +
  theme_classic(base_size = 15)+
  geom_smooth(method = 'lm')
CO2_flux_v_OM_plot #weak positive relationship -- no significant

CO2_flux_v_BioM_plot <- ggplot(data = BD_Veg_flux_data, aes(x = Biomass_g_cm2, y = CO2_Flux_umol_m.2_h)) +
  geom_point() +
  theme_classic(base_size = 15)+
  geom_smooth(method = 'lm')
CO2_flux_v_BioM_plot #v strong positive relationship -- SIGNIFICANT! p = 0.003

CO2_flux_v_BD_lm <- lm(CO2_Flux_umol_m.2_h ~ Core_BD_g_cm3, data = BD_Veg_flux_data)
summary(CO2_flux_v_BD_lm) # significant p = 0.03
CO2_flux_v_OM_lm <- lm(CO2_Flux_umol_m.2_h ~ Perc_OM, data = BD_Veg_flux_data)
summary(CO2_flux_v_OM_lm) #not sig
CO2_flux_v_BioM_lm <- lm(CO2_Flux_umol_m.2_h ~ Biomass_g_cm2, data = BD_Veg_flux_data)
summary(CO2_flux_v_BioM_lm) #SIGNIFICANT! p = 0.003
CO2_flux_v_BioM_BD_lm <- lm(CO2_Flux_umol_m.2_h ~ (Biomass_g_cm2+Core_BD_g_cm3), data = BD_Veg_flux_data)
summary(CO2_flux_v_BioM_BD_lm) #SIGNIFICANT! p = 0.022
CO2_flux_v_BioM_BD_lm2 <- lm(CO2_Flux_umol_m.2_h ~ (Biomass_g_cm2*Core_BD_g_cm3), data = BD_Veg_flux_data)
summary(CO2_flux_v_BioM_BD_lm2) #SIGNIFICANT! p = 0.001

ggPredict(CO2_flux_v_BioM_BD_lm,interactive=TRUE)


#### Below-ground GHG conc data ######################################################################################################################
#no significant differences among groups for CH4, CO2 or N2O: -- Not really a valid test sice this is pseudoreplication....
kruskal.test(BD_Veg_flux_data$CH4_conc_6in_PPM, g = BD_Veg_flux_data$Location) #p = 0.1266
kruskal.test(BD_Veg_flux_data$C02_conc_6in_PPM, g = BD_Veg_flux_data$Location) #p = 0.1266
kruskal.test(BD_Veg_flux_data$N2O_conc_6in_PPM, g = BD_Veg_flux_data$Location) #p = 0.5127

#boxplots of GHGs -- use points instead
CH4_below_pt <- ggplot(data=subset(BD_Veg_flux_data, Location=='Control'|Location=='Experimental'), aes(x=Location, y=CH4_conc_6in_PPM)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15) +
  theme(text = element_text(family = "Lato")) +
  labs(y=expression(CH[4]~concentration~(PPM)))
CH4_below_pt

CO2_below_pt <- ggplot(data=subset(BD_Veg_flux_data, Location=='Control'|Location=='Experimental'), aes(x=Location, y=C02_conc_6in_PPM)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15) +
  theme(text = element_text(family = "Lato")) +
  labs(y=expression(CO[2]~concentration~(PPM))) +
  ggtitle("Concentration")
CO2_below_pt

N2O_below_pt <- ggplot(data=subset(BD_Veg_flux_data, Location=='Control'|Location=='Experimental'), aes(x=Location, y=N2O_conc_6in_PPM)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15) +
  theme(text = element_text(family = "Lato")) +
  labs(y=expression(N[2]*O~concentration~(PPM)))
N2O_below_pt
GHGs_below_pt <- ggarrange(CO2_below_pt, CH4_below_pt, N2O_below_pt, ncol = 1, align = "v")
GHGs_below_pt

ggsave("200107 GHG below Pts.png", plot = GHGs_pt, device = "png", path = NULL,
       scale = 1, width = 4, height = 7.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200107 GHG below Pts.pdf", plot = GHGs_pt, device = cairo_pdf, 
       scale = 1, width = 6.5, height = 7.5, units = "in",
       dpi = 600)

Flux_v_belowground_GHG <- ggarrange(CO2_pt, CO2_below_pt, CH4_pt, CH4_below_pt, 
                                    N2O_pt, N2O_below_pt, ncol = 2, nrow = 3, align = "v")
Flux_v_belowground_GHG ################################################################################################ Manuscript Fig 3

ggsave("200512 GHG above v below Pts.png", plot = Flux_v_belowground_GHG, device = "png", 
       scale = 1, width = 12, height = 8.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200512 GHG above v below Pts.pdf", plot = Flux_v_belowground_GHG, device = cairo_pdf, 
       scale = 1, width = 12, height = 8.5, units = "in",
       dpi = 600)

#### Rarifying #########################################################################################################


## RARIFY ##
    #pmoA: #Plot the rarefaction curves using vegan function rarecurve():of ps
rarecurve(t(otu_table(pmoA_phylo_raw)), ylab = "pmoA ASVs", label = FALSE)#it appears that all samples plateau

rarecurve(t(otu_table(pmoA_phylo_raw_rm3s)), ylab = "ASVs", label = FALSE, title = "pmoA") #it appears that all samples plateau

#Quality Filtering Decisions: pmoA: KEEP all taxa (not just with >3 reads ("pmoA_phylo_raw_rm3s") (makes 2 samples have 0 ASVs)); rarify @ 10 ASVs (drops 6 samples)
#pmoA_phylo_3s_rare = rarefy_even_depth(pmoA_phylo_raw_rm3s, rngseed = 16, sample.size = 10) # removes: AHC08AHC13AHC24AHC27AHC31 & 54OTUs
#taxa_sums(pmoA_phylo_3s_rare)
pmoA_phylo_rare = rarefy_even_depth(pmoA_phylo_raw, rngseed = 16, sample.size = 10) # removes: AHC13AHC24AHC27AHC31AHC37 & 69OTUs
taxa_sums(pmoA_phylo_rare) #use this as it more conservatively preserves ASVs...
#RENAME RARIFIED OBJECT as "pmoA_phylo" to use existing code below:
pmoA_phylo <- pmoA_phylo_rare 

    #nosZ:
rarecurve(t(otu_table(nosZ_phylo_raw)), ylab = "nosZ ASVs", label = FALSE) #it appears that all samples plateau
rarecurve(t(otu_table(nosZ_phylo_raw_rm10s)), ylab = "nosZ ASVs", label = FALSE) #it appears that all samples plateau
##Quality Filtering Decisions: nosZ: KEEP taxa with >10 reads ("nosZ_phylo_raw_rm10s"); rarify @4510 ASVs (drops 2 samples)
nosZ_phylo_10s_rare = rarefy_even_depth(nosZ_phylo_raw_rm10s, rngseed = 16, sample.size = 4510) # removes: AHC94AHC95 & 93 OTUs
#RENAME RARIFIED OBJECT as "nosZ_phylo" to use existing code below:
nosZ_phylo <- nosZ_phylo_10s_rare 

############### S T A T I S T I C A L    C O M M U N I T Y    A N A L Y S I S ####################

#### Alpha Diversity ####
#pmoA:
#Alpha Diversity:
pmoA_alphadiv <- estimate_richness(pmoA_phylo, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")) #Can't estimate fisher's...
pmoA_alpha_plots <- plot_richness(pmoA_phylo, measures = c("Observed", "Shannon", "Simpson"), x = "Sampling_ID", nrow = 3) + xlab("Sample") + ggtitle(expression(italic("pmoA")))
plot_richness(pmoA_phylo, measures = c("Observed", "Shannon", "Simpson"), x = "Depth") + scale_color_viridis(discrete = T) + 
  geom_boxplot(aes(group = Depth)) +
  theme(legend.position = "bottom")

#Richness statistics:
pairwise.wilcox.test(pmoA_alphadiv$Observed, sample_data(pmoA_phylo)$Depth) #not significant
pairwise.wilcox.test(pmoA_alphadiv$Shannon, sample_data(pmoA_phylo)$Depth) #not significant
pairwise.wilcox.test(pmoA_alphadiv$Simpson, sample_data(pmoA_phylo)$Depth) #not significant
pairwise.wilcox.test(pmoA_alphadiv$Observed, sample_data(pmoA_phylo)$Location) #not significant
pairwise.wilcox.test(pmoA_alphadiv$Shannon, sample_data(pmoA_phylo)$Location) #not significant
pairwise.wilcox.test(pmoA_alphadiv$Simpson, sample_data(pmoA_phylo)$Location) #not significant


#stacked Taxonomy bar plots
  #total abundance - based on rarification set to 10
#pmoA_tax_tot_abund <- plot_bar(pmoA_phylo, fill = "Genus", x = "Sampling_ID") + 
#  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + theme(legend.position = "top") +
#  labs(x = "Sample")
#pmoA_tax_tot_abund_depth <- plot_bar(pmoA_phylo, fill = "Genus", x = "Sampling_ID", color = "Genus") + facet_grid(~Depth, scales = "free") + 
#  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + theme(legend.position = "top") +
#  labs(x = "Sample")

  #relative abundance - based on rarified data
pmoA_rel_abund_tax_bar <- plot_bar(pmoA_phylo_pct, fill = "Genus", x = "Sampling_ID") +
  geom_bar(stat="identity", position="stack") +
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + 
  theme(legend.position = "top", legend.text = element_text(size = 9, family = "Lato")) +
  labs(x = "Sample", y = "Abundance (%)")
pmoA_rel_abund_tax_bar_depth <- plot_bar(pmoA_phylo_pct, fill = "Genus", x = "Sampling_ID") + facet_grid(~Depth, scales = "free") + 
  geom_bar(stat="identity", position="stack") +
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + 
  theme(legend.position = "top", legend.text = element_text(size = 9, family = "Lato")) +
  labs(x = "Sample", y = "Abundance (%)")
pmoA_rel_abund_tax_bar_depth ################################################################################################ Manuscript Fig 5
ggsave("200421 pmoA rel abund Tax Bar.png", plot = pmoA_rel_abund_tax_bar, device = "png", 
       scale = 1, width = 8.5, height = 6.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200512 pmoA rel abund Tax Bar by depth.png", plot = pmoA_rel_abund_tax_bar_depth, device = "png", 
       scale = 1, width = 8.5, height = 6.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("Fig 5 200512 pmoA rel abund Tax Bar by depth.pdf", plot = pmoA_rel_abund_tax_bar_depth, device = cairo_pdf, 
       scale = 1, width = 12, height = 8.5, units = "in",
       dpi = 600)


#nosZ:
#Alpha Diversity: -- based on rarified data
nosZ_alphadiv <- estimate_richness(nosZ_phylo)
nosZ_alpha_plots <- plot_richness(nosZ_phylo, measures = c("Observed", "Shannon", "Simpson"), x = "Sampling_ID", nrow = 3) + xlab("Sample") + ggtitle(expression(italic("nosZ")))
all_alpha_plots <- ggarrange(pmoA_alpha_plots, nosZ_alpha_plots, nrow = 2)
ggsave("200421 all alpha diversity plots.png", plot = all_alpha_plots, device = "png", 
       scale = 1, width = 7.5, height = 10, units = "in",
       dpi = 600, type = "cairo")

plot_richness(nosZ_phylo, measures = c("Observed", "Shannon", "Simpson"), x = "Sample_Type", 
              color = "Location") + scale_color_viridis(discrete = T) + 
  geom_boxplot(aes(group = Sample_Type)) +
  theme(legend.position = "bottom")

plot_richness(nosZ_phylo, measures = c("Observed", "Shannon", "Simpson"), x = "Depth", 
              color = "Location") + scale_color_viridis(discrete = T) + 
  geom_boxplot(aes(group = Depth)) +
  theme(legend.position = "bottom")

#Richness statistics: -- no significant differences
pairwise.wilcox.test(nosZ_alphadiv$Observed, sample_data(nosZ_phylo)$Depth) #not significant
pairwise.wilcox.test(nosZ_alphadiv$Shannon, sample_data(nosZ_phylo)$Depth) #not significant
pairwise.wilcox.test(nosZ_alphadiv$Simpson, sample_data(nosZ_phylo)$Depth) #not significant
pairwise.wilcox.test(nosZ_alphadiv$Observed, sample_data(nosZ_phylo)$Location) #not significant
pairwise.wilcox.test(nosZ_alphadiv$Shannon, sample_data(nosZ_phylo)$Location) #not significant
pairwise.wilcox.test(nosZ_alphadiv$Simpson, sample_data(nosZ_phylo)$Location) #not significant
pairwise.wilcox.test(nosZ_alphadiv$Observed, sample_data(nosZ_phylo)$Sample_Type) #not significant
pairwise.wilcox.test(nosZ_alphadiv$Shannon, sample_data(nosZ_phylo)$Sample_Type) #not significant
pairwise.wilcox.test(nosZ_alphadiv$Simpson, sample_data(nosZ_phylo)$Sample_Type) #not significant


#stacked bar plots
  #total abundance
#plot_bar(nosZ_phylo, fill = "Genus", x = "Sampling_ID") + facet_grid(~Sample_Type, scales = "free") + 
#  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 3)) + theme(legend.position = "top") +
#  labs(x = "Sample")

#relative abundance
nosZ_rel_abund_tax_bar <- plot_bar(nosZ_phylo_pct, fill = "Genus", x = "Sampling_ID") + 
  geom_bar(stat="identity", position="stack") +
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 3)) + 
  theme(legend.position = "top", legend.text = element_text(size = 9, family = "Lato")) +
  labs(x = "Sample", y = "Abundance (%)")  
nosZ_rel_abund_tax_bar ############################################################################################### Manuscript Fig 6 TOP
nosZ_rel_abund_tax_bar_type <- plot_bar(nosZ_phylo_pct, fill = "Genus", x = "Sampling_ID") + facet_grid(~Sample_Type, scales = "free") + 
  geom_bar(stat="identity", position="stack") +
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 3)) + 
  theme(legend.position = "top", legend.text = element_text(size = 9, family = "Lato")) +
  labs(x = "Sample", y = "Abundance (%)")
nosZ_rel_abund_tax_bar_depth <- plot_bar(nosZ_phylo_pct, fill = "Genus", x = "Sampling_ID") + facet_grid(~Depth, scales = "free") + 
  geom_bar(stat="identity", position="stack") +
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 3)) + 
  theme(legend.position = "top", legend.text = element_text(size = 9, family = "Lato")) +
  labs(x = "Sample", y = "Abundance (%)")
nosZ_rel_abund_tax_bar_depth_type <- plot_bar(nosZ_phylo_pct, fill = "Genus", x = "Sampling_ID") + facet_grid(Sample_Type~Depth, scales = "free") + 
  geom_bar(stat="identity", position="stack") +
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 3)) + 
  theme(legend.position = "top", legend.text = element_text(size = 9, family = "Lato")) +
  labs(x = "Sample", y = "Abundance (%)")

ggsave("200512 nosZ rel abund Tax Bar.png", plot = nosZ_rel_abund_tax_bar, device = "png", 
       scale = 1, width = 8.5, height = 6.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200512 nosZ rel abund Tax Bar.pdf", plot = nosZ_rel_abund_tax_bar, device = cairo_pdf, 
scale = 1, width = 12, height = 8.5, units = "in",
dpi = 600)

ggsave("200421 nosZ rel abund Tax Bar by type.png", plot = nosZ_rel_abund_tax_bar_type, device = "png", 
       scale = 1, width = 8.5, height = 6.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200421 nosZ rel abund Tax Bar by depth.png", plot = nosZ_rel_abund_tax_bar_depth, device = "png", 
       scale = 1, width = 8.5, height = 6.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200421 nosZ rel abund Tax Bar by depth Type.png", plot = nosZ_rel_abund_tax_bar_depth_type, device = "png", 
       scale = 1, width = 8.5, height = 8.5, units = "in",
       dpi = 600, type = "cairo")

#subset to lowest 10% taxon abundance...
low10_nosZ_phylo_pct_taxa <- nosZ_phylo_pct %>%  
  psmelt() %>%
  group_by(Genus)  %>%
  summarise(Abundance = sum(Abundance)) %>%
  filter(Abundance > 0)  %>%
  filter(Abundance <= 10)  %>% #retain taxa that account for less than 10% of all assigned nosZ sequences
  as.data.frame()
# create dataframe with lowest freq 10% of assigned seqs and list abundance as a % of "original" total
low10_nosZ_phylo_pct <- nosZ_phylo_pct %>%  
  psmelt() %>%
  filter(Genus %in% low10_nosZ_phylo_pct_taxa$Genus)  %>%
  filter(Abundance > 0)  %>%
  as.data.frame()

nosZ_low_abund_taxa_bar <- ggplot(data=low10_nosZ_phylo_pct, aes(x=Sampling_ID, y=Abundance/sum(Abundance)*100, fill=Genus)) + 
  geom_bar(aes(), stat="identity", position="stack",color = "black") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 3)) + theme(legend.position = "bottom",
                                                                           axis.text.x = element_text(angle = 270),
                                                                           legend.text = element_text(size = 9)) +
  labs(x = "Sample", y = "Abundance (%)")
nosZ_low_abund_taxa_bar ########################################################################################### Mansucript Fig 6 Bottom

nosZ_low_abund_taxa_bar_depth <- ggplot(data=low10_nosZ_phylo_pct, aes(x=Sampling_ID, y=Abundance/sum(Abundance)*100, fill=Genus)) + 
  geom_bar(aes(), stat="identity", position="stack",color = "black") + # facet_grid(~Sample_Type, scales = "free") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 3)) + theme(legend.position = "top",
                                                                           axis.text.x = element_text(angle = 270),
                                                                           legend.text = element_text(size = 9)) +
  labs(x = "Sample", y = "Abundance in Sample (%)") + facet_grid(~Depth)
nosZ_low_abund_taxa_bar_depth

nosZ_low_abund_taxa_bar_depth_type <- ggplot(data=low10_nosZ_phylo_pct, aes(x=Sampling_ID, y=Abundance/sum(Abundance)*100, fill=Genus)) + 
  geom_bar(aes(), stat="identity", position="stack",color = "black") + # facet_grid(~Sample_Type, scales = "free") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 3)) + theme(legend.position = "top",
                                                                           axis.text.x = element_text(angle = 270),
                                                                           legend.text = element_text(size = 9)) +
  labs(x = "Sample", y = "Abundance in Sample (%)") + facet_grid(Sample_Type~Depth, scales = "free_x")
nosZ_low_abund_taxa_bar_depth_type

ggsave("200421 nosZ low abund Tax Bar.png", plot = nosZ_low_abund_taxa_bar, device = "png", 
       scale = 1, width = 8.5, height = 6.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200421 nosZ low abund Tax Bar depth type.png", plot = nosZ_low_abund_taxa_bar_depth_type, device = "png", 
       scale = 1, width = 8.5, height = 6.5, units = "in",
       dpi = 600, type = "cairo")

#to get legends to NOT compress in ggarrange:
leg1 <- get_legend(nosZ_rel_abund_tax_bar)
nosZ_rel_abund_tax_bar <- nosZ_rel_abund_tax_bar + theme(legend.position = "none")
leg2 <- get_legend(nosZ_low_abund_taxa_bar)
nosZ_low_abund_taxa_bar <- nosZ_low_abund_taxa_bar + theme(legend.position = "none")


nosZ_all_low_tax_barplot <- ggarrange(leg1, nosZ_rel_abund_tax_bar + rremove("xlab"), nosZ_low_abund_taxa_bar, leg2, 
                                      ncol = 1, labels = c("", "A", "B", ""), heights = c(1,2,2,1))
nosZ_all_low_tax_barplot ####################################################################################################### Manuscript Fig 6
ggsave("2005121 nosZ all low abund Tax Bar.png", plot = nosZ_all_low_tax_barplot, device = "png", 
       scale = 1, width = 8.5, height = 12, units = "in",
       dpi = 600, type = "cairo")
ggsave("Fig 6 2005121 nosZ all low abund Tax Bar.pdf", plot = nosZ_all_low_tax_barplot, device = cairo_pdf, 
       scale = 1, width = 8.5, height = 12, units = "in", dpi = 600)

#### Beta diversity ##################################################################################################################################

#Beta diversity --- phyloseq & vegan
beta_nosZ_bray_dist <- distance(nosZ_phylo, method = "bray")
nosZ_bray_mds <- metaMDS(beta_nosZ_bray_dist)
nosZ_bray_mds_data <- as.data.frame(nosZ_bray_mds$points)

#To use the sample data in the plotting, we can combine the coordinate data with the sample data table: -- need to add sequencing ID labels to each to make join
nosZ_bray_mds_data$SequencingID <- rownames(nosZ_bray_mds_data)
gene_metadata$SequencingID <- rownames(gene_metadata)

nosZ_bray_mds_data <- dplyr::left_join(nosZ_bray_mds_data, gene_metadata)

#Plot MDS ordination: No good clustering based on either...
ggplot(nosZ_bray_mds_data, aes(x = MDS1, y = MDS2, color = Depth)) +
  geom_point()
ggplot(nosZ_bray_mds_data, aes(x = MDS1, y = MDS2, color = Sample_Type)) +
  geom_point()
#same thing using just phyloseq:
nosZ_ps_NMDS <- ordinate(nosZ_phylo, method = "NMDS", distance = "bray")
plot_ordination(nosZ_phylo, nosZ_ps_NMDS, type = "samples", color = "Depth", shape = "Sample_Type") #no good clusters using bray or (w)unifrac distances

nosZ_ps_PCoA <- ordinate(nosZ_phylo, "PCoA", "bray")
plot_ordination(nosZ_phylo, nosZ_ps_PCoA, color = "Depth", shape = "Sample_Type") #sep(ish) by sample type using unifrac, sep(ish) by depth with wunifrac (axes describe 45% variability)

nosZ_bray_PCoA <- ordinate(nosZ_phylo, "PCoA", "bray")
nosZ_bray_PCoA_plot <- plot_ordination(nosZ_phylo, nosZ_bray_PCoA, shape = "Sample_Type", color = "Depth") #plot ordination -- no good separations
mytitle <- expression(paste(italic("nosZ")))
basic_nosZ_PCoA <- nosZ_bray_PCoA_plot +
  ggtitle(mytitle)+
  geom_point(aes(size=3), show.legend = FALSE)+
  geom_text_repel(mapping = aes(label=Sampling_ID), size=4, show.legend = FALSE, color ='black') +
  scale_color_viridis(end = .6, discrete = TRUE) + scale_shape_discrete (name = "Sample Type")
basic_nosZ_PCoA

#exploring super important taxa...
nosZ_ASVs = as(otu_table(nosZ_phylo), "matrix") #export from phyloseq
tnosZ_ASVs<- t(nosZ_ASVs) #transpose

nosZ_BC_nmds = vegan::metaMDS(tnosZ_ASVs, distance="bray", k=2, trymax=1000) #cal distance

#ID taxa driving differences
nosZ_fit<- vegan::envfit(nosZ_BC_nmds, tnosZ_ASVs)
nosZ_fit
nosZ_df <- data.frame((nosZ_fit$vectors)$arrows, (nosZ_fit$vectors)$r, (nosZ_fit$vectors)$pvals)
write.csv(nosZ_df, "200421 taxa_Fit_nosZ.csv")
#filter out significant taxa:
nosZ_df_p0.05 <- nosZ_df %>% filter(X.nosZ_fit.vectors..pvals < 0.05)


#join ASVs to taxonomy string:
nosZ_df$ASV <- rownames(nosZ_df)
nosZ_taxonomy_join <- nosZ_taxonomy_in %>% rename(ASV = Feature.ID)
nosZ_df_tax <- left_join(nosZ_df, nosZ_taxonomy_join)
write.csv(nosZ_df_tax, "200427 ASV_taxa_Fit_nosZ.csv")
nosZ_df_tax_sig <- left_join(nosZ_df_p0.05, nosZ_taxonomy_join)
write.csv(nosZ_df_tax_sig, "200427 sig_ASV_taxa_Fit_nosZ.csv")

#Create arrows on constrained ordination.... 
bray_nosZ_dist <- phyloseq::distance(physeq = nosZ_phylo, method = "bray")
#Significance testing:
#nosZ dataframe:
nosZ_df = data.frame(sample_data(nosZ_phylo))
adonis(bray_nosZ_dist ~ Depth+Sample_Type, data = nosZ_df) #depth p = 0.001; R2 = 0.06; sample type p = 0,001, R2 = 0.05 => conclude depths & sample type have different centroids
# Homogeneity of dispersion test
nosZ_beta <- betadisper(bray_nosZ_dist, nosZ_df$Depth)
permutest(nosZ_beta) #not significant so can conclude that groups prob don't have same dispersions, lending credence to adonis output


# CAP ordinate -- doesn't work with samples that have NA in metadata -- don't use Below-ground GHG []
nosZ_cap_ord <- ordinate(
  physeq = nosZ_phylo, 
  method = "CAP",
  distance = bray_nosZ_dist,
  formula = ~ Depth+Sample_Type)

print(anova(nosZ_cap_ord)) #p = 0.001 for depth + sample type


# CAP plot
nosZ_cap_plot <- plot_ordination(physeq = nosZ_phylo, 
  ordination = nosZ_cap_ord, 
  color = "Depth") + 
  aes(shape = Sample_Type) + 
  geom_point(aes(colour = Depth), size = 4) + 
  scale_color_viridis(discrete = T, end = 0.7)
Final_nosZ_cap_plot <- nosZ_cap_plot +
  scale_shape_discrete(name = "Sample Type") + ggtitle(expression(paste(italic("nosZ"))))
Final_nosZ_cap_plot
# Now add the environmental variables as arrows (https://github.com/ccsosa/DADA2-pipeline/blob/master/Main_code_chapter2/004_2_Cluster_Beta.R)-- doesn't work 


#Beta diversity --- phyloseq & vegan for pmoA
#calculate NDMS and plot it using just phyloseq:
pmoA_ps_NMDS <- ordinate(pmoA_phylo, method = "NMDS", distance = "unifrac") #try distance "bray" (dropped values + no good clustering) "wunifrac" (no good clustering)
plot_ordination(pmoA_phylo, pmoA_ps_NMDS, type = "samples", color = "Depth") #no good clustering

#calculate PCoA with bray
pmoA_bray_PCoA <- ordinate(pmoA_phylo, "PCoA", "bray") #no good clustering
pmoA_bray_PCoA_plot <- plot_ordination(pmoA_phylo, pmoA_bray_PCoA, shape = "Sample_Type", color = "Depth") #plot ordination -- separations not explained by ANY variable
mytitle <- expression(paste(italic("pmoA")))
basic_pmoA_PCoA <- pmoA_bray_PCoA_plot +
  ggtitle(mytitle)+
  geom_point(aes(size=3), show.legend = FALSE)+
  geom_text_repel(mapping = aes(label=Sampling_ID), size=4, show.legend = FALSE, color ='black') +
  scale_color_viridis(end = .6, discrete = TRUE) + scale_shape_discrete (name = "Sample Type")
basic_pmoA_PCoA 

pmoA_nosZ_BC_PCoA <- ggarrange(basic_pmoA_PCoA, basic_nosZ_PCoA, common.legend = T)
pmoA_nosZ_BC_PCoA
ggsave("200421 pmoA nosZ BC PCoA.png", plot = pmoA_nosZ_BC_PCoA, device = "png", 
       scale = 1, width = 8.5, height = 5, units = "in",
       dpi = 600, type = "cairo")

#exploring super important taxa...
pmoA_ASVs = as(otu_table(pmoA_phylo), "matrix") #export from phyloseq
tpmoA_ASVs<- t(pmoA_ASVs) #transpose

pmoA_BC_nmds = vegan::metaMDS(tpmoA_ASVs, distance="bray", k=2, trymax=1000) #cal distance

#ID taxa driving differences
pmoA_fit<- vegan::envfit(pmoA_BC_nmds, tpmoA_ASVs)
pmoA_fit

pmoA_df <- data.frame((pmoA_fit$vectors)$arrows, (pmoA_fit$vectors)$r, (pmoA_fit$vectors)$pvals)
write.csv(pmoA_df, "200421 taxa_Fit_pmoA.csv")
#filter out significant taxa:
pmoA_df_p0.05 <- pmoA_df %>% filter(X.pmoA_fit.vectors..pvals < 0.05)

#join ASVs to taxonomy string:
pmoA_df$ASV <- rownames(pmoA_df)
pmoA_taxonomy_join <- pmoA_taxonomy_in %>% rename(ASV = Feature.ID)
pmoA_df_tax <- left_join(pmoA_df, pmoA_taxonomy_join)
write.csv(pmoA_df_tax, "200427 ASV_taxa_Fit_pmoA.csv")
pmoA_df_tax_sig <- left_join(pmoA_df_p0.05, pmoA_taxonomy_join)
write.csv(pmoA_df_tax_sig, "200427 sig_ASV_taxa_Fit_pmoA.csv")



#constrained ordination... can't seem to add arrows. Canonincal anlaysis of Principle Coordinates (CAP) = constrain the ordination axes to linear combinations of environmental variables.
      #http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
bray_pmoA_dist <- phyloseq::distance(physeq = pmoA_phylo, method = "bray")
#Significance testing:
#pmoA dataframe:
pmoA_df = data.frame(sample_data(pmoA_phylo))
adonis(bray_pmoA_dist ~ Depth+Sample_Type, data = pmoA_df) #depth p = 0.015; R2 = 0.083;  => conclude depths have different centroids
# Homogeneity of dispersion test
pmoA_beta <- betadisper(bray_pmoA_dist, pmoA_df$Depth)
permutest(pmoA_beta) #not significant so can conclude that groups prob don't have same dispersions, lending credence to adonis output


# CAP ordinate -- doesn't work with samples that have NA in metadata -- don't use Below-ground GHG []
pmoA_cap_ord <- ordinate(
  physeq = pmoA_phylo, 
  method = "CAP",
  distance = bray_pmoA_dist,
  formula = ~ Depth+Sample_Type)
print(anova(pmoA_cap_ord)) #p = 0.01 for depth + sample type
scor = vegan::scores(pmoA_cap_ord, display=c("sp", "cn", "bp"), scaling=2) 

# CAP plot
pmoA_cap_plot <- plot_ordination(physeq = pmoA_phylo, 
  ordination = pmoA_cap_ord, color = "Depth") + 
  aes(shape = Sample_Type) + 
  geom_point(aes(colour = Depth), size = 4) +
  scale_color_viridis(discrete = T, end = 0.7) 
Final_pmoA_cap_plot <- pmoA_cap_plot +
  scale_shape_discrete(name = "Sample Type") + ggtitle(expression(paste(italic("pmoA"))))
Final_pmoA_cap_plot

# Now add the environmental variables as arrows (http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html)
pmoA_arrowmat <- vegan::scores(pmoA_cap_ord, display = "bp") #SKIP b/c problem: creates 3 categories, rather than 2 vectors with just depth or sample type; "Cn" represents centroids


pmoA_nosZ_cap <- ggarrange(Final_pmoA_cap_plot, Final_nosZ_cap_plot, common.legend = T)
pmoA_nosZ_cap 
ggsave("200421 pmoA nosZ cap plot.png", plot = pmoA_nosZ_cap, device = "png", 
       scale = 1, width = 8.5, height = 5, units = "in",
       dpi = 600, type = "cairo")


#### Dendrograms: Clustering samples based on bray-curtis dissimilarity ####
#Need transposed ASV tables (nosZ_ASVs_t & pmoA_ASVs_t; grouping info is in metadata)
  #nosZ:
nosZ_betad<-vegdist(nosZ_ASVs_t,method="bray")

# Use Adonis to test for overall differences --- how to look at?!
nosZ_res_adonis <- adonis(nosZ_betad ~ Depth, gene_metadata) #999 permutations
print(nosZ_res_adonis) #sig diff among Sample_Type (p = 0.001; R2 = .04), Depth (p = 0.003; R2 = .06) ... 6% of variation described by sample depth, 4% by sample type
#no sig diffs by sampling_ID, Location, BD, Biomass, OM, GHG flux

#Cluster the samples
nosZ_hc <- hclust(nosZ_betad, method = "average")
nosZ_hc_d <- dendro_data(nosZ_hc)
#add Sampling_ID to "labels" in "Type" Column
nosZ_hc_d$labels$Sampling_ID<-gene_metadata[as.character(nosZ_hc_d$labels$label),1]
nosZ_hc_d$labels$Location<-gene_metadata[as.character(nosZ_hc_d$labels$label),2]
nosZ_hc_d$labels$Depth<-gene_metadata[as.character(nosZ_hc_d$labels$label),4]
nosZ_hc_d$labels$Sample_Type<-gene_metadata[as.character(nosZ_hc_d$labels$label),5]

  #pmoA:
pmoA_betad<-vegdist(pmoA_ASVs_t,method="bray")

# Use Adonis to test for overall differences
pmoA_res_adonis <- adonis(pmoA_betad ~ Depth, gene_metadata[c(1:35,37:41,43),]) #999 permutations
print(pmoA_res_adonis) #sig diff among Depth (p = 0.039; R2 = .07)
#no sig diffs by sampling_ID, Sample_Type, Location, BD, Biomass, OM, GHG flux

#Cluster the samples
pmoA_hc <- hclust(pmoA_betad, method = "average")
pmoA_hc_d <- dendro_data(pmoA_hc)
#add Sampling_ID to "labels" in "Type" Column
pmoA_hc_d$labels$Sampling_ID<-gene_metadata[as.character(pmoA_hc_d$labels$label),1]
pmoA_hc_d$labels$Location<-gene_metadata[as.character(pmoA_hc_d$labels$label),2]
pmoA_hc_d$labels$Depth<-gene_metadata[as.character(pmoA_hc_d$labels$label),4]
pmoA_hc_d$labels$Sample_Type<-gene_metadata[as.character(pmoA_hc_d$labels$label),5]

## Plot clusters (rarified)

nosZ_dendro <- ggplot(data = segment(nosZ_hc_d)) +
  geom_segment(data = segment(nosZ_hc_d), aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = nosZ_hc_d$labels, aes(x, y, label = Sampling_ID, color = Sample_Type), hjust = 1) +
  geom_point(data = nosZ_hc_d$labels, aes(x, y, color = Sample_Type), shape = "") + 
  coord_flip() +
  ylab("Distance (Bray-Curtis)") +
  ggtitle(expression(italic("nosZ"))) +
  theme_classic(base_size = 15) + ylim(-.1,1) +
  theme(text = element_text(family = "Lato"), axis.title.y = element_blank(), legend.position = "bottom", axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = "Soil Type") +
  scale_color_viridis(end = .6, discrete = TRUE)+
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15)))
nosZ_dendro  ##################################################################################################### Manuscript Fig 4 Right

pmoA_dendro <- ggplot(data = segment(pmoA_hc_d)) +
  geom_segment(data = segment(pmoA_hc_d), aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = pmoA_hc_d$labels, aes(x, y, label = Sampling_ID, color = Depth), hjust = 1) +
  geom_point(data = pmoA_hc_d$labels, aes(x, y, color = Depth), shape = "") +
  coord_flip() +
  ylab("Distance (Bray-Curtis)") +
  ggtitle(expression(italic("pmoA"))) +
  theme_classic(base_size = 15) + ylim(-.1,1) +
  theme(text = element_text(family = "Lato"), axis.title.y = element_blank(), legend.position = "bottom", axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.text.y = element_blank()) +
  scale_color_viridis(begin = 0.25, end = .74, discrete = TRUE) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15)))
pmoA_dendro   ##################################################################################################### Manuscript Fig 4 Left

#to get legends to NOT compress in ggarrange:
leg3 <- get_legend(pmoA_dendro)
pmoA_dendro <- pmoA_dendro + theme(legend.position = "none")
leg4 <- get_legend(nosZ_dendro)
nosZ_dendro <- nosZ_dendro + theme(legend.position = "none")

dendros <- ggarrange(pmoA_dendro, nosZ_dendro, leg3, leg4, ncol = 2, nrow = 2, heights = c(7,1))  
dendros ####################################################################################################################### Manuscript Fig 4
ggsave("200512 BC dendrograms.png", plot = dendros, device = "png", 
       scale = 1, width = 7, height = 8.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("Fig 4 200512 BC dendrograms.pdf", plot = dendros, device = cairo_pdf, 
       scale = 1, width = 7.1, height = 8.5, units = "in", dpi = 600)


#### LEfSe Analysis: #### 
#-- in order for GALAXY to read file into LEfSe format to run analysis, need to make sure file does not contain sample IDs.
#LDA Effect Size (LEfSe) (Segata et. al 2010) is An algorithm for High-Dimensional biomarker discovery and explanation that identifies genomic features (genes, pathways, or taxa) characterizing the differences between two or more biological conditions (or classes, see figure below). It emphasizes both statistical significance and biological relevance, allowing researchers to identify differentially abundant features that are also consistent with biologically meaningful categories (subclasses). LEfSe first robustly identifies features that are statistically different among biological classes. It then performs additional tests to assess whether these differences are consistent with respect to expected biological behavior.
      #Specifically, we first use the non-parametric factorial Kruskal-Wallis (KW) sum-rank test to detect features with significant differential abundance with respect to the class of interest; biological significance is subsequently investigated using a set of pairwise tests among subclasses using the (unpaired) Wilcoxon rank-sum test. As a last step, LEfSe uses Linear Discriminant Analysis to estimate the effect size of each differentially abundant feature and, if desired by the investigator, to perform dimension reduction.
## pmoA
# this script converts phyloseq object into lefse recognized file format https://wipperman.github.io/TBRU/Treatment/#part-3-analyze-data-with-lefse
sample.data <- phyloseq::sample_data(pmoA_phylo) %>% data.frame(stringsAsFactors = FALSE)
sample.data$sample <- rownames(sample.data)
#sample.data$Sample_Type <- gsub(" ", "_", sample.data$Sample_Type) #remove spaces because LEfSe galaxy site can't read DFs with special characters?
#
keepvars <- c("Depth","Sample_Type") #don't include sample ID rows...?! "sample","Sampling_ID"
keepvars <- unique(keepvars[!is.na(keepvars)])
lefse.samp <- sample.data[, keepvars]
#
sample0 <- t(lefse.samp) %>% as.matrix()
colnames(sample0) <- sample0[1,]
sample0 <- as.data.frame(sample0) 
#
data0 <- otu_table(pmoA_phylo) %>% as.data.frame(stringsAsFactors = FALSE)
data1 <- data0 %>% as.data.frame(keep.rownames=T, stringsAsFactors = FALSE) %>% t() %>% as.data.frame(keep.rownames=T)
sample1 <- sample0 %>% as.data.frame(keep.rownames=T, stringsAsFactors = FALSE) %>% t() %>% as.data.frame(keep.rownames=T)
pre.lefse <- bind_cols(sample1, data1) %>% t() %>% as.data.frame() %>% rownames_to_column()

#writes table for LEfSe to run on...
write.table(pre.lefse,file ="200423 pmoA_lefse.txt",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #can be uploaded to https://huttenhower.sph.harvard.edu/galaxy/
#ERROR: "Traceback (most recent call last):
# File "/shed_tools/testtoolshed.g2.bx.psu.edu/repos/george-weingart/lefse/a6284ef17bf3/lefse/format_input.py", line 435, in <module>
#  feats = numerical_values(feats,params['norm_v'])
# File "/shed_tools/testtoolshed --- resolved by dropping Sampling_ID row

#OUTPUT: no significant features found, apparently : Number of significantly discriminative features: 0 ( 1 ) before internal wilcoxon
    #No features with significant differences between the two classes
    #Number of discriminative features with abs LDA score > 2.0 : 0

## nosZ
# this script converts phyloseq object into lefse recognized file format https://wipperman.github.io/TBRU/Treatment/#part-3-analyze-data-with-lefse
sample.data <- phyloseq::sample_data(nosZ_phylo) %>% data.frame(stringsAsFactors = FALSE)
sample.data$sample <- rownames(sample.data)
sample.data$Sample_Type <- gsub(" ", "_", sample.data$Sample_Type) #remove spaces because LEfSe galaxy site can't read DFs with special characters?
#
keepvars <- c("Depth","Sample_Type") #don't include sample ID rows...?! "sample","Sampling_ID"
keepvars <- unique(keepvars[!is.na(keepvars)])
lefse.samp <- sample.data[, keepvars]
#
sample0 <- t(lefse.samp) %>% as.matrix()
colnames(sample0) <- sample0[1,]
sample0 <- as.data.frame(sample0) 
#
data0 <- otu_table(nosZ_phylo) %>% as.data.frame(stringsAsFactors = FALSE)
data1 <- data0 %>% as.data.frame(keep.rownames=T, stringsAsFactors = FALSE) %>% t() %>% as.data.frame(keep.rownames=T)
sample1 <- sample0 %>% as.data.frame(keep.rownames=T, stringsAsFactors = FALSE) %>% t() %>% as.data.frame(keep.rownames=T)
pre.lefse <- bind_cols(sample1, data1) %>% t() %>% as.data.frame() %>% rownames_to_column()

#writes table for LEfSe to run on...
write.table(pre.lefse,file ="200423 nosZ_lefse.txt",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #can be uploaded to https://huttenhower.sph.harvard.edu/galaxy/
#output: no significant features found: Number of significantly discriminative features: 0 ( 41 ) before internal wilcoxon
    #No features with significant differences between the two classes
    #Number of discriminative features with abs LDA score > 2.0 : 0


#### SPIEC-EASI Network Analysis #### 
#https://github.com/zdk123/SpiecEasi
#unrarified pmoA
pmoA_raw_se.mb <- spiec.easi(pmoA_phylo_raw, method='mb', lambda.min.ratio=1e-2,
                           nlambda=50, pulsar.params=list(rep.num=40))
pmoA_raw_ig2.mb <- adj2igraph(getRefit(pmoA_raw_se.mb),  vertex.attr=list(name=taxa_names(pmoA_phylo_raw)))
plot_network(pmoA_raw_ig2.mb, pmoA_phylo_raw, type='taxa', color="Genus") #color network dots by genus ... nothing jumps out
plot_network(pmoA_raw_ig2.mb, pmoA_phylo_raw, color="Depth") #some weirdo ggplot error Must request at least one colour from a hue palette
# examine ROC over lambda path and PR over the stars index for the selected graph
getStability(pmoA_raw_se.mb)#0.048 #since we have only a rough, discrete sampling of networks along the lambda path, we should check how far we are from the target stability threshold (0.05).
sum(getRefit(pmoA_raw_se.mb))/2 #149...nodes?

#pmoA_raw_se.gl <- spiec.easi(pmoA_phylo_raw, method='glasso', lambda.min.ratio=1e-2,
#                             nlambda=50, pulsar.params=list(rep.num=40)) #seems to take MUCH longer than the 'mb' method
#pmoA_raw_ig2.gl <- adj2igraph(getRefit(pmoA_raw_se.gl),  vertex.attr=list(name=taxa_names(pmoA_phylo_raw)))
#plot_network(pmoA_raw_ig2.gl, pmoA_phylo_raw, type='taxa', color="Genus") #color network dots by genus ... nothing jumps out
#getStability(pmoA_raw_se.gl)#0.048 #since we have only a rough, discrete sampling of networks along the lambda path, we should check how far we are from the target stability threshold (0.05).
#sum(getRefit(pmoA_raw_se.gl))/2 #149...nodes?

#rarified pmoA -- on difference
pmoA_rare_se.mb <- spiec.easi(pmoA_phylo, method='mb', lambda.min.ratio=1e-2,
                             nlambda=20, pulsar.params=list(rep.num=40))
pmoA_rare_ig2.mb <- adj2igraph(getRefit(pmoA_rare_se.mb),  vertex.attr=list(name=taxa_names(pmoA_phylo)))
plot_network(pmoA_rare_ig2.mb, pmoA_phylo, type='taxa', color="Genus")

#unrarified nosZ
nosZ_raw_se.mb <- spiec.easi(nosZ_phylo_raw, method='mb', lambda.min.ratio=1e-2,
                             nlambda=20, pulsar.params=list(rep.num=20)) # Takes FOREVER to run - crashes before it completes
# examine ROC over lambda path and PR over the stars index for the selected graph
getStability(nosZ_raw_se.mb)#0.048 # check how far we are from the target stability threshold (0.05).
sum(getRefit(nosZ_raw_se.mb))/2 #149...nodes?

nosZ_raw_ig2.mb <- adj2igraph(getRefit(nosZ_raw_se.mb),  vertex.attr=list(name=taxa_names(nosZ_phylo_raw)))
plot_network(nosZ_raw_ig2.mb, nosZ_phylo_raw, type='taxa', color="Genus") #color network dots by genus ... nothing jumps out


