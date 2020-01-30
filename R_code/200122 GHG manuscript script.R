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
library(vegan) # for calculating distances for dendrogram etc.
library(ggdendro) # for making dendrogram
library(phyloseq)


sessionInfo()
#R version 3.6.2 (2019-12-12)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows >= 8 x64 (build 9200)
#[1] microbiome_1.8.0   phyloseq_1.30.0    ggpubr_0.2.4       magrittr_1.5       viridis_0.5.1      viridisLite_0.3.0 
#[7] ggiraphExtra_0.2.9 forcats_0.4.0      stringr_1.4.0      dplyr_0.8.3        purrr_0.3.3        readr_1.3.1       
#[13] tidyr_1.0.0        tibble_2.1.3       ggplot2_3.2.1      tidyverse_1.3.0    psych_1.9.12.31 

#Read in files:
BD_Veg_flux_data <- read.csv("raw_data/190913 GHG Soil Summary.csv")

nosZ_taxonomy_table <- read.csv("raw_data/200107 nosZ taxonomy table.csv") #samples v taxonomy...cleaned up taxonomy from original QIIME2 files
pmoA_taxonomy_table <- read.csv("raw_data/200122 pmoA taxonomy table.csv") #samples v taxonomy...cleaned up taxonomy from original QIIME2 files -- including replacing "streptomyces" with "Bacteria"
nosZ_ASVs_in <- read.csv("raw_data/200107 nosZ ASV table.csv")
pmoA_ASVs_in <- read.csv("raw_data/200107 pmoA ASV table.csv") #filtered from from original dada2 outputs files to only keep samples with >=1 feature
nosZ_taxonomy_in <- read.delim("raw_data/200109 nosZ_taxonomy.tsv", sep = "\t") #cleaned up taxonomy from original QIIME2 files -- removed extra ;; 
pmoA_taxonomy_in <- read.delim("raw_data/200108 pmoA_taxonomy.tsv", sep = "\t") #cleaned up taxonomy from original QIIME2 files -- removed extra ;; 
  #replace streptomyces with "Bacteria" Taxonomic assignment: Bacteria;Actinobacteria;Streptomycetales;Streptomycetaceae;Streptomyces
pmoA_taxonomy_in$Taxon <- sapply(pmoA_taxonomy_in$Taxon,function(x) {x <- gsub("Bacteria;Actinobacteria;Streptomycetales;Streptomycetaceae;Streptomyces", "Bacteria",x)})

nosZ_rtree <- read_tree("raw_data/nosZ_rooted_tree.nwk")
pmoA_rtree <- read_tree("raw_data/pmoA_rooted_tree.nwk")
gene_metadata <- read.delim("raw_data/200113_AHC_sequencing_sample_GHG_metadata.txt", sep = "\t", row.names = 1)


#change ASV & taxonomy tables to be phyloseq-compatible (numeric matix)
nosZ_ASVs <- nosZ_ASVs_in %>%
  select(-OTU_ID_nosZ) %>%
  as.matrix
rownames(nosZ_ASVs) <- nosZ_ASVs_in$OTU_ID_nosZ
nosZ_taxonomy <- nosZ_taxonomy_in %>%
  select(Taxon) %>%
  as.matrix
rownames(nosZ_taxonomy) <- nosZ_taxonomy_in$Feature.ID

pmoA_ASVs <- pmoA_ASVs_in %>%
  select(-OTU_ID_pmoA) %>%
  as.matrix
rownames(pmoA_ASVs) <- pmoA_ASVs_in$OTU_ID_pmoA
pmoA_taxonomy <- pmoA_taxonomy_in %>%
  select(Taxon) %>%
  as.matrix
rownames(pmoA_taxonomy) <- pmoA_taxonomy_in$Feature.ID

# for dendrogram / non-phyloseq analyses
nosZ_ASVs_t <- t(nosZ_ASVs)
pmoA_ASVs_t <- t(pmoA_ASVs)

#use phyloseq() to combine text files (otu tables, taxonomy, sample_data etc.) into phyloseq object... https://joey711.github.io/phyloseq/import-data.html
  #structure: # phyloseq(otu_table(GP), phy_tree(GP), tax_table(GP), sample_data(GP))
nosZ_phylo <- phyloseq(otu_table(nosZ_ASVs, taxa_are_rows = T), phy_tree(nosZ_rtree), tax_table(nosZ_taxonomy), sample_data(gene_metadata))
pmoA_phylo <- phyloseq(otu_table(pmoA_ASVs, taxa_are_rows = T), phy_tree(pmoA_rtree), tax_table(pmoA_taxonomy), sample_data(gene_metadata))
nosZ_phylo_pct <- transform_sample_counts(nosZ_phylo, function(x) 100 * x/sum(x))
pmoA_phylo_pct <- transform_sample_counts(pmoA_phylo, function(x) 100 * x/sum(x))

#do some basic ggplot themesettings:
theme_set(theme_classic(base_size = 15))

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
  labs(y=expression(paste(CH[4]," flux (", mu,mol/m^2/h,")"))) +
  geom_hline(yintercept = 0, color = 'grey', linetype = 2)
CH4_pt

CO2_pt <- ggplot(data=BD_Veg_flux_data, aes(x=Location, y=CO2_Flux_umol_m.2_h)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15)+
  labs(y=expression(paste(CO[2]," flux (", mu,mol/m^2/h,")"))) +
  ggtitle("Flux")
CO2_pt

N2O_pt <- ggplot(data=BD_Veg_flux_data, aes(x=Location, y=N2O_Flux_umol_m.2_h)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15)+
  labs(y=expression(paste(N[2]*O," flux (", mu,mol/m^2/h,")")))
N2O_pt
GHGs_pt <- ggarrange(CO2_pt, CH4_pt, N2O_pt, ncol = 1, align = "v")
GHGs_pt

ggsave("200121 GHG Pts.png", plot = GHGs_pt, device = "png", path = NULL,
       scale = 1, width = 4, height = 7.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200121 GHG Pts.pdf", plot = GHGs_pt, device = cairo_pdf, 
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
  labs(y=expression(CH[4]~concentration~(PPM)))
CH4_below_pt

CO2_below_pt <- ggplot(data=subset(BD_Veg_flux_data, Location=='Control'|Location=='Experimental'), aes(x=Location, y=C02_conc_6in_PPM)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15)+
  labs(y=expression(CO[2]~concentration~(PPM))) +
  ggtitle("Concentration")
CO2_below_pt

N2O_below_pt <- ggplot(data=subset(BD_Veg_flux_data, Location=='Control'|Location=='Experimental'), aes(x=Location, y=N2O_conc_6in_PPM)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15)+
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
Flux_v_belowground_GHG

ggsave("200121 GHG above v below Pts.png", plot = Flux_v_belowground_GHG, device = "png", 
       scale = 1, width = 12, height = 8.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200121 GHG above v below Pts.pdf", plot = Flux_v_belowground_GHG, device = cairo_pdf, 
       scale = 1, width = 12, height = 8.5, units = "in",
       dpi = 600)

#### Phyoseq alpha diversity, taxonomy & other stuff #########################################################################################################

#pmoA:
#Alpha Diversity:
pmoA_alphadiv <- estimate_richness(pmoA_phylo)
pmoA_alpha_plots <- plot_richness(pmoA_phylo, measures = c("Observed", "Shannon", "Simpson"), x = "Sampling_ID", nrow = 3) + xlab("Sample") + ggtitle(expression(italic("pmoA")))
plot_richness(pmoA_phylo, measures = c("Observed", "Shannon", "Simpson"), x = "Depth", 
              color = "Location") + scale_color_viridis(discrete = T) + 
  geom_boxplot(aes(group = Depth)) +
  theme(legend.position = "bottom")
#stacked bar plots
  #total abundance
plot_bar(pmoA_phylo, fill = "Taxon", x = "Sampling_ID") + facet_grid(~Depth, scales = "free") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + theme(legend.position = "top") +
  labs(x = "Sample")
  #relative abundance
pmoA_rel_abund_tax_bar <- plot_bar(pmoA_phylo_pct, fill = "Taxon", x = "Sampling_ID") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + 
  theme(legend.position = "top", legend.text = element_text(size = 9)) +
  labs(x = "Sample", y = "Abundance (%)")
pmoA_rel_abund_tax_bar_depth <- plot_bar(pmoA_phylo_pct, fill = "Taxon", x = "Sampling_ID") + facet_grid(~Depth, scales = "free") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + 
  theme(legend.position = "top", legend.text = element_text(size = 9)) +
  labs(x = "Sample", y = "Abundance (%)")
ggsave("200122 pmoA rel abund Tax Bar.png", plot = pmoA_rel_abund_tax_bar, device = "png", 
       scale = 1, width = 8.5, height = 6.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200122 pmoA rel abund Tax Bar by depth.png", plot = pmoA_rel_abund_tax_bar_depth, device = "png", 
       scale = 1, width = 8.5, height = 6.5, units = "in",
       dpi = 600, type = "cairo")



#nosZ:
#Alpha Diversity:
nosZ_alphadiv <- estimate_richness(nosZ_phylo)
nosZ_alpha_plots <- plot_richness(nosZ_phylo, measures = c("Observed", "Shannon", "Simpson"), x = "Sampling_ID", nrow = 3) + xlab("Sample") + ggtitle(expression(italic("nosZ")))
all_alpha_plots <- ggarrange(pmoA_alpha_plots, nosZ_alpha_plots, nrow = 2)
ggsave("200127 all alpha diversity plots.png", plot = all_alpha_plots, device = "png", 
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

#stacked bar plots
#total abundance
plot_bar(nosZ_phylo, fill = "Taxon", x = "Sampling_ID") + facet_grid(~Sample_Type, scales = "free") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + theme(legend.position = "top") +
  labs(x = "Sample")
#relative abundance
nosZ_rel_abund_tax_bar <- plot_bar(nosZ_phylo_pct, fill = "Taxon", x = "Sampling_ID") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + 
  theme(legend.position = "top", legend.text = element_text(size = 9)) +
  labs(x = "Sample", y = "Abundance (%)")
nosZ_rel_abund_tax_bar_type <- plot_bar(nosZ_phylo_pct, fill = "Taxon", x = "Sampling_ID") + facet_grid(~Sample_Type, scales = "free") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + 
  theme(legend.position = "top", legend.text = element_text(size = 9)) +
  labs(x = "Sample", y = "Abundance (%)")
nosZ_rel_abund_tax_bar_depth <- plot_bar(nosZ_phylo_pct, fill = "Taxon", x = "Sampling_ID") + facet_grid(~Depth, scales = "free") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + 
  theme(legend.position = "top", legend.text = element_text(size = 9)) +
  labs(x = "Sample", y = "Abundance (%)")
nosZ_rel_abund_tax_bar_depth_type <- plot_bar(nosZ_phylo_pct, fill = "Taxon", x = "Sampling_ID") + facet_grid(Sample_Type~Depth, scales = "free") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + 
  theme(legend.position = "top", legend.text = element_text(size = 9)) +
  labs(x = "Sample", y = "Abundance (%)")

ggsave("200122 nosZ rel abund Tax Bar.png", plot = nosZ_rel_abund_tax_bar, device = "png", 
       scale = 1, width = 9.5, height = 10.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200122 nosZ rel abund Tax Bar by type.png", plot = nosZ_rel_abund_tax_bar_type, device = "png", 
       scale = 1, width = 9.5, height = 10.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200122 nosZ rel abund Tax Bar by depth.png", plot = nosZ_rel_abund_tax_bar_depth, device = "png", 
       scale = 1, width = 9.5, height = 10.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200122 nosZ rel abund Tax Bar by depth Type.png", plot = nosZ_rel_abund_tax_bar_depth_type, device = "png", 
       scale = 1, width = 9.5, height = 12.5, units = "in",
       dpi = 600, type = "cairo")

#subset to lowest 10% taxon abundance...
low10_nosZ_phylo_pct_taxa <- nosZ_phylo_pct %>%  
  psmelt() %>%
  group_by(Taxon)  %>%
  summarise(Abundance = sum(Abundance)) %>%
  filter(Abundance > 0)  %>%
  filter(Abundance <= 10)  %>% #retain taxa that account for less than 10% of all assigned nosZ sequences
  as.data.frame()
# create dataframe with lowest freq 10% of assigned seqs and list abundance as a % of "new" total
low10_nosZ_phylo_pct <- nosZ_phylo_pct %>%  
  psmelt() %>%
  filter(Taxon %in% low10_nosZ_phylo_pct_taxa$Taxon)  %>%
  filter(Abundance > 0)  %>%
  as.data.frame()

nosZ_low_abund_taxa_bar <- ggplot(data=low10_nosZ_phylo_pct, aes(x=Sampling_ID, y=Abundance/sum(Abundance)*100, fill=Taxon)) + 
  geom_bar(aes(), stat="identity", position="stack",color = "black") + # facet_grid(~Sample_Type, scales = "free") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + theme(legend.position = "top",
                                                                           axis.text.x = element_text(angle = 270),
                                                                           legend.text = element_text(size = 9)) +
  labs(x = "Sample", y = "Abundance (%)")
nosZ_low_abund_taxa_bar

nosZ_low_abund_taxa_bar_depth <- ggplot(data=low10_nosZ_phylo_pct, aes(x=Sampling_ID, y=Abundance/sum(Abundance)*100, fill=Taxon)) + 
  geom_bar(aes(), stat="identity", position="stack",color = "black") + # facet_grid(~Sample_Type, scales = "free") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + theme(legend.position = "top",
                                                                           axis.text.x = element_text(angle = 270),
                                                                           legend.text = element_text(size = 9)) +
  labs(x = "Sample", y = "Abundance (%)") + facet_grid(~Depth)
nosZ_low_abund_taxa_bar_depth

nosZ_low_abund_taxa_bar_depth_type <- ggplot(data=low10_nosZ_phylo_pct, aes(x=Sampling_ID, y=Abundance/sum(Abundance)*100, fill=Taxon)) + 
  geom_bar(aes(), stat="identity", position="stack",color = "black") + # facet_grid(~Sample_Type, scales = "free") + 
  scale_fill_viridis(discrete = T, guide = guide_legend(ncol = 1)) + theme(legend.position = "top",
                                                                           axis.text.x = element_text(angle = 270),
                                                                           legend.text = element_text(size = 9)) +
  labs(x = "Sample", y = "Abundance (%)") + facet_grid(Sample_Type~Depth, scales = "free_x")
nosZ_low_abund_taxa_bar_depth_type

ggsave("200122 nosZ low abund Tax Bar.png", plot = nosZ_low_abund_taxa_bar, device = "png", 
       scale = 1, width = 9.5, height = 10.5, units = "in",
       dpi = 600, type = "cairo")
ggsave("200122 nosZ low abund Tax Bar depth type.png", plot = nosZ_low_abund_taxa_bar_depth_type, device = "png", 
       scale = 1, width = 9.5, height = 10.5, units = "in",
       dpi = 600, type = "cairo")



plot_tree(nosZ_phylo, color = "Depth", shape = "Location", label.tips = "Sampling_ID", 
          ladderize = TRUE, nodelabf=nodeplotblank) # To make circular: + coord_polar(theta = "y")

nosZ_network <- make_network(nosZ_phylo, type = "samples", distance = "unifrac", max.dist = 0.7)
plot_network(nosZ_network)


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

#Beta diversity --- phyloseq & vegan for pmoA
#calculate NDMS and plot it using just phyloseq:
pmoA_ps_NMDS <- ordinate(pmoA_phylo, method = "NMDS", distance = "unifrac") #try distance "bray" (dropped values + no good clustering) "wunifrac" (no good clustering)
plot_ordination(pmoA_phylo, pmoA_ps_NMDS, type = "samples", color = "Depth") #no good clustering
pmoA_ps_PCoA <- ordinate(pmoA_phylo, "PCoA", "bray")  #try distance "unifrac" (prob interpeting cuz extra rows) "wunifrac" (prob interpreting)
plot_ordination(pmoA_phylo, pmoA_ps_PCoA, color = "Depth") #no good clustering

#### Dendrograms: Clustering samples based on bray-curtis dissimilarity ####
#Need transposed ASV tables (nosZ_ASVs_t & pmoA_ASVs_t; grouping info is in metadata)
  #nosZ:
nosZ_betad<-vegdist(nosZ_ASVs_t,method="bray")

# Use Adonis to test for overall differences --- how to look at?!
nosZ_res_adonis <- adonis(nosZ_betad ~ Sample_Type, gene_metadata) #999 permutations
print(nosZ_res_adonis) #sig diff among Sample_Type (p = 0.001; R2 = .04), Depth (p = 0.003; R2 = .06) ... 6% of variablion described by sample depth, 4% by sample type
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
print(pmoA_res_adonis) #sig diff among Depth (p = 0.028; R2 = .07)
#no sig diffs by sampling_ID, Sample_Type, Location, BD, Biomass, OM, GHG flux

#Cluster the samples
pmoA_hc <- hclust(pmoA_betad, method = "average")
pmoA_hc_d <- dendro_data(pmoA_hc)
#add Sampling_ID to "labels" in "Type" Column
pmoA_hc_d$labels$Sampling_ID<-gene_metadata[as.character(pmoA_hc_d$labels$label),1]
pmoA_hc_d$labels$Location<-gene_metadata[as.character(pmoA_hc_d$labels$label),2]
pmoA_hc_d$labels$Depth<-gene_metadata[as.character(pmoA_hc_d$labels$label),4]
pmoA_hc_d$labels$Sample_Type<-gene_metadata[as.character(pmoA_hc_d$labels$label),5]

## Plot clusters

nosZ_dendro <- ggplot(data = segment(nosZ_hc_d)) +
  geom_segment(data = segment(nosZ_hc_d), aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = nosZ_hc_d$labels, aes(x, y, label = Sampling_ID, color = Sample_Type), hjust = 1) +
  coord_flip() +
  ylab("Distance (Bray-Curtis)") +
  ggtitle(expression(italic("nosZ"))) +
  theme_classic(base_size = 15) + ylim(-.1,1) +
  theme(axis.title.y = element_blank(), legend.position = "bottom", axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.text.y = element_blank()) +
  labs(color = "Soil Type") +
  scale_color_viridis(end = .6, discrete = TRUE)
nosZ_dendro  

pmoA_dendro <- ggplot(data = segment(pmoA_hc_d)) +
  geom_segment(data = segment(pmoA_hc_d), aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = pmoA_hc_d$labels, aes(x, y, label = Sampling_ID, color = Depth), hjust = 1) +
  coord_flip() +
  ylab("Distance (Bray-Curtis)") +
  ggtitle(expression(italic("pmoA"))) +
  theme_classic(base_size = 15) + ylim(-.1,1) +
  theme(axis.title.y = element_blank(), legend.position = "bottom", axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.text.y = element_blank()) +
  scale_color_viridis(end = .65, discrete = TRUE)
pmoA_dendro  

dendros <- ggarrange(pmoA_dendro, nosZ_dendro)  
dendros
ggsave("200122 BC dendrograms.png", plot = dendros, device = "png", 
       scale = 1, width = 7, height = 8.5, units = "in",
       dpi = 600, type = "cairo")

