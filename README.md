# Gut microbiota in human melioidosis

This is the code used for the main analyses in "Microbial Resilience in Human Melioidosis: Understanding Gut Microbiome and Resistome Dynamics in Diagnosis to Recovery" (submitted). For questions: Bob Kullberg at r.f.j.kullberg@amsterdamumc.nl or Chiranjay Mukhopadhyay at chiranjay.m@manipal.edu

Data protection regulations do not allow public sharing of individual patient data. Please refer to the data availability statement in the manuscript for information on how to gain data access. 


## Step 1 - Load libraries

```
library(tidyverse)
library(yingtools2)
library(phyloseq)
library(rlang)
library(microbiome)
library(vegan)
library(glmnet)
library(Maaslin2)
library(RColorBrewer) 
library(cowplot)
library(ggpubr) 
library(ggrepel)
library(scales)
library(readxl)
library(ggalluvial)
library(rstatix)
library(ComplexHeatmap)
```

## Step 2 - Load data

Microbiota sequence data is already preprocessed and a count table is produced, which is integrated with the taxonomy and a phylogenetic tree using the phyloseq package (details described in the manuscript). 

```
df <- read_csv("~/metadata.csv") # metadata
phyloseq <- readRDS("~/phyloseq.RDS") # phyloseq file
amr <- read_tsv("~/amr_assemblies_melioidosis.tsv") %>%
  splitstackshape::concat.split(split.col = "#FILE", sep = "assemblies/") %>%
  splitstackshape::concat.split(split.col = "#FILE_2", sep = "_final") %>%
  dplyr::rename(sample = `#FILE_2_1`) %>%
  right_join(df) %>%
  dplyr::select(sample, group, group_time, reads, SEQUENCE, GENE, `%COVERAGE`, `%IDENTITY`, PRODUCT, RESISTANCE) %>%
  mutate(RESISTANCE = if_else(RESISTANCE == "CEPHALOSPORIN", "BETA-LACTAM", RESISTANCE)) # these categories were combined given their overlap
```

The main column is “group_time”, which classifies samples into controls and melioidosis patients at the different timepoints (diagnosis, day 14/discharge, 6 months)
```
table(df$group_time)
```


## Step 3 - Microbiota composition (figure 1)
Assess beta diversity:
```
set.seed(88)
bray <- phyloseq::distance(phyloseq, method = "bray") 
df.bray <- as.matrix(bray)
df.bray <- as.data.frame(row.names(df.bray)) %>%
  mutate(sample = `row.names(df.bray)`) %>%
  left_join(df)

#plot_ordination(phyloseq, ord, type="samples", color="group_time")

bray.plot <- plot_ordination(phyloseq, ord, type = "samples", color = "group_time", justDF = T) %>% # use 'justDF = T' to get a dataframe with the  coordinates for PCoA plot
  rownames_to_column(var = "sample") %>%
  dplyr::select(Axis.1, Axis.2, sample, group_time)

centroids <- aggregate(cbind(Axis.1, Axis.2)~group_time, data=bray.plot, mean) # calculate the centroids per group
bray.plot <- merge(bray.plot, centroids, by="group_time", suffixes=c("", ".centroid")) # merge the centroid data with PCoA data

ggplot(bray.plot, aes(Axis.1, Axis.2, color = group_time)) +
  geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2, color= group_time), alpha = 0.5)+ 
  geom_point(data=bray.plot, aes(color=group_time),size=4.5,alpha=1) + # Create points for each participant
  geom_label_repel(data = centroids, aes(x=Axis.1, y=Axis.2, label=c("Control", "Melioidosis - diagnosis", "Melioidosis - day 14", "Melioidosis - 6 months")), size=4.5, fill = ggplot2::alpha(c("white"),0.92)) +
  theme_bw() +
  xlab("22.0% variance") + #Label X-axis
  ylab("16.9% variance") + #Label Y-axis
  scale_color_manual(values=c("#91c3e1", "#660033", "#FC8D59", "#FEE08B")) +
  theme(legend.position = "none")
```

To assess significance of differences between timepoints (e.g. controls vs. patients at diagnosis), permutational ANOVA was used:
```
phyloseq1 <- phyloseq
df1 <- df %>%
  filter(group_time == "control day 1" | group_time == "case day 1") # Controls vs patients at diagnosis
sample_data(phyloseq1) <- set.samp(df1)
set.seed(88) 
bray <- phyloseq::distance(phyloseq1, method = "bray") 
df.bray <- as.matrix(bray)
df.bray <- as.data.frame(row.names(df.bray)) %>%
  mutate(sample = `row.names(df.bray)`) %>%
  left_join(df1)
adonis2((bray ~ group_time), by = 'margin', data = df.bray, permutations =9999)
rm(phyloseq1, df1, bray, df.bray)
```
This analysis was repeated for the other comparisons (e.g. controls vs. patients at 6 months) using similar code.   

Next, we assessed determinants of interindividual dissimilarities in patients at diagnosis and community controls:
```
phyloseq1 <- phyloseq
df1 <- df %>%
  filter(group_time == "control day 1" | group_time == "case day 1")
sample_data(phyloseq1) <- set.samp(df1)

set.seed(88) 
bray <- phyloseq::distance(phyloseq1, method = "bray") 
df.bray <- as.matrix(bray)
df.bray <- as.data.frame(row.names(df.bray)) %>%
  mutate(sample = `row.names(df.bray)`) %>%
  left_join(df1)

variables <-c("group", "Age", "Sex", "BMI", "Any comorbidity", "Diabetes", "Hypertension", "COPD")

set.seed(88)
R2 <- list((bray ~ group),
           (bray ~ Age),
           (bray ~ Sex),
           (bray ~ BMI),
           (bray ~ Comorbidities),
           (bray ~ Diabetes),
           (bray ~ Hypertension),
           (bray ~ COPD)) %>%
  lapply(adonis2, data = df.bray, by = 'margin', permutations = 9999) %>%
  lapply(function(x) x[["R2"]][[1]]) %>%
  do.call(c, .)

Pvalues <- list((bray ~ group),
                (bray ~ Age),
                (bray ~ Sex),
                (bray ~ BMI),
                (bray ~ Comorbidities),
                (bray ~ Diabetes),
                (bray ~ Hypertension),
                (bray ~ COPD)) %>%
  lapply(adonis2, data = df.bray, by = 'margin', permutations = 9999) %>%
  lapply(function(x) x[["Pr(>F)"]][[1]]) %>%
  do.call(c, .)

dissim <- as.data.frame(cbind(variables, R2, Pvalues))
dissim <- dissim %>%
  mutate(R2 = as.numeric(as.character(R2))) %>%
  mutate(variables = fct_relevel(variables, 
                                 "group",  "Sex",  "Any comorbidity", "Diabetes",
                                 "COPD", "BMI", "Hypertension", "Age"))

ggplot(dissim, aes(x=R2, y=variables)) + 
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = scales::percent_format()) +
  xlab("Contribution to interindividual dissimilarities in microbiota (R2)") +
  ylab("") +
  theme_bw(base_size = 12)
```

Shannon diversity at the species level was calculated: 
```
P.alpha <- microbiome::aggregate_taxa(phyloseq, level = "Species") 
alpha <- estimate_richness(P.alpha)
alpha$sample <- row.names(alpha)
alpha$sample <- gsub("X","",as.character(alpha$sample))
alpha <- left_join(df, alpha)
  
rm(P.alpha)

comparisons <- list(c("control day 1", "case day 1"), c("control day 1", "case day 14"), c("control day 1", "case 6 months"))
ggplot(alpha, aes(x= group_time, y = `sample_sums(phyloseq)`, fill = group_time)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "black", pch = 21, alpha =.75, size = 3.5, width = 0.2)+
  scale_color_manual(values=c("#91c3e1", "#660033", "#FC8D59", "#FEE08B")) +
  scale_fill_manual(values=c("#91c3e1", "#660033", "#FC8D59", "#FEE08B")) +
  stat_compare_means(comparisons = comparisons) +
  theme_cowplot(font_size = 12)+
  xlab("")+
  ylab("Shannon diversity")+
  theme(legend.position = "none")
```

The relative abundance of  butyrate-producing bacteria was calculated by combining the abundances of bacteria that are known to be the most abundant drivers of butyrate production:
```
P.butyrate <- microbiome::transform(phyloseq, "clr")
butyrate.producers <- c("Butyricimonas", "Odoribacter", "Anaerostipes", "Anaerobutyricum", "Agathobacter", 
                        "Butyrivibrio_A","Coprococcus", "Coprococcus_A", "Roseburia","Shuttleworthia","Butyricicoccus",
                        "Flavonifractor", "Pseudoflavonifractor", 
                        "Faecalibacterium","Oscillibacter", "Gemmiger") 
butyrate.clr <- get.otu.melt(P.butyrate, filter.zero = F) %>%
  mutate(Genus = gsub("g__","", as.character(Genus))) %>% # remove the "g__" addition if necessary
  mutate(butyrate_producers = if_else(Genus %in% butyrate.producers, Genus, "Other", "Other")) %>%
  filter(butyrate_producers != "Other") %>%
  group_by(sample)%>%
  dplyr::summarize(butyrate.clr = sum(numseqs))

df <- left_join(df, butyrate.clr)
rm(butyrate.clr)

comparisons <- list(c("control day 1", "case day 1"), c("control day 1", "case day 14"), c("control day 1", "case 6 months"))
ggplot(df, aes(x= group_time, y = butyrate.clr, fill = group_time)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "black", pch = 21, alpha =.75, size = 3.5, width = 0.2)+
  scale_color_manual(values=c("#91c3e1", "#660033", "#FC8D59", "#FEE08B")) +
  scale_fill_manual(values=c("#91c3e1", "#660033", "#FC8D59", "#FEE08B")) +
  stat_compare_means(comparisons = comparisons) +
  theme_cowplot(font_size = 12)+
  xlab("")+
  ylab("Butyrate-producing bacteria")+
  theme(legend.position = "none")
```

Alluvial plots were used to visualise microbiota differences between the different timepoints:
```
P.phylum <- microbiome::aggregate_taxa(phyloseq, level = "Phylum")

most_abundant_phyla <- as.data.frame(sort(taxa_sums(P.phylum), TRUE)[1:5])
most_abundant_phyla$otu <- row.names(most_abundant_phyla)
composition_phylum <- get.otu.melt(P.phylum, filter.zero =F)  %>%
  mutate(Phylum = if_else(otu %in% most_abundant_phyla$otu, Phylum, "Others")) %>%
  mutate(Phylum = fct_relevel(Phylum, "p__Actinobacteriota", "p__Bacteroidota", "p__Firmicutes", "p__Firmicutes_A", "p__Proteobacteria", "Others")) 
comp_time_phy <- composition_phylum %>%
  group_by(group_time, Phylum) %>%
  summarise(pctseqs = mean(pctseqs)) %>%
  mutate(Phylum = fct_relevel(Phylum, "p__Actinobacteriota", "p__Bacteroidota", "p__Firmicutes", "p__Firmicutes_A", "p__Proteobacteria", "Others")) 

comp_time_phylum <- comp_time_phy %>%
  group_by(group_time) %>%
  mutate(pctseqs = pctseqs* 100/sum(pctseqs)) %>%
  ggplot(aes(y = pctseqs, x = group_time, fill = Phylum)) +
  geom_flow(aes(alluvium = Phylum), alpha= .4, color = "white",
            curve_type = "linear", 
            width = .5) +
  geom_col(width = .5, color = "white") +
  scale_fill_manual(values=c("#FFC107", "#004D40", '#88CCEE', "#003366", "#cc6677", "#808080")) +
  scale_y_continuous(NULL, expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw()+
  theme(axis.title.x=element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "right") 
```
Similar code was used for the alluvial plot at genus level. 

MaAsLin2 was used for differential abundance analysis:
```
ps.da <- microbiome::aggregate_taxa(phyloseq, level = "Species") # analysis is performed at species level
ps.da <- core(ps.da, detection=1, prevalence=10/100, include.lowest=T) # analysis was limited to species present in at least 10 percent of the samples

features <- t(get.otu(ps.da))
metadata <- get.samp(ps.da) %>%
  dplyr::select(group, group_time, sample) %>%
  column_to_rownames(var = "sample")

maaslin <- Maaslin2(
  input_data = features, 
  input_metadata = metadata, 
  output = "maaslin_output",
  fixed_effects = "group_time", 
  reference = "group_time,control day 1")

control.d1 <- as.data.frame(maaslin[["results"]])  %>%
  filter(value == "case day 1") %>%
  filter(qval < 0.05) %>%
  mutate(group=ifelse(coef < 0, "Higher in controls", "Higher in cases"))

ggplot(control.d1, aes(x=reorder(feature, coef), y=coef, fill=group)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  ylab("log 2-Fold change") +
  xlab("") + 
  theme_bw(base_size = 11) +
  scale_fill_manual(values=c("#660033", "#91c3e1")) +
  scale_y_continuous(limits = c(-2.45, 3.4)) + 
  theme(legend.position = "none")
```
Similar code was used to compare other timepoints (day 14, 6 months) to controls. 


## 4 - Resistome analysis (figure 2)

Compare the total number of detected fecal antimicrobial resistance (AMR)-associated genes:
```
amr.count <- amr %>%
  dplyr::group_by(sample, group, group_time, reads) %>%
  dplyr::summarise(count = length(GENE)) %>%
  ungroup() %>%
  mutate(count_corrected = (count/reads)*1000000) # counts per million reads

comparisons <- list(c("control day 1", "case day 1"), c("control day 1", "case day 14"), c("control day 1", "case 6 months"))

ggplot(amr.count, aes(x= group_time, y = count_corrected, fill = group_time)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "black", pch = 21, alpha =.75, size = 3.5, width = 0.2)+
  scale_color_manual(values=c("#91c3e1", "#660033", "#FC8D59", "#FEE08B")) +
  scale_fill_manual(values=c("#91c3e1", "#660033", "#FC8D59", "#FEE08B")) +
  stat_compare_means(comparisons = comparisons) +
  theme_cowplot(font_size = 12)+
  scale_y_continuous(trans=scales::pseudo_log_trans(), breaks = c(0, 1, 5, 10, 20)) +
  xlab("")+
  ylab("Number of AMR genes")+
  theme(legend.position = "none")
```

PCoA was used to show overall resistome profiles
```
amr.pca <- amr %>%
  dplyr::select(sample, GENE, value) %>%
  reshape2::dcast(sample ~ GENE, value.var = "value",  sum) %>%
  column_to_rownames(var = "sample")
amr.pca <- as.matrix(amr.pca)
amr.otu <- otu_table(amr.pca, taxa_are_rows = F)

bray <- phyloseq::distance(amr.otu, method = "bray") 
df.bray <- as.matrix(bray)
df.bray <- as.data.frame(row.names(df.bray)) %>%
  mutate(sample = `row.names(df.bray)`) %>%
  left_join(df)

ord <- ordinate(amr.otu, method = "PCoA", distance = bray)

bray.plot <- plot_ordination(phyloseq, ord, type = "samples", color = "group_time", justDF = T) %>%
  rownames_to_column(var = "sample") %>%
  dplyr::select(Axis.1, Axis.2, sample, group_time)

centroids <- aggregate(cbind(Axis.1, Axis.2)~group_time, data=bray.plot, mean) # calculate the centroids per group
bray.plot <- merge(bray.plot, centroids, by="group_time", suffixes=c("", ".centroid")) # merge the centroid data with PCoA data

ggplot(bray.plot, aes(Axis.1, Axis.2, color = group_time)) +
  geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2, color= group_time), alpha = 0.5)+ 
  geom_point(data=bray.plot, aes(color=group_time),size=4.5,alpha=1) + # Create points for each participant
  geom_label_repel(data = centroids, aes(x=Axis.1, y=Axis.2, label=c("Control", "Melioidosis - diagnosis", "Melioidosis - day 14", "Melioidosis - 6 months")), size=4.5, fill = ggplot2::alpha(c("white"),0.92)) +
  theme_bw() +
  scale_color_manual(values=c("#91c3e1", "#660033", "#FC8D59", "#FEE08B")) +
  theme(legend.position = "none")
```

Bar plot showing the AMR genes per class of antibiotics:
```
colours <- c("#b1004e", #Aminoglycosides 
            "#00488c", #BETA-LACTAM
            "#ffff4d", #CARBAPENEM
            "#8b6500", # CHLORAMPHENICOL
            "#BF40BF", #MACROLIDES 
            "#bcbcbc", #OTHER
            "#ffa293", #QUINOLONE
            "#00f8fb", #STREPTOMYCIN
            "#69932d",  #SULFONAMIDE
            "#00959f", #TETRACYCLINES 
            "#00c06d") #TRIMETHOPRIM
            
amrcompplot <- amr[amr$RESISTANCE %in% names(sort(table(amr$RESISTANCE), decreasing = TRUE)[1:10]), ]
amrcompplot <- amr %>%
  mutate(resistance_filtered = if_else(RESISTANCE %in% selection, RESISTANCE, "Other")) %>%
  mutate(resistance_filtered = as.factor(resistance_filtered)) %>%
  mutate(group_time = as.character(group_time)) %>%
  dplyr::group_by(sample, group_time, reads, resistance_filtered, .drop=F) %>%
  dplyr::summarise(count = length(GENE)) %>%
  ungroup() %>%
  mutate(count_corrected = (count/reads)*1000000) %>% # counts per million reads
  mutate(group_time = as.factor(group_time)) %>%
  mutate(group_time = fct_relevel(group_time, "control day 1", "case day 1", "case day 14", "case 6 months")) 

amrcompplot %>%
  ggplot(aes(x=reorder(sample, -count), y=count, fill=resistance_filtered)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=colours) +
  theme_bw()+
  ylab("Number of AMR genes") +
  facet_wrap(~group_time, scales = "free_x", nrow = 1) +
  scale_y_continuous(expand = c(0,0,0,1),) +
  theme(axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")
```

Finally, we compared the number of AMR genes per class of antibiotics between timepoints:
```
selection <- c("AMINOGLYCOSIDE", "BETA-LACTAM", 
               "CARBAPENEM", "CHLORAMPHENICOL", "MACROLIDE", 
               "QUINOLONE",  "STREPTOMYCIN", "SULFONAMIDE",
               "TETRACYCLINE",  "TRIMETHOPRIM") # Classes that were most often present were selected
comparisons <- list(c("control day 1", "case day 1"), c("control day 1", "case day 14"), c("control day 1", "case 6 months"))
amr.abx <- amr %>%
  filter(RESISTANCE %in% selection) %>%
  mutate(RESISTANCE = as.factor(RESISTANCE)) %>%
  mutate(group_time = as.character(group_time)) %>%
  dplyr::group_by(sample, group_time, reads, RESISTANCE, .drop=F) %>%
  dplyr::summarise(count = length(GENE)) %>%
  ungroup() %>%
  mutate(count_corrected = (count/reads)*1000000) %>% # counts per million reads
  mutate(group_time = as.factor(group_time)) %>%
  mutate(group_time = fct_relevel(group_time, "control day 1", "case day 1", "case day 14", "case 6 months")) 

amr.abx %>%
  ggplot(aes(x=group_time, y=count_corrected, fill=group_time)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "black", pch = 21, alpha =.75, size = 3.5, width = 0.2)+
  scale_fill_manual(values=c("#91c3e1", "#660033", "#FC8D59", "#FEE08B")) +
  theme_bw()+
  facet_wrap(~RESISTANCE, nrow = 2) +
  ylab("Number of AMR genes") +
  stat_compare_means(comparisons = comparisons) +
  theme(axis.title.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
```

## 5 - Linking antibiotic exposure to gut microbiota composition and resistome (figure 3)
First, we assessed the contribution of clinical variables to dissimilarities in microbiota composition at 6 months:
```
ps.abxm6 <- phyloseq
df.abxm6 <- df %>%
  filter(time == "6")
sample_data(ps.abxm6) <- set.samp(df.abxm6)

set.seed(88) 
bray <- phyloseq::distance(ps.abxm6, method = "bray") 
df.bray <- as.matrix(bray)
df.bray <- as.data.frame(row.names(df.bray)) %>%
  mutate(sample = `row.names(df.bray)`) %>%
  left_join(df.abxm6)

variables <-c("Age", "Sex", "BMI", "Any comorbidity", "Diabetes", "Hypertension",
              "Chronic kidney disease", "Chronic liver disaease", "Neutrophils", "Lymphocytes", 
              "SOFA at admission", "ICU admission", "Length of hospital stay", "Disseminated/Localised", 
              "Cephalosporins", "Meropenem", "Penicillins with B-lac inhib", "Co-trimoxazol", "Number of antibiotics", "Days of Therapy")

set.seed(88)
R2 <- list((bray ~ Age),
           (bray ~ Sex),
           (bray ~ BMI),
           (bray ~ Comorbidities),
           (bray ~ Diabetes),
           (bray ~ Hypertension),
           (bray ~ `Chronic kidney disease`), 
           (bray ~ `Chronic liver disease`), 
           (bray ~ Neutrophil),
           (bray ~ Lymphocyte),
           (bray ~ sofa_total), 
           (bray ~ `ICU stayed at any time point`),
           (bray ~ `Hospital stay (days)`),
           (bray ~ disseminated), 
           (bray ~ cephalosporins_pre_sample3), 
           (bray ~ meropenem_pre_sample3), 
           (bray ~ penicillin_pre_sample3), 
           (bray ~ trimethoprim_pre_sample3), 
           (bray ~ no_antibiotics_pre_sample3), 
           (bray ~ days_of_therapy_pre_sample3)) %>%
  lapply(adonis2, data = df.bray, by = 'margin', permutations = 9999) %>%
  lapply(function(x) x[["R2"]][[1]]) %>%
  do.call(c, .)

Pvalues <- list((bray ~ Age),
                (bray ~ Sex),
                (bray ~ BMI),
                (bray ~ Comorbidities),
                (bray ~ Diabetes),
                (bray ~ Hypertension),
                (bray ~ `Chronic kidney disease`), 
                (bray ~ `Chronic liver disease`), 
                (bray ~ Neutrophil),
                (bray ~ Lymphocyte),
                (bray ~ sofa_total), 
                (bray ~ `ICU stayed at any time point`),
                (bray ~ `Hospital stay (days)`),
                (bray ~ disseminated), 
                (bray ~ cephalosporins_pre_sample3), 
                (bray ~ meropenem_pre_sample3), 
                (bray ~ penicillin_pre_sample3), 
                (bray ~ trimethoprim_pre_sample3), 
                (bray ~ no_antibiotics_pre_sample3), 
                (bray ~ days_of_therapy_pre_sample3)) %>%
  lapply(adonis2, data = df.bray, by = 'margin', permutations = 9999) %>%
  lapply(function(x) x[["Pr(>F)"]][[1]]) %>%
  do.call(c, .)

dissim <- as.data.frame(cbind(variables, R2, Pvalues))
dissim <- dissim %>%
  mutate(R2 = as.numeric(as.character(R2))) %>%
  mutate(Pvalues = as.numeric(as.character(Pvalues)))

ggplot(dissim, aes(x=R2, y=reorder(variables, -R2))) + 
  geom_bar(stat = "identity", fill = "#FEE08B") +
  #geom_text(aes(x = -.005, y = reorder(variables, -R2), label = paste("p=", round(Pvalues,3)))) +
  scale_x_continuous(labels = scales::percent_format()) +
  xlab("Contribution to interindividual dissimilarities in microbiota (R2)") +
  ggtitle("Melioidosis - 6 months") +
  ylab("") +
  theme_bw(base_size = 10)

rm(dissim, Pvalues, R2, bray, df.bray)
```

Next, we assessed changes in microbiota associated with antibiotic exposure: 
```
features <- t(get.otu(ps.abxm6))
metadata <- get.samp(ps.abxm6) %>%
  dplyr::select(cephalosporins_pre_sample3, trimethoprim_pre_sample3, days_of_therapy_pre_sample3, `Hospital.stay..days.`, 
                sample) %>%
  column_to_rownames(var = "sample")

cephalosporins_pre_sample3 <- as.data.frame(Maaslin2(
  input_data = features, 
  input_metadata = metadata, 
  output = "maaslin_output",
  fixed_effects = "cephalosporins_pre_sample3")[["results"]])

trimethoprim_pre_sample3 <- as.data.frame(Maaslin2(
  input_data = features, 
  input_metadata = metadata, 
  output = "maaslin_output",
  fixed_effects = "trimethoprim_pre_sample3")[["results"]])

dot_pre_sample3 <- as.data.frame(Maaslin2(
  input_data = features, 
  input_metadata = metadata, 
  output = "maaslin_output",
  fixed_effects = "days_of_therapy_pre_sample3")[["results"]]) 

los_pre_sample3 <- as.data.frame(Maaslin2(
  input_data = features, 
  input_metadata = metadata, 
  output = "maaslin_output",
  fixed_effects = "Hospital.stay..days.")[["results"]]) 

res_pre_samp3 <- rbind(cephalosporins_pre_sample3, trimethoprim_pre_sample3, dot_pre_sample3, los_pre_sample3) %>% 
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  mutate(group = if_else(coef > 0, "Higher", "Lower")) %>%
  mutate(metadata = fct_relevel(metadata, 'trimethoprim_pre_sample3', 'cephalosporins_pre_sample3',
                                'days_of_therapy_pre_sample3', 'Hospital.stay..days.')) %>%
  filter(qval < 0.05)

ggplot(res_pre_samp3, aes(x=reorder(feature, coef), y=coef, fill=group)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  ylab("Coefficient") +
  xlab("") + 
  theme_bw(base_size = 12) +
  facet_wrap(~metadata, scales = "free_y") +
  ylim(-5.7,2) +
  scale_fill_manual(values=c("#882255", "#FC8D59")) +
  theme(legend.position = "none")
```

Next, we assessed associations between the length of hospital stay, days of antibiotic therapy or antiboitics and gut resistome profiles. Here, the code for the associations with length of hospital stay and with meropenem are provided:
```
vars <- c("Aminoglycoside", "Beta-lactam", "Carbapenem", "Chloramphenicol",
          "Macrolide", "Quinolone", "Streptomycin", "Sulfonamide", "Tetracycline", "Trimethoprim")
estimate <- rbind(
  cor.test(amr.abx$AMINOGLYCOSIDE, amr.abx$hospstay, data = amr.abx), 
  cor.test(amr.abx$`BETA-LACTAM`, amr.abx$hospstay, data = amr.abx),
  cor.test(amr.abx$CARBAPENEM, amr.abx$hospstay, data = amr.abx),
  cor.test(amr.abx$CHLORAMPHENICOL, amr.abx$hospstay, data = amr.abx),
  cor.test(amr.abx$MACROLIDE, amr.abx$hospstay, data = amr.abx),
  cor.test(amr.abx$QUINOLONE, amr.abx$hospstay, data = amr.abx),
  cor.test(amr.abx$STREPTOMYCIN, amr.abx$hospstay, data = amr.abx),
  cor.test(amr.abx$SULFONAMIDE, amr.abx$hospstay, data = amr.abx),
  cor.test(amr.abx$TETRACYCLINE, amr.abx$hospstay, data = amr.abx),
  cor.test(amr.abx$TRIMETHOPRIM, amr.abx$hospstay, data = amr.abx)) %>%
  data.frame() %>%
  mutate(estimate = as.numeric(estimate)) %>%
  cbind(vars) %>%
  dplyr::select(vars, estimate) %>%
  mutate(hospstay_effsize = ifelse(estimate < 0.1, 0, estimate),
         effsize_cat = cut(estimate,
                           breaks = c(-1, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 1)) %>% as.factor()) %>%
  select(vars, hospstay_effsize) %>%
  column_to_rownames("vars") %>%
  as.matrix()

col_fun = circlize::colorRamp2(breaks=c(-.6, 0, .6), c( "#003366", "#EEEEEE", "#660033"))

Heatmap(estimate,
        col = col_fun,
        show_column_names = T, 
        show_row_names = T, 
        row_names_side = "left",
        cluster_rows = F,
        cluster_columns = F,
        border=T,
        cluster_row_slices = F,
        heatmap_legend_param = list(at = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5),
                                    labels = c("<-0.5", "<-0.3", "<-0.1", ">-0.1 & <0.1", ">0.1", ">0.3", ">0.5"),
                                    title = "Pearson r", color_bar = "discrete", fontsize = 8),
        rect_gp = gpar(col = "black", lwd = 0.5),
        border_gp = gpar(col = "black", lwd = 1))
```

```
estimate <- rbind(
  cohens_d(AMINOGLYCOSIDE ~ meropenem_pre_sample, data = amr.abx, hedges.correction = T, var.equal = F), 
  cohens_d(`BETA-LACTAM` ~ meropenem_pre_sample, data = amr.abx, hedges.correction = T, var.equal = F), 
  cohens_d(CARBAPENEM ~ meropenem_pre_sample, data = amr.abx, hedges.correction = T, var.equal = F), 
  cohens_d(CHLORAMPHENICOL ~ meropenem_pre_sample, data = amr.abx, hedges.correction = T, var.equal = F), 
  cohens_d(MACROLIDE ~ meropenem_pre_sample, data = amr.abx, hedges.correction = T, var.equal = F), 
  cohens_d(QUINOLONE ~ meropenem_pre_sample, data = amr.abx, hedges.correction = T, var.equal = F), 
  cohens_d(STREPTOMYCIN ~ meropenem_pre_sample, data = amr.abx, hedges.correction = T, var.equal = F), 
  cohens_d(SULFONAMIDE ~ meropenem_pre_sample, data = amr.abx, hedges.correction = T, var.equal = F), 
  cohens_d(TETRACYCLINE ~ meropenem_pre_sample, data = amr.abx, hedges.correction = T, var.equal = F), 
  cohens_d(TRIMETHOPRIM ~ meropenem_pre_sample, data = amr.abx, hedges.correction = T, var.equal = F)) %>%
  mutate(effsize = -effsize) %>% # transform so usage of abx gives higher value to aid interpretation
  mutate(meropenem_effsize = ifelse(effsize > -0.2 & effsize < 0.2, 0, effsize),
         effsize_cat = cut(effsize,
                           breaks = c(-Inf, -0.8, -0.5, -0.2, 0.2, 0.5, 0.8, Inf)) %>% as.factor()) %>%
  select(`.y.`, meropenem_effsize) %>%
  column_to_rownames(".y.") %>%
  as.matrix()
  
col_fun = circlize::colorRamp2(breaks=c(-1, 0, 1), c( "#003366", "#EEEEEE", "#660033"))

Heatmap(estimate,
        col = col_fun,
        show_column_names = T, 
        show_row_names = T, 
        row_names_side = "left",
        cluster_rows = F,
        cluster_columns = F,
        border=T,
        cluster_row_slices = F,
        heatmap_legend_param = list(at = c(-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8),
                                    labels = c("<-0.8", "<-0.5", "<-0.2", ">-0.2 & <0.2", ">0.2", ">0.5", ">0.8"),
                                    title = "Hedges' g", color_bar = "discrete", fontsize = 8),
        rect_gp = gpar(col = "black", lwd = 0.5),
        border_gp = gpar(col = "black", lwd = 1)) #4x4
```

```
amr.merop <- amr %>%
  filter(group_time != "control day 1") %>%
  dplyr::select(sample, GENE, meropenem_pre_sample) %>%
  mutate(value = 1)  %>%
  reshape2::dcast(sample + meropenem_pre_sample ~ GENE, value.var = "value")

res.merop = data.frame(variable = "x",
                      pvalue=2, padj=2, 
                      presence_users = 999, 
                      presence_nonusers = 999,
                      stringsAsFactors = F)

for(i in 3:202) {
  amr.merop[,i] = if_else(amr.merop[,i] > 1, 1, amr.merop[,i])
  pvalue = fisher.test(table(amr.merop$meropenem_pre_sample, amr.merop[,i]))[["p.value"]]
  res.merop[i, "variable"] = colnames(amr.merop[i])
  res.merop[i, "pvalue"] = pvalue
  res.merop[i, "presence_users"] = count(amr.merop[amr.merop$meropenem_pre_sample==1,][,i]==1)
  res.merop[i, "presence_nonusers"] = count(amr.merop[amr.merop$meropenem_pre_sample==0,][,i]==1)
}

res.merop <- res.merop %>%
  filter(!is.na(variable)) %>%
  filter(variable != "x") %>%
  mutate(presence_users_relative = (presence_users/32)*100) %>%
  mutate(presence_nonusers_relative = (presence_nonusers/38)*100) %>%
  mutate(percentage_change = presence_users_relative-presence_nonusers_relative) %>%
  mutate(label = if_else(pvalue<0.05, variable, "")) %>%
  mutate(group = if_else(percentage_change > 0, "Higher in users", "Higher in nonusers"))

ggplot(res.merop) +
  geom_point(aes(x = percentage_change, y = -log10(pvalue), fill=group),shape=21, size = 5, alpha = 1, colour="black") +
  geom_text_repel(aes(x = percentage_change, y = -log10(pvalue), label = label), min.segment.length = unit(0, 'lines'),
                  nudge_x = 0, nudge_y=0.1, segment.alpha = 0.3, force=2) +
  xlab("Change in proportion of participants with AMR gene") +
  ylab("-log10 p-value") +
  theme_classic() +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  scale_fill_manual(values=c("#003366", "#660033")) +
  theme(legend.position = "none")
```

## 6 - Microbiota and nosocomial infections with Klebsiella pneumoniae (figure 4)

Given the higher abundances of Klebsiella pneumoniae during melioidosis and because K. pneumoniae was the most common cause of secondary infection (Supplementary Figure 3), we explored the relationship between fecal K. pneumoniae and secondary K. pneumoniae infections:

```
df.c1 <- df %>%
  filter(group_time == "case day 1" | group_time == "case day 14") # select samples during hospital stay

phyloseq.c1 <- phyloseq
sample_data(phyloseq.c1) <- set.samp(df.c1)

set.seed(88) 
bray <- phyloseq::distance(phyloseq.c1, method = "bray") 
df.bray <- as.matrix(bray)
df.bray <- as.data.frame(row.names(df.bray)) %>%
  mutate(sample = `row.names(df.bray)`) %>%
  left_join(df.c1)
adonis2((bray ~ `sec_infection_klebsiella`), by = 'margin', data = df.bray, permutations =9999) 

ord <- ordinate(phyloseq.c1, method = "PCoA", distance = bray)
plot_ordination(phyloseq.c1, ord, type="samples", color="sec_infection_klebsiella")
bray.plot <- plot_ordination(phyloseq.c1, ord, type = "samples", color = "sec_infection_klebsiella", justDF = T) %>% #
  rownames_to_column(var = "sample") %>%
  left_join(df.c1)
centroids <- aggregate(cbind(Axis.1, Axis.2)~`sec_infection_klebsiella`, data=bray.plot, mean) # calculate the centroids per group
bray.plot <- merge(bray.plot, centroids, by="sec_infection_klebsiella", suffixes=c("", ".centroid")) # merge the centroid data with PCoA data

ggplot(bray.plot, aes(Axis.1, Axis.2, color = `sec_infection_klebsiella`)) +
  geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2, color= `sec_infection_klebsiella`), alpha = 0.3)+ 
  geom_point(data=bray.plot, aes(color=`sec_infection_klebsiella`, shape = group_time),size=3.5,alpha=1.0) + # Create points for each participant
  geom_label_repel(data = centroids, aes(x=Axis.1, y=Axis.2, label=c("No K. pneumoniae infection", "K. pneumoniae infection")), size=5, fill = ggplot2::alpha(c("white"),0.76)) +
  theme_bw() +
  xlab("22.8% variance") + #Label X-axis
  ylab("19.3% variance") + #Label Y-axis
  scale_color_manual(values=c("#004D40", "#FFC107")) +
  theme(legend.position = "bottom")
```

```
klebs <- get.otu.melt(phyloseq.c1, filter.zero = F) %>%
  filter(Species == "s__Klebsiella pneumoniae") %>%
  dplyr::group_by(sample, reads, group_time,day_of_diagnose, day_sec_sample, sec_infection_klebsiella, sec_infection_klebsiella_day) %>%
  summarise(sum = sum(pctseqs)) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = sampleday_relativeklebs, y = sum)) +
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "grey36", size = 0.4) + 
  geom_point() +
  theme_bw() +
  geom_smooth() +
  xlab("Day relative to first positive culture with K. pneumoniae") +
  ylab("Klebsiella pneumoniae (%)") +
  scale_y_continuous(limits = c(-0.15,1.1))
```
```
klebs <- get.otu.melt(phyloseq, filter.zero = F) %>%
  filter(Species == "s__Klebsiella pneumoniae") %>%
  dplyr::group_by(sample, group_time, sec_infection_klebsiella) %>%
  summarise(sum = sum(pctseqs)) %>%
  dplyr::ungroup() 

klebs %>%
  mutate(sec_infection_klebsiella_controls = if_else(group_time == "control day 1", "control", sec_infection_klebsiella)) %>%
  filter(group_time != "case 6 months") %>%
  mutate(sec_infection_klebsiella_controls = fct_relevel(sec_infection_klebsiella_controls, "control", "0", "1")) %>%
  ggplot(aes(x= sec_infection_klebsiella_controls, y = sum, fill = sec_infection_klebsiella_controls)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(shape = group_time), color = "black", alpha =.75, size = 2, width = 0.2)+
  stat_compare_means(comparisons = list(c("0", "1"), c("1", "control"), c("0", "control"))) +
  theme_cowplot()+
  xlab("")+
  ylab("Klebsiella pneumoniae (%)")+
  scale_fill_manual(values=c("#000000", "#004D40", "#FFC107")) +
  theme(legend.position = "none")
```

