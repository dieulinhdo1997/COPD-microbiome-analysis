#load BLAST data
blast <- read.delim("COPD_sherlock1_result.txt", header = TRUE, stringsAsFactors = FALSE)
#check data structure
head(blast)
#add column names
colnames(blast) <- c("query ID", "Organism", "% identity", "alignment length", "mismatch", "gap opens", "q.start", "q.end", "s,start", "s.end", "E_Value", "Bit_score")

install.packages("tidyverse")
library(tidyverse)
library(dplyr)

min_identity <- 95
filtered_microbes <- blast %>% 
  filter(`% identity` >= min_identity) %>% #Filter for identity >= 95
  distinct() #remove duplicates
print(filtered_microbes)

#select the row with the highest Bit score for each query ID
best_hits <- filtered_microbes %>%
  group_by(`query ID`) %>%
  slice_max(`Bit_score`, with_ties = TRUE) %>%
  ungroup()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager:: install("Biostrings")
library(Biostrings)

database_file <- readDNAStringSet("C:\\Users\\Admin\\Downloads\\all_seqs.fa\\all_seqs.fa")
headers <- names(database_file)
taxonomy_db <- data.frame(organism_code = sub(" .*", "", headers),  # Extract bacterial code (before first space)
                          organism_name = sub("^[^ ]+ ", "", headers),  # Extract organism name (after first space)
                          stringsAsFactors = FALSE)
taxonomy_db <- taxonomy_db %>% rename("organism_code" = "Organism")
write.table(taxonomy_db, "taxonomy.txt", sep = "\t", row.names = FALSE, quote = FALSE)
rm(database_file)

best_hits_with_taxa <- best_hits %>%
  left_join(taxonomy_db, by = "Organism")

write.table(best_hits_with_taxa, "best_hits_with_taxa.txt", sep = "\t", row.names = FALSE, quote = FALSE)

best_hits_with_taxa <- best_hits_with_taxa %>%
  mutate(organism_name = as.character(organism_name))

best_hits_final <- best_hits_with_taxa %>%
  group_by(`query ID`) %>%
  arrange(organism_name) %>%  # Sort alphabetically
  filter(row_number() == 1) %>%  # Select the first row
  ungroup()

filtered_organisms <- best_hits_final %>%
  select(Organism, organism_name) %>%
  distinct()

filtered_organisms_1 <- filtered_organisms %>%
  mutate(organism_name = str_extract(organism_name, "^\\S+\\s\\S+"))

filtered_organisms_1 <- filtered_organisms_1 %>%
  distinct(organism_name, .keep_all = TRUE)
write.table(filtered_organisms_1, "reference_bacteria_sherlock1.txt", sep = "\t", row.names = FALSE, quote = FALSE)

best_hits_final <- best_hits_final %>% 
  mutate(organism_name_2 = str_extract(organism_name, "^\\S+\\s\\S+"))

write.table(best_hits_final, "best_hits_final_sherlock1.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Count the number of contigs for each organism
contig_counts <- best_hits_final %>% 
  distinct(`query ID`, organism_name_2) %>% #ensure uniqueness
  group_by(organism_name_2) %>%
  summarise(contig_count = n(), .groups = "drop")
write.table(contig_counts, "contig_counts_sherlock1.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Count the number of bases for each organism
fasta_file <- readDNAStringSet("C:\\Users\\Admin\\Downloads\\COPD\\20250401_sherlock1_BB_trinity.Trinity.fasta")
names(fasta_file) <- sub(" .*", "", names(fasta_file))
fasta_df <- data.frame(sequence_id = names(fasta_file), #extract sequence IDs
                       sequence_length = width(fasta_file) #calculate the length of eachsequence
)
merged_data <- best_hits_final %>%
  inner_join(fasta_df, by = c("query ID" = "sequence_id"))
base_counts <- merged_data %>%
  group_by(organism_name_2) %>%
  summarise(total_bases = sum(sequence_length)) %>%
  ungroup()
write.table(base_counts, "base_counts_sherlock1.txt", sep = "\t", row.names = FALSE, quote = FALSE)

install.packages("ggplot2")
install.packages("pheatmap")
library(ggplot2)
library(pheatmap)
#load bacterial read count data (cpm file) and metadata
read_count_table <- read.delim("20250414_Sherlock1_bronchial_brushes_bact_exp.txt", row.names=1, check.names= FALSE)
metadata <- read.delim("20250414_Sherlock1_bronchial_brushes_samples.txt", row.names = 1, check.names =  FALSE)
row.names(read_count_table) <- read_count_table[,1]
read_count_table_1 <- read_count_table[,-1]

#differential analysis (edgeR)
library(ggplot2)
library(ggpubr) 
library(ggrepel)
#differential analysis for disease group
library(edgeR)
group_diseasegroup <- factor(metadata_diseasegroup$Disease, levels = c("Control", "Moderate COPD", "Severe COPD")) #Set up group variable
dge_diseasegroup <- DGEList(counts = raw_count_table_diseasegroup, group = group_diseasegroup) #Create DGEList
dge_diseasegroup <- calcNormFactors(dge_diseasegroup) #TMM normalization
design_diseasegroup <- model.matrix(~ 0 + group_diseasegroup) # Design matrix
colnames(design_diseasegroup) <- make.names(colnames(design_diseasegroup))
colnames(design_diseasegroup)
dge_diseasegroup <- estimateDisp(dge_diseasegroup, design_diseasegroup) # Estimate dispersions
fit_diseasegroup <- glmFit(dge_diseasegroup, design_diseasegroup) # Fit GLM model
#Create contrasts for pairwise comparisons
contrast_matrix <- makeContrasts(
  Mod_vs_Ctrl = group_diseasegroupModerate.COPD - group_diseasegroupControl,
  Sev_vs_Ctrl = group_diseasegroupSevere.COPD - group_diseasegroupControl,
  Sev_vs_Mod  = group_diseasegroupSevere.COPD - group_diseasegroupModerate.COPD,
  levels = design_diseasegroup
) # Define contrast matrix using valid names
#Run glmLRT for each pair
lrt_mod_ctrl <- glmLRT(fit_diseasegroup, contrast = contrast_matrix[,"Mod_vs_Ctrl"])
lrt_sev_ctrl <- glmLRT(fit_diseasegroup, contrast = contrast_matrix[,"Sev_vs_Ctrl"])
lrt_sev_mod  <- glmLRT(fit_diseasegroup, contrast = contrast_matrix[,"Sev_vs_Mod"])
#Extract result tables
res_mod_ctrl <- topTags(lrt_mod_ctrl, n = Inf)$table
res_mod_ctrl$FDR <- p.adjust(res_mod_ctrl$PValue, method = "fdr")

res_sev_ctrl <- topTags(lrt_sev_ctrl, n = Inf)$table
res_sev_ctrl$FDR <- p.adjust(res_sev_ctrl$PValue, method = "fdr")

res_sev_mod <- topTags(lrt_sev_mod, n = Inf)$table
res_sev_mod$FDR <- p.adjust(res_sev_mod$PValue, method = "fdr")
#View number of significant taxa in each pairwise comparison
cat("Significant taxa (Moderate vs Control):", sum(res_mod_ctrl$FDR < 0.05), "\n")
cat("Significant taxa (Severe vs Control):", sum(res_sev_ctrl$FDR < 0.05), "\n")
cat("Significant taxa (Severe vs Moderate):", sum(res_sev_mod$FDR < 0.05), "\n")

sig_taxa_mod_ctrl <- rownames(res_mod_ctrl)[res_mod_ctrl$FDR < 0.05]
sig_taxa_sev_ctrl <- rownames(res_sev_ctrl)[res_sev_ctrl$FDR < 0.05]
sig_taxa_sev_mod  <- rownames(res_sev_mod)[res_sev_mod$FDR < 0.05]

res_mod_ctrl$Species <- rownames(res_mod_ctrl)
res_sev_ctrl$Species <- rownames(res_sev_ctrl)
res_sev_mod$Species  <- rownames(res_sev_mod)

res_mod_ctrl_sig <- res_mod_ctrl %>% filter(FDR<0.05)
res_sev_ctrl_sig <- res_sev_ctrl %>% filter(FDR<0.05)
res_sev_mod_sig <- res_sev_mod %>% filter(FDR<0.05)
#top 50 moderate vs control
top50_mod_ctrl <- res_mod_ctrl %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
ggplot(top50_mod_ctrl, aes(x = reorder(Species, logFC), y = logFC, fill = logFC > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "pink", "FALSE" = "green"),
    labels = c("FALSE" = "Control", "TRUE" = "Moderate COPD"),
    name = "Higher in"
  ) +
  labs(title = "Top 50 Differentially Abundant Species (Moderate vs Control)",
       x = "Species",
       y = "log2 Fold Change (Moderate - Control)") +
  theme_minimal()
#volcano plot
res_mod_ctrl <- res_mod_ctrl %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_mod_ctrl, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_mod_ctrl %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: Moderate vs Control",
    x = "log2 Fold Change (Moderate - Control)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#top 50 severe vs control
top50_sev_ctrl <- res_sev_ctrl %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
ggplot(top50_sev_ctrl, aes(x = reorder(Species, logFC), y = logFC, fill = logFC > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "blue", "FALSE" = "green"),
    labels = c("FALSE" = "Control", "TRUE" = "Severe COPD"),
    name = "Higher in"
  ) +
  labs(title = "Top 50 Differentially Abundant Species (Severe vs Control)",
       x = "Species",
       y = "log2 Fold Change (Severe - Control)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6),  # 👈 adjust this size as needed
  )
#volcano plot
res_sev_ctrl <- res_sev_ctrl %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_sev_ctrl, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_sev_ctrl %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: Severe vs Control",
    x = "log2 Fold Change (Severe - Control)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#top 50 severe vs moderate
top50_sev_mod <- res_sev_mod %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
ggplot(top50_sev_mod, aes(x = reorder(Species, logFC), y = logFC, fill = logFC > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "blue", "FALSE" = "pink"),
    labels = c("FALSE" = "Moderate COPD", "TRUE" = "Severe COPD"),
    name = "Higher in"
  ) +
  labs(title = "Top 50 Differentially Abundant Species (Severe vs Modrate)",
       x = "Species",
       y = "log2 Fold Change (Severe - Moderate)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6), 
  )
#volcano plot
res_sev_mod <- res_sev_mod %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_sev_mod, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_sev_mod %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: Severe vs Moderate",
    x = "log2 Fold Change (Severe - Moderate)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
ggplot(res_sev_mod, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7) +
  geom_text_repel(
    data = res_sev_mod %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  labs(title = "Volcano Plot: Severe vs Moderate",
       x = "log2 Fold Change (Severe - Moderate)",
       y = "-log10(FDR)",
       color = "Significant (FDR < 0.05)") +
  theme_minimal()
#heatmap
lrt_diseasegroup <- glmLRT(fit_diseasegroup)
res_diseasegroup <- topTags(lrt_diseasegroup, n = Inf)$table
res_diseasegroup$FDR <- p.adjust(res_diseasegroup$PValue, method = "fdr")
res_diseasegroup$Species <- rownames(res_diseasegroup)
res_diseasegroup_sig <- res_diseasegroup %>% filter(FDR<0.05)
cpm_diseasegroup <- cpm(dge_diseasegroup, log = FALSE)
sig_species_diseasegroup <- res_diseasegroup_sig$Species
print(sig_species_diseasegroup)
cpm_diseasegroup_sig <- cpm_diseasegroup[sig_species_diseasegroup, , drop = FALSE]
setdiff(sig_species_diseasegroup, rownames(cpm_diseasegroup))
scaled_cpm_diseasegroup <- t(scale(t(cpm_diseasegroup)))
scaled_cpm_diseasegroup[!is.finite(scaled_cpm_diseasegroup)] <- 0
annotation_col_diseasegroup <- data.frame(
  Disease = metadata_diseasegroup[colnames(cpm_diseasegroup), "Disease"]
)
rownames(annotation_col_diseasegroup) <- colnames(cpm_diseasegroup)
sample_ordered_disease <- rownames(annotation_col_diseasegroup)[
  order(annotation_col_diseasegroup$Disease)
]
scaled_cpm_diseasegroup_ordered <- scaled_cpm_diseasegroup[, sample_ordered_disease]
annotation_col_diseasegroup_ordered <- annotation_col_diseasegroup[sample_ordered_disease, , drop = FALSE]
annotation_colors_diseasegroup <- list(
  Disease = c("Control" = "green", 
              "Moderate COPD" = "pink", 
              "Severe COPD" = "blue")
)
pheatmap(scaled_cpm_diseasegroup_ordered,  
         annotation_col = annotation_col_diseasegroup_ordered,
         annotation_colors = annotation_colors_diseasegroup,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 3,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species (Disease Group)")
#top50
top50_diseasegroup <- rownames(res_diseasegroup[order(res_diseasegroup$FDR), ])[1:50]
cpm_diseasegroup_top50 <- cpm_diseasegroup[top50_diseasegroup, , drop = FALSE]
scaled_cpm_diseasegroup_top50 <- t(scale(t(cpm_diseasegroup_top50)))
scaled_cpm_diseasegroup_top50[!is.finite(scaled_cpm_diseasegroup_top50)] <- 0  # handle NaNs
annotation_col_diseasegroup_top50 <- data.frame(
  Disease = metadata_diseasegroup[colnames(cpm_diseasegroup_top50), "Disease"]
)
rownames(annotation_col_diseasegroup_top50) <- colnames(cpm_diseasegroup_top50)
sample_ordered_disease_top50 <- rownames(annotation_col_diseasegroup_top50)[
  order(annotation_col_diseasegroup_top50$Disease)
]
scaled_cpm_diseasegroup_top50_ordered <- scaled_cpm_diseasegroup_top50[, sample_ordered_disease_top50]
annotation_col_diseasegroup_top50_ordered <- annotation_col_diseasegroup_top50[sample_ordered_disease_top50, , drop = FALSE]
abs_max <- max(abs(scaled_cpm_diseasegroup_top50_ordered), na.rm = TRUE)
breaks <- seq(-abs_max, abs_max, length.out = 101)
pheatmap(scaled_cpm_diseasegroup_top50_ordered,  
         annotation_col = annotation_col_diseasegroup_top50_ordered,
         annotation_colors = annotation_colors_diseasegroup,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 8,
         color = custom_colors,
         breaks = breaks,
         main = "Top 50 Significant Species (Disease Group)")

#differential analysis for gender
group_gender <- factor(metadata$Gender, levels = c("male", "female")) #create group factor
dge_gender <- DGEList(counts = raw_count_table, group = group_gender)
dge_gender <- calcNormFactors(dge_gender)
design_gender <- model.matrix(~ 0 + group_gender)
dge_gender <- estimateDisp(dge_gender, design_gender)
fit_gender <- glmFit(dge_gender, design_gender)
contrast_gender <- makeContrasts(Female_vs_Male = group_genderfemale - group_gendermale, levels = design_gender)
lrt_gender <- glmLRT(fit_gender, contrast = contrast_gender)
res_gender <- topTags(lrt_gender, n= Inf)$table
res_gender$FDR <- p.adjust(res_gender$PValue, method = "fdr")
res_gender$Species <- rownames(res_gender)
res_gender <- res_gender %>%
  filter(FDR < 0.05)
#bar chart
res_gender$Higher_in <- ifelse(res_gender$logFC > 0, "Female", "Male")
ggplot(res_gender, aes(x = reorder(Species, logFC), y = logFC, fill = Higher_in)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Female" = "red", "Male" = "blue"),
    name = "Higher in"
  ) +
  labs(title = "Top Differentially Abundant Species (Female vs Male)",
       x = "Species",
       y = "log2 Fold Change (Female - Male)") +
  theme_minimal()
#volcano plot
res_gender <- res_gender %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_gender, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_gender %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: Female vs Male",
    x = "log2 Fold Change (Female - Male)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#extract cpm and make the heatmap
cpm_gender <- cpm(dge_gender, log = FALSE)
sig_species_gender <- res_gender$Species
cpm_gender_sig <- cpm_gender[sig_species_gender, , drop = FALSE]
scaled_cpm_gender <- t(scale(t(cpm_gender_sig)))
annotation_col_gender <- data.frame(
  Gender = metadata[colnames(cpm_gender_sig), "Gender"]
)
rownames(annotation_col_gender) <- colnames(cpm_gender_sig)
annotation_colors_gender <- list(
  Gender = c("male" = "blue", "female" = "red")
)
gender_order <- annotation_col_gender$Gender
sample_ordered <- rownames(annotation_col_gender)[order(annotation_col_gender$Gender)]
scaled_cpm_gender_ordered <- scaled_cpm_gender[, sample_ordered]
annotation_col_gender_ordered <- annotation_col_gender[sample_ordered, , drop = FALSE]

# Calculate symmetric color scale around 0
max_val <- max(scaled_cpm_gender, na.rm = TRUE)
min_val <- min(scaled_cpm_gender, na.rm = TRUE)
abs_max <- max(abs(c(max_val, min_val)))

breaks <- seq(-abs_max, abs_max, length.out = 101)

# Define matching color palette
custom_colors <- colorRampPalette(c("darkblue", "lightblue", "white", "red", "darkred"))(100)

# Plot heatmap
pheatmap(scaled_cpm_gender_ordered,  
         annotation_col = annotation_col_gender_ordered,
         annotation_colors = annotation_colors_gender,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE, 
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 5,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species (Gender)")


#differential analysis for age group
metadata_age <- metadata %>%
  filter(!is.na(Age_Group))
raw_count_table_age <- raw_count_table[, rownames(metadata_age)]
group_age <- factor(metadata_age$Age_Group, levels = c("Young", "Old"))
dge_age <- DGEList(counts = raw_count_table_age, group = group_age)
dge_age <- calcNormFactors(dge_age)
design_age <- model.matrix(~ 0 + group_age)
dge_age <- estimateDisp(dge_age, design_age)
fit_age <- glmFit(dge_age, design_age)
contrast_age <- makeContrasts(Old_vs_Young = group_ageOld - group_ageYoung, levels = design_age)
lrt_age <- glmLRT(fit_age, contrast = contrast_age)
res_age <- topTags(lrt_age, n = Inf)$table
res_age$FDR <- p.adjust(res_age$PValue, method = "fdr")
res_age$Species <- rownames(res_age)
res_age_sig <- res_age %>%
  filter(FDR < 0.05)
top50_age_sig <- res_age_sig %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
top50_age_sig$Higher_in <- ifelse(top50_age_sig$logFC > 0, "Old", "Young")
ggplot(top50_age_sig, aes(x = reorder(Species, logFC), y = logFC, fill = Higher_in)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Old" = "red", "Young" = "blue"),
    name = "Higher in"
  ) +
  labs(title = "Top Differentially Abundant Species (Old vs Young)",
       x = "Species",
       y = "log2 Fold Change (Old - Young)") +
  theme_minimal()
#volcano plot
res_age <- res_age %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_age, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_age %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: Old vs Young",
    x = "log2 Fold Change (Old - Young)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#heatmap
cpm_age <- cpm(dge_age, log = FALSE)
sig_species_age <- res_age$Species[res_age$FDR < 0.05]
cpm_age_sig <- cpm_age[sig_species_age, , drop = FALSE]
scaled_cpm_age <- t(scale(t(cpm_age_sig)))
annotation_col_age <- data.frame(
  Age_Group = metadata_age[colnames(cpm_age_sig), "Age_Group"]
)
rownames(annotation_col_age) <- colnames(cpm_age_sig)
annotation_colors_age <- list(
  Age_Group = c("Young" = "blue", "Old" = "red")
)
sample_ordered_age <- rownames(annotation_col_age)[order(annotation_col_age$Age_Group)]
scaled_cpm_age_ordered <- scaled_cpm_age[, sample_ordered_age]
annotation_col_age_ordered <- annotation_col_age[sample_ordered_age, , drop = FALSE]
max_val <- max(scaled_cpm_age, na.rm = TRUE)
min_val <- min(scaled_cpm_age, na.rm = TRUE)
abs_max <- max(abs(c(max_val, min_val)))
breaks <- seq(-abs_max, abs_max, length.out = 101)
custom_colors <- colorRampPalette(c("darkblue", "lightblue", "white", "red", "darkred"))(100)
pheatmap(scaled_cpm_age_ordered,  
         annotation_col = annotation_col_age_ordered,
         annotation_colors = annotation_colors_age,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 5,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species (Age Group)")


#differential analysis for CMH
metadata_CMH <- metadata %>%
  filter(!is.na(CMH))
raw_count_table_CMH <- raw_count_table[, rownames(metadata_CMH)]
group_CMH <- factor(metadata_CMH$CMH, levels = c("FALSE", "TRUE"))
dge_CMH <- DGEList(counts = raw_count_table_CMH, group = group_CMH)
dge_CMH <- calcNormFactors(dge_CMH)
design_CMH <- model.matrix(~ 0 + group_CMH)
dge_CMH <- estimateDisp(dge_CMH, design_CMH)
fit_CMH <- glmFit(dge_CMH, design_CMH)
contrast_CMH <- makeContrasts(
  CMH_TRUE_vs_FALSE = group_CMHTRUE - group_CMHFALSE,
  levels = design_CMH
)
lrt_CMH <- glmLRT(fit_CMH, contrast = contrast_CMH)
res_CMH <- topTags(lrt_CMH, n = Inf)$table
res_CMH$FDR <- p.adjust(res_CMH$PValue, method = "fdr")
res_CMH$Species <- rownames(res_CMH)
res_CMH_sig <- res_CMH %>% 
  filter(FDR < 0.05)
top50_CMH_sig <- res_CMH %>%
  arrange(FDR) %>%
  slice(1:50)
top50_CMH_sig$Higher_in <- ifelse(top50_CMH_sig$logFC > 0, "CMH-Positive", "CMH-Negative")
#bar chart
ggplot(top50_CMH_sig, aes(x = reorder(Species, logFC), y = logFC, fill = Higher_in)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("CMH-Positive" = "red", "CMH-Negative" = "blue"),
    name = "Higher in"
  ) +
  labs(title = "Top Differentially Abundant Species (CMH Positive vs Negative)",
       x = "Species",
       y = "log2 Fold Change (CMH Positive - Negative)") +
  theme_minimal()
#volcano plot
res_CMH <- res_CMH %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_CMH, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_CMH %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: CMH Positive vs Negative",
    x = "log2 Fold Change (CMH Positive - Negative)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#heatmap
cpm_CMH <- cpm(dge_CMH, log = FALSE)
sig_species_CMH <- res_CMH$Species[res_CMH$FDR < 0.05]
cpm_CMH_sig <- cpm_CMH[sig_species_CMH, , drop = FALSE]
scaled_cpm_CMH <- t(scale(t(cpm_CMH_sig)))
scaled_cpm_CMH[!is.finite(scaled_cpm_CMH)] <- 0
annotation_col_CMH <- data.frame(
  CMH = metadata_CMH[colnames(cpm_CMH_sig), "CMH"]
)
rownames(annotation_col_CMH) <- colnames(cpm_CMH_sig)
annotation_colors_CMH <- list(
  CMH = c("FALSE" = "blue", "TRUE" = "red")
)
sample_ordered_CMH <- rownames(annotation_col_CMH)[order(annotation_col_CMH$CMH)]
scaled_cpm_CMH_ordered <- scaled_cpm_CMH[, sample_ordered_CMH]
annotation_col_CMH_ordered <- annotation_col_CMH[sample_ordered_CMH, , drop = FALSE]
max_val <- max(scaled_cpm_CMH_ordered, na.rm = TRUE)
min_val <- min(scaled_cpm_CMH_ordered, na.rm = TRUE)
abs_max <- max(abs(c(max_val, min_val)))
breaks <- seq(-abs_max, abs_max, length.out = 101)
custom_colors <- colorRampPalette(c("darkblue", "lightblue", "white", "red", "darkred"))(100)
scaled_cpm_CMH_ordered <- as.matrix(scaled_cpm_CMH_ordered)
annotation_col_CMH_ordered$CMH <- as.factor(annotation_col_CMH_ordered$CMH)
pheatmap(scaled_cpm_CMH_ordered,  
         annotation_col = annotation_col_CMH_ordered,
         annotation_colors = annotation_colors_CMH,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 5,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species (CMH)")


#calculate diversity indices
install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages("ggpubr")
library(phyloseq)
library(vegan)
library(ggpubr)

raw_count_table <- read.delim("20250414_Sherlock1_bronchial_brushes_bact_readcounts.txt", header=TRUE, stringsAsFactors = FALSE)
rownames(raw_count_table) <- raw_count_table$Name
raw_count_table <- raw_count_table[,-1]
str(raw_count_table)
otu1 <- otu_table(raw_count_table, taxa_are_rows = TRUE)
meta <- sample_data(metadata_diseasegroup)

ps1 <- phyloseq(otu1, meta)
ps1
#alpha diversity index: Shannon
shannon_diversity <- estimate_richness(ps1, measures = "Shannon")
shannon_df <- data.frame(Sample = row.names(shannon_diversity), Shannon = shannon_diversity$Shannon)
shannon_df <- merge(shannon_df, sample_data(ps1), by = "Sample")
row.names(shannon_df) <- shannon_df$Sample
merged_shannon <- cbind(metadata[row.names(shannon_df), ], shannon_df[, "Shannon"])
print(colnames(merged_shannon))
colnames(merged_shannon)[colnames(merged_shannon) == 'shannon_df[, "Shannon"]'] <- "Shannon"
# gender
ggboxplot(
  merged_shannon,
  x = "Gender",
  y = "Shannon",
  fill = "Gender",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("male", "female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
# disease
ggboxplot(
  merged_shannon,
  x = "COPD",
  y = "Shannon",
  fill = "COPD",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#age group
ggboxplot(
  merged_shannon,
  x = "Age_Group",
  y = "Shannon",
  fill = "Age_Group",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Young", "Old")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#disease classification
ggboxplot(
  merged_shannon,
  x = "GOLD",
  y = "Shannon",
  fill = "GOLD",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_shannon$Shannon, na.rm = TRUE) + 8  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_shannon$GOLD), 2, simplify = FALSE),
    label = "p.signif",       
    hide.ns = FALSE            
  ) +
  theme_minimal() +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 7) 
  )
#for disease group 
merged_shannon_clean <- merged_shannon %>%
  filter(!is.na(Disease), !is.na(Shannon))
ggboxplot(
  merged_shannon_clean,
  x = "Disease",
  y = "Shannon",
  fill = "Disease",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_shannon_clean$Shannon, na.rm = TRUE) + 8
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_shannon_clean$Disease), 2, simplify = FALSE),
    label = "p.signif",
    hide.ns = FALSE
  ) +
  theme_minimal() +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 7)
  )
#CMH
merged_shannon_clean_CMH <- merged_shannon %>%
  filter(!is.na(CMH), !is.na(Shannon))
ggboxplot(
  merged_shannon_clean_CMH,
  x = "CMH",
  y = "Shannon",
  fill = "CMH",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
 #calculate alpha diversity indices 
alpha_df <- estimate_richness(
  ps1,
  measures = c("Observed", "Shannon", "Simpson", "Chao1")
)

merged_alpha <- cbind(metadata[rownames(alpha_df), ], alpha_df)

#Observe
#for gender
ggboxplot(
  merged_alpha,
  x = "Gender",
  y = "Observed",
  fill = "Gender",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("male", "female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
# disease
ggboxplot(
  merged_alpha,
  x = "COPD",
  y = "Observed",
  fill = "COPD",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#age group
ggboxplot(
  merged_alpha,
  x = "Age_Group",
  y = "Observed",
  fill = "Age_Group",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Young", "Old")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#disease group
ggboxplot(
  merged_alpha,
  x = "Disease",
  y = "Observed",
  fill = "Disease",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha$Observed, na.rm = TRUE) + 100  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha$Disease), 2, simplify = FALSE),
    label = "p.signif",     
    hide.ns = FALSE, 
    exact = FALSE
  ) +
  theme_minimal() + theme_pub +
  theme(
    axis.text.x = element_text(size = 7) 
  )
#CMH
merged_alpha_clean_CMH <- merged_alpha %>%
  filter(!is.na(CMH))
ggboxplot(
  merged_alpha_clean_CMH,
  x = "CMH",
  y = "Observed",
  fill = "CMH",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#Chao1
#for gender
ggboxplot(
  merged_alpha,
  x = "Gender",
  y = "Chao1",
  fill = "Gender",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("male", "female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
# disease
ggboxplot(
  merged_alpha,
  x = "COPD",
  y = "Chao1",
  fill = "COPD",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#age group
ggboxplot(
  merged_alpha,
  x = "Age_Group",
  y = "Chao1",
  fill = "Age_Group",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Young", "Old")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#disease group
ggboxplot(
  merged_alpha,
  x = "Disease",
  y = "Chao1",
  fill = "Disease",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha$Chao1, na.rm = TRUE) + 150  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha$Disease), 2, simplify = FALSE),
    label = "p.signif",       
    hide.ns = FALSE, 
    exact = FALSE
  ) +
  theme_minimal() + theme_pub +
  theme(
    axis.text.x = element_text(size = 7) 
  )
#CMH
ggboxplot(
  merged_alpha_clean_CMH,
  x = "CMH",
  y = "Chao1",
  fill = "CMH",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#Simpson
#for gender
ggboxplot(
  merged_alpha,
  x = "Gender",
  y = "Simpson",
  fill = "Gender",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("male", "female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() +theme_pub
# disease
ggboxplot(
  merged_alpha,
  x = "COPD",
  y = "Simpson",
  fill = "COPD",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#age group
ggboxplot(
  merged_alpha,
  x = "Age_Group",
  y = "Simpson",
  fill = "Age_Group",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Young", "Old")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#disease group
ggboxplot(
  merged_alpha,
  x = "Disease",
  y = "Simpson",
  fill = "Disease",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha$Simpson, na.rm = TRUE) + 2  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha$Disease), 2, simplify = FALSE),
    label = "p.signif",       
    hide.ns = FALSE,
    exact = FALSE
  ) +
  theme_minimal() + theme_pub +
  theme(
    axis.text.x = element_text(size = 7) 
  )
#CMH
ggboxplot(
  merged_alpha_clean_CMH,
  x = "CMH",
  y = "Simpson",
  fill = "CMH",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#calculate Bray-Curtis beta diversity index
otu2 <- otu_table(read_count_table_1, taxa_are_rows = TRUE)
meta <- sample_data(metadata_diseasegroup)
ps2 <- phyloseq(otu2, meta)
ps2
bray_curtis <- distance(ps2, method = "bray")
print(bray_curtis)
pcoa <- ordinate(ps2, method = "PCoA", distance = bray_curtis)
print(pcoa)
# Get percent variance explained by PC1 and PC2
var_exp <- round(pcoa$values$Relative_eig[1:2] * 100, 1)
# Create data frame for plotting
plot_data <- plot_ordination(ps2, pcoa, type = "samples", justDF = TRUE)
# Add grouping variable (gender)
plot_data$Gender <- sample_data(ps2)$Gender
# Run PERMANOVA (replace "Gender" with your grouping variable)
theme_pub <- theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )
perm_test_gender <- adonis2(bray_curtis ~ Gender, data = metadata_diseasegroup)
pval_gender <- signif(perm_test_gender$`Pr(>F)`[1], 3)
ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = Gender)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) by Gender"
  ) +
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_gender), 
           hjust = 1, size = 3.5) + theme_pub

# Add grouping variable (disease)
plot_data$COPD <- sample_data(ps2)$COPD
perm_test_disease <- adonis2(bray_curtis ~ COPD, data = metadata_diseasegroup)
pval_disease <- signif(perm_test_disease$`Pr(>F)`[1], 3)
ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = COPD)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +  # 95% confidence ellipses
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) (disease)"
  ) +
  annotate("text", 
           x = 0.5,  
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_disease), 
           hjust = 1, size = 3.5) + theme_pub

# Add grouping variable (age group)
plot_data$Age_Group <- sample_data(ps2)$Age_Group
perm_test_ageGroup <- adonis2(bray_curtis ~ Age_Group, data = metadata_diseasegroup)
pval_ageGroup <- signif(perm_test_ageGroup$`Pr(>F)`[1], 3)
ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = Age_Group)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +  # 95% confidence ellipses
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) (age group)"
  ) +
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_ageGroup), 
           hjust = 1, size = 3.5) + theme_pub


# Add grouping variable (disease group)
plot_data$Disease <- sample_data(ps2)$Disease
perm_test_diseaseGroup <- adonis2(bray_curtis ~ Disease, data = metadata_diseasegroup)
pval_diseaseGroup <- signif(perm_test_diseaseGroup$`Pr(>F)`[1], 3)
ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = Disease)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +  # 95% confidence ellipses
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) (Disease group)"
  ) +
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_diseaseGroup), 
           hjust = 1, size = 3.5) + theme_pub

#pairwise permanova
# Set the grouping variable
group_levels <- unique(metadata_diseasegroup$Disease)
# Store results
pairwise_results <- list()
# Loop through all pairwise combinations
for (i in 1:(length(group_levels)-1)) {
  for (j in (i+1):length(group_levels)) {
    # Define the groups
    group1 <- group_levels[i]
    group2 <- group_levels[j]
    # Subset metadata for just those two groups
    subset_metadata <- metadata_diseasegroup[metadata_diseasegroup$Disease %in% c(group1, group2), , drop = FALSE]
    # Subset the distance matrix accordingly
    subset_dist <- as.dist(as.matrix(bray_curtis)[rownames(subset_metadata), rownames(subset_metadata)])
    # Run adonis2
    result <- adonis2(subset_dist ~ Disease, data = subset_metadata)
    pairwise_results[[paste(group1, "vs", group2)]] <- result
  }
}
pairwise_results
# Add grouping variable (CMH)
plot_data$CMH <- sample_data(ps2)$CMH
plot_data_clean <- plot_data[!is.na(plot_data$CMH), ]
# Convert sample_data to a data frame
meta1 <- as(sample_data(ps2), "data.frame")
# Keep only samples with non-NA CMH
meta_clean1 <- meta1[!is.na(meta1$CMH), , drop = FALSE]
# Get sample names
common_samples_B <- rownames(meta_clean)
# Subset distance matrix
bray_curtis_clean <- as.matrix(bray_curtis)[common_samples_B, common_samples_B]
perm_test_CMH_B <- adonis2(bray_curtis_clean ~ CMH, data = meta_clean1)
pval_CMH_B <- signif(perm_test_CMH_B$`Pr(>F)`[1], 3)
print(perm_test_CMH_B)
ggplot(plot_data_clean, aes(x = Axis.1, y = Axis.2, color = CMH)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) in bronchial samples (CMH)"
  ) + 
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_CMH_B), 
           hjust = 1, vjust = 0, size = 4) +
  theme_pub

#taxSEA
BiocManager::install("TaxSEA")
library(TaxSEA)
#gender
taxsea_data_gender <- setNames(res_gender$logFC, res_gender$Species)
taxsea_results_gender <- TaxSEA(taxon_ranks = taxsea_data_gender)
metabolite_gender <- taxsea_results_gender$Metabolite_producers
metabolite_gender <- metabolite_gender %>%
  mutate(logP = -log10(PValue),
         direction = ifelse(median_rank_of_set_members > 0, "Enriched", "Depleted"),
         signed_logP = ifelse(direction == "Enriched", logP, -logP)) %>%
  filter(FDR < 0.05) #nothing is statistically significant
health_gender <- taxsea_results_gender$Health_associations
health_gender <- health_gender %>%
  mutate(logP = -log10(PValue),
         direction = ifelse(median_rank_of_set_members > 0, "Enriched", "Depleted"),
         signed_logP = ifelse(direction == "Enriched", logP, -logP)) %>%
  filter(FDR < 0.05) #nothing is statistically significant

#age group
taxsea_data_agegroup <- setNames(res_age_sig$logFC, res_age_sig$Species)
taxsea_results_agegroup <- TaxSEA(taxon_ranks = taxsea_data_agegroup)
metabolite_agegroup <- taxsea_results_agegroup$Metabolite_producers
health_agegroup <- taxsea_results_agegroup$Health_associations
bugsig_agegroup <- taxsea_results_agegroup$BugSigDB
metabolite_agegroup <- metabolite_agegroup %>%
  mutate(logP = -log10(PValue),
         direction = ifelse(median_rank_of_set_members > 0, "Enriched", "Depleted"),
         signed_logP = ifelse(direction == "Enriched", logP, -logP)) %>%
  filter(FDR < 0.05) #nothing is statistically significant

health_agegroup <- health_agegroup %>%
  mutate(logP = -log10(PValue),
         direction = ifelse(median_rank_of_set_members > 0, "Enriched", "Depleted"),
         signed_logP = ifelse(direction == "Enriched", logP, -logP)) %>%
  filter(FDR < 0.05) #nothing is statistically significant

bugsig_agegroup <- bugsig_agegroup %>%
  mutate(logP = -log10(PValue),
         direction = ifelse(median_rank_of_set_members > 0, "Enriched", "Depleted"),
         signed_logP = ifelse(direction == "Enriched", logP, -logP)) %>%
  filter(FDR < 0.05) #nothing is statistically significant


#for CMH 
taxsea_data_CMH <- setNames(res_CMH_sig$logFC, res_CMH_sig$Species)
taxsea_results_CMH <- TaxSEA(taxon_ranks = taxsea_data_CMH)
metabolite_CMH <- taxsea_results_CMH$Metabolite_producers #nothing is statistically significant
health_CMH <- taxsea_results_CMH$Health_associations #nothing is statistically significant
bugsig_CMH <- taxsea_results_CMH$BugSigDB #nothing is statistically significant

#for disease group
#control vs moderate COPD
taxsea_data_MvsC <- setNames(res_mod_ctrl_sig$logFC, res_mod_ctrl_sig$Species)
taxsea_results_MvsC <- TaxSEA(taxon_ranks = taxsea_data_MvsC)
metabolite_MvsC <- taxsea_results_MvsC$Metabolite_producers #FDR > 0.05

#control vs severe COPD
taxsea_data_SvsC <- setNames(res_sev_ctrl_sig$logFC, res_sev_ctrl_sig$Species)
taxsea_results_SvsC <- TaxSEA(taxon_ranks = taxsea_data_SvsC) 
metabolite_SvsC <- taxsea_results_SvsC$Metabolite_producers #FDR > 0.05
health_SvsC <- taxsea_results_SvsC$Health_associations 
bugsig_SvsC <- taxsea_results_SvsC$BugSigDB 

health_SvsC <- health_SvsC %>%
  mutate(logP = -log10(PValue),
         direction = ifelse(median_rank_of_set_members > 0, "Enriched", "Depleted"),
         signed_logP = ifelse(direction == "Enriched", logP, -logP)) %>%
  filter(FDR < 0.05)
print(health_SvsC$TaxonSet)
ggplot(health_SvsC, aes(x = signed_logP, y = reorder(taxonSetName, signed_logP), fill = direction)) +
  geom_col(width = 0.6) + 
  scale_fill_manual(values = c("Enriched" = "red", "Depleted" = "blue")) +
  labs(
    x = expression(-log[10](PValue)),
    y = "Taxon Set Name",
    title = "Taxon Set Enrichment Analysis (Health Association) (Severe COPD vs Control)"
  ) +
  theme_minimal() + 
  theme_pub

print(health_SvsC$TaxonSet)
library(stringr)
bugsig_SvsC <- bugsig_SvsC %>%
  mutate(logP = -log10(PValue),
         direction = ifelse(median_rank_of_set_members > 0, "Enriched", "Depleted"),
         signed_logP = ifelse(direction == "Enriched", logP, -logP)) %>%
  filter(FDR < 0.05)
ggplot(bugsig_SvsC, aes(x = signed_logP, y = reorder(taxonSetName, signed_logP), fill = direction)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Enriched" = "red", "Depleted" = "blue")) +
  labs(
    x = expression(-log[10](PValue)),
    y = "Taxon Set Name",
    title = "Taxon Set Enrichment Analysis (BugSigDB) (Severe COPD vs control)"
  ) + scale_y_discrete(labels = function(x) gsub("_", "_\n", x)) +
  theme_minimal() + theme_pub
print(bugsig_SvsC$TaxonSet)
#severe vs moderate COPD
taxsea_data_SvsM <- setNames(res_sev_mod_sig$logFC, res_sev_mod_sig$Species)
taxsea_results_SvsM <- TaxSEA(taxon_ranks = taxsea_data_SvsM)
metabolite_SvsM <- taxsea_results_SvsM$Metabolite_producers #FDR > 0.05
health_SvsM <- taxsea_results_SvsM$Health_associations 
bugsig_SvsM <- taxsea_results_SvsM$BugSigDB 

health_SvsM <- health_SvsM %>%
  mutate(logP = -log10(PValue),
         direction = ifelse(median_rank_of_set_members > 0, "Enriched", "Depleted"),
         signed_logP = ifelse(direction == "Enriched", logP, -logP)) %>%
  filter(FDR < 0.05)
ggplot(health_SvsM, aes(x = signed_logP, y = reorder(taxonSetName, signed_logP), fill = direction)) +
  geom_col(width=0.5) +
  scale_fill_manual(values = c("Enriched" = "red", "Depleted" = "blue")) +
  labs(
    x = expression(-log[10](PValue)),
    y = "Taxon Set Name",
    title = "Taxon Set Enrichment Analysis (Health Association) (Severe vs Moderate COPD)"
  ) +
  theme_minimal() + theme_pub

bugsig_SvsM <- bugsig_SvsM %>%
  mutate(logP = -log10(PValue),
         direction = ifelse(median_rank_of_set_members > 0, "Enriched", "Depleted"),
         signed_logP = ifelse(direction == "Enriched", logP, -logP)) %>%
  filter(FDR < 0.05)
ggplot(bugsig_SvsM, aes(x = signed_logP, y = reorder(taxonSetName, signed_logP), fill = direction)) +
  geom_col(width=0.5) +
  scale_fill_manual(values = c("Enriched" = "red", "Depleted" = "blue")) +
  labs(
    x = expression(-log[10](PValue)),
    y = "Taxon Set Name",
    title = "Taxon Set Enrichment Analysis (BugSigDB) (Severe vs moderate COPD)"
  ) + scale_y_discrete(labels = function(x) gsub("_", "_\n", x)) +
  theme_minimal() + theme_pub


XB_readcounts <- read.delim("20250423_Sherlock1_XB_cpms.txt", header = TRUE, stringsAsFactors = FALSE)
XB_rawcounts <- read.delim("20250423_Sherlock1_XB_readcounts.txt", header = TRUE, stringsAsFactors = FALSE)
XB_metadata <- read.delim("20250423_Sherlock1_XB_samples.txt", header = TRUE, stringsAsFactors = FALSE)
unique(XB_metadata$Tissue)

row.names(XB_readcounts) <- XB_readcounts$Name
#Nasal brush sample
NB_metadata <- XB_metadata %>%
  filter(Tissue == "Nasal Brush")
NB_sample_ids <- row.names(NB_metadata)
NB_readcounts <- XB_readcounts[, NB_sample_ids]
NB_rawcounts <- XB_rawcounts[,NB_sample_ids]
NB_metadata$Disease <- case_when(
  NB_metadata$GOLD %in% c("GOLD 1", "GOLD 2") ~ "Moderate COPD",
  NB_metadata$GOLD %in% c("GOLD 3", "GOLD 4") ~ "Severe COPD",
  NB_metadata$GOLD == "Not applicable" ~ "Control",
  TRUE ~ NA_character_  # if some GOLD values don't match any of the above
)
#differential analysis
#differential analysis for disease group
library(edgeR)
group_diseasegroup_NB <- factor(metadata_diseasegroup_NB$Disease, levels = c("Control", "Moderate COPD", "Severe COPD")) #Set up group variable
dge_diseasegroup_NB <- DGEList(counts = NB_rawcounts_diseasegroup, group = group_diseasegroup_NB) #Create DGEList
dge_diseasegroup_NB <- calcNormFactors(dge_diseasegroup_NB) #TMM normalization
design_diseasegroup_NB <- model.matrix(~ 0 + group_diseasegroup_NB) # Design matrix
colnames(design_diseasegroup_NB) <- make.names(colnames(design_diseasegroup_NB))
colnames(design_diseasegroup_NB)
dge_diseasegroup_NB <- estimateDisp(dge_diseasegroup_NB, design_diseasegroup_NB) # Estimate dispersions
fit_NB <- glmFit(dge_diseasegroup_NB, design_diseasegroup_NB) # Fit GLM model
#Create contrasts for pairwise comparisons
contrast_matrix_NB <- makeContrasts(
  Mod_vs_Ctrl_NB = group_diseasegroup_NBModerate.COPD - group_diseasegroup_NBControl,
  Sev_vs_Ctrl_NB = group_diseasegroup_NBSevere.COPD - group_diseasegroup_NBControl,
  Sev_vs_Mod_NB  = group_diseasegroup_NBSevere.COPD - group_diseasegroup_NBModerate.COPD,
  levels = design_diseasegroup_NB
) # Define contrast matrix using valid names
#Run glmLRT for each pair
lrt_mod_ctrl_NB <- glmLRT(fit_NB, contrast = contrast_matrix_NB[,"Mod_vs_Ctrl_NB"])
lrt_sev_ctrl_NB <- glmLRT(fit_NB, contrast = contrast_matrix_NB[,"Sev_vs_Ctrl_NB"])
lrt_sev_mod_NB  <- glmLRT(fit_NB, contrast = contrast_matrix_NB[,"Sev_vs_Mod_NB"])
#Extract result tables
res_mod_ctrl_NB <- topTags(lrt_mod_ctrl_NB, n = Inf)$table
res_mod_ctrl_NB$FDR <- p.adjust(res_mod_ctrl_NB$PValue, method = "fdr")

res_sev_ctrl_NB <- topTags(lrt_sev_ctrl_NB, n = Inf)$table
res_sev_ctrl_NB$FDR <- p.adjust(res_sev_ctrl_NB$PValue, method = "fdr")

res_sev_mod_NB <- topTags(lrt_sev_mod_NB, n = Inf)$table
res_sev_mod_NB$FDR <- p.adjust(res_sev_mod_NB$PValue, method = "fdr")

sig_taxa_mod_ctrl_NB <- rownames(res_mod_ctrl_NB)[res_mod_ctrl_NB$FDR < 0.05]
sig_taxa_sev_ctrl_NB <- rownames(res_sev_ctrl_NB)[res_sev_ctrl_NB$FDR < 0.05]
sig_taxa_sev_mod_NB  <- rownames(res_sev_mod_NB)[res_sev_mod_NB$FDR < 0.05]

res_mod_ctrl_NB$Species <- rownames(res_mod_ctrl_NB)
res_sev_ctrl_NB$Species <- rownames(res_sev_ctrl_NB)
res_sev_mod_NB$Species  <- rownames(res_sev_mod_NB)

res_mod_ctrl_NB_sig <- res_mod_ctrl_NB %>% filter(FDR<0.05)
res_sev_ctrl_NB_sig <- res_sev_ctrl_NB %>% filter(FDR<0.05)
res_sev_mod_NB_sig <- res_sev_mod_NB %>% filter(FDR<0.05)
#top 50 moderate vs control
top50_mod_ctrl_NB <- res_mod_ctrl_NB %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
ggplot(top50_mod_ctrl_NB, aes(x = reorder(Species, logFC), y = logFC, fill = logFC > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "pink", "FALSE" = "green"),
    labels = c("FALSE" = "Control", "TRUE" = "Moderate COPD"),
    name = "Higher in"
  ) +
  labs(title = "Top 50 Differentially Abundant Species in Nasal samples (Moderate vs Control)",
       x = "Species",
       y = "log2 Fold Change (Moderate - Control)") +
  theme_minimal()
#volcano plot
res_mod_ctrl_NB <- res_mod_ctrl_NB %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_mod_ctrl_NB, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +
  geom_text_repel(
    data = res_mod_ctrl_NB %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot in nasal samples: Moderate vs Control",
    x = "log2 Fold Change (Moderate - Control)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#top 50 severe vs control
top50_sev_ctrl_NB <- res_sev_ctrl_NB %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
ggplot(top50_sev_ctrl_NB, aes(x = reorder(Species, logFC), y = logFC, fill = logFC > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "blue", "FALSE" = "green"),
    labels = c("FALSE" = "Control", "TRUE" = "Severe COPD"),
    name = "Higher in"
  ) +
  labs(title = "Top 50 Differentially Abundant Species in nasal samples (Severe vs Control)",
       x = "Species",
       y = "log2 Fold Change (Severe - Control)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6), 
  )
#volcano plot
res_sev_ctrl_NB <- res_sev_ctrl_NB %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_sev_ctrl_NB, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +
  geom_text_repel(
    data = res_sev_ctrl_NB %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot in nasal samples: Severe vs Control",
    x = "log2 Fold Change (Severe - Control)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#top 50 severe vs moderate
top50_sev_mod_NB <- res_sev_mod_NB %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
ggplot(top50_sev_mod_NB, aes(x = reorder(Species, logFC), y = logFC, fill = logFC > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "blue", "FALSE" = "pink"),
    labels = c("FALSE" = "Moderate COPD", "TRUE" = "Severe COPD"),
    name = "Higher in"
  ) +
  labs(title = "Top 50 Differentially Abundant Species in nasal samples (Severe vs Modrate)",
       x = "Species",
       y = "log2 Fold Change (Severe - Moderate)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6), 
  )
#volcano plot
res_sev_mod_NB <- res_sev_mod_NB %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_sev_mod_NB, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +
  geom_text_repel(
    data = res_sev_mod_NB %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot in nasal samples: Severe vs Moderate ",
    x = "log2 Fold Change (Severe - Moderate)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#heatmap
lrt_diseasegroup_NB <- glmLRT(fit_NB)
res_diseasegroup_NB <- topTags(lrt_diseasegroup_NB, n = Inf)$table
res_diseasegroup_NB$FDR <- p.adjust(res_diseasegroup_NB$PValue, method = "fdr")
res_diseasegroup_NB$Species <- rownames(res_diseasegroup_NB)
res_diseasegroup_NB_sig <- res_diseasegroup_NB %>% filter(FDR < 0.05)
sig_species_diseasegroup_NB <- intersect(res_diseasegroup_NB_sig$Species,
                                         rownames(dge_diseasegroup_NB))
cpm_diseasegroup_NB <- cpm(dge_diseasegroup_NB, log = FALSE)
cpm_diseasegroup_sig_NB <- cpm_diseasegroup_NB[sig_species_diseasegroup_NB, , drop = FALSE]
scaled_cpm_diseasegroup_NB <- t(scale(t(cpm_diseasegroup_sig_NB)))
scaled_cpm_diseasegroup_NB[!is.finite(scaled_cpm_diseasegroup_NB)] <- 0
annotation_col_diseasegroup_NB <- data.frame(
  Disease = metadata_diseasegroup_NB[colnames(cpm_diseasegroup_NB), "Disease"]
)
rownames(annotation_col_diseasegroup_NB) <- colnames(cpm_diseasegroup_NB)
sample_ordered_disease_NB <- rownames(annotation_col_diseasegroup_NB)[
  order(annotation_col_diseasegroup_NB$Disease)
]
scaled_cpm_diseasegroup_NB_ordered <- scaled_cpm_diseasegroup_NB[, sample_ordered_disease_NB]
annotation_col_diseasegroup_NB_ordered <- annotation_col_diseasegroup_NB[sample_ordered_disease_NB, , drop = FALSE]
annotation_colors_diseasegroup <- list(
  Disease = c("Control" = "green", 
              "Moderate COPD" = "pink", 
              "Severe COPD" = "blue")
)
abs_max <- max(abs(scaled_cpm_diseasegroup_NB_ordered), na.rm = TRUE)
breaks <- seq(-abs_max, abs_max, length.out = 101)
pheatmap(scaled_cpm_diseasegroup_NB_ordered,  
         annotation_col = annotation_col_diseasegroup_NB_ordered,
         annotation_colors = annotation_colors_diseasegroup,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 2,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species in nasal samples (Disease Group)")

#Top 50
top50_diseasegroup_NB <- rownames(res_diseasegroup_NB[order(res_diseasegroup_NB$FDR), ])[1:50]
top50_diseasegroup_NB <- intersect(top50_diseasegroup_NB, rownames(cpm_diseasegroup_NB))
cpm_diseasegroup_top50_NB <- cpm_diseasegroup_NB[top50_diseasegroup_NB, , drop = FALSE]
scaled_cpm_diseasegroup_top50_NB <- t(scale(t(cpm_diseasegroup_top50_NB)))
scaled_cpm_diseasegroup_top50_NB[!is.finite(scaled_cpm_diseasegroup_top50_NB)] <- 0
annotation_col_diseasegroup_top50_NB <- data.frame(
  Disease = metadata_diseasegroup_NB[colnames(cpm_diseasegroup_top50_NB), "Disease"]
)
rownames(annotation_col_diseasegroup_top50_NB) <- colnames(cpm_diseasegroup_top50_NB)
sample_ordered_disease_top50_NB <- rownames(annotation_col_diseasegroup_top50_NB)[
  order(annotation_col_diseasegroup_top50_NB$Disease)
]
scaled_cpm_diseasegroup_top50_NB_ordered <- scaled_cpm_diseasegroup_top50_NB[, sample_ordered_disease_top50_NB]
annotation_col_diseasegroup_top50_NB_ordered <- annotation_col_diseasegroup_top50_NB[sample_ordered_disease_top50_NB, , drop = FALSE]
pheatmap(scaled_cpm_diseasegroup_top50_NB_ordered,  
         annotation_col = annotation_col_diseasegroup_top50_NB_ordered,
         annotation_colors = annotation_colors_diseasegroup,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 8,
         color = custom_colors,
         breaks = breaks,
         main = "Top 50 Significant Species in nasal samples(Disease Group)")


#differential analysis for gender
group_gender_NB <- factor(NB_metadata$Gender, levels = c("male", "female")) #create group factor
dge_gender_NB <- DGEList(counts = NB_rawcounts, group = group_gender_NB)
dge_gender_NB <- calcNormFactors(dge_gender_NB)
design_gender_NB <- model.matrix(~ 0 + group_gender_NB)
dge_gender_NB <- estimateDisp(dge_gender_NB, design_gender_NB)
fit_gender_NB <- glmFit(dge_gender_NB, design_gender_NB)
contrast_gender_NB <- makeContrasts(Female_vs_Male_NB = group_gender_NBfemale - group_gender_NBmale, levels = design_gender_NB)
lrt_gender_NB <- glmLRT(fit_gender_NB, contrast = contrast_gender_NB)
res_gender_NB <- topTags(lrt_gender_NB, n= Inf)$table
res_gender_NB$FDR <- p.adjust(res_gender_NB$PValue, method = "fdr")
res_gender_NB$Species <- rownames(res_gender_NB)
res_gender_NB_sig <- res_gender_NB %>%
  filter(FDR < 0.05)
#bar chart
res_gender_NB_sig$Higher_in <- ifelse(res_gender_NB_sig$logFC > 0, "Female", "Male")
ggplot(res_gender_NB_sig, aes(x = reorder(Species, logFC), y = logFC, fill = Higher_in)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Female" = "red", "Male" = "blue"),
    name = "Higher in"
  ) +
  labs(title = "Top Differentially Abundant Species in nasal samples (Female vs Male)",
       x = "Species",
       y = "log2 Fold Change (Female - Male)") +
  theme_minimal()
#volcano plot
res_gender_NB <- res_gender_NB %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_gender_NB, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +
  geom_text_repel(
    data = res_gender_NB %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot in nasal samples: Female vs Male",
    x = "log2 Fold Change (Female - Male)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#extract cpm and make the heatmap
cpm_gender_NB <- cpm(dge_gender_NB, log = FALSE)
sig_species_gender_NB <- res_gender_NB_sig$Species
cpm_gender_sig_NB <- cpm_gender_NB[sig_species_gender_NB, , drop = FALSE]
scaled_cpm_gender_NB <- t(scale(t(cpm_gender_sig_NB)))
annotation_col_gender_NB <- data.frame(
  Gender = NB_metadata[colnames(cpm_gender_sig_NB), "Gender"]
)
rownames(annotation_col_gender_NB) <- colnames(cpm_gender_sig_NB)
annotation_colors_gender_NB <- list(
  Gender = c("male" = "blue", "female" = "red")
)
gender_order_NB <- annotation_col_gender_NB$Gender
sample_ordered_NB <- rownames(annotation_col_gender_NB)[order(annotation_col_gender_NB$Gender)]
scaled_cpm_gender_ordered_NB <- scaled_cpm_gender_NB[, sample_ordered_NB]
annotation_col_gender_ordered_NB <- annotation_col_gender_NB[sample_ordered_NB, , drop = FALSE]

# Calculate symmetric color scale around 0
max_val <- max(scaled_cpm_gender_NB, na.rm = TRUE)
min_val <- min(scaled_cpm_gender_NB, na.rm = TRUE)
abs_max <- max(abs(c(max_val, min_val)))

breaks <- seq(-abs_max, abs_max, length.out = 101)

# Plot heatmap
pheatmap(scaled_cpm_gender_ordered_NB,  
         annotation_col = annotation_col_gender_ordered_NB,
         annotation_colors = annotation_colors_gender_NB,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE, 
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 8,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species in nasal samples (Gender)")


#differential analysis for age group
metadata_age_NB <- NB_metadata %>%
  filter(!is.na(Age_Group))
NB_rawcounts_age <- NB_rawcounts[, rownames(metadata_age_NB)]
group_age_NB <- factor(metadata_age_NB$Age_Group, levels = c("Young", "Old"))
dge_age_NB <- DGEList(counts = NB_rawcounts_age, group = group_age_NB)
dge_age_NB <- calcNormFactors(dge_age_NB)
design_age_NB <- model.matrix(~ 0 + group_age_NB)
dge_age_NB <- estimateDisp(dge_age_NB, design_age_NB)
fit_age_NB <- glmFit(dge_age_NB, design_age_NB)
contrast_age_NB <- makeContrasts(Old_vs_Young_NB = group_age_NBOld - group_age_NBYoung, levels = design_age_NB)
lrt_age_NB <- glmLRT(fit_age_NB, contrast = contrast_age_NB)
res_age_NB <- topTags(lrt_age_NB, n = Inf)$table
res_age_NB$FDR <- p.adjust(res_age_NB$PValue, method = "fdr")
res_age_NB$Species <- rownames(res_age_NB)
res_age_NB_sig <- res_age_NB %>%
  filter(FDR < 0.05)
top50_age_NB_sig <- res_age_NB_sig %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
top50_age_NB_sig$Higher_in <- ifelse(top50_age_NB_sig$logFC > 0, "Old", "Young")
ggplot(top50_age_NB_sig, aes(x = reorder(Species, logFC), y = logFC, fill = Higher_in)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Old" = "red", "Young" = "blue"),
    name = "Higher in"
  ) +
  labs(title = "Top Differentially Abundant Species in nasal samples (Old vs Young)",
       x = "Species",
       y = "log2 Fold Change (Old - Young)") +
  theme_minimal()
#volcano plot
res_age_NB <- res_age_NB %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_age_NB, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +
  geom_text_repel(
    data = res_age_NB %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
  labs(
    title = "Volcano Plot: Old vs Young (nasal samples)",
    x = "log2 Fold Change (Old - Young)",
    y = "-log10(FDR)",
    color = "Significance"
  ) +
  theme_minimal()

#heatmap
cpm_age_NB <- cpm(dge_age_NB, log = FALSE)
sig_species_age_NB <- res_age_NB$Species[res_age_NB$FDR < 0.05]
cpm_age_sig_NB <- cpm_age_NB[sig_species_age_NB, , drop = FALSE]
scaled_cpm_age_NB <- t(scale(t(cpm_age_sig_NB)))
annotation_col_age_NB <- data.frame(
  Age_Group_NB = metadata_age_NB[colnames(cpm_age_sig_NB), "Age_Group"]
)
rownames(annotation_col_age_NB) <- colnames(cpm_age_sig_NB)
annotation_colors_age_NB <- list(
  Age_Group_NB = c("Young" = "blue", "Old" = "red")
)
sample_ordered_age_NB <- rownames(annotation_col_age_NB)[order(annotation_col_age_NB$Age_Group_NB)]
scaled_cpm_age_ordered_NB <- scaled_cpm_age_NB[, sample_ordered_age_NB]
annotation_col_age_ordered_NB <- annotation_col_age_NB[sample_ordered_age_NB, , drop = FALSE]
max_val <- max(scaled_cpm_age_NB, na.rm = TRUE)
min_val <- min(scaled_cpm_age_NB, na.rm = TRUE)
abs_max <- max(abs(c(max_val, min_val)))
breaks <- seq(-abs_max, abs_max, length.out = 101)
custom_colors <- colorRampPalette(c("darkblue", "lightblue", "white", "red", "darkred"))(100)
pheatmap(scaled_cpm_age_ordered_NB,  
         annotation_col = annotation_col_age_ordered_NB,
         annotation_colors = annotation_colors_age_NB,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 7,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species in nasal samples (Age Group)")


#differential analysis for CMH
metadata_CMH_NB <- NB_metadata %>%
  filter(!is.na(CMH))
NB_rawcounts_CMH <- NB_rawcounts[, rownames(metadata_CMH_NB)]
group_CMH_NB <- factor(metadata_CMH_NB$CMH, levels = c("FALSE", "TRUE"))
dge_CMH_NB <- DGEList(counts = NB_rawcounts_CMH, group = group_CMH_NB)
dge_CMH_NB <- calcNormFactors(dge_CMH_NB)
design_CMH_NB <- model.matrix(~ 0 + group_CMH_NB)
dge_CMH_NB <- estimateDisp(dge_CMH_NB, design_CMH_NB)
fit_CMH_NB <- glmFit(dge_CMH_NB, design_CMH_NB)
contrast_CMH_NB <- makeContrasts(
  CMH_TRUE_vs_FALSE_NB = group_CMH_NBTRUE - group_CMH_NBFALSE,
  levels = design_CMH_NB
)
lrt_CMH_NB <- glmLRT(fit_CMH_NB, contrast = contrast_CMH_NB)
res_CMH_NB <- topTags(lrt_CMH_NB, n = Inf)$table
res_CMH_NB$FDR <- p.adjust(res_CMH_NB$PValue, method = "fdr")
res_CMH_NB$Species <- rownames(res_CMH_NB)
res_CMH_NB_sig <- res_CMH_NB %>% 
  filter(FDR < 0.05)
top50_CMH_NB_sig <- res_CMH_NB_sig %>%
  arrange(FDR) %>%
  slice(1:50)
top50_CMH_NB_sig$Higher_in <- ifelse(top50_CMH_NB_sig$logFC > 0, "CMH-Positive", "CMH-Negative")
#bar chart
ggplot(top50_CMH_NB_sig, aes(x = reorder(Species, logFC), y = logFC, fill = Higher_in)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("CMH-Positive" = "red", "CMH-Negative" = "blue"),
    name = "Higher in"
  ) +
  labs(title = "Top Differentially Abundant Species in nasal species (CMH Positive vs Negative)",
       x = "Species",
       y = "log2 Fold Change (CMH Positive - Negative)") +
  theme_minimal()
#volcano plot
res_CMH_NB <- res_CMH_NB %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_CMH_NB, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_CMH_NB %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: CMH Positive vs Negative",
    x = "log2 Fold Change (CMH Positive - Negative)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#heatmap
cpm_CMH_NB <- cpm(dge_CMH_NB, log = FALSE)
sig_species_CMH_NB <- res_CMH_NB$Species[res_CMH_NB$FDR < 0.05]
cpm_CMH_NB_sig <- cpm_CMH_NB[sig_species_CMH_NB, , drop = FALSE]
scaled_cpm_CMH_NB <- t(scale(t(cpm_CMH_NB_sig)))
scaled_cpm_CMH_NB[!is.finite(scaled_cpm_CMH_NB)] <- 0
annotation_col_CMH_NB <- data.frame(
  CMH = metadata_CMH_NB[colnames(cpm_CMH_NB_sig), "CMH"]
)
rownames(annotation_col_CMH_NB) <- colnames(cpm_CMH_NB_sig)
annotation_colors_CMH_NB <- list(
  CMH = c("FALSE" = "blue", "TRUE" = "red")
)
sample_ordered_CMH_NB <- rownames(annotation_col_CMH_NB)[order(annotation_col_CMH_NB$CMH)]
scaled_cpm_CMH_ordered_NB <- scaled_cpm_CMH_NB[, sample_ordered_CMH_NB]
annotation_col_CMH_ordered_NB <- annotation_col_CMH_NB[sample_ordered_CMH_NB, , drop = FALSE]
max_val <- max(scaled_cpm_CMH_ordered_NB, na.rm = TRUE)
min_val <- min(scaled_cpm_CMH_ordered_NB, na.rm = TRUE)
abs_max <- max(abs(c(max_val, min_val)))
breaks <- seq(-abs_max, abs_max, length.out = 101)
custom_colors <- colorRampPalette(c("darkblue", "lightblue", "white", "red", "darkred"))(100)
scaled_cpm_CMH_ordered_NB <- as.matrix(scaled_cpm_CMH_ordered_NB)
annotation_col_CMH_ordered_NB$CMH <- as.factor(annotation_col_CMH_ordered_NB$CMH)
pheatmap(scaled_cpm_CMH_ordered_NB,  
         annotation_col = annotation_col_CMH_ordered_NB,
         annotation_colors = annotation_colors_CMH_NB,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 8,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species in nasal samples (CMH)")

#diversity indices for nasal samples
otu_NB <- otu_table(NB_rawcounts, taxa_are_rows = TRUE)
meta_NB <- sample_data(metadata_diseasegroup_NB)
ps_NB <- phyloseq(otu_NB, meta_NB)
ps_NB
#alpha diversity 
alpha_df_NB <- estimate_richness(ps_NB, measures = c("Observed", "Shannon", "Simpson", "Chao1"))
merged_alpha_NB <- cbind(
  metadata_diseasegroup_NB[rownames(alpha_df_NB), , drop = FALSE],
  alpha_df_NB
)
#Observed
#for gender
ggboxplot(
  merged_alpha_NB,
  x = "Gender",
  y = "Observed",
  fill = "Gender",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("male", "female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#for CMH
# Filter out rows where CMH or Observed is NA
filtered_data_NB <- merged_alpha_NB %>%
  filter(!is.na(CMH))
ggboxplot(
  filtered_data_NB,
  x = "CMH",
  y = "Observed",
  fill = "CMH",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#age group
ggboxplot(
  merged_alpha_NB,
  x = "Age_Group",
  y = "Observed",
  fill = "Age_Group",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Young", "Old")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#disease group
ggboxplot(
  merged_alpha_NB,
  x = "Disease",
  y = "Observed",
  fill = "Disease",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha_NB$Observed, na.rm = TRUE) + 100  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha_NB$Disease), 2, simplify = FALSE),
    label = "p.signif",       
    hide.ns = FALSE, 
    exact = FALSE
  ) +
  theme_minimal() + theme_pub +
  theme(
    axis.text.x = element_text(size = 7) 
  )
#Chao1
#for gender
ggboxplot(
  merged_alpha_NB,
  x = "Gender",
  y = "Chao1",
  fill = "Gender",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("male", "female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#age group
ggboxplot(
  merged_alpha_NB,
  x = "Age_Group",
  y = "Chao1",
  fill = "Age_Group",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Young", "Old")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#CMH
ggboxplot(
  filtered_data_NB,
  x = "CMH",
  y = "Chao1",
  fill = "CMH",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#disease group
ggboxplot(
  merged_alpha_NB,
  x = "Disease",
  y = "Chao1",
  fill = "Disease",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha_NB$Chao1, na.rm = TRUE) + 150  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha_NB$Disease), 2, simplify = FALSE),
    label = "p.signif",       
    hide.ns = FALSE,
    exact = FALSE
  ) +
  theme_minimal() + theme_pub +
  theme(
    axis.text.x = element_text(size = 7) 
  )
#Simpson
#for gender
ggboxplot(
  merged_alpha_NB,
  x = "Gender",
  y = "Simpson",
  fill = "Gender",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("male", "female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() +theme_pub
#age group
ggboxplot(
  merged_alpha_NB,
  x = "Age_Group",
  y = "Simpson",
  fill = "Age_Group",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Young", "Old")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#disease group
ggboxplot(
  merged_alpha_NB,
  x = "Disease",
  y = "Simpson",
  fill = "Disease",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha_NB$Simpson, na.rm = TRUE) + 2  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha_NB$Disease), 2, simplify = FALSE),
    label = "p.signif",       
    hide.ns = FALSE, 
    exact = FALSE
  ) +
  theme_minimal() + theme_pub +
  theme(
    axis.text.x = element_text(size = 7) 
  )
#CMH
ggboxplot(
  filtered_data_NB,
  x = "CMH",
  y = "Simpson",
  fill = "CMH",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#Shannon
#for gender
ggboxplot(
  merged_alpha_NB,
  x = "Gender",
  y = "Shannon",
  fill = "Gender",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("male", "female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() +theme_pub
#age group
ggboxplot(
  merged_alpha_NB,
  x = "Age_Group",
  y = "Shannon",
  fill = "Age_Group",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Young", "Old")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#disease group
ggboxplot(
  merged_alpha_NB,
  x = "Disease",
  y = "Shannon",
  fill = "Disease",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha_NB$Simpson, na.rm = TRUE) + 5  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha_NB$Disease), 2, simplify = FALSE),
    label = "p.signif",      
    hide.ns = FALSE, 
    exact = FALSE
  ) +
  theme_minimal() + theme_pub +
  theme(
    axis.text.x = element_text(size = 7) 
  )
#CMH
ggboxplot(
  filtered_data_NB,
  x = "CMH",
  y = "Shannon",
  fill = "CMH",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#beta diversity
otu_NB2 <- otu_table(NB_readcounts_diseasegroup, taxa_are_rows = TRUE)
ps_NB2 <- phyloseq(otu_NB2, meta_NB)
ps_NB2

bray_curtis_NB <- distance(ps_NB2, method="bray")
print(bray_curtis_NB)
pcoa_NB <- ordinate(ps_NB2, method = "PCoA", distance = bray_curtis_NB)
print(pcoa_NB)
# Get percent variance explained by PC1 and PC2
var_exp_NB <- round(pcoa_NB$values$Relative_eig[1:2] * 100, 1)
# Create data frame for plotting
plot_data_NB <- plot_ordination(ps_NB2, pcoa_NB, type = "samples", justDF = TRUE)
# Add grouping variable (gender)
plot_data_NB$Gender <- sample_data(ps_NB2)$Gender
# Run PERMANOVA (replace "Gender" with your grouping variable)

perm_test_gender_NB <- adonis2(bray_curtis_NB ~ Gender, data = metadata_diseasegroup_NB)
pval_gender_NB <- signif(perm_test_gender_NB$`Pr(>F)`[1], 3)
ggplot(plot_data_NB, aes(x = Axis.1, y = Axis.2, color = Gender)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) in nasal samples by Gender"
  ) +
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_gender_NB), 
           hjust = 1, size = 3.5) + theme_pub

# Add grouping variable (age group)
plot_data_NB$Age_Group <- sample_data(ps_NB2)$Age_Group
perm_test_ageGroup_NB <- adonis2(bray_curtis_NB ~ Age_Group, data = metadata_diseasegroup_NB)
pval_ageGroup_NB <- signif(perm_test_ageGroup_NB$`Pr(>F)`[1], 3)
ggplot(plot_data_NB, aes(x = Axis.1, y = Axis.2, color = Age_Group)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +  # 95% confidence ellipses
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) in nasal samples (age group)"
  ) +
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_ageGroup_NB), 
           hjust = 1, size = 3.5) + theme_pub


# Add grouping variable (disease group)
plot_data_NB$Disease <- sample_data(ps_NB2)$Disease
perm_test_diseaseGroup_NB <- adonis2(bray_curtis_NB ~ Disease, data = metadata_diseasegroup_NB)
pval_diseaseGroup_NB <- signif(perm_test_diseaseGroup_NB$`Pr(>F)`[1], 3)
ggplot(plot_data_NB, aes(x = Axis.1, y = Axis.2, color = Disease)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +  # 95% confidence ellipses
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) in nasal samples (Disease group)"
  ) +
  annotate("text", 
           x = 0.5,  # Change based on your plot scale
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_diseaseGroup_NB), 
           hjust = 1, size = 3.5) + theme_pub
#pairwise permanova
# Set the grouping variable
group_levels_NB <- unique(metadata_diseasegroup_NB$Disease)
# Store results
pairwise_results_NB <- list()
# Loop through all pairwise combinations
for (i in 1:(length(group_levels_NB)-1)) {
  for (j in (i+1):length(group_levels_NB)) {
    # Define the groups
    group1 <- group_levels_NB[i]
    group2 <- group_levels_NB[j]
    # Subset metadata for just those two groups
    subset_metadata_NB <- metadata_diseasegroup_NB[metadata_diseasegroup_NB$Disease %in% c(group1, group2), , drop = FALSE]
    # Subset the distance matrix accordingly
    subset_dist_NB <- as.dist(as.matrix(bray_curtis_NB)[rownames(subset_metadata_NB), rownames(subset_metadata_NB)])
    # Run adonis2
    result <- adonis2(subset_dist_NB ~ Disease, data = subset_metadata_NB)
    pairwise_results_NB[[paste(group1, "vs", group2)]] <- result
  }
}
# Print all results
pairwise_results_NB
# Add grouping variable (CMH)
plot_data_NB$CMH <- sample_data(ps_NB2)$CMH
lot_data_NB_clean <- plot_data_NB[!is.na(plot_data_NB$CMH), ]
# Convert sample_data to a data frame
meta_NB <- as(sample_data(ps_NB2), "data.frame")
# Keep only samples with non-NA CMH
meta_NB_clean <- meta_NB[!is.na(meta_NB$CMH), , drop = FALSE]
# Get sample names
common_samples <- rownames(meta_NB_clean)
# Subset distance matrix
bray_curtis_NB_clean <- as.matrix(bray_curtis_NB)[common_samples, common_samples]
perm_test_CMH <- adonis2(bray_curtis_NB_clean ~ CMH, data = meta_NB_clean)
pval_CMH <- signif(perm_test_CMH$`Pr(>F)`[1], 3)
print(perm_test_CMH)
ggplot(lot_data_NB_clean, aes(x = Axis.1, y = Axis.2, color = CMH)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) in nasal samples (CMH)"
  ) + 
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_CMH), 
           hjust = 1, vjust = 0, size = 4) +
  theme_pub

#taxSEA 
#gender
taxsea_data_gender_NB <- setNames(res_gender_NB_sig$logFC, res_gender_NB_sig$Species)
taxsea_results_gender_NB <- TaxSEA(taxon_ranks = taxsea_data_gender_NB)
metabolite_gender_NB <- taxsea_results_gender_NB$Metabolite_producers #nothing is statistically significant

#age group
taxsea_data_agegroup_NB <- setNames(res_age_NB_sig$logFC, res_age_NB_sig$Species)
taxsea_results_agegroup_NB <- TaxSEA(taxon_ranks = taxsea_data_agegroup_NB)
metabolite_agegroup_NB <- taxsea_results_agegroup_NB$Metabolite_producers #nothing is statisically significant
health_agegroup_NB <- taxsea_results_agegroup_NB$Health_associations #nothing is statisically significant
bugsig_agegroup_NB <- taxsea_results_agegroup_NB$BugSigDB #nothing is statisically significant

#for CMH 
taxsea_data_CMH_NB <- setNames(res_CMH_NB_sig$logFC, res_CMH_NB_sig$Species)
taxsea_results_CMH_NB <- TaxSEA(taxon_ranks = taxsea_data_CMH_NB)
metabolite_CMH_NB <- taxsea_results_CMH_NB$Metabolite_producers #nothing is statistically significant

#for disease group
#control vs moderate COPD
taxsea_data_MvsC_NB <- setNames(res_mod_ctrl_NB_sig$logFC, res_mod_ctrl_NB_sig$Species)
taxsea_results_MvsC_NB <- TaxSEA(taxon_ranks = taxsea_data_MvsC_NB)

#control vs severe COPD
taxsea_data_SvsC_NB <- setNames(res_sev_ctrl_NB_sig$logFC, res_sev_ctrl_NB_sig$Species)
taxsea_results_SvsC_NB <- TaxSEA(taxon_ranks = taxsea_data_SvsC_NB) 
metabolite_SvsC_NB <- taxsea_results_SvsC_NB$Metabolite_producers #FDR > 0.05
health_SvsC_NB <- taxsea_results_SvsC_NB$Health_associations #FDR > 0.05
bugsig_SvsC_NB <- taxsea_results_SvsC_NB$BugSigDB #FDR > 0.05

#severe vs moderate COPD
taxsea_data_SvsM_NB <- setNames(res_sev_mod_NB_sig$logFC, res_sev_mod_NB_sig$Species)
taxsea_results_SvsM_NB <- TaxSEA(taxon_ranks = taxsea_data_SvsM_NB)
metabolite_SvsM_NB <- taxsea_results_SvsM_NB$Metabolite_producers #FDR > 0.05
health_SvsM_NB <- taxsea_results_SvsM_NB$Health_associations #FDR > 0.05
bugsig_SvsM_NB <- taxsea_results_SvsM_NB$BugSigDB #FDR > 0.05

#nasal brush vs bronchial brush
XB_readcounts <- XB_readcounts[,-1]
row.names(XB_readcounts) <- XB_readcounts[,1]
XB_readcounts_1 <- XB_readcounts[,-1]
row.names(XB_metadata) <- XB_metadata[,1]
XB_metadata <- XB_metadata[,-1]
rank_test_xb <- apply(XB_readcounts_1, 1, function(x) {
  group <- XB_metadata$Tissue
  test_result <- wilcox.test(x~group, exact = FALSE)
  return(test_result$p.value)
})
print(rank_test_xb)
p_adj_xb <- p.adjust(rank_test_xb, method = "fdr")
rank_test_results_xb <- data.frame(Species=rownames(XB_readcounts_1),
                                       p_value = rank_test_xb, p_adj = p_adj_xb)
print(rank_test_results_xb)
significant_species_xb <- rank_test_results_xb$Species[rank_test_results_xb$p_adj < threshold] #filter species with p value < 0.05
print(significant_species_xb)
filtered_data_xb <- XB_readcounts_1[significant_species_xb,,drop=FALSE] #subset to keep only significant species
scaled_data_xb <- t(scale(t(filtered_data_xb)))
ordered_samples_xb <- order(XB_metadata$Tissue) #sort by tissue
print(ordered_samples_xb)
sample_names_ordered_xb <- rownames(XB_metadata)[ordered_samples_xb] #get ordered sample names
setdiff(sample_names_ordered_xb, colnames(scaled_data_xb))
#reorder scaled data and metadata
scaled_data_xb <- scaled_data_xb[, sample_names_ordered_xb] #reorder columns in data 
metadata_ordered_xb <- XB_metadata[sample_names_ordered_xb,,drop=FALSE] #reorder metadata rows
#generation annotation data after ordering
annotation_col_xb <- data.frame(Tissue = metadata_ordered_xb$Tissue, row.names=colnames(scaled_data_xb))
annotation_colors_xb <- list(Tissue = c("Lung Brush" = "red", "Nasal Brush" = "blue"))
custom_colors <- colorRampPalette(c("darkblue", "lightblue", "white", "red", "darkred"))(100)
#generate heatmap
pheatmap(scaled_data_xb,
         annotation_col = annotation_col_xb,
         annotation_colors = annotation_colors_xb,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         scale = "row",
         cluster_cols = FALSE,
         color = custom_colors,
         fontsize_row = 2,
         fontsize_col = 1,
         labels_col = NA, 
         main = "Heatmap of Significant Bacterial Abundance (Lung Brush vs. Nasal Brush)")

## select top 50 most significant based on FDR
top50_species_xb <- rank_test_results_xb %>%
  arrange(p_adj) %>%
  filter(!is.na(p_adj)) %>%  # remove NA if any
  slice_head(n = 50) %>%
  pull(Species)
print(top50_species_xb)
filtered_data_xb_50 <- XB_readcounts_1[top50_species_xb, , drop = FALSE]
scaled_data_xb_50 <- t(scale(t(filtered_data_xb_50))) 
ordered_samples_xb_50 <- order(XB_metadata$Tissue) 
sample_names_ordered_xb_50 <- rownames(XB_metadata)[ordered_samples_xb_50] #get ordered sample names
#reorder scaled data and metadata
scaled_data_xb_50 <- scaled_data_xb_50[, sample_names_ordered_xb_50] #reorder columns in data 
metadata_ordered_xb_50 <- XB_metadata[sample_names_ordered_xb_50,,drop=FALSE] #reorder metadata rows
#generation annotation data after ordering
annotation_col_xb_50 <- data.frame(Tissue = as.factor(metadata_ordered_xb_50$Tissue), 
                                             row.names = colnames(scaled_data_xb_50))
annotation_colors_xb <- list(Tissue = c(
  "Lung Brush" = "red",
  "Nasal Brush" = "blue"
))
custom_colors <- colorRampPalette(c("darkblue", "lightblue", "white", "red", "darkred"))(100)
#generate heatmap
pheatmap(scaled_data_xb_50,
         annotation_col = annotation_col_xb_50,
         annotation_colors = annotation_colors_xb,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         scale = "row",
         cluster_cols = FALSE,
         color = custom_colors,
         fontsize_row = 8,
         main = "Heatmap of Significant Bacterial Abundance (Lung Brush vs Nasal Brush)",
         show_colnames = FALSE) 


#differential analysis
group_xb <- factor(XB_metadata$Tissue, levels = c("Nasal Brush", "Lung Brush")) #create group factor
dge_xb <- DGEList(counts = XB_rawcounts, group = group_xb)
dge_xb <- calcNormFactors(dge_xb)
design_xb <- model.matrix(~ 0 + group_xb)
dge_xb <- estimateDisp(dge_xb, design_xb)
fit_xb <- glmFit(dge_xb, design_xb)
colnames(design_xb) <- make.names(colnames(design_xb))
colnames(design_xb)
contrast_xb <- makeContrasts(Nasal_vs_Lung = group_xbNasal.Brush - group_xbLung.Brush, levels = design_xb)
lrt_xb <- glmLRT(fit_xb, contrast = contrast_xb)
res_xb <- topTags(lrt_xb, n= Inf)$table
res_xb$FDR <- p.adjust(res_xb$PValue, method = "fdr")
res_xb$Species <- rownames(res_xb)
res_xb_sig <- res_xb %>%
  filter(FDR < 0.05)
top50_xb <- res_xb_sig %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
top50_xb$Higher_in <- ifelse(top50_xb$logFC > 0, "Nasal Brush", "Lung Brush")
ggplot(top50_xb, aes(x = reorder(Species, logFC), y = logFC, fill = Higher_in)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Nasal Brush" = "red", "Lung Brush" = "blue"),
    name = "Higher in"
  ) +
  labs(title = "Top Differentially Abundant Species (Nasal Brush vs Lung Brush)",
       x = "Species",
       y = "log2 Fold Change (Nasal Brush - Lung Brush)") +
  theme_minimal()
#volcano plot
res_xb <- res_xb %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_xb, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_xb %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: Nasal vs Lung",
    x = "log2 Fold Change (Nasal Brush - Lung Brush)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#heatmap of unpaired analysis
cpm_xb <- cpm(dge_xb, log = FALSE)
sig_species_xb <- res_xb$Species[res_xb$FDR < 0.05]
cpm_xb_sig <- cpm_xb[sig_species_xb, , drop = FALSE]
scaled_cpm_xb <- t(scale(t(cpm_xb_sig)))
scaled_cpm_xb[!is.finite(scaled_cpm_xb)] <- 0
annotation_col_xb <- data.frame(
  Tissue = XB_metadata[colnames(cpm_xb_sig), "Tissue"]
)
rownames(annotation_col_xb) <- colnames(cpm_xb_sig)
annotation_colors_CMH <- list(
  Tissue = c("Nasal Brush" = "blue", "Lung Brush" = "red")
)
sample_ordered_xb <- rownames(annotation_col_xb)[order(annotation_col_xb$Tissue)]
scaled_cpm_xb_ordered <- scaled_cpm_xb[, sample_ordered_xb]
annotation_col_xb_ordered <- annotation_col_xb[sample_ordered_xb, , drop = FALSE]
max_val <- max(scaled_cpm_xb_ordered, na.rm = TRUE)
min_val <- min(scaled_cpm_xb_ordered, na.rm = TRUE)
abs_max <- max(abs(c(max_val, min_val)))
breaks <- seq(-abs_max, abs_max, length.out = 101)
custom_colors <- colorRampPalette(c("darkblue", "lightblue", "white", "red", "darkred"))(100)
scaled_cpm_xb_ordered <- as.matrix(scaled_cpm_CMH_ordered)
annotation_col_xb_ordered$CMH <- as.factor(annotation_col_CMH_ordered$CMH)
pheatmap(scaled_cpm_xb_ordered,  
         annotation_col = annotation_col_xb_ordered,
         annotation_colors = annotation_colors_xb,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 2,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species (Nasal Brush vs Lung Brush)")
#top50 

top50_xb <- rownames(res_xb[order(res_xb$FDR), ])[1:50]
top50_xb <- cpm_xb[top50_xb, , drop = FALSE]
scaled_cpm_xb_top50 <- t(scale(t(top50_xb)))
scaled_cpm_xb_top50[!is.finite(scaled_cpm_xb_top50)] <- 0
annotation_col_nasalvslung_top50 <- data.frame(
  Tissue = XB_metadata[colnames(scaled_cpm_xb_top50), "Tissue"]
)
rownames(annotation_col_nasalvslung_top50) <- colnames(scaled_cpm_xb_top50)
sample_order <- rownames(annotation_col_nasalvslung_top50)[order(annotation_col_nasalvslung_top50$Tissue)]
scaled_cpm_xb_top50 <- scaled_cpm_xb_top50[, sample_order]
annotation_col_nasalvslung_top50 <- annotation_col_nasalvslung_top50[sample_order, , drop = FALSE]
annotation_colors_xb <- list(
  Tissue = c(
    "Lung Brush" = "red",
    "Nasal Brush" = "blue"
  )
)
max_val_top50 <- max(scaled_cpm_xb_top50, na.rm = TRUE)
min_val_top50 <- min(scaled_cpm_xb_top50, na.rm = TRUE)
abs_max_top50 <- max(abs(c(max_val_top50, min_val_top50)))
breaks <- seq(-abs_max_top50, abs_max_top50, length.out = 101)
pheatmap(scaled_cpm_xb_top50,
         annotation_col = annotation_col_nasalvslung_top50,
         annotation_colors = annotation_colors_xb,  
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         color = colorRampPalette(c("darkblue", "lightblue", "white", "red", "darkred"))(100),
         breaks = breaks,
         show_colnames = FALSE,
         fontsize_row = 8,
         main = "Top 50 Significant Taxa (Nasal vs Lung)")

#calculate diversity indices
row.names(XB_rawcounts) <- XB_rawcounts[,2]
XB_rawcounts <- XB_rawcounts[,-1] #x2
otu_xb <- otu_table(XB_rawcounts, taxa_are_rows = TRUE)
row.names(XB_metadata) <- XB_metadata[,1]
XB_metadata <- XB_metadata[,-1]
meta_xb <- sample_data(XB_metadata)
ps3 <- phyloseq(otu_xb, meta_xb)
ps3
alpha_xb <- estimate_richness(ps3, measures = c("Shannon", "Observed", "Simpson", "Chao1"))
merged_alpha_xb <- cbind(XB_metadata[rownames(alpha_xb), ], alpha_xb)
#Observed
ggboxplot(
  merged_alpha_xb,
  x = "Tissue",
  y = "Observed",
  fill = "Tissue",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Lung Brush", "Nasal Brush")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() +theme_pub
#Chao1
ggboxplot(
  merged_alpha_xb,
  x = "Tissue",
  y = "Chao1",
  fill = "Tissue",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Lung Brush", "Nasal Brush")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() +theme_pub
#Shannon
ggboxplot(
  merged_alpha_xb,
  x = "Tissue",
  y = "Shannon",
  fill = "Tissue",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Lung Brush", "Nasal Brush")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() +theme_pub
#Simpson
ggboxplot(
  merged_alpha_xb,
  x = "Tissue",
  y = "Simpson",
  fill = "Tissue",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Lung Brush", "Nasal Brush")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() +theme_pub
#calculate Bray-Curtis beta diversity index
otu_xb_2 <- otu_table(XB_readcounts_1, taxa_are_rows = TRUE)
ps4 <- phyloseq(otu_xb_2, meta_xb)
ps4
bray_curtis_xb <- distance(ps4, method = "bray")
print(bray_curtis_xb)
pcoa_xb <- ordinate(ps4, method = "PCoA", distance = bray_curtis_xb)
print(pcoa_xb)
# Get percent variance explained by PC1 and PC2
var_exp_xb <- round(pcoa_xb$values$Relative_eig[1:2] * 100, 1)
# Create data frame for plotting
plot_data_xb <- plot_ordination(ps4, pcoa_xb, type = "samples", justDF = TRUE)
# Add grouping variable (gender)
plot_data_xb$Tissue <- sample_data(ps4)$Tissue
# Run PERMANOVA (replace "Gender" with your grouping variable)
perm_test_xb <- adonis2(bray_curtis_xb ~ Tissue, data = XB_metadata)
pval_xb <- signif(perm_test_xb$`Pr(>F)`[1], 3)
ggplot(plot_data_xb, aes(x = Axis.1, y = Axis.2, color = Tissue)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) by Tissue"
  ) +
  annotate("text", 
           x = 0.5,  # Change based on your plot scale
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_xb), 
           hjust = 1, size = 3.5) + theme_pub

#taxSEA
taxsea_data_xb <- setNames(res_xb_sig$logFC, res_xb_sig$Species)
taxsea_results_xb <- TaxSEA(taxon_ranks = taxsea_data_xb)
metabolite_xb <- taxsea_results_xb$Metabolite_producers
metabolite_xb <- metabolite_xb %>%
  mutate(logP = -log10(PValue),
         direction = ifelse(median_rank_of_set_members > 0, "Enriched", "Depleted"),
         signed_logP = ifelse(direction == "Enriched", logP, -logP)) %>%
  filter(FDR < 0.05) #FDR > 0.05
health_xb <- taxsea_results_xb$Health_associations
health_xb <- health_xb %>%
  mutate(logP = -log10(PValue),
         direction = ifelse(median_rank_of_set_members > 0, "Enriched", "Depleted"),
         signed_logP = ifelse(direction == "Enriched", logP, -logP)) %>%
  filter(FDR < 0.05) #FDR > 0.05
bugsig_xb <- taxsea_results_xb$BugSigDB
bugsig_xb <- bugsig_xb %>%
  mutate(logP = -log10(PValue),
         direction = ifelse(median_rank_of_set_members > 0, "Enriched", "Depleted"),
         signed_logP = ifelse(direction == "Enriched", logP, -logP)) %>%
  filter(FDR < 0.05) #FDR > 0.05


#paired test
paired_samples <- XB_metadata %>%
  group_by(Patient) %>%
  summarise(tissue_count = n_distinct(Tissue)) %>%
  filter(tissue_count >= 2) %>%  # Keep only patients with BOTH nasal and bronchial
  pull(Patient)
metadata_paired <- XB_metadata %>%
  filter(Patient %in% paired_samples)
row.names(XB_readcounts) <- XB_readcounts[,2]
nasal_samples <- metadata_paired %>% filter(Tissue == "Nasal Brush") %>% arrange(Patient)
lung_samples  <- metadata_paired %>% filter(Tissue == "Lung Brush")  %>% arrange(Patient)

install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)
# Create group variable (Tissue) and blocking variable (Patient)
group <- factor(metadata_paired$Tissue)
patient <- factor(metadata_paired$Patient)
# Create DGEList from raw counts
dge <- DGEList(counts = XB_rawcounts_1)
dge <- calcNormFactors(dge)  # TMM normalization
paired_analysis_mat <- model.matrix(~ patient + group)
dge <- estimateDisp(dge, paired_analysis_mat)
fit <- glmFit(dge, paired_analysis_mat)
# Compare Lung vs Nasal (last coefficient)
lrt <- glmLRT(fit, coef = "groupNasal Brush")
res <- topTags(lrt, n = Inf)$table
res$FDR <- p.adjust(res$PValue, method = "fdr")
sig_taxa <- rownames(res)[res$FDR < 0.05]
length(sig_taxa)  # see how many significant taxa you have
# Extract CPM values
CPM <- cpm(dge, log = FALSE)
# Subset logCPM matrix to keep only significant taxa
CPM_sig <- CPM[sig_taxa,,drop = FALSE ]

common_samples <- intersect(colnames(CPM_sig), rownames(metadata_paired))
CPM_sig <- CPM_sig[, common_samples]
nasal_samples <- metadata_paired %>% filter(Tissue == "Nasal Brush") %>% arrange(Patient)
lung_samples  <- metadata_paired %>% filter(Tissue == "Lung Brush")  %>% arrange(Patient)
nasal_abund <- CPM_sig[, row.names(nasal_samples)]
lung_abund  <- CPM_sig[, row.names(lung_samples)]
# Add small value to avoid log(0)
epsilon <- 1e-6
logFC_per_patient <- log2((nasal_abund + epsilon) / (lung_abund + epsilon))
abs_max <- max(abs(logFC_per_patient), na.rm = TRUE)
breaks <- seq(-abs_max, abs_max, length.out = 101)
pheatmap(logFC_per_patient,
         cluster_rows = TRUE,         
         cluster_cols = TRUE,         
         color = custom_colors,
         breaks = breaks,
         fontsize_row = 2,
         fontsize_col = 8,
         show_colnames = FALSE,       
         show_rownames = TRUE,
         main = "log2FC (Nasal vs Lung) per Patient")
#top50
top50_taxa <- rownames(res[order(res$FDR), ])[1:50]
print(top50_taxa)
CPM_top50 <- CPM[top50_taxa, , drop = FALSE]
nasal_abund_top50 <- CPM_top50[, rownames(nasal_samples)]
lung_abund_top50  <- CPM_top50[, rownames(lung_samples)]
logFC_top50 <- log2((nasal_abund_top50 + epsilon) / (lung_abund_top50 + epsilon))
abs_max <- max(abs(logFC_top50), na.rm = TRUE)
breaks <- seq(-abs_max, abs_max, length.out = 101)
custom_colors <- colorRampPalette(c("darkblue","lightblue", "white", "red","darkred"))(100)
pheatmap(logFC_top50,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = custom_colors,
         breaks = breaks,
         fontsize_row = 8,
         show_colnames = FALSE,
         show_rownames = TRUE,
         main = "Top 50 Taxa: log2FC (Nasal vs Lung) per Patient")


#load 16S data
readcount_16S <- read.delim("20250726_16S_readcounts2.txt", header = TRUE, stringsAsFactors = FALSE)
sample_16S <- read.delim("20250726_16S_samples.txt", header = TRUE, stringsAsFactors = FALSE)

#species level
readcount_16S_1 <- readcount_16S[,-1]
readcount_16S_1 <- readcount_16S %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")
species_16S <- readcount_16S_1$Species
species_xb <- XB_readcounts$Name
#overlapping species
overlap_species <- intersect(species_16S, species_xb)
print(overlap_species)
unique_to_16S <- setdiff(species_16S, species_xb)
print(unique_to_16S)
unique_to_xb <- setdiff(species_xb, species_16S)
print(unique_to_xb)
# Write overlapping species
write.table(overlap_species, "overlap_species.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# Write species only in 16S
write.table(unique_to_16S, "unique_to_16S.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# Write species only in XB
write.table(unique_to_xb, "unique_to_xb.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

install.packages("VennDiagram")
library(VennDiagram)
grid.newpage() 
venn.plot_species <- draw.pairwise.venn(
  area1 = length(unique_to_16S) + length(overlap_species),
  area2 = length(unique_to_xb) + length(overlap_species),
  cross.area = length(overlap_species),
  category = c("16S rRNA sequencing", "Metatranscriptomic sequencing"),
  fill = c("pink", "lightgreen"),
  alpha = 0.7,
  cex = 2,
  cat.cex = 1.5,
  cat.col = c("black", "black"),
  cat.pos = c(-30, 30),        
  cat.dist = c(-0.05, -0.08)   
)
packageVersion("VennDiagram")
#genus level
genus_16S <- word(species_16S, 1)
print(genus_16S)
genus_xb <- word(species_xb, 1)
print(genus_xb)
unique_genus_16S <- unique(genus_16S)
print(unique_genus_16S)
unique_genus_xb <- unique(genus_xb)
print(unique_genus_xb)
#overlap genus
overlap_genus <- intersect(unique_genus_16S, unique_genus_xb)
unique_to_16S_genus <- setdiff(unique_genus_16S, unique_genus_xb)
print(unique_to_16S)
unique_to_xb_genus <- setdiff(unique_genus_xb, unique_genus_16S)
print(unique_to_xb)
# Write overlapping genus
write.table(overlap_genus, "overlap_genus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# Write species only in 16S
write.table(unique_to_16S_genus, "unique_to_16S_genus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# Write species only in XB
write.table(unique_to_xb_genus, "unique_to_xb_genus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
venn.plot_genus <- draw.pairwise.venn(
  area1 = length(unique_to_16S_genus) + length(overlap_genus),
  area2 = length(unique_to_xb_genus) + length(overlap_genus),
  cross.area = length(overlap_genus),
  category = c("16S rRNA sequencing", "Metatranscriptomic sequencing"),
  fill = c("pink", "lightgreen"),
  alpha = 0.7,
  cex = 2,
  cat.cex = 1.5,
  cat.col = c("black", "black"),
  cat.pos = c(-30, 30),       
  cat.dist = c(-0.05, -0.08)   
)
#number of taxa resolved to species level in 16S data
resolved_species_16S <- species_16S[!grepl("sp$", species_16S)]
print(resolved_species_16S)
#number of taxa resolved to species level in metagenomic data
resolved_species_xb <- species_xb[!grepl("sp$", species_xb)]
print(resolved_species_xb)

#calculate proportion of reads assign to species level in 16S data
readcount_16S_1 <- as.data.frame(readcount_16S_1)
rownames(readcount_16S_1) <- readcount_16S_1[[1]]
readcount_16S_1 <- readcount_16S_1[,-1]
# Assumes species-level names include a space (e.g., "Escherichia coli") and not "sp"
species_level_16S <- grepl(" ", rownames(readcount_16S_1)) & !grepl(" sp$", rownames(readcount_16S_1))
# Calculate total and species-level reads per sample
total_reads_per_sample_16S <- colSums(readcount_16S_1)
species_reads_per_sample_16S <- colSums(readcount_16S_1[species_level_16S, ])
# Calculate proportion per sample
proportion_species_level_16S <- species_reads_per_sample_16S / total_reads_per_sample_16S * 100  # in percentage
#Summary stats
mean_prop_16S <- mean(proportion_species_level_16S)
sd_prop_16S <- sd(proportion_species_level_16S)
cat("Proportion of reads resolved to species level per sample in 16S data (%):\n")
print(round(proportion_species_level_16S, 2))
cat("\nMean proportion:", round(mean_prop_16S, 2), "%\n")
cat("Standard deviation:", round(sd_prop_16S, 2), "%\n")

#calculate proportion of reads assign to species level in metatranscriptomic data
species_level_xb <- grepl(" ", rownames(XB_rawcounts_1)) & !grepl(" sp$", rownames(XB_rawcounts))
# Calculate total and species-level reads per sample
total_reads_per_sample_xb <- colSums(XB_rawcounts)
species_reads_per_sample_xb <- colSums(XB_rawcounts[species_level_xb, ])
# Calculate proportion per sample
proportion_species_level_xb <- species_reads_per_sample_xb / total_reads_per_sample_xb * 100  # in percentage
#Summary stats
mean_prop_xb <- mean(proportion_species_level_xb)
sd_prop_xb <- sd(proportion_species_level_xb)
cat("Proportion of reads resolved to species level per sample in metatranscriptomic data (%):\n")
print(round(proportion_species_level_xb, 2))
cat("\nMean proportion:", round(mean_prop_xb, 2), "%\n")
cat("Standard deviation:", round(sd_prop_xb, 2), "%\n")

metadata_16S <- read.delim("20250726_16S_samples2.txt", header = TRUE, stringsAsFactors = FALSE)
metadata_16S_1 <- na.omit(metadta_16S)
metadata_16S_1 <- metadata_16S_1 %>% mutate(Disease = case_when(
  GOLD == "Not applicable" ~ "Control",
  GOLD %in% c("GOLD 1", "GOLD 2") ~ "Moderate COPD",
  GOLD %in% c("GOLD 3", "GOLD 4") ~ "Severe COPD"
))
metadata_16S_1 <- metadata_16S_1 %>% mutate(Age_Group = case_when(
  Age <= 63 ~ "Young", 
  Age > 63 ~ "Old"
))
readcount_16S_filtered <- readcount_16S_1[, colnames(readcount_16S_1) %in% metadata_16S_1$Sample]
library(phyloseq)
#diversity indices for 16S data
rownames(metadata_16S_1) <- metadata_16S_1[[1]]
metadata_16S_1 <- metadata_16S_1[, -1]
otu_16S <- otu_table(readcount_16S_filtered, taxa_are_rows = TRUE)
meta_16S <- sample_data(metadata_16S_1)
ps_16S <- phyloseq(otu_16S, meta_16S)
ps_16S
#alpha diversity 
alpha_df_16S <- estimate_richness(ps_16S, measures = c("Observed", "Shannon", "Simpson", "Chao1"))
merged_alpha_16S <- cbind(
  metadata_16S_1[rownames(alpha_df_16S), , drop = FALSE],
  alpha_df_16S
)
#Observed
library(ggpubr)
#for gender
ggboxplot(
  merged_alpha_16S,
  x = "Sex",
  y = "Observed",
  fill = "Sex",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("male", "female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#for CMH
ggboxplot(
  merged_alpha_16S,
  x = "CMH",
  y = "Observed",
  fill = "CMH",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#age group
ggboxplot(
  merged_alpha_16S,
  x = "Age_Group",
  y = "Observed",
  fill = "Age_Group",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Young", "Old")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#disease group
ggboxplot(
  merged_alpha_16S,
  x = "Disease",
  y = "Observed",
  fill = "Disease",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha_16S$Observed, na.rm = TRUE) + 150  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha_16S$Disease), 2, simplify = FALSE),
    label = "p.signif",       
    hide.ns = FALSE, 
    exact = FALSE
  ) +
  theme_minimal() + theme_pub +
  theme(
    axis.text.x = element_text(size = 7) 
  )
#Chao1
#for gender
ggboxplot(
  merged_alpha_16S,
  x = "Sex",
  y = "Chao1",
  fill = "Sex",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("male", "female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#for CMH
ggboxplot(
  merged_alpha_16S,
  x = "CMH",
  y = "Chao1",
  fill = "CMH",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#age group
ggboxplot(
  merged_alpha_16S,
  x = "Age_Group",
  y = "Chao1",
  fill = "Age_Group",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Young", "Old")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#disease group
ggboxplot(
  merged_alpha_16S,
  x = "Disease",
  y = "Chao1",
  fill = "Disease",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha_16S$Observed, na.rm = TRUE) + 150  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha_16S$Disease), 2, simplify = FALSE),
    label = "p.signif",       
    hide.ns = FALSE, 
    exact = FALSE
  ) +
  theme_minimal() + theme_pub +
  theme(
    axis.text.x = element_text(size = 7) 
  )
#Shannon
#for gender
ggboxplot(
  merged_alpha_16S,
  x = "Sex",
  y = "Shannon",
  fill = "Sex",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("male", "female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#for CMH
ggboxplot(
  merged_alpha_16S,
  x = "CMH",
  y = "Shannon",
  fill = "CMH",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#age group
ggboxplot(
  merged_alpha_16S,
  x = "Age_Group",
  y = "Shannon",
  fill = "Age_Group",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Young", "Old")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#disease group
ggboxplot(
  merged_alpha_16S,
  x = "Disease",
  y = "Shannon",
  fill = "Disease",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha_16S$Shannon, na.rm = TRUE) + 3  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha_16S$Disease), 2, simplify = FALSE),
    label = "p.signif",       
    hide.ns = FALSE, 
    exact = FALSE
  ) +
  theme_minimal() + theme_pub +
  theme(
    axis.text.x = element_text(size = 7) 
  )
#Simpson
#for gender
ggboxplot(
  merged_alpha_16S,
  x = "Sex",
  y = "Simpson",
  fill = "Sex",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("male", "female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#for CMH
ggboxplot(
  merged_alpha_16S,
  x = "CMH",
  y = "Simpson",
  fill = "CMH",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("TRUE", "FALSE")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#age group
ggboxplot(
  merged_alpha_16S,
  x = "Age_Group",
  y = "Simpson",
  fill = "Age_Group",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Young", "Old")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#disease group
ggboxplot(
  merged_alpha_16S,
  x = "Disease",
  y = "Simpson",
  fill = "Disease",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha_16S$Shannon, na.rm = TRUE) + 0.1  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha_16S$Disease), 2, simplify = FALSE),
    label = "p.signif",       
    hide.ns = FALSE, 
    exact = FALSE
  ) +
  theme_minimal() + theme_pub +
  theme(
    axis.text.x = element_text(size = 7) 
  )
#beta diversity
bray_curtis_16S <- distance(ps_16S, method="bray")
library(vegan)
print(bray_curtis_16S)
pcoa_16S <- ordinate(ps_16S, method = "PCoA", distance = bray_curtis_16S)
print(pcoa_16S)
# Get percent variance explained by PC1 and PC2
var_exp_16S <- round(pcoa_16S$values$Relative_eig[1:2] * 100, 1)
# Create data frame for plotting
plot_data_16S <- plot_ordination(ps_16S, pcoa_16S, type = "samples", justDF = TRUE)
# Add grouping variable (gender)
plot_data_16S$Sex <- sample_data(ps_16S)$Sex
# Run PERMANOVA (replace "Gender" with your grouping variable)

perm_test_gender_16S <- adonis2(bray_curtis_16S ~ Sex, data = metadata_16S_1)
pval_gender_16S <- signif(perm_test_gender_16S$`Pr(>F)`[1], 3)
ggplot(plot_data_16S, aes(x = Axis.1, y = Axis.2, color = Sex)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) in 16S data by Gender"
  ) +
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_gender_16S), 
           hjust = 1, size = 3.5) + theme_pub
#CMH
plot_data_16S$CMH <- sample_data(ps_16S)$CMH
perm_test_CMH_16S <- adonis2(bray_curtis_16S ~ CMH, data = metadata_16S_1)
pval_CMH_16S <- signif(perm_test_CMH_16S$`Pr(>F)`[1], 3)
ggplot(plot_data_16S, aes(x = Axis.1, y = Axis.2, color = CMH)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) in 16S data by CMH"
  ) +
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_CMH_16S), 
           hjust = 1, size = 3.5) + theme_pub
#Age group
plot_data_16S$Age_Group <- sample_data(ps_16S)$Age_Group
perm_test_age_16S <- adonis2(bray_curtis_16S ~ Age_Group, data = metadata_16S_1)
pval_age_16S <- signif(perm_test_age_16S$`Pr(>F)`[1], 3)
ggplot(plot_data_16S, aes(x = Axis.1, y = Axis.2, color = Age_Group)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) in 16S data by age group"
  ) +
  annotate("text", 
           x = 0.5, 
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_age_16S), 
           hjust = 1, size = 3.5) + theme_pub
# Add grouping variable (disease group)
plot_data_16S$Disease <- sample_data(ps_16S)$Disease
perm_test_disease_16S <- adonis2(bray_curtis_16S ~ Disease, data = metadata_16S_1)
pval_disease_16S <- signif(perm_test_disease_16S$`Pr(>F)`[1], 3)
ggplot(plot_data_16S, aes(x = Axis.1, y = Axis.2, color = Disease)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +  # 95% confidence ellipses
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) in 16S data by Disease "
  ) +
  annotate("text", 
           x = 0.5,  # Change based on your plot scale
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_disease_16S), 
           hjust = 1, size = 3.5) + theme_pub

#pairwise permanova
# Set the grouping variable
group_levels_16S <- unique(metadata_16S_1$Disease)
# Store results
pairwise_results_16S <- list()
# Loop through all pairwise combinations
for (i in 1:(length(group_levels_16S)-1)) {
  for (j in (i+1):length(group_levels_16S)) {
    # Define the groups
    group1 <- group_levels_16S[i]
    group2 <- group_levels_16S[j]
    # Subset metadata for just those two groups
    subset_metadata <- metadata_16S_1[metadata_16S_1$Disease %in% c(group1, group2), , drop = FALSE]
    # Subset the distance matrix accordingly
    subset_dist <- as.dist(as.matrix(bray_curtis_16S)[rownames(subset_metadata), rownames(subset_metadata)])
    # Run adonis2
    result <- adonis2(subset_dist ~ Disease, data = subset_metadata)
    pairwise_results_16S[[paste(group1, "vs", group2)]] <- result
  }
}
# Print all results
pairwise_results_16S

#correlation test between 16S and nasal metatranscriptomic data
#species
readcount_16S_common <- readcount_16S_1[overlap_species, ]
readcount_meta_common_nasal <- NB_rawcounts[overlap_species, ]
# Calculate total counts per species
readcount_16S_common <- as.data.frame(readcount_16S_common) %>%
  mutate(total_count = rowSums(across(where(is.numeric))))
total_read_species_16S <- readcount_16S_common$total_count
readcount_meta_common_nasal <- readcount_meta_common_nasal %>%
  mutate(total_count = rowSums(across(where(is.numeric))))
total_read_species_meta_nasal <- readcount_meta_common_nasal$total_count
# Calculate Spearman correlation
correlation_nasal <- cor.test(total_read_species_16S, total_read_species_meta_nasal, method = "spearman", exact = FALSE)
print(correlation_nasal)
# Prepare log-transformed values for plotting
log_16S_species <- log10(total_read_species_16S + 1)
log_meta_nasal_species <- log10(total_read_species_meta_nasal + 1)
# Scatter plot
plot(log_16S_species, log_meta_nasal_species,
     xlab = "16S total counts (log10)",
     ylab = "Metatranscriptomic total counts (log10)",
     main = "Correlation between 16S and nasal metatranscriptomic data",
     pch = 19, col = rgb(0, 0, 1, 0.5))
# Fit linear regression line and add it in blue
fit <- lm(log_meta_nasal_species ~ log_16S_species)
abline(fit, col = "blue", lwd = 2)
# Add correlation text (Spearman)
text(x = max(log_16S_species) * 0.7,
     y = max(log_meta_nasal_species) * 0.9,
     labels = paste("Spearman rho =", round(correlation_nasal$estimate, 3), "\n", "p value =", round(correlation_nasal$p.value, 3) ),
     cex = 1.2, col = "black")

#genus
readcount_16S_1_genus <- readcount_16S_1 %>%
  rownames_to_column("Species") %>%
  mutate(Genus = sub(" .*", "", Species), .after = Species)
readcount_16S_genus <- readcount_16S_1_genus %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")
NB_rawcounts_genus <- NB_rawcounts %>%
  rownames_to_column("Species") %>%
  mutate(Genus = sub(" .*", "", Species), .after = Species)
readcount_NB_genus <- NB_rawcounts_genus %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")
# Filter to overlapping genera
readcount_16S_genus_common <- filter(readcount_16S_genus, Genus %in% overlap_genus)
readcount_NB_genus_common  <- filter(readcount_NB_genus, Genus %in% overlap_genus)
# Order the metatranscriptomic data to match the 16S genus order
readcount_NB_genus_common <- readcount_NB_genus_common %>%
  arrange(match(Genus, readcount_16S_genus_common$Genus))
# Calculate total counts per genus (sum across numeric columns)
readcount_16S_genus_common <- readcount_16S_genus_common %>%
  mutate(total_count = rowSums(across(where(is.numeric))))
readcount_NB_genus_common <- readcount_NB_genus_common %>%
  mutate(total_count = rowSums(across(where(is.numeric))))
# Extract total counts vectors aligned by genus
total_read_genus_16S <- readcount_16S_genus_common$total_count
total_read_genus_NB  <- readcount_NB_genus_common$total_count
# Calculate Spearman correlation
correlation_nasal_genus <- cor.test(total_read_genus_16S, total_read_genus_NB, method = "spearman", exact = FALSE)
print(correlation_nasal_genus)
# Prepare log-transformed values for plotting
log_16S_genus <- log10(total_read_genus_16S + 1)
log_NB_genus <- log10(total_read_genus_NB + 1)
# Plot correlation
plot(log_16S_genus, log_NB_genus,
     xlab = "16S total counts (log10)",
     ylab = "Metatranscriptomic total counts (log10)",
     main = "Correlation between 16S and nasal metatranscriptomic data (genus level)",
     pch = 19, col = rgb(0, 0, 1, 0.5))
# Fit and add regression line (blue)
fit_genus <- lm(log_NB_genus ~ log_16S_genus)
abline(fit_genus, col = "blue", lwd = 2)
# Add correlation text with Spearman rho label
text(x = max(log_16S_genus) * 0.7,
     y = max(log_NB_genus) * 0.9,
     labels = paste("Spearman rho =", round(correlation_nasal_genus$estimate, 3), "\n", "p value =", correlation_nasal_genus$p.value),
     cex = 1.2, col = "black")

#correlation between 16S and bronchial metatranscriptomic data
BB_metadata <- XB_metadata %>%
  filter(Tissue == "Lung Brush")
BB_sample_ids <- row.names(BB_metadata)
BB_rawcounts <- XB_rawcounts[,BB_sample_ids]
#species
readcount_16S_common <- readcount_16S_1[overlap_species, ]
readcount_meta_common_bronchial <- BB_rawcounts[overlap_species,]
readcount_16S_common <- as.data.frame(readcount_16S_common) %>%
  mutate(total_count = rowSums(across(where(is.numeric))))
total_read_species_16S <- readcount_16S_common$total_count
readcount_meta_common_bronchial <- readcount_meta_common_bronchial %>%
  mutate(total_count = rowSums(across(where(is.numeric))))
total_read_species_meta_bronchial <- readcount_meta_common_bronchial$total_count
correlation_bronchial <- cor.test(total_read_species_16S, total_read_species_meta_bronchial, method = "spearman", exact = FALSE)
print(correlation_bronchial)
# Prepare log-transformed values for plotting
log_16S_species <- log10(total_read_species_16S + 1)
log_BB_species <- log10(total_read_species_meta_bronchial + 1)
plot(log_16S_species, log_BB_species,
     xlab = "16S total counts (log10)",
     ylab = "Metatranscriptomic total counts (log10)",
     main = "Correlation between 16S and bronchial metatranscriptomic data",
     pch = 19, col = rgb(0, 0, 1, 0.5))
# Fit and add regression line (blue)
fit_bronchial_species <- lm(log_BB_species ~ log_16S_species)
abline(fit_bronchial_species, col = "blue", lwd = 2)
# Add correlation text with Spearman rho label
text(x = max(log_16S_species) * 0.7,
     y = max(log_BB_species) * 0.9,
     labels = paste("Spearman rho =", round(correlation_bronchial$estimate, 3), "\n", "p value =", round(correlation_bronchial$p.value, 3) ),
     cex = 1.2, col = "black")


#genus
readcount_16S_genus <- readcount_16S_1 %>%
  rownames_to_column("Species") %>%
  mutate(Genus = sub(" .*", "", Species), .after = Species)
readcount_16S_genus <- readcount_16S_1 %>%
  rownames_to_column("Species") %>%
  mutate(Genus = sub(" .*", "", Species), .after = Species) %>%  # extract genus
  select(-Species) %>%  # remove Species column
  group_by(Genus) %>% 
  summarise(across(everything(), sum)) %>%
  ungroup()
BB_rawcounts_genus <- BB_rawcounts %>%
  rownames_to_column("Species") %>%
  mutate(Genus = sub(" .*", "", Species), .after = Species) %>%  # extract genus
  select(-Species) %>%  # remove Species column
  group_by(Genus) %>% 
  summarise(across(everything(), sum)) %>%
  ungroup()
readcount_16S_genus_common <- readcount_16S_genus %>%
  filter(Genus %in% overlap_genus) %>%
  arrange(Genus)
readcount_BB_genus_common <- BB_rawcounts_genus %>%
  filter(Genus %in% overlap_genus) %>%
  arrange(Genus)
# Calculate total counts
readcount_16S_genus_common <- readcount_16S_genus_common %>%
  mutate(total_count = rowSums(across(where(is.numeric))))
readcount_BB_genus_common <- readcount_BB_genus_common %>%
  mutate(total_count = rowSums(across(where(is.numeric))))
# Extract aligned vectors
total_read_genus_16S <- readcount_16S_genus_common$total_count
print(total_read_genus_16S)
total_read_genus_BB <- readcount_BB_genus_common$total_count
print(total_read_genus_BB) 
# Now correlation will work
correlation_bronchial_genus <- cor.test(total_read_genus_16S, total_read_genus_BB, method = "spearman", exact = FALSE)
print(correlation_bronchial_genus)
# Prepare log-transformed values for plotting
log_16S_genus <- log10(total_read_genus_16S + 1)
log_BB_genus <- log10(total_read_genus_BB + 1)
plot(log_16S_genus, log_BB_genus,
     xlab = "16S total counts (log10)",
     ylab = "Metatranscriptomic total counts (log10)",
     main = "Correlation between 16S and bronchial metatranscriptomic data (genus level)",
     pch = 19, col = rgb(0, 0, 1, 0.5))
# Fit and add regression line (blue)
fit_bronchial_genus <- lm(log_BB_genus ~ log_16S_genus)
abline(fit_bronchial_genus, col = "blue", lwd = 2)
# Add correlation text with Spearman rho label
text(x = max(log_16S_genus) * 0.7,
     y = max(log_BB_genus) * 0.9,
     labels = paste("Spearman rho =", round(correlation_bronchial_genus$estimate, 3), "\n", "p value =", round(correlation_bronchial_genus$p.value, 3)),
     cex = 1.2, col = "black")




#make histogram
readcount_16S_df <- readcount_16S_1 %>%
  as.data.frame() %>%
  rownames_to_column("Species") %>%
  mutate(Genus = sub(" .*", "", Species), .after =  Species)
meta_counts_df <- XB_rawcounts %>%
  as.data.frame() %>%
  rownames_to_column("Species") %>%
  mutate(Genus = sub(" .*", "", Species), .after = Species)
genus_16S <- readcount_16S_df %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() %>%
  mutate(total_count = rowSums(across(where(is.numeric))))
genus_meta <- meta_counts_df %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() %>%
  mutate(total_count = rowSums(across(where(is.numeric))))
overlap_counts <- genus_16S %>%
  filter(Genus %in% overlap_genus) %>%
  pull(total_count)
only_16S_counts <- genus_16S %>%
  filter(Genus %in% unique_to_16S_genus) %>%
  pull(total_count)
max_count <- max(c(overlap_counts, only_16S_counts))

#make histogram
# Set plotting layout: 1 row, 2 columns
par(mfrow = c(1, 2))
# Histogram for overlapping genera
hist(log10(overlap_counts + 1),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     main = "genera found in 16S rRNA 
   and metatranscriptomic data",
     xlab = "Log10 Total read counts per genus",
     ylab = "Frequency",
     ylim = c(0, 25),
     xlim = c(0, 7))

# Histogram for 16S-only genera
hist(log10(only_16S_counts + 1),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     main = "genera found only in 16S rRNA data ",
     xlab = "Log10 Total read counts per genus",
     ylab = "Frequency",
     ylim = c(0, 25),
     xlim = c(0, 7))

#analysis for 16S
group_diseasegroup_16S <- factor(metadata_16S_1$Disease, levels = c("Control", "Moderate COPD", "Severe COPD")) #Set up group variable
dge_diseasegroup_16S <- DGEList(counts = readcount_16S_filtered, group = group_diseasegroup_16S) #Create DGEList
dge_diseasegroup_16S <- calcNormFactors(dge_diseasegroup_16S) #TMM normalization
design_diseasegroup_16S <- model.matrix(~ 0 + group_diseasegroup_16S) # Design matrix
colnames(design_diseasegroup_16S) <- make.names(colnames(design_diseasegroup_16S))
colnames(design_diseasegroup_16S)
dge_diseasegroup_16S <- estimateDisp(dge_diseasegroup_16S, design_diseasegroup_16S) # Estimate dispersions
fit_diseasegroup_16S <- glmFit(dge_diseasegroup_16S, design_diseasegroup_16S) # Fit GLM model
#Create contrasts for pairwise comparisons
contrast_matrix_16S <- makeContrasts(
  Mod_vs_Ctrl = group_diseasegroup_16SModerate.COPD - group_diseasegroup_16SControl,
  Sev_vs_Ctrl = group_diseasegroup_16SSevere.COPD - group_diseasegroup_16SControl,
  Sev_vs_Mod  = group_diseasegroup_16SSevere.COPD - group_diseasegroup_16SModerate.COPD,
  levels = design_diseasegroup_16S
) # Define contrast matrix using valid names
#Run glmLRT for each pair
lrt_mod_ctrl_16S <- glmLRT(fit_diseasegroup_16S, contrast = contrast_matrix_16S[,"Mod_vs_Ctrl"])
lrt_sev_ctrl_16S <- glmLRT(fit_diseasegroup_16S, contrast = contrast_matrix_16S[,"Sev_vs_Ctrl"])
lrt_sev_mod_16S  <- glmLRT(fit_diseasegroup_16S, contrast = contrast_matrix_16S[,"Sev_vs_Mod"])
#Extract result tables
res_mod_ctrl_16S <- topTags(lrt_mod_ctrl_16S, n = Inf)$table
res_mod_ctrl_16S$FDR <- p.adjust(res_mod_ctrl_16S$PValue, method = "fdr")

res_sev_ctrl_16S <- topTags(lrt_sev_ctrl_16S, n = Inf)$table
res_sev_ctrl_16S$FDR <- p.adjust(res_sev_ctrl_16S$PValue, method = "fdr")

res_sev_mod_16S <- topTags(lrt_sev_mod_16S, n = Inf)$table
res_sev_mod_16S$FDR <- p.adjust(res_sev_mod_16S$PValue, method = "fdr")
#View number of significant taxa in each pairwise comparison
cat("Significant taxa (Moderate vs Control):", sum(res_mod_ctrl_16S$FDR < 0.05), "\n")
cat("Significant taxa (Severe vs Control):", sum(res_sev_ctrl_16S$FDR < 0.05), "\n")
cat("Significant taxa (Severe vs Moderate):", sum(res_sev_mod_16S$FDR < 0.05), "\n")

sig_taxa_mod_ctrl_16S <- rownames(res_mod_ctrl_16S)[res_mod_ctrl_16S$FDR < 0.05]
sig_taxa_sev_ctrl_16S <- rownames(res_sev_ctrl_16S)[res_sev_ctrl_16S$FDR < 0.05]
sig_taxa_sev_mod_16S  <- rownames(res_sev_mod_16S)[res_sev_mod_16S$FDR < 0.05]

res_mod_ctrl_16S$Species <- rownames(res_mod_ctrl_16S)
res_sev_ctrl_16S$Species <- rownames(res_sev_ctrl_16S)
res_sev_mod_16S$Species  <- rownames(res_sev_mod_16S)

res_mod_ctrl_sig_16S <- res_mod_ctrl_16S %>% filter(FDR<0.05) #nothing is statistically significant
res_sev_ctrl_sig_16S <- res_sev_ctrl_16S %>% filter(FDR<0.05)
res_sev_mod_sig_16S <- res_sev_mod_16S %>% filter(FDR<0.05)

#top 50 severe vs control
top50_sev_ctrl_16S <- res_sev_ctrl_16S %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
ggplot(top50_sev_ctrl_16S, aes(x = reorder(Species, logFC), y = logFC, fill = logFC > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "blue", "FALSE" = "green"),
    labels = c("FALSE" = "Control", "TRUE" = "Severe COPD"),
    name = "Higher in"
  ) +
  labs(title = "Differentially Abundant Species (Severe vs Control) in 16S rRNA sequencing",
       x = "Species",
       y = "log2 Fold Change (Severe - Control)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # 👈 adjust this size as needed
  )
#volcano plot
res_sev_ctrl_16S <- res_sev_ctrl_16S %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_sev_ctrl_16S, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_sev_ctrl_16S %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: Severe vs Control",
    x = "log2 Fold Change (Severe - Control)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#top 50 severe vs moderate
top50_sev_mod_16S <- res_sev_mod_16S %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
ggplot(top50_sev_mod_16S, aes(x = reorder(Species, logFC), y = logFC, fill = logFC > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "blue", "FALSE" = "pink"),
    labels = c("FALSE" = "Moderate COPD", "TRUE" = "Severe COPD"),
    name = "Higher in"
  ) +
  labs(title = "Differentially Abundant Species (Severe vs Modrate) in 16S rRNA sequencing",
       x = "Species",
       y = "log2 Fold Change (Severe - Moderate)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6), 
  )
#volcano plot
res_sev_mod_16S <- res_sev_mod_16S %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_sev_mod_16S, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_sev_mod_16S %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: Severe vs Moderate",
    x = "log2 Fold Change (Severe - Moderate)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()

#heatmap
lrt_diseasegroup_16S <- glmLRT(fit_diseasegroup_16S)
res_diseasegroup_16S <- topTags(lrt_diseasegroup_16S, n = Inf)$table
res_diseasegroup_16S$FDR <- p.adjust(res_diseasegroup_16S$PValue, method = "fdr")
res_diseasegroup_16S$Species <- rownames(res_diseasegroup_16S)
res_diseasegroup_sig_16S <- res_diseasegroup_16S %>% filter(FDR<0.05)
cpm_diseasegroup_16S <- cpm(dge_diseasegroup_16S, log = FALSE)
sig_species_diseasegroup_16S <- res_diseasegroup_sig_16S$Species
print(sig_species_diseasegroup_16S)
cpm_diseasegroup_sig_16S <- cpm_diseasegroup_16S[sig_species_diseasegroup_16S, , drop = FALSE]
setdiff(sig_species_diseasegroup_16S, rownames(cpm_diseasegroup_16S))
scaled_cpm_diseasegroup_16S <- t(scale(t(cpm_diseasegroup_16S)))
scaled_cpm_diseasegroup_16S[!is.finite(scaled_cpm_diseasegroup_16S)] <- 0
row.names(metadata_16S_1) <- metadata_16S_1[,1]
annotation_col_diseasegroup_16S <- data.frame(
  Disease = metadata_16S_1[colnames(cpm_diseasegroup_16S), "Disease"]
)
rownames(annotation_col_diseasegroup_16S) <- colnames(cpm_diseasegroup_16S)
sample_ordered_disease_16S <- rownames(annotation_col_diseasegroup_16S)[
  order(annotation_col_diseasegroup_16S$Disease)
]
scaled_cpm_diseasegroup_ordered_16S <- scaled_cpm_diseasegroup_16S[, sample_ordered_disease_16S]
annotation_col_diseasegroup_ordered_16S <- annotation_col_diseasegroup_16S[sample_ordered_disease_16S, , drop = FALSE]
annotation_colors_diseasegroup_16S <- list(
  Disease = c("Control" = "green", 
              "Moderate COPD" = "pink", 
              "Severe COPD" = "blue")
)
library(pheatmap)
pheatmap(scaled_cpm_diseasegroup_ordered_16S,  
         annotation_col = annotation_col_diseasegroup_ordered_16S,
         annotation_colors = annotation_colors_diseasegroup_16S,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 3,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species (Disease Group)")
#top50
top50_diseasegroup_16S <- rownames(res_diseasegroup_16S[order(res_diseasegroup_16S$FDR), ])[1:50]
cpm_diseasegroup_top50_16S <- cpm_diseasegroup_16S[top50_diseasegroup_16S, , drop = FALSE]
scaled_cpm_diseasegroup_top50_16S <- t(scale(t(cpm_diseasegroup_top50_16S)))
scaled_cpm_diseasegroup_top50_16S[!is.finite(scaled_cpm_diseasegroup_top50_16S)] <- 0  # handle NaNs
annotation_col_diseasegroup_top50_16S <- data.frame(
  Disease = metadata_16S_1[colnames(cpm_diseasegroup_top50_16S), "Disease"]
)
rownames(annotation_col_diseasegroup_top50_16S) <- colnames(cpm_diseasegroup_top50_16S)
sample_ordered_disease_top50_16S <- rownames(annotation_col_diseasegroup_top50_16S)[
  order(annotation_col_diseasegroup_top50_16S$Disease)
]
scaled_cpm_diseasegroup_top50_ordered_16S <- scaled_cpm_diseasegroup_top50_16S[, sample_ordered_disease_top50_16S]
annotation_col_diseasegroup_top50_ordered_16S <- annotation_col_diseasegroup_top50_16S[sample_ordered_disease_top50_16S, , drop = FALSE]
abs_max <- max(abs(scaled_cpm_diseasegroup_top50_ordered_16S), na.rm = TRUE)
breaks <- seq(-abs_max, abs_max, length.out = 101)
pheatmap(scaled_cpm_diseasegroup_top50_ordered_16S,  
         annotation_col = annotation_col_diseasegroup_top50_ordered_16S,
         annotation_colors = annotation_colors_diseasegroup_16S,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 8,
         color = custom_colors,
         breaks = breaks,
         main = "Top 50 Significant Species (Disease Group) in 16S rRNA data")

#differential analysis for gender
group_gender_16S <- factor(metadata_16S_1$Sex, levels = c("male", "female")) #create group factor
dge_gender_16S <- DGEList(counts = readcount_16S_filtered, group = group_gender_16S)
dge_gender_16S <- calcNormFactors(dge_gender_16S)
design_gender_16S <- model.matrix(~ 0 + group_gender_16S)
dge_gender_16S <- estimateDisp(dge_gender_16S, design_gender_16S)
fit_gender_16S <- glmFit(dge_gender_16S, design_gender_16S)
contrast_gender_16S <- makeContrasts(Female_vs_Male = group_gender_16Sfemale - group_gender_16Smale, levels = design_gender_16S)
lrt_gender_16S <- glmLRT(fit_gender_16S, contrast = contrast_gender_16S)
res_gender_16S <- topTags(lrt_gender_16S, n= Inf)$table
res_gender_16S$FDR <- p.adjust(res_gender_16S$PValue, method = "fdr")
res_gender_16S$Species <- rownames(res_gender_16S)
res_gender_16S <- res_gender_16S %>%
  filter(FDR < 0.05)
#bar chart
res_gender_16S$Higher_in <- ifelse(res_gender_16S$logFC > 0, "Female", "Male")
ggplot(res_gender_16S, aes(x = reorder(Species, logFC), y = logFC, fill = Higher_in)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Female" = "red", "Male" = "blue"),
    name = "Higher in"
  ) +
  labs(title = "Top Differentially Abundant Species (Female vs Male)",
       x = "Species",
       y = "log2 Fold Change (Female - Male)") +
  theme_minimal()
#volcano plot
res_gender_16S <- res_gender_16S %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_gender_16S, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_gender_16S %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: Female vs Male",
    x = "log2 Fold Change (Female - Male)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#extract cpm and make the heatmap
cpm_gender_16S <- cpm(dge_gender_16S, log = FALSE)
sig_species_gender_16S <- res_gender_16S$Species
cpm_gender_sig_16S <- cpm_gender_16S[sig_species_gender_16S, , drop = FALSE]
scaled_cpm_gender_16S <- t(scale(t(cpm_gender_sig_16S)))
annotation_col_gender_16S <- data.frame(
  Sex = metadata_16S_1[colnames(cpm_gender_sig_16S), "Sex"]
)
rownames(annotation_col_gender_16S) <- colnames(cpm_gender_sig_16S)
annotation_colors_gender_16S <- list(
  Sex = c("male" = "blue", "female" = "red")
)
gender_order_16S <- annotation_col_gender_16S$Sex
sample_ordered_16S <- rownames(annotation_col_gender_16S)[order(annotation_col_gender_16S$Sex)]
print(sample_ordered_16S)
scaled_cpm_gender_ordered_16S <- scaled_cpm_gender_16S[, sample_ordered_16S]
annotation_col_gender_ordered_16S <- annotation_col_gender_16S[sample_ordered_16S, , drop = FALSE]

# Calculate symmetric color scale around 0
max_val <- max(scaled_cpm_gender_16S, na.rm = TRUE)
min_val <- min(scaled_cpm_gender_16S, na.rm = TRUE)
abs_max <- max(abs(c(max_val, min_val)))

breaks <- seq(-abs_max, abs_max, length.out = 101)

# Define matching color palette
custom_colors <- colorRampPalette(c("darkblue", "lightblue", "white", "red", "darkred"))(100)

# Plot heatmap
pheatmap(scaled_cpm_gender_ordered_16S,  
         annotation_col = annotation_col_gender_ordered_16S,
         annotation_colors = annotation_colors_gender_16S,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE, 
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 8,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species (Gender)")


#differential analysis for age group
group_age_16S <- factor(metadata_16S_1$Age_Group, levels = c("Young", "Old"))
dge_age_16S <- DGEList(counts = readcount_16S_filtered, group = group_age_16S)
dge_age_16S <- calcNormFactors(dge_age_16S)
design_age_16S <- model.matrix(~ 0 + group_age_16S)
dge_age_16S <- estimateDisp(dge_age_16S, design_age_16S)
fit_age_16S <- glmFit(dge_age_16S, design_age_16S)
contrast_age_16S <- makeContrasts(Old_vs_Young = group_age_16SOld - group_age_16SYoung, levels = design_age_16S)
lrt_age_16S <- glmLRT(fit_age_16S, contrast = contrast_age_16S)
res_age_16S <- topTags(lrt_age_16S, n = Inf)$table
res_age_16S$FDR <- p.adjust(res_age_16S$PValue, method = "fdr")
res_age_16S$Species <- rownames(res_age_16S)
res_age_sig_16S <- res_age_16S %>%
  filter(FDR < 0.05)
top50_age_sig_16S <- res_age_sig_16S %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
top50_age_sig_16S$Higher_in <- ifelse(top50_age_sig_16S$logFC > 0, "Old", "Young")
ggplot(top50_age_sig_16S, aes(x = reorder(Species, logFC), y = logFC, fill = Higher_in)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Old" = "red", "Young" = "blue"),
    name = "Higher in"
  ) +
  labs(title = "Top Differentially Abundant Species (Old vs Young)",
       x = "Species",
       y = "log2 Fold Change (Old - Young)") +
  theme_minimal()
#volcano plot
res_age_16S <- res_age_16S %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_age_16S, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_age_16S %>% arrange(FDR) %>% slice(1:4),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: Old vs Young",
    x = "log2 Fold Change (Old - Young)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#heatmap
cpm_age_16S <- cpm(dge_age_16S, log = FALSE)
sig_species_age_16S <- res_age_16S$Species[res_age_16S$FDR < 0.05]
cpm_age_sig_16S <- cpm_age_16S[sig_species_age_16S, , drop = FALSE]
scaled_cpm_age_16S <- t(scale(t(cpm_age_sig_16S)))
annotation_col_age_16S <- data.frame(
  Age_Group = metadata_16S_1[colnames(cpm_age_sig_16S), "Age_Group"]
)
rownames(annotation_col_age_16S) <- colnames(cpm_age_sig_16S)
annotation_colors_age_16S <- list(
  Age_Group = c("Young" = "blue", "Old" = "red")
)
sample_ordered_age_16S <- rownames(annotation_col_age_16S)[order(annotation_col_age_16S$Age_Group)]
scaled_cpm_age_ordered_16S <- scaled_cpm_age_16S[, sample_ordered_age_16S]
annotation_col_age_ordered_16S <- annotation_col_age_16S[sample_ordered_age_16S, , drop = FALSE]
max_val <- max(scaled_cpm_age_16S, na.rm = TRUE)
min_val <- min(scaled_cpm_age_16S, na.rm = TRUE)
abs_max <- max(abs(c(max_val, min_val)))
breaks <- seq(-abs_max, abs_max, length.out = 101)
custom_colors <- colorRampPalette(c("darkblue", "lightblue", "white", "red", "darkred"))(100)
pheatmap(scaled_cpm_age_ordered_16S,  
         annotation_col = annotation_col_age_ordered_16S,
         annotation_colors = annotation_colors_age_16S,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 8,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species (Age Group)")


#differential analysis for CMH
group_CMH_16S <- factor(metadata_16S_1$CMH, levels = c("FALSE", "TRUE"))
dge_CMH_16S <- DGEList(counts = readcount_16S_filtered, group = group_CMH_16S)
dge_CMH_16S <- calcNormFactors(dge_CMH_16S)
design_CMH_16S <- model.matrix(~ 0 + group_CMH_16S)
dge_CMH_16S <- estimateDisp(dge_CMH_16S, design_CMH_16S)
fit_CMH_16S <- glmFit(dge_CMH_16S, design_CMH_16S)
contrast_CMH_16S <- makeContrasts(
  CMH_TRUE_vs_FALSE = group_CMH_16STRUE - group_CMH_16SFALSE,
  levels = design_CMH_16S
)
lrt_CMH_16S <- glmLRT(fit_CMH_16S, contrast = contrast_CMH_16S)
res_CMH_16S <- topTags(lrt_CMH_16S, n = Inf)$table
res_CMH_16S$FDR <- p.adjust(res_CMH_16S$PValue, method = "fdr")
res_CMH_16S$Species <- rownames(res_CMH_16S)
res_CMH_sig_16S <- res_CMH_16S %>% 
  filter(FDR < 0.05)  #only Staphylococcus haemolyticus is significant

#taxSEA for 16S
#gender
taxsea_data_gender_16S <- setNames(res_gender_16S$logFC, res_gender_16S$Species)
taxsea_results_gender_16S <- TaxSEA(taxon_ranks = taxsea_data_gender_16S)
metabolite_gender_16S <- taxsea_results_gender_16S$Metabolite_producers #nothing is statistically significant

#age group
taxsea_data_agegroup_16S <- setNames(res_age_16S$logFC, res_age_16S$Species)
taxsea_results_agegroup_16S <- TaxSEA(taxon_ranks = taxsea_data_agegroup_16S)
metabolite_agegroup_16S <- taxsea_results_agegroup_16S$Metabolite_producers #nothing is statistically significant
health_agegroup_16S <- taxsea_results_agegroup_16S$Health_associations #nothing is statistically significant
bugsig_agegroup_16S <- taxsea_results_agegroup_16S$BugSigDB #nothing is statistically significant

#for CMH 
taxsea_data_CMH_16S <- setNames(res_CMH_sig_16S$logFC, res_CMH_sig_16S$Species)
taxsea_results_CMH_16S <- TaxSEA(taxon_ranks = taxsea_data_CMH_16S)

#for disease group
#control vs severe COPD
taxsea_data_SvsC_16S <- setNames(res_sev_ctrl_sig_16S$logFC, res_sev_ctrl_sig_16S$Species)
taxsea_results_SvsC_16S <- TaxSEA(taxon_ranks = taxsea_data_SvsC_16S) 
metabolite_SvsC_16S <- taxsea_results_SvsC_16S$Metabolite_producers #FDR > 0.05

#severe vs moderate COPD
taxsea_data_SvsM_16S <- setNames(res_sev_mod_sig_16S$logFC, res_sev_mod_sig_16S$Species)
taxsea_results_SvsM_16S <- TaxSEA(taxon_ranks = taxsea_data_SvsM_16S)
metabolite_SvsM_16S <- taxsea_results_SvsM_16S$Metabolite_producers #FDR > 0.05
