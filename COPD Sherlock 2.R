# Load BLAST results
blast <- read.delim("COPD_result.txt", header = TRUE, stringsAsFactors = FALSE)

#check data structure
head(blast)

colnames(blast) <- c("query ID", "Organism", "% identity", "alignment length", "mismatch", "gap opens", "q.start", "q.end", "s.start", "s.end", "E_Value", "Bit Score")

install.packages("tidyverse")
library(tidyverse)
library(dplyr)

#set mininal identity threshold
min_identity <- 95

filtered_microbes <- blast %>% 
  filter(`% identity` >= min_identity) %>% #filter for identity >= 95
  distinct() #remove duplicates

head(filtered_microbes)

# select rows with highest Bit Score
best_hits <- filtered_microbes %>%
  group_by(`query ID`) %>%
  slice_max(`Bit Score`, with_ties = TRUE) %>%
  ungroup()

print(best_hits)

# Install Bioconductor manager (needed for Biostrings)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Install Biostrings from Bioconductor
BiocManager::install("Biostrings")
library(Biostrings)

database_file <- readDNAStringSet("C:\\Users\\Admin\\Downloads\\all_seqs.fa\\all_seqs.fa")
headers <- names(database_file)
taxonomy_db <- data.frame(
  organism_code = sub(" .*", "", headers),  # Extract bacterial code (before first space)
  organism_name = sub("^[^ ]+ ", "", headers),  # Extract organism name (after first space)
  stringsAsFactors = FALSE
)
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
write.table(filtered_organisms_1, "reference_bacteria.txt", sep = "\t", row.names = FALSE, quote = FALSE)

best_hits_final <- best_hits_final %>% 
  mutate(organism_name_2 = str_extract(organism_name, "^\\S+\\s\\S+"))
write.table(best_hits_final, "best_hits_final", sep = "\t", row.names = FALSE, quote = FALSE)

# Count the number of contigs for each organism
contig_counts <- best_hits_final %>% 
  distinct(`query ID`, organism_name_2) %>% #ensure uniqueness
  group_by(organism_name_2) %>%
  summarise(contig_count = n(), .groups = "drop")
write.table(contig_counts, "contig_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Count the number of bases for each organism
fasta_file <- readDNAStringSet("C:\\Users\\Admin\\Downloads\\COPD\\20250311_sherlock2_brush_trinity.Trinity.fasta", format = "fasta")
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
write.table(base_counts, "base_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

install.packages("ggplot2")
install.packages("pheatmap")
library(ggplot2)
library(pheatmap)

# Load bacterial read count data and metadata
read_count_table <- read.delim("20250325_Sherlock2_brushes_bact_exp.txt", row.names = 1, check.names =  FALSE)
metadata <- read.delim("20250325_Sherlock2_brushes_samples.txt", row.names = 1, check.names = FALSE)

head(read_count_table)
rownames(read_count_table) <- read_count_table[,1]
read_count_table_1 <- read_count_table[,-1]
read_count_table_1 <- read_count_table_1[, colnames(read_count_table_1) %in% rownames(metadata)]

#differential analysis
if (!requireNamespace("edgeR", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("edgeR")
}
library(edgeR)
library(ggpubr) 
library(ggrepel)
#disease group
colnames(raw_count_table) <- gsub("\\.", "-", colnames(raw_count_table))
metadata$Group[metadata$Group == "MildModerateCOPD"] <- "Moderate COPD"
metadata$Group[metadata$Group == "SevereCOPD"] <- "Severe COPD"
group_diseasegroup <- factor(metadata$Group, levels = c("Control", "Moderate COPD", "Severe COPD")) #Set up group variable
dge_diseasegroup <- DGEList(counts = raw_count_table, group = group_diseasegroup) #Create DGEList
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
  levels = design_diseasegroup) # Define contrast matrix using valid names
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
  Group = metadata[colnames(cpm_diseasegroup), "Group"]
)
rownames(annotation_col_diseasegroup) <- colnames(cpm_diseasegroup)
sample_ordered_disease <- rownames(annotation_col_diseasegroup)[
  order(annotation_col_diseasegroup$Group)
]
scaled_cpm_diseasegroup_ordered <- scaled_cpm_diseasegroup[, sample_ordered_disease]
annotation_col_diseasegroup_ordered <- annotation_col_diseasegroup[sample_ordered_disease, , drop = FALSE]
annotation_colors_diseasegroup <- list(
  Group = c("Control" = "green", 
              "Moderate COPD" = "pink", 
              "Severe COPD" = "blue")
)
abs_max <- max(abs(scaled_cpm_diseasegroup_ordered), na.rm = TRUE)
breaks <- seq(-abs_max, abs_max, length.out = 101)
pheatmap(scaled_cpm_diseasegroup_ordered,  
         annotation_col = annotation_col_diseasegroup_ordered,
         annotation_colors = annotation_colors_diseasegroup,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 1,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species (Disease Group)")
#top50
top50_diseasegroup <- rownames(res_diseasegroup[order(res_diseasegroup$FDR), ])[1:50]
cpm_diseasegroup_top50 <- cpm_diseasegroup[top50_diseasegroup, , drop = FALSE]
scaled_cpm_diseasegroup_top50 <- t(scale(t(cpm_diseasegroup_top50)))
scaled_cpm_diseasegroup_top50[!is.finite(scaled_cpm_diseasegroup_top50)] <- 0  # handle NaNs
annotation_col_diseasegroup_top50 <- data.frame(
  Group = metadata[colnames(cpm_diseasegroup_top50), "Group"]
)
rownames(annotation_col_diseasegroup_top50) <- colnames(cpm_diseasegroup_top50)
sample_ordered_disease_top50 <- rownames(annotation_col_diseasegroup_top50)[
  order(annotation_col_diseasegroup_top50$Group)
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
group_gender <- factor(metadata$SEX, levels = c("Male", "Female")) #create group factor
dge_gender <- DGEList(counts = raw_count_table, group = group_gender)
dge_gender <- calcNormFactors(dge_gender)
design_gender <- model.matrix(~ 0 + group_gender)
dge_gender <- estimateDisp(dge_gender, design_gender)
fit_gender <- glmFit(dge_gender, design_gender)
contrast_gender <- makeContrasts(Female_vs_Male = group_genderFemale - group_genderMale, levels = design_gender)
lrt_gender <- glmLRT(fit_gender, contrast = contrast_gender)
res_gender <- topTags(lrt_gender, n= Inf)$table
res_gender$FDR <- p.adjust(res_gender$PValue, method = "fdr")
res_gender$Species <- rownames(res_gender)
res_gender_sig <- res_gender %>%
  filter(FDR < 0.05)
#bar chart
res_gender_sig$Higher_in <- ifelse(res_gender_sig$logFC > 0, "Female", "Male")
ggplot(res_gender_sig, aes(x = reorder(Species, logFC), y = logFC, fill = Higher_in)) +
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
sig_species_gender <- res_gender_sig$Species
cpm_gender_sig <- cpm_gender[sig_species_gender, , drop = FALSE]
scaled_cpm_gender <- t(scale(t(cpm_gender_sig)))
annotation_col_gender <- data.frame(
  SEX = metadata[colnames(cpm_gender_sig), "SEX"]
)
rownames(annotation_col_gender) <- colnames(cpm_gender_sig)
annotation_colors_gender <- list(
  SEX = c("Male" = "blue", "Female" = "red")
)
gender_order <- annotation_col_gender$SEX
sample_ordered <- rownames(annotation_col_gender)[order(annotation_col_gender$SEX)]
scaled_cpm_gender_ordered <- scaled_cpm_gender[, sample_ordered]
annotation_col_gender_ordered <- annotation_col_gender[sample_ordered, , drop = FALSE]

# Calculate symmetric color scale around 0
max_val <- max(scaled_cpm_gender, na.rm = TRUE)
min_val <- min(scaled_cpm_gender, na.rm = TRUE)
abs_max <- max(abs(c(max_val, min_val)))
breaks <- seq(-abs_max, abs_max, length.out = 101)
# Plot heatmap
pheatmap(scaled_cpm_gender_ordered,  
         annotation_col = annotation_col_gender_ordered,
         annotation_colors = annotation_colors_gender,
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
group_age <- factor(metadata$Age_Group, levels = c("Young", "Old"))
dge_age <- DGEList(counts = raw_count_table, group = group_age)
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
  Age_Group = metadata[colnames(cpm_age_sig), "Age_Group"]
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
         fontsize_row = 8,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species (Age Group)")
#differential analysis for smoking
group_smoking <- factor(metadata$Smoking_Status, levels = c("Ex_smoker", "Current_smoker"))
dge_smoking <- DGEList(counts = raw_count_table, group = group_smoking)
dge_smoking <- calcNormFactors(dge_smoking)
design_smoking <- model.matrix(~ 0 + group_smoking)
dge_smoking <- estimateDisp(dge_smoking, design_smoking)
fit_smoking <- glmFit(dge_smoking, design_smoking)
contrast_smoking <- makeContrasts(ExSmoker_vs_CurrentSmoker = group_smokingEx_smoker - group_smokingCurrent_smoker, levels = design_smoking)
lrt_smoking <- glmLRT(fit_smoking, contrast = contrast_smoking)
res_smoking <- topTags(lrt_smoking, n = Inf)$table
res_smoking$FDR <- p.adjust(res_smoking$PValue, method = "fdr")
res_smoking$Species <- rownames(res_smoking)
res_smoking_sig <- res_smoking %>%
  filter(FDR < 0.05)
top50_smoking_sig <- res_smoking_sig %>%
  arrange(FDR) %>%
  slice(1:50)
#bar chart
top50_smoking_sig$Higher_in <- ifelse(top50_smoking_sig$logFC > 0, "Ex_smoker", "Current_smoker")
ggplot(top50_smoking_sig, aes(x = reorder(Species, logFC), y = logFC, fill = Higher_in)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Ex_smoker" = "red", "Current_smoker" = "blue"),
    name = "Higher in"
  ) +
  labs(title = "Top Differentially Abundant Species (Ex smoker vs Current smoker)",
       x = "Species",
       y = "log2 Fold Change (Ex smoker - Current smoker)") +
  theme_minimal()
#volcano plot
res_smoking <- res_smoking %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up",
      FDR < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )
ggplot(res_smoking, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.7) +  
  geom_text_repel(
    data = res_smoking %>% arrange(FDR) %>% slice(1:10),
    aes(label = Species),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: Ex smoker vs Current smoker",
    x = "log2 Fold Change (Ex smoker - Current smoker)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()
#heatmap
cpm_smoking <- cpm(dge_smoking, log = FALSE)
sig_species_smoking <- res_smoking$Species[res_smoking$FDR < 0.05]
cpm_smoking_sig <- cpm_smoking[sig_species_smoking, , drop = FALSE]
scaled_cpm_smoking <- t(scale(t(cpm_age_smoking)))
annotation_col_smoking <- data.frame(
  Smoking_Status = metadata[colnames(cpm_smoking_sig), "Smoking_Status"]
)
rownames(annotation_col_smoking) <- colnames(cpm_smoking_sig)
annotation_colors_smoking <- list(
  Smoking_Status = c("Current_smoker" = "blue", "Ex_smoker" = "red")
)
sample_ordered_smoking <- rownames(annotation_col_smoking)[order(annotation_col_smoking$Smoking_Status)]
scaled_cpm_smoking_ordered <- scaled_cpm_smoking[, sample_ordered_smoking]
annotation_col_smoking_ordered <- annotation_col_smoking[sample_ordered_smoking, , drop = FALSE]
max_val <- max(scaled_cpm_smoking, na.rm = TRUE)
min_val <- min(scaled_cpm_smoking, na.rm = TRUE)
abs_max <- max(abs(c(max_val, min_val)))
breaks <- seq(-abs_max, abs_max, length.out = 101)
custom_colors <- colorRampPalette(c("darkblue", "lightblue", "white", "red", "darkred"))(100)
pheatmap(scaled_cpm_smoking_ordered,  
         annotation_col = annotation_col_smoking_ordered,
         annotation_colors = annotation_colors_smoking,
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         clustering_method = "complete",
         scale = "none",
         show_colnames = FALSE,
         fontsize_row = 8,
         color = custom_colors,
         breaks = breaks,
         main = "Heatmap of Significant Species (Smoking Status)")


#calculate diversity indices
install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages("ggplot2")
library(phyloseq)
library(ggplot2)
library(vegan)
library(ape)
install.packages("ggpubr")
library(ggpubr)
# OTU table (read counts)
otu <- otu_table(as.matrix(read_count_table_1), taxa_are_rows = TRUE)
#create phyloseq object
sample_data <- sample_data(metadata)
physeq <- phyloseq(otu, sample_data)
physeq
#beta diversity
#calculate Bray-Curtis dissimilarity
bray_dist <- phyloseq::distance(physeq, method = "bray")
print(bray_dist)
#perfrom PCoA
ord <- ordinate(physeq, method = "PCoA", distance = bray_dist)
print(ord)
# Get percent variance explained by PC1 and PC2
var_exp <- round(ord$values$Relative_eig[1:2] * 100, 1)
# Create data frame for plotting
plot_data <- plot_ordination(physeq, ord, type = "samples", justDF = TRUE)
# Add grouping variable 
plot_data$Group <- sample_data(physeq)$Group
# Plot with ggplot2
permanova_disease <- adonis2(bray_dist ~ Group, data = sample_df)
pval_disease <- signif(permanova_disease$`Pr(>F)`[1], 3)
ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +  # 95% confidence ellipses
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) (disease status)"
  ) +
  annotate("text", 
           x = 0.5,  # Change based on your plot scale
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_disease), 
           hjust = 1, size = 3.5) + theme_pub
group_levels <- unique(sample_df$Group)
# Pairwise comparisons
pairwise_results <- list()
for(i in 1:(length(group_levels)-1)) {
  for(j in (i+1):length(group_levels)) {
    subset_df <- sample_df[sample_df$Group %in% c(group_levels[i], group_levels[j]), ]
    subset_dist <- as.dist(as.matrix(bray_dist)[rownames(subset_df), rownames(subset_df)])
    result <- adonis2(subset_dist ~ Group, data = subset_df)
    pairwise_results[[paste(group_levels[i], "vs", group_levels[j])]] <- result
  }
}
pairwise_results
# smoking status
# Add grouping variable 
plot_data$Smoking_Status <- sample_data(physeq)$Smoking_Status
# Plot with ggplot2
permanova_smoking <- adonis2(bray_dist ~ Smoking_Status, data = sample_df)
pval_smoking <- signif(permanova_smoking$`Pr(>F)`[1], 3)
ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = Smoking_Status)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +  # 95% confidence ellipses
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) (smoking status)"
  )  +
  annotate("text", 
           x = 0.5,  # Change based on your plot scale
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_smoking), 
           hjust = 1, size = 3.5) + theme_pub

# age group
# Add grouping variable 
plot_data$Age_Group <- sample_data(physeq)$Age_Group
# Plot with ggplot2
permanova_ageGroup <- adonis2(bray_dist ~ Age_Group, data = sample_df)
pval_ageGroup <- signif(permanova_ageGroup$`Pr(>F)`[1], 3)
ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = Age_Group)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +  # 95% confidence ellipses
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) (Age Group)"
  ) +
  annotate("text", 
           x = 0.5,  # Change based on your plot scale
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_ageGroup), 
           hjust = 1, size = 3.5) + theme_pub

# sex
# Add grouping variable (change `Group` to your actual grouping variable)
plot_data$SEX <- sample_data(physeq)$SEX
# Plot with ggplot2
permanova_sex <- adonis2(bray_dist ~ SEX, data = sample_df)
pval_gender <- signif(permanova_sex$`Pr(>F)`[1], 3)
ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = SEX)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, type = "t") +  # 95% confidence ellipses
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "Beta Diversity (Bray-Curtis) (Gender)"
  ) +
  annotate("text", 
           x = 0.5,  # Change based on your plot scale
           y = 0.5, 
           label = paste0("PERMANOVA p = ", pval_gender), 
           hjust = 1, size = 3.5) + theme_pub

#alpha diversity
raw_count_table <- read.delim("20250412_Sherlock2_brushes_bact_readcounts.txt", header = TRUE, stringsAsFactors = FALSE)
raw_count_table <- raw_count_table[, -1] #do this twice 
rownames(raw_count_table) <- raw_count_table$Name
colnames(raw_count_table) <- gsub("^X", "", colnames(raw_count_table))
raw_count <- otu_table(as.matrix(raw_count_table), taxa_are_rows = TRUE)

sample_names(raw_count)
sample_names(sample_data)
sample_names(raw_count) <- gsub("\\.", "-", sample_names(raw_count))
setdiff(sample_names(raw_count),sample_names(sample_data))

physeq_1 <- phyloseq(raw_count, sample_data)
physeq_1
#calculate alpha diversity indices 
#observed
observed <- estimate_richness(physeq_1, measures = c("Observed"))
rownames(observed) <- gsub("^X", "", rownames(observed))
rownames(observed) <- gsub("\\.", "-", rownames(observed))
merged_observed <- cbind(metadata[rownames(observed),], observed)
merged_observed$Group[merged_observed$Group == "MildModerateCOPD"] <- "Moderate COPD"
merged_observed$Group[merged_observed$Group == "SevereCOPD"] <- "Severe COPD"

#gender
theme_pub <- theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )
ggboxplot(
  merged_observed,
  x = "SEX",
  y = "Observed",
  fill = "SEX",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Male", "Female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
# smoking status
ggboxplot(
  merged_observed,
  x = "Smoking_Status",
  y = "Observed",
  fill = "Smoking_Status",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Ex_smoker", "Current_smoker")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#age group
ggboxplot(
  merged_observed,
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
#disease 
ggboxplot(
  merged_observed,
  x = "Group",
  y = "Observed",
  fill = "Group",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_observed$Observed, na.rm = TRUE) + 100  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_observed$Group), 2, simplify = FALSE),
    label = "p.signif",       
    hide.ns = FALSE            
  ) +
  theme_minimal() +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 7)  # Adjust as needed
  )
#Shannon
shannon <- estimate_richness(physeq_1, measures = c("Shannon"))
merged_alpha <- cbind(merged_observed[rownames(shannon),], shannon)
#gender
ggboxplot(
  merged_alpha,
  x = "SEX",
  y = "Shannon",
  fill = "SEX",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Male", "Female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
# smoking status
ggboxplot(
  merged_alpha,
  x = "Smoking_Status",
  y = "Shannon",
  fill = "Smoking_Status",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Ex_smoker", "Current_smoker")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
#age group
ggboxplot(
  merged_alpha,
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
#disease 
ggboxplot(
  merged_alpha,
  x = "Group",
  y = "Shannon",
  fill = "Group",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha$Shannon, na.rm = TRUE) + 4  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha$Group), 2, simplify = FALSE),
    label = "p.signif",      
    hide.ns = FALSE            
  ) +
  theme_minimal() +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 8)  # Adjust as needed
  )

#chao1
chao1 <- estimate_richness(physeq_1, measures = c("Chao1"))
rownames(chao1) <- gsub("^X", "", rownames(chao1))
rownames(chao1) <- gsub("\\.", "-", rownames(chao1))
merged_alpha <- cbind(merged_alpha[rownames(chao1),], chao1)
#gender
ggboxplot(
  merged_alpha,
  x = "SEX",
  y = "Chao1",
  fill = "SEX",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Male", "Female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
# smoking status
ggboxplot(
  merged_alpha,
  x = "Smoking_Status",
  y = "Chao1",
  fill = "Smoking_Status",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Ex_smoker", "Current_smoker")),  # specify the pair to compare
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
#disease 
ggboxplot(
  merged_alpha,
  x = "Group",
  y = "Chao1",
  fill = "Group",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha$Chao1, na.rm = TRUE) + 100  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha$Group), 2, simplify = FALSE),
    label = "p.signif",       
    hide.ns = FALSE            
  ) +
  theme_minimal() +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 8)  # Adjust as needed
  )
#simpson
simpson <- estimate_richness(physeq_1, measures = c("Simpson"))
merged_alpha <- cbind(merged_alpha[rownames(simpson),], simpson)
#gender
ggboxplot(
  merged_alpha,
  x = "SEX",
  y = "Simpson",
  fill = "SEX",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Male", "Female")),  # specify the pair to compare
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() + theme_pub
# smoking status
ggboxplot(
  merged_alpha,
  x = "Smoking_Status",
  y = "Simpson",
  fill = "Smoking_Status",
  palette = "jco"
) +
  stat_compare_means(
    comparisons = list(c("Ex_smoker", "Current_smoker")),  # specify the pair to compare
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
#disease 
ggboxplot(
  merged_alpha,
  x = "Group",
  y = "Simpson",
  fill = "Group",
  palette = "jco"
) +
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(merged_alpha$Simpson, na.rm = TRUE) + 0.5  # Kruskal-Wallis p-value at the top
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = combn(unique(merged_alpha$Group), 2, simplify = FALSE),
    label = "p.signif",     
    hide.ns = FALSE            
  ) +
  theme_minimal() +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 8)  
  )

#taxSEA
BiocManager::install("TaxSEA")
library(TaxSEA)
data("TaxSEA_test_data")  # Loads example data
str(TaxSEA_test_data)
metabolite <- taxsea_result$Metabolite_producers

#for smoking
taxsea_data_smoking <- setNames(res_smoking_sig$logFC, res_smoking_sig$Species)
print(taxsea_data_smoking)
taxsea_result_smoking <- TaxSEA(taxon_ranks = taxsea_data_smoking)
metabolites_smoking <- taxsea_result_smoking$Metabolite_producers #FDR>0.05
health_smoking <- taxsea_result_smoking$Health_associations #FDR>0.05
bugsig_smoking <- taxsea_result_smoking$BugSigDB #FDR>0.05

#for gender 
taxsea_data_gender <- setNames(res_gender_sig$logFC, res_gender_sig$Species)
taxsea_result_gender <- TaxSEA(taxon_ranks = taxsea_data_gender) 
metabolites_gender <- taxsea_result_gender$Metabolite_producers #FDR>0.05
health_gender <- taxsea_result_gender$Health_associations #FDR>0.05
bugsig_gender <- taxsea_result_gender$BugSigDB #FDR>0.05

#for age group 
taxsea_data_ageGroup <- setNames(res_age_sig$logFC, res_age_sig$Species)
taxsea_result_ageGroup <- TaxSEA(taxon_ranks = taxsea_data_ageGroup)
metabolites_ageGroup <- taxsea_result_ageGroup$Metabolite_producers #FDR>0.05

#disease control vs moderate COPD
taxsea_data_diseaseMvsC <- setNames(res_mod_ctrl_sig$logFC, res_mod_ctrl_sig$Species)
taxsea_result_diseaseMvsC <- TaxSEA(taxon_ranks = taxsea_data_diseaseMvsC) #FDR>0.05
metabolites_diseaseMvsC <- taxsea_result_diseaseMvsC$Metabolite_producers

#disease control vs severe COPD
taxsea_data_diseaseSvsC <- setNames(res_sev_ctrl_sig$logFC, res_sev_ctrl_sig$Species)
taxsea_result_diseaseSvsC <- TaxSEA(taxon_ranks = taxsea_data_diseaseSvsC)
metabolites_diseaseSvsC <- taxsea_result_diseaseSvsC$Metabolite_producers #FDR>0.05

#disease moderate vs severe COPD
taxsea_data_diseaseSvsM <- setNames(res_sev_mod_sig$logFC, res_sev_mod_sig$Species)
taxsea_result_diseaseSvsM <- TaxSEA(taxon_ranks = taxsea_data_diseaseSvsM)
metabolites_diseaseSvsM <- taxsea_result_diseaseSvsM$Metabolite_producers #FDR>0.05
health_diseaseSvsM <- taxsea_result_diseaseSvsM$Health_associations #FDR>0.05
bugsig_diseaseSvsM <- taxsea_result_diseaseSvsM$BugSigDB #FDR>0.05

