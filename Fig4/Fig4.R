################################################################################
#
# Figure 4: K13 mutant and wildtype parasites relatedness
# 
# Author: George Tollefson
# Date: June 2, 2025
#
################################################################################

################################################################################
#
# General Notes for all figure scripts: 
#
# Necessary raw data input files for Figures 1,4, and 5 are stored in the github 
# repo under /data and include:
# variant genotyping data
# sample metadata
# complexity of infection (COI) estimations
# IBD estimations 
# To run the scripts for Figures 1,4, and 5 you can change the loading paths in
# the Figure scripts 1,4 and 5 to match your local installation of the 
# repository with the raw data.

# Scripts for Figures 2 and 3 are packaged with the processed data which may be 
# loaded directly into the figure scripts provided.
#
################################################################################


################################################################################
#
# Data processing, custom function definitions. SHARED ACROSS FIG 1, 4, and 5
# Fig 1 specific code starts line 205
#
# Note: This must be run for Fig 1,4, and 5 and is included at head of each 
# figure script in the repo.
#
################################################################################

# Custom color definitions to be used throughout paper main figures
custom_palette <- c(
  "Gondar Zuria" = "#31A354",        # Green
  "Tach Armachiho" = "#ffc425",     # Yellow
  "RDT Negative" = "blue",             # Blue
  "RDT Positive" = "grey60",           # Grey
  "Mutant" = "red",             # Red
  "Wildtype" = "grey60",           # Grey
  "RDT_NegativityPercent" = "blue", # Blue
  "622I_prevalence" = "red",      # Red
  "622I_mut_and_RDT_neg" = "purple" # Purple
)

# Define output directory
figure_output_directory <- "./figures"

# Define unified figure save function
save_plot_as_images <- function(ggplot_object, filename_prefix, output_dir) {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define file paths
  svg_file <- file.path(output_dir, paste0(filename_prefix, ".svg"))
  pdf_file <- file.path(output_dir, paste0(filename_prefix, ".pdf"))
  png_file <- file.path(output_dir, paste0(filename_prefix, ".png"))
  
  # Save in SVG format
  ggsave(svg_file, ggplot_object, dpi = 600, width = 8, height = 6, units = "in", device = "svg")
  
  # Save in PDF format
  ggsave(pdf_file, ggplot_object, dpi = 600, width = 8, height = 6, units = "in", device = "pdf")
  
  # Save in PNG format
  ggsave(png_file, ggplot_object, dpi = 600, width = 8, height = 6, units = "in", device = "png")
  
  message("Plots saved as SVG, PDF, and PNG in the directory: ", output_dir)
}

# Define unified figure save function
save_plot_as_images_custom_ibd_network_two_legend_cols <- function(ggplot_object, filename_prefix, output_dir) {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define file paths
  svg_file <- file.path(output_dir, paste0(filename_prefix, ".svg"))
  pdf_file <- file.path(output_dir, paste0(filename_prefix, ".pdf"))
  png_file <- file.path(output_dir, paste0(filename_prefix, ".png"))
  
  # Save in SVG format
  ggsave(svg_file, ggplot_object, dpi = 600, width = 9, height = 6, units = "in", device = "svg")
  
  # Save in PDF format
  ggsave(pdf_file, ggplot_object, dpi = 600, width = 9, height = 6, units = "in", device = "pdf")
  
  # Save in PNG format
  ggsave(png_file, ggplot_object, dpi = 600, width = 9, height = 6, units = "in", device = "png")
  
  message("Plots saved as SVG, PDF, and PNG in the directory: ", output_dir)
}

################################################################################
#
# Prepare DR input data and metadata
#
################################################################################

# Load libraries
library(miplicorn)
library(dplyr)
library(readr)
library(rstatix)
library(ggplot2)
library(ggsignif)

# Define input and output directories

# Define data directory with input files
data_dir <- "./data"

# Load miptools tables
ref_file <- ("reference_AA_table.csv")
alt_file <- ("alternate_AA_table.csv")
cov_file <- ("coverage_AA_table.csv")

# put your table paths here:
ref_file <-read_tbl_reference(file.path(data_dir,ref_file))
alt_file <- read_tbl_alternate(file.path(data_dir,alt_file))
cov_file <- read_tbl_coverage(file.path(data_dir,cov_file))

# assemble table with ref, alt, and coverage data combined into one table
variant_table <-data.frame(sample_id=ref_file$sample,gene_id=ref_file$gene_id,gene_name=ref_file$gene,
                           mutation_name=ref_file$mutation_name,
                           mutation=ref_file$aa_change,umi_ref=ref_file$ref_umi_count,
                           umi_alt=alt_file$alt_umi_count, umi_cov=cov_file$coverage)

# make new column with total coverage for each locus and sample by adding UMI ref + UMI alt
variant_table<-mutate(variant_table,total_cov=umi_ref+umi_alt)

# Filtering based on UMI (read) coverage equal to or greater than 5
variant_table_filtered <- filter(variant_table, total_cov >= 5)

# Optionally write the filtered variant table kf1 to comma separated file to read in future instead of regenerating
# write.table(x=kf1,file = "filtered_variants_coverage_greater_than_5.csv",row.names=FALSE,col.names=TRUE)
# kf1 <- read.table("filtered_variants_coverage_greater_than_5.csv", header=TRUE, sep=",")

# Load metadata
met <- read.table("./data/MetadataAJZ_k13_status.tsv", header = T, sep = ",")

# For easier viewing during programming, subset metadata to columns we will use including: district, kebele, gametocyte, RDT result, season, occupation, lat and longitude
metadata <- met[,c("sample_id","month","season","seasons","year","district","kebele","longitude","latitude","gametocyte","age_group","sex","rdt_result", "treatment_history", "travel_history","parasitemia", "parasitemia", "age_group" )]

# Clean your kf1 table sample names to match metadata by removing the "-Gondar-1" from the kf1 table
variant_table_filtered$sample_id <- gsub("-Gondar-1", "", variant_table_filtered$sample_id)

# Merge metadata with variant table
colnames(metadata)[1] <- "sample_id"

variant_table_filtered_with_metadata <- merge(x=metadata, y=variant_table_filtered, by = "sample_id", all.x = T)

# Prepare prevalence table input for histograms and boxplots
variant_summary_stats <- variant_table_filtered_with_metadata %>%
  group_by(mutation_name) %>%
  get_summary_stats(total_cov,type = "full")

# extract denominator from kf4
den <- data.frame(mutation_name=variant_summary_stats$mutation_name,denominator=variant_summary_stats$n)

# add the denominator (total denom per variant in full sample set) to the variant table (this denominator is recalculated by metadata grouping variable in each figure)
variant_table_filtered_with_metadata_and_den <- merge(x=variant_table_filtered_with_metadata,y=den,by = "mutation_name", type=full)

################################################################################
# Load and preprocess COI data 
################################################################################

# load COI tables and add a new column indicating if the samepl_id has COI of 1 or greater than 1 to the existing dataframe of sample id key with 622I mutant status 
coi_data_pre <- read.table("./data/COI_output/IBC2FULL_PF_genomics_Ethiopia_DR2_IBC_merged.COI_dist_plot.samplemiss-25.sitemiss-20_COIfromGDS.tsv_summary.txt",sep="\t",header=T)
# Filter rows where Corp is "C"
filtered_coi_data_pre <- subset(coi_data_pre, CorP == "C")

# Create a dataframe with sample ids and COI values
coi_data <- data.frame(sample_ids = filtered_coi_data_pre$name, coi_values = filtered_coi_data_pre$mean)

# clean coi_samplenames
coi_data$sample_ids <- gsub("-Gondar-1", "", coi_data$sample_ids)

# Merge COI data with K13 mutant status
coi_data_status <- coi_data %>%
  left_join(variant_table_filtered_with_metadata %>%
              select(sample_id, gene_name, mutation, umi_alt) %>%
              filter(gene_name == "k13") %>%
              mutate(k13_status = ifelse(umi_alt > 1, "Mutant", "Wildtype")),
            by = c("sample_ids" = "sample_id")) %>%
  mutate(k13_status = ifelse(is.na(k13_status), "Wildtype", k13_status))

coi_data_status <- coi_data_status %>%
  filter(mutation == "Arg622Ile")

################################################################################
# Figure 4 - K13 mutant and wildtype parasites relatedness
################################################################################

# 4A) Pairwise IBD sharing between within monoclonal K13 R622I mutant and wildtype parasites within kebeles compared by study district

# Load required libraries
library(tidyverse)
library(ggpubr)

################################################################################
# Step 1: Load IBD Results 
################################################################################

# Load IBD results
ibd_results_path <- "./data/IBD_results/IBC2FULL_PF_genomics_Ethiopia.site20.sample_25.hmmidb.hmm_fract.txt"
ibd <- read.table(ibd_results_path, sep = "\t", header = TRUE)

################################################################################
# Step 2: Process IBD Data and Merge with Metadata
################################################################################

# Load COI data and filter for COI = 1
coi_one_samples <- coi_data %>%
  filter(coi_values == 1) %>%
  pull(sample_ids)

# Clean IBD sample names and filter for COI = 1
ibd1_prefilter <- ibd %>%
  select(sample1, sample2, fract_sites_IBD, N_informative_sites) %>%
  mutate(
    sample1 = gsub("-Gondar-1", "", sample1),
    sample2 = gsub("-Gondar-1", "", sample2)
  )

ibd1 <- ibd1_prefilter %>%
  filter(sample1 %in% coi_one_samples & sample2 %in% coi_one_samples) %>%
  filter(fract_sites_IBD >= 0.95 & N_informative_sites >= 50)

# Extract unique sample-to-district and kebele mapping
district_kebele_mapping <- variant_table_filtered_with_metadata %>%
  select(sample_id, district, kebele) %>%
  distinct()

# Extract unique K13 status per sample
k13_status_lookup <- variant_table_filtered_with_metadata %>%
  filter(gene_name == "k13") %>%
  mutate(K13_status = ifelse(mutation == "Arg622Ile" & umi_alt > 1, "Mutant", "Wildtype")) %>%
  select(sample_id, K13_status) %>%
  distinct()

# Merge district and kebele information for `sample1`
ibd1 <- ibd1 %>%
  left_join(district_kebele_mapping, by = c("sample1" = "sample_id")) %>%
  rename(district_sample1 = district, kebele_sample1 = kebele)

# Merge district and kebele information for `sample2`
ibd1 <- ibd1 %>%
  left_join(district_kebele_mapping, by = c("sample2" = "sample_id")) %>%
  rename(district_sample2 = district, kebele_sample2 = kebele)

# Merge K13 status for `sample1`
ibd1 <- ibd1 %>%
  left_join(k13_status_lookup, by = c("sample1" = "sample_id")) %>%
  rename(k13_status_sample1 = K13_status)

# Merge K13 status for `sample2`
ibd1 <- ibd1 %>%
  left_join(k13_status_lookup, by = c("sample2" = "sample_id")) %>%
  rename(k13_status_sample2 = K13_status)

# Filter for pairs where both samples have the same K13 status
filtered_ibd_coi_one_with_k13 <- ibd1 %>%
  filter(k13_status_sample1 == k13_status_sample2) %>%
  mutate(K13_status = k13_status_sample1)  # Create a single K13 status column

# Convert fract_sites_IBD to percentage
filtered_ibd_coi_one_with_k13 <- filtered_ibd_coi_one_with_k13 %>%
  mutate(fract_sites_IBD_percent = fract_sites_IBD * 100)

################################################################################
# Step 3: Aggregate Mean IBD Sharing at the Kebele Level
################################################################################

ibd_kebele_summary <- filtered_ibd_coi_one_with_k13 %>%
  group_by(district_sample1, kebele_sample1, K13_status) %>%
  summarise(
    mean_ibd_kebele = mean(fract_sites_IBD_percent, na.rm = TRUE),
    sd_ibd_kebele = sd(fract_sites_IBD_percent, na.rm = TRUE),
    .groups = "drop"
  )

# perform statistical tests (t-tests between within-kebele IBD means between districts)
pvalues_kebele <- ibd_kebele_summary %>%
  group_by(district_sample1) %>%
  summarise(
    p_value = t.test(mean_ibd_kebele ~ K13_status, var.equal = TRUE)$p.value
  ) %>%
  mutate(
    p_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Compute bracket positions for significance annotation at the district level
bracket_positions <- ibd_kebele_summary %>%
  group_by(district_sample1, K13_status) %>%
  summarise(ymax = max(mean_ibd_kebele, na.rm = TRUE), .groups = "drop") %>%
  spread(K13_status, ymax) %>%
  mutate(
    y_position = max(Mutant, Wildtype, na.rm = TRUE) + 5  # Adjust bracket above max value
  ) %>%
  left_join(pvalues_kebele, by = "district_sample1")

# Plot boxplot with kebele dots and p-values
ibd_boxplot <- ggplot(ibd_kebele_summary, aes(x = K13_status, y = mean_ibd_kebele)) +
  geom_boxplot(aes(fill = K13_status), outlier.shape = NA) +  # Boxplots by K13 status
  geom_jitter(aes(color = K13_status), width = 0.2, size = 2.5, shape = 21, fill = "white", stroke = 0.7) +  # Kebele dots
  scale_fill_manual(values = custom_palette[c("Mutant", "Wildtype")]) +
  scale_color_manual(values = custom_palette[c("Mutant", "Wildtype")]) +
  
  # Adding brackets for p-values (district-level significance)
  geom_segment(data = bracket_positions, aes(x = 1, xend = 2, y = y_position, yend = y_position), color = "black", linewidth = 0.5) +
  geom_text(data = bracket_positions, aes(x = 1.5, y = y_position + 2, label = p_label), size = 8, fontface = "italic") +
  
  # Formatting
  theme_minimal() +
  labs(
    title = "IBD Sharing by K13 Status and District",
    x = "K13 Status",
    y = "Mean IBD Sharing per Kebele (%)"
  ) +
  theme(
    text = element_text(size = 12),  # Ensures all text is at least size 12
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),  # X-axis title size 14
    axis.title.y = element_text(size = 14),  # Y-axis title size 14
    plot.title = element_blank(),  # Main title size 14
    strip.text = element_text(size = 14, face = "bold"),  # Facet-wrapped (district) title size 14
    legend.title = element_blank(),  # Removes legend title
    legend.text = element_text(size = 12),  # Removes legend title
    legend.position = "none"  # Keeps the legend in place
  )  +
  facet_wrap(~ district_sample1) +  # Facet by district
  coord_cartesian(ylim = c(0, max(ibd_kebele_summary$mean_ibd_kebele, na.rm = TRUE) + 10))

# Print the plot
print(ibd_boxplot)

filename <- "4A_mean_IBD_sharing_COI1_622I mutant_vs_wildtype_by_district"
save_plot_as_images(ibd_boxplot, filename,figure_output_directory)

#4B) Network analysis of highly related parasite pairs  (IBD≥0.95) by RDT test results shaped by K13 mutation status

# Load required libraries
library(ggraph)
library(igraph)
library(tidygraph)
library(dplyr)
library(ggplot2)

# Updated Plot Network Function
plot_network_two_cats_4B <- function(met, coi_data, ibd, figure_output_directory, 
                                  color_var, shape_var, color_title, shape_title, custom_palette) {
  
  # Make a temporary copy of the met data to avoid altering the original
  met_temp <- met
  
  # Adjust RDT status format only within this function (without altering the original data)
  met_temp <- met_temp %>%
    mutate(Result.RDT = case_when(
      Result.RDT == "RDT positive" ~ "RDT Positive",
      Result.RDT == "RDT negative" ~ "RDT Negative",
      TRUE ~ Result.RDT
    ))
  
  # Extract COI = 1 samples
  coi_one_samples <- coi_data %>%
    filter(coi_values == 1) %>%
    pull(sample_ids)
  
  # Clean IBD sample names and filter for COI = 1
  ibd1 <- ibd %>%
    select(sample1, sample2, fract_sites_IBD, N_informative_sites) %>%
    mutate(
      sample1 = gsub("-Gondar-1", "", sample1),
      sample2 = gsub("-Gondar-1", "", sample2)
    )
  
  ibd_filtered <- ibd1 %>%
    filter(sample1 %in% coi_one_samples & sample2 %in% coi_one_samples) %>%
    filter(fract_sites_IBD >= 0.95 & N_informative_sites >= 50)
  
  # Prepare metadata for plotting
  metadata_clean <- met_temp %>%
    mutate(
      K13_status = case_when(
        k13_allele == "Arg622Ile" ~ "622I",
        k13_allele == "Cys580Tyr" ~ "580Y",
        TRUE ~ "Wildtype"
      )
    ) %>%
    select(sample_id, Result.RDT, K13_status) %>%
    distinct(sample_id, .keep_all = TRUE)
  
  # Ensure edges have matching samples in metadata
  edges <- ibd_filtered %>%
    select(sample1, sample2, fract_sites_IBD) %>%
    rename(from = sample1, to = sample2, weight = fract_sites_IBD) %>%
    filter(from %in% metadata_clean$sample_id & to %in% metadata_clean$sample_id)
  
  # Create nodes with valid samples only
  nodes <- metadata_clean %>%
    filter(sample_id %in% unique(c(edges$from, edges$to))) %>%
    rename(name = sample_id)
  
  # Build the network graph
  g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  tidy_g <- as_tbl_graph(g)
  
  # Define shape palette for K13 mutations
  shape_palette <- c("622I" = 17, "580Y" = 15, "Wildtype" = 16)  # Triangle, Square, Circle
  
  # Generate the network plot
  network_plot <- ggraph(tidy_g, layout = "fr") +
    geom_edge_link(alpha = 0.5, color = "gray70") +
    geom_node_point(aes(color = Result.RDT, shape = K13_status), size = 4) +
    scale_color_manual(values = custom_palette) +
    scale_shape_manual(values = shape_palette) +
    labs(
      color = color_title,
      shape = shape_title
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    )
  
  return(network_plot)
}

Fig4B_network_plot <- plot_network_two_cats_4B(
  met = met,
  coi_data = coi_data,
  ibd = ibd,
  figure_output_directory = figure_output_directory,
  color_var = "Result.RDT",
  shape_var = "K13_status",
  color_title = "RDT Status",
  shape_title = "K13 Mutation",
  custom_palette = custom_palette
)

filename <- "4B_RDT_and_k13_status_network"
save_plot_as_images(Fig4B_network_plot,filename,figure_output_directory)

# 4C) Network analysis of highly related parasite pairs  (IBD≥0.95) by district shaped by K13 mutation status.

# Plot Network Function for District and K13 Status
plot_network_two_cats_4C <- function(met, coi_data, ibd, figure_output_directory, 
                                  color_var, shape_var, color_title, shape_title, custom_palette) {
  
  # Make a temporary copy of the met data to avoid altering the original
  met_temp <- met
  
  # Adjust RDT status format only within this function (without altering the original data)
  met_temp <- met_temp %>%
    mutate(Result.RDT = case_when(
      Result.RDT == "RDT positive" ~ "RDT Positive",
      Result.RDT == "RDT negative" ~ "RDT Negative",
      TRUE ~ Result.RDT
    ))
  
  # Extract COI = 1 samples
  coi_one_samples <- coi_data %>%
    filter(coi_values == 1) %>%
    pull(sample_ids)
  
  # Clean IBD sample names and filter for COI = 1
  ibd1 <- ibd %>%
    select(sample1, sample2, fract_sites_IBD, N_informative_sites) %>%
    mutate(
      sample1 = gsub("-Gondar-1", "", sample1),
      sample2 = gsub("-Gondar-1", "", sample2)
    )
  
  ibd_filtered <- ibd1 %>%
    filter(sample1 %in% coi_one_samples & sample2 %in% coi_one_samples) %>%
    filter(fract_sites_IBD >= 0.95 & N_informative_sites >= 50)
  
  # Prepare metadata for plotting
  metadata_clean <- met_temp %>%
    mutate(
      K13_status = case_when(
        k13_allele == "Arg622Ile" ~ "622I",
        k13_allele == "Cys580Tyr" ~ "580Y",
        TRUE ~ "Wildtype"
      )
    ) %>%
    select(sample_id, District, K13_status) %>%
    distinct(sample_id, .keep_all = TRUE)
  
  # Ensure edges have matching samples in metadata
  edges <- ibd_filtered %>%
    select(sample1, sample2, fract_sites_IBD) %>%
    rename(from = sample1, to = sample2, weight = fract_sites_IBD) %>%
    filter(from %in% metadata_clean$sample_id & to %in% metadata_clean$sample_id)
  
  # Create nodes with valid samples only
  nodes <- metadata_clean %>%
    filter(sample_id %in% unique(c(edges$from, edges$to))) %>%
    rename(name = sample_id)
  
  # Build the network graph
  g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  tidy_g <- as_tbl_graph(g)
  
  # Define shape palette for K13 mutations
  shape_palette <- c("622I" = 17, "580Y" = 15, "Wildtype" = 16)  # Triangle, Square, Circle
  
  # Generate the network plot
  network_plot <- ggraph(tidy_g, layout = "fr") +
    geom_edge_link(alpha = 0.5, color = "gray70") +
    geom_node_point(aes(color = District, shape = K13_status), size = 4) +
    scale_color_manual(values = custom_palette) +
    scale_shape_manual(values = shape_palette) +
    labs(
      color = color_title,
      shape = shape_title
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    )
  
  return(network_plot)
}

# Plot 4C - District and k13 status
Fig4C_network_plot <- plot_network_two_cats_4C(
  met = met,
  coi_data = coi_data,
  ibd = ibd,
  figure_output_directory = figure_output_directory,
  color_var = "District",
  shape_var = "K13_status",
  color_title = "District",
  shape_title = "K13 Mutation",
  custom_palette = custom_palette[c("Gondar Zuria", "Tach Armachiho")]
)

# Save the plot
filename <- "4C_District_and_k13_status_network"
save_plot_as_images(Fig4C_network_plot, filename, figure_output_directory)

print(Fig4C_network_plot)

#4D) Network analysis of highly related parasites (IBD≥0.95) at Kebele (village) level shaped by K13 mutation status

# Kebele custom colors
kebele_colors <- c(
  "dodgerblue2", "#E31A1C",  # red
  "green4", "#6A3D9A",       # purple
  "#FF7F00",                 # orange
  "black", "gold1",
  "skyblue2", "#FB9A99",     # lt pink
  "palegreen2", "#CAB2D6",   # lt purple
  "#FDBF6F",                 # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "black"
)
kebele_levels <- unique(met$Kebele)
names(kebele_colors) <- kebele_levels[seq_along(kebele_colors)]


plot_network_two_cats_4D <- function(met, coi_data, ibd, figure_output_directory, 
                                  color_var, shape_var, color_title, shape_title, custom_palette) {
  
  ################################################################################
  # Filter IBD Data for COI = 1 and High IBD Sharing
  ################################################################################
  
  # Extract COI = 1 samples
  coi_one_samples <- coi_data %>%
    filter(coi_values == 1) %>%
    pull(sample_ids)
  
  # Process IBD file (remove "-Gondar-1" from sample IDs)
  ibd1 <- ibd %>%
    select(sample1, sample2, fract_sites_IBD, N_informative_sites) %>%
    mutate(sample1 = gsub("-Gondar-1", "", sample1),
           sample2 = gsub("-Gondar-1", "", sample2))
  
  # Filter for COI = 1 samples and high IBD sharing
  ibd_filtered <- ibd1 %>%
    filter(sample1 %in% coi_one_samples & sample2 %in% coi_one_samples) %>%
    filter(fract_sites_IBD >= 0.95 & N_informative_sites >= 50)
  
  # Create edges from IBD data
  edges <- ibd_filtered %>%
    select(sample1, sample2, fract_sites_IBD) %>%
    rename(from = sample1, to = sample2, weight = fract_sites_IBD)
  
  ################################################################################
  # Prepare Metadata for Network
  ################################################################################
  # First, create K13_status column to ensure it exists
  met <- met %>%
    mutate(
      K13_status = case_when(
        k13_allele == "Arg622Ile" ~ "622I",
        k13_allele == "Cys580Tyr" ~ "580Y",
        TRUE ~ "Wildtype"
      )
    )
  
  # Select relevant metadata based on user input
  metadata_clean <- met %>%
    select(sample_id, all_of(color_var), all_of(shape_var), K13_status) %>%
    distinct(sample_id, .keep_all = TRUE)  # Remove duplicates
  
  
  # Create nodes from metadata (only include nodes present in edges)
  nodes <- metadata_clean %>%
    filter(sample_id %in% unique(c(edges$from, edges$to))) %>%
    rename(name = sample_id)
  
  # Build network graph
  g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  
  # Convert to tidygraph format (all metadata is already in g, no need to merge again)
  tidy_g <- as_tbl_graph(g)
  
  ################################################################################
  # Define Aesthetics (Color and Shape)
  ################################################################################
  
  # Define shape palette for K13 mutations
  shape_palette <- c("622I" = 17, "580Y" = 15, "Wildtype" = 16)  # Triangle, Square, Circle
  
  ################################################################################
  # Generate Network Plot
  ################################################################################
  
  # Generate network plot with dynamic metadata categories and custom legend titles
  network_plot <- ggraph(tidy_g, layout = "fr") +
    geom_edge_link(alpha = 0.5, color = "gray70") +
    geom_node_point(aes_string(color = color_var, shape = shape_var), size = 4) +
    scale_color_manual(values = custom_palette) +
    scale_shape_manual(values = shape_palette) +
    labs(
      color = color_title,   # Custom title for color legend
      shape = shape_title    # Custom title for shape legend
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    )
  
  
  return(network_plot)
}

# Plot 4D - Kebele and k13 status
Fig4D_network_plot <- plot_network_two_cats_4D(
  met = met,
  coi_data = coi_data,
  ibd = ibd,
  figure_output_directory = figure_output_directory,
  color_var = "Kebele",
  shape_var = "K13_status",
  color_title = "Kebele",
  shape_title = "K13 Mutation",
  custom_palette = kebele_colors
)

filename <- "4D_Kebele_k13_status_network"
save_plot_as_images_custom_ibd_network_two_legend_cols(Fig4D_network_plot,filename,figure_output_directory)
