################################################################################
#
# Figure 5: Monthly prevalence of the K13 622I mutation, RDT negativity, P. falciparum positivity and genetic relatedness in P. falciparum parasites
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
# Figure 5 - Monthly prevalence of the K13 622I mutation, RDT negativity, 
# P. falciparum positivity and genetic relatedness in P. falciparum parasites
################################################################################

#5A) Stacked barplot of RDT and 622I with case cound and rainfall trendlines

# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(terra)
library(stringr)
library(viridis)
library(grDevices)

month_order <- c("November", "December", "January", "February", "March", "April", 
                 "May", "June", "July", "August", "September", "October")

# load the rainfall data generated from CHIRPS data as described in the paper (at kebele level, by year and month, averaged to district level by month in appropriate)
rainfall_for_plot <- read.table("data/misc/rainfall/rainfall_data_from_CHIRPS_for_districts_averaged_by_kebele.csv", 
                                sep = ",", 
                                header = TRUE, 
                                stringsAsFactors = FALSE)

# Prepare rainfall data for line overlay: 
# Note: This code was used to prepare the rainfall data file uploaded to the github repo for this project. Commenting out since can not upload the ti ffiles for rainfall to github due to size limits.
# This produces the file in ./data/misc/rainfall/rainfall_data_from_CHIRPS_for_districts_averaged_by_kebele.csv
#
# # define site coordinates
# sites <- data.frame(
#   Site = c("Gondar Zuria", "Tach Armachiho"),
#   Latitude = c(12.4, 13.0),
#   Longitude = c(37.5, 37.3)
# )
# 
# # Convert to SpatVector
# sites_vect <- vect(sites, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
# 
# # === 2. List all .tif files ===
# tif_dir <- "./data/misc/rainfall"
# tif_files <- list.files(path = tif_dir, pattern = "chirps-v3.0.*\\.tif$", full.names = TRUE)
# 
# # Extract rainfall for each site from each raster
# all_data <- lapply(tif_files, function(tif) {
#   # Load raster
#   r <- rast(tif)
#   
#   # Extract rainfall
#   values <- terra::extract(r, sites_vect)[, 2]
#   
#   # Parse year and month from filename
#   month_str <- str_extract(tif, "\\d{4}\\.\\d{2}")
#   
#   # Combine with site info
#   data.frame(
#     Site = sites$Site,
#     Latitude = sites$Latitude,
#     Longitude = sites$Longitude,
#     Month = month_str,
#     Rainfall_mm = values
#   )
# }) %>% bind_rows()
# 
# rainfall_data <- all_data %>%
#   mutate(
#     Month = as.numeric(as.character(Month)),  # Ensure numeric
#     year = floor(Month),
#     month_num = round((Month - year) * 100),
#     month_name = month.name[month_num],
#     month_name = factor(month_name, levels = month_order)
#   ) %>%
#   filter(month_name %in% month_order) %>%
#   rename(district = Site) %>%
#   select(district, month = month_name, Rainfall_mm) %>%
#   mutate(month = factor(month, levels = month_order))
# 
# # --- Prepare rainfall data for line overlay ---
# rainfall_for_plot <- rainfall_data %>%
#   group_by(district, month) %>%
#   summarise(Value = mean(Rainfall_mm, na.rm = TRUE), .groups = "drop") %>%
#   mutate(
#     Metric = "Rainfall (mm)"
#   )

group_palette <- c(
  "622I Mutation" = unname(custom_palette["Mutant"]),
  "RDT Negative" = unname(custom_palette["RDT Negative"])
)

variant_table_filtered_with_metadata_and_den <- variant_table_filtered_with_metadata_and_den %>%
  mutate(month = factor(month, levels = month_order))

# 622I prevalence (only among k13/622I sequenced samples)
mut_622I_data <- variant_table_filtered_with_metadata_and_den %>%
  filter(gene_name == "k13", mutation == "Arg622Ile", !is.na(umi_alt), !is.na(month), !is.na(district)) %>%
  group_by(district, month) %>%
  summarise(
    Proportion = mean(umi_alt >= 2, na.rm = TRUE),
    Metric = "622I Mutation",
    .groups = "drop"
  )

# RDT negativity prevalence
rdt_data <- variant_table_filtered_with_metadata_and_den %>%
  filter(!is.na(rdt_result), !is.na(month), !is.na(district)) %>%
  mutate(is_rdt_negative = rdt_result == "RDT negative") %>%
  group_by(district, month) %>%
  summarise(
    Proportion = mean(is_rdt_negative, na.rm = TRUE),
    Metric = "RDT Negative",
    .groups = "drop"
  )

# Combine into one dataset
prevalence_data <- bind_rows(mut_622I_data, rdt_data) %>%
  mutate(
    month = factor(month, levels = month_order),
    Metric = factor(Metric, levels = c("622I Mutation", "RDT Negative"))
  )

# Temporary custom color palette for this plot
temp_palette <- c(
  "622I Mutation" = unname(custom_palette["Mutant"]),      # Red
  "RDT Negative" = unname(custom_palette["RDT Negative"])  # Blue
)

# Define line colors for Rainfall and Case Count
line_colors <- c(
  "Rainfall (mm)" = "deepskyblue",
  "Case Count" = "#36454F"
)

# Line types for rainfall and case count
line_types <- c(
  "Rainfall (mm)" = "solid",
  "Case Count" = "solid"
)

# Prepare case count data for legend
# add monthly case counts
# load and clean case data
gondar_cases <- read_csv("./data/misc/case_counts/Montly_Pf_case_Counts_Gondar_Zuria.csv") %>%
  mutate(district = "Gondar Zuria")

tach_cases <- read_csv("./data/misc/case_counts/Monthly_Pf_case_Count_Tach_Armachiho.csv") %>%
  mutate(district = "Tach Armachiho")

# Combine and clean
case_counts <- bind_rows(gondar_cases, tach_cases) %>%
  rename(
    month = Month,
    cases = `Number of P. falciparum infected cases`
  ) %>%
  mutate(
    month = factor(month, levels = month_order)
  )

case_for_plot <- case_counts %>%
  mutate(
    Value = cases,
    Metric = "Case Count"
  )

# Combine for plotting
line_data <- bind_rows(rainfall_for_plot, case_for_plot)

# Find the max value for scaling
max_rainfall <- max(rainfall_for_plot$Value, na.rm = TRUE)
max_cases <- max(case_for_plot$Value, na.rm = TRUE)
max_value <- max(max_rainfall, max_cases)

# --- Plot: Grouped bars + rainfall and case count lines ---
p <- ggplot() +
  # Grouped bars for 622I and RDT (scaled to percentage)
  geom_bar(
    data = prevalence_data,
    aes(x = month, y = Proportion * 100, fill = Metric),  # Multiply by 100 for percentage
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  # Combined line for Rainfall and Case Count
  geom_line(
    data = line_data,
    aes(
      x = month,
      y = (Value / max_value) * 100,  # Scale the line values to match percentage
      color = Metric,
      linetype = Metric,
      group = Metric
    ),
    size = 1.2
  ) +
  # Point overlay for rainfall and case counts
  geom_point(
    data = line_data,
    aes(
      x = month,
      y = (Value / max_value) * 100,  # Scale the point values to match percentage
      color = Metric,
      shape = Metric
    ),
    size = 3,
    fill = "white"
  ) +
  facet_wrap(~ district) +
  scale_fill_manual(values = temp_palette) +
  scale_color_manual(values = line_colors) +
  scale_linetype_manual(values = line_types) +
  scale_shape_manual(values = c("Rainfall (mm)" = 21, "Case Count" = 22)) +
  scale_y_continuous(
    limits = c(0, 110),
    labels = scales::label_number(scale = 1),  # Display as whole numbers (percent)
    name = "Prevalence",
    sec.axis = sec_axis(
      trans = ~ . * (max_value / 100),  # Scale back to original line values
      name = "Line Metric Count"
    )
  ) +
  labs(
    x = "Month",
    y = "Prevalence (%)",   # Updated y-axis label
    fill = "Bar Metric",
    color = "Line Metric",
    linetype = "Line Metric",
    shape = "Line Metric"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y.right = element_text(color = "blue", size = 13),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "top"
  )

# --- Save the plot ---
filename <- "5A_grouped_barplot_with_case_counts_and_rainfall_overlay_percentage_fixed"
save_plot_as_images(p, filename, figure_output_directory)

#5B) Network diagram of 622I by month

plot_network_two_cats_month <- function(met, coi_data, ibd, figure_output_directory, 
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
  
  month_order <- c(
    "November", "December", "January", "February", "March", "April",
    "May", "June", "July", "August", "September", "October"
  )
  
  # Generate network plot with dynamic metadata categories and custom legend titles
  network_plot <- ggraph(tidy_g, layout = "fr") +
    geom_edge_link(alpha = 0.5, color = "gray70") +
    geom_node_point(aes_string(color = color_var, shape = shape_var), size = 4) +
    scale_color_manual(
      values = custom_palette,
      breaks = month_order,  # ensures legend follows Nov–Oct
      name = "Month"
    ) +
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

# Reorder months to start from November
met$Month <- factor(met$Month, levels = c(
  "November", "December", "January", "February", "March", "April", 
  "May", "June", "July", "August", "September", "October"
))

# Define color gradient palette following this new order
custom_palette_months <- setNames(
  grDevices::colorRampPalette(c("#8B0000", "#FF6347", "#FFC0CB", "#ADD8E6", "#00008B"))(12),
  levels(met$Month)
)

# Plot 5C - Month and k13 status
Fig5B_network_plot <- plot_network_two_cats_month(
  met = met,
  coi_data = coi_data,
  ibd = ibd,
  figure_output_directory = figure_output_directory,
  color_var = "Month",
  shape_var = "K13_status",
  color_title = "Month",
  shape_title = "K13 Mutation",
  custom_palette = custom_palette_months
)

filename <- "5B_month_and_k13_status_network"
save_plot_as_images(Fig5B_network_plot,filename,figure_output_directory)

#5C)
# Load necessary libraries
library(dplyr)
library(ggplot2)

# Set month order: November to October
month_order <- c("November", "December", "January", "February", "March", "April", 
                 "May", "June", "July", "August", "September", "October")

# Create lookup table for sample -> month
month_lookup <- met %>%
  select(sample_id, Month) %>%
  filter(!is.na(Month)) %>%
  mutate(Month = factor(Month, levels = month_order))

# Merge month info with both samples in IBD pairs
ibd_month <- filtered_ibd_coi_one_with_k13 %>%
  left_join(month_lookup, by = c("sample1" = "sample_id")) %>%
  rename(month1 = Month) %>%
  left_join(month_lookup, by = c("sample2" = "sample_id")) %>%
  rename(month2 = Month) %>%
  filter(month1 == month2) %>%
  mutate(Month = month1)

# Summarise mean IBD by month and K13 status
monthly_ibd_summary <- ibd_month %>%
  group_by(K13_status, Month) %>%
  summarise(
    mean_ibd = mean(fract_sites_IBD_percent, na.rm = TRUE),
    sd_ibd = sd(fract_sites_IBD_percent, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Month = factor(Month, levels = month_order)
  )

################################################################################
# Updated: Monthly Mean IBD Barplot Grouped by Month (No Faceting)
################################################################################

# Grouped barplot (Mutant/Wildtype side-by-side per month)
monthly_ibd_grouped_plot <- ggplot(monthly_ibd_summary, aes(x = Month, y = mean_ibd, fill = K13_status)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_errorbar(
    aes(ymin = mean_ibd - sd_ibd, ymax = mean_ibd + sd_ibd),
    position = position_dodge(width = 0.8), width = 0.2
  ) +
  scale_fill_manual(values = custom_palette[c("Mutant", "Wildtype")]) +
  labs(
    title = "Monthly Mean IBD Sharing by K13 Genotype",
    x = "Month",
    y = "Mean IBD Sharing (%)",
    fill = "K13 Status"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )

filename <- "5C_monthly_mean_IBD_sharing_grouped_barplot"
save_plot_as_images(monthly_ibd_grouped_plot, filename, figure_output_directory)

print("Grouped barplot of monthly mean IBD sharing (grouped by month) saved successfully!")

# 5D)

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(forcats)

# Ensure month is ordered
month_order <- c(
  "November", "December", "January", "February", "March", "April", 
  "May", "June", "July", "August", "September", "October"
)

# Define group color palette with desired order
group_palette <- c(
  "Mutant - RDT Positive" = unname(custom_palette["Mutant"]),          # Red (Top)
  "Mutant - RDT Negative" = unname(custom_palette["622I_mut_and_RDT_neg"]),  # Purple (Middle)
  "Wildtype - RDT Negative" = unname(custom_palette["RDT Negative"]),  # Blue (Second from Bottom)
  "Wildtype - RDT Positive" = unname(custom_palette["Wildtype"])       # Grey (Bottom)
)

# Prepare stacked data by month and district
stacked_month_data <- variant_table_filtered_with_metadata_and_den %>%
  filter(gene_name == "k13", mutation == "Arg622Ile", 
         !is.na(month), !is.na(district), !is.na(rdt_result)) %>%
  mutate(
    carrier = ifelse(umi_alt >= 2, "Mutant", "Wildtype"),
    rdt_status = ifelse(rdt_result == "RDT negative", "RDT Negative", "RDT Positive"),
    group = paste(carrier, rdt_status, sep = " - "),
    month = factor(month, levels = month_order)
  ) %>%
  group_by(district, month, group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(district, month) %>%
  mutate(prevalence = (n / sum(n)) * 100) %>%
  ungroup() %>%
  mutate(
    group = factor(str_trim(group), levels = names(group_palette))
  )

# Plot
p <- ggplot(stacked_month_data, aes(x = month, y = prevalence, fill = group)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  facet_wrap(~ district, scales = "free_x") +
  scale_fill_manual(values = group_palette, name = "622I & RDT Status") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Show as percentage
  labs(
    x = "Month",
    y = "Prevalence (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

# Save
filename <- "5D_stacked_barplot_622I_RDT_by_month_percent"
save_plot_as_images(p, filename, figure_output_directory)

print("Stacked barplot with prevalence (percentage) saved successfully!")


