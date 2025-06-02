################################################################################
#
# Figure 1: Prevalences of K13 622I mutation and RDT negativity in P. falciparum 
# infections in Northwest Ethiopia
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
# Figure 1 - Prevalences  of K13 622I mutation and RDT negativity in P. falciparum infections in Northwest Ethiopia.
################################################################################

# Load necessary libraries
library(ggplot2)
library(sf)
library(stringr)

# 1A) Nested maps with zoomed in collection map to show kebele within districts

# 1A.1 whole country

# Read the GeoJSON file for districts
geojson_file <- './data/mapping/Ethiopia_AdminBoundaries_woredas.geojson'
district_geojson <- st_read(geojson_file)

# Ensure CRS is consistent
district_geojson <- st_transform(district_geojson, crs = st_crs(4326))

# Filter for the districts of interest
districts_of_interest <- c("Tach Armacho", "Gonder Zuria")
filtered_geojson <- district_geojson %>%
  filter(WOREDANAME %in% districts_of_interest)

# Add a column for categorical district assignment
filtered_geojson <- filtered_geojson %>%
  mutate(district_category = case_when(
    WOREDANAME == "Tach Armacho" ~ "Tach Armachiho",
    WOREDANAME == "Gonder Zuria" ~ "Gondar Zuria"
  ))

bbox <- st_bbox(filtered_geojson)

# Plot the map
p <- ggplot() +
  geom_sf(data = district_geojson, fill = "grey90", color = "grey60", size = 0.2) +  # Background map
  geom_sf(data = filtered_geojson, aes(fill = district_category), color = "black", size = 0.3) +  # Highlight districts
  scale_fill_manual(values = custom_palette, name = "Districts") +  # Custom color palette
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.position = "none",
    plot.margin = margin(5.5, 40, 5.5, 5.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) 

filename <- "1A_whole_country_ethiopian_woredas_collection_sites"
save_plot_as_images(p,filename,figure_output_directory)

# 1A.2 Zoomed in collection sites map for inlay in full country map to show kebele

# Define a mapping of district names
district_name_mapping <- c(
  "Tach Armacho" = "Tach Armachiho",
  "Gonder Zuria" = "Gondar Zuria"
)

# Replace district names in the GeoJSON data with the corrected names
filtered_geojson <- filtered_geojson %>%
  mutate(WOREDANAME = recode(WOREDANAME, !!!district_name_mapping))

# Adjust the metadata to filter based on corrected district names and include kebele
unique_locations <- met %>%
  filter(district %in% district_name_mapping) %>%  # Filter for districts in metadata
  distinct(longitude, latitude, district, kebele)

# Remove Kurdi since it is in Lay Armachiho
unique_locations_filtered <- unique_locations %>%
  filter(kebele != "Kurbi")  

# Plot the map with sample locations and district names
p <- ggplot() +
  geom_sf(data = district_geojson, fill = "grey90", color = "grey60", size = 0.2) +  # Background map
  geom_sf(data = filtered_geojson, aes(fill = district_category), color = "black", size = 0.3) +  # Highlight districts
  scale_fill_manual(values = custom_palette, name = "Districts") +  # Custom color palette
  # Add red dots for sample locations
  geom_point(
    data = unique_locations_filtered,
    aes(x = longitude, y = latitude),
    color = "red",
    size = 2,
    shape = 21,
    fill = "red"
  ) +
  geom_sf_label(
    data = filtered_geojson,
    aes(label = WOREDANAME),
    size = 3,
    color = "black",
    fill = "white",
    label.padding = unit(0.1, "lines"),
    label.size = 0.2
  ) +
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"])
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.position = "none",
    plot.margin = margin(5.5, 40, 5.5, 5.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) 

filename <- "1A_ZOOM_ethiopian_woredas_collection_sites_zoom_with_woredas_labels"
save_plot_as_images(p,filename,figure_output_directory)

# Plot the map with sample locations WITHOUT district names
p <- ggplot() +
  geom_sf(data = district_geojson, fill = "grey90", color = "grey60", size = 0.2) +  # Background map
  geom_sf(data = filtered_geojson, aes(fill = district_category), color = "black", size = 0.3) +  # Highlight districts
  scale_fill_manual(values = custom_palette, name = "Districts") +  # Custom color palette
  # Add red dots for sample locations
  geom_point(
    data = unique_locations_filtered,
    aes(x = longitude, y = latitude),
    color = "red",
    size = 2,
    shape = 21,
    fill = "red"
  ) +
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"])
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.position = "none",
    plot.margin = margin(5.5, 40, 5.5, 5.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) 

filename <- "1A_ZOOM_ethiopian_woredas_collection_sites_zoom_without_woredas_labels"
save_plot_as_images(p,filename,figure_output_directory)

#1B) Grouped barplot by district of 622I prev and RDT negativity prev
# --- 1. Filter and prepare data ---
k13_samples <- variant_table_filtered_with_metadata_and_den %>%
  filter(gene_name == "k13", mutation == "Arg622Ile", !is.na(umi_alt), !is.na(district)) %>%
  mutate(
    district = factor(district),
    value = (umi_alt >= 2) * 100,  # Convert to percentage
    Metric = "622I_prevalence"
  ) %>%
  select(district, value, Metric)

rdt_samples <- variant_table_filtered_with_metadata_and_den %>%
  filter(gene_name == "k13", mutation == "Arg622Ile", !is.na(rdt_result), !is.na(district)) %>%
  mutate(
    district = factor(district),
    value = (rdt_result == "RDT negative") * 100,  # Convert to percentage
    Metric = "RDT_NegativityPercent"
  ) %>%
  select(district, value, Metric)

all_samples <- bind_rows(k13_samples, rdt_samples)

# --- 2. Run t-tests (Print to Table) ---
ttest_results <- all_samples %>%
  group_by(Metric) %>%
  summarise(
    p = t.test(value ~ district, var.equal = TRUE)$p.value,
    .groups = "drop"
  )

# Print t-test results as a table
print(ttest_results)

# --- 3. Summary stats for plotting (percentage values) ---
bar_data <- all_samples %>%
  group_by(district, Metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    lower = mean(value, na.rm = TRUE) - 1.96 * sd(value, na.rm = TRUE) / sqrt(n()),
    upper = mean(value, na.rm = TRUE) + 1.96 * sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# --- 4. Base plot with the correct color mapping ---
p <- ggplot(bar_data, aes(x = district, y = mean, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  scale_fill_manual(values = custom_palette) +
  labs(
    x = "District",
    y = "Prevalence (%)",  # Updated to reflect percentage
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "top"
  )

# --- 5. Extend y-axis to accommodate percentages ---
p <- p + coord_cartesian(ylim = c(0, 100), clip = "off")

# --- 6. Save plot ---
filename <- "1B_Grouped_Barplot_ByDistrict_622I_RDT"
#ggsave(paste0(filename, ".png"), plot = p, width = 8, height = 6, dpi = 300)

save_plot_as_images(p,filename,figure_output_directory)

#1C) Stacked barplot of 622I prevalence co-existing with and without RDT negativity

# prepare the data for plotting
stacked_data <- variant_table_filtered_with_metadata_and_den %>%
  filter(gene_name == "k13", mutation == "Arg622Ile", !is.na(kebele), !is.na(rdt_result), !is.na(district)) %>%
  mutate(
    carrier = ifelse(umi_alt >= 2, "Mutant", "Wildtype"),
    rdt_status = ifelse(rdt_result == "RDT negative", "RDT Negative", "RDT Positive"),
    group = paste(carrier, rdt_status, sep = " - ")
  ) %>%
  group_by(district, kebele, group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(district, kebele) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    group = str_trim(group),
    kebele = factor(kebele, levels = sort(unique(kebele)))  # Alphabetically sort kebeles
  )

# alphabetically sort the kebeles
stacked_data <- stacked_data %>%
  mutate(kebele = factor(kebele, levels = sort(unique(kebele))))

# define colors based on unified color pallete defined above
group_palette <- c(
  "Mutant - RDT Positive" = unname(custom_palette["Mutant"]),                # Red
  "Mutant - RDT Negative" = unname(custom_palette["622I_mut_and_RDT_neg"]),  # Purple
  "Wildtype - RDT Negative" = unname(custom_palette["RDT Negative"]),        # Blue
  "Wildtype - RDT Positive" = unname(custom_palette["Wildtype"])             # Grey
)

# Convert 'group' to factor with matching levels for color mapping
stacked_data <- stacked_data %>%
  mutate(group = factor(group, levels = names(group_palette)))

# plot the stacked barplot
p <- ggplot(stacked_data, aes(x = kebele, y = prop, fill = group)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ district, scales = "free_x") +  # Free x-axis to avoid repeating kebeles
  scale_fill_manual(values = group_palette, name = "622I & RDT Status") +
  labs(
    x = "Kebele",
    y = "Proportion of Samples",
    title = "Stacked Barplot of 622I & RDT Status by Kebele"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

filename <- "1C_stacked_barplot_622I_RDT_by_kebele"
save_plot_as_images(p,filename,figure_output_directory)
print(paste("Stacked barplot saved as", filename, ".png"))

# 1D) Boxplot of prevalence of pfhrp2-based RDT negativity by K13 622I status by kebele across districts.

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# Create a temporary color palette for this plot only
temp_palette <- c(
  "Mutant" = unname(custom_palette["622I_mut_and_RDT_neg"]),  # Purple for mutant
  "Wildtype" = unname(custom_palette["RDT Negative"])             # Grey for wildtype
)

# Filter for k13 mutation 622I only and relevant columns
rdt_data <- variant_table_filtered_with_metadata_and_den %>%
  filter(gene_name == "k13", mutation == "Arg622Ile", !is.na(kebele), !is.na(district)) %>%
  mutate(
    k13_status = case_when(
      umi_alt >= 2 ~ "Mutant",    # Define Mutant if umi_alt >= 2
      TRUE ~ "Wildtype"           # Otherwise, it's Wildtype
    ),
    rdt_negative = ifelse(rdt_result == "RDT negative", 1, 0)  # Convert RDT result to binary
  ) %>%
  group_by(district, kebele, k13_status) %>%
  summarise(
    rdt_neg_percentage = mean(rdt_negative, na.rm = TRUE) * 100,  # Convert to percentage
    sample_size = n(),
    .groups = "drop"
  )

# Ensure k13_status is a factor for proper grouping
rdt_data$k13_status <- factor(rdt_data$k13_status, levels = c("Mutant", "Wildtype"))

# Define the comparison for statistical testing (Mutant vs Wildtype)
comparisons <- list(c("Mutant", "Wildtype"))

# Generate the boxplot with facets by district
p <- ggplot(rdt_data, aes(x = k13_status, y = rdt_neg_percentage, fill = k13_status)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.5, color = "black") +  # Jittered points
  stat_compare_means(
    comparisons = comparisons, 
    method = "t.test", 
    label = "p.signif", 
    size = 6,                 # Increase p-value text size
    bracket.size = 1.2        # Increase thickness of the significance bracket
  ) +  
  facet_wrap(~ district) +  # Facet by district
  labs(
    x = "K13 622I Status",
    y = "RDT Negativity Percentage"
  ) +
  theme_minimal() +
  scale_fill_manual(values = temp_palette) +  # Apply the temporary custom palette
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label size
    legend.position = "none"  # Remove legend
  )

# Save the plot
filename <- "1D_boxplot_RDT_negativity_by_mutation_status_faceted"
save_plot_as_images(p, filename, figure_output_directory)