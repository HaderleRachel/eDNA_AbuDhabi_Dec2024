####Carto####
library(sf)
library(ggplot2)
library(terra)       # for get_tiles()
library(scales)      # for rgb()

# Load your transects from the KML (drop the 3rd dimension if there is one)
transects <- st_read("sites.kml", options = "Z=drop")
transects
# Define the area around Abu Dhabi (WGS84)
xmin <- 54.25
xmax <- 54.5
ymin <- 24.35
ymax <- 24.5
bbox <- st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), crs = 4326)

# Download the satellite basemap
library(maptiles)
basemap <- get_tiles(bbox, provider = "Esri.WorldImagery", crop = TRUE, zoom = 15)

# Convert to data.frame for ggplot
basemap_df <- as.data.frame(basemap, xy = TRUE)
basemap_df$fill <- rgb(
  basemap_df$red,
  basemap_df$green,
  basemap_df$blue,
  maxColorValue = 255
)
basemap_df$fill <- rgb(
  basemap_df$Esri.WorldImagery_15_21321_14082_1,
  basemap_df$Esri.WorldImagery_15_21321_14082_2,
  basemap_df$Esri.WorldImagery_15_21321_14082_3,
  maxColorValue = 255
)

# Plot
ggplot() +
  geom_raster(data = basemap_df, aes(x = x, y = y, fill = fill)) +
  scale_fill_identity() +
  geom_sf(data = transects, color = "red", size = 1) +
  coord_sf(
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    expand = FALSE
  ) +
  theme_minimal()



transects_centroids <- st_centroid(transects)
transects_centroids$id <- seq_len(nrow(transects_centroids))

# Add the names manually
transects_centroids$name <- c("Mangrove", "Al Hudayriat", "Seagrass", "Bateen")

# Map with custom names
ggplot() +
  geom_raster(data = basemap_df, aes(x = x, y = y, fill = fill)) +
  scale_fill_identity() +
  geom_sf(data = transects, color = "red", size = 1) +
  geom_sf_text(
    data = transects_centroids,
    aes(label = name),
    color = "white",
    size = 4,
    fontface = "bold"
  ) +
  coord_sf(
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    expand = FALSE
  ) +
  labs(x = "Longitude", y = "Latitude")+
  theme_minimal()


library(sf)
library(dplyr)

# Compute the distance matrix (in meters if geometry is WGS84)
dist_matrix <- st_distance(transects_centroids)

# Convert to long data.frame
dist_df <- as.data.frame(as.table(as.matrix(dist_matrix)))
colnames(dist_df) <- c("site1_index", "site2_index", "distance")

# Add site names
dist_df$site1 <- transects_centroids$name[dist_df$site1_index]
dist_df$site2 <- transects_centroids$name[dist_df$site2_index]

# Keep only each unique pair (site1 < site2)
dist_df <- dist_df %>%
  filter(site1 != site2) %>%          # remove 0 distances
  rowwise() %>%
  mutate(tmp = list(sort(c(site1, site2)))) %>%
  mutate(site1 = tmp[[1]], site2 = tmp[[2]]) %>%
  dplyr::select(-tmp) %>%
  distinct()

dist_df


####eDNA Data####
data_abu <- read.csv2(file="input_data.csv")


# Relevant columns
spy_cols <- c("SPY2405886", "SPY2405893", 
              "SPY2405895", "SPY2405896",
              "SPY2405887", "SPY2405892",
              "SPY2405890", "SPY2405891")

# Robust cleaning function
clean_numeric <- function(x) {
  x <- gsub("[\u00A0\u202F ]", "", x)  # remove normal, non-breaking, and narrow non-breaking spaces
  x[x == ""] <- NA                    # replace empty strings with NA
  as.numeric(x)
}

# Apply to all SPY* columns
data_abu[spy_cols] <- lapply(data_abu[spy_cols], clean_numeric)


# Step 2: Create grouped columns
data_abu$Al_Hudayriat <- rowSums(data_abu[, c("SPY2405886", "SPY2405893")], na.rm = TRUE)
data_abu$Bateen <- rowSums(data_abu[, c("SPY2405895", "SPY2405896")], na.rm = TRUE)
data_abu$Mangrove <- rowSums(data_abu[, c("SPY2405887", "SPY2405892")], na.rm = TRUE)
data_abu$Seagrass <- rowSums(data_abu[, c("SPY2405890", "SPY2405891")], na.rm = TRUE)

# Step 3: Convert to presence/absence
data_abu$Al_Hudayriat <- ifelse(data_abu$Al_Hudayriat > 0, 1, 0)
data_abu$Bateen <- ifelse(data_abu$Bateen > 0, 1, 0)

data_abu$Mangrove <- ifelse(data_abu$Mangrove > 0, 1, 0)
data_abu$Seagrass <- ifelse(data_abu$Seagrass > 0, 1, 0)

data_abu <- data_abu[, !(names(data_abu) %in% spy_cols)]

colSums(data_abu[, (ncol(data_abu)-3):ncol(data_abu)] == 1)

data_abu
write.xlsx(data_abu, "tableau taxons motus.xlsx")

###Bargraph & Venn diagram###
library(dplyr)
library(tidyr)
library(ggplot2)

library(ggvenn)
library(patchwork)

# ----Barplot by class by site ----
data_abu <- read.csv2("tableau taxons motus.csv")

library(dplyr)

n_unique_species <- data_abu %>%
  filter(!is.na(specificEpithet) & specificEpithet != "") %>%  # keep non-null species
  distinct(TAXON, .keep_all = TRUE) %>%       # keep one unique line per TAXON
  nrow()                                       # count number of rows

n_unique_species


taxon_counts <- data_abu %>%
  dplyr::select(TAXON, Al_Hudayriat, Bateen, Mangrove, Seagrass) %>%
  pivot_longer(cols = c(Al_Hudayriat, Bateen, Mangrove, Seagrass), names_to = "site", values_to = "presence") %>%
  filter(presence == 1) %>%
  distinct(site, TAXON) %>%
  count(site, name = "n_taxa")

print(taxon_counts)


total_taxa <- data_abu %>%
  filter(Al_Hudayriat == 1 | Bateen == 1 | Mangrove == 1 | Seagrass == 1) %>%
  distinct(TAXON) %>%
  count()


print(total_taxa)

# Convert to long format
data_long <- data_abu %>%
  dplyr::select(TAXON, class, Al_Hudayriat, Bateen, Mangrove, Seagrass) %>%
  pivot_longer(cols = c(Al_Hudayriat, Bateen, Mangrove, Seagrass), names_to = "site", values_to = "presence") %>%
  filter(presence == 1) %>%
  distinct(site, class, TAXON)

# Count unique MOTUs per site and per class
motu_counts <- data_long %>%
  group_by(site, class) %>%
  summarise(n_motus = n(), .groups = "drop")

# Ensure all site x class combinations exist (even with 0)
all_combos <- expand.grid(
  site = unique(motu_counts$site),
  class = unique(motu_counts$class)
)

motu_counts <- left_join(all_combos, motu_counts, by = c("site", "class")) %>%
  mutate(n_motus = ifelse(is.na(n_motus), 0, n_motus))


# Barplot
# Relevel the factor to define the order of sites
motu_counts$site <- factor(motu_counts$site,
                           levels = c("Al_Hudayriat", "Bateen", "Mangrove", "Seagrass"))

# Barplot
barplot <- ggplot(motu_counts, aes(x = site, y = n_motus, fill = class)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.5) +
  geom_text(aes(label = n_motus),
            position = position_dodge(width = 0.6),
            vjust = -0.3, size = 6) +
  labs(
    x = "Site",
    y = "Number of taxa",
    fill = "Class",
    title = NULL
  ) +
  theme_minimal(base_size = 16)

library(ggplot2)

barplot <- ggplot(motu_counts, aes(x = site, y = n_motus, fill = class)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.7), 
           width = 0.6, color = "black") +
  geom_text(aes(label = n_motus),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 4.5, family = "sans") +
  labs(
    x = "Site",
    y = "Number of taxa",
    fill = "Class"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.x = element_blank(),   # remove vertical gridlines
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",                # move legend above plot
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 14),
    plot.margin = margin(10, 10, 10, 10)    # breathing space
  ) +
  scale_fill_brewer(palette = "Set2")       # nice colorblind-friendly palette


library(ggplot2)

barplot <- ggplot(motu_counts, aes(x = site, y = n_motus, fill = class)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.7), 
           width = 0.6, color = "black") +
  geom_text(aes(label = n_motus),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 4.5, family = "sans") +
  labs(
    x = "Site",
    y = "Number of taxa",
    fill = "Class"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.x = element_blank(),   # remove vertical gridlines
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",                # move legend above plot
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 14),
    plot.margin = margin(10, 10, 10, 10)    # breathing space
  ) +
  scale_fill_brewer(palette = "Set2")       # nice colorblind-friendly palette

barplot


library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
data_abu_fb <- read.xlsx("tableau taxons publi.xlsx")

# Site columns
site_cols <- names(data_abu_fb)[13:16]

# Data preparation
plot_data <- data_abu_fb %>%
  filter(class == "Teleostei", !is.na(DemersPelag)) %>%
  pivot_longer(
    cols = all_of(site_cols),
    names_to = "site",
    values_to = "present"
  ) %>%
  filter(present == 1) %>%
  group_by(site, DemersPelag) %>%
  summarise(count = n(), .groups = "drop")


# ---- 2. Venn diagram with ggvenn ----
library(dplyr)
library(ggvenn)
library(patchwork)
library(ggplot2)

# Filter data_abu to keep only Teleostei with non-NA Class
data_abu_filtered <- data_abu %>%
  filter(class == "Teleostei")

data_abu_filtered <- data_abu %>%
  filter(class == "Teleostei") %>%
  group_by(TAXON) %>%
  summarise(across(tail(names(.), 4), ~ sum(.), .names = "{col}")) %>%
  ungroup() %>%
  mutate(across(tail(names(.), 4), ~ ifelse(. > 1, 1, .)))

##############
# 1. Venn diagram by TAXA
venn_list_taxa <- list(
  Al_Hudayriat = data_abu_filtered$TAXON[data_abu_filtered$Al_Hudayriat == 1],
  Bateen = data_abu_filtered$TAXON[data_abu_filtered$Bateen == 1],
  Mangrove = data_abu_filtered$TAXON[data_abu_filtered$Mangrove == 1],
  Seagrass = data_abu_filtered$TAXON[data_abu_filtered$Seagrass == 1]
)

venn_plot_taxa <- ggvenn(venn_list_taxa, show_percentage = FALSE,
                         fill_color = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                         stroke_size = 0.3, set_name_size = 4, text_size = 5) + theme(plot.margin = margin(20, 20, 20, 20))+
  ggtitle("TAXA")

# Filter data_abu to keep only Teleostei with non-NA Class
data_abu_filtered <- data_abu %>%
  filter(!is.na(class), class == "Teleostei")

data_abu_filtered <- data_abu %>%
  filter(!is.na(class), class == "Teleostei")%>%
  group_by(TAXON, family) %>%
  summarise(across(tail(names(.), 4), ~ sum(.), .names = "{col}")) %>%
  ungroup()%>%
  mutate(across(tail(names(.), 4), ~ ifelse(. > 1, 1, .)))

############

# 2. Venn diagram by FAMILY
venn_list_family <- list(
  Al_Hudayriat = data_abu_filtered$family[data_abu_filtered$Al_Hudayriat == 1],
  Bateen = data_abu_filtered$family[data_abu_filtered$Bateen == 1],
  Mangrove = data_abu_filtered$family[data_abu_filtered$Mangrove == 1],
  Seagrass = data_abu_filtered$family[data_abu_filtered$Seagrass == 1]
)


venn_plot_family <- ggvenn(venn_list_family, show_percentage = FALSE,
                           fill_color = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                           stroke_size = 0.3, set_name_size = 4, text_size = 5) + theme(plot.margin = margin(20, 20, 20, 20))+
  
  ggtitle("FAMILIES")

# 3. Combination with centered main title
final_plot <-  venn_plot_family + venn_plot_taxa +
  plot_layout(ncol = 2, widths = c(1.1, 1)) +
  plot_annotation(
    title = "Comparison of Shared Families and Taxa Across Sites",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# Display
print(final_plot)


library(dplyr)
library(tibble)
library(purrr)
library(combinat)

get_venn_details <- function(venn_list) {
  sites <- names(venn_list)
  all_taxa <- unique(unlist(venn_list))
  
  # For each taxon/family, determine in how many (and which) sites it appears
  venn_df <- map_dfr(all_taxa, function(taxon) {
    present_in <- names(venn_list)[map_lgl(venn_list, ~ taxon %in% .x)]
    tibble(
      Element = taxon,
      Sites = paste(sort(present_in), collapse = ","),
      N_sites = length(present_in)
    )
  }) %>% arrange(desc(N_sites))
  
  return(venn_df)
}

# Details for TAXA
venn_details_taxa <- get_venn_details(venn_list_taxa)
cat("=== TAXA DETAILS ===\n")
print(venn_details_taxa)

# Details for FAMILIES
venn_details_family <- get_venn_details(venn_list_family)
cat("=== FAMILY DETAILS ===\n")


print(venn_details_family)

# Export details to CSV files
write.csv(venn_details_taxa, "venn_details_taxa.csv", row.names = FALSE)
write.csv(venn_details_family, "venn_details_family.csv", row.names = FALSE)

# Number of elements per exact site combination
summary_taxa <- venn_details_taxa %>%
  count(Sites, name = "N_taxa") %>%
  arrange(desc(N_taxa))

print(summary_taxa)

summary_family <- venn_details_family %>%
  count(Sites, name = "N_families") %>%
  arrange(desc(N_families))

print(summary_family)



####taxonomic beta-div####
library(dplyr)
library(tidyr)
library(adespatial)

# Step 1: keep only the necessary columns
data_subset <- data_abu %>%
  dplyr::select(TAXON, Al_Hudayriat, Bateen, Mangrove, Seagrass)

# Step 2: identify duplicates and rename them
#data_subset <- data_subset %>%
#  group_by(scientificName) %>%
#  mutate(scientificName = make.unique(scientificName, sep = "_")) %>%
#  ungroup()

# Step 3: pivot to long format (sites in rows)
data_long <- data_subset %>%
  pivot_longer(cols = Al_Hudayriat:Seagrass, names_to = "Site", values_to = "Presence")

# Step 4: pivot to wide format (one scientific name per column)
library(dplyr)
library(tidyr)

data_long_unique <- data_long %>%
  group_by(Site, TAXON) %>%
  summarise(Presence = max(Presence), .groups = "drop")

data_wide <- data_long_unique %>%
  pivot_wider(
    names_from = TAXON,
    values_from = Presence,
    values_fill = 0
  )


data_wide[is.na(data_wide)] <- 0

data_matrix_wide <- as.data.frame(data_wide)
rownames(data_matrix_wide) <- data_matrix_wide$Site
data_matrix_wide$Site <- NULL


jaccard<-data_matrix_wide
jaccard[,] <- lapply(jaccard[,], as.numeric)

jaccard[,] <- lapply(jaccard[,], as.numeric)
jaccard[is.na(jaccard)] <- 0

JAC <- beta.div.comp(jaccard, quant = FALSE, save.abc = FALSE)
JAC

JACD <- JAC$D
JACD

JAC_sim <- 1 - JACD
JACturn <- JAC$repl
JACnest <- JAC$rich

####phylogenetic beta-div####
# -------------------------------
# 1. Load packages
# -------------------------------
library(dplyr)
library(tidyr)
library(ape)
library(U.PhyloMaker)
library(ggtree)
library(phangorn)

# -------------------------------
# 2. Prepare the data
# -------------------------------
# Site columns
sites <- c("Al_Hudayriat", "Bateen", "Mangrove", "Seagrass")

# Long pivot and occurrence filtering
data_long <- data_abu %>%
  pivot_longer(
    cols = all_of(sites),
    names_to = "site",
    values_to = "occurrence"
  ) %>%
  filter(!is.na(occurrence) & occurrence != 0) %>%
  distinct(TAXON, site, .keep_all = TRUE)

# Presence summary by site
data_site <- data_long %>%
  group_by(TAXON, site) %>%
  summarise(
    Presence = as.numeric(any(occurrence > 0)),
    across(where(~ !is.list(.)), ~ first(.)),
    .groups = "drop"
  )

# Species list with taxonomy
taxa_list <- tibble::tibble(
  species = data_site$TAXON,
  genus = data_site$genus,
  family = data_site$family,
  order = data_site$order
)

sp.list <- taxa_list %>%
  dplyr::select(species, genus, family, order) %>%
  mutate(species.relative = NA, genus.relative = NA)

# -------------------------------
# 3. Load megatrees
# -------------------------------
megatree_fish <- read.tree("phylogeny/fish_megatree.tre")
megatree_mm <- read.tree("phylogeny/mammal_megatree.tre")
megatree_oi <- read.tree("phylogeny/bird_megatree.tre")
megatree_shark <- read.tree("phylogeny/shark_megatree_consensus.tre")

# Create a dummy root
combined_tree <- bind.tree(megatree_mm, megatree_oi, where = "root", position = 0)
combined_tree <- bind.tree(combined_tree, megatree_fish, where = "root", position = 0)
combined_tree <- bind.tree(combined_tree, megatree_shark, where = "root", position = 0)


plot(combined_tree)


# -------------------------------
# 4. Load genus lists
# -------------------------------
gen.list_fish <- read.csv("phylogeny/fish_genus_list.csv")
gen.list_mm   <- read.csv("phylogeny/mammal_genus_list.csv")
gen.list_oi   <- read.csv("phylogeny/bird_genus_list.csv")
gen.list_shark <- read.csv("phylogeny/Species.list.csv")

# Merge genus lists
gen.list <- rbind(gen.list_fish, gen.list_mm, gen.list_oi)

# -------------------------------
# 5. Generate the phylogeny with U.PhyloMaker
# -------------------------------
result <- U.PhyloMaker::phylo.maker(
  sp.list,
  combined_tree,
  gen.list,
  nodes.type = 1,
  scenario = 3
)

phylo <- result$phylo
phylo
plot(phylo)
phylo$edge.length[phylo$edge.length < 0] <- 1e-6
plot(phylo)

library(phangorn)

phylo_rooted <- midpoint(phylo)
is.rooted(phylo_rooted)

plot(phylo_rooted)
library(ggtree)

# -------------------------------
# 6. Visualize the phylogeny
# -------------------------------
# Rectangular
ggtree(phylo_rooted, layout = "rectangular") +
  geom_tiplab(size = 2)


colnames(jaccard) <- gsub(" ", "_", colnames(jaccard))

# Names in your dissimilarity matrix
colnames(jaccard)

# Tip names in your tree
phylo$tip.label

library(picante)
unifrac <- unifrac(jaccard, phylo_rooted)
unifrac

unifrac_sim <- 1 - unifrac

library(betapart)
phylobeta <- phylo.beta.pair(jaccard, phylo_rooted, index.family = "jaccard")

phylobeta

phylobetaturn <- phylobeta$phylo.beta.jtu
phylobetanest<- phylobeta$phylo.beta.jne

# Global sums for turnover, nestedness, and total dissimilarity
total_turnover <- sum(phylobetaturn)  # Total turnover contribution
total_nestedness <- sum(phylobetanest)  # Total nestedness contribution
total_dissimilarity <- total_turnover + total_nestedness  # Total dissimilarity

# Calculate percentages
percent_turnover <- (total_turnover / total_dissimilarity) * 100
percent_nestedness <- (total_nestedness / total_dissimilarity) * 100

# Display results
cat("Global percentage of turnover (phylo.beta.jtu):", round(percent_turnover, 2), "%\n")
cat("Global percentage of nestedness (phylo.beta.jne):", round(percent_nestedness, 2), "%\n")


####PCA beta-div####
# Function to calculate mean pairwise β-diversity for each station
mean_pairwise_beta <- function(mat) {
  mat <- as.matrix(mat)
  station_means <- rowMeans(mat, na.rm = TRUE)  # Mean of each row
  return(station_means)
}

# Apply the function to each β-diversity matrix
Bjne_Phyl <- mean_pairwise_beta(phylobetanest)

Bjtu_Phyl <- mean_pairwise_beta(phylobetaturn)


Bjtu_Tax <- mean_pairwise_beta(JACturn)


Bjne_Tax <- mean_pairwise_beta(JACnest)

# Combine the mean β-diversity matrices into a single matrix
combined_matrix <- cbind(Bjne_Phyl, Bjtu_Phyl, Bjtu_Tax, Bjne_Tax)
combined_matrix
rownames(combined_matrix) <- c(
  "Al_Hudayriat",
  "Bateen",
  "Mangrove",
  "Seagrass"
)


#### Turnover/nestedness matrix ####
library(pheatmap)

# Example: JAC$repl and JAC$rich are your triangular matrices

# 1. Convert to matrices
jac_turn <- as.matrix(JAC$repl)
jac_nest <- as.matrix(JAC$rich)

# 2. Complete each matrix to make them symmetric
complete_symmetric <- function(mat) {
  full <- mat
  full[upper.tri(full)] <- t(mat)[upper.tri(mat)]
  diag(full) <- 0  # Optional: or NA if preferred
  return(full)
}

jac_turn_full <- complete_symmetric(jac_turn)
jac_nest_full <- complete_symmetric(jac_nest)

# 3. Function to combine turnover and nestedness into a single matrix
combine_beta_matrices <- function(turn_mat, nest_mat) {
  stopifnot(all(dim(turn_mat) == dim(nest_mat)))
  stopifnot(all(rownames(turn_mat) == rownames(nest_mat)))
  n <- nrow(turn_mat)
  combined <- matrix(NA, n, n)
  rownames(combined) <- rownames(turn_mat)
  colnames(combined) <- colnames(turn_mat)
  
  combined[upper.tri(combined)] <- turn_mat[upper.tri(turn_mat)]
  combined[lower.tri(combined)] <- nest_mat[lower.tri(nest_mat)]
  diag(combined) <- NA  # or 0 if you want to display the diagonal
  
  return(combined)
}

# 4. Apply the function
beta_taxo_combined <- combine_beta_matrices(jac_turn_full, jac_nest_full)

# 5. Visualization (optional but publication-ready)
pheatmap(
  beta_taxo_combined,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Taxonomic β-diversity",
  fontsize = 12
)

phyl_turn_full <- complete_symmetric(as.matrix(phylobetaturn))
phyl_nest_full <- complete_symmetric(as.matrix(phylobetanest))

beta_phylo_combined <- combine_beta_matrices(phyl_turn_full, phyl_nest_full)

pheatmap(
  beta_phylo_combined,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Phylogenetic β-diversity",
  fontsize = 12
)

# Station names
statnam <- c("Al_Hudayriat", "Bateen", "Mangrove", "Seagrass")

res.pca <- prcomp(combined_matrix)
library(factoextra)
RES.PCA <- fviz_pca_biplot(res.pca, habillage = statnam, labelsize = 5, repel = TRUE, legend.title = "Sample", title = NULL) +
  theme(text = element_text(size = 15))
RES.PCA

library(ggplot2)
library(factoextra)

# Compute % of variance
var_pct <- round(100 * (res.pca$sdev^2) / sum(res.pca$sdev^2), 1)



RES.PCA <- fviz_pca_biplot(
  res.pca,
  habillage = statnam,        # Group colors
  label = "all",              # Show all labels (set to "ind" for individuals)
  repel = TRUE,               # Avoid label overlap
  addEllipses = TRUE,         # Add 95% ellipses
  ellipse.level = 0.95,
  pointsize = 2,              # Point size (if individuals are shown)
  labelsize = 4,  
  shape.ind = 16,
  legend.title = "Sample",    # Legend title
  title = NULL
) +
  scale_color_viridis_d(option = "plasma", begin = 0, end = 0.85) +  # Clean colors
  theme_minimal(base_size = 16) +         # Clean theme
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    axis.line = element_line(color = "black")
  )

RES.PCA


####Faith's PD and SES####
####Calculating PD####
colnames(jaccard) <- gsub(" ", "_", colnames(jaccard))
#rownames(jaccard_par_annee) <- jaccard_par_annee$annee

pd.result <- pd(jaccard, phylo_rooted, include.root = TRUE)

# Check the results
print(pd.result)


phydist <- cophenetic(phylo)


# Extract common species
species_common <- intersect(colnames(jaccard), phylo$tip.label)

# Filter jaccard_par_annee and phydist to keep only these species
jaccard_par_annee_filt <- jaccard[, species_common]
phydist_filt <- cophenetic(phylo)[species_common, species_common]



# Re-run the calculation
ses.mpd.result <- ses.mpd(
  jaccard_par_annee_filt,
  phydist_filt,
  null.model = "taxa.labels",
  abundance.weighted = FALSE,
  runs = 1000
)
ses.mpd.result

ses.mntd.result <- ses.mntd(
  jaccard_par_annee_filt,
  phydist_filt,
  null.model = "taxa.labels",
  abundance.weighted = FALSE,
  runs = 1000
)
ses.mntd.result

# Calculate the standardized effect size of phylogenetic diversity for fish
ses.pd.result <- ses.pd(jaccard, phylo_rooted, null.model = "taxa.labels", run = 1000)
ses.pd.result



library(lirrr)
lsf.str("package:lirrr")

mvpd.result <- mvpd(jaccard_par_annee_filt,
                    phydist_filt,
                    abundance.weighted = FALSE)

mvpd.result


ses.vpd <- function(samp, dis, null.model = c("taxa.labels", "richness", "frequency", "sample.pool", 
                                              "phylogeny.pool", "independentswap", "trialswap"), 
                    abundance.weighted = FALSE, runs = 999, iterations = 1000) {
  dis <- as.matrix(dis)
  
  # Calculate observed VPD values by extracting the `vpd` column from mvpd results
  vpd.obs <- mvpd(samp, dis, abundance.weighted = abundance.weighted)$vpd
  null.model <- match.arg(null.model)
  
  # Calculate VPD values for each null model
  vpd.rand <- switch(null.model,
                     taxa.labels = t(replicate(runs, mvpd(samp, taxaShuffle(dis), abundance.weighted=abundance.weighted)$vpd)),
                     richness = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="richness"), dis, abundance.weighted)$vpd)),
                     frequency = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="frequency"), dis, abundance.weighted)$vpd)),
                     sample.pool = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="richness"), dis, abundance.weighted)$vpd)),
                     phylogeny.pool = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="richness"), 
                                                             taxaShuffle(dis), abundance.weighted)$vpd)),
                     independentswap = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="independentswap", iterations), dis, abundance.weighted)$vpd)),
                     trialswap = t(replicate(runs, mvpd(randomizeMatrix(samp, null.model="trialswap", iterations), dis, abundance.weighted)$vpd))
  )
  
  # Calculate mean and standard deviation of simulated values
  vpd.rand.mean <- apply(X = vpd.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
  vpd.rand.sd <- apply(X = vpd.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
  
  # Calculate SES values for VPD
  vpd.obs.z <- (vpd.obs - vpd.rand.mean) / vpd.rand.sd
  vpd.obs.rank <- apply(X = rbind(vpd.obs, vpd.rand), MARGIN = 2, FUN = rank)[1, ]
  vpd.obs.rank <- ifelse(is.na(vpd.rand.mean), NA, vpd.obs.rank)
  
  data.frame(ntaxa = specnumber(samp), vpd.obs, vpd.rand.mean, vpd.rand.sd, vpd.obs.rank, 
             vpd.obs.z, vpd.obs.p = vpd.obs.rank / (runs + 1), runs = runs, row.names = row.names(samp))
}

ses.vpd.result <- ses.vpd(jaccard_par_annee_filt, phydist_filt, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)

ses.vpd.result



library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# 1. Extract PD SES
pd_df <- ses.pd.result %>%
  rownames_to_column(var = "site") %>%
  dplyr::select(site, value = pd.obs.z, p = pd.obs.p) %>%
  mutate(metric = "PD")

# 2. Extract MPD SES
mpd_df <- ses.mpd.result %>%
  rownames_to_column(var = "site") %>%
  dplyr::select(site, value = mpd.obs.z, p = mpd.obs.p) %>%
  mutate(metric = "MPD")

# 3. Extract VPD SES
vpd_df <- ses.vpd.result %>%
  rownames_to_column(var = "site") %>%
  dplyr::select(site, value = vpd.obs.z, p = vpd.obs.p) %>%
  mutate(metric = "VPD")

# 4. Merge
metrics_df <- bind_rows(pd_df, mpd_df, vpd_df) %>%
  mutate(sig = ifelse(p < 0.05 | p > 0.95, "*", ""))

metrics_df <- metrics_df %>%
  mutate(site = gsub("_", " ", site))

ggplot(metrics_df, aes(x = site, y = value, fill = metric)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.8), 
           color = "black", width = 0.7) +
  geom_text(aes(label = sig),
            position = position_dodge(width = 0.8),
            vjust = -0.8, size = 5, color = "black", fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.6) +
  scale_fill_manual(values = c("MPD" = "#E41A1C", 
                               "PD" = "#4DAF4A", 
                               "VPD" = "#377EB8")) +
  labs(x = "Site", y = "SES (obs.z)", fill = "Metric") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
