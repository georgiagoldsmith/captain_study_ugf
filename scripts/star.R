library(rredlist)
library(readxl)
library(here)
library(tidyverse)
library(sf)
library(terra)
library(readr)

###################################################################
# SECTION 1: Assign IUCN status to UGF birds
# Queries the IUCN Red List API for each species in the bird list
# and retrieves their conservation status and population trend
###################################################################

# set IUCN_REDLIST_KEY in .Renviron rather than hardcoding it here
api_key <- Sys.getenv("IUCN_REDLIST_KEY")

# Load bird list — resume from checkpoint if available to avoid
# re-querying species already processed in a previous run
checkpoint_path <- "birds_progress.csv"

if (file.exists(checkpoint_path)) {
  birds <- read.csv(checkpoint_path)
  cat("Resuming from checkpoint:", checkpoint_path, "\n")
} else {
  birds <- read_excel(here("data", "birds_ugf.xlsx"))
  # split scientific name into genus and species for API queries
  birds <- birds |>
    separate(scientific_name, into = c("genus", "species"),
             sep = " ", remove = FALSE)
  birds$IUCN           <- NA
  birds$population_trend <- NA
}

# Query IUCN API for each species
# skips species already filled in (safe resume after interruption)
for (i in seq_len(nrow(birds))) {
  if (!is.na(birds$IUCN[i])) next
  cat("Species", i, "of", nrow(birds), "-", birds$scientific_name[i], "\n")
  
  out <- tryCatch(
    rl_species(genus = birds$genus[i], species = birds$species[i], key = api_key),
    error = function(e) NULL
  )
  
  if (!is.null(out) && nrow(out$assessments) > 0) {
    # filter to most recent global assessment only
    global_latest <- out$assessments[
      out$assessments$latest == TRUE & grepl("Global", out$assessments$scopes), 
    ]
    if (nrow(global_latest) > 0) {
      birds$IUCN[i] <- global_latest$red_list_category_code[1]
      # get population trend from assessment details
      assess_detail <- tryCatch(
        rl_assessment(id = global_latest$assessment_id[1], key = api_key),
        error = function(e) NULL
      )
      if (!is.null(assess_detail)) {
        birds$population_trend[i] <- assess_detail$population_trend$description$en
      }
    }
  }
  
  Sys.sleep(0.5)  # rate limiting 
  
  # save checkpoint every 25 species in case of interruption
  if (i %% 25 == 0) {
    write.csv(birds, checkpoint_path, row.names = FALSE)
    cat("  Progress saved at species", i, "\n")
  }
}

write.csv(birds, "birds_species_status.csv", row.names = FALSE)

###################################################################
# SECTION 2: Load all data
###################################################################

birds_ugf <- read_csv(here("data/Birds w species status/birds_species_status.csv"))
aoh_list <- read_csv(here("data/birds_aoh/Birds_list_AOH.csv"))
ugf_boundary <- st_read(here("data/UGF_gp.shp/UGF_gp.shp"))

# reproject boundary to WGS84 to match AOH rasters
ugf_boundary_wgs84 <- st_transform(ugf_boundary, 4326)
ugf_vect <- vect(ugf_boundary_wgs84)  # convert to terra SpatVector for raster operations

# path to clipped AOH rasters
clipped_dir <- here("data/birds_aoh/clipped")

###################################################################
# SECTION 3: Match bird species to AOH list
# The AOH list (Lumbierres et al.) uses BirdLife taxonomy
# which differs slightly from our bird list taxonomy
###################################################################

# fix taxonomy mismatch — Bleda eximia -> Bleda eximius
birds_ugf <- birds_ugf |>
  mutate(scientific_name = recode(scientific_name,
                                  "Bleda eximia" = "Bleda eximius"))

# inner join keeps only species present in both datasets
matched <- birds_ugf |>
  inner_join(aoh_list, by = c("scientific_name" = "BINOMIAL"))

cat("Matched:", nrow(matched), "of", nrow(birds_ugf), "species\n")

# check which species didn't match
unmatched <- birds_ugf |>
  anti_join(aoh_list, by = c("scientific_name" = "BINOMIAL"))

cat("Unmatched:", nrow(unmatched), "species\n")

###################################################################
# SECTION 4: Extract relevant AOH rasters from downloaded folders
# The AOH dataset (Lumbierres et al.) is organised by order and
# migratory status. This function copies only the rasters for
# our bird species from each downloaded folder into a single
# extracted/ folder, preferring _R (resident) over _B (breeding)
###################################################################

process_aoh_folder <- function(folder_path, species_list, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  matched_files <- c()
  for (sp in species_list) {
    sp_file <- gsub(" ", "_", sp)
    r_file <- file.path(folder_path, paste0(sp_file, "_R.tif"))  # resident
    b_file <- file.path(folder_path, paste0(sp_file, "_B.tif"))  # breeding
    no_suffix <- file.path(folder_path, paste0(sp_file, ".tif")) # non-migratory
    if (file.exists(r_file)) matched_files <- c(matched_files, r_file)
    else if (file.exists(b_file)) matched_files <- c(matched_files, b_file)
    else if (file.exists(no_suffix)) matched_files <- c(matched_files, no_suffix)
  }
  file.copy(matched_files, output_dir)
  cat("Copied", length(matched_files), "species rasters\n")
}

# list of downloaded AOH folders to process
folders <- c(
  "Birds_migratory_pterocliformes",
  "Birds_non_migratory_pelecaniformes",
  "Birds_non_migratory_psittaciformes",
  "Birds_non_migratory_pterocliformes",
  "Birds_non_migratory_struthioniformes",
  "Birds_non_migratory_strigiformes",
  "Birds_non_migratory_trogoniformes",
  "Birds_migratory_pelecaniformes",
  "Birds_non_migratory_piciformes"
)

for (folder in folders) {
  cat("Processing", folder, "\n")
  process_aoh_folder(
    folder_path = here("data/birds_aoh", folder),
    species_list = matched$scientific_name,
    output_dir = here("data/birds_aoh/extracted")
  )
}

###################################################################
# SECTION 5: Calculate % of global AOH within the UGF
# For each species:
#   1. Calculate total global AOH (sum of all pixels x pixel area)
#   2. Clip raster to UGF boundary and calculate local AOH
#   3. Express local AOH as % of global AOH (= Ps,i in STAR formula)
#   4. Save clipped raster and delete global raster to save space
###################################################################

process_species_aoh <- function(tif_path) {
  sp_name <- gsub("_R$|_B$", "", tools::file_path_sans_ext(basename(tif_path)))
  sp_name <- gsub("_", " ", sp_name)
  
  r <- rast(tif_path)
  
  # approximate pixel area in km² (faster than cellSize() for large rasters)
  pixel_area_km2 <- prod(res(r)) * (111.32^2)
  
  # total global AOH — sum of all habitat pixels x pixel area
  total_aoh_km2 <- global(r, "sum", na.rm = TRUE)$sum * pixel_area_km2
  
  # clip to UGF bounding box first, then mask to exact boundary
  r_clipped <- tryCatch({
    crop(r, ugf_vect)
  }, error = function(e) {
    cat("No overlap for", sp_name, "- skipping clip\n")
    return(NULL)
  })
  
  # if species has no habitat in UGF, record 0 and move on
  if (is.null(r_clipped)) {
    file.remove(tif_path)
    return(data.frame(
      scientific_name = sp_name,
      total_aoh_km2 = total_aoh_km2,
      ugf_aoh_km2 = 0,
      pct_aoh_in_ugf = 0
    ))
  }
  
  r_clipped <- mask(r_clipped, ugf_vect)
  ugf_aoh_km2 <- global(r_clipped, "sum", na.rm = TRUE)$sum * pixel_area_km2
  
  # save clipped raster before deleting global
  out_path <- file.path(here("data/birds_aoh/clipped"), basename(tif_path))
  writeRaster(r_clipped, out_path, overwrite = TRUE, datatype = "INT1U")
  
  # delete global raster to free up storage
  file.remove(tif_path)
  gc()
  
  return(data.frame(
    scientific_name = sp_name,
    total_aoh_km2 = total_aoh_km2,
    ugf_aoh_km2 = ugf_aoh_km2,
    pct_aoh_in_ugf = ugf_aoh_km2 / total_aoh_km2 * 100  # Ps,i for STAR formula
  ))
}

# run for all extracted rasters
tif_files <- list.files(here("data/birds_aoh/extracted"), pattern = "\\.tif$", full.names = TRUE)

results <- list()
for (i in seq_along(tif_files)) {
  cat("Processing", i, "of", length(tif_files), ":", basename(tif_files[i]), "\n")
  results[[i]] <- process_species_aoh(tif_files[i])
  gc()
}

# append new results to existing aoh_table (avoids overwriting previous batches)
new_results <- bind_rows(results)
aoh_table <- bind_rows(aoh_table, new_results)
write_csv(aoh_table, here("outputs/aoh_table.csv"))

###################################################################
# SECTION 6: Pull threat data from IUCN Red List API
# For each species, retrieves threat scope and severity from the
# most recent global assessment. These are used to calculate C
# (relative contribution of each threat to extinction risk)
# following the STAR methodology (Mair et al. 2021)
###################################################################

threat_list <- list()

for (i in 1:nrow(birds_star)) {
  sp <- birds_star$scientific_name[i]
  genus <- word(sp, 1)
  species <- word(sp, 2)
  
  result <- tryCatch({
    assessment <- rl_species(genus, species, key = api_key)
    # get most recent global assessment
    assessment_id <- assessment$assessments$assessment_id[
      assessment$assessments$latest == TRUE
    ][1]
    threats <- rl_assessment(assessment_id, key = api_key)$threats
    if (!is.null(threats) && nrow(threats) > 0) {
      threats$scientific_name <- sp
      threats
    } else NULL
  }, error = function(e) NULL)
  
  threat_list[[i]] <- result
  Sys.sleep(0.5)
  cat("Done", i, "of", nrow(birds_star), ":", sp, "\n")
}

threats_df <- bind_rows(threat_list)

# decline values from STAR Supplementary Table 2
# each cell = expected % population decline from that scope x severity combination
decline_lookup <- data.frame(
  scope = c(
    rep("Whole (>90%)", 6),
    rep("Majority (50-90%)", 6),
    rep("Minority (<50%)", 6)
  ),
  severity = rep(c(
    "Very Rapid Declines",
    "Rapid Declines", 
    "Slow, Significant Declines",
    "Negligible declines",
    "No decline",
    "Causing/Could cause fluctuations"
  ), 3),
  decline = c(
    63, 24, 10, 1, 0, 10,  # Whole (>90%)
    52, 18,  9, 0, 0,  9,  # Majority (50-90%)
    24,  7,  5, 0, 0,  5   # Minority (<50%)
  )
)

# calculate C = decline from this threat / total decline from all threats
# unknown severity treated as 0 (excluded from calculation)
threats_df <- threats_df |>
  select(-any_of(c("decline", "C", "scope_num", "severity_num"))) |>
  left_join(decline_lookup, by = c("scope", "severity")) |>
  filter(timing != "Past, Unlikely to Return") |>  # exclude historical threats
  mutate(decline = replace_na(decline, 0)) |>       # unknown severity = 0
  group_by(scientific_name) |>
  mutate(
    total_decline = sum(decline, na.rm = TRUE),
    C = ifelse(total_decline == 0, 0, decline / total_decline)
  ) |>
  ungroup()

write_csv(threats_df, here("outputs/threats_df.csv"))

###################################################################
# SECTION 7: Calculate STAR-T scores
# STAR-T formula (Mair et al. 2021):
#   T(t,i) = Σ P(s,i) x Ws x C(s,t)
# where:
#   P(s,i) = % of global AOH within UGF (pct_aoh_in_ugf / 100)
#   Ws     = IUCN Red List weight (NT=1, VU=2, EN=3, CR=4)
#   C(s,t) = relative contribution of threat t to species s extinction risk
###################################################################

# join IUCN status — remove duplicate columns first if join has been run before
aoh_table <- aoh_table |>
  left_join(birds_ugf |> select(scientific_name, IUCN, population_trend),
            by = "scientific_name")

# filter to threatened species with habitat in UGF
# LC and DD excluded per STAR methodology (Mair et al. 2021)
birds_star <- aoh_table |>
  filter(pct_aoh_in_ugf > 0) |>
  filter(IUCN %in% c("NT", "VU", "EN", "CR")) |>
  mutate(
    # assign IUCN weight
    Ws = case_when(
      IUCN == "NT" ~ 1,
      IUCN == "VU" ~ 2,
      IUCN == "EN" ~ 3,
      IUCN == "CR" ~ 4
    ),
    # species-level STAR-T = P(s,i) x Ws
    # (summing C across all threats = 1, so simplifies to this)
    star_t = (pct_aoh_in_ugf / 100) * Ws
  )

# calculate STAR-T per species per threat for threat breakdown
star_data <- threats_df |>
  left_join(birds_star |> select(scientific_name, pct_aoh_in_ugf, Ws),
            by = "scientific_name") |>
  mutate(star_t = (pct_aoh_in_ugf / 100) * Ws * C)

# aggregate STAR-T by threat type to identify priority threats in UGF
star_by_threat <- star_data |>
  group_by(en, code) |>
  summarise(star_t = sum(star_t, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(star_t))

star_by_threat |> print(n = Inf)

# total STAR-T for the UGF
total_star_t <- sum(star_by_threat$star_t)
cat("Total STAR-T for UGF:", total_star_t, "\n")

write_csv(birds_star, here("outputs/birds_star_t.csv"))

###################################################################
# SECTION 8: Build priority raster
# Each cell receives the sum of STAR-T scores for all species
# whose habitat (AOH) overlaps that cell. Higher values indicate
# areas where conserving habitat would contribute more to reducing
# extinction risk of threatened bird species in the UGF.
###################################################################

# create blank template raster covering full UGF extent
# resolution matches original AOH rasters (~100m)

template <- rast(ext(ugf_vect), resolution = 0.0009920722, crs = "EPSG:4326")
values(template) <- 0
priority_raster <- template

for (i in 1:nrow(birds_star_t)) {
  sp <- birds_star_t$scientific_name[i]
  star_score <- birds_star_t$star_t[i]
  
  # find clipped raster for this species
  sp_file <- file.path(clipped_dir, paste0(gsub(" ", "_", sp), "_R.tif"))
  if (!file.exists(sp_file)) {
    sp_file <- file.path(clipped_dir, paste0(gsub(" ", "_", sp), ".tif"))
  }
  
  if (file.exists(sp_file)) {
    cat("Processing", i, "of", nrow(birds_star_t), ":", sp, "\n")
    r <- rast(sp_file)
    # reproject to match template grid
    r_proj <- project(r, template, method = "near")
    # force binary — 1 = habitat present, 0 = absent
    r_proj <- ifel(is.na(r_proj), 0, 1)
    # add weighted habitat layer to accumulating priority raster
    priority_raster <- priority_raster + (r_proj * star_score)
    gc()
  }
}

# mask to UGF boundary — areas outside boundary set to NA
priority_raster <- mask(priority_raster, ugf_vect)

# normalize to 0-1 for easier interpretation
priority_norm <- (priority_raster - global(priority_raster, "min", na.rm = TRUE)$min) / 
  (global(priority_raster, "max", na.rm = TRUE)$max - global(priority_raster, "min", na.rm = TRUE)$min)

# save outputs
writeRaster(priority_raster, here("outputs/bird_priority_raster.tif"), overwrite = TRUE)
writeRaster(priority_norm, here("outputs/bird_priority_normalized.tif"), overwrite = TRUE)

png(here("outputs/bird_priority_raster.png"), width = 2000, height = 1200, res = 150)
plot(priority_raster)
dev.off()

png(here("outputs/bird_priority_normalized.png"), width = 2000, height = 1200, res = 150)
plot(priority_norm)
dev.off()

###################################################################
# SECTION 9: Plot individual species AOH maps
# Shows the clipped AOH raster for each threatened bird species
# in the UGF — useful for visualizing spatial coverage
###################################################################

png(here("outputs/species_aoh_plots.png"), width = 3000, height = 1500, res = 150)
par(mfrow = c(3, 6))  # 3 rows x 6 columns for 18 species

for (i in 1:nrow(birds_star)) {
  sp <- birds_star$scientific_name[i]
  sp_file <- file.path(clipped_dir, paste0(gsub(" ", "_", sp), "_R.tif"))
  if (!file.exists(sp_file)) {
    sp_file <- file.path(clipped_dir, paste0(gsub(" ", "_", sp), ".tif"))
  }
  if (file.exists(sp_file)) {
    r <- rast(sp_file)
    coltab(r) <- NULL  # remove embedded color table to allow custom colors
    plot(r, main = sp, cex.main = 0.6, legend = FALSE, col = "purple")
  }
}

par(mfrow = c(1, 1))  # reset plot layout
dev.off()