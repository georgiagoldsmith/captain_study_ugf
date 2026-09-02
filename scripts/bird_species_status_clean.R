
library(rredlist)
library(readxl)
library(here)
library(tidyverse)
library(sf)

###################################################################

# Assign IUCN status to UGF birds

Sys.setenv(IUCN_REDLIST_KEY = "EXJCZgc5cxkksBJ9fqKHnR6aEqqugamF8ZrK") 

# Load data 
checkpoint_path <- "birds_progress.csv"

if (file.exists(checkpoint_path)) {
  birds <- read.csv(checkpoint_path)
  cat("Resuming from checkpoint:", checkpoint_path, "\n")
} else {
  birds <- read_excel(here("data", "birds_ugf.xlsx"))
  birds <- birds |>
    separate(scientific_name, into = c("genus", "species"),
             sep = " ", remove = FALSE)
  birds$IUCN           <- NA
  birds$population_trend <- NA
}

# Main loop 
for (i in seq_len(nrow(birds))) {
  
  # Skip rows already filled (allows safe resume)
  if (!is.na(birds$IUCN[i])) next
  
  cat("Species", i, "of", nrow(birds), "-", birds$scientific_name[i], "\n")
  
  out <- tryCatch(
    rl_species(
      genus  = birds$genus[i],
      species = birds$species[i],
      key    = Sys.getenv("IUCN_REDLIST_KEY")
    ),
    error = function(e) {
      cat("  API error for", birds$scientific_name[i], ":", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (!is.null(out) && nrow(out$assessments) > 0) {
    
    # Filter to most recent global assessment
    # Check the scopes string matches what the API actually returns
    global_latest <- out$assessments[
      out$assessments$latest == TRUE & grepl("Global", out$assessments$scopes), 
    ]
    
    if (nrow(global_latest) > 0) {
      birds$IUCN[i] <- global_latest$red_list_category_code[1]
      
      assess_detail <- tryCatch(
        rl_assessment(
          id  = global_latest$assessment_id[1],
          key = Sys.getenv("IUCN_REDLIST_KEY")
        ),
        error = function(e) {
          cat("  Assessment detail error for", birds$scientific_name[i], "\n")
          NULL
        }
      )
      
      if (!is.null(assess_detail)) {
        birds$population_trend[i] <- assess_detail$population_trend$description$en
      }
      
    } else {
      cat("  No global assessment found for", birds$scientific_name[i], "\n")
    }
    
  } else {
    cat("  No results found for", birds$scientific_name[i], "\n")
  }
  
  # Rate limiting — be polite to the API
  Sys.sleep(0.5)
  
  # Save checkpoint every 25 species
  if (i %% 25 == 0) {
    write.csv(birds, checkpoint_path, row.names = FALSE)
    cat("  Progress saved at species", i, "\n")
  }
}

# Final save 
write.csv(birds, "birds_species_status.csv", row.names = FALSE)
cat("Done! Results saved to birds_species_status.csv\n")

# Check how many were matched 
cat("IUCN status assigned:", sum(!is.na(birds$IUCN)), "of", nrow(birds), "species\n")
missing <- birds[is.na(birds$IUCN), "scientific_name"]
if (length(missing) > 0) {
  cat("Missing matches — check these names against IUCN taxonomy:\n")
  print(missing)
}

install.packages("writexl")
library(writexl)

write_xlsx(birds, here("outputs/bird_species_status.xlsx"))

api_key <- "EXJCZgc5cxkksBJ9fqKHnR6aEqqugamF8ZrK"

rl_version(key = api_key)
rl_species("Accipiter", "badius", key = api_key)

# test
rl_threats("Accipiter", "badius", key = api_key)

################################################################################################

# AOH list


library(readr)

birds_ugf <- read_csv(here("data/Birds w species status/birds_species_status.csv"))
aoh_list <- read_csv(here("data/birds_aoh/Birds_list_AOH.csv"))

# check what columns are available
names(aoh_list)
names(birds_ugf)

# check what orders are present
unique(aoh_list$order)  # column name may differ, check names() first

# join birds in the AOH list to the list of CI SDMs w IUCN status 
matched <- birds_ugf |>
  inner_join(aoh_list, by = c("scientific_name" = "BINOMIAL"))  

# how many matched
nrow(matched)

# which didn't match
unmatched <- birds |>
  anti_join(aoh_list, by = c("scientific_name" = "BINOMIAL"))

nrow(unmatched)

unmatched$scientific_name

# check for different names
iucn_lookup <- lapply(1:nrow(unmatched), function(i) {
  result <- tryCatch(
    rl_species(unmatched$genus[i], unmatched$species[i], key = api_key),
    error = function(e) NULL
  )
  if (!is.null(result)) {
    data.frame(
      original_name = unmatched$scientific_name[i],
      iucn_name = result$taxon$scientific_name,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(original_name = unmatched$scientific_name[i], iucn_name = NA)
  }
})

iucn_lookup_df <- bind_rows(iucn_lookup)
print(iucn_lookup_df)

genus_lookup <- lapply(1:nrow(unmatched), function(i) {
  result <- tryCatch(
    rl_species(unmatched$genus[i], key = api_key),
    error = function(e) NULL
  )
  if (!is.null(result)) {
    data.frame(
      original_name = unmatched$scientific_name[i],
      iucn_matches = paste(result$assessments$taxon_scientific_name, collapse = "; "),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(original_name = unmatched$scientific_name[i], iucn_matches = NA)
  }
})

genus_lookup_df <- bind_rows(genus_lookup)
print(genus_lookup_df)

synonym_lookup <- lapply(1:nrow(unmatched), function(i) {
  result <- tryCatch(
    rl_synonyms(unmatched$genus[i], unmatched$species[i], key = api_key),
    error = function(e) NULL
  )
  if (!is.null(result)) {
    data.frame(
      original_name = unmatched$scientific_name[i],
      accepted_name = result$result$accepted_name,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(original_name = unmatched$scientific_name[i], accepted_name = NA)
  }
})

synonym_lookup_df <- bind_rows(synonym_lookup)
print(synonym_lookup_df)

excluded <- read_csv(here("data/birds_aoh/Birds_list_excluded.csv"))
unmatched$scientific_name %in% excluded$BINOMIAL

# check if any appear under slightly different names in aoh_list
aoh_list$BINOMIAL[grep("Alethe", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Anthreptes", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Bleda", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Bubo", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Cossypha", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Canirallus", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Caprimulgus", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Cinnyricinclus", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Coracina", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Ispidina", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Muscicapa", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Petrochelidon", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Phyllastrephus", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Prinia", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Serinus", aoh_list$BINOMIAL)]
aoh_list$BINOMIAL[grep("Zoothera", aoh_list$BINOMIAL)]

# fix the name mismatch
birds_ugf <- birds_ugf |>
  mutate(scientific_name = recode(scientific_name,
                                  "Bleda eximia" = "Bleda eximius"))

# re-run the join
matched <- birds_ugf |>
  inner_join(aoh_list, by = c("scientific_name" = "BINOMIAL"))

matched |>
  count(Order_) |>
  arrange(desc(n))

matched |>
  count(Order_, is.migratory) |>
  arrange(Order_, is.migratory) |>
  print( n = Inf)

#################################################################################################################

# Extract AOH rasters


library(terra)

# save only tifs that match the birds ugf dataset
process_aoh_folder <- function(folder_path, species_list, output_dir) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # list all tifs in the folder
  all_tifs <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)
  
  matched_files <- c()
  for (sp in species_list) {
    sp_file <- gsub(" ", "_", sp)
    r_file <- file.path(folder_path, paste0(sp_file, "_R.tif"))
    b_file <- file.path(folder_path, paste0(sp_file, "_B.tif"))
    no_suffix <- file.path(folder_path, paste0(sp_file, ".tif"))
    
    if (file.exists(r_file)) matched_files <- c(matched_files, r_file)
    else if (file.exists(b_file)) matched_files <- c(matched_files, b_file)
    else if (file.exists(no_suffix)) matched_files <- c(matched_files, no_suffix)
  }
  
  # copy matched files to output
  file.copy(matched_files, output_dir)
  cat("Copied", length(matched_files), "species rasters\n")
}

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

###################################################################################################

# % habitat in UGF


## folder the extracted rasters are
accipiter_badius_aoh <- rast(here("data/birds_aoh/extracted/Accipiter_badius_R.tif"))
crs(accipter_badius_aoh)
st_crs(ugf_boundary)

## reproject ugf boundary to match AOH raster CRS (WGS84)
ugf_boundary <- st_transform(ugf_boundary, 4326)
ugf_vect <- vect(ugf_boundary)

## list all tifs
tif_files <- list.files(here("data/birds_aoh/extracted"), pattern = "\\.tif$", full.names = TRUE)
length(tif_files)
tif_files[1]

## function to process one raster
process_species_aoh <- function(tif_path) {
  
  sp_name <- gsub("_R$|_B$", "", tools::file_path_sans_ext(basename(tif_path)))
  sp_name <- gsub("_", " ", sp_name)
  
  r <- rast(tif_path)
  pixel_areas <- cellSize(r, unit = "km")
  
  # total global AOH
  total_aoh_km2 <- global(r * pixel_areas, "sum", na.rm = TRUE)$sum
  
  # clip to UGF
  r_clipped <- tryCatch({
    crop(r, ugf_vect)
  }, error = function(e) {
    cat("No overlap for", sp_name, "- skipping clip\n")
    return(NULL)
  })
  
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
  
  # UGF AOH
  ugf_aoh_km2 <- global(r_clipped * crop(pixel_areas, ugf_vect), "sum", na.rm = TRUE)$sum
  
  # save clipped raster FIRST
  out_path <- file.path(here("data/birds_aoh/clipped"), basename(tif_path))
  writeRaster(r_clipped, out_path, overwrite = TRUE)
  
  # THEN delete global raster
  file.remove(tif_path)
  
  return(data.frame(
    scientific_name = sp_name,
    total_aoh_km2 = total_aoh_km2,
    ugf_aoh_km2 = ugf_aoh_km2,
    pct_aoh_in_ugf = ugf_aoh_km2 / total_aoh_km2 * 100
  ))
}

test <- process_species_aoh(tif_files[1])
print(test)

## run for all species
tif_files <- list.files(here("data/birds_aoh/extracted"), pattern = "\\.tif$", full.names = TRUE)

results <- list()
for (i in seq_along(tif_files)) {
  cat("Processing", i, "of", length(tif_files), ":", basename(tif_files[i]), "\n")
  results[[i]] <- process_species_aoh(tif_files[i])
  gc()
}

new_results <- bind_rows(results)

# add to existing table instead of overwriting
aoh_table <- bind_rows(aoh_table, new_results)

# save 
write_csv(aoh_table, here("outputs/aoh_table.csv"))

aoh_table |>
  left_join(matched |> select(scientific_name, Order_, is.migratory),
            by = "scientific_name") |>
  count(Order_, is.migratory)


length(list.files(here("data/birds_aoh/clipped")))

##############################################################################

# Get threats for all species in birds_star



# get assessment ID for a species first
assessment <- rl_species("Accipiter", "badius", key = api_key)
assessment_id <- assessment$assessments$assessment_id[1]  # get latest

# then get threats for that assessment
rl_assessment(assessment_id, key = api_key)

# pull threats for all birds_star species
threat_list <- list()

for (i in 1:nrow(birds_star)) {
  sp <- birds_star$scientific_name[i]
  genus <- word(sp, 1)
  species <- word(sp, 2)
  
  result <- tryCatch({
    assessment <- rl_species(genus, species, key = api_key)
    
    # get latest global assessment id
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

names(threats_df)
head(threats_df)

# check for species missing threat data
birds_star |>
  filter(!scientific_name %in% threats_df$scientific_name) |>
  select(scientific_name, IUCN)

unique(threats_df$scope)
unique(threats_df$severity)

# values from STAR Supplementary Table 2
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
    63, 24, 10, 1, 0, 10,  # Whole
    52, 18,  9, 0, 0,  9,  # Majority
    24,  7,  5, 0, 0,  5   # Minority
  )
)

# join to threats_df
threats_df <- threats_df |>
  select(-any_of(c("decline", "C", "scope_num", "severity_num"))) |>
  left_join(decline_lookup, by = c("scope", "severity")) |>
  filter(timing != "Past, Unlikely to Return") |>
  group_by(scientific_name) |>
  mutate(
    total_decline = sum(decline, na.rm = TRUE),
    C = decline / total_decline
  ) |>
  ungroup()

# which species are missing
missing_threats <- birds_star |>
  filter(!scientific_name %in% threats_df$scientific_name) |>
  pull(scientific_name)

length(missing_threats)
head(missing_threats)

threat_list <- list()

for (i in seq_along(missing_threats)) {
  sp <- missing_threats[i]
  genus <- word(sp, 1)
  species <- word(sp, 2)
  
  result <- tryCatch({
    assessment <- rl_species(genus, species, key = api_key)
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
  cat("Done", i, "of", length(missing_threats), ":", sp, "\n")
}

new_threats <- bind_rows(threat_list)
threats_df <- bind_rows(threats_df, new_threats)

# Assign zero to threats with unknown severity 
threats_df |> filter(is.na(decline)) |> count(scope, severity)

threats_df <- threats_df |>
  mutate(decline = replace_na(decline, 0)) |>
  group_by(scientific_name) |>
  mutate(
    total_decline = sum(decline, na.rm = TRUE),
    C = ifelse(total_decline == 0, 0, decline / total_decline)
  ) |>
  ungroup()

# save immediately
write_csv(threats_df, here("outputs/threats_df.csv"))

################################################################################

# Calculating  STAR

####################
## IUCN score
## join species status
aoh_table <- aoh_table |>
  left_join(birds_ugf |> select(scientific_name, IUCN, population_trend),
            by = "scientific_name")


## filter out birds with no overlap
no_overlap_ugf <- aoh_table |>
  filter(pct_aoh_in_ugf == 0)

birds_aoh_ugf <- aoh_table |>
  filter(pct_aoh_in_ugf > 0)

# filter out birds of least concern
birds_star <- aoh_table |>
  filter(pct_aoh_in_ugf > 0) |>
  filter(IUCN %in% c("NT", "VU", "EN", "CR"))

# add weights for IUCN status
birds_star <- birds_star |>
  mutate(Ws = case_when(
    IUCN == "NT" ~ 1,
    IUCN == "VU" ~ 2,
    IUCN == "EN" ~ 3,
    IUCN == "CR" ~ 4,
  ))

###############
## Threat score
## join threat data to birds_star
star_data <- threats_df |>
  left_join(birds_star |> select(scientific_name, pct_aoh_in_ugf, Ws),
            by = "scientific_name")

###############
## STAR-T
## calculate STAR-T per species per threat
star_data <- star_data |>
  mutate(
    star_t = (pct_aoh_in_ugf / 100) * Ws * C
  )

## total STAR-T for UGF by threat
star_by_threat <- star_data |>
  group_by(description$en, code) |>
  summarise(star_t = sum(star_t, na.rm = TRUE)) |>
  arrange(desc(star_t))
star_by_threat |>
  arrange(desc(star_t)) |>
  print(n = Inf)

# total STAR-T for UGF overall
total_star_t <- sum(star_by_threat$star_t)
total_star_t

birds_star <- star_data |>
  group_by(scientific_name) |>
  summarise(star_t = sum(star_t, na.rm = TRUE)) |>
  left_join(birds_star |> select(-any_of("star_t")), by = "scientific_name") |>
  mutate(star_t = round(star_t, 6))

write.csv(birds_star, here("outputs/birds_star_t.csv"))

################################################################################
#some checks
# species with no UGF overlap
no_overlap <- aoh_table |>
  filter(ugf_aoh_km2 == 0) |>
  left_join(matched |> select(scientific_name, Order_, is.migratory),
            by = "scientific_name")

no_overlap |> select(scientific_name, Order_, is.migratory)
not_yet <- matched |>
  filter(!scientific_name %in% aoh_table$scientific_name) |>
  left_join(matched |> select(scientific_name, Order_, is.migratory),
            by = "scientific_name")

not_yet |>
  filter(Order_.x == "ANSERIFORMES")

length(unique(aoh_table$scientific_name))

# which matched species are still missing
matched |>
  filter(!scientific_name %in% aoh_table$scientific_name) |>
  nrow()

# duplicates?
sum(duplicated(aoh_table$scientific_name))

# species in aoh_table with ugf_aoh_km2 > 0
aoh_table |>
  filter(ugf_aoh_km2 > 0) |>
  nrow()

# check for NAs in ugf_aoh_km2
sum(is.na(aoh_table$ugf_aoh_km2))

# see the range of values
summary(aoh_table$ugf_aoh_km2)

# species with ugf_aoh_km2 = 0
aoh_table |>
  filter(ugf_aoh_km2 == 0) |>
  nrow()







# how many unique species in threats_df
length(unique(threats_df$scientific_name))

# which birds_star species are missing from threats_df
birds_star |>
  filter(!scientific_name %in% threats_df$scientific_name) |>
  select(scientific_name, IUCN)





summary(birds_star_t$pct_aoh_in_ugf)

birds_star_t <- birds_star_t |>
  select(-any_of(c("star_t.x", "star_t.y"))) |>
  mutate(star_t = (pct_aoh_in_ugf / 100) * Ws)

birds_star_t |> 
  select(scientific_name, IUCN, pct_aoh_in_ugf, Ws, star_t) |> 
  arrange(desc(star_t))


















