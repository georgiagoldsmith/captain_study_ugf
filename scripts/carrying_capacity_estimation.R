# Calculating carrying capacities
# G Goldsmith
# August 2026
###############################################################################

# Carrying-capacity scalar based on Damuth's Law:
# Fit an order-stratified density~body-mass model using real empirical
# density records (TetraDENSITY 2.0, Santini et al. 2024) and body mass
# (AVONET, Tobias et al. 2022), then predict K for the 18 UGF species --
# using TetraDENSITY density estimates directly where a species has a match,
# and the fitted model only for the rest.

# Common-slope, taxonomic order-specific-intercept model (ANCOVA): 
# log10(density) ~ log10(mass) + Order, fit on ALL ~1850 bird species with both 
# a density record and a body mass (not just the 18 target species). 

library(readr)
library(readxl)
library(dplyr)
library(stringr)

# raw inputs: real density records (all tetrapods, worldwide), body mass
# (all birds, worldwide), and our 18 UGF species + their taxonomic orders
td <- read_csv("/tmp/tetradensity.csv", show_col_types = FALSE)
avonet <- read_excel("/tmp/avonet.xlsx", sheet = "AVONET1_BirdLife")
traits <- read_csv(here::here("data/Birds species status and traits/species_traits.csv"), show_col_types = FALSE)
orders <- read_csv(here::here("outputs/bird_orders.csv"), show_col_types = FALSE)

# --- STEP 1: build the training set ---

# keep birds only; standardize the two density units TetraDENSITY reports
# birds in (pairs count 2 individuals) so every row is directly comparable
birds <- td |> filter(Class == "Aves") |>
  mutate(density_ind_km2 = ifelse(Density_unit == "pairs/km2", Density * 2, Density))

# a species can have many separate density records (different studies/sites)
# collapse to one row per species using the median across all of them
sp_density <- birds |> group_by(Species_rep) |>
  summarise(density_ind_km2 = median(density_ind_km2, na.rm = TRUE), n_records = n(), .groups = "drop")

avonet_sm <- avonet |> select(Species1, Order1, Mass) |> filter(!is.na(Mass))

# inner_join keeps only species present in BOTH datasets (i.e. species we
# both have a real measured density for AND know the body mass of)
train <- sp_density |>
  inner_join(avonet_sm, by = c("Species_rep" = "Species1")) |>
  filter(density_ind_km2 > 0)

cat("Training set: n =", nrow(train), "bird species with density + mass\n")

# --- STEP 2: fit ONE shared line across all orders (common slope), letting
# each order shift the intercept ---
train$logD <- log10(train$density_ind_km2)
train$logM <- log10(train$Mass)
train$Order1 <- factor(train$Order1)

fit <- lm(logD ~ logM + Order1, data = train)
cat("\nGlobal slope (log-log):", round(coef(fit)["logM"], 3),
    "(textbook Damuth exponent is -0.75, for comparison)\n")

# --- STEP 3: predict density for our 18 species from mass + order alone,
# using the line fit in step 2 (none of these 18 were in the training set) ---
sp18 <- traits |>
  select(species, biomass) |>
  left_join(orders |> select(scientific_name, order), by = c("species" = "scientific_name")) |>
  rename(Mass = biomass, Order1 = order)

sp18$Order1 <- str_to_title(tolower(sp18$Order1))  # match AVONET's "Passeriformes" casing
sp18$Order1 <- factor(sp18$Order1, levels = levels(train$Order1))

missing_order <- sp18$Order1[is.na(sp18$Order1)]
if (length(missing_order) > 0) cat("\nWARNING: order not in training data for some species -- will fall back to global slope only\n")

sp18$logM <- log10(sp18$Mass)
pred <- predict(fit, newdata = sp18, se.fit = TRUE)
sp18$density_pred_ind_km2 <- 10 ^ pred$fit

# --- STEP 4: exception -- for any of our 18 that DOES have a real
# TetraDENSITY record, use that measured value instead of the prediction ---
direct <- train |> filter(Species_rep %in% sp18$species) |>
  select(species = Species_rep, density_direct = density_ind_km2, n_records)

sp18 <- sp18 |> left_join(direct, by = "species")
sp18$density_final <- ifelse(!is.na(sp18$density_direct), sp18$density_direct, sp18$density_pred_ind_km2)
sp18$source <- ifelse(!is.na(sp18$density_direct), "TetraDENSITY direct match", "order-stratified model")

# --- STEP 5: convert density (individuals/km2) into individuals per
# CAPTAIN grid cell. Cells are 3km x 3km = 9 km2, so density is multiplied by 9 
# for per cell count ---
sp18$K_tetradensity <- round(sp18$density_final * 9, 1)

out <- sp18 |> select(species, Order1, Mass, density_final, source, K_tetradensity) |>
  arrange(desc(K_tetradensity))

write_csv(out, here::here("outputs/K_tetradensity.csv"))
print(out, n = Inf)
