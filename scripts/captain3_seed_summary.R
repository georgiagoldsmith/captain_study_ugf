# Cross-seed summaries of the corrected 3km sweep: contact sheet, 10-policy
# selection frequency, and top-3 agreement.
#
#   Rscript scripts/captain3_seed_summary.R                  # all three
#   Rscript scripts/captain3_seed_summary.R frequency top3   # named subset
#
# Merged from captain3_seed_grid_corrected.R,
# captain3_selection_frequency_corrected.R and captain3_top3_agreement.R.
# All three re-declared the SEEDS vector and repeated the same "read
# selected_grid.csv, drop the locked-in PAs, accumulate" loop -- adding a seed
# to the sweep used to mean editing three files. Output filenames unchanged.
#
# Locked-in protected areas are held out of every count below: each policy
# contains them by construction, so including them would swamp the gradients.
suppressMessages({library(here);library(terra);library(sf)})

DD    <- "/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/captain2-main/ugf_data_3km_corrected"
SEEDS <- c(42,675,12345,41898,176843,367181179,915315400,232373165,3245678920,8052026)
run_dir <- function(s) paste0("outputs/captain3_seed", s, "_corrected")

b  <- vect(st_transform(st_read(here("data/UGF_gp.shp","UGF_gp.shp"),quiet=TRUE),3857))
co <- mask(rast(file.path(DD,"environmental_layers","costs.tif")), b, touches=TRUE)
pa <- as.matrix(rast(file.path(DD,"environmental_layers","protected_areas.tif")), wide=TRUE); pa[is.na(pa)] <- 0

as_rast <- function(m) mask(rast(m, extent=ext(co), crs=crs(co)), co)
grid_of <- function(s) as.matrix(read.csv(here(run_dir(s),"selected_grid.csv"), header=FALSE))
new_of  <- function(s) (grid_of(s) > 0) & (pa != 1)          # new selections only

# final-epoch avg_reward from a seed's RL training log
reward_of <- function(s) {
  L <- readLines(here(run_dir(s),"training_log.tsv"))
  as.numeric(strsplit(L[length(L)],"\t")[[1]][3])
}
# per-cell count of how many of `seeds` selected the cell
accumulate <- function(seeds) Reduce(`+`, lapply(seeds, function(s) new_of(s) * 1))

FIGS <- list(

  # 2 cols x 5 rows, matching the UGF's ~2.5:1 aspect ratio
  grid = function() {
    png(here("outputs","captain3_seed_grid_corrected.png"), width=2200, height=1950, res=150)
    layout(matrix(1:12, nrow=6, byrow=TRUE), heights=c(rep(1,5),0.34))
    for (i in seq_along(SEEDS)) {
      s <- SEEDS[i]
      L <- readLines(here(run_dir(s),"training_log.tsv")); f <- strsplit(L[length(L)],"\t")[[1]]
      r <- as_rast(ifelse(pa == 1, 2, ifelse(grid_of(s) > 0, 1, 0)))
      par(mar=c(0.1,0.3,1.6,0.3))
      # NB panels are labelled by ORDINAL (1..10), not by seed value -- carried
      # over from captain3_seed_grid_corrected.R unchanged. Swap `i` for `s`
      # below to label each panel with its actual seed.
      plot(r, col=c("grey92","#1f78b4","grey40"), breaks=c(-0.5,0.5,1.5,2.5),
           main=sprintf("Seed %d   reward %.3f   LC %.2f", i, as.numeric(f[3]), as.numeric(f[12])),
           cex.main=1.25, legend=FALSE, axes=FALSE, mar=c(0.1,0.3,1.6,0.3))
    }
    par(mar=c(0,0,0,0)); plot.new()
    legend("center", legend=c("Not selected","Newly selected","Protected Areas"),
           fill=c("grey92","#1f78b4","grey40"), ncol=3, bty="n", cex=1.35, xpd=TRUE)
    par(mar=c(0,0,0,0)); plot.new()
    dev.off(); cat("contact sheet written\n")
  },

  # each cell coloured by how many of the 10 policies selected it, 1..10
  frequency = function() {
    acc <- accumulate(SEEDS)
    writeRaster(as_rast(acc), here("outputs","captain3_selection_frequency_corrected.tif"), overwrite=TRUE)

    # display: PAs into a sentinel below zero so they read as their own flat class
    disp <- as_rast(ifelse(pa == 1, -1, acc))
    grad <- hcl.colors(10, "YlGnBu", rev=TRUE)
    cols <- c("grey40","grey93", grad)
    brks <- c(-1.5,-0.5, seq(0.5, 10.5, by=1))

    png(here("outputs","captain3_selection_frequency_corrected.png"), width=2200, height=1250, res=150)
    # tighter legend panel, and a taller top margin so the title sits well clear of the map
    layout(matrix(1:2, nrow=2), heights=c(6,1.15))
    par(mar=c(0.2,2.2,6.5,2.2))
    plot(disp, col=cols, breaks=brks, legend=FALSE, main="", cex.main=1.9)
    title(main="Selection frequency across 10 policies (corrected 3km)", cex.main=1.9, line=4.2)
    par(mar=c(0,0,0,0)); plot.new()
    # title drawn separately -- as a legend() title it sits above the panel and clips
    text(0.5, 0.93, "number of policies selecting the cell", cex=1.3, xpd=TRUE)
    legend(x=0.5, y=0.60, xjust=0.5, yjust=0.5,
           legend=c("never", as.character(1:10)), fill=c("grey93", grad),
           ncol=11, bty="n", cex=1.3, xpd=TRUE)
    legend(x=0.5, y=0.16, xjust=0.5, yjust=0.5,
           legend="Protected Areas (locked in, in every policy)", fill="grey40",
           ncol=1, bty="n", cex=1.4, xpd=TRUE)
    dev.off()

    tb <- table(factor(acc[acc > 0 & pa != 1], levels=1:10))
    cat("cells by number of policies selecting them:\n")
    for (i in 1:10) cat(sprintf("  %2d/10  %6d\n", i, tb[[i]]))
    cat(sprintf("\nunion %d | all-10 %d (%.0f%%)\n", sum(acc > 0 & pa != 1), tb[[10]],
                100 * tb[[10]] / sum(acc > 0 & pa != 1)))
  },

  # agreement among the three best-rewarded policies
  top3 = function() {
    rew  <- vapply(SEEDS, reward_of, 0)
    top3 <- SEEDS[order(-rew)][1:3]
    cat("top 3 by reward:", paste(sprintf("%s (%.3f)", top3, sort(rew, decreasing=TRUE)[1:3]),
        collapse=", "), "\n")

    acc <- accumulate(top3)
    writeRaster(as_rast(acc), here("outputs","captain3_top3_agreement.tif"), overwrite=TRUE)

    disp <- as_rast(ifelse(pa == 1, -1, acc))
    cols <- c("grey40","grey93","#fdd49e","#41b6c4","#0b3d91")     # PA, never, 1, 2, 3
    brks <- c(-1.5,-0.5,0.5,1.5,2.5,3.5)
    n1 <- sum(acc == 1 & pa != 1); n2 <- sum(acc == 2 & pa != 1); n3 <- sum(acc == 3 & pa != 1)
    uni <- n1 + n2 + n3
    # size of a single plan, for the "% of one plan" figure -- computed from the
    # top-3 runs rather than the 4827 that used to be hardcoded here
    plan <- round(mean(vapply(top3, function(s) sum(new_of(s)), 0)))

    png(here("outputs","captain3_top3_agreement.png"), width=2200, height=1250, res=150)
    layout(matrix(1:2, nrow=2), heights=c(6,1.15))
    par(mar=c(0.2,2.2,6.5,2.2))
    plot(disp, col=cols, breaks=brks, legend=FALSE, main="", axes=TRUE)
    title(main=sprintf("Agreement among the 3 best-rewarded policies (IoU = %.2f)", n3/uni),
          cex.main=1.9, line=4.2)
    par(mar=c(0,0,0,0)); plot.new()
    text(0.5, 0.93, "number of the 3 policies selecting the cell", cex=1.3, xpd=TRUE)
    legend(x=0.5, y=0.60, xjust=0.5, yjust=0.5,
           legend=c("never", sprintf("1  (%s)", format(n1, big.mark=",")),
                    sprintf("2  (%s)", format(n2, big.mark=",")),
                    sprintf("3  (%s)", format(n3, big.mark=","))),
           fill=cols[-1], ncol=4, bty="n", cex=1.4, xpd=TRUE)
    legend(x=0.5, y=0.16, xjust=0.5, yjust=0.5,
           legend="Protected Areas (locked in, in every policy)", fill="grey40",
           ncol=1, bty="n", cex=1.4, xpd=TRUE)
    dev.off()
    cat(sprintf("1/3 %d | 2/3 %d | 3/3 %d | union %d -> all-three %.0f%% of union, %.0f%% of one plan\n",
                n1, n2, n3, uni, 100*n3/uni, 100*n3/plan))
  }
)

sel <- commandArgs(trailingOnly = TRUE)
if (!length(sel)) sel <- names(FIGS)
unknown <- setdiff(sel, names(FIGS))
if (length(unknown)) stop("unknown section(s): ", paste(unknown, collapse=", "),
                          "\navailable: ", paste(names(FIGS), collapse=", "))
for (nm in sel) { cat("---", nm, "\n"); FIGS[[nm]]() }
