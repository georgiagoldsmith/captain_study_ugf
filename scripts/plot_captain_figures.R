# Plot CAPTAIN 3 solution
# G Goldsmith 
# August 2026
###############################################################################

suppressMessages({library(here);library(terra);library(sf)})

DD   <- "/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/captain2-main/ugf_data_3km_corrected/environmental_layers"
SEED <- 232373165
RUN_PRESENT <- "outputs/captain3_seed232373165_corrected"
RUN_FUTURE  <- "outputs/captain3_seed232373165_corrected_future"
RUN_STAGED  <- "outputs/captain3_seed232373165_corrected_future_staged"

ugf <- vect(st_transform(st_read(here("data/UGF_gp.shp","UGF_gp.shp"),quiet=TRUE),3857))
co  <- mask(rast(file.path(DD,"costs.tif")), ugf, touches=TRUE)
pa  <- as.matrix(rast(file.path(DD,"protected_areas.tif")), wide=TRUE); pa[is.na(pa)] <- 0

SOL_COLS  <- c("grey93","#1f78b4","grey40")
PAIR_COLS <- c("grey93","#e66101","#2c7bb6","#1a9641","grey40")
fmt <- function(x) format(x, big.mark=",")

grid_of <- function(d) as.matrix(read.csv(here(d,"selected_grid.csv"), header=FALSE))
# new selections only -- PAs are locked in for every run, so they are not a choice
new_of  <- function(d) (grid_of(d) > 0) & (pa != 1)
as_rast <- function(m) mask(rast(m, extent=ext(co), crs=crs(co)), co)

# --- one solution on its own ---
# 0 not selected, 1 newly selected, 2 pre-existing protected area. PAs are
# tested FIRST so they read as their own class -- they are all selected, being
# locked in, but they are not something the model chose.
plot_solution <- function(run, outfile, title_text) {
  g <- grid_of(run)
  n <- sum(g > 0 & pa != 1)
  r <- as_rast(ifelse(pa == 1, 2, ifelse(g > 0, 1, 0)))
  png(here("outputs", outfile), width=2200, height=1150, res=150)
  layout(matrix(1:2, nrow=2), heights=c(6,0.85))
  par(mar=c(0.2,2.2,5.2,2.2))
  plot(r, col=SOL_COLS, breaks=c(-0.5,0.5,1.5,2.5), legend=FALSE, main="", axes=TRUE)
  title(main=title_text, cex.main=1.85, line=3.2)
  par(mar=c(0,0,0,0)); plot.new()
  legend(x=0.5, y=0.55, xjust=0.5, yjust=0.5,
         legend=c("Not selected", sprintf("Newly selected (%s)", fmt(n)),
                  "Protected Areas (locked in)"),
         fill=SOL_COLS, ncol=3, bty="n", cex=1.35, xpd=TRUE)
  dev.off()
  invisible(n)
}

# --- two solutions compared ---
# 0 neither, 1 A only, 2 B only, 3 both, 4 locked-in protected area.
# title_fmt takes one %.2f, the IoU over new selections.
plot_pair <- function(run_a, run_b, outfile, title_fmt, lab_a, lab_b) {
  a <- new_of(run_a); b <- new_of(run_b)
  n_a <- sum(a & !b); n_b <- sum(b & !a); n_both <- sum(a & b)
  iou <- n_both / sum(a | b)
  r <- as_rast(ifelse(pa == 1, 4, ifelse(a & b, 3, ifelse(b, 2, ifelse(a, 1, 0)))))
  png(here("outputs", outfile), width=2200, height=1250, res=150)
  layout(matrix(1:2, nrow=2), heights=c(6,1.15))
  par(mar=c(0.2,2.2,6.5,2.2))
  plot(r, col=PAIR_COLS, breaks=c(-0.5,0.5,1.5,2.5,3.5,4.5), legend=FALSE, main="", axes=TRUE)
  title(main=sprintf(title_fmt, iou), cex.main=1.75, line=4.2)
  par(mar=c(0,0,0,0)); plot.new()
  # legend() sizes every column to the WIDEST label, so a short entry beside a
  # long one leaves a large gap. Keep the four labels comparable in length and
  # put the explanation in the title instead.
  legend(x=0.5, y=0.66, xjust=0.5, yjust=0.5,
         legend=c("neither", sprintf("%s (%s)", lab_a, fmt(n_a)),
                  sprintf("%s (%s)", lab_b, fmt(n_b)), sprintf("both (%s)", fmt(n_both))),
         fill=PAIR_COLS[1:4], ncol=4, bty="n", cex=1.3, xpd=TRUE)
  legend(x=0.5, y=0.18, xjust=0.5, yjust=0.5,
         legend="Protected Areas (locked in)", fill="grey40", ncol=1, bty="n", cex=1.3, xpd=TRUE)
  dev.off()
  invisible(c(a_only=n_a, b_only=n_b, both=n_both, iou=iou))
}

FIGS <- list(

  present = function() {
    n <- plot_solution(RUN_PRESENT, "best_solution_present_seed232373165.png",
           sprintf("CAPTAIN 3 solution, present-day conditions (seed %d)", SEED))
    L <- readLines(here(RUN_PRESENT,"training_log.tsv")); f <- strsplit(L[length(L)],"\t")[[1]]
    cat(sprintf("seed %d | reward %.3f | %s new cells | LC %s\n", SEED, as.numeric(f[3]), fmt(n),
        sub(".*threat_0\\s+","",grep("threat_0",readLines(here(RUN_PRESENT,"summary.txt")),value=TRUE))))
  },

  future = function() {
    n <- plot_solution(RUN_FUTURE, "future_solution_seed232373165.png",
           sprintf("CAPTAIN 3 solution with future disturbance (seed %d)", SEED))
    cat(sprintf("future solution | %s new cells\n", fmt(n)))
  },

  future_vs_static = function() {
    r <- plot_pair(RUN_PRESENT, RUN_FUTURE, "future_vs_static_seed232373165.png",
           sprintf("Static vs future-disturbance solution, seed %d (IoU = %%.2f)", SEED),
           "dropped", "added")
    cat(sprintf("static only %d | future only %d | both %d | IoU %.3f\n", r[1], r[2], r[3], r[4]))
  },

  staged = function() {
    r <- plot_pair(RUN_FUTURE, RUN_STAGED, "staged_vs_oneshot_seed232373165.png",
           sprintf("One-shot vs staged protection under future disturbance, seed %d (IoU = %%.2f)", SEED),
           "one-shot only", "staged only")
    cat(sprintf("one-shot only %d | staged only %d | both %d | IoU %.3f\n", r[1], r[2], r[3], r[4]))
  },

  # 82% of the study area is unchanged, so a diverging ramp centred on zero
  # shows the change far better than plotting the two layers side by side.
  disturbance_delta = function() {
    p <- mask(rast(file.path(DD,"disturbance.tif")), ugf, touches=TRUE)
    f <- mask(rast(file.path(DD,"disturbance_future_2075.tif")), ugf, touches=TRUE)
    d <- f - p
    v <- values(d); v <- v[!is.na(v)]
    same <- 100 * mean(abs(v) < 0.001)
    up   <- 100 * mean(v >=  0.001)
    down <- 100 * mean(v <= -0.001)
    brk  <- c(-0.9,-0.3,-0.1,-0.001,0.001,0.1,0.3,0.5,0.9)
    cols <- c("#2166ac","#67a9cf","#d1e5f0","grey92","#fddbc7","#ef8a62","#d6604d","#b2182b")
    png(here("outputs","disturbance_delta_2075.png"), width=2200, height=1150, res=150)
    layout(matrix(1:2, nrow=2), heights=c(6,0.9))
    par(mar=c(0.2,2.2,5.2,2.2))
    plot(d, col=cols, breaks=brk, legend=FALSE, main="", axes=TRUE)
    title(main="Change in disturbance, present -> 2075 (SSP4-RCP3.4)", cex.main=1.85, line=3.2)
    par(mar=c(0,0,0,0)); plot.new()
    legend(x=0.5, y=0.62, xjust=0.5, yjust=0.5, ncol=8, bty="n", cex=1.15, xpd=TRUE,
           legend=c("<= -0.3","-0.3 to -0.1","-0.1 to 0","no change","0 to 0.1",
                    "0.1 to 0.3","0.3 to 0.5","> 0.5"), fill=cols)
    mtext(sprintf("%.0f%% of cells unchanged  |  %.1f%% more disturbed, %.1f%% less  |  mean %+.3f",
                  same, up, down, mean(v)), side=1, line=-1.4, cex=1.15)
    dev.off()
    cat(sprintf("disturbance delta | %.0f%% unchanged | %.1f%% up | %.1f%% down | mean %+.3f\n",
        same, up, down, mean(v)))
  }
)

sel <- commandArgs(trailingOnly = TRUE)
if (!length(sel)) sel <- names(FIGS)
unknown <- setdiff(sel, names(FIGS))
if (length(unknown)) stop("unknown figure(s): ", paste(unknown, collapse=", "),
                          "\navailable: ", paste(names(FIGS), collapse=", "))
for (nm in sel) { cat("---", nm, "\n"); FIGS[[nm]]() }
