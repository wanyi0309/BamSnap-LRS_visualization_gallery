library(scatterplot3d)

# Read Highlight data
dat <- read.delim(
  "Highlight_3D_plot.txt",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)
# Convert to numeric
dat$Plot_interval_bp <- as.numeric(
  gsub(",", "", dat[["Plot interval (bp)"]])
)
dat$Displayed_records <- as.numeric(
  gsub(",", "", dat[["total_reads"]])
)
dat$Runtime_s <- as.numeric(dat$Runtime_s)
dat$Maximum_RSS_GiB <- as.numeric(
  dat$Maximum_RSS_GiB
)
# Check data
sum(is.na(dat$Plot_interval_bp))
sum(is.na(dat$Displayed_records))
sum(is.na(dat$Runtime_s))
sum(is.na(dat$Maximum_RSS_GiB))


# Plot A: Runtime
pdf(
  "Highlight_Runtime_3D.pdf",
  width = 8,
  height = 7
)
s3d <- scatterplot3d(
  x = dat$Plot_interval_bp,
  y = dat$Displayed_records,
  z = dat$Runtime_s,
  type = "h",
  pch = 19,
  color = "grey65",
  angle = 45,
  grid = TRUE,
  box = TRUE,
  xlab = "Plot interval (bp)",
  ylab = "Displayed alignment records",
  zlab = "Runtime (s)",
  main = "Highlight: Runtime",
  cex.axis = 0.9,
  cex.lab = 1.1
)

# Add solid points
s3d$points3d(
  dat$Plot_interval_bp,
  dat$Displayed_records,
  dat$Runtime_s,
  pch = 19,
  col = "#2878B5",
  cex = 1
)
dev.off()

# Plot B: Memory
pdf(
  "Highlight_Memory_3D.pdf",
  width = 8,
  height = 7
)
G3d <- scatterplot3d(
  x = dat$Plot_interval_bp,
  y = dat$Displayed_records,
  z = dat$Maximum_RSS_GiB,
  type = "h",
  pch = 19,
  color = "grey65",
  angle = 45,
  grid = TRUE,
  box = TRUE,
  xlab = "Plot interval (bp)",
  ylab = "Displayed alignment records",
  zlab = "Maximum RSS (GiB)",
  main = "Highlight: Memory usage",
  cex.axis = 0.9,
  cex.lab = 1.1
)

# Add solid points
G3d$points3d(
  dat$Plot_interval_bp,
  dat$Displayed_records,
  dat$Maximum_RSS_GiB,
  pch = 19,
  col = "#D97932",
  cex = 1
)
dev.off()
