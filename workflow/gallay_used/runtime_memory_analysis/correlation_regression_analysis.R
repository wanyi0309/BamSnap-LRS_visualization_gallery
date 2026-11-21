
# ==========================================
# 1. Read data
# ==========================================

dup <- read.delim(
  "DUP_3D_plot.txt",
  sep = "\t",
  check.names = FALSE
)

highlight <- read.delim(
  "Highlight_3D_plot.txt",
  sep = "\t",
  check.names = FALSE
)


# ==========================================
# 2. Prepare DUP data
# ==========================================

dup$Region_size <- as.numeric(
  gsub(",", "", dup[["Plot interval (bp)"]])
)

dup$Records <- as.numeric(
  gsub(",", "", dup[["all_displayed_records"]])
)

dup$Runtime <- as.numeric(dup$Runtime_s)

dup$Memory <- as.numeric(dup$Maximum_RSS_GiB)


# ==========================================
# 3. Prepare Highlight data
# ==========================================

highlight$Region_size <- as.numeric(
  gsub(",", "", highlight[["Plot interval (bp)"]])
)

highlight$Records <- as.numeric(
  gsub(",", "", highlight[["total_reads"]])
)

highlight$Runtime <- as.numeric(highlight$Runtime_s)

highlight$Memory <- as.numeric(
  highlight$Maximum_RSS_GiB
)


# ==========================================
# 4. Check data
# ==========================================

stopifnot(
  nrow(dup) == 50,
  nrow(highlight) == 50,
  all(complete.cases(dup[, c(
    "Region_size", "Records", "Runtime", "Memory"
  )])),
  all(complete.cases(highlight[, c(
    "Region_size", "Records", "Runtime", "Memory"
  )])),
  all(unlist(dup[, c(
    "Region_size", "Records", "Runtime", "Memory"
  )]) > 0),
  all(unlist(highlight[, c(
    "Region_size", "Records", "Runtime", "Memory"
  )]) > 0)
)


# ==========================================
# 5. Multiple linear regression
# ==========================================

# DUP
dup_runtime <- lm(
  log(Runtime) ~ log(Region_size) + log(Records),
  data = dup
)

dup_memory <- lm(
  log(Memory) ~ log(Region_size) + log(Records),
  data = dup
)


# Highlight
highlight_runtime <- lm(
  log(Runtime) ~ log(Region_size) + log(Records),
  data = highlight
)

highlight_memory <- lm(
  log(Memory) ~ log(Region_size) + log(Records),
  data = highlight
)


# ==========================================
# 6. Display results
# ==========================================

cat("\n========== DUP: Runtime ==========\n")
print(summary(dup_runtime))

cat("\n========== DUP: Memory ==========\n")
print(summary(dup_memory))

cat("\n========== Highlight: Runtime ==========\n")
print(summary(highlight_runtime))

cat("\n========== Highlight: Memory ==========\n")
print(summary(highlight_memory))


# ==========================================
# 7. Export regression results
# ==========================================

get_results <- function(model, mode, outcome) {
  
  s <- summary(model)
  ci <- confint(model)
  
  data.frame(
    Mode = mode,
    Outcome = outcome,
    Variable = rownames(coef(s)),
    Coefficient = coef(s)[, 1],
    Std_Error = coef(s)[, 2],
    t_value = coef(s)[, 3],
    P_value = coef(s)[, 4],
    CI_lower = ci[, 1],
    CI_upper = ci[, 2],
    Adjusted_R2 = s$adj.r.squared
  )
}

results <- rbind(
  get_results(dup_runtime, "DUP", "Runtime"),
  get_results(dup_memory, "DUP", "Memory"),
  get_results(highlight_runtime, "Highlight", "Runtime"),
  get_results(highlight_memory, "Highlight", "Memory")
)

#write.table(
  results,
  "regression_results.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\nRegression results saved!\n")



# Correlation between two predictors
r <- cor(
  log(highlight$Region_size),
  log(highlight$Records)
)

# Calculate VIF
VIF <- 1 / (1 - r^2)

cat("Correlation:", r, "\n")
cat("VIF:", VIF, "\n")


library(car)
vif(highlight_runtime)
vif(highlight_memory)


# Region size vs Memory
cor.test(
  log(highlight$Region_size),
  log(highlight$Memory),
  method = "pearson"
)

# Records vs Memory
cor.test(
  log(highlight$Records),
  log(highlight$Memory),
  method = "pearson"
)
