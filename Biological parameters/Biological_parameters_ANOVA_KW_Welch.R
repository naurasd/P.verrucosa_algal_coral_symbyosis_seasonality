############################################################
# BIOLOGICAL PARAMETERS — automatic transformations + test choice
# - Reads: Biological_Parameters.xlsx
# - First column: Season
# - Other columns: biological parameters (names exactly as in Excel)
#
# Pipeline per variable:
# 1) Try transformations (raw / sqrt / log(x+1) / log(-x+1) / logit if %)
# 2) Check Shapiro–Wilk normality per Season (on transformed values)
# 3) If all seasons normal:
#       Levene’s test
#       - if variances unequal: Welch ANOVA + Games–Howell
#       - else: classical ANOVA + Tukey
#    Else:
#       Kruskal–Wallis + Dunn (Bonferroni)
############################################################

############################################################
# PACKAGES
############################################################
library(readxl)
library(dunn.test)
library(rstatix)   # games_howell_test
library(dplyr)
library(car)

############################################################
# 0. READ DATA
############################################################

# set working directory to the directory the script is saved in

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

bio <- read_xlsx("Biological Data for ANOVA  Kruskal-Wallis Welch.xlsx", na = "NA")

# Ensure Season exists and has correct order
if (!"Season" %in% colnames(bio)) stop("Column 'Season' not found.")
bio <- bio %>% filter(!is.na(Season))
bio$Season <- factor(bio$Season, levels = c("Spring", "Summer", "Fall", "Winter"))

# Biological variable names (everything except Season)
bio_vars <- setdiff(colnames(bio), "Season")

cat("\nBiological variables detected:\n")
print(bio_vars)

############################################################
# 1. HELPERS
############################################################

# Shapiro per group
check_norm <- function(x, group) {
  out <- tapply(x, group, function(z) {
    z <- z[!is.na(z)]
    if (length(z) < 3) return(NA_real_)
    shapiro.test(z)$p.value
  })
  unlist(out)
}

# Levene
levene_p <- function(df, var) {
  lev <- leveneTest(df[[var]] ~ df$Season)
  as.numeric(lev[1, "Pr(>F)"])
}

# Welch ANOVA
welch_anova_test <- function(df, var, group_var = "Season") {
  form <- as.formula(paste0("`", var, "` ~ ", group_var))
  oneway.test(form, data = df, var.equal = FALSE)
}

# Games–Howell post-hoc
games_howell_test <- function(df, var, group_var = "Season") {
  form <- as.formula(paste0("`", var, "` ~ ", group_var))
  df %>% rstatix::games_howell_test(formula = form)
}

# Robust “is this % data?” detector:
# - either name suggests percent
# - or values look like 0–100 (or 0–1) with no negatives
looks_like_percent <- function(x, var_name) {
  z <- x[!is.na(x)]
  if (length(z) < 3) return(FALSE)
  if (grepl("%|percent|percentage", var_name, ignore.case = TRUE)) return(TRUE)
  if (min(z) >= 0 && max(z) <= 1) return(TRUE)
  if (min(z) >= 0 && max(z) <= 100) return(TRUE)
  FALSE
}

# Logit transform helper (expects proportions 0..1, avoids 0/1 exactly)
logit_transform <- function(p) {
  eps <- 1e-6
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}

# Create transformed vector given a method
apply_transform <- function(x, method, var_name) {
  if (method == "raw") return(x)
  
  if (method == "sqrt") {
    if (any(x < 0, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(sqrt(x))
  }
  
  if (method == "log(x+1)") {
    if (any(x <= -1, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(log(x + 1))
  }
  
  if (method == "log(-x+1)") {
    if (any(x >= 1, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(log(-x + 1))
  }
  
  if (method == "logit") {
    z <- x
    # if it looks like 0–100, convert to proportion
    if (max(z, na.rm = TRUE) > 1) z <- z / 100
    if (any(z < 0 | z > 1, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(logit_transform(z))
  }
  
  rep(NA_real_, length(x))
}

# Pick “best” transformation:
# - Try a set; choose the one that maximizes the minimum Shapiro p across seasons
# - If any transformation yields ALL seasons p > 0.05, take the one with highest min-p
choose_best_transform <- function(df, var) {
  x <- df[[var]]
  
  methods <- c("raw", "sqrt", "log(x+1)", "log(-x+1)")
  if (looks_like_percent(x, var)) methods <- c("logit", methods)
  
  scores <- data.frame(
    method = methods,
    min_p = NA_real_,
    all_normal = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(methods)) {
    m <- methods[i]
    xt <- apply_transform(x, m, var)
    if (all(is.na(xt))) {
      scores$min_p[i] <- NA_real_
      scores$all_normal[i] <- FALSE
      next
    }
    pvec <- check_norm(xt, df$Season)
    scores$min_p[i] <- suppressWarnings(min(pvec, na.rm = TRUE))
    scores$all_normal[i] <- all(pvec > 0.05, na.rm = TRUE)
  }
  
  # Prefer those that make all seasons normal; otherwise take highest min_p
  if (any(scores$all_normal, na.rm = TRUE)) {
    cand <- scores %>% filter(all_normal == TRUE)
    best <- cand$method[which.max(cand$min_p)]
  } else {
    best <- scores$method[which.max(scores$min_p)]
  }
  
  list(best_method = best, score_table = scores)
}

# Pretty printers (so you can copy stats into text)
print_anova_stats <- function(aov_obj, var_name) {
  s <- summary(aov_obj)[[1]]
  F_value <- s["Season", "F value"]
  p_value <- s["Season", "Pr(>F)"]
  df1 <- s["Season", "Df"]
  df2 <- s["Residuals", "Df"]
  cat("\nANOVA (classical) for ", var_name, ": F(", df1, ", ", df2, ") = ",
      round(F_value, 3), ", p = ", signif(p_value, 3), "\n", sep = "")
}

print_welch_stats <- function(welch_obj, var_name) {
  cat("\nWelch ANOVA for ", var_name, ": F(", round(welch_obj$parameter[1], 3), ", ",
      round(welch_obj$parameter[2], 3), ") = ", round(welch_obj$statistic, 3),
      ", p = ", signif(welch_obj$p.value, 3), "\n", sep = "")
}

print_kw_stats <- function(kw_obj, var_name) {
  H_val <- unname(kw_obj$statistic)
  dfv <- unname(kw_obj$parameter)
  cat("\nKruskal–Wallis for ", var_name, ": H = ", round(H_val, 3),
      ", df = ", dfv, ", p = ", signif(kw_obj$p.value, 3), "\n", sep = "")
}

############################################################
# 2. MAIN ANALYSIS FUNCTION (BIO)
############################################################
analyze_bio <- function(df, var) {
  
  cat("\n==============================\n")
  cat("ANALYSIS FOR", var, "\n")
  cat("==============================\n")
  
  if (!var %in% colnames(df)) stop("Column ", var, " not found.")
  
  # Work on a clean subset
  d <- df %>% filter(!is.na(Season), !is.na(.data[[var]]))
  
  # Ensure numeric (Excel sometimes reads weirdly)
  d[[var]] <- suppressWarnings(as.numeric(d[[var]]))
  
  # Choose best transformation automatically
  tr <- choose_best_transform(d, var)
  best_method <- tr$best_method
  
  d[[paste0(var, "_tr")]] <- apply_transform(d[[var]], best_method, var)
  
  # Drop NAs created by transform (if any)
  d <- d %>% filter(!is.na(.data[[paste0(var, "_tr")]]))
  
  cat("\nChosen transformation: ", best_method, "\n", sep = "")
  cat("\nTransformation candidates (min Shapiro p across seasons):\n")
  print(tr$score_table)
  
  # Normality per season on transformed data
  norm_p <- check_norm(d[[paste0(var, "_tr")]], d$Season)
  cat("\nShapiro–Wilk p-values per season (TRANSFORMED):\n")
  print(norm_p)
  
  normal_all <- all(norm_p > 0.05, na.rm = TRUE)
  
  if (normal_all) {
    # Levene on transformed data
    lev <- leveneTest(d[[paste0(var, "_tr")]] ~ d$Season)
    cat("\nLevene's test (TRANSFORMED):\n")
    print(lev)
    
    if (as.numeric(lev[1, "Pr(>F)"]) < 0.05) {
      cat("\nDecision: normal but heteroscedastic → Welch ANOVA + Games–Howell\n")
      
      # Welch ANOVA
      # (We use a temporary column name without spaces for formula safety)
      tmp_name <- "Y_tmp"
      d[[tmp_name]] <- d[[paste0(var, "_tr")]]
      
      wel <- oneway.test(Y_tmp ~ Season, data = d, var.equal = FALSE)
      print_welch_stats(wel, var)
      
      cat("\nGames–Howell post-hoc (on transformed values):\n")
      gh <- d %>% rstatix::games_howell_test(Y_tmp ~ Season)
      print(gh)
      
      return(invisible(list(
        var = var,
        transform = best_method,
        test = "Welch ANOVA + Games–Howell",
        norm_p = norm_p,
        levene = lev,
        welch = wel,
        posthoc = gh
      )))
      
    } else {
      cat("\nDecision: normal + homoscedastic → classical ANOVA + Tukey\n")
      
      tmp_name <- "Y_tmp"
      d[[tmp_name]] <- d[[paste0(var, "_tr")]]
      
      aov_obj <- aov(Y_tmp ~ Season, data = d)
      print_anova_stats(aov_obj, var)
      
      cat("\nTukey post-hoc (on transformed values):\n")
      tk <- TukeyHSD(aov_obj)
      print(tk)
      
      return(invisible(list(
        var = var,
        transform = best_method,
        test = "Classical ANOVA + Tukey",
        norm_p = norm_p,
        levene = lev,
        aov = aov_obj,
        posthoc = tk
      )))
    }
    
  } else {
    cat("\nDecision: not normal → Kruskal–Wallis + Dunn (Bonferroni)\n")
    
    # Kruskal on transformed values (still OK; KW uses ranks)
    tmp_name <- "Y_tmp"
    d[[tmp_name]] <- d[[paste0(var, "_tr")]]
    
    kw <- kruskal.test(Y_tmp ~ Season, data = d)
    print_kw_stats(kw, var)
    
    cat("\nDunn post-hoc (Bonferroni) (on transformed values):\n")
    dn <- dunn.test(d[[tmp_name]], d$Season, method = "bonferroni")
    print(dn)
    
    return(invisible(list(
      var = var,
      transform = best_method,
      test = "Kruskal–Wallis + Dunn (Bonferroni)",
      norm_p = norm_p,
      kw = kw,
      posthoc = dn
    )))
  }
}

############################################################
# 3. RUN FOR EACH BIOLOGICAL VARIABLE
############################################################
cat("\n=========== BIOLOGICAL VARIABLES ===========\n")

# Run all
bio_results <- list()
for (v in bio_vars) {
  bio_results[[v]] <- analyze_bio(bio, v)
}

############################################################
# 4. QUICK SUMMARY (so you can copy the ‘# …’ comments style)
############################################################
cat("\n=========== SUMMARY (copy/paste helper) ===========\n")
for (v in bio_vars) {
  res <- bio_results[[v]]
  if (is.null(res)) next
  cat('analyze_bio(bio, "', v, '")   # ', res$test,
      " (transform: ", res$transform, ")\n", sep = "")
}

###############
analyze_bio(bio, "C host tissue")             # Classical ANOVA + Tukey (transform: log(x+1))
analyze_bio(bio, "N host tissue")             # Welch ANOVA + Games–Howell (transform: logit)
analyze_bio(bio, "C:N ratio host tissue")     # Welch ANOVA + Games–Howell (transform: raw)
analyze_bio(bio, "C algal symbiont")          # Kruskal–Wallis + Dunn (Bonferroni) (transform: raw)
analyze_bio(bio, "N algal symbiont")          # Kruskal–Wallis + Dunn (Bonferroni) (transform: raw)
analyze_bio(bio, "C:N ratio algal symbiont")  # Kruskal–Wallis + Dunn (Bonferroni) (transform: raw)
analyze_bio(bio, "δ13C host tissue")          # Classical ANOVA + Tukey (transform: raw)
analyze_bio(bio, "δ15N host tissue")          # Classical ANOVA + Tukey (transform: logit)
analyze_bio(bio, "δ13C algal symbiont")       # Kruskal–Wallis + Dunn (Bonferroni) (transform: raw)
analyze_bio(bio, "δ15N algal symbiont")       # Kruskal–Wallis + Dunn (Bonferroni) (transform: logit)
analyze_bio(bio, "δ15N host-algal symbiont")  # Classical ANOVA + Tukey (transform: log(x+1))
analyze_bio(bio, "δ13C host-algal symbiont")  # Welch ANOVA + Games–Howell (transform: raw)
analyze_bio(bio, "Algal symbiont density")    # Classical ANOVA + Tukey (transform: raw)
analyze_bio(bio, "Mitotic index")             # Welch ANOVA + Games–Howell (transform: log(-x+1))
analyze_bio(bio, "Chlorophyll a")             # Classical ANOVA + Tukey (transform: raw)
analyze_bio(bio, "nifH (Ct values)")          # Classical ANOVA + Tukey (transform: raw)

######Table

## ============================================================
## SUPPLEMENT TABLE (BIO): SS, DF, test statistic, p-value
## - Uses your existing `bio_results` list
## - SS is only defined for classical ANOVA; Welch/KW -> NA
## ============================================================

bio_supp_table <- bind_rows(lapply(names(bio_results), function(param) {
  
  res <- bio_results[[param]]
  
  # Default fields
  transform_used <- res$transform
  test_used      <- res$test
  
  # ---- Classical ANOVA ----
  if (grepl("Classical ANOVA", test_used, fixed = TRUE)) {
    
    s <- summary(res$aov)[[1]]
    
    SS_season <- unname(s["Season", "Sum Sq"])
    DF1       <- unname(s["Season", "Df"])
    DF2       <- unname(s["Residuals", "Df"])
    stat      <- unname(s["Season", "F value"])
    pval      <- unname(s["Season", "Pr(>F)"])
    
    return(data.frame(
      Parameter      = param,
      Model          = "Classical one-way ANOVA",
      Transformation = transform_used,
      SS             = SS_season,
      DF1            = DF1,
      DF2            = DF2,
      Test_statistic = stat,
      Statistic_type = "F",
      p_value        = pval,
      stringsAsFactors = FALSE
    ))
  }
  
  # ---- Welch ANOVA ----
  if (grepl("Welch ANOVA", test_used, fixed = TRUE)) {
    
    wel <- res$welch
    
    DF1  <- unname(wel$parameter[1])  # numerator df
    DF2  <- unname(wel$parameter[2])  # denominator df
    stat <- unname(wel$statistic)
    pval <- unname(wel$p.value)
    
    return(data.frame(
      Parameter      = param,
      Model          = "Welch one-way ANOVA",
      Transformation = transform_used,
      SS             = NA_real_,
      DF1            = DF1,
      DF2            = DF2,
      Test_statistic = stat,
      Statistic_type = "F",
      p_value        = pval,
      stringsAsFactors = FALSE
    ))
  }
  
  # ---- Kruskal–Wallis ----
  kw <- res$kw
  
  return(data.frame(
    Parameter      = param,
    Model          = "Kruskal–Wallis",
    Transformation = transform_used,
    SS             = NA_real_,
    DF1            = unname(kw$parameter),  # df
    DF2            = NA_real_,
    Test_statistic = unname(kw$statistic),
    Statistic_type = "H",
    p_value        = unname(kw$p.value),
    stringsAsFactors = FALSE
  ))
})) %>%
  # Pretty formatting for df and p (optional, but Word-friendly)
  mutate(
    DF1 = ifelse(is.na(DF1), NA, round(DF1, 3)),
    DF2 = ifelse(is.na(DF2), NA, round(DF2, 3)),
    SS  = ifelse(is.na(SS), NA, round(SS, 5)),
    Test_statistic = round(Test_statistic, 3),
    p_value = signif(p_value, 3)
  )

## 1) Print as tab-delimited text
write.table(bio_supp_table,sep = "\t", row.names = FALSE, quote = FALSE)

## 2) Or print as a markdown table 
print(bio_supp_table)
