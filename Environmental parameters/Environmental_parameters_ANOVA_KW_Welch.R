############################################################
# ENVIRONMENTAL PARAMETERS — automatic transformations + test choice
#
# Pipeline per variable:
# 1) Try transformations (raw / sqrt / log(x+1) / log(-x+1))
# 2) Check Shapiro–Wilk normality per Season (on transformed values)
#    - If a Season has > 5000 observations, a random subset of 5000 is used
#      (Shapiro–Wilk limit in R).
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
library(car)
library(rstatix)   # games_howell_test
library(dplyr)

############################################################
# 0. READ DATA
############################################################

# set working directory to the directory the script is saved in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ---- MONTHLY RECORDED PARAMETERS ----
monthly <- read_xlsx("Environmental Data file Tilstra et al.xlsx",
                     sheet = "all_monthly", na = "NA") %>%
  na.omit()
# This omits all rows where Season is NA (measurements that do not belong to coral sampling seasons)

# Ensure Season exists and has correct order
if (!"Season" %in% colnames(monthly)) stop("Column 'Season' not found in monthly sheet.")
monthly <- monthly %>% filter(!is.na(Season))
monthly$Season <- factor(monthly$Season, levels = c("Spring", "Summer", "Fall", "Winter"))

# remove date column

monthly<-monthly[-2]

# Simplify some parameter names

colnames(monthly)[4]  <- "PAR"
colnames(monthly)[5]  <- "Nitrate"
colnames(monthly)[6]  <- "Nitrite"
colnames(monthly)[7]  <- "Ammonia"
colnames(monthly)[8]  <- "DIN"
colnames(monthly)[9]  <- "Phosphate"
colnames(monthly)[10] <- "DIN_DIP"

monthly_vars <- setdiff(colnames(monthly), "Season")

cat("\nMonthly variables detected:\n")
print(monthly_vars)

# ---- TEMPERATURE ----
temp <- read_xlsx("Environmental Data file Tilstra et al.xlsx", sheet = "Temperature",
                  col_types = c("text", "text", "numeric"),
                  na = "NA") %>%
  na.omit()

# Ensure Season exists and has correct order
if (!"Season" %in% colnames(temp)) stop("Column 'Season' not found in Temperature sheet.")
temp <- temp %>% filter(!is.na(Season))
temp$Season <- factor(temp$Season, levels = c("Spring", "Summer", "Fall", "Winter"))

# remove date column

temp<-temp[-2]

# ---- BIMONTHLY RECORDED PARAMETERS ----
bimonthly <- read_xlsx("Environmental Data file Tilstra et al.xlsx",
                     sheet = "all_bimonthly", na = "NA") %>%
  na.omit()

# Ensure Season exists and has correct order
if (!"Season" %in% colnames(bimonthly)) stop("Column 'Season' not found in bimonthly sheet.")
bimonthly <- bimonthly %>% filter(!is.na(Season))
bimonthly$Season <- factor(bimonthly$Season, levels = c("Spring", "Summer", "Fall", "Winter"))

# remove date column

bimonthly<-bimonthly[-2]

# Simplify parameter names

colnames(bimonthly)[3]  <- "DO"
colnames(bimonthly)[4]  <- "DOC"
colnames(bimonthly)[6]  <- "DON"
colnames(bimonthly)[7]  <- "DOC_DON"

# Remove TDN parameter which will not be further analysed
bimonthly<-bimonthly %>% select(-5)

bimonthly_vars <- setdiff(colnames(bimonthly), "Season")

cat("\nBimonthly variables detected:\n")
print(bimonthly_vars)

############################################################
# 1. HELPERS
############################################################

# Shapiro per group
# - Uses random subset of 5000 if a season has > 5000 values (Shapiro limit)
check_norm <- function(x, group) {
  out <- tapply(x, group, function(z) {
    z <- z[!is.na(z)]
    if (length(z) < 3) return(NA_real_)
    if (length(z) > 5000) {
      set.seed(1)
      z <- sample(z, 5000)
    }
    shapiro.test(z)$p.value
  })
  unlist(out)
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

# Robust “is this % data?” detector
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
    if (max(z, na.rm = TRUE) > 1) z <- z / 100
    if (any(z < 0 | z > 1, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(logit_transform(z))
  }
  
  rep(NA_real_, length(x))
}

# Pick “best” transformation:
# - Choose the one that maximizes the minimum Shapiro p across seasons
# - Prefer transforms that make ALL seasons p > 0.05
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
  
  if (any(scores$all_normal, na.rm = TRUE)) {
    cand <- scores %>% filter(all_normal == TRUE)
    best <- cand$method[which.max(cand$min_p)]
  } else {
    best <- scores$method[which.max(scores$min_p)]
  }
  
  list(best_method = best, score_table = scores)
}

# Pretty printers (copy-friendly)
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
# 2. MAIN ANALYSIS FUNCTION (ENV) — mirrors analyze_bio()
############################################################
analyze_env <- function(df, var) {
  
  cat("\n==============================\n")
  cat("ANALYSIS FOR", var, "\n")
  cat("==============================\n")
  
  if (!var %in% colnames(df)) stop("Column ", var, " not found.")
  
  d <- df %>% filter(!is.na(Season), !is.na(.data[[var]]))
  d[[var]] <- suppressWarnings(as.numeric(d[[var]]))
  
  tr <- choose_best_transform(d, var)
  best_method <- tr$best_method
  
  d[[paste0(var, "_tr")]] <- apply_transform(d[[var]], best_method, var)
  d <- d %>% filter(!is.na(.data[[paste0(var, "_tr")]]))
  
  cat("\nChosen transformation: ", best_method, "\n", sep = "")
  cat("\nTransformation candidates (min Shapiro p across seasons):\n")
  print(tr$score_table)
  
  norm_p <- check_norm(d[[paste0(var, "_tr")]], d$Season)
  cat("\nShapiro–Wilk p-values per season (TRANSFORMED):\n")
  print(norm_p)
  
  normal_all <- all(norm_p > 0.05, na.rm = TRUE)
  
  if (normal_all) {
    lev <- leveneTest(d[[paste0(var, "_tr")]] ~ d$Season)
    cat("\nLevene's test (TRANSFORMED):\n")
    print(lev)
    
    if (as.numeric(lev[1, "Pr(>F)"]) < 0.05) {
      cat("\nDecision: normal but heteroscedastic → Welch ANOVA + Games–Howell\n")
      
      d$Y_tmp <- d[[paste0(var, "_tr")]]
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
      
      d$Y_tmp <- d[[paste0(var, "_tr")]]
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
    
    d$Y_tmp <- d[[paste0(var, "_tr")]]
    kw <- kruskal.test(Y_tmp ~ Season, data = d)
    print_kw_stats(kw, var)
    
    cat("\nDunn post-hoc (Bonferroni) (on transformed values):\n")
    dn <- dunn.test(d$Y_tmp, d$Season, method = "bonferroni")
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
# 3. RUN FOR EACH ENV VARIABLE
############################################################
cat("\n=========== MONTHLY VARIABLES ===========\n")
monthly_results <- list()
for (v in monthly_vars) monthly_results[[v]] <- analyze_env(monthly, v)

cat("\n=========== TEMPERATURE ===========\n")
temp_results <- list()
temp_results[["Temperature"]] <- analyze_env(temp, "Temperature")

cat("\n=========== BIMONTHLY VARIABLES ===========\n")
bimonthly_results <- list()
for (v in bimonthly_vars) bimonthly_results[[v]] <- analyze_env(bimonthly, v)

############################################################
# 4. QUICK SUMMARY (copy/paste helper)
############################################################
cat("\n=========== SUMMARY (copy/paste helper) ===========\n")

cat("\n# MONTHLY\n")
for (v in monthly_vars) {
  res <- monthly_results[[v]]
  if (!is.null(res)) cat('analyze_env(monthly, "', v, '")   # ', res$test,
                         " (transform: ", res$transform, ")\n", sep = "")
}

cat("\n# TEMPERATURE\n")
resT <- temp_results[["Temperature"]]
cat('analyze_env(temp, "Temperature")   # ', resT$test,
    " (transform: ", resT$transform, ")\n", sep = "")

cat("\n# BIMONTHLY\n")
for (v in bimonthly_vars) {
  res <- bimonthly_results[[v]]
  if (!is.null(res)) cat('analyze_env(bimonthly, "', v, '")   # ', res$test,
                         " (transform: ", res$transform, ")\n", sep = "")
}

############################################################
# 5. OPTIONAL: SUPPLEMENT TABLES 
############################################################
make_supp_table <- function(results_list) {
  bind_rows(lapply(names(results_list), function(param) {
    
    res <- results_list[[param]]
    if (is.null(res)) return(NULL)
    
    transform_used <- res$transform
    test_used      <- res$test
    
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
    
    if (grepl("Welch ANOVA", test_used, fixed = TRUE)) {
      wel <- res$welch
      DF1  <- unname(wel$parameter[1])
      DF2  <- unname(wel$parameter[2])
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
    
    kw <- res$kw
    data.frame(
      Parameter      = param,
      Model          = "Kruskal–Wallis",
      Transformation = transform_used,
      SS             = NA_real_,
      DF1            = unname(kw$parameter),
      DF2            = NA_real_,
      Test_statistic = unname(kw$statistic),
      Statistic_type = "H",
      p_value        = unname(kw$p.value),
      stringsAsFactors = FALSE
    )
  })) %>%
    mutate(
      DF1 = ifelse(is.na(DF1), NA, round(DF1, 3)),
      DF2 = ifelse(is.na(DF2), NA, round(DF2, 3)),
      SS  = ifelse(is.na(SS), NA, round(SS, 5)),
      Test_statistic = round(Test_statistic, 3),
      p_value = signif(p_value, 3)
    )
}

monthly_supp_table <- make_supp_table(monthly_results)
temp_supp_table    <- make_supp_table(temp_results)
bimonthly_supp_table      <- make_supp_table(bimonthly_results)

cat("\n=========== SUPPLEMENT TABLE: MONTHLY ===========\n")
write.table(monthly_supp_table, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=========== SUPPLEMENT TABLE: TEMPERATURE ===========\n")
write.table(temp_supp_table, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=========== SUPPLEMENT TABLE: BIMONTHLY ===========\n")
write.table(bimonthly_supp_table, sep = "\t", row.names = FALSE, quote = FALSE)