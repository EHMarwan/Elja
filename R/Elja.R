
utils::globalVariables(c("odd_ratio", "coefficients"))

#' Linear regression for EnvWAS/EWAS analysis
#'
#' @description
#' A tool for Environment-Wide Association Studies (EnvWAS / EWAS) which are repeated analysis. This function is espacially for linear regressions and allows the addition of adjustment variables.
#'
#'
#' @param var A categorical and binary variable. It is generally your outcome.
#' @param var_adjust A vector containing the names of the fixed adjustment variables for all the models.
#' @param data A dataframe containing all the variables needed for the analysis.
#' @param manplot Generate a Manhattan plot of the results of the analysis.
#' @param nbvalmanplot The number of variables to include in each Manhattan plot.
#' @param Bonferroni Add a dashed bar to the Manhattan plot showing the Bonferroni significance level.
#' @param FDR Add a dashed bar to the Manhattan plot showing the False Discovery Rate (Benjamini-Hochberg method) significance threshold. NA if all p-values > FDR corrected p-values.
#' @param manplotsign Generates a Manhattan plot with only significant results (p<0.05).
#'
#' @import ggplot2 dplyr devtools
#' @rawNamespace import(stats, except = c("filter","lag"))
#' @rawNamespace import(MASS, except = "select")
#'
#' @return A Dataframe with results for each variable of the model.
#'
#' @export
#'
#' @examples
#' ### Creating a dataframe with random variables
#'
#'
#'exposure1 <- as.factor(sample(c("Always", "Often","Never"), size = 400, replace = TRUE))
#'exposure2 <- as.factor(sample(c("Always", "Often","Never"), size = 400, replace = TRUE))
#'exposure3 <- as.factor(sample(c("Exposed", "Not exposed"), size = 400, replace = TRUE))
#'exposure4 <- as.factor(sample(c("Yes", "No"), size = 400, replace = TRUE))
#'exposure5 <- as.factor(sample(c("Yes", "No"), size = 400, replace = TRUE))
#'exposure6 <- as.factor(sample(c("Yes", "No"), size = 400, replace = TRUE))
#'exposure7 <- as.numeric(sample(0:400, size = 400, rep = TRUE))
#'exposure8 <- as.numeric(sample(0:400, size = 400, rep = TRUE))
#'exposure9 <- as.numeric(sample(0:0.1, size = 400, rep = TRUE))
#'outcome <- as.numeric(sample(0:100, size = 400, rep = TRUE))
#'
#'data <- data.frame(exposure1,exposure2,exposure3,exposure4,exposure5,
#'                   exposure6,exposure7,exposure8,exposure9, outcome)
#'
#'
#'### Launch of the EnvWAS analysis
#'
#'ELJAlinear(var = 'outcome',data = data, var_adjust = c('exposure5'))
#'
#'@references
#'1. Dunn OJ. Multiple Comparisons Among Means. Journal of the American Statistical Association. 1961;56(293):52‑64.
#'2. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society: Series B (Methodological). 1995;57(1):289‑300.
#'
ELJAlinear <- function(var, var_adjust = NULL, data,
                      manplot = TRUE, nbvalmanplot = 100, Bonferroni = FALSE, FDR = FALSE, manplotsign = FALSE) {


  # Definition de y et x ( x = variable explicative dans la boucle + variable d'ajustements fixe s'il y en a)
  y <- data[, var]

  if (!is.null(var_adjust)) {
    x_adjust <- data[, var_adjust, drop = FALSE]
    x <- cbind(data[, -which(names(data) %in% c(var, var_adjust))], x_adjust)
    x_formula <- paste0(paste0(names(x_adjust), collapse = " + "), collapse = " + ")
    x_formula <- ifelse(length(x_formula) > 0, paste0(x_formula, " + "), "")
  } else {
    x <- data[, -which(names(data) == var)]
    x_formula <- ""
  }

  n <- length(x)

  # Analyses de regressions pour chaque variable explicative

  x_cols <- colnames(x)
  lm_results <- list()


  for (col in x_cols) {
    if(col %in% var_adjust) {
      next # passer à la variable explicative suivante si elle fait partie des variables d'ajustement
    }
    fit <- lm(as.formula(paste0(var, " ~ ", x_formula, col)), data = data)
    summary_fit <- summary(fit)
    coef_names <- rownames(summary_fit$coefficients)[-1] # Exclure l'intercept

    for (name in coef_names) {
      level <- gsub("x\\[.*\\]", " ", name) # Extraire la modalite de la variable explicative
      coeff <- summary_fit$coefficients[name, 1]
      ci <- suppressMessages(confint(fit, parm = name, level = 0.95))
      p_value <- summary_fit$coefficients[name, 4]
      n <- nobs(fit)
      aic <- AIC(fit)
      lm_results[[paste0(col, "_", level)]] <- data.frame(
        variable = col,
        level = level,
        coefficients = coeff,
        ci_low = ci[1],
        ci_high = ci[2],
        p_value = p_value,
        n = n,
        AIC = aic
      )
    }
  }

  # Assemblage des resultats dans un dataframe unique
  results <- do.call(rbind, lm_results)
  results$variable_level <- rownames(results)
  results$variable <- NULL



  # Suppressions des Betas associes aux variables d'ajustement
  if (!is.null(var_adjust)) {
    results <- results[!grepl(paste(var_adjust, collapse="|"), results$variable_level),]
    rownames(results) <- results$level
  }

  # Stockage du tableau des resultats
  results$variable_level <- NULL
  results <<- results

  nbvar <- as.numeric(length(results$p_value))

  # Calcul du FDR

  results_test <- results %>% arrange(p_value) %>% mutate(ligne = row_number())

  results_test$q <- (results_test$ligne / length(results_test$p_value)) * 0.05

  results_test$test <- results_test$p_value < results_test$q

  FDRcalc <- results_test %>% dplyr::filter(test == TRUE) %>% slice_tail(n = 1) %>% select(q) %>% as.numeric()

  # Creation du Manhattan Plot
  if (manplot == TRUE) {

    # Creer un nouveau tableau de resultats avec autant de multiples de nbvalmanplot
    n_new <- ceiling(nrow(results) / nbvalmanplot) * nbvalmanplot
    results_new <- results[1:n_new, ]

    # Ajout de valeurs manquantes pour completer la derniere sous-liste
    n_missing <- n_new - nrow(results_new)
    if (n_missing > 0) {
      results_new[(nrow(results_new)+1):(n_new),] <- NA
    }

    # Diviser le nouveau tableau de resultats en sous groupes de nbvalmanplot lignes
    results_list <- split(results_new, ceiling(seq_along(results_new$p_value)/nbvalmanplot))

    # Creation d'un manhattan plot pour chaque sous groupe de donnees
    for (i in seq_along(results_list)) {
      sub_results <- results_list[[i]]

      # Verification de si le sous groupe contient des valeurs manquantes
      if (any(is.na(sub_results$p_value))) {
        sub_results <- na.omit(sub_results)
        n_sub <- nrow(sub_results)
      } else {
        n_sub <- nbvalmanplot
      }

      # Creer le manhattan plot pour les sous groupes
      Manplot <- ggplot(sub_results, aes(x = level, y = -log10(p_value),
                                        color = ifelse(coefficients > 0, "> 0", "< 0"))) +
        geom_point() +
        scale_color_manual(values = c("> 0" = "red", "< 0" = "green")) +
        coord_flip() +
        labs(color = "coefficients") +
        guides(color = guide_legend(title = "Coefficients"))

      if (Bonferroni == TRUE & FDR == TRUE) {
        Manplot <- Manplot +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(0.05/nbvar), linetype = paste0( "Bonferroni : ",round(0.05/nbvar,5))), colour = "black",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(FDRcalc), linetype = paste0( "FDR : ",round(FDRcalc,5))), colour = "purple",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed","dashed","dashed"),
                                guide = guide_legend(override.aes = list(color = c("black","purple","#ffd800"))))
      } else if (Bonferroni == TRUE & FDR == FALSE) {
          Manplot <- Manplot +
            geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
            geom_hline(aes(yintercept = -log10(0.05/nbvar), linetype = paste0( "Bonferroni : ",round(0.05/nbvar,5))), colour = "black",linewidth = 0.9) +
            scale_linetype_manual(name = "Threshold",values = c("dashed","dashed"),
                                  guide = guide_legend(override.aes = list(color = c("black","#ffd800"))))
      } else if (Bonferroni == FALSE & FDR == TRUE) {
        Manplot <- Manplot +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(FDRcalc), linetype = paste0( "FDR : ",round(FDRcalc,5))), colour = "purple",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed","dashed"),
                                guide = guide_legend(override.aes = list(color = c("purple","#ffd800"))))
      } else {
        Manplot <- Manplot +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed"),
                                guide = guide_legend(override.aes = list(color = c("#ffd800"))))
      }

      Manplot <<- Manplot
      print(Manplot)
    }

  } else {
    cat("The manhattan plot is not shown.\n")
  }

  # Creation du tableau avec les valeurs significatives uniquement

  if (manplotsign == TRUE) {

    # Creer un nouveau tableau de resultats avec autant de multiples de nbvalmanplot
    resultssign <- results %>% dplyr::filter(p_value < 0.05)
    n_new_sign <- ceiling(nrow(resultssign) / nbvalmanplot) * nbvalmanplot
    resultssign_new <- resultssign[1:n_new_sign, ]

    # Ajouter des valeurs manquantes pour completer la derniere sous liste
    n_missing_sign <- n_new_sign - nrow(resultssign_new)
    if (n_missing_sign > 0) {
      resultssign_new[(nrow(resultssign_new)+1):(n_new_sign),] <- NA
    }

    # Diviser le nouveau tableau de resultats en sous groupes de 30 lignes
    results_list <- split(resultssign_new, ceiling(seq_along(resultssign_new$p_value)/nbvalmanplot))

    # Creer un manhattan plot pour chaque sous-groupe de donnees
    for (i in seq_along(results_list)) {
      sub_results <- results_list[[i]]

      # Verification de si le sous-groupe contient des valeurs manquantes
      if (any(is.na(sub_results$p_value))) {
        sub_results <- na.omit(sub_results)
        n_sub <- nrow(sub_results)
      } else {
        n_sub <- nbvalmanplot
      }

      # Creation d'un manhattan plot pour chaque sous-groupe de donnees
      Manplotsign <- ggplot(sub_results, aes(x = level, y = -log10(p_value),
                                            color = ifelse(coefficients > 0, "> 0", "< 0"))) +
        geom_point() +
        scale_color_manual(values = c("> 0" = "red", "< 0" = "green")) +
        coord_flip() +
        labs(color = "coefficients") +
        guides(color = guide_legend(title = "Coefficients"))

      if (Bonferroni == TRUE & FDR == TRUE) {
        Manplotsign <- Manplotsign +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(0.05/nbvar), linetype = paste0( "Bonferroni : ",round(0.05/nbvar,5))), colour = "black",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(FDRcalc), linetype = paste0( "FDR : ",round(FDRcalc,5))), colour = "purple",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed","dashed","dashed"),
                                guide = guide_legend(override.aes = list(color = c("black","purple","#ffd800"))))
      } else if (Bonferroni == TRUE & FDR == FALSE) {
        Manplotsign <- Manplotsign +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(0.05/nbvar), linetype = paste0( "Bonferroni : ",round(0.05/nbvar,5))), colour = "black",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed","dashed"),
                                guide = guide_legend(override.aes = list(color = c("black","#ffd800"))))
      } else if (Bonferroni == FALSE & FDR == TRUE) {
        Manplotsign <- Manplotsign +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(FDRcalc), linetype = paste0( "FDR : ",round(FDRcalc,5))), colour = "purple",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed","dashed"),
                                guide = guide_legend(override.aes = list(color = c("purple","#ffd800"))))
      } else {
        Manplotsign <- Manplotsign +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed"),
                                guide = guide_legend(override.aes = list(color = c("#ffd800"))))
      }

      Manplotsign <<- Manplotsign
      print(Manplotsign)
    }

  } else {
    cat("The manhattan plot (manplotsign) is not shown.\n")
  }
}



#' Logistic regression tool for EnvWAS/EWAS analysis
#'
#' @description
#' A tool for Environment-Wide Association Studies (EnvWAS / EWAS) which are repeated analysis. This function is espacially for logistic regression based on the glm function with a binomial family with a logit link and allows the addition of adjustment variables.
#'
#'
#' @param var A categorical and binary variable. It is generally your outcome.
#' @param var_adjust A vector containing the names of the fixed adjustment variables for all the models.
#' @param data A dataframe containing all the variables needed for the analysis.
#' @param manplot Generate a Manhattan plot of the results of the analysis.
#' @param nbvalmanplot The number of variables to include in each Manhattan plot.
#' @param Bonferroni Add a dashed bar to the Manhattan plot showing the Bonferroni significance level.
#' @param FDR Add a dashed bar to the Manhattan plot showing the False Discovery Rate (Benjamini-Hochberg method) significance threshold. NA if all p-values > FDR corrected p-values.
#' @param manplotsign Generates a Manhattan plot with only significant results (p<0.05).
#'
#' @import ggplot2 dplyr devtools
#' @rawNamespace import(stats, except = c("filter","lag"))
#' @rawNamespace import(MASS, except = "select")
#'
#' @return A Dataframe with results for each variable of the model.
#' @export
#'
#' @examples
#'#' ### Creating a dataframe with random variables
#'
#'
#' exposure1 <- as.factor(sample(c("Always", "Often","Never"), size = 400, replace = TRUE))
#' exposure2 <- as.factor(sample(c("Always", "Often","Never"), size = 400, replace = TRUE))
#' exposure3 <- as.factor(sample(c("Exposed", "Not exposed"), size = 400, replace = TRUE))
#' exposure4 <- as.factor(sample(c("Yes", "No"), size = 400, replace = TRUE))
#' exposure5 <- as.factor(sample(c("Yes", "No"), size = 400, replace = TRUE))
#' exposure6 <- as.factor(sample(c("Yes", "No"), size = 400, replace = TRUE))
#' exposure7 <- as.numeric(sample(0:400, 400, rep = TRUE))
#' exposure8 <- as.numeric(sample(0:400, 400, rep = TRUE))
#' exposure9 <- as.numeric(sample(0:0.1, 400, rep = TRUE))
#' outcome <- as.factor(sample(c("Yes", "No"), size = 400, replace = TRUE))
#'
#' data <- data.frame(exposure1,exposure2,exposure3,exposure4,exposure5,
#'                    exposure6,exposure7,exposure8,exposure9, outcome)
#'
#'
#' ### Launch of the EnvWAS analysis
#'
#' ELJAlogistic(var = 'outcome',data = data, var_adjust = c('exposure5'))
#'
#'@references
#'1. Dunn OJ. Multiple Comparisons Among Means. Journal of the American Statistical Association. 1961;56(293):52‑64.
#'2. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society: Series B (Methodological). 1995;57(1):289‑300.
#'
#'
ELJAlogistic <- function(var, var_adjust = NULL, data,
                        manplot = TRUE, nbvalmanplot = 100, Bonferroni = FALSE, FDR = FALSE, manplotsign = FALSE) {

  # Definition de y et x ( x = variable explicative dans la boucle + variable d'ajustements fixe s'il y en a)

  y <- data[, var]

  if (!is.null(var_adjust)) {
    x_adjust <- data[, var_adjust, drop = FALSE]
    x <- cbind(data[, -which(names(data) %in% c(var, var_adjust))], x_adjust)
    x_formula <- paste0(paste0(names(x_adjust), collapse = " + "), collapse = " + ")
    x_formula <- ifelse(length(x_formula) > 0, paste0(x_formula, " + "), "")
  } else {
    x <- data[, -which(names(data) == var)]
    x_formula <- ""
  }

  n <- length(x)


  # Analyses de regressions pour chaque variable explicative

  x_cols <- colnames(x)
  glm_results <- list()


  for (col in x_cols) {
    if(col %in% var_adjust) {
      next # passer à la variable explicative suivante si elle fait partie des variables d'ajustement
    }
    fit <- glm(as.formula(paste0(var, " ~ ", x_formula, col)), data = data, family = binomial(link = "logit"))
    summary_fit <- summary(fit)
    coef_names <- rownames(summary_fit$coefficients)[-1] # Exclure l'intercept

    for (name in coef_names) {
      level <- gsub("x\\[.*\\]", " ", name) # Extraire la modalite de la variable explicative
      odds <- exp(summary_fit$coefficients[name, 1])
      ci <- exp(suppressMessages(confint(fit, parm = name, level = 0.95)))
      p_value <- summary_fit$coefficients[name, 4]
      n <- nobs(fit)
      aic <- AIC(fit)
      glm_results[[paste0(col, "_", level)]] <- data.frame(
        variable = col,
        level = level,
        odd_ratio = odds,
        ci_low = ci[1],
        ci_high = ci[2],
        p_value = p_value,
        n = n,
        AIC = aic
      )
    }
  }

  # Assemblage des resultats dans un dataframe unique

  results <- do.call(rbind, glm_results)
  results$variable_level <- rownames(results)
  results$variable <- NULL


  # Suppression des OR associes aux variables d'ajustement

  if (!is.null(var_adjust)) {
    results <- results[!grepl(paste(var_adjust, collapse="|"), results$variable_level),]
    rownames(results) <- results$level
  }

  # Stockage du tableau des resultats et nettoyage de celui-ci

  results$variable_level <- NULL
  results <<- results

  nbvar <- as.numeric(length(results$p_value))

  # Calcul du FDR

  results_test <- results %>% arrange(p_value) %>% mutate(ligne = row_number())

  results_test$q <- (results_test$ligne / length(results_test$p_value)) * 0.05

  results_test$test <- results_test$p_value < results_test$q

  FDRcalc <- results_test %>% dplyr::filter(test == TRUE) %>% slice_tail(n = 1) %>% select(q) %>% as.numeric()

  # Traçage des Manhattan plots

  if (manplot == TRUE) {

    # Creation d'autant de Manhattan plot que de multiple de nbvalmanplot

    n_new <- ceiling(nrow(results) / nbvalmanplot) * nbvalmanplot
    results_new <- results[1:n_new, ]

    # Evite l'affichage de NA lorsque le nombre de Manhattan plot n'est pas un mutliple de nbvalmanplot
    n_missing <- n_new - nrow(results_new)
    if (n_missing > 0) {
      results_new[(nrow(results_new)+1):(n_new),] <- NA
    }

    results_list <- split(results_new, ceiling(seq_along(results_new$p_value)/nbvalmanplot))

    # Creer un manhattan plot pour chaque sous-groupe de donnees
    for (i in seq_along(results_list)) {
      sub_results <- results_list[[i]]

      # Verifier si le sous-groupe contient des valeurs manquantes
      if (any(is.na(sub_results$p_value))) {
        sub_results <- na.omit(sub_results)
        n_sub <- nrow(sub_results)
      } else {
        n_sub <- nbvalmanplot
      }

      # Creer le manhattan plot pour le sous-groupe
      Manplot <- ggplot(sub_results, aes(x = level, y = -log10(p_value),
                                        color = ifelse(odd_ratio > 1, "> 1", "< 1"))) +
        geom_point() +
        scale_color_manual(values = c("> 1" = "red", "< 1" = "green")) +
        coord_flip() +
        labs(color = "OR") +
        guides(color = guide_legend(title = "Odds ratio"))

      if (Bonferroni == TRUE & FDR == TRUE) {
          Manplot <- Manplot +
            geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
            geom_hline(aes(yintercept = -log10(0.05/nbvar), linetype = paste0( "Bonferroni : ",round(0.05/nbvar,5))), colour = "black",linewidth = 0.9) +
            geom_hline(aes(yintercept = -log10(FDRcalc), linetype = paste0( "FDR : ",round(FDRcalc,5))), colour = "purple",linewidth = 0.9) +
            scale_linetype_manual(name = "Threshold",values = c("dashed","dashed","dashed"),
                                  guide = guide_legend(override.aes = list(color = c("black","purple","#ffd800"))))
        } else if (Bonferroni == TRUE & FDR == FALSE) {
          Manplot <- Manplot +
            geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
            geom_hline(aes(yintercept = -log10(0.05/nbvar), linetype = paste0( "Bonferroni : ",round(0.05/nbvar,5))), colour = "black",linewidth = 0.9) +
            scale_linetype_manual(name = "Threshold",values = c("dashed","dashed"),
                                  guide = guide_legend(override.aes = list(color = c("black","#ffd800"))))
        } else if (Bonferroni == FALSE & FDR == TRUE) {
          Manplot <- Manplot +
            geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
            geom_hline(aes(yintercept = -log10(FDRcalc), linetype = paste0( "FDR : ",round(FDRcalc,5))), colour = "purple",linewidth = 0.9) +
            scale_linetype_manual(name = "Threshold",values = c("dashed","dashed"),
                                  guide = guide_legend(override.aes = list(color = c("purple","#ffd800"))))
        } else {
          Manplot <- Manplot +
            geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
            scale_linetype_manual(name = "Threshold",values = c("dashed"),
                                  guide = guide_legend(override.aes = list(color = c("#ffd800"))))
        }

      Manplot <<- Manplot
      print(Manplot)
    }

  } else {
    cat("The manhattan plot is not shown.\n")
  }

  if (manplotsign == TRUE) {

    # Creer un nouveau tableau de resultats avec les résultats uniquement significatifs
    resultssign <- results %>% dplyr::filter(p_value < 0.05)
    n_new_sign <- ceiling(nrow(resultssign) / nbvalmanplot) * nbvalmanplot
    resultssign_new <- resultssign[1:n_new_sign, ]

    # Ajouter des valeurs manquantes pour completer la derniere sous-liste
    n_missing_sign <- n_new_sign - nrow(resultssign_new)
    if (n_missing_sign > 0) {
      resultssign_new[(nrow(resultssign_new)+1):(n_new_sign),] <- NA
    }

    # Diviser le nouveau tableau de resultats en sous-groupes de 30 lignes
    results_list <- split(resultssign_new, ceiling(seq_along(resultssign_new$p_value)/nbvalmanplot))

    # Creer un manhattan plot pour chaque sous-groupe de donnees
    for (i in seq_along(results_list)) {
      sub_results <- results_list[[i]]

      # Verifier si le sous-groupe contient des valeurs manquantes
      if (any(is.na(sub_results$p_value))) {
        sub_results <- na.omit(sub_results)
        n_sub <- nrow(sub_results)
      } else {
        n_sub <- nbvalmanplot
      }

      # Creer le manhattan plot pour le sous-groupe
      Manplotsign <- ggplot(sub_results, aes(x = level, y = -log10(p_value),
                                            color = ifelse(odd_ratio > 1, "> 1", "< 1"))) +
        geom_point() +
        scale_color_manual(values = c("> 1" = "red", "< 1" = "green")) +
        coord_flip() +
        labs(color = "OR") +
        guides(color = guide_legend(title = "Odds ratio"))


      if (Bonferroni == TRUE & FDR == TRUE) {
        Manplotsign <- Manplotsign +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(0.05/nbvar), linetype = paste0( "Bonferroni : ",round(0.05/nbvar,5))), colour = "black",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(FDRcalc), linetype = paste0( "FDR : ",round(FDRcalc,5))), colour = "purple",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed","dashed","dashed"),
                                guide = guide_legend(override.aes = list(color = c("black","purple","#ffd800"))))
      } else if (Bonferroni == TRUE & FDR == FALSE) {
        Manplotsign <- Manplotsign +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(0.05/nbvar), linetype = paste0( "Bonferroni : ",round(0.05/nbvar,5))), colour = "black",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed","dashed"),
                                guide = guide_legend(override.aes = list(color = c("black","#ffd800"))))
      } else if (Bonferroni == FALSE & FDR == TRUE) {
        Manplotsign <- Manplotsign +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(FDRcalc), linetype = paste0( "FDR : ",round(FDRcalc,5))), colour = "purple",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed","dashed"),
                                guide = guide_legend(override.aes = list(color = c("purple","#ffd800"))))
      } else {
        Manplotsign <- Manplotsign +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed"),
                                guide = guide_legend(override.aes = list(color = c("#ffd800"))))
      }

      Manplotsign <<- Manplotsign
      print(Manplotsign)
    }

  } else {
    cat("The manhattan plot (manplotsign) is not shown.\n")
  }
}




#' Generalized Linear Models regression for EnvWAS/EWAS analysis
#'
#' @description
#' A tool for Environment-Wide Association Studies (EnvWAS / EWAS) which are repeated analysis. This function is espacially for generalized linear models 'glm' and allows the addition of adjustment variables.
#'
#'
#' @param var A categorical and binary variable. It is generally your outcome.
#' @param var_adjust A vector containing the names of the fixed adjustment variables for all the models.
#' @param family The family and the link use for the glm function.
#' @param data A dataframe containing all the variables needed for the analysis.
#' @param manplot Generate a Manhattan plot of the results of the analysis.
#' @param nbvalmanplot The number of variables to include in each Manhattan plot.
#' @param Bonferroni Add a dashed bar to the Manhattan plot showing the Bonferroni significance threshold.
#' @param FDR Add a dashed bar to the Manhattan plot showing the False Discovery Rate (Benjamini-Hochberg method) significance threshold. NA if all p-values > FDR corrected p-values.
#' @param manplotsign Generates a Manhattan plot with only significant results (p<0.05).
#'
#' @import ggplot2 dplyr devtools
#' @rawNamespace import(stats, except = c("filter","lag"))
#' @rawNamespace import(MASS, except = "select")
#'
#' @return A Dataframe with results for each variable of the model.
#' @export
#'
#' @examples
#' ### Creating a dataframe with random variables
#'
#'
#' exposure1 <- as.factor(sample(c("Always", "Often","Never"), size = 400, replace = TRUE))
#' exposure2 <- as.factor(sample(c("Always", "Often","Never"), size = 400, replace = TRUE))
#' exposure3 <- as.factor(sample(c("Exposed", "Not exposed"), size = 400, replace = TRUE))
#' exposure4 <- as.factor(sample(c("Yes", "No"), size = 400, replace = TRUE))
#' exposure5 <- as.factor(sample(c("Yes", "No"), size = 400, replace = TRUE))
#' exposure6 <- as.factor(sample(c("Yes", "No"), size = 400, replace = TRUE))
#' exposure7 <- as.numeric(sample(0:400, 400, rep = TRUE))
#' exposure8 <- as.numeric(sample(0:400, 400, rep = TRUE))
#' exposure9 <- as.numeric(sample(0:0.1, 400, rep = TRUE))
#' outcome <- as.factor(sample(c("Yes", "No"), size = 400, replace = TRUE))
#'
#' data <- data.frame(exposure1,exposure2,exposure3,exposure4,exposure5,
#'                    exposure6,exposure7,exposure8,exposure9, outcome)
#'
#'
#' ### Launch of the EnvWAS analysis
#'
#' ELJAglm(var = 'outcome',data = data, var_adjust = c('exposure5'))
#'
#'@references
#'1. Dunn OJ. Multiple Comparisons Among Means. Journal of the American Statistical Association. 1961;56(293):52‑64.
#'2. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society: Series B (Methodological). 1995;57(1):289‑300.
#'
#'
ELJAglm <- function(var, var_adjust = NULL, family = binomial(link = "logit"), data,
                        manplot = TRUE, nbvalmanplot = 100, Bonferroni = FALSE, FDR = FALSE, manplotsign = FALSE) {

  # Definition de y et x ( x = variable explicative dans la boucle + variable d'ajustements fixe s'il y en a)

  y <- data[, var]

  if (!is.null(var_adjust)) {
    x_adjust <- data[, var_adjust, drop = FALSE]
    x <- cbind(data[, -which(names(data) %in% c(var, var_adjust))], x_adjust)
    x_formula <- paste0(paste0(names(x_adjust), collapse = " + "), collapse = " + ")
    x_formula <- ifelse(length(x_formula) > 0, paste0(x_formula, " + "), "")
  } else {
    x <- data[, -which(names(data) == var)]
    x_formula <- ""
  }

  n <- length(x)


  # Analyses de regressions pour chaque variable explicative

  x_cols <- colnames(x)
  glm_results <- list()


  for (col in x_cols) {
    if(col %in% var_adjust) {
      next # passer à la variable explicative suivante si elle fait partie des variables d'ajustement
    }
    fit <- glm(as.formula(paste0(var, " ~ ", x_formula, col)), data = data, family = family)
    summary_fit <- summary(fit)
    coef_names <- rownames(summary_fit$coefficients)[-1] # Exclure l'intercept

    for (name in coef_names) {
      level <- gsub("x\\[.*\\]", " ", name) # Extraire la modalite de la variable explicative
      coef <- summary_fit$coefficients[name, 1]
      ci <- suppressMessages(confint(fit, parm = name, level = 0.95))
      p_value <- summary_fit$coefficients[name, 4]
      n <- nobs(fit)
      aic <- AIC(fit)
      glm_results[[paste0(col, "_", level)]] <- data.frame(
        variable = col,
        level = level,
        coefficients = coef,
        ci_low = ci[1],
        ci_high = ci[2],
        p_value = p_value,
        n = n,
        AIC = aic
      )
    }
  }

  # Assemblage des resultats dans un dataframe unique

  results <- do.call(rbind, glm_results)
  results$variable_level <- rownames(results)
  results$variable <- NULL


  # Suppression des coefficients associes aux variables d'ajustement

  if (!is.null(var_adjust)) {
    results <- results[!grepl(paste(var_adjust, collapse="|"), results$variable_level),]
    rownames(results) <- results$level
  }

  # Stockage du tableau des resultats et nettoyage de celui-ci

  results$variable_level <- NULL
  results <<- results

  nbvar <- as.numeric(length(results$p_value))

  # Calcul du FDR

  results_test <- results %>% arrange(p_value) %>% mutate(ligne = row_number())

  results_test$q <- (results_test$ligne / length(results_test$p_value)) * 0.05

  results_test$test <- results_test$p_value < results_test$q

  FDRcalc <- results_test %>% dplyr::filter(test == TRUE) %>% slice_tail(n = 1) %>% select(q) %>% as.numeric()

  # Traçage des Manhattan plots

  if (manplot == TRUE) {

    # Creation d'autant de Manhattan plot que de multiple de nbvalmanplot

    n_new <- ceiling(nrow(results) / nbvalmanplot) * nbvalmanplot
    results_new <- results[1:n_new, ]

    # Evite l'affichage de NA lorsque le nombre de Manhattan plot n'est pas un mutliple de nbvalmanplot
    n_missing <- n_new - nrow(results_new)
    if (n_missing > 0) {
      results_new[(nrow(results_new)+1):(n_new),] <- NA
    }

    results_list <- split(results_new, ceiling(seq_along(results_new$p_value)/nbvalmanplot))

    # Creer un manhattan plot pour chaque sous-groupe de donnees
    for (i in seq_along(results_list)) {
      sub_results <- results_list[[i]]

      # Verifier si le sous-groupe contient des valeurs manquantes
      if (any(is.na(sub_results$p_value))) {
        sub_results <- na.omit(sub_results)
        n_sub <- nrow(sub_results)
      } else {
        n_sub <- nbvalmanplot
      }

      # Creer le manhattan plot pour le sous-groupe
      Manplot <- ggplot(sub_results, aes(x = level, y = -log10(p_value),
                                        color = ifelse(coefficients > 0, "> 0", "< 0"))) +
        geom_point() +
        scale_color_manual(values = c("> 0" = "red", "< 0" = "green")) +
        coord_flip() +
        labs(color = "coefficients") +
        guides(color = guide_legend(title = "Coefficients"))

      if (Bonferroni == TRUE & FDR == TRUE) {
          Manplot <- Manplot +
            geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
            geom_hline(aes(yintercept = -log10(0.05/nbvar), linetype = paste0( "Bonferroni : ",round(0.05/nbvar,5))), colour = "black",linewidth = 0.9) +
            geom_hline(aes(yintercept = -log10(FDRcalc), linetype = paste0( "FDR : ",round(FDRcalc,5))), colour = "purple",linewidth = 0.9) +
            scale_linetype_manual(name = "Threshold",values = c("dashed","dashed","dashed"),
                                  guide = guide_legend(override.aes = list(color = c("black","purple","#ffd800"))))
        } else if (Bonferroni == TRUE & FDR == FALSE) {
          Manplot <- Manplot +
            geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
            geom_hline(aes(yintercept = -log10(0.05/nbvar), linetype = paste0( "Bonferroni : ",round(0.05/nbvar,5))), colour = "black",linewidth = 0.9) +
            scale_linetype_manual(name = "Threshold",values = c("dashed","dashed"),
                                  guide = guide_legend(override.aes = list(color = c("black","#ffd800"))))
        } else if (Bonferroni == FALSE & FDR == TRUE) {
          Manplot <- Manplot +
            geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
            geom_hline(aes(yintercept = -log10(FDRcalc), linetype = paste0( "FDR : ",round(FDRcalc,5))), colour = "purple",linewidth = 0.9) +
            scale_linetype_manual(name = "Threshold",values = c("dashed","dashed"),
                                  guide = guide_legend(override.aes = list(color = c("purple","#ffd800"))))
        } else {
          Manplot <- Manplot +
            geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
            scale_linetype_manual(name = "Threshold",values = c("dashed"),
                                  guide = guide_legend(override.aes = list(color = c("#ffd800"))))
        }

      Manplot <<- Manplot
      print(Manplot)
    }

  } else {
    cat("The manhattan plot is not shown.\n")
  }

  if (manplotsign == TRUE) {

    # Creer un nouveau tableau de resultats avec autant de multiples de 30
    resultssign <- results %>% dplyr::filter(p_value < 0.05)
    n_new_sign <- ceiling(nrow(resultssign) / nbvalmanplot) * nbvalmanplot
    resultssign_new <- resultssign[1:n_new_sign, ]

    # Ajouter des valeurs manquantes pour completer la derniere sous-liste
    n_missing_sign <- n_new_sign - nrow(resultssign_new)
    if (n_missing_sign > 0) {
      resultssign_new[(nrow(resultssign_new)+1):(n_new_sign),] <- NA
    }

    # Diviser le nouveau tableau de resultats en sous-groupes de 30 lignes
    results_list <- split(resultssign_new, ceiling(seq_along(resultssign_new$p_value)/nbvalmanplot))

    # Creer un manhattan plot pour chaque sous-groupe de donnees
    for (i in seq_along(results_list)) {
      sub_results <- results_list[[i]]

      # Verifier si le sous-groupe contient des valeurs manquantes
      if (any(is.na(sub_results$p_value))) {
        sub_results <- na.omit(sub_results)
        n_sub <- nrow(sub_results)
      } else {
        n_sub <- nbvalmanplot
      }

      # Creer le manhattan plot pour le sous-groupe
      Manplotsign <- ggplot(sub_results, aes(x = level, y = -log10(p_value),
                                            color = ifelse(coefficients > 0, "> 0", "< 0"))) +
        geom_point() +
        scale_color_manual(values = c("> 0" = "red", "< 0" = "green")) +
        coord_flip() +
        labs(color = "Coefficients") +
        guides(color = guide_legend(title = "Coefficients"))


      if (Bonferroni == TRUE & FDR == TRUE) {
        Manplotsign <- Manplotsign +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(0.05/nbvar), linetype = paste0( "Bonferroni : ",round(0.05/nbvar,5))), colour = "black",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(FDRcalc), linetype = paste0( "FDR : ",round(FDRcalc,5))), colour = "purple",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed","dashed","dashed"),
                                guide = guide_legend(override.aes = list(color = c("black","purple","#ffd800"))))
      } else if (Bonferroni == TRUE & FDR == FALSE) {
        Manplotsign <- Manplotsign +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(0.05/nbvar), linetype = paste0( "Bonferroni : ",round(0.05/nbvar,5))), colour = "black",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed","dashed"),
                                guide = guide_legend(override.aes = list(color = c("black","#ffd800"))))
      } else if (Bonferroni == FALSE & FDR == TRUE) {
        Manplotsign <- Manplotsign +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          geom_hline(aes(yintercept = -log10(FDRcalc), linetype = paste0( "FDR : ",round(FDRcalc,5))), colour = "purple",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed","dashed"),
                                guide = guide_legend(override.aes = list(color = c("purple","#ffd800"))))
      } else {
        Manplotsign <- Manplotsign +
          geom_hline(aes(yintercept = -log10(0.05), linetype = paste0( "Unadjusted : 0.05")), colour = "#ffd800",linewidth = 0.9) +
          scale_linetype_manual(name = "Threshold",values = c("dashed"),
                                guide = guide_legend(override.aes = list(color = c("#ffd800"))))
      }

      Manplotsign <<- Manplotsign
      print(Manplotsign)
    }

  } else {
    cat("The manhattan plot (manplotsign) is not shown.\n")
  }
}


