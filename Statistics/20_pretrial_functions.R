library(car)
library(emmeans)
library(multcomp)
library(here)
library(tidyverse)
library(car)
library(lmerTest)
library(LambertW)
library(rcompanion)
library(msm)
library(gt)

# ANOVA functions

test_models <- function(df, model_selection = "fixed") {
  if (model_selection == "AIC") {
    model_1 <- aov(response ~ bacteria, data = df) # easiest model
    model_2 <- lmer(response ~ bacteria + (1|block), data = df)
    actual_AIC = 9999
    for (i in c(1:2)) { ###################################################
      model <- get(paste0("model_", i))
      print(paste("Model", i, round(AIC(model),2)))
      singular = FALSE
      if(any(class(model) %in% "lmerModLmerTest")) {
        print(isSingular(model))
        singular = isSingular(model)
      } # evaluates singularity of model
      if (AIC(model) < actual_AIC & isFALSE(singular)) {
        actual_AIC <- AIC(model)
        best_model <- model
      }
    }
  } else if (model_selection == "fixed") {
    model <- lmer(response ~ bacteria + (1|block), data = df)
    if (VarCorr(model)["block"] == 0) {
      best_model <- aov(response ~ bacteria + block, data = df)
    } else {
      best_model <- model
    }
  }
  print(best_model)
  return(best_model)
}

create_n_table <- function(df) {
  t <- as.data.frame(table(df$bacteria)) %>%
    dplyr::rename("bacteria" = Var1, "n" = Freq)
  return(t)
}

evaluate_model <- function(model) {
  if(any(class(model) %in% "lmerModLmerTest")) {
    plot(model)
    #ranova(model)
  }
  plot(residuals(model))
  abline(a = 0, b = 0)
  plot(abs(residuals(model)))
  qqPlot(resid(model))
  print(shapiro.test(resid(model)))
  #print(ks.test(resid(model), "pnorm"))
}

do_anova <- function(model) {
  if(any(class(model) %in% "lmerModLmerTest")) {
    anova <- anova(model, type = "III", ddf="Kenward-Roger")
  } else {
    anova <- Anova(model, type = "III")
  }
  print(anova)
  return(anova)
}

pairwise_comparisons <- function(model) {
  emmeans <- emmeans(model, ~bacteria)
  cld <- cld(emmeans, Letters = letters, details = F, reversed = T)
  return(cld)
}

pairwise_comparisons_two_factorial <- function(model) {
  emmeans <- emmeans(model, ~bacteria+enzyme+bacteria:enzyme)
  cld <- cld(emmeans, Letters = letters, details = F, reversed = T)
  return(cld)
}

########################### own transform tukey function

transformTukey2 <- function(df, start = -10, end = 10, int = 0.25) {
  result_frame <- data.frame(transformation = seq(start, end, int),
                             W = 0)
  for (i in 1:nrow(result_frame)) {
    transformation_factor <- result_frame$transformation[i]
    transformed_df <- df
    # test transformations
    if(transformation_factor > 0) {
      transformed <- df$response^(transformation_factor)
      transformed_df$response <- transformed
    } else if (transformation_factor == 0) {
      transformed <- log(df$response)
      transformed_df$response <- transformed
    } else if (transformation_factor < 0) {
      transformed <- -1 * df$response^transformation_factor
      transformed_df$response <- transformed
    }
    if (any(is.infinite(transformed_df$response)) == FALSE & any(is.nan(transformed_df$response)) == FALSE) { # only validate if no NA or infinities in transformed data
      x <- suppressMessages(suppressWarnings(
        { capture.output(model <- test_models(transformed_df)) }
      ))
      test_statistic <- shapiro.test(resid(model))
      if (deviance(model) < sqrt(.Machine$double.eps)) { # if precision gets to low (similar to test in car::Anova)
        result_frame$W[i] <- 0 # do not use this transformation
      } else {
        result_frame$W[i] <- test_statistic$statistic
      }
    }
  }
  plot(x = result_frame$transformation, y = result_frame$W)
  W <- result_frame$transformation[which.max(result_frame$W)] 
  print(W)
  # generate output
  if(W > 0) {
    transformed <- df$response^(W)
  } else if (W == 0) {
    transformed <- log(df$response)
  } else if (W < 0) {
    transformed <- -1 * df$response^W
  }
  transformed_list <- list(transformed, W)
  return(transformed_list)
}

# function to perform transformations on concentration column

transform_data <- function(df, transformation = c("gauss", "tukey", "tukey2", "log", "logit"), type = c("h", "hh", "s")) {
  if (transformation == "gauss") {
    transformed <- Gaussianize(df$response, type = type, return.tau.mat = T)
    
  } else if (transformation == "tukey") {
    
    transformed1 <- transformTukey(df$response)
    transformed2 <- transformTukey(df$response, returnLambda = T)
    transformed <- list(transformed1, transformed2)
  } else if (transformation == "tukey2") { # leads to choosing the worst model
    transformed <- transformTukey2(df)
  } else if (transformation == "log") {
    transformed1 <- log(df$response)
    transformed2 <- 1
    transformed <- list(transformed1, transformed2)
  } else if (transformation == "logit") {
    transformed1 <- log((df$response/100)/(1-(df$response/100)))
    transformed2 <- 1
    transformed <- list(transformed1, transformed2)
  }
  df$response <- transformed[[1]][1:length(transformed[[1]])]
  transformation_factor <- transformed[[2]]
  transformation <- list(df, transformation_factor)
  return(transformation)
}

# function for backtransformation of means

backtransform_data <- function(means, transformation = c("gauss", "tukey", "tukey2", "log", "logit"), transformation_factor) {
  if (transformation == "gauss") {
    means$emmean <- Gaussianize(means$emmean, inverse = TRUE, tau.mat = transformation_factor)[,1]
    means$lower.CL <- Gaussianize(means$lower.CL, inverse = TRUE, tau.mat = transformation_factor)[,1]
    means$upper.CL <- Gaussianize(means$upper.CL, inverse = TRUE, tau.mat = transformation_factor)[,1]
    means$SE <- (means$upper.CL - means$lower.CL) / (2 * 1.96)
  } else if (transformation %in% c("tukey", "tukey2")) {
    if(transformation_factor > 0) {
      means$emmean <- means$emmean^(1/transformation_factor)
      means$lower.CL <- means$lower.CL^(1/transformation_factor)
      means$upper.CL <- means$upper.CL^(1/transformation_factor)
      means$SE <- (means$upper.CL - means$lower.CL) / (2 * 1.96)
    } else if (transformation_factor == 0) {
      means$emmean <- exp(means$emmean)
      means$lower.CL <- exp(means$lower.CL)
      means$upper.CL <- exp(means$upper.CL)
      means$SE <- (means$upper.CL - means$lower.CL) / (2 * 1.96)
    } else if (transformation_factor < 0) {
      means$emmean <- (-1 * means$emmean)^(1/transformation_factor)
      means$lower.CL <- (-1 * means$lower.CL)^(1/transformation_factor)
      means$upper.CL <- (-1 * means$upper.CL)^(1/transformation_factor)
      means$SE <- (means$upper.CL - means$lower.CL) / (2 * 1.96)
    }
  } else if (transformation == "log") {
    means$emmean <- exp(means$emmean)
    means$lower.CL <- exp(means$lower.CL)
    means$upper.CL <- exp(means$upper.CL)
    means$SE <- (means$upper.CL - means$lower.CL) / (2 * 1.96)
  } else if (transformation == "logit") {
    means$emmean <- 100*(exp(means$emmean)/(exp(means$emmean) + 1))
    means$lower.CL <- 100*(exp(means$lower.CL)/(exp(means$lower.CL) + 1))
    means$upper.CL <- 100*(exp(means$upper.CL)/(exp(means$upper.CL) + 1))
    means$SE <- (means$upper.CL - means$lower.CL) / (2 * 1.96)
  }
  return(means)
}

# function to test different transformations and find the best one according to shapiro wilk test

test_transformations <- function(df) {
  result_frame <- data.frame(transformation = c("none", "tukey", "tukey2", "log", "logit"),
                             W = 0)
  for (i in 1:nrow(result_frame)) {
    transformation <- result_frame$transformation[i]
    if (transformation == "none") {
      transformed_df <- df
    } else {
      transformed_df <- transform_data(df, transformation = transformation)[[1]]
    }
    
    if (NaN %in% transformed_df$response) {
      result_frame$W[i] <- 0
    } else {
      model <- test_models(transformed_df)
      test_statistic <- shapiro.test(resid(model))
      result_frame$W[i] <- test_statistic$statistic
    }
    
  }
  best_transformation <- result_frame$transformation[which.max(result_frame$W)]
  return(best_transformation)
}


combined_comparison <- function(df,transformation = "none", 
                                type = c("h", "hh", "s")) { #, correction_filter
  if (transformation == "test") {
    transformation = test_transformations(df)
    print(paste("Best transformation:", transformation))
  }
  
  if (transformation == "none") {
    
    model <- test_models(df)
    
    n <- create_n_table(df)
    
    evaluate_model(model)
    
    ano <- do_anova(model)
    
    cld <- pairwise_comparisons(model)
    cld2 <- cld
    
  } else {
    
    transformed <- transform_data(df, transformation = transformation, type = type)
    transformed_df <- transformed[[1]]
    transformation_factor <- transformed[[2]]
    
    model <- test_models(transformed_df)
    
    n <- create_n_table(df)
    
    evaluate_model(model)
    
    ano <- do_anova(model)
    
    cld <- pairwise_comparisons(model)
    
    cld <- backtransform_data(cld, transformation = transformation, transformation_factor = transformation_factor)
    
    model_untransformed <- test_models(df)
    cld2 <- pairwise_comparisons(model_untransformed)
    
  }
  
  return(list(ano = ano, cld = cld, cld2 = cld2, n = n))
}


create_results_table <- function(input_df, comparison_object, response_name = "response", digits = 1) {
  P <- comparison_object$ano[which(rownames(comparison_object$ano) == "bacteria"), 
                             which(colnames(comparison_object$ano) == "Pr(>F)")]
  Print <- ifelse(P < 0.001, paste("< 0.001"), paste(format(round(P,3), nsmall = 3)))
  n <- comparison_object$n
  cld <- comparison_object$cld
  cld2 <- comparison_object$cld2 %>%
    inner_join(n, by = "bacteria")
  pSEM <- paste(format(round(sum(cld2$SE*(cld2$n-1)) / (sum(n$n)-length(n$n)),digits+1), nsmall = digits+1))
  table <- dplyr::select(cld2, bacteria, emmean, SE) %>%
    inner_join(dplyr::select(cld, bacteria, .group), by = "bacteria") %>%
    arrange(bacteria) %>%
    inner_join(distinct(dplyr::select(input_df, bacteria, species)), by = "bacteria") %>%
    mutate(!!sym(response_name) := paste0(format(round(emmean, digits), nsmall = digits),
                                          if (length(unique(cld$.group))==1) {""} else {str_trim(.group)})) %>%
    #dplyr::select(diet = description, !!sym(response_name)) %>%
    dplyr::select(bacteria, !!sym(response_name)) %>%
    add_row(bacteria = "Pooled SEM", !!sym(response_name) := pSEM) %>%
    add_row(bacteria = "P-value", !!sym(response_name) := Print)
  return(table)
}

create_results_plot <- function(input_df, comparison_object, response_name = "response", y_axis = "value") {
  P <- comparison_object$ano[which(rownames(comparison_object$ano) == "bacteria"), 
                             which(colnames(comparison_object$ano) == "Pr(>F)")]
  Print <- ifelse(P < 0.001, paste("< 0.001"), paste("=", round(P,3)))
  n <- comparison_object$n
  cld <- comparison_object$cld
  cld2 <- comparison_object$cld2 %>%
    inner_join(n, by = "bacteria")
  pSEM <- sum(cld2$SE*(cld2$n-1)) / (sum(n$n)-length(n$n))
  table <- dplyr::select(cld2, bacteria, emmean, SE) %>%
    inner_join(dplyr::select(cld, bacteria, .group), by = "bacteria") %>%
    arrange(bacteria) %>%
    inner_join(distinct(dplyr::select(input_df, bacteria, species)), by = "bacteria") %>%
    #mutate(description = factor(description, levels = c("Pea3", "Pea6", "Pea7", "Pea11"))) %>%
    mutate(lower = emmean - SE, upper = emmean + SE)
  p <- ggplot(table, aes(x = bacteria)) +
    geom_bar(aes(y = emmean), stat = "identity", width = 0.4) +
    geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.1, size = 1) +
    annotate("label", label = paste("REML\n P-value", Print), x = 7, y = max(table$upper)/10, size = 6, fill = "snow2") +
    {if(P < 0.05)geom_text(aes(y = (emmean + SE)*1.05, label = str_trim(.group)), size = 5, position = position_nudge(x = 0))}+
    scale_y_continuous(limits = c(min(0, table$lower),1.1*max(table$upper)), expand = expansion(mult = c(0, .1))) +
    labs(x = "bacteria", y = y_axis, title = response_name) +
    guides(x=guide_axis(angle = 30))
  print(p)
}

save_results_table <- function(table, title, name) {
  gt <- gt(table, rowname_col = "bacteria") %>%
    tab_header(title = title) %>%
    opt_align_table_header(align = "left") %>%
    tab_options(table.border.top.color = "white",
                heading.title.font.size = px(20),
                column_labels.border.top.width = 3,
                column_labels.border.top.color = "black", 
                column_labels.border.bottom.width = 2,
                column_labels.border.bottom.color = "black",
                table_body.border.bottom.color = "black",
                table_body.border.bottom.width = 3,
                table.border.bottom.color = "white",
                table.width = pct(100),
                table.background.color = "white") %>%
    tab_style(style = cell_borders(sides = c("top", "bottom"),color = "white"),
              locations = cells_body()) %>%
    tab_style(style = cell_borders(sides = c("top", "bottom", "right"),color = "white"),
              locations = cells_stub(rows = TRUE)) %>%
    cols_align(align="center", columns = -bacteria) %>%
    opt_table_font(font = google_font("Merriweather"))
  gtsave(gt, filename = paste0(name, ".png"), "plots")
  return(gt)
}

create_results_plot2 <- function(cld_object, response_name = "response", y_axis = "value") {
  P = 0
  cld <- cld_object
  pSEM <- mean(cld$SE)
  table <- dplyr::select(cld, dm, inoculation, emmean, SE, .group) %>%
    mutate(lower = emmean - SE, upper = emmean + SE)
  p <- ggplot(table, aes(x = dm, fill = inoculation)) +
    geom_bar(aes(y = emmean), position = position_dodge(width = 0.5), stat = "identity", width = 0.4) +
    geom_errorbar(aes(ymax = upper, ymin = lower), position = position_dodge(width = 0.5), width = 0.1, size = 1) +
    {if(P < 0.05)geom_text(aes(y = (emmean + SE)*1.05, label = str_trim(.group)), position = position_dodge(width = 0.5), size = 5)}+
    scale_y_continuous(limits = c(min(0, table$lower),1.1*max(table$upper)), expand = expansion(mult = c(0, .1))) +
    labs(x = "pea ratio", y = y_axis, title = response_name) +
    scale_fill_manual(values = c("grey70", "grey50", "grey30"))
  print(p)
}

# select column from combined frame

select_response <- function(df, selected_response) {
  filtered_df <- df %>%
    dplyr::select(sampleid, bacteria, block,
                  response = !!sym(selected_response))
}

select_response_two_factorial <- function(df, selected_response) {
  filtered_df <- df %>%
    dplyr::select(sampleid, bacteria, enzyme, block, box,
                  response = !!sym(selected_response))
}

# plotting function

plot_pairwise <- function(filtered_df, output, selected_response, save = F) { #, correction_filter
  cld <- output$cld
  P <- output$ano[which(rownames(output$ano) == "bacteria"), which(colnames(output$ano) == "Pr(>F)")]
  Print <- ifelse(P < 0.001, paste("< 0.001"), paste("=", round(P,3)))
  p <- ggplot(cld, aes(x = bacteria)) +
    geom_boxplot(data = filtered_df, aes(y = response), width = 0.2, position = position_nudge(x = -0.1)) +
    geom_point(data = filtered_df, aes(y = response, shape = block), 
               size = 2, position = position_nudge(x = -0.1), stroke = 2) +
    #geom_point(aes(y = emmean), size = 2, position = position_nudge(x = 0)) +
    #geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0, lwd = 1, position = position_nudge(x = 0)) +
    {if(P < 0.05)geom_text(aes(y = emmean, label = str_trim(.group)), size = 5, position = position_nudge(x = 0.2))}+
    labs(x = "bacteria", y = "value", title = paste(selected_response)) + #, correction_filter
    annotate("label", label = paste("bacteria\n P", Print), x = 8, y = max(filtered_df$response), size = 5)
  print(p)
}

# function combining everything

combined_comparison_nutrition <- function(df, selected_response, transformation = "none", 
                                          type = c("h", "hh", "s")) { #, correction_filter
  filtered_df <- select_response(df, selected_response = selected_response) #, correction_filter
  
  out <- combined_comparison(df = filtered_df, transformation = transformation, type = type)
  
  plot_pairwise(filtered_df, out, selected_response, save = save) #, correction_filter
  
  return(out)
}

# one liner function for each response

combined_comparison_and_results <- function(df, selected_response, transformation = "none",
                                            type = c("h", "hh", "s"), response_name, y_axis,
                                            save_name, save, digits = 1) {
  script_path <- rstudioapi::getSourceEditorContext()$path
  script_name <- basename(script_path)
  number <- str_extract(script_name, "[0-9]+")
  
  out <- combined_comparison_nutrition(df = df, selected_response = selected_response, 
                                       transformation = transformation, type = type)
  if (isTRUE(save)) {
    save_big(paste0(number, "_", save_name))
  }
  
  table <- create_results_table(input_df = df, comparison_object = out, response_name = response_name, digits = digits)
  
  create_results_plot(input_df = df, comparison_object = out, response_name = response_name, 
                      y_axis = y_axis)
  if (isTRUE(save)) {
    save_big(paste0(number, "_result_", save_name))
  }
  return(table)
}



save_results_table <- function(table, title, name) {
  gt <- gt(table, rowname_col = "bacteria") %>%
    tab_header(title = title) %>%
    opt_align_table_header(align = "left") %>%
    tab_options(table.border.top.color = "white",
                heading.title.font.size = px(20),
                column_labels.border.top.width = 3,
                column_labels.border.top.color = "black", 
                column_labels.border.bottom.width = 2,
                column_labels.border.bottom.color = "black",
                table_body.border.bottom.color = "black",
                table_body.border.bottom.width = 3,
                table.border.bottom.color = "white",
                table.width = pct(100),
                table.background.color = "white") %>%
    tab_style(style = cell_borders(sides = c("top", "bottom"),color = "white"),
              locations = cells_body()) %>%
    tab_style(style = cell_borders(sides = c("top", "bottom", "right"),color = "white"),
              locations = cells_stub(rows = TRUE)) %>%
    cols_align(align="center", columns = -bacteria) %>%
    opt_table_font(font = google_font("Merriweather"))
  gtsave(gt, filename = paste0(name, ".png"), "plots")
  return(gt)
}
