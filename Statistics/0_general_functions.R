library(here)
library(tidyverse)
library(ggbeeswarm)
library(emmeans)
library(multcomp)
library(car)
library(lmerTest)
library(LambertW)
library(rcompanion)
library(msm)
library(gt)
library(venn)


source("E:/R/source/ggplot2_theme_bw.R")

# save plots as jpeg

save_big <- function(name, width = 20, height = 15) {
  ggsave(filename= paste0(name, ".jpeg"),
         plot = last_plot(),
         device= "jpeg", 
         path = "plots", 
         units = "cm", 
         width = width,
         height = height,
         scale=1,
         dpi="print")
}

# create venn diagramm

save_venn <- function(input_df, save_name, save = FALSE){
  if (isTRUE(save)) jpeg(filename = paste0("plots/", save_name, ".jpeg"), width = 7.7, height = 7.7, unit="cm", res = 1000)
  venn(input_df, 
       zcolor = colors,
       ilabels = "counts",
       ilcs = 0.6,
       sncs = 0.8,
       box = FALSE)
  if (isTRUE(save)) dev.off()
}


# save output as html

save_output <- function() {
  script_path <- rstudioapi::getSourceEditorContext()$path
  script_name <- basename(script_path)
  rmarkdown::render(input = script_name, 
                    output_file = paste0(str_remove(script_name, "\\.R"), 
                                         str_extract(Sys.time(), "\\d+.\\d+.\\d+"), ".html"), 
                    output_dir = "output")
}

# function to retrieve number of current script

get_script_number <- function() {
  script_path <- rstudioapi::getSourceEditorContext()$path
  script_name <- basename(script_path)
  number <- str_extract(script_name, "[0-9]+")
  return(number)
}

# function taking a P value and format it for printing
# smaller than 0.001 is transformed into < 0.001, otherwise rounded to 3 digits
# < or = gets included

round_p_value <- function(P) {
  Print <- ifelse(P < 0.001, paste("< 0.001"), paste("=", format(round(P,3), nsmall = 3)))
  return(Print)
}

# read in function for all nutrition tsv's

read_in_nutrition <- function(path) {
  tibble <- read_tsv(path) %>%
    rename_all(.funs = ~ gsub("\\s+|-+|\\,+|\\(|\\)","_", .) %>% gsub("_$","", .) %>% gsub("_+","_",.) %>% tolower) %>%
    mutate(across(c(animal, period, square, diet), as.character))
  return(tibble)
}

# select column from combined frame

select_response <- function(df, selected_response) {
  filtered_df <- df %>%
    dplyr::select(animal, period, square, diet,
                  response = !!sym(selected_response))
  return(filtered_df)
}

# separate ileal and faecal samples

select_matrix <- function(df, selected_matrix) {
  filtered_df <- df %>%
    filter(matrix == selected_matrix) %>%
    dplyr::select(-matrix)
  return(filtered_df)
}

# ANOVA functions

test_models <- function(df, model_style = "klein", model_selection = "fixed") {
  if (model_selection == "AIC") {
    model_1 <- aov(response ~ diet, data = df) # easiest model
    if (model_style == "klein") {
      model_2 <- lmer(response ~ diet + (1|animal) + (1|period), data = df)
      model_3 <- lmer(response ~ diet + (1|animal), data = df) # reduced model
      model_4 <- lmer(response ~ diet + (1|period), data = df) # reduced model
    } else if (model_style == "heyer") {
      # heyer models
      model_2 <- lmer(response ~ diet + (1|square) + (1|square:animal) + (1|square:period), data = df)
      model_3 <- lmer(response ~ diet + (1|square/animal), data = df) # reduced model
      model_4 <- lmer(response ~ diet + (1|square/period), data = df) # reduced model
    } else if (model_style == "fixed") {
      model_1 <- aov(response ~ diet + animal + period, data = df)
      model_2 <- aov(response ~ diet + animal + period, data = df)
      model_3 <- aov(response ~ diet + animal + period, data = df)
      model_4 <- aov(response ~ diet + animal + period, data = df)
    }
    actual_AIC = 9999
    for (i in c(1:4)) { ###################################################
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
    if (model_style == "klein") {
      model <- lmer(response ~ diet + (1|animal) + (1|period), data = df)
      if (VarCorr(model)["animal"] == 0 & VarCorr(model)["period"] == 0) {
        best_model <- aov(response ~ diet + animal + period, data = df)
      } else if (VarCorr(model)["animal"] == 0) {
        best_model = lmer(response ~ diet + animal + (1|period), data = df)
      } else if (VarCorr(model)["period"] == 0) {
        best_model = lmer(response ~ diet + (1|animal) + period, data = df)
      } else
        best_model <- model
    } else if (model_style == "heyer") {
      # heyer models
      best_model <- lmer(response ~ diet + (1|square) + (1|square:animal) + (1|square:period), data = df)
    } else if (model_style == "fixed") {
      best_model <- aov(response ~ diet + animal + period, data = df)
    } else if (model_style == "ferment") {
      model <- lmer(response ~ diet + (1|period), data = df)
      if (VarCorr(model)["period"] == 0) {
        best_model <- aov(response ~ diet + period, data = df)
      } else {
        best_model <- model
      }
    }
  }
  print(best_model)
  return(best_model)
}

create_n_table <- function(df) {
  t <- as.data.frame(table(df$diet)) %>%
    dplyr::rename("diet" = Var1, "n" = Freq)
  return(t)
}

evaluate_model <- function(model) {
  if(any(class(model) %in% "lmerModLmerTest")) {
    plot(model)
    #ranova(model)
  }
  plot(residuals(model))
  abline(a = 0, b = 0)
  # do residuals have constant variance?
  plot(fitted(model), residuals(model))
  abline(a = 0, b = 0)
  # linearity of predictors (not relevant for categorical predictors?)
  plot(model.response(model.frame(model)), residuals(model))
  abline(a = 0, b = 0)
  plot(abs(residuals(model)))
  plot(rstudent(model))
  abline(a = 0, b = 0)
  qqPlot(resid(model))
  shapiro <- shapiro.test(resid(model))
  #print(ks.test(resid(model), "pnorm"))
  return(shapiro)
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
  emmeans <- emmeans(model, ~diet)
  contrasts <- contrast(emmeans, "pairwise")
  cld <- cld(emmeans, Letters = letters, details = F, reversed = T)
  print(contrasts)
  return(cld)
}

########################### own transform tukey function

transformTukey2 <- function(df, start = -10, end = 10, int = 0.25, model_style) {
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
        { capture.output(model <- test_models(transformed_df, model_style = model_style)) }
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

transform_data <- function(df, transformation = c("gauss", "tukey", "tukey2", "log", "logit"), 
                           type = c("h", "hh", "s"), model_style) {
  if (transformation == "gauss") {
    transformed <- Gaussianize(df$response, type = type, return.tau.mat = T)
    
  } else if (transformation == "tukey") {
    
    transformed1 <- transformTukey(df$response)
    transformed2 <- transformTukey(df$response, returnLambda = T)
    transformed <- list(transformed1, transformed2)
  } else if (transformation == "tukey2") { # leads to choosing the worst model
    transformed <- transformTukey2(df, model_style = model_style)
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

test_transformations <- function(df, model_style) {
  result_frame <- data.frame(transformation = c("none", "tukey", "tukey2", "log", "logit", "gauss", "gauss", "gauss"),
                             type = c(NA, NA, NA, NA, NA, "h", "hh", "s"),
                             W = 0)
  for (i in 1:nrow(result_frame)) {
    transformation <- result_frame$transformation[i]
    type <- result_frame$type[i]
    if (transformation == "none") {
      transformed_df <- df
    } else {
      transformed_df <- transform_data(df, transformation = transformation, 
                                       type = type, model_style = model_style)[[1]]
    }
    
    if (any(is.infinite(transformed_df$response)) == FALSE & any(is.nan(transformed_df$response)) == FALSE) {
      x <- suppressMessages(suppressWarnings(
        { capture.output(model <- test_models(transformed_df, model_style = model_style)) }
      ))
      test_statistic <- shapiro.test(resid(model))
      if (deviance(model) < sqrt(.Machine$double.eps)) { # if precision gets to low (similar to test in car::Anova)
        result_frame$W[i] <- 0 # do not use this transformation
      } else {
        result_frame$W[i] <- test_statistic$statistic
      }
    }
    
  }
  best_transformation <- result_frame$transformation[which.max(result_frame$W)] # returning only first 
  best_type <- result_frame$type[which.max(result_frame$W)]
  return(list(best_transformation = best_transformation, best_type = best_type))
}


combined_comparison <- function(df,transformation = "none", 
                                type = c("h", "hh", "s"), model_style = "klein") { #, correction_filter
  if (transformation == "test") {
    best_transformation = test_transformations(df, model_style = model_style)
    transformation <- best_transformation$best_transformation
    type <- best_transformation$best_type
    print(paste("Best transformation:", transformation))
  }
  
  if (transformation == "none") {
    transformation_factor = NA # compatible for output
    
    model <- test_models(df, model_style = model_style)
    
    n <- create_n_table(df)
    
    shapiro <- evaluate_model(model)
    print(shapiro)
    
    ano <- do_anova(model)
    
    cld <- pairwise_comparisons(model)
    
    cld2 <- cld # for compatibility with transformation (original means)
    
  } else {
    
    transformed <- transform_data(df, transformation = transformation, type = type, model_style = model_style)
    transformed_df <- transformed[[1]]
    transformation_factor <- transformed[[2]]
    
    model <- test_models(transformed_df, model_style = model_style)
    
    n <- create_n_table(df)
    
    shapiro <- evaluate_model(model)
    print(shapiro)
    
    ano <- do_anova(model)
    
    cld <- pairwise_comparisons(model)
    
    cld <- backtransform_data(cld, transformation = transformation, transformation_factor = transformation_factor)
    
    model_untransformed <- test_models(df, model_style = model_style)
    cld2 <- pairwise_comparisons(model_untransformed)
    
  }
  
  if (ano[which(rownames(ano) == "diet"), which(colnames(ano) == "Pr(>F)")] > .05) {
    cld <- cld %>%
      mutate(.group = "")
  } # if not significant, remove contrasts (if calculated)
  
  return(list(ano = ano, cld = cld, cld2 = cld2, n = n, 
              transformation = c(transformation, transformation_factor, type),
              shapiro = shapiro))
}

create_results_table <- function(input_df, comparison_object, response_name = "response", digits = 1) {
  P <- comparison_object$ano[which(rownames(comparison_object$ano) == "diet"), 
                             which(colnames(comparison_object$ano) == "Pr(>F)")]
  Print <- str_remove(round_p_value(P), "= ")
  n <- comparison_object$n
  cld <- comparison_object$cld
  cld2 <- comparison_object$cld2 %>%
    inner_join(n, by = "diet")
  pSEM <- paste(format(round(sum(cld2$SE*(cld2$n-1)) / (sum(n$n)-length(n$n)),digits+1), nsmall = digits+1))
  table <- dplyr::select(cld2, diet, emmean, SE) %>%
    inner_join(dplyr::select(cld, diet, .group), by = "diet") %>%
    arrange(diet) %>%
    inner_join(distinct(dplyr::select(input_df, diet)), by = "diet") %>%
    mutate(!!sym(response_name) := paste0(format(round(emmean, digits), nsmall = digits),
                                          if (length(unique(cld$.group))==1) {""} else {str_trim(.group)})) %>%
    #dplyr::select(diet = description, !!sym(response_name)) %>%
    dplyr::select(diet, !!sym(response_name)) %>%
    add_row(diet = "Pooled SEM", !!sym(response_name) := pSEM) %>%
    add_row(diet = "P-value", !!sym(response_name) := Print)
  return(table)
}

# plotting function

plot_pairwise <- function(filtered_df, output, selected_response, save = F) { #, correction_filter
  cld <- output$cld
  P <- output$ano[which(rownames(output$ano) == "diet"), which(colnames(output$ano) == "Pr(>F)")]
  Print <- round_p_value(P)
  p <- ggplot(cld, aes(x = diet)) +
    geom_boxplot(data = filtered_df, aes(y = response), width = 0.2, position = position_nudge(x = -0.1)) +
    geom_point(data = filtered_df, aes(y = response, shape = animal, color = period), 
               size = 2, position = position_nudge(x = -0.1), stroke = 2) +
    #geom_point(aes(y = emmean), size = 2, position = position_nudge(x = 0)) +
    #geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0, lwd = 1, position = position_nudge(x = 0)) +
    {if(P < 0.05)geom_text(aes(y = emmean, label = str_trim(.group)), size = 5, position = position_nudge(x = 0.2))}+
    labs(x = "diet", y = "value", title = paste(selected_response)) + #, correction_filter
    annotate("label", label = paste("diet\n P", Print), x = 4.25, y = max(filtered_df$response), size = 5) +
    #scale_color_manual(values = c("black", "grey20", "grey40", "grey60")) +
    scale_color_manual(values = c("red", "blue", "green", "yellow")) +
    scale_shape_manual(values = c(16,17,15,18,21,24,22,23))
  print(p)
}

create_results_plot <- function(input_df, comparison_object, response_name = "response", y_axis = "value") {
  P <- comparison_object$ano[which(rownames(comparison_object$ano) == "diet"), 
                             which(colnames(comparison_object$ano) == "Pr(>F)")]
  Print <- round_p_value(P)
  n <- comparison_object$n
  cld <- comparison_object$cld
  cld2 <- comparison_object$cld2 %>%
    inner_join(n, by = "diet")
  pSEM <- sum(cld2$SE*(cld2$n-1)) / (sum(n$n)-length(n$n))
  table <- dplyr::select(cld2, diet, emmean, SE) %>%
    inner_join(dplyr::select(cld, diet, .group), by = "diet") %>%
    arrange(diet) %>%
    inner_join(distinct(dplyr::select(input_df, diet)), by = "diet") %>%
    mutate(lower = emmean - SE, upper = emmean + SE)
  p <- ggplot(table, aes(x = diet)) +
    geom_bar(aes(y = emmean), stat = "identity", width = 0.4) +
    geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.1, size = 1) +
    annotate("label", label = paste("F-Test \n P-value", Print), x = nrow(n), y = max(table$upper)/10, size = 6, fill = "snow2") +
    {if(P < 0.05)geom_text(aes(y = (2*emmean + SE)/2, label = str_trim(.group)), size = 5, position = position_nudge(x = 0.1))}+
    scale_y_continuous(limits = c(min(0, table$lower),max(table$upper)), expand = expansion(mult = c(0, .1))) +
    labs(x = "diet", y = y_axis, title = response_name)
  print(p)
}

save_results_table <- function(table, title, name) {
  gt <- gt(table, rowname_col = "diet") %>%
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
    cols_align(align="center", columns = -diet) %>%
    opt_table_font(font = google_font("Merriweather"))
  gtsave(gt, filename = paste0(name, ".png"), "plots")
  return(gt)
}

# to split eg metaproteomic data into ileum and feces

filter_ileum <- function(input_df) {
  output <- input_df %>%
    filter(str_detect(sampleid, "^1"))
  return(output)
}

filter_faeces <- function(input_df) {
  output <- input_df %>%
    filter(str_detect(sampleid, "^3"))
  return(output)
}

filter_ferment <- function(input_df) {
  output <- input_df %>%
    filter(str_detect(sampleid, "^4"))
  return(output)
}
