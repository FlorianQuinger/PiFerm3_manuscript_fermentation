library(here)

source("0_general_functions.R")
source("20_pretrial_functions.R")

# load data

ph <- read_tsv("data/2_pretrial1.txt") %>%
  mutate(inoculation = as.character(inoculation),
         dm = as.character(dm),
         sampleid = as.character(sampleid))

ph_long <- ph %>%
  pivot_longer(-c(sampleid, dm, inoculation), values_to = "ph", names_to="time") %>%
  mutate(time = as.numeric(str_remove(time, "h")))

# to test different transformations


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
    if (any(is.infinite(transformed_df$response)) == FALSE & any(is.nan(transformed_df$response)) == FALSE) { 
      model <- aov(transformed_df$response ~ dm + inoculation+ dm:inoculation, data = df)
      
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


# plot by dm 

create_plots_for_dm <- function(dm_filter) {
  p <- ph_long %>%
    filter(dm == dm_filter) %>%
    ggplot(aes(x = time, y = ph, color = inoculation)) +
    geom_point(aes(shape = inoculation), size = 3) +
    geom_line(aes(group = sampleid), lwd= 1) +
    scale_color_manual(values = colors) +
    labs(y = "pH", x = "time (h)", shape = "conc. (µl)", color = "conc. (µl)", 
         title = paste("Pea ratio", dm_filter)) +
    ylim(3.5, 6.5)
 return(p) 
}

for (i in c("20", "30", "40", "50")) {
  create_plots_for_dm(dm = i)
  save_big(paste0("21_ph_curve_", i))
}

ph_long %>%
  ggplot(aes(x = time, y = ph, color = inoculation)) +
  geom_point(aes(shape = inoculation), size = 3) +
  geom_line(aes(group = sampleid, linetype = dm), lwd= 1) +
  scale_color_manual(values = colors) +
  labs(y = "pH", x = "time (h)", shape = "conc. (µl)", color = "conc. (µl)") +
  ylim(3.5, 6.5)

# statistial comparison of factors

# at t0

ph0 <- ph_long %>%
  filter(time == 0)

model_0 <- aov(ph ~ dm+ inoculation + dm:inoculation, data = ph0)
plot(model_0)
shapiro.test(resid(model_0))
Anova(model_0, type = "III")

cld_0dm <- cld(emmeans(model_0, ~dm), 
               Letters = letters, details = F, reversed = T)
cld_0inoc <- cld(emmeans(model_0, ~inoculation), 
                 Letters = letters, details = F, reversed = T)
cld_0 <- cld(emmeans(model_0, ~ dm + inoculation + dm:inoculation), 
             Letters = letters, details = F, reversed = T)

# at t2

ph2 <- ph_long %>%
  filter(time == 2)

model_2 <- aov(ph ~ dm + inoculation + dm:inoculation, data = ph2)
plot(model_2)
shapiro.test(resid(model_2))
Anova(model_2, type = "III")

cld_2dm <- cld(emmeans(model_2, ~dm), 
               Letters = letters, details = F, reversed = T)
cld_2inoc <- cld(emmeans(model_2, ~inoculation), 
                 Letters = letters, details = F, reversed = T)
cld_2 <- cld(emmeans(model_2, ~ dm + inoculation + dm:inoculation), 
             Letters = letters, details = F, reversed = T)

create_results_plot2(cld_2, response_name = "pH 2h", y_axis = "pH")
save_big("21_ph_2h")

# at t4

ph4 <- ph_long %>%
  filter(time == 4)

model_4 <- aov(ph ~ dm + inoculation+ dm:inoculation, data = ph4)
plot(model_4)
shapiro.test(resid(model_4))
model_4t <- aov(transformTukey(ph4$ph) ~ dm + inoculation+ dm:inoculation, data = ph4)
model_4t <- aov(Gaussianize(ph4$ph, type = "s") ~ dm + inoculation+ dm:inoculation, data = ph4)
model_4t <- aov((10^-ph4$ph)*10e5 ~ dm + inoculation+ dm:inoculation, data = ph4)
plot(model_4t)
shapiro.test(resid(model_4t))
Anova(model_4t, type = "III")

cld_4dm <- cld(emmeans(model_4, ~dm), 
               Letters = letters, details = F, reversed = T)
cld_4inoc <- cld(emmeans(model_4, ~inoculation), 
                 Letters = letters, details = F, reversed = T)
cld_4 <- cld(emmeans(model_4, ~ dm + inoculation + dm:inoculation), 
             Letters = letters, details = F, reversed = T)

cld_4tdm <- cld(emmeans(model_4t, ~dm), 
               Letters = letters, details = F, reversed = T)
cld_4tinoc <- cld(emmeans(model_4t, ~inoculation), 
                 Letters = letters, details = F, reversed = T)

create_results_plot2(cld_4, response_name = "pH 4h", y_axis = "pH")
save_big("21_ph_4h")

# at t6

ph6 <- ph_long %>%
  filter(time == 6)

model_6 <- aov(ph ~ dm + inoculation+ dm:inoculation, data = ph6)
plot(model_6)
shapiro.test(resid(model_6))
Anova(model_6, type = "III")

cld_6dm <- cld(emmeans(model_6, ~dm), 
               Letters = letters, details = F, reversed = T)
cld_6inoc <- cld(emmeans(model_6, ~inoculation), 
                 Letters = letters, details = F, reversed = T)
cld_6 <- cld(emmeans(model_6, ~ dm + inoculation + dm:inoculation), 
             Letters = letters, details = F, reversed = T)

create_results_plot2(cld_6, response_name = "pH 6h", y_axis = "pH")
save_big("21_ph_6h")

# at t8

ph8 <- ph_long %>%
  filter(time == 8)

model_8 <- aov(ph ~ dm + inoculation+ dm:inoculation, data = ph8)
plot(model_8)
shapiro.test(resid(model_8))
Anova(model_8, type = "III")

# cld_8dm <- cld(emmeans(model_8, ~dm), 
#                Letters = letters, details = F, reversed = T)
# cld_8inoc <- cld(emmeans(model_8, ~inoculation), 
#                  Letters = letters, details = F, reversed = T)
cld_8 <- cld(emmeans(model_8, ~ dm + inoculation + dm:inoculation), 
             Letters = letters, details = F, reversed = T)

create_results_plot2(cld_8, response_name = "pH 8h", y_axis = "pH")
save_big("21_ph_8h")

# at t12

ph12 <- ph_long %>%
  filter(time == 12)

model_12 <- aov(ph ~ dm + inoculation+ dm:inoculation, data = ph12)
plot(model_12)
shapiro.test(resid(model_12))
model_12t <- aov(transformTukey(ph12$ph) ~ dm + inoculation+ dm:inoculation, data = ph12)
model_12t <- aov(Gaussianize(ph12$ph, type = "hh") ~ dm + inoculation+ dm:inoculation, data = ph12)
model_12t <- aov((10^-ph12$ph)*10e5 ~ dm + inoculation+ dm:inoculation, data = ph12)
plot(model_12t)
shapiro.test(resid(model_12t))
Anova(model_12t, type = "III")

# cld_12dm <- cld(emmeans(model_12, ~dm), 
#                Letters = letters, details = F, reversed = T)
# cld_12inoc <- cld(emmeans(model_12, ~inoculation), 
#                  Letters = letters, details = F, reversed = T)
cld_12 <- cld(emmeans(model_12, ~ dm + inoculation + dm:inoculation), 
             Letters = letters, details = F, reversed = T)

cld_12t <- cld(emmeans(model_12t, ~ dm + inoculation + dm:inoculation), 
              Letters = letters, details = F, reversed = T)

create_results_plot2(cld_12t, response_name = "pH 12h", y_axis = "pH")
save_big("21_ph_12h")

# at t24

ph24 <- ph_long %>%
  filter(time == 24)

model_24 <- aov(ph ~ dm + inoculation+ dm:inoculation, data = ph24)
plot(model_24)
shapiro.test(resid(model_24))
model_24t <- aov(transformTukey(ph24$ph) ~ dm + inoculation+ dm:inoculation, data = ph24)
model_24t <- aov(Gaussianize(ph24$ph, type = "s") ~ dm + inoculation+ dm:inoculation, data = ph24)
#model_24t <- aov((10^-ph24$ph)*10e5 ~ dm + inoculation+ dm:inoculation, data = ph24)
plot(model_24t)
shapiro.test(resid(model_24t))
Anova(model_24t, type = "III")

cld_24dm <- cld(emmeans(model_24, ~dm), 
                Letters = letters, details = F, reversed = T)
cld_24inoc <- cld(emmeans(model_24, ~inoculation), 
                  Letters = letters, details = F, reversed = T)
cld_24 <- cld(emmeans(model_24, ~ dm + inoculation + dm:inoculation), 
              Letters = letters, details = F, reversed = T)

cld_24tdm <- cld(emmeans(model_24t, ~dm), 
                Letters = letters, details = F, reversed = T)
cld_24tinoc <- cld(emmeans(model_24t, ~inoculation), 
                  Letters = letters, details = F, reversed = T)

create_results_plot2(cld_24, response_name = "pH 24h", y_axis = "pH")
save_big("21_ph_24h")

# at t36

ph36 <- ph_long %>%
  filter(time == 36)

model_36 <- aov(ph ~ dm + inoculation+ dm:inoculation, data = ph36)
plot(model_36)
shapiro.test(resid(model_36))
model_36t <- aov(transformTukey(ph36$ph) ~ dm + inoculation+ dm:inoculation, data = ph36)
model_36t <- aov(Gaussianize(ph36$ph, type = "h") ~ dm + inoculation+ dm:inoculation, data = ph36)
#model_36t <- aov((10^-ph36$ph)*10e5 ~ dm + inoculation+ dm:inoculation, data = ph36)
#model_36t <- aov((ph36$ph^2) ~ dm + inoculation+ dm:inoculation, data = ph36)
#model_36t <- aov(transformTukey2(dplyr::rename(ph36, response = ph))[[1]] ~ dm + inoculation+ dm:inoculation, data = ph36)
plot(model_36t)
shapiro.test(resid(model_36t))
Anova(model_36t, type = "III")

cld_36dm <- cld(emmeans(model_36, ~dm), 
                Letters = letters, details = F, reversed = T)
cld_36inoc <- cld(emmeans(model_36, ~inoculation), 
                  Letters = letters, details = F, reversed = T)
cld_36 <- cld(emmeans(model_36, ~ dm + inoculation + dm:inoculation), 
              Letters = letters, details = F, reversed = T)

cld_36tdm <- cld(emmeans(model_36t, ~dm), 
                 Letters = letters, details = F, reversed = T)
cld_36tinoc <- cld(emmeans(model_36t, ~inoculation), 
                   Letters = letters, details = F, reversed = T)

create_results_plot2(cld_36, response_name = "pH 36h", y_axis = "pH")
save_big("21_ph_36h")


# at t48

ph48 <- ph_long %>%
  filter(time == 48)

model_48 <- aov(ph ~ dm + inoculation+ dm:inoculation, data = ph48)
plot(model_48)
shapiro.test(resid(model_48))
model_48t <- aov(transformTukey(ph48$ph) ~ dm + inoculation+ dm:inoculation, data = ph48)
model_48t <- aov(Gaussianize(ph48$ph, type = "hh") ~ dm + inoculation+ dm:inoculation, data = ph48)
#model_48t <- aov((10^-ph48$ph)*10e5 ~ dm + inoculation+ dm:inoculation, data = ph48)
#model_48t <- aov((ph48$ph^2) ~ dm + inoculation+ dm:inoculation, data = ph48)
#model_48t <- aov(transformTukey2(dplyr::rename(ph48, response = ph))[[1]] ~ dm + inoculation+ dm:inoculation, data = ph48)
plot(model_48t)
shapiro.test(resid(model_48t))
Anova(model_48t, type = "III")

cld_48dm <- cld(emmeans(model_48, ~dm), 
                Letters = letters, details = F, reversed = T)
cld_48inoc <- cld(emmeans(model_48, ~inoculation), 
                  Letters = letters, details = F, reversed = T)
cld_48 <- cld(emmeans(model_48, ~ dm + inoculation + dm:inoculation), 
              Letters = letters, details = F, reversed = T)

cld_48tdm <- cld(emmeans(model_48t, ~dm), 
                Letters = letters, details = F, reversed = T)

create_results_plot2(cld_48, response_name = "pH 48h", y_axis = "pH")
save_big("21_ph_48h")
