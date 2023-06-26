## Code for: White (in review) Deceptive pollinator lures benefit from physical 
## and perceptual proximity to flowers

## --- Clear out --- ##
  rm(list = ls())
  
## --- Libraries --- ##
  library(patchwork)
  library(performance)
  library(pavo)
  library(tidyverse)
  library(interplot)
  library(MuMIn)
  library(easystats)
  library(emmeans)
  
## --- Data --- ##
  caps_obs <- read.csv('../data/dat_captures_tidy.csv')
  caps_exp <- read.csv('../data/dat_exp_tidy.csv')
  spec_flowers <- procspec(as.rspec(read.csv('../data/spec_flowers.csv'), lim = c(300, 700)), opt = 'smooth', fixneg = 'zero')
  spec_spiders <- procspec(as.rspec(read.csv('../data/spec_spiders.csv'), lim = c(300, 700)), opt = 'smooth', fixneg = 'zero')
  spec_all <- merge(spec_spiders, spec_flowers)
  
## --- Visual modelling --- ##
  
  # Create a new spider/flower ID variable for ease of merging with
  # colour-distance data
  caps_obs$comb_id <- paste0(caps_obs$spider_id, "-", caps_obs$flower_id)
  vmod <-  vismodel(merge(spec_spiders, spec_flowers), visual = 'drosophila')
  vdist <- coldist(colspace(vmod, space = 'categorical'))  
  
  # Estimate colour distances in Troje (1993) space using Drosophila 
  # visual phenotype
  vdist <- 
    spec_all |> 
    vismodel(visual = 'drosophila') |> 
    colspace(space = 'categorical') |> 
    coldist() |> 
    mutate(comb_id = paste0(patch1, "-", patch2))
  
  # Merge into primary dataset  
  caps_obs <- left_join(caps_obs, dplyr::select(vdist, comb_id, dS))
  
## --- Statistical modelling --- ##

  # Examine distributions for observational data
  hist(caps_obs$dist_cm)
  quantile(caps_obs$dist_cm, probs = seq(0.05, 0.95, 0.05))  # physical distance
  
  hist(caps_obs$dS)
  quantile(caps_obs$dS, probs = seq(0.05, 0.95, 0.05))  # colour distance
  
## Observational  
  # Full model
  obs_mod_full <- glm(caps_hr ~ dS + poly(dS, 2) + poly(dS, 3) + dist_cm + 
                      dS:dist_cm + poly(dS, 2):dist_cm + poly(dS, 3):dist_cm, 
                      data = caps_obs,
                      na.action = 'na.fail')
  # Model reduction
  dredge(obs_mod_full)
  
  # Leading model (rename vars for convenient display)
  obs_mod <- glm(caps_hr ~ dS * dist_cm, 
                 data = caps_obs,
                 na.action = 'na.fail')
  
  # Inspect diagnostic plots
  check_model(obs_mod)
  
  # Model summary & fit
  summary(obs_mod)
  r2(obs_mod)
  
  # Tidy model results table 
  model_parameters(obs_mod)

## Experimental
  # Examine distribution
  hist(caps_exp$intercept_rate)
  
  # Model
  exp_mod <- glm(intercept_rate ~ proximity * colour, data = caps_exp)
  
  # Inspect diagnostic plots
  check_model(exp_mod)
  
  # Model summary & fit
  summary(exp_mod)
  r2(exp_mod)
  
  # Post-hoc test
  caps_exp$treatment <- factor(caps_exp$treatment)
  summary(glht(glm(intercept_rate ~ treatment, data = caps_exp), linfct = mcp(treatment = "Tukey")))
  pairs(emmeans(glm(intercept_rate ~ treatment, data = caps_exp), "treatment"), adjust = 'none')
  
  # Tidy model results table 
  model_parameters(exp_mod)
  ggplot(caps_exp, aes(x = treatment, y = intercept_rate)) +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.4) +
    ylab('Interceptions / hour') +
    xlab('Treatment') +
    scale_x_discrete(labels = c("Far & colour-match", "Far & colour-mismatch", "Near & colour-match", 'Near & colour-mismatch')) +
    theme_classic()
  
## --- Plots --- ##
  
  # Spectra
  png('../figs/fig_specs.png', width = 26, height = 14, units = 'cm', res = 300)
  par(mfrow = c(1, 2))
    plot(spec_spiders, 
         ylim = c(0, 100))
    text(305, 100, '(a)', cex = 1.4)
    plot(spec_flowers, 
         ylab = '',
         ylim = c(0, 100))
    text(305, 100, '(b)', cex = 1.4)
  dev.off()
  
## Observational
  
  # Colour-distance distribution
  (plot_coldist <- 
      ggplot(caps_obs, aes(x = dS)) +
      geom_histogram(color = "black", fill = 'grey') +
      geom_vline(aes(xintercept = 0.02), linetype = "dashed", size = 1, alpha = 0.5) +
      geom_vline(aes(xintercept = 0.26), linetype = "dashed", size = 1, alpha = 0.5) +
      xlab('Colour distance (dS)') +
      ylab('Count') +
      ylim(c(0, 7)) +
      theme_classic())
  
  # Physical distance distribution
  (plot_physdist <- 
      ggplot(caps_obs, aes(x = dist_cm)) +
      geom_histogram(color = "black", fill = 'grey', bins = 20) +
      geom_vline(aes(xintercept = 60), linetype = "dashed", size = 1, alpha = 0.5) +
      geom_vline(aes(xintercept = 430), linetype = "dashed", size = 1, alpha = 0.5) +
      ylim(c(0, 7)) +
      xlab('Physical distance (cm)') +
      ylab(' ') +
      theme_classic())
  
  # Combine into a single graphic
  plot_coldist + plot_physdist
  
  # Save
  ggsave('../figs/fig_obsdists.tiff', height = 5)
  ggsave('../figs/fig_obsdists.png', height = 5)
  
  # Interception rate ~ colour distance
  (plot_ds <- 
      ggplot(caps_obs, aes(x = dS, y = caps_hr)) +
      geom_point() +
      geom_smooth(method = 'lm', linewidth = 0.5, color = 'black') +
      ylim(0, 7) +
      xlab('Colour distance (dS)') +
      ylab('Interceptions / hour') +
      theme_classic())
  
  # Interception rate ~ physical distance
  (plot_dist <- 
      ggplot(caps_obs, aes(x = dist_cm, y = caps_hr)) +
      geom_point() +
      geom_smooth(method = 'lm', linewidth = 0.5, color = 'black') +
      xlab('Physical distance (cm)') + 
      ylim(0, 7) +
      ylab(' ') +
      theme_classic())
  
  # Conditional coefficients of colour distance ~ physical distance 
  # (Interaction plot)
  a <- interplot(m = obs_mod, var1 = "dS", var2 = "dist_cm", plot = FALSE)
  (plot_int <- 
    interplot(m = obs_mod, var1 = "dS", var2 = "dist_cm", plot = TRUE) +
    #geom_point(data = a, aes(x = dist_cm, y = coef)) +
    xlab('Physial distance (cm)') +
    ylab('Coefficient for \n colour distance (dS)') +
    theme_classic())
  
  # Combine into a single graphic
  plot_int / (plot_ds | plot_dist)
  
  # Save
  ggsave('../figs/fig_models.tiff', height = 9)
  ggsave('../figs/fig_models.png', height = 9)
  
## Experimental
  # Treatment effects
  ggplot(caps_exp, aes(x = treatment, y = intercept_rate)) +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.4) +
    ylab('Interceptions / hour') +
    xlab('Treatment') +
    scale_x_discrete(labels = c("Far & colour-match", "Far & colour-mismatch", "Near & colour-match", 'Near & colour-mismatch')) +
    theme_classic() +
    annotate('text', x = 1, y = 6.2, label = 'c', size = 5) +
    annotate('text', x = 2, y = 6.2, label = 'c', size = 5) +
    annotate('text', x = 3, y = 6.2, label = 'a', size = 5) +
    annotate('text', x = 4, y = 6.2, label = 'b', size = 5)
  
  # Save
  ggsave('../figs/fig_exp.tiff')
  ggsave('../figs/fig_exp.png')
  
  
  