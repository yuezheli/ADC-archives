# date: 1/15/2026 
# author: Yuezhe Li 
# purpose of this code: to generate figures for the repo

library(dplyr)
library(tidyverse)
library(readr)
library(ggplot2)
library(here)

# Setup -------------------------------------------------------------------
figDir <- here::here("deliv", "figure")

options(mrg.script = "script/pk-summary-ints-pl.R", 
        pmplots.path.type = "proj",
        mrggsave.dir = figDir,
        mrggsave.dev = "pdf")

setwd(here::here())

theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))

# PL -------------------------------------------------------------------
tdm1_adc_ints <- read_csv(here::here("data/t-dm1-tissue-ints.csv"))
tdm1_tissue_ints <- read_csv(here::here("data/t-dm1-pl-tissue-ints.csv"))
tdm1_tissue_endo <- read_csv(here::here("data/t-dm1-pl-tissue-endo.csv"))

tdxd_adc_ints <- read_csv(here::here("data/t-dxd-tissue-ints.csv"))
tdxd_tissue_ints <- read_csv(here::here("data/t-dxd-pl-tissue-ints.csv"))
tdxd_tissue_endo <- read_csv(here::here("data/t-dxd-pl-tissue-endo.csv"))

cm_adc_ints <- read_csv(here::here("data/CanAg-dm1-tissue-ints.csv"))
cm_tissue_ints <- read_csv(here::here("data/CanAg-dm1-pl-tissue-ints.csv"))
cm_tissue_endo <- read_csv(here::here("data/CanAg-dm1-pl-tissue-endo.csv"))

# visualization -------------------------------------------------------------------

# ADC PK comparison, cantuzumab mersantine vs T-DM1 
plt_adc_cm_tdm1 <- ggplot( 
  data = rbind(tdm1_adc_ints %>% mutate(drug = "T-DM1"), 
               cm_adc_ints %>% mutate(drug = "Cantuzumab Mertansine")) %>% 
    rename(`large intestine` = "la_int", 
           `small intestine` = "sm_int", 
           `aqueous humor` = ah, 
           `vitreous humor` = vh, 
           `iris-ciliary body` = icb, 
           `bone marrow` = marrow, 
           `lymph nodes` = ln) %>%
    filter(Dose %in% c(0.3, 0.6, 1.2, 1.8, 2.4, 3.6, 4.8)) %>%
    pivot_longer(cols = -c(time_hr, Dose, drug), names_to = c("organ"), values_to = c("conc")), 
  aes(x = time_hr/24, y = conc, col = as.factor(Dose), linetype = as.factor(drug) )) +
  geom_line() + 
  labs(x="Time (hr)", y = "tissue interstitial ADC concentration (uM)", col = "Dose (mg/kg)", linetype = "Drug") +
  scale_y_continuous(trans = 'log10', limits = c(1E-7, 1E1)) +
  scale_x_continuous(breaks = c(0, 21, 42, 63, 84), limits = c(0, 84)) + 
  facet_wrap(~organ)

# ADC PK prediction, T-DM1
plt_adc_tdm1 <- ggplot( 
  data = tdm1_adc_ints %>% 
    rename(`large intestine` = "la_int", 
           `small intestine` = "sm_int", 
           `aqueous humor` = ah, 
           `vitreous humor` = vh, 
           `iris-ciliary body` = icb, 
           `bone marrow` = marrow, 
           `lymph nodes` = ln) %>%
    select( -c("cornea" )) %>%
    filter(Dose %in% c(0.3, 0.6, 1.2, 1.8, 2.4, 3.6, 4.8)) %>%
    pivot_longer(cols = -c(time_hr, Dose), names_to = c("organ"), values_to = c("conc")), 
  aes(x = time_hr/24, y = conc, col = as.factor(Dose) )) +
  geom_line() + 
  labs(x="Time (hr)", y = "tissue interstitial ADC concentration (uM)", col = "Dose (mg/kg)") +
  scale_y_continuous(trans = 'log10', limits = c(1E-7, 1E1)) +
  scale_x_continuous(breaks = c(0, 21, 42, 63, 84), limits = c(0, 84)) + 
  facet_wrap(~organ)

# ADC PK prediction, T-Dxd
plt_adc_tdxd <- ggplot( 
  data = tdxd_adc_ints %>% 
    rename(`large intestine` = "la_int", 
           `small intestine` = "sm_int", 
           `aqueous humor` = ah, 
           `vitreous humor` = vh, 
           `iris-ciliary body` = icb, 
           `bone marrow` = marrow, 
           `lymph nodes` = ln) %>%
    select( -c("cornea" )) %>%
    filter(Dose %in% c(0.8, 1.6, 3.2, 5.4, 6.4, 8.0)) %>%
    pivot_longer(cols = -c(time_hr, Dose), names_to = c("organ"), values_to = c("conc")), 
  aes(x = time_hr/24, y = conc, col = as.factor(Dose) )) +
  geom_line() + 
  labs(x="Time (hr)", y = "tissue interstitial ADC concentration (uM)", col = "Dose (mg/kg)") +
  scale_y_continuous(trans = 'log10', limits = c(1E-7, 1E1)) +
  scale_x_continuous(breaks = c(0, 21, 42, 63, 84), limits = c(0, 84)) + 
  facet_wrap(~organ)


# PL PK comparison, cantuzumab mersantine vs T-DM1 
plt_endo_pl_cm_tdm1 <- ggplot( 
  data = rbind(tdm1_tissue_endo %>% mutate(drug = "T-DM1"), 
               cm_tissue_endo %>% mutate(drug = "Cantuzumab Mertansine")) %>% 
    #select( -c("la_int", "sm_int", "ah", "vh", "icb") ) %>% 
    filter(Dose %in% c(0.3, 0.6, 1.2, 1.8, 2.4, 3.6, 4.8)) %>%
    pivot_longer(cols = -c(time_hr, Dose, drug), names_to = c("organ"), values_to = c("conc")), 
  aes(x = time_hr/24, y = conc, col = as.factor(Dose), linetype = as.factor(drug) )) +
  geom_line() + 
  labs(x="Time (hr)", y = "PL endothelial concentration (uM)", col = "Dose (mg/kg)", linetype = "Drug") +
  scale_y_continuous(trans = 'log10', limits = c(1E-7, 1E1)) +
  scale_x_continuous(breaks = c(0, 21, 42, 63, 84), limits = c(0, 84)) + 
  facet_wrap(~organ)

# PL PK comparison, T-DM1 vs T-Dxd
plt_endo_pl_tdm1_tdxd <- ggplot( 
  data = rbind(tdm1_tissue_endo %>% mutate(drug = "T-DM1"), 
               tdxd_tissue_endo %>% mutate(drug = "T-Dxd")) %>% 
    rename(`large intestine` = "la_int", 
           `small intestine` = "sm_int", 
           `aqueous humor` = ah, 
           `vitreous humor` = vh, 
           `iris-ciliary body` = icb, 
           `bone marrow` = marrow) %>%
    select( -c("cornea" )) %>%
    filter(Dose %in% c(3.6, 5.4)) %>%
    pivot_longer(cols = -c(time_hr, Dose, drug), names_to = c("organ"), values_to = c("conc")), 
  aes(x = time_hr/24, y = conc, col = as.factor(Dose), linetype = as.factor(drug) )) +
  geom_line() + 
  labs(x="Time (hr)", y = "PL endothelial concentration (uM)", col = "Dose (mg/kg)", linetype = "Drug") +
  scale_y_continuous(trans = 'log10', limits = c(1E-5, 1E1)) +
  scale_x_continuous(breaks = c(0, 21, 42, 63, 84), limits = c(0, 84)) + 
  facet_wrap(~organ)


# PL PK prediction, T-Dxd
plt_pl_tdxd <- ggplot( 
  data = rbind(tdxd_tissue_ints %>% mutate(Type = "Interstitium"), 
               tdxd_tissue_endo %>% mutate(Type = "Endothelial")
               ) %>% 
    rename(`large intestine` = "la_int", 
           `small intestine` = "sm_int", 
           `iris-ciliary body` = icb, 
           `bone marrow` = marrow
          ) %>%
    select( -c("cornea", "ah", "vh", "plasma" )) %>%
    filter(Dose %in% c(0.8, 1.6, 3.2, 5.4, 6.4, 8.0)) %>%
    pivot_longer(cols = -c(time_hr, Dose, Type), names_to = c("organ"), values_to = c("conc")), 
  aes(x = time_hr/24, y = conc, col = as.factor(Dose), linetype = as.factor(Type) )) +
  geom_line() + 
  labs(x="Time (hr)", y = "tissue payload concentration (uM)", col = "Dose (mg/kg)") +
  scale_y_continuous(trans = 'log10', limits = c(1E-9, 1E-1), breaks = c(1E-9, 1E-7, 1E-5, 1E-3, 1E-1)) +
  scale_x_continuous(breaks = c(0, 21, 42, 63, 84), limits = c(0, 84)) + 
  facet_wrap(~organ)

# save figrues -------------------------------------------------------------------
ggsave(here::here("deliv/figure/ADC-T-DM1-CM.png"), plt_adc_cm_tdm1, width = 12, height = 6, units = c("in"), dpi = 300)
ggsave(here::here("deliv/figure/PL-T-DM1-CM.png"), plt_endo_pl_cm_tdm1, width = 12, height = 6, units = c("in"), dpi = 300)
ggsave(here::here("deliv/figure/ADC-T-DM1.png"), plt_adc_tdm1, width = 10, height = 6, units = c("in"), dpi = 300)
ggsave(here::here("deliv/figure/ADC-T-DXD.png"), plt_adc_tdxd, width = 10, height = 6, units = c("in"), dpi = 300)
ggsave(here::here("deliv/figure/PL-T-DM1-T-DXD.png"), plt_endo_pl_tdm1_tdxd, width = 10, height = 6, units = c("in"), dpi = 300)
ggsave(here::here("deliv/figure/PL-T-DXD-INTS-ENDO.png"), plt_pl_tdxd, width = 10, height = 6, units = c("in"), dpi = 300)

