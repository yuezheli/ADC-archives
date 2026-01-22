# date: 7/2/2025 
# author: Yuezhe Li 
# purpose of this code: analysis on payload dynamics in liver endothelial cells 

rm(list=ls()) 
gc()

library(here)
library(data.table)
library(dplyr)
library(tidyverse)
library(pmtables)
library(ggplot2)

free_he_pl <- function(df, timespan_d = 84){
  df_hep <- df %>% 
    filter(time_hr <= timespan_d * 24) %>%
    group_by(dose_mgkg) %>% 
    mutate(
      Cmax_he_nM = sig(max(pl_liver_uM) * 1E3 , digits = 3) , 
      Cavg_he_nM = sig(mean(pl_liver_uM) * 1E3, digits = 3) 
    ) %>%
    slice(n()) %>%
    ungroup() %>% 
    select(-c("adc_plasma_uM", "pl_liver_uM"))
  return(df_hep)
}

##======= free payload comparison between T-DM1 and T-Dxd =======##
tdm1_q3w <- read.csv(here("data/sim/trastuzumab-emtansine-q3w.csv")) 
tdxd_q3w <- read.csv(here("data/sim/trastuzumab-deruxtecan-q3w.csv")) 

tdm1_q3w_he_pl <- free_he_pl(tdm1_q3w)
tdxd_q3w_he_pl <- free_he_pl(tdxd_q3w)

rbind(tdm1_q3w_he_pl %>% mutate(ADC = "T-DM1"), 
      tdxd_q3w_he_pl %>% mutate(ADC = "T-Dxd")) %>% 
  select(`Dose (mg/kg)` = dose_mgkg, 
         `Cmax (nM)` = Cmax_he_nM, 
         `Cavg (nM)` = Cavg_he_nM, 
         ADC
         ) %>% 
  st_new() %>%
  stable( r_file = "script/payload-analysis.R", 
          output_file = file.path("deliv/tab/", paste0("trastuzumab-based-adc", ".tex")), 
          panel = "ADC") %>% 
  st2report(stem="trastuzumab-based-adc",output_dir = "deliv/tab/", caption = "Free payload in liver endothelial cells; simulated for 12 weeks")

plt_trastuzumab <- ggplot(
  data = rbind(
    tdm1_q3w %>% filter(dose_mgkg == 3.6) %>% mutate(adc = "T-DM1"), 
    tdxd_q3w %>% filter(dose_mgkg == 5.4) %>% mutate(adc = "T-Dxd")) %>% 
    mutate(Dose = paste0(dose_mgkg, " mg/kg ", adc)) ,  
  aes(x = time_hr/24, y = pl_liver_uM * 1E3, col = Dose, linetype = Dose)
) + 
  geom_line(alpha = 0.7, size = 1.5) + 
  scale_x_continuous(breaks = seq(0, 84, 7)) + 
  scale_y_continuous(limits = c(0.1, 1E2), trans = "log10") + 
  labs(x = "Time (Day)", y = "Free payload concentration (nM)") +
  theme_bw() + theme(legend.position="bottom", panel.grid.minor = element_blank())

plt_trastuzumab

ggsave(here("deliv/figure/pl/trastuzumab-free-pl.png"), plt_trastuzumab, width = 140, height = 100, units = c("mm"))

##======= free payload comparison of cantuzumab mertansine =======##
cm_q3w <- read.csv(here("data/sim/cantuzumab-mertansine-q3w.csv")) 
cm_qw <- read.csv(here("data/sim/cantuzumab-mertansine-qw.csv")) 
cm_tiw <- read.csv(here("data/sim/cantuzumab-mertansine-tiw.csv")) 

cm_q3w_he_pl <- free_he_pl(cm_q3w)
cm_qw_he_pl <- free_he_pl(cm_qw)
cm_tiw_he_pl <- free_he_pl(cm_tiw)

rbind(
  cm_q3w_he_pl %>% mutate(Schedule = "Q3W"), 
  cm_qw_he_pl %>% mutate(Schedule = "QW"), 
  cm_tiw_he_pl %>% mutate(Schedule = "TIW")
) %>% 
  select(`Dose (mg/m2)` = dose_mgkg, 
         `Cmax (nM)` = Cmax_he_nM, 
         `Cavg (nM)` = Cavg_he_nM, 
         Schedule
) %>% 
  st_new() %>%
  stable( r_file = "script/payload-analysis.R", 
          output_file = file.path("deliv/tab/", paste0("cantuzumab-mertansine", ".tex")), 
          panel = "Schedule") %>% 
  st2report(stem="cantuzumab-mertansine",output_dir = "deliv/tab/", caption = "Free payload in liver endothelial cells; simulated for 12 weeks")

# visual comparison between payload dynamics of DM1 in different cantuzumab mertansine dosing scheme 

plt_cm_schedule <- ggplot(
  data = rbind(
  cm_q3w %>% filter(dose_mgkg == 235) %>% mutate(schedule = "Q3W"), 
  cm_qw %>% filter(dose_mgkg == 115) %>% mutate(schedule = "QW"), 
  cm_tiw %>% filter(dose_mgkg == 45) %>% mutate(schedule = "TIW")) %>% 
      mutate(Dose = paste0(dose_mgkg, " mg/m2 ", schedule)) ,  
  aes(x = time_hr/24, y = pl_liver_uM * 1E3, col = Dose, linetype = Dose)
  ) + 
  geom_line(alpha = 0.7, size = 1.5) + 
  scale_x_continuous(breaks = seq(0, 84, 7)) + 
  labs(x = "Time (Day)", y = "Free DM1 concentration (nM)") +
  theme_bw() + theme(legend.position="bottom", panel.grid.minor = element_blank())

plt_cm_schedule

ggsave(here("deliv/figure/pl/cantuzumab-mertansine-free-pl.png"), plt_cm_schedule, width = 140, height = 100, units = c("mm"))

##======= free payload comparison, MMAE =======##
bv_q3w <- read.csv(here("data/sim/brentuximab-vedotin-q3w.csv")) 
pola_q3w <- read.csv(here("data/sim/polatuzumab-vedotin-q3w.csv")) 

bv_q3w_he_pl <- free_he_pl(bv_q3w)
pola_q3w_he_pl <- free_he_pl(pola_q3w)

rbind(bv_q3w_he_pl %>% mutate(ADC = "Brentuximab vedotin"), 
      pola_q3w_he_pl %>% mutate(ADC = "Polatuzumab vedotin")) %>% 
  select(`Dose (mg/kg)` = dose_mgkg, 
         `Cmax (nM)` = Cmax_he_nM, 
         `Cavg (nM)` = Cavg_he_nM, 
         ADC
  ) %>% 
  st_new() %>%
  stable( r_file = "script/payload-analysis.R", 
          output_file = file.path("deliv/tab/", paste0("mmae-based-adc", ".tex")), 
          panel = "ADC") %>% 
  st2report(stem="mmae-based-adc",output_dir = "deliv/tab/", caption = "Free payload in liver endothelial cells; simulated for 12 weeks")


##======= free payload comparison, PBD =======##
vt_q3w <- read.csv(here("data/sim/vadastuximab-talirine-q3w.csv")) 
vt_q4w <- read.csv(here("data/sim/vadastuximab-talirine-q4w.csv")) 


rbind(
  free_he_pl(vt_q3w) %>% mutate(Schedule = "Q3W"), 
  free_he_pl(vt_q4w) %>% mutate(Schedule = "Q4W")
)  %>% 
  select(`Dose (mg/kg)` = dose_mgkg, 
         `Cmax (nM)` = Cmax_he_nM, 
         `Cavg (nM)` = Cavg_he_nM, 
         Schedule
  ) %>% 
  st_new() %>%
  stable( r_file = "script/payload-analysis.R", 
          output_file = file.path("deliv/tab/", paste0("sgn-33a-pbd", ".tex")), 
          panel = "Schedule") %>% 
  st2report(stem="sgn-33a-pbd",output_dir = "deliv/tab/", caption = "Free payload in liver endothelial cells; simulated for 12 weeks")

