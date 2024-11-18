# date: 4/5/24
# author: Yuezhe Li 
# purpose of this code: to test a PBPK model for Dox 

rm(list=ls())

library(dplyr)
library(ggplot2)
library(gridExtra)
library(mrgsolve)

# read in observed data 
# obs <- read.csv("data/dox_pk_obs_2.csv")
obs <- read.csv("data/dox_pk_obs.csv") %>% filter(Species == "DOX")

# molecular weight 
MW = 543 # [g/mol]

source("dox-pbpk/PartitionCoef.R")

Kp_PT <- calcKp_PT(logP= 1.27, pKa = c(9.53, 8.94), fup = 0.25, BP=1.15, type=6)

modA <- mread("dox-pbpk/small_mol_PBPK") %>% param(lapply(Kp_PT,"*",5.3119)) %>% init(URINE = 0)

# single IV bolus dose, 60mg/m2
sim1 <- modA %>% ev(amt=60*1.73*1E3/MW, cmt="VEN", addl=0, rate=0) %>% mrgsim(delta = 0.1, end = 168) %>% as.data.frame() %>% filter(time > 0)
sim2 <- modA %>% ev(amt=50*1.73*1E3/MW, cmt="VEN", addl=0, rate=0) %>% mrgsim(delta = 0.1, end = 168) %>% as.data.frame() %>% filter(time > 0)
sim3 <- modA %>% ev(amt=75*1.73*1E3/MW, cmt="VEN", addl=0, rate=0) %>% mrgsim(delta = 0.1, end = 168) %>% as.data.frame() %>% filter(time > 0)

sim_combined <- rbind(sim1 %>% mutate(Dose = "60mg/m2"), 
                      sim2 %>% mutate(Dose = "50mg/m2"),
                      sim3 %>% mutate(Dose = "75mg/m2")) %>% 
  select(time, Cvenous, Dose)

#p_plasma <- ggplot(data = sim1, aes(x = time, y = Cvenous*1E3*MW)) + 
#  geom_line() + 
#  geom_point(data = obs, aes(x = Time_h, y = Conc_mgL)) + 
#  coord_cartesian(ylim = c(0.001, 1)) + scale_y_log10() + theme_bw()

p_plasma <- ggplot(data = sim_combined, aes(x = time, y = Cvenous*0.25/1.15, group = Dose, col = Dose)) + 
  geom_line() + 
  geom_point(data = obs, aes(x = Time_min/60, y = Conc_uM, col = Dose)) + scale_y_log10()

p_plasma

