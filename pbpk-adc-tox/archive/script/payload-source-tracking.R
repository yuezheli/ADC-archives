# date: 7/2/2025 
# author: Yuezhe Li 
# purpose of this code: visualization of bar plot of payload source 

rm(list=ls()) 
gc()

library(here)
library(data.table)
library(dplyr)
library(tidyverse)
library(pmtables)
library(ggplot2)

# read in simulation data 
pl <- read.csv(here("data/sim/free-pl-tracking.csv")) %>% 
  arrange(desc(deg_mem_perc)) %>%
  rename(`nonspecific` = deg_mem_perc, 
         `FcRn` = deg_fcrn_perc)
  

pl_stack <- melt(pl, id.vars = "Organ") 

# visualization 
p_deg_pl <- ggplot(pl_stack, aes(fill = variable, x = Organ, y = value) ) + 
  geom_bar(position = position_dodge(), stat = "identity", width = 0.7, alpha = 0.6) + 
  labs(x = "", y = "% of payload released from ADC degredation \n in endothelial cells", fill = "Source") +
  facet_wrap(~variable, scales = "free_y") + 
  scale_x_discrete(limits = pl$Organ) +
  theme_bw() + theme(panel.grid.minor = element_blank(), legend.position="bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

p_deg_pl

# save figure 
ggsave(here("deliv/figure/pl/deg-pl-distribution.png"), p_deg_pl, width = 120, height = 100, units = c("mm"))
