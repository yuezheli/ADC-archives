# date: 4/5/24
# author: Yuezhe Li 
# purpose: functions to calculate partition coefficients 
# https://dmd.aspetjournals.org/content/48/10/903

library(dplyr)

# all default params were updated to doxorubicin 
# https://pubmed.ncbi.nlm.nih.gov/35335919/

calcKp_PT <- function(logP= 1.27, pKa = c(9.53, 8.94), fup = 0.25, BP=1.15, type=6){

  dat <- read.csv( paste(getwd(), "dox-pbpk/tissue_comp_P_T.csv", sep = "/") ) 
  
  dat_all <- dat %>% filter(!tissue %in% c("Plasma","Adipose","RBCs"))
  
  n <- length(dat$tissue)
  Kp_all <- vector(mode = "numeric", length = n)
  
  Vwp <- dat$f_water[dat$tissue == "Plasma"]
  Vnlp <- dat$f_n_l[dat$tissue == "Plasma"]
  Vphp <- dat$f_pl[dat$tissue == "Plasma"]
  
  dat2 <- dat %>% filter(!tissue %in% c("Plasma","RBCs"))
  
  Vwt <- dat2$f_water[dat2$tissue != "Adipose"]
  Vwad <- dat2$f_water[dat2$tissue == "Adipose"]
  Vnlt <- dat2$f_n_l[dat2$tissue != "Adipose"]
  Vnlad <- dat2$f_n_l[dat2$tissue == "Adipose"]
  Vpht <- dat2$f_pl[dat2$tissue != "Adipose"]
  Vphad <- dat2$f_pl[dat2$tissue == "Adipose"]
  
  pH <- dat$pH[dat$tissue == "Adipose"]
  logD <- 1.115*logP-1.35 #logD is the olive oil:buffer(water) partition coefficient of nonionized species
  
  logD_star <- switch(type,
                      #1-neutral
                      logD,   
                      #2-monoprotic acid
                      logD-log10(1+10^(pH-pKa)),
                      #3-monoprotic base
                      logD-log10(1+10^(pKa-pH)), 
                      #4-diprotic acid
                      logD-log10(1+10^(2*pH-pKa[1]-pKa[2])),
                      #5-diprotic base
                      logD-log10(1+10^(pKa[1]+pKa[2]-2*pH)), 
                      #6-monoprotic acid monoprotic base (acid comes first)
                      logD-log10(1+10^(pKa[2]-pKa[1])),  
                      #7-triprotic acid
                      logD-log10(1+10^(3*pH-pKa[1]-pKa[2]-pKa[3])),  
                      #8-triprotic base
                      logD-log10(1+10^(pKa[1]+pKa[2]+pKa[3]-3*pH)),  
                      #9-diprotic acid monoprotic base (first two are acid)
                      logD-log10(1+10^(pH-pKa[1]-pKa[2]+pKa[3])), 
                      #10-diprotic base monoprotic acid (first one is acid)
                      logD-log10(1+10^(pKa[2]+pKa[3]-pKa[1]-pH)))       
  
  D_star <- 10^logD_star   
  Kpad <- ((D_star*(Vnlad+0.3*Vphad)+(1*(Vwad+0.7*Vphad)))/(D_star*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))) * fup
  
  P <- 10^logP
  fut <- 1/(1+((1-fup)/fup)*0.5)
  Kpt <- ((P*(Vnlt+0.3*Vpht)+(1*(Vwt+0.7*Vpht)))/(P*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))) * (fup/fut)
  
  nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  nms_all <- paste("Kp", nms_all, sep="")
  nms <- c("Kpad",nms_all)
  Kp <- as.list(c(Kpad,Kpt))
  names(Kp) <- nms
  
  return(Kp)
  
}

calcKp_Schmitt <- function(logP= 1.27, pKa = c(9.53, 8.94), fup = 0.25, BP=1.15, type=6){
  dat <- read.csv( paste(getwd(), "dox-pbpk/tissue_comp_Schmitt.csv", sep = "/") )
  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Plasma")) 
  logMA <- logP  #in case we don't have a direct logMA
  K_n_pl <- 10^logMA    #neutral phospholipids:water partition coefficient
  #K_protein <- 0.163+0.0221*K_n_pl    #protein:water partition; Schmitt, Walter (2008)
  K_protein <- ((0.81 + 0.11 * K_n_pl)/24.92)*5
  pH <- dat_all$pH
  alpha <- 1e-3  #ratio between distribution coefficient at given pH (D) and that in neutral form (D0)
  W <- switch(type,
              #1-neutral
              0,
              #2-monoprotic acid
              10^(pH-pKa),
              #3-monoprotic base
              10^(pKa-pH),
              #4-diprotic acid
              10^(pH-pKa[1])+10^(2*pH-pKa[1]-pKa[2]), 
              #5-diprotic base
              10^(pKa[2]-pH)+10^(pKa[1]+pKa[2]-2*pH), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH)+10^(pH-pKa[1]),
              #7-triprotic acid
              10^(pH-pKa[1])+10^(2*pH-pKa[1]-pKa[2])+10^(3*pH-pKa[1]-pKa[2]-pKa[3]))
  
  if(type==1 | type==2 | type==4 | type==7)
  { # neutral, monoprotic acid, diprotic acid, triprotic acid
    K_n_l <- K_n_pl*(((1-alpha)/(1+W))+alpha)
    K_a_pl <- K_n_pl*((1/(1+W))+0.05*(1-(1/(1+W)))) 
  }
  else
  {
    if(type==3)
      { # monoprotic base
        K_n_l <- K_n_pl*(((1-alpha)/(1+W))+alpha) 
        K_a_pl <- K_n_pl*((1/(1+W))+20*(1-(1/(1+W)))) 
    }
    if(type==5)
      { # diprotic base 
        F1 <- (1/(1+10^(pKa[1]-pH)))
        F2 <- (1/(1+10^(pKa[2]-pH)))
        K_n_l <- K_n_pl*(F1*F2 + alpha*((1-F1)*F2 + F1*(1-F2)) + (1-F1)*(1-F2))
        K_a_pl <- K_n_pl*(F1*F2 + 20*((1-F1)*F2 + F1*(1-F2)) + (1-F1)*(1-F2))
      }
    if(type==6)
      { # monoprotic acid monoprotic base (acid comes first)
        F1 <- (1/(1+10^(pH-pKa[1])))
        F2 <- (1/(1+10^(pKa[2]-pH)))
        K_n_l <- K_n_pl*(F1*F2 + alpha*((1-F1)*F2 + F1*(1-F2)) + (1-F1)*(1-F2))
        K_a_pl <- K_n_pl*(F1*F2 + 0.05*(1-F1)*F2 + 20*(1-(F2))*F1 + (1-F1)*(1-F2))
      }
  }
  kp <- (dat_all$f_water+(K_n_l*dat_all$f_n_l)+(K_n_pl*dat_all$f_n_pl)+(K_a_pl*dat_all$f_a_pl)+(K_protein*dat_all$f_proteins))*fup
  
  dat2 <- data.frame(tissue=dat_all$tissue, Kp=kp)
  name <- dat2$tissue %>% substr(1,2) %>% tolower()
  name <- paste("Kp", name, sep="")
  Kp <- as.list(dat2$Kp)
  names(Kp) <- name
  return(Kp)
}









