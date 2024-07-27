
rm(list = ls())
graphics.off()
# ---------------------------------------------
currentDir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currentDir)

libs <- c("readxl", "readxl", "dplyr", "vars",
          'gridExtra','ggpubr',
          "sandwich", "tseries", "forecast","stargazer",
          "stringr","sandwich","tidyverse",
          "xtable","tidyr","lubridate",'hpfilter')
lapply(libs, require, character.only = T)
# -----------------------------------------------

# OWN FUNS---------------------------
source('stock_pred_fun.R')



# load data----------------------------------
inflation_yoy <- read_excel("PredictorData2022.xlsx", sheet = "inflace_yoy")


inflation_yoy <- inflation_yoy %>%
  pivot_longer(-rok,values_to='cpi_yoy') %>% 
  mutate(datum_help=paste0(rok,"-",name)) %>% 
  mutate(datum=ceiling_date(ym(datum_help), 'month') - days(1),.before = 1) %>%
  dplyr::select(datum,cpi_yoy ) %>% 
  dplyr::arrange(datum) %>% 
  mutate(cpi_yoy_trend=ewma(cpi_yoy,0.987),
         cpi_yoy_detr=cpi_yoy-cpi_yoy_trend) %>% 
  mutate(cpi_yoy_detr_lag=dplyr::lag(cpi_yoy_detr),
         cpi_yoy_trend_lag=dplyr::lag(cpi_yoy_trend),
         cpi_yoy_lag=dplyr::lag(cpi_yoy)) %>%
  dplyr::select(-c('cpi_yoy_detr'))




monthly <- read_excel("PredictorData2022.xlsx", 
                      sheet = "Monthly", col_types = c("text", 
                                                       "numeric", "numeric", "numeric", 
                                                       "numeric", "numeric", "numeric", 
                                                       "numeric", "numeric", "numeric", 
                                                       "numeric", "numeric", "numeric", 
                                                       "numeric", "numeric", "numeric", 
                                                       "numeric", "numeric",'numeric','numeric','numeric','numeric'))


# data manipulation
monthly <- monthly %>% mutate(year=str_sub(yyyymm,1,4), month=str_sub(yyyymm,5,6)) %>% 
  mutate(datum_help=paste0(year,"-",month)) %>% 
  mutate(datum=ceiling_date(ym(datum_help), 'month') - days(1),.before = 1) %>% 
  dplyr::select(-c(month,year,datum_help,yyyymm)) %>% 
  mutate(return_m=log(Index)-lag(log(Index))) %>% 
  filter(row_number()!=1) %>% 
  mutate(rolling_avg_old = cumsum(return_m) / row_number()) %>% 
  mutate(Index_trend=ewma(log(Index),0.98),Index_detr=log(Index)-Index_trend) %>% 
  # Rfree is end-of-month tbl divided by 12 (so applies to the next month)
  # equity premium and log equity premium
  mutate(equity_premium=CRSP_SPvw-Rfree,
         log_equity_premium=log(CRSP_SPvw+1)-log(1+Rfree)) %>%
  mutate(equity_premium_lead=lead(equity_premium)) %>% 
  mutate(log_equity_premium3m = rollapplyr(log_equity_premium, width = 3, FUN = sum, partial = TRUE)) %>%
  mutate(log_equity_premium_lead=lead(log_equity_premium)) %>% 
  # mutate(infl_lag=lag(infl)) %>% 
  filter(year(datum) >1921) %>% 
  
  mutate(TBILL_trend=ewma(tbl,0.98),
         tbl_detr=tbl-TBILL_trend,
         lty_detr=lty-TBILL_trend) %>% 
  mutate(EP=log(E12)-log(Index),
         DP=log(D12)-log(Index)) %>% 
  mutate(tms=tbl-lty) %>% 
  filter(year(datum) >1926) %>% 
  mutate(rolling_avg = cumsum(CRSP_SPvw) / row_number()) %>% 
  mutate(rolling_avg_1m = cumsum(log_equity_premium) / row_number()) %>% 
  mutate(rolling_avg_3m = cumsum(log_equity_premium3m) / row_number()) %>% 
  # remove last row
  filter(row_number() <= n()-1) %>%
  mutate(VAR1=stats::filter(abs(CRSP_SPvw), rep(1,12), sides = 1) * sqrt(pi / 2) *sqrt(12) / 1200) %>% 
  left_join(inflation_yoy, by='datum') %>% 
  mutate(real_rate=tbl-cpi_yoy) %>% 
  mutate(real_rate2=tbl-lag(cpi_yoy)) %>% 
  mutate(real_rate_detr=real_rate-ewma(real_rate,0.98),
         real_rate_detr2=real_rate2-ewma(real_rate,0.98),
         real_rate_tr=ewma(real_rate,0.999)) %>% 
  mutate(tbl_detr2=tbl-cpi_yoy_trend_lag) %>% 
  mutate(tbl_detr3=tbl-ewma(real_rate,0.999)-cpi_yoy_trend_lag) %>% 
  mutate(picovina=tbl-cpi_yoy_trend_lag-real_rate_detr2) %>% 
  mutate(lty_detr1=lty-cpi_yoy_trend_lag) %>% 
  mutate(lty_detr2=lty-ewma(real_rate,0.999)-cpi_yoy_trend_lag) %>% 
  mutate(FED=exp(EP)-lty) %>% 
  mutate(FED_detr=FED-cpi_yoy_trend_lag-pull(hpfilter::hp1(as.data.frame(FED-cpi_yoy_trend_lag),lambda = 128800)))%>% 
  mutate(dfr=corpr-100*ltr ) %>% 
  mutate(dfy=BAA-lty) %>% 
  mutate(dfy_cpi=BAA-lty-cpi_yoy_trend_lag) %>% 
  dplyr::select(-c(Index,D12,E12,CRSP_SPvwx,return_m))%>% 
  # remove column with NA
  filter(year(datum) >1949) %>% 
  mutate(UN_one_sided=UN-pull(hpfilter::hp1(as.data.frame(UN),lambda = 128800))) %>% 
  mutate(UN=lag(UN),
         UN_one_sided=lag(UN_one_sided)) %>% 
  mutate(IP_HP_two_sided=mFilter::hpfilter(INPRO, freq = 128800)$cycle) %>% 
  mutate(IP_HP_two_sided=lag(IP_HP_two_sided)) %>% 
  mutate(IP_HP_one_sided=INPRO-pull(hpfilter::hp1(as.data.frame(INPRO),lambda = 128800))) %>% 
  # interesting FACT !!!  one-sided HP do not forecast ERP
  # TWO sided HP (two sided filter ) forecast ERP
  mutate(IP_HP_one_sided=lag(IP_HP_one_sided)) %>% 
  mutate(dfr_detr=pull(hp1(as.data.frame(dfr,lambda = 128800)))) %>% 
  mutate(dfy_detr=pull(hp1(as.data.frame(dfy_cpi,lambda = 128800)))) %>% 
  filter(year(datum) >1949) %>% 
  # remove column with NA
  select_if(~ !any(is.na(.))) %>% 
  filter(year(datum) >1950) %>% 
  mutate(log_equity_premium3 = rollapplyr(log_equity_premium_lead, width = 3, FUN = sum, partial = TRUE)) %>%
  mutate(log_equity_premium3 = dplyr::lead(log_equity_premium3,2))



#### 1 or 3 month ####
#### chose one !!!!! ####
#### option A
monthly_selected <- monthly %>% dplyr::select(log_equity_premium_lead,cpi_yoy_detr_lag,cpi_yoy_lag,infl,tbl_detr2,tbl_detr3,tbl,lty_detr1,lty_detr2,lty)
horizon <- 1
monthly$rolling_avg_fin <- monthly$rolling_avg_1m
# -----------------
#### option B
# monthly_selected <- monthly %>% dplyr::select(log_equity_premium3,cpi_yoy_detr_lag,cpi_yoy_lag,infl,tbl_detr2,tbl_detr3,tbl,lty_detr1,lty_detr2,lty)
# horizon <- 3
# monthly_selected <- monthly_selected %>% rename(log_equity_premium_lead=log_equity_premium3)
# monthly$rolling_avg_fin <- monthly$rolling_avg_3m
# ------------


# in sample test
model_list <- list()
pval_list <- list()
std_error_list <- list()

replication_n <- 10000
for (i in 2:dim(monthly_selected)[2]) {
  
  pomocny_dataset <- monthly_selected %>% dplyr::select(1,i) %>% filter(complete.cases(.))
  model <- lm(log_equity_premium_lead~ ., data=pomocny_dataset)
  std_error <- coeftest(model, vcov.=NeweyWest(model, lag=horizon, prewhite=FALSE, adjust=TRUE, verbose=TRUE))
  print(summary(model))
  pval <- pval_boot_fun(x=10*as.numeric(pull(pomocny_dataset[,2])), y=10*pomocny_dataset$log_equity_premium_lead,replication=replication_n)
  model_list[[(i-1)]] <- model
  pval_list[[(i-1)]] <- pval
  std_error_list[[(i-1)]] <- std_error
}




# model summary
stargazer(model_list[[1]],model_list[[2]],model_list[[3]],model_list[[4]],model_list[[5]],model_list[[6]],model_list[[7]],model_list[[8]],model_list[[9]]
          
)


stargazer(model_list[[1]],model_list[[2]],model_list[[3]],model_list[[4]],model_list[[5]],
          model_list[[6]],model_list[[7]],model_list[[8]],model_list[[9]],
          se = list(std_error_list[[1]][,"Std. Error"], std_error_list[[2]][,"Std. Error"]
                    , std_error_list[[3]][,"Std. Error"], std_error_list[[4]][,"Std. Error"]
                    , std_error_list[[5]][,"Std. Error"], std_error_list[[6]][,"Std. Error"]
                    , std_error_list[[7]][,"Std. Error"]
                    , std_error_list[[8]][,"Std. Error"]
                    , std_error_list[[9]][,"Std. Error"])
          
)


pval_list


# --------------------------------



# Rolling R squared ------------------------------------- (not used)
R_result <- c()
coef_result <- c()
for (ii in 24:dim(monthly)[1]) {
  holder_temp <- summary(lm(monthly$equity_premium_lead[1:ii]~ monthly$tbl_detr2[1:ii]))
  holder_temp$adj.r.squared
  R_result <- append(R_result,holder_temp$adj.r.squared)
  
  coef_result <- append(coef_result,coef(holder_temp)[2])
  
  
}
ts.plot(R_result)
abline(h=0)

# --------------------------




# out of sample test
R2_OS <- list()
CW <- list()

##### CHANGE HERE !!
# chose one
# cambell_restriction  <-  TRUE
cambell_restriction  <-  FALSE


for (i in 2:dim(monthly_selected)[2]) {
  
  pomocny_dataset=monthly_selected %>% dplyr::select(1,i)
  vv1 <- c()
  vv2 <- c()
  f1 <- c()
  f2 <- c()
  
  coeff=c()
  
  for (ii in 72:(dim(pomocny_dataset)[1]-horizon)) {
    holder_temp <- coef(lm(pomocny_dataset$log_equity_premium_lead[1:ii]~ pull(pomocny_dataset[1:ii,2])))
    
    if(cambell_restriction==TRUE){
      # cambell thomson restriction 1
      holder_temp[2] <- ifelse(holder_temp[2]>0,0,holder_temp[2])
    }
    
    coeff <- append(coeff,holder_temp[2])
    
    
    
    result <- holder_temp[2]*pull(pomocny_dataset[(ii+1),2])+holder_temp[1]
    
    if(cambell_restriction==TRUE){
      # cambell thomson restriction 2
      result <- ifelse(result<0,0,result)
    }
    
    
    v1 <- as.numeric(pomocny_dataset$log_equity_premium_lead[(ii+1)]-result)
    
    v2 <- as.numeric(pomocny_dataset$log_equity_premium_lead[(ii+1)]-monthly$rolling_avg_fin[ii])
    
    vv1 <- append(vv1,v1)
    vv2 <- append(vv2,v2)
    
    ff1 <- as.numeric(result)
    ff2 <- as.numeric(monthly$rolling_avg_fin[ii])
    
    f1 <- append(f1,ff1)
    f2 <- append(f2,ff2)
    
    
    
    
  }
  ts.plot(coeff)
  abline(h=0)
  
  R2_OS[[i]] <- round((1-THEIL_U_RMSE_R(vv1,vv2)),3)
  print(colnames(pomocny_dataset))
  
  print(  round((1-THEIL_U_RMSE_R(vv1,vv2))*100,3))
  
  CW[[i]] <- round(clark_west(monthly$equity_premium_lead[73:length(monthly$equity_premium_lead)], f2, f1,1)$pvalue,3)
  
  print(  clark_west(monthly$equity_premium_lead[73:length(monthly$equity_premium_lead)], f2, f1,1))
  print(  cw(vv2,vv1,f2,f1,1) )
  
  
  
}


unlist(R2_OS)
unlist(CW)






