setwd("C:/Users/mauro/Downloads/final553")

library(readxl)
library(stats)
library(tidyverse)
library(dampack)
library(plotly)

# 1. Optimization ########################################################################
us_states <- read_excel("us-states.xlsx")
ct <- subset(us_states,state=='Connecticut')

ct$cases_1 <- lag(ct$cases,1,NA)
ct$new_cases <- ct$cases-ct$cases_1
ct$new_cases[1] <- 0

rho <- 1/14
res <- c()

f <- function(beta){
  for (i in 1:100){
    res <- append(res,exp((beta - rho) * 3 * i) - ct$cases[i])
  }
  sum(res**2)
}

beta_opt <- optimize(f,interval = c(0, 1),tol = 0.0001)

# 2.Plot ########################################################################
simu0 <- function(test_freq,test_se,R){
  
  initial_susceptible <- 4990
  initial_infected <- 10
  R0 <- R
  cycles.per.day <- 3
  
  days_to_incubation <- 3
  time_to_recovery <- 14
  
  pct_advancing_to_symptoms <- 30
  symptom_case_fatality_ratio <- 0.0005
  
  
  test_sensitivity <- test_se
  test_specificity <- ifelse(R==1.5,0.997,ifelse(R==2.5|R==3.5,0.98,NA))
  time_to_return_fps <- 1
  
  new_infections_per_shock <- ifelse(R==1.5,5,ifelse(R==2.5,10,ifelse(R==3.5, 25,NA)))
  num.exogenous.shocks <- 1
  frequency_exogenous_shocks <- 7
  frequency.exogenous.shocks <- cycles.per.day*frequency_exogenous_shocks
  
  cycles.per.test <- test_freq*cycles.per.day
  
  rho <- 1/(time_to_recovery*cycles.per.day)
  sigma <- rho*(pct_advancing_to_symptoms/100/(1-pct_advancing_to_symptoms/100))
  beta <- R0*(rho+sigma)
  delta <- (symptom_case_fatality_ratio/(1-symptom_case_fatality_ratio))*rho
  theta <- 1/(days_to_incubation*cycles.per.day)
  mu <- 1/(cycles.per.day*time_to_return_fps)
  
  n.cycle <- 240
  
  mat <- matrix(c(0,initial_susceptible,0,0,initial_infected,0,0,0,0), nrow = 1)
  mat <- rbind(mat,c(1,
                     max(0,mat[1,2]*(1-beta*(mat[1,5]/(mat[1,2]+mat[1,5]+mat[1,4])))+mat[1,3]*mu),
                     max(0,mat[1,3]*(1-mu)),
                     max(0,mat[1,4]*(1-theta)+ beta*(mat[1,2]*mat[1,5]/(mat[1,2]+mat[1,5]+mat[1,4]))),
                     max(0,mat[1,5]*(1-sigma-rho)+mat[1,4]*theta),
                     max(0,mat[1,6]*(1-delta-rho)+(mat[1,5]+mat[1,7])*sigma),
                     0,
                     max(0,mat[1,8]+(mat[1,5]+mat[1,6]+mat[1,7])*rho),
                     max(0,delta*mat[1,6]+mat[1,9])))
  
  superspreader.event <- 0
  superspreader.event <- c(superspreader.event, 
                           (1:n.cycle %% frequency.exogenous.shocks == 0)*num.exogenous.shocks)
  
  for(i in 2:n.cycle) {
    mat <- 
      rbind(
        mat,
        c(i,
          max(0,mat[i,2]*(1-beta*(mat[i,5]/(mat[i,2]+mat[i,5]+mat[i,4])))+mat[i,3]*mu-mat[i-1,2]*(1-test_specificity)/cycles.per.test-superspreader.event[i+1]*new_infections_per_shock),
          max(0,mat[i,3]*(1-mu)+mat[i-1,2]*(1-test_specificity)/cycles.per.test),
          max(0,mat[i,4]*(1-theta)+beta*(mat[i,2]*mat[i,5]/(mat[i,2]+mat[i,5]+mat[i,4]))+superspreader.event[i+1]*new_infections_per_shock),
          max(0,mat[i,5]*(1-sigma-rho)+mat[i,4]*theta-mat[i-1,5]*test_sensitivity/cycles.per.test),
          max(0,mat[i,6]*(1-delta-rho)+(mat[i,5]+mat[i,7])*sigma),
          max(0,mat[i,7]*(1-sigma-rho)+mat[i-1,5]*test_sensitivity/cycles.per.test),
          max(0,mat[i,8]+(mat[i,5]+mat[i,6]+mat[i,7])*rho),
          max(0,delta*mat[i,6]+mat[i,9]))
      )
  }
  
  mat <- cbind(mat, superspreader.event)
  mat[241,c(3,6,7)]
  
  names.df <- c("Cycle","Susceptible","FP","Exposed","Asympt","Symptoms","TP","Recovered","Dead","Superspreader Event")
  
  
  df <- 
    mat %>% 
    as_tibble() %>% 
    rename_all(~names.df) %>% 
    mutate(`Persons Tested` = (lag(Susceptible,1,NA)+lag(Exposed,1,NA)+lag(Asympt,1,NA))/cycles.per.test,
           `Total TPs` = lag(Asympt,2,NA)*test_sensitivity/cycles.per.test,
           `Total FPs` = lag(Susceptible,2,NA)*(1-test_specificity)/cycles.per.test,
           `Total TNs` = lag(Susceptible,2,NA)*test_specificity/cycles.per.test,
           `Total FNs` = lag(Exposed,2,NA)+lag(Asympt,2,NA)*(1-test_sensitivity)/cycles.per.test) %>% 
    mutate(Day = Cycle/cycles.per.day,
           `True Positive` = TP,
           Symptoms = Symptoms,
           `False Positive` = FP,
           Total = TP+Symptoms+FP) %>% 
    mutate(`New Infections` = lag(Asympt,1,NA)*beta*lag(Susceptible,1,NA)/(lag(Susceptible,1,NA)+lag(Exposed,1,NA)+lag(Asympt,1,NA)),
           `New Infections` = ifelse(Cycle>1,
                                     `New Infections`+pmin(`Superspreader Event`*new_infections_per_shock,lag(Susceptible,1,NA)),
                                     `New Infections`),
           `New Infections` = ifelse(is.na(`New Infections`),0,`New Infections`),
           `Cumulative Infections` = cumsum(`New Infections`),
           `%Cumulative Infections` = `Cumulative Infections`/initial_susceptible)
  return(df)
}

simu <- function(test_freq,test_se,beta_calc,R=2.5){
  
  options(scipen = 200)
  
  initial_susceptible <- 4990
  initial_infected <- 10
  R0 <- R
  cycles.per.day <- 3
  
  days_to_incubation <- 3
  time_to_recovery <- 14
  
  pct_advancing_to_symptoms <- 30
  symptom_case_fatality_ratio <- 0.0005
  
  test_sensitivity <- test_se
  test_specificity <- 0.98
  time_to_return_fps <- 1
  
  new_infections_per_shock <- 10
  num.exogenous.shocks <- 1
  frequency_exogenous_shocks <- 7
  frequency.exogenous.shocks <- cycles.per.day*frequency_exogenous_shocks
  
  cycles.per.test <- test_freq*cycles.per.day
  
  rho <- 1/(time_to_recovery*cycles.per.day)
  sigma <- rho*(pct_advancing_to_symptoms/100/(1-pct_advancing_to_symptoms/100))
  
  beta <- ifelse(beta_calc==0,R0*(rho+sigma),ifelse(beta_calc==1,beta_opt$minimum,NA))
  delta <- (symptom_case_fatality_ratio/(1-symptom_case_fatality_ratio))*rho
  theta <- 1/(days_to_incubation*cycles.per.day)
  mu <- 1/(cycles.per.day*time_to_return_fps)
  
  n.cycle <- 240
  
  mat <- matrix(c(0,initial_susceptible,0,0,initial_infected,0,0,0,0,0), nrow = 1)
  mat <- rbind(mat,c(1,
                     max(0,mat[1,2]*(1-beta*(mat[1,5]/(mat[1,2]+mat[1,5]+mat[1,4])))+mat[1,3]*mu),
                     max(0,mat[1,3]*(1-mu)),
                     max(0,mat[1,4]*(1-theta)+ beta*(mat[1,2]*mat[1,5]/(mat[1,2]+mat[1,5]+mat[1,4]))),
                     max(0,mat[1,5]*(1-sigma-rho)+mat[1,4]*theta),
                     max(0,mat[1,6]*(1-delta-rho)+(mat[1,5]+mat[1,7])*sigma),
                     0,
                     max(0,mat[1,8]+(mat[1,5]+mat[1,6]+mat[1,7])*rho),
                     max(0,delta*mat[1,6]+mat[1,9]),
                     0))
  
  superspreader.event <- 0
  superspreader.event <- c(superspreader.event, 
                           (1:n.cycle %% frequency.exogenous.shocks == 0)*num.exogenous.shocks)
  
  for(i in 2:n.cycle) {
    mat <- 
      rbind(
        mat,
        c(i,
          max(0,mat[i,2]*(1-beta*((mat[i,5]+mat[i,10])/(mat[i,2]+mat[i,5]+mat[i,4]+mat[i,10])))+mat[i,3]*mu-mat[i-1,2]*(1-test_specificity)/cycles.per.test-superspreader.event[i+1]*new_infections_per_shock),
          max(0,mat[i,3]*(1-mu)+mat[i-1,2]*(1-test_specificity)/cycles.per.test),
          max(0,mat[i,4]*(1-theta)+beta*(mat[i,2]*(mat[i,5]+mat[i,10])/(mat[i,2]+mat[i,5]+mat[i,4]+mat[i,10]))+superspreader.event[i+1]*new_infections_per_shock),
          max(0,mat[i,5]*(1-sigma-rho)+mat[i,4]*theta-mat[i-1,5]*test_sensitivity/cycles.per.test),
          max(0,mat[i,6]*(1-delta-rho)+(mat[i,5]+mat[i,7])*sigma-mat[i-1,6]*test_sensitivity/cycles.per.test+mat[i,10]*mu),
          max(0,mat[i,7]*(1-sigma-rho)+mat[i-1,5]*test_sensitivity/cycles.per.test),
          max(0,mat[i,8]+(mat[i,5]+mat[i,6]+mat[i,7])*rho),
          max(0,mat[i,9]+delta*mat[i,6]),
          max(0,mat[i,10]*(1-mu)+mat[i-1,6]*test_sensitivity/cycles.per.test)
      ))
  }
  
  mat <- cbind(mat, superspreader.event)
  mat[241,c(3,6,7,10)]
  
  names.df <- c("Cycle","Susceptible","FP","Exposed","Asympt","Symptoms","TP","Recovered","Dead","False Negative","Superspreader Event")
  
  df_stats <- 
    mat %>% 
    as_tibble() %>% 
    rename_all(~names.df) %>% 
    mutate(`Persons Tested` = (lag(Susceptible,1,NA)+lag(Exposed,1,NA)+lag(Asympt,1,NA))/cycles.per.test,
           `Total TPs` = lag(Asympt,2,NA)*test_sensitivity/cycles.per.test,
           `Total FPs` = lag(Susceptible,2,NA)*(1-test_specificity)/cycles.per.test,
           `Total TNs` = lag(Susceptible,2,NA)*test_specificity/cycles.per.test,
           `Total FNs` = lag(Exposed,2,NA)+lag(Asympt,2,NA)*(1-test_sensitivity)/cycles.per.test) %>% 
    mutate(Day = Cycle/cycles.per.day,
           `True Positive` = TP,
           Symptoms = Symptoms,
           `False Positive` = FP,
           Total = TP+Symptoms+FP) %>% 
    mutate(`New Infections` = (lag(Asympt,1,NA)+lag(`False Negative`,1,NA))*beta*lag(Susceptible,1,NA)
           /(lag(Susceptible,1,NA)+lag(Exposed,1,NA)+lag(Asympt,1,NA)+lag(`False Negative`,1,NA)),
           `New Infections` = ifelse(Cycle>1,
                                     `New Infections`+pmin(`Superspreader Event`*new_infections_per_shock,lag(Susceptible,1,NA)),
                                     `New Infections`),
           `New Infections` = ifelse(is.na(`New Infections`),0,`New Infections`),
           `Cumulative Infections` = cumsum(`New Infections`),
           `%Cumulative Infections` = `Cumulative Infections`/initial_susceptible)
  return(df_stats)
}

iso_plot <- function(df_stats,main){
  df_stats %>%
    select(Day, `True Positive`, Symptoms, `False Positive`) %>%
    pivot_longer(`True Positive`:`False Positive`, names_to = "Group", values_to = "Value") %>%
    mutate(Group = as.factor(Group),
           Group = forcats::fct_relevel(Group, levels = c("True Positive", "Symptoms", "False Positive")),
           Group = forcats::fct_recode(Group,
                                       "Asymptomatic (TP)" = "True Positive",
                                       "Symptomatic" = "Symptoms",
                                       "Uninfected (FP)" = "False Positive")) %>%
    group_by(Day) %>%
    arrange(Group) %>%
    mutate(`New Students` = sum(Value),
           Students = cumsum(Value)) %>%
    plot_ly(x = ~Day,
            y = ~Students,
            color = ~Group,
            colors = RColorBrewer::brewer.pal(9,"YlOrRd")[c(3,6,9)],
            alpha = 0.7,
            type = "scatter",
            mode = "lines",
            fill = 'tonexty',
            hoverinfo = "text") %>%
    layout(title = "Composition of Isolation Pool") %>%
    layout(yaxis = list(title = "Number of Students")) %>%
    layout(autosize = TRUE,
           margin = list(l = 75,
                         r = 75,
                         b = 75,
                         t = 75,
                         pad = 10)) %>%
    config(displaylogo = FALSE)
}




# 3. CEA Analysis ########################################################################

cea_stats <- function(cost,df_stats){
 
  test_cost <- cost
  confirmatory_test_cost <-100
  
  sum.stat <-
    df_stats %>%
    slice(2:n()) %>%
    summarize(`Total Persons Tested in 80 days` = sum(`Persons Tested`, na.rm = TRUE),
              `Total Confirmatory Tests Performed` = sum(`Total TPs`, na.rm = TRUE) + sum(`Total FPs`, na.rm = TRUE),
              `Average Isolation Unit Census` = mean(`Total`, na.rm = TRUE),
              `Average %TP in Isolation` = 1-(mean(`False Positive`, na.rm = TRUE)/mean(`Total`, na.rm = TRUE)),
              `Total testing cost` = `Total Persons Tested in 80 days`*test_cost+`Total Confirmatory Tests Performed`*confirmatory_test_cost,
              `Total Infections` = last(`Cumulative Infections`))

  sum.stat <- list(
    ## Expected outputs
    number_tested = sum.stat$`Total Persons Tested in 80 days`,
    number_confirmatory_tests = sum.stat$`Total Confirmatory Tests Performed`,
    average_iu_census = sum.stat$`Average Isolation Unit Census`,
    average_pct_isolated = sum.stat$`Average %TP in Isolation`,
    testing_cost = sum.stat$`Total testing cost`,
    infections = sum.stat$`Total Infections`
  )

  sum.stat
  return(sum.stat)
}


cea_output <- function(beta_calc, R=2.5){
  cea_table <- matrix(nrow=13,ncol=5,NA)
  
  test_freq_name <- c('Daily','2 day','3 day','Weekly')
  test_se_name <- c('70', '80','90')
  test_freq_ <- c(1, 2, 3, 7)
  test_se_ <- c(0.7, 0.8, 0.9)
  test_cost <- c(10,25,50)
  
  for (i in 1:4){
    for (k in 1:3){
      cea_table[(i-1)*3+k+1,1] <- test_freq_name[i]
      cea_table[(i-1)*3+k+1,2] <- test_se_name[k]
      cea_table[(i-1)*3+k+1,3] <- cea_stats(test_cost[k],simu(test_freq_[i], test_se_[k], beta_calc, R))$testing_cost
      cea_table[(i-1)*3+k+1,4] <- cea_stats(test_cost[k],simu(test_freq_[i], test_se_[k], beta_calc, R))$infections
    }
  }
  
  cea_table[1,1] <- 'Symptom Based'
  cea_table[1,3] <- 0
  cea_table[1,4] <- cea_stats(0, simu(9999999999999999999999, 0.7, beta_calc, R))$infections
  
  cea_table <- cea_table[order(-as.numeric(cea_table[,4])),]
  
  cea_table[1,3]<-0
  
  icer <- calculate_icers(cost=as.numeric(cea_table[1:13,3]), 
                          effect=-as.numeric(cea_table[1:13,4]), 
                          strategies=c(1:13))
  
  as.numeric(icer$Strategy)
  ND <- as.numeric(subset(icer, Status=='ND')$Strategy)
  
  cea_table[1,5] <- NA
  cea_table[2:13,5] <- 'Dominated'
  
  cea_table[ND[1],5] <- (as.numeric(cea_table[ND[1],3])-as.numeric(cea_table[1,3]))/(as.numeric(cea_table[1,4])-as.numeric(cea_table[ND[1],4]))
  for (i in 2:length(ND)){
    cea_table[ND[i],5] <- (as.numeric(cea_table[ND[i],3])-as.numeric(cea_table[ND[i-1],3]))/(as.numeric(cea_table[ND[i-1],4])-as.numeric(cea_table[ND[i],4]))
  }
  cea_table[1,3]<-NA
  cea_table <- as.data.frame(cea_table)
  colnames(cea_table) <- c('Frequency','Sensitivity','Cost','Infection','CEA')
  return(cea_table)
  
}

iso_plot(simu(test_freq=1, test_se=0.7, beta_calc=0, R=2.5))
iso_plot(simu(test_freq=2, test_se=0.7, beta_calc=0, R=2.5))
iso_plot(simu(test_freq=7, test_se=0.7, beta_calc=0, R=2.5))
iso_plot(simu(test_freq=99999999999999, test_se=0.7, beta_calc=0, R=2.5))

iso_plot(simu(test_freq=1, test_se=0.7, beta_calc=1))
iso_plot(simu(test_freq=2, test_se=0.7, beta_calc=1))
iso_plot(simu(test_freq=7, test_se=0.7, beta_calc=1))
iso_plot(simu(test_freq=99999999999999, test_se=0.7, beta_calc=1))

cea_output(beta_calc=0, R=2.5)
cea_output(beta_calc=0, R=3.5)
cea_output(beta_calc=0, R=1.5)

cea_output(beta_calc=1)
#write.csv(cea_table,'cea_scenario0.csv')
write.csv(simu0(test_freq=2, test_se=0.7, R=2.5),'table0.csv')

write.csv(simu(test_freq=2, test_se=0.7, beta_calc=0, R=2.5),'table1.csv')

par(mfrow=c(2,2))
iso_plot(simu(test_freq=1, test_se=0.7, beta_calc=0, R=2.5))  #daily testing 
iso_plot(simu(test_freq=2, test_se=0.7, beta_calc=0, R=2.5))  #testing per 2 days
iso_plot(simu(test_freq=7, test_se=0.7, beta_calc=0, R=2.5))  #weekly testing 
iso_plot(simu(test_freq=99999999999999, test_se=0.7, beta_calc=0, R=2.5))  #symptom-based testing



write.csv(simu(test_freq=2, test_se=0.7, beta_calc=1),'table2.csv')


iso_plot(simu(test_freq=1, test_se=0.7, beta_calc=1))
iso_plot(simu(test_freq=2, test_se=0.7, beta_calc=1))
iso_plot(simu(test_freq=7, test_se=0.7, beta_calc=1))
iso_plot(simu(test_freq=99999999999999, test_se=0.7, beta_calc=1))

cea_output(beta_calc=0, R=2.5)
write.csv(cea_output(beta_calc=0, R=2.5),'cea_base.csv')
write.csv(cea_output(beta_calc=0, R=3.5),'cea_worst.csv')
write.csv(cea_output(beta_calc=0, R=1.5),'cea_best.csv')
write.csv(cea_output(beta_calc=1),'cea_Con.csv')

cea_output(beta_calc=1)  #Connecticut scenario
