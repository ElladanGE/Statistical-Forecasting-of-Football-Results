# References :

# Data taken from : https://fbref.com/en/comps/9/schedule/Premier-League-Scores-and-Fixtures

# Modelling association football scores by M. J. MAHER

# Analysis of Sports Data by Using Bivariate Poisson Models
# Author(s): Dimitris Karlis and Ioannis Ntzoufras
# Source: Journal of the Royal Statistical Society. Series D (The Statistician) , 2003, Vol.
# 52, No. 3 (2003), pp. 381-393

# Modelling Association Football Scores and Inefficiencies in the Football Betting Market
# Mark J. Dixon,Stuart G. Coles

library(ggplot2)
library(tidyr)
library(dplyr)
library(alabama)
library(tidyverse)

DC_model =function(data1,true_data,week,xi1,strt_date){
  
  # Read Premier League 2021/2022 Data
  data = read.csv(data1)
  
  # Looking at our data
  head(data)
  
  # Looking at the mean of our Home and Away Goals
  mean(data$Home.Goals)
  mean(data$Away.Goals)
  
  # Plotting our Home and Away Goals
  hist(data$Home.Goals)
  hist(data$Away.Goals)
  
  # Tau dependance function
  tau = function(x, y, lambda, mu, rho){
    if (x == 0 & y == 0){
      return(1 - (lambda*mu*rho))
    } 
    else if (x == 0 & y == 1){
      return(1 + (lambda*rho))
    } 
    else if (x == 1 & y == 0){
      return(1 + (mu*rho))
    } 
    else if (x == 1 & y == 1){
      return(1 - rho)
    } 
    else {
      return(1)
    }
  }
  
  
  # Our log likelihood function
  pois_log_lik = function(parameters, home_goals, away_goals, home_teams, away_teams, weights = NULL)
  {
    param = relist(parameters, equal_parameters)
    home_log_lik = away_log_lik = tau_log_lik = c()
    att_pow = param$att_pow
    def_pow = param$def_pow
    home_adv = param$home_adv
    rho = param$rho
    
    for (i in 1:length(home_goals))
    {
      exp_goals_home = exp(att_pow[home_teams[i]] + def_pow[away_teams[i]] + home_adv)
      exp_goals_away = exp(att_pow[away_teams[i]] + def_pow[home_teams[i]])
      
      tau_log_lik[i] = log(tau(home_goals[i], away_goals[i], exp_goals_home, exp_goals_away, rho))
      home_log_lik[i] = dpois(home_goals[i], exp_goals_home, log = TRUE)
      away_log_lik[i] = dpois(away_goals[i], exp_goals_away, log = TRUE)
    }
    if (is.null(weights)){
      return (-sum(home_log_lik+away_log_lik))
    }
    else{
      return (-sum((home_log_lik+away_log_lik+tau_log_lik)*weights))
    }
  }
  
  # Attacking Constraint Function
  att_constr = function(parameters, home_goals, away_goals, home_teams, away_teams, weights = NULL)
  {
    n = length(equal_parameters$att_pow)
    param = relist(parameters, equal_parameters)
    attack = param$att_pow
    return((sum(attack)/n) - 1)
  }
  
  
  
  
  # Exponential Weighting Day
  Exp_Weights_Day = function(dates, curr_day, xi = 0){
    date_diff = as.Date(dates, "%d/%m/%Y") - as.Date(curr_day,"%d/%m/%Y")
    date_diff = as.numeric(date_diff*-1)
    phi = exp(-xi*date_diff)
    phi[date_diff <= 0] = 0
    return(phi)
  }
  
  dates = data$Date
  
  
  # Init. our data
  
  home_goals = data$Home.Goals
  away_goals = data$Away.Goals
  teams = unique(data$Home)
  home_teams = data$Home
  away_teams = data$Away
  
  # Our initial parameters (0 for att/def pow and 2 for home adv) 
  equal_parameters <- list(
    att_pow = rep(0.01, length(teams)) %>% `names<-`(teams[1:length(teams)]),
    def_pow = rep(-0.08, length(teams)) %>% `names<-`(teams[1:length(teams)]),
    home_adv = 0.06,
    rho = 0.03
  )
  xi = xi1
  weights = Exp_Weights_Day(dates, strt_date, xi)
  # Calculating MLE, optim will search for optimatl parameters
  
  # mle_fit = optim(par = unlist(equal_parameters),
  #                  fn = pois_log_lik,
  #                  home_goals = home_goals, away_goals = away_goals,
  #                  home_teams = home_teams, away_teams = away_teams, weights = weights, method = "BFGS")
  
  # Converting our parameters back from log() base
  
  # true_mle_fit = exp(mle_fit$par)
  # true_mle = relist(true_mle_fit, equal_parameters)
  
  aug_mle_fit = auglag(par = unlist(equal_parameters), fn = pois_log_lik, heq = att_constr, home_goals = home_goals, away_goals = away_goals,
                       home_teams = home_teams, away_teams = away_teams, weights = weights)
  aug_true_mle = exp(aug_mle_fit$par)
  true_mle = relist(aug_true_mle, equal_parameters)
  # Win Prediction
  
  # Prediction match results using double Poisson
  comp_data = read.csv(true_data)
  games_predict = comp_data %>% filter(Wk == week)
  len = length(games_predict$Home)
  probability_matrix = data.frame(Match = 1:len, Home_Win = len, Draw = len, Away_Win = len, True_Outcome = len)
  for (i in seq_len(len)){
    temp = matrix(ncol=10,nrow = 10)
    atth = true_mle$att_pow[games_predict$Home[i]]
    defh = true_mle$def_pow[games_predict$Home[i]]
    atta = true_mle$att_pow[games_predict$Away[i]]
    defa = true_mle$def_pow[games_predict$Away[i]]
    home_adv = true_mle$home_adv
    rho = log(true_mle$rho)
    for (k in 0:9)
    {
      for (l in 0:9)
        
      {
        temp[k+1,l+1] = dpois(k, atth*defa*home_adv)*dpois(l,atta*defh)
      }
    }
    if (i == 1){
      write.table(temp, file = "Matrix.txt", append = TRUE, sep = ",", dec = ".",
                  row.names = TRUE, col.names = TRUE)
    }
   
    temp[1,1] = temp[1,1]*tau(0,0,atth*defa*home_adv, atta*defh, rho)
    temp[2,1] = temp[2,1]*tau(1,0,atth*defa*home_adv, atta*defh, rho)
    temp[1,2] = temp[1,2]*tau(0,1,atth*defa*home_adv, atta*defh, rho)
    temp[2,2] = temp[2,2]*tau(1,1,atth*defa*home_adv, atta*defh, rho)
    probability_matrix$Match[i] = paste(games_predict$Home[i],"-",games_predict$Away[i])
    probability_matrix$Home_Win[i] = sum(temp[lower.tri(temp)], FALSE)
    probability_matrix$Away_Win[i] = sum(temp[upper.tri(temp)], FALSE)
    probability_matrix$Draw[i] = sum(diag(temp))
    if (games_predict$Home.Goals[i] > games_predict$Away.Goals[i]){
      probability_matrix$True_Outcome[i] = "HOME WIN"
    }
    else if (games_predict$Home.Goals[i] < games_predict$Away.Goals[i]){
      probability_matrix$True_Outcome[i] = "AWAY WIN"
    }
    else{
      probability_matrix$True_Outcome[i] = "DRAW"
    }
  }
  
  
  # Log score
  log_score = 0
  for (i in 1:len){
    if (probability_matrix$True_Outcome[i] == "HOME WIN"){
      log_score = log_score + log(probability_matrix$Home_Win[i])
    }
    if (probability_matrix$True_Outcome[i] == "AWAY WIN"){
      log_score = log_score + log(probability_matrix$Away_Win[i])
    }
    if (probability_matrix$True_Outcome[i] == "DRAW"){
      log_score = log_score + log(probability_matrix$Draw[i])
    }
  }
  
  
  # Xi plot
  library(reshape2)
  date_plt = strt_date
  df_w = data.frame(Date = as.Date(dates, "%d/%m/%Y"), "xi" = weights)
  df_w$"0.0001" = Exp_Weights_Day(dates, date_plt, 0.0001)
  df_w$"0.01" = Exp_Weights_Day(dates, date_plt, 0.01)
  df_w$"0" = Exp_Weights_Day(dates, date_plt, 0)
  df_test = melt(df_w, id.vars = "Date", variable.name = "Weights")
  ggplot(df_test, aes(x = Date, y = value)) + geom_line(aes(color = Weights))
  
  
  
  write.table(probability_matrix, file = "Results4.csv", append = TRUE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
  
  write.table(log_score, file = "log_score4.csv", append = TRUE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
  
  my_list = list("Attack" = true_mle$att_pow, "Defense" = true_mle$def_pow, 
                 "Home" = true_mle$home_adv, "Rho" = true_mle$rho)
  return(my_list)
}

week1 = DC_model("Cross_Valid.csv","18-19.csv",1,0.00325,"10/08/2018")
week2 = DC_model("Cross_Valid.csv","18-19.csv",2,0.004,"18/08/2018")
week3 = DC_model("Cross_Valid.csv","18-19.csv",3,0.004,"25/08/2018")
week4 = DC_model("Cross_Valid.csv","18-19.csv",4,0.004,"01/09/2018")
week5 = DC_model("Cross_Valid.csv","18-19.csv",5,0.004,"15/09/2018")
week6 = DC_model("Cross_Valid.csv","18-19.csv",6,0.004,"22/09/2018")
week7 = DC_model("Cross_Valid.csv","18-19.csv",7,0.004,"29/09/2018")
week8 = DC_model("Cross_Valid.csv","18-19.csv",8,0.004,"05/10/2018")
week9 = DC_model("Cross_Valid.csv","18-19.csv",9,0.004,"20/10/2018")
week10 = DC_model("Cross_Valid.csv","18-19.csv",10,0.004,"27/10/2018")
week11 = DC_model("Cross_Valid.csv","18-19.csv",11,0.004,"03/11/2018")
week12 = DC_model("Cross_Valid.csv","18-19.csv",12,0.004,"10/11/2018")
week13 = DC_model("Cross_Valid.csv","18-19.csv",13,0.004,"24/11/2018")
week14 = DC_model("Cross_Valid.csv","18-19.csv",14,0.00325,"30/11/2018")
week15 = DC_model("Cross_Valid.csv","18-19.csv",15,0.004,"04/12/2018")
week16 = DC_model("Cross_Valid.csv","18-19.csv",16,0.004,"08/12/2018")
week17 = DC_model("Cross_Valid.csv","18-19.csv",17,0.004,"15/12/2018")
week18 = DC_model("Cross_Valid.csv","18-19.csv",18,0.004,"21/12/2018")
week19 = DC_model("Cross_Valid.csv","18-19.csv",19,0.004,"26/12/2018")




# Plotting Teams Parameters
require("ggrepel")
teams = read.csv('18-19.csv')
teams = unique(teams$Home)
team_param = data.frame(Attack = 1:20, Defense = 1:20)
row.names(team_param) = teams
ctr = 0
for (i in teams){
  ctr = ctr + 1
  team_param$Attack[ctr] = week19$Attack[i]
  team_param$Defensive[ctr] = week19$Defense[i]
}
ggplot(data = team_param, aes(x = Defensive, y = Attack)) + geom_point(alpha = 0.5) +
  geom_text_repel(label = teams ) + theme_minimal() + ggtitle("Team Attack and Defensive Parameters") + 
  xlab("Defense Parameter") + ylab("Attack Parameter")

# Sum of log score
log_data = read.csv('log_score4.csv')
log_score = sum(as.numeric(log_data$log))

a = read.csv('Results3.csv')
