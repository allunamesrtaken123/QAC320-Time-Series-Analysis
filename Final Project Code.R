library(readr)
library(dplyr)


Matches_Odds <- read_csv("~/COURSES/QAC320/Matches_Odds.csv", col_types = 
                           cols(date_created = col_datetime(format ="%m/%d/%Y %H:%M")))

obs_counts <- Matches_Odds %>% 
  count(match_id) %>% 
  filter(n >= 43) %>% 
  select(match_id)

moreMatches <- Matches_Odds %>% 
  filter(match_id %in% obs_counts$match_id) %>% 
  select(match_id,
         home_team_odd,
         away_team_odd)
moreMatchesFull <- Matches_Odds %>% filter(match_id %in% obs_counts$match_id)

# Get data into proper format to be clustered
widen <- function(d,col_name,data_col){
  n <- 0
  match_ids <- unique(d %>% select(col_name))
  for (e in 1:nrow(match_ids)){
    n <- max(n,nrow(d %>% filter(d[col_name] == match_ids[[1]][e])))
  }
  widen_help <- function(d,col_name,data_col,n){
    i <- unique(d %>% select(col_name))[[1]][1]
    t <- unlist(d %>% filter(d[col_name] == i) %>% select(data_col))
    nas <- rep(NA,n - length(t))
    t <- data.frame(c(t,nas))
    colnames(t) <- paste0(data_col,i)
    if (i == max(unique(d %>% select(col_name)))){return(t)
    }else{return(cbind(t,widen_help(d %>% filter(d[col_name] != i),col_name,data_col,n)))}
  }
  return(widen_help(d,col_name,data_col,n))
}

moreHomeWidened <- widen(moreMatches,'match_id','home_team_odd')
moreAwayWidened <- widen(moreMatches,'match_id','away_team_odd')
moreWidened <- cbind(moreHomeWidened,moreAwayWidened)
moreWidenedT <- t(moreWidened)

library(cluster)
library(factoextra)
# Pairwise distances of time series
dist2 <- cluster::daisy(moreWidenedT,metric = "gower")

# Silhouette Plot for Optimal k
fviz_nbclust(moreWidenedT, pam, method = "silhouette") + theme_classic()

# Cluster!
k2 <- pam(x = dist2, k = 2)
k4 <- pam(x = dist2, k = 4)


# Test Type-Cluster Overlap ----

empty <- data.frame(matrix(ncol = 8,nrow=0))
colnames(empty) <- c('match_id','home','away','max','min',
                     'type','max_val','min_val')
f2 <- function(d,acc){
  result <- data.frame(matrix(ncol = 8,nrow=4))
  colnames(result) <- c('match_id','home','away','max','min',
                        'type','max_val','min_val')
  n = nrow(d)
  
  result$away <- c(0,0,1,1)
  result$home <- c(1,1,0,0)
  
  result$match_id <- head(d$match_id,1)
  
  max_home_index <- 1
  min_home_index <- 1
  max_away_index <- 1
  min_away_index <- 1
  
  max_home <- d[1,2]
  min_home <- d[1,2]
  max_away <- d[1,3]
  min_away <- d[1,3]
  
  for (i in 1:n){
    if (d[i,2] > max_home){max_home <- d[i,2] ; max_home_index <- i
    }else if (d[i,2] < min_home){min_home <- d[i,2] ; min_home_index <- i
    }else if (d[i,3] > max_away){max_away <- d[i,3] ; max_away_index <- i
    }else if (d[i,3] < min_away){min_away <- d[i,3] ; min_away_index <- i}
  }
  
  result$max <- c(if(max_home_index != 1 & max_home_index != n){1}else{0},
                  if(max_home_index != 1 & max_home_index != n){1}else{0},
                  if(max_away_index != 1 & max_away_index != n){1}else{0},
                  if(max_away_index != 1 & max_away_index != n){1}else{0})
  
  result$min <- c(if(min_home_index != 1 & min_home_index != n){1}else{0},
                  if(min_home_index != 1 & min_home_index != n){1}else{0},
                  if(min_away_index != 1 & min_away_index != n){1}else{0},
                  if(min_away_index != 1 & min_away_index != n){1}else{0})
  
  result$max_val <- c(max_home,max_home,max_away,max_away)
  result$min_val <- c(min_home,min_home,min_away,min_away)
  
  for (i in 1:4){
    if (result[i,]$max == 1){
      if (result[i,]$min == 1){result[i,]$type = "type 1"}else{result[i,]$type = "type 2"}
    }else{
      if (result[i,]$min == 1){result[i,]$type = "type 3"}else{result[i,]$type = "type 4"}
    }
  }
  
  return(rbind(acc,result))
}
multi_f2 <- function(d,result){
  for (ids in unique(d$match_id)){
    result <- f2(d %>% filter(match_id == ids),result)
  }
  return(result)
}


clus2 <- unname(k2$clustering)
clus4 <- unname(k4$clustering)
more_type <- multi_f2(moreMatches,empty)
type_vec <- rbind(more_type %>% filter(home==1) %>% 
                    distinct(match_id,.keep_all = T) %>% select(type),
                  more_type %>% filter(away==1) %>% 
                    distinct(match_id,.keep_all = T) %>% select(type))

name_vec <- c()
for (i in more_type$match_id){name_vec <- c(name_vec,paste("home",i))}
for (i in more_type$match_id){name_vec <- c(name_vec,paste("away",i))}

type_cluster2 <- cbind(name_vec,clus2,type_vec)
type_cluster4 <- cbind(name_vec,clus4,type_vec)

type_cluster2 <- cbind(name_vec,clus2,type_vec)
table(type_cluster2$clus2,type_cluster2$type)

type_cluster4 <- cbind(name_vec,clus4,type_vec)
table(type_cluster4$clus4,type_cluster4$type)


# Medoids (look at k2$medoids and k4$medoids to find the medoids) ----
c1mk2 <- Matches_Odds %>% filter(match_id == 700309) %>% select(home_team_odd)
c2mk2 <- Matches_Odds %>% filter(match_id == 659942) %>% select(away_team_odd)

c1mk4 <- Matches_Odds %>% filter(match_id == 687824) %>% select(away_team_odd)
c2mk4 <- Matches_Odds %>% filter(match_id == 690485) %>% select(home_team_odd)
c3mk4 <- Matches_Odds %>% filter(match_id == 125609) %>% select(away_team_odd)
c4mk4 <- Matches_Odds %>% filter(match_id == 690559) %>% select(away_team_odd)

# Find (And Report) the Elements Furthest From the Medoid for Each Cluster ----
# Done for k = 2,4

# k = 2
clustsk2 <- data.frame(k2$clustering)
c <- 1
for (medoid in k2$medoids){
  ids <- clustsk2 %>% filter(k2.clustering == c) %>% row.names()
  max_dist <- 0
  min_dist <- 100000000000000
  max_test = ""
  min_test = ""
  dists <- as.matrix(dist2)[medoid,]
  for (j in ids){
    if (max_dist < dists[j]){max_dist <- dists[j] ; max_test <- j}
    if (min_dist > dists[j] && dists[j] > 0){min_dist <- dists[j] ; min_test <- j}
  }
  print(paste(paste("Cluster:",c),paste("Furthest Object:",max_test),sep = "    " ))
  print(paste(paste("Cluster:",c),paste("Closest Object:",min_test),sep = "    " ))
  c <- c + 1
}

c1fk2 <- Matches_Odds %>% filter(match_id == 688364) %>% select(home_team_odd)
c2fk2 <- Matches_Odds %>% filter(match_id == 666823) %>% select(away_team_odd)

c1ck2 <- Matches_Odds %>% filter(match_id == 691952) %>% select(home_team_odd)
c2ck2 <- Matches_Odds %>% filter(match_id == 692618) %>% select(away_team_odd)

# k = 4
clustsk4 <- data.frame(k4$clustering)
c <- 1
for (medoid in k4$medoids){
  ids <- clustsk4 %>% filter(k4.clustering == c) %>% row.names()
  max_dist <- 0
  min_dist <- 10000000000000
  max_test = ""
  min_test = ""
  dists <- as.matrix(dist2)[medoid,]
  for (j in ids){
    if (max_dist < dists[j]){max_dist <- dists[j] ; max_test <- j}
    if (min_dist > dists[j] && dists[j] > 0){min_dist <- dists[j] ; min_test <- j}
  }
  print(paste(paste("Cluster:",c),paste("Furthest Object:",max_test),sep = "    " ))
  print(paste(paste("Cluster:",c),paste("Closest Object:",min_test),sep = "    " ))
  c <- c + 1
}

c1fk4 <- Matches_Odds %>% filter(match_id == 688364) %>% select(home_team_odd)
c2fk4 <- Matches_Odds %>% filter(match_id == 666823) %>% select(home_team_odd)
c3fk4 <- Matches_Odds %>% filter(match_id == 33576) %>% select(home_team_odd)
c4fk4 <- Matches_Odds %>% filter(match_id == 666823) %>% select(away_team_odd)

c1ck4 <- Matches_Odds %>% filter(match_id == 688525) %>% select(away_team_odd)
c2ck4 <- Matches_Odds %>% filter(match_id == 691328) %>% select(home_team_odd)
c3ck4 <- Matches_Odds %>% filter(match_id == 659942) %>% select(away_team_odd)
c4ck4 <- Matches_Odds %>% filter(match_id == 665256) %>% select(away_team_odd)



max_obs <- moreMatchesFull %>% count(match_id, sort = T) %>% select(n) %>% head(1) %>% pull()
# Naive Methods For k = 2 ------------------------------------------------------
library(forecast)
# Naive Methods for Cluster One of k=2
meanmc1mk2 <- meanf(c1mk2$home_team_odd, h = max_obs)
naivemc1mk2 <- naive(c1mk2$home_team_odd, h = max_obs)
driftmc1mk2 <- rwf(c1mk2$home_team_odd, h = max_obs, drift = T)

sprintf("RMSE of Mean Model (cluster 1 of k=2) on Far Point: %f",accuracy(meanmc1mk2,c1fk2$home_team_odd)[2,2])
sprintf("RMSE of Mean Model (cluster 1 of k=2) on Close Point: %f",accuracy(meanmc1mk2,c1ck2$home_team_odd)[2,2])
sprintf("RMSE of Naive Model (cluster 1 of k=2) on Far Point: %f",accuracy(naivemc1mk2,c1fk2$home_team_odd)[2,2])
sprintf("RMSE of Naive Model (cluster 1 of k=2) on Close Point: %f",accuracy(naivemc1mk2,c1ck2$home_team_odd)[2,2])
sprintf("RMSE of Drift Model (cluster 1 of k=2) on Far Point: %f",accuracy(driftmc1mk2,c1fk2$home_team_odd)[2,2])
sprintf("RMSE of Drift Model (cluster 1 of k=2) on Close Point: %f",accuracy(driftmc1mk2,c1ck2$home_team_odd)[2,2])

# Naive Methods for Cluster Two of k=2
meanmc2mk2 <- meanf(c2mk2$away_team_odd, h = max_obs)
naivemc2mk2 <- naive(c2mk2$away_team_odd, h = max_obs)
driftmc2mk2 <- rwf(c2mk2$away_team_odd, h = max_obs, drift = T)

sprintf("RMSE of Mean Model (cluster 2 of k=2) on Far Point: %f",accuracy(meanmc2mk2,c2fk2$away_team_odd)[2,2])
sprintf("RMSE of Mean Model (cluster 2 of k=2) on Close Point: %f",accuracy(meanmc2mk2,c2ck2$away_team_odd)[2,2])
sprintf("RMSE of Naive Model (cluster 2 of k=2) on Far Point: %f",accuracy(naivemc2mk2,c2fk2$away_team_odd)[2,2])
sprintf("RMSE of Naive Model (cluster 2 of k=2) on Close Point: %f",accuracy(naivemc2mk2,c2ck2$away_team_odd)[2,2])
sprintf("RMSE of Drift Model (cluster 2 of k=2) on Far Point: %f",accuracy(driftmc2mk2,c2fk2$away_team_odd)[2,2])
sprintf("RMSE of Drift Model (cluster 2 of k=2) on Close Point: %f",accuracy(driftmc2mk2,c2ck2$away_team_odd)[2,2])

# Naive Methods For k = 4 ------------------------------------------------------
# Naive Methods for Cluster One of k=4
meanmc1mk4 <- meanf(c1mk4$away_team_odd, h = max_obs)
naivemc1mk4 <- naive(c1mk4$away_team_odd, h = max_obs)
driftmc1mk4 <- rwf(c1mk4$away_team_odd, h = max_obs, drift = T)

# Naive Methods for Cluster Two of k=4
meanmc2mk4 <- meanf(c2mk4$home_team_odd, h = max_obs)
naivemc2mk4 <- naive(c2mk4$home_team_odd, h = max_obs)
driftmc2mk4 <- rwf(c2mk4$home_team_odd, h = max_obs, drift = T)

# Naive Methods for Cluster Three of k=4
meanmc3mk4 <- meanf(c3mk4$away_team_odd, h = max_obs)
naivemc3mk4 <- naive(c3mk4$away_team_odd, h = max_obs)
driftmc3mk4 <- rwf(c3mk4$away_team_odd, h = max_obs, drift = T)

# Naive Methods for Cluster Four of k=4
meanmc4mk4 <- meanf(c4mk4$away_team_odd, h = max_obs)
naivemc4mk4 <- naive(c4mk4$away_team_odd, h = max_obs)
driftmc4mk4 <- rwf(c4mk4$away_team_odd, h = max_obs, drift = T)

# Spline for k = 2 Medoids -----------------------------------------------------

# Cluster 1
fit.splinefc1mk2 <- splinef(c1mk2$home_team_odd, h = max_obs)

sprintf("RMSE of Spline Model (cluster 1 of k=2) on Far Point: %f",accuracy(fit.splinefc1mk2,c1fk2$home_team_odd)[2,2])
sprintf("RMSE of Spline Model (cluster 1 of k=2) on Close Point: %f",accuracy(fit.splinefc1mk2,c1ck2$home_team_odd)[2,2])



# Cluster 2
fit.splinefc2mk2 <- splinef(c2mk2$away_team_odd, h = max_obs)
casts.splfc2mk2 <- forecast(fit.splinefc2mk2)
plot(casts.splfc2mk2)

sprintf("RMSE of Spline Model (cluster 2 of k=2) on Far Point: %f",forecast::accuracy(fit.splinefc2mk2,c2fk2$away_team_odd)[2,2])
sprintf("RMSE of Spline Model (cluster 2 of k=2) on Close Point: %f",forecast::accuracy(fit.splinefc2mk2,c2ck2$away_team_odd)[2,2])

# Spline for k = 4 Medoids------------------------------------------------------

# Cluster 1
fit.splinefc1mk4 <- splinef(c1mk4, h = max_obs)
casts.splfc1mk4 <- forecast(fit.splinefc1mk4)
plot(casts.splfc1mk4)

# Cluster 2
fit.splinefc2mk4 <- splinef(c2mk4, h = max_obs)
casts.splfc2mk4 <- forecast(fit.splinefc2mk4)
plot(casts.splfc2mk4)

# Cluster 3
fit.splinefc3mk4 <- splinef(c3mk4, h = max_obs)
casts.splfc3mk4 <- forecast(fit.splinefc3mk4)
plot(casts.splfc3mk4)

# Cluster 4
fit.splinefc4mk4 <- splinef(c4mk4, h = max_obs)
casts.splfc4mk4 <- forecast(fit.splinefc4mk4)
plot(casts.splfc4mk4)


# Getting each Medoid to be a Time Series
c1mk2ts <- ts(c1mk2)
c2mk2ts <- ts(c2mk2)

c1mk4ts <- ts(c1mk4)
c2mk4ts <- ts(c2mk4)
c3mk4ts <- ts(c3mk4)
c4mk4ts <- ts(c4mk4)


# NNETAR For k = 2 Medoids -----------------------------------------------------
library(forecast)
nnetar.c1mk2 <- nnetar(diff(c1mk2$home_team_odd))
nnetar.c2mk2 <- nnetar(diff(c2mk2$away_team_odd))

accuracy(nnetar(c1fk2$home_team_odd,model= nnetar.c1mk2))
accuracy(nnetar(diff(c1fk2$home_team_odd),model= nnetar.c1mk2))

accuracy(nnetar(diff(c1ck2$home_team_odd),model= nnetar.c1mk2))

accuracy(nnetar(diff(c2fk2$away_team_odd),model= nnetar.c2mk2))
accuracy(nnetar(diff(c2ck2$away_team_odd),model= nnetar.c2mk2))
# NNETAR For k = 4 Medoids -----------------------------------------------------
nnetar.c1mk4 <- nnetar(diff(c1mk4$away_team_odd))
accuracy(nnetar(diff(c1fk4$home_team_odd),model= nnetar.c1mk4))
accuracy(nnetar(diff(c1ck4$away_team_odd),model= nnetar.c1mk4))


nnetar.c2mk4 <- nnetar(diff(c2mk4$home_team_odd))
accuracy(nnetar(diff(c2fk4$home_team_odd),model= nnetar.c2mk4))
accuracy(nnetar(diff(c2ck4$home_team_odd),model= nnetar.c2mk4))

nnetar.c3mk4 <- nnetar(diff(c3mk4$away_team_odd))
accuracy(nnetar(diff(c3fk4$home_team_odd),model= nnetar.c3mk4))
accuracy(nnetar(diff(c3ck4$away_team_odd),model= nnetar.c3mk4))

nnetar.c4mk4 <- nnetar(diff(c4mk4$away_team_odd))
accuracy(nnetar(diff(c4fk4$away_team_odd),model= nnetar.c4mk4))
accuracy(nnetar(diff(c4ck4$away_team_odd),model= nnetar.c4mk4))



# ARIMA/ARFIMA For k = 2 Medoids-------------------------------------------------------
library(tseries)
library(forecast)
# Cluster 1
adf.test(c1mk2ts)
BoxCox.lambda(c1mk2$home_team_odd)
c1mk2tsL <- BoxCox(c1mk2$home_team_odd,BoxCox.lambda(c1mk2$home_team_odd))


c1mk2ts010 <- arima(c1mk2tsL, order = c(0,1,0))
c1mk2ts010
adf.test(c1mk2ts010$residuals) # Differencing makes it stationary
acf(c1mk2ts010$residuals)
pacf(c1mk2ts010$residuals)

c1mk2ts011 <- arima(c1mk2tsL, order = c(0,1,1))
c1mk2ts011
acf(c1mk2ts011$residuals)
pacf(c1mk2ts011$residuals)

c1mk2ts012 <- arima(c1mk2tsL, order = c(0,1,2))
c1mk2ts012
# Turns out order c(0,1,2) is the best (AIC is 18 lower than (0,1,1))
auto.arima(c1mk2tsL) 
AIC(arfima(c1mk2tsL))

# Cluster 2 (Looks like ARIMA is NOT the way to go on this one)
adf.test(c2mk2ts) # Not stationary (p > .05)
BoxCox.lambda(c2mk2ts, upper = 10) # Lambda is 4.8

#Box Cox doesn't make it stationary but closer
c2mk2tsL <- BoxCox(c2mk2ts,BoxCox.lambda(c2mk2ts, upper = 10))

acf(c2mk2tsL) # Shows trend
pacf(c2mk2tsL) # 1 significant lag

c2mk2tsL010 <- arima(c2mk2tsL, order = c(0,1,0))
acf(c2mk2tsL010$residuals)
pacf(c2mk2tsL010$residuals)

adf.test(c2mk2tsL010$residuals) #It's stationary but AIC is high
auto.arima(c2mk2tsL) # This says c(0,1,0) is in fact the best you can do
AIC(arfima(c2mk2tsL))

# ARIMA For k = 4 Medoids-------------------------------------------------------

# Cluster 1
adf.test(c1mk4$away_team_odd)
BoxCox.lambda(c1mk4$away_team_odd, lower = -10) # Lambda is -6.725ish
c1mk4tsL <- BoxCox(c1mk4$away_team_odd, BoxCox.lambda(c1mk4$away_team_odd, lower = -10))
adf.test(c1mk4tsL)

acf(c1mk4tsL)
pacf(c1mk4tsL)
c1mk4tsL010 <- arima(c1mk4tsL, order = c(0,1,0))
acf(c1mk4tsL010$residuals)
pacf(c1mk4tsL010$residuals)

c1mk4tsL011 <- arima(c1mk4tsL, order = c(0,1,1))
c1mk4tsL011
acf(c1mk4tsL011$residuals)
pacf(c1mk4tsL011$residuals)

c1mk4tsL110 <- arima(c1mk4tsL, order = c(1,1,0))
c1mk4tsL110
acf(c1mk4tsL110$residuals)
pacf(c1mk4tsL110$residuals)

c1mk4tsL011100 <- arima(c1mk4tsL, order = c(0,1,1), 
                        seasonal = list(order = c(1,0,0), period = 10))
acf(c1mk4tsL011100$residuals)
pacf(c1mk4tsL011100$residuals)
c1mk4tsL011100
auto.arima(c1mk4tsL) # Turns out that order c(0,1,1) with no seasonality is best

AIC(arfima(c1mk4tsL)) # ARFIMA BETTER!
arfima(c1mk4tsL)
library(forecast)


accuracy(forecast::Arima(c1fk4$home_team_odd,model = c1mk4tsL011))
accuracy(forecast::Arima(c1ck4$away_team_odd,model = c1mk4tsL011))

# Cluster 2
adf.test(c2mk4$home_team_odd)
BoxCox.lambda(c2mk4$home_team_odd)
c2mk4tsL <- BoxCox(c2mk4$home_team_odd, BoxCox.lambda(c2mk4$home_team_odd))
adf.test(c2mk4tsL)

acf(c2mk4tsL)
pacf(c2mk4tsL)
c2mk4tsL010 <- arima(c2mk4tsL, order = c(0,1,0))
c2mk4tsL010
acf(c2mk4tsL010$residuals)
pacf(c2mk4tsL010$residuals)

c2mk4tsL011 <- arima(c2mk4tsL, order = c(0,1,1))
c2mk4tsL011
acf(c2mk4tsL011$residuals)
pacf(c2mk4tsL011$residuals)

c2mk4tsL110 <- arima(c2mk4tsL, order = c(1,1,0))
c2mk4tsL110
acf(c2mk4tsL110$residuals)
pacf(c2mk4tsL110$residuals)

c2mk4tsL111 <- arima(c2mk4tsL, order = c(1,1,1)) #This one does the best
c2mk4tsL111
acf(c2mk4tsL111$residuals)
pacf(c2mk4tsL111$residuals)

accuracy(forecast::Arima(c2ck4$home_team_odd,model = c2mk4tsL111))
accuracy(forecast::Arima(c2fk4$home_team_odd,model = c2mk4tsL111))

AIC(arfima(c2mk4tsL))

# Cluster 3
adf.test(c3mk4$away_team_odd) # Not stationary
BoxCox.lambda(c3mk4$away_team_odd)
c3mk4tsL <- BoxCox(c3mk4$away_team_odd,BoxCox.lambda(c3mk4$away_team_odd))
adf.test(c3mk4tsL) #Still not stationary

acf(c3mk4tsL)
pacf(c3mk4tsL)

c3mk4tsL010 <- arima(c3mk4tsL, order = c(0,1,0))
c3mk4tsL010
acf(c3mk4tsL010$residuals)
pacf(c3mk4tsL010$residuals)

c3mk4tsL011 <- arima(c3mk4tsL, order = c(0,1,1)) #BEST
c3mk4tsL110 <- arima(c3mk4tsL, order = c(1,1,0))
c3mk4tsL111 <- arima(c3mk4tsL, order = c(1,1,1))
c3mk4tsL012 <- arima(c3mk4tsL, order = c(0,1,2))
AIC(c3mk4tsL011)
AIC(c3mk4tsL110)
AIC(c3mk4tsL111)
AIC(c3mk4tsL012)
auto.arima(c3mk4tsL)
AIC(arfima(c3mk4tsL))

accuracy(forecast::Arima(c3ck4$away_team_odd,model = c3mk4tsL011))
accuracy(forecast::Arima(c3fk4$home_team_odd,model = c3mk4tsL011))

# Cluster 4
adf.test(c4mk4$away_team_odd)
BoxCox.lambda(c4mk4$away_team_odd, lower = -5)
c4mk4tsL <- BoxCox(c4mk4$away_team_odd,BoxCox.lambda(c4mk4$away_team_odd, lower = -5))

acf(c4mk4tsL)
pacf(c4mk4tsL)
c4mk4tsL010 <- arima(c4mk4tsL, order = c(0,1,0))

c4mk4tsL011 <- arima(c4mk4tsL, order = c(0,1,1))
c4mk4tsL110 <- arima(c4mk4tsL, order = c(1,1,0))
c4mk4tsL111 <- arima(c4mk4tsL, order = c(1,1,1))
c4mk4tsL012 <- arima(c4mk4tsL, order = c(0,1,2))
c4mk4tsL013 <- arima(c4mk4tsL, order = c(0,1,3)) #Best ARIMA
c4mk4tsL014 <- arima(c4mk4tsL, order = c(0,1,4))
AIC(c4mk4tsL011)
AIC(c4mk4tsL110)
AIC(c4mk4tsL111)
AIC(c4mk4tsL012)
AIC(c4mk4tsL013)
AIC(c4mk4tsL014)
auto.arima(c2mk4tsL)

accuracy(forecast::Arima(c4ck4$away_team_odd,model = c4mk4tsL013))
accuracy(forecast::Arima(c4fk4$away_team_odd,model = c4mk4tsL013))

AIC(arfima(c4mk4tsL))
arfima(c4mk4tsL)
arfima
#ARFIMA ----

# Cluster 1 k = 2
arfima.c1mk2 <- arfima(c1mk2$home_team_odd)
accuracy(arfima(c1ck2$home_team_odd,model = arfima.c1mk2))
accuracy(arfima(c1fk2$home_team_odd,model = arfima.c1mk2))

# Cluster 2 k = 2
arfima.c2mk2 <- arfima(c2mk2$away_team_odd)
accuracy(arfima(c2ck2$away_team_odd,model = arfima.c2mk2))
accuracy(arfima(c2fk2$away_team_odd,model = arfima.c2mk2))

# Cluster 1 k = 4
arfima.c1mk4 <- arfima(c1mk4$away_team_odd)
accuracy(arfima(c1ck4$away_team_odd,model = arfima.c1mk4))
accuracy(arfima(c1fk4$home_team_odd,model = arfima.c1mk4))

# Cluster 2 k = 4
arfima.c2mk4 <- arfima(c2mk4$home_team_odd)
accuracy(arfima(c2ck4$home_team_odd,model = arfima.c2mk4))
accuracy(arfima(c2fk4$home_team_odd,model = arfima.c2mk4))

# Cluster 3 k = 4
arfima.c3mk4 <- arfima(c3mk4$away_team_odd)
accuracy(arfima(c3ck4$away_team_odd,model = arfima.c3mk4))
accuracy(arfima(c3fk4$home_team_odd,model = arfima.c3mk4))

# Cluster 4 k = 4
arfima.c4mk4 <- arfima(c4mk4$away_team_odd)
accuracy(arfima(c4ck4$away_team_odd,model = arfima.c4mk4))
accuracy(arfima(c4fk4$away_team_odd,model = arfima.c4mk4))

# Evaluating Which is Best Model -----------------------------------------------
# Should I select a few members from each cluster and take average accuracy score?
# Test for ARCH effect in the residuals of the best models!







# Home vs. Away Distribution in Clusters ---------------------------------------


for (i in c(1,2)){
  print(as.data.frame(data.frame(k2$clustering) %>% filter(k2.clustering == i) %>%
                        row.names() %>% substr(1,4) %>% table()))
}

hak2dist<- matrix(c(376,325,22,73),nrow = 2, ncol = 2, 
       dimnames = list(c("home","away"),c("c1","c2")))


for (i in c(1:4)){
  print(as.data.frame(data.frame(k4$clustering) %>% filter(k4.clustering == i) %>%
                        row.names() %>% substr(1,4) %>% table()))
}

hak4dist<- matrix(c(229,193,97,183,59,19,13,3),nrow = 2, ncol = 4, 
                  dimnames = list(c("home","away"),c("c1","c2","c3","c4")))

# Shows most popular leagues in the data set----
distinct(moreMatchesFull,competition_name,match_id) %>% count(competition_name,sort = T)

europa <- moreMatchesFull %>% filter(competition_name == "Europa League") %>%
  distinct(match_id)

for (id in europa$match_id){
  print(paste("home_team_odd",id,sep = ""))
  print(k2$clustering[paste("home_team_odd",id,sep = "")])
  
  print(paste("away_team_odd",id,sep = ""))
  print(k2$clustering[paste("away_team_odd",id,sep = "")])
}

for (id in europa$match_id){
  print(paste("home_team_odd",id,sep = ""))
  print(k4$clustering[paste("home_team_odd",id,sep = "")])
  
  print(paste("away_team_odd",id,sep = ""))
  print(k4$clustering[paste("away_team_odd",id,sep = "")])
}

#RMSE Shit ---------------------------------------------------------------------
library(Metrics)
# Cluster 1 k = 4----
sprintf("RMSE of Mean Model (cluster 1 of k=4) on Far Point: %f",accuracy(meanmc1mk4,c1fk4$home_team_odd)[2,2])
sprintf("RMSE of Mean Model (cluster 1 of k=4) on Close Point: %f",accuracy(meanmc1mk4,c1ck4$away_team_odd)[2,2])
sprintf("RMSE of Naive Model (cluster 1 of k=4) on Far Point: %f",accuracy(naivemc1mk4,c1fk4$home_team_odd)[2,2])
sprintf("RMSE of Naive Model (cluster 1 of k=4) on Close Point: %f",accuracy(naivemc1mk4,c1ck4$away_team_odd)[2,2])
sprintf("RMSE of Drift Model (cluster 1 of k=4) on Far Point: %f",accuracy(driftmc1mk4,c1fk4$home_team_odd)[2,2])
sprintf("RMSE of Drift Model (cluster 1 of k=4) on Close Point: %f",accuracy(driftmc1mk4,c1ck4$away_team_odd)[2,2])
sprintf("RMSE of Spline Model (cluster 1 of k=4) on Far Point: %f",accuracy(fit.splinefc1mk4,c1fk4$home_team_odd)[2,2])
sprintf("RMSE of Spline Model (cluster 1 of k=4) on Close Point: %f",accuracy(fit.splinefc1mk4,c1ck4$away_team_odd)[2,2])
# NNETAR
rmse(c1fk4$home_team_odd,forecast(c1fk4$home_team_odd,h = count(c1fk4)$n,model = nnetar.c1mk4)$mean)
rmse(c1ck4$away_team_odd,forecast(c1ck4$away_team_odd,h = count(c1ck4)$n,model = nnetar.c1mk4)$mean)
#ARIMA
rmse(c1fk4$home_team_odd,forecast(c1fk4$home_team_odd,h = count(c1fk4)$n,model = c1mk4tsL011)$fitted)
rmse(c1ck4$away_team_odd,forecast(c1ck4$away_team_odd,h = count(c1ck4)$n,model = c1mk4tsL011)$fitted)
accuracy(c1mk4tsL011, c1fk4$away_team_odd)
forecast(diff(c1fk2$home_team_odd),h = count(c1fk2)$n,model = nnetar.c1mk2) 
# Gives same values when I do differenced or non-differenced

c1fk2$home_team_odd
#ARFIMA
rmse(c1fk4$home_team_odd,forecast(object = c1fk4$home_team_odd,h = count(c1fk4)$n,model = arfima(c1mk4$away_team_odd))$fitted)
rmse(c1ck4$away_team_odd,forecast(c1ck4$away_team_odd,h = count(c1ck4)$n,model = arfima(c1mk4$away_team_odd))$fitted)
# Cluster 2 k = 4----
sprintf("RMSE of Mean Model (cluster 2 of k=4) on Far Point: %f",accuracy(meanmc2mk4,c2fk4$home_team_odd)[2,2])
sprintf("RMSE of Mean Model (cluster 2 of k=4) on Close Point: %f",accuracy(meanmc2mk4,c2ck4$home_team_odd)[2,2])
sprintf("RMSE of Naive Model (cluster 2 of k=4) on Far Point: %f",accuracy(naivemc2mk4,c2fk4$home_team_odd)[2,2])
sprintf("RMSE of Naive Model (cluster 2 of k=4) on Close Point: %f",accuracy(naivemc2mk4,c2ck4$home_team_odd)[2,2])
sprintf("RMSE of Drift Model (cluster 2 of k=4) on Far Point: %f",accuracy(driftmc2mk4,c2fk4$home_team_odd)[2,2])
sprintf("RMSE of Drift Model (cluster 2 of k=4) on Close Point: %f",accuracy(driftmc2mk4,c2ck4$home_team_odd)[2,2])
sprintf("RMSE of Spline Model (cluster 2 of k=4) on Far Point: %f",accuracy(fit.splinefc2mk4,c2fk4$home_team_odd)[2,2])
sprintf("RMSE of Spline Model (cluster 2 of k=4) on Close Point: %f",accuracy(fit.splinefc2mk4,c2ck4$home_team_odd)[2,2])
# NNETAR
rmse(c2fk4$home_team_odd,forecast(c2fk4$home_team_odd,h = count(c2fk4)$n,model = nnetar.c2mk4)$mean)
rmse(c2ck4$home_team_odd,forecast(c2ck4$home_team_odd,h = count(c2ck4)$n,model = nnetar.c2mk4)$mean)
#ARIMA
rmse(c2fk4$home_team_odd,forecast(c2fk4$home_team_odd,h = count(c2fk4)$n,model = c2mk4tsL111)$fitted)
rmse(c2ck4$home_team_odd,forecast(c2ck4$home_team_odd,h = count(c2ck4)$n,model = c2mk4tsL111)$fitted)

# Cluster 3 k = 4----
sprintf("RMSE of Mean Model (cluster 3 of k=4) on Far Point: %f",accuracy(meanmc3mk4,c3fk4$home_team_odd)[2,2])
sprintf("RMSE of Mean Model (cluster 3 of k=4) on Close Point: %f",accuracy(meanmc3mk4,c3ck4$away_team_odd)[2,2])
sprintf("RMSE of Naive Model (cluster 3 of k=4) on Far Point: %f",accuracy(naivemc3mk4,c3fk4$home_team_odd)[2,2])
sprintf("RMSE of Naive Model (cluster 3 of k=4) on Close Point: %f",accuracy(naivemc3mk4,c3ck4$away_team_odd)[2,2])
sprintf("RMSE of Drift Model (cluster 3 of k=4) on Far Point: %f",accuracy(driftmc3mk4,c3fk4$home_team_odd)[2,2])
sprintf("RMSE of Drift Model (cluster 3 of k=4) on Close Point: %f",accuracy(driftmc3mk4,c3ck4$away_team_odd)[2,2])
sprintf("RMSE of Spline Model (cluster 3 of k=4) on Far Point: %f",accuracy(fit.splinefc3mk4,c3fk4$home_team_odd)[2,2])
sprintf("RMSE of Spline Model (cluster 3 of k=4) on Close Point: %f",accuracy(fit.splinefc3mk4,c3ck4$away_team_odd)[2,2])

# Cluster 4 k = 4----
sprintf("RMSE of Mean Model (cluster 4 of k=4) on Far Point: %f",accuracy(meanmc4mk4,c4fk4$away_team_odd)[2,2])
sprintf("RMSE of Mean Model (cluster 4 of k=4) on Close Point: %f",accuracy(meanmc4mk4,c4ck4$away_team_odd)[2,2])
sprintf("RMSE of Naive Model (cluster 4 of k=4) on Far Point: %f",accuracy(naivemc4mk4,c4fk4$away_team_odd)[2,2])
sprintf("RMSE of Naive Model (cluster 4 of k=4) on Close Point: %f",accuracy(naivemc4mk4,c4ck4$away_team_odd)[2,2])
sprintf("RMSE of Drift Model (cluster 4 of k=4) on Far Point: %f",accuracy(driftmc4mk4,c4fk4$away_team_odd)[2,2])
sprintf("RMSE of Drift Model (cluster 4 of k=4) on Close Point: %f",accuracy(driftmc4mk4,c4ck4$away_team_odd)[2,2])
sprintf("RMSE of Spline Model (cluster 4 of k=4) on Far Point: %f",accuracy(fit.splinefc4mk4,c4fk4$away_team_odd)[2,2])
sprintf("RMSE of Spline Model (cluster 4 of k=4) on Close Point: %f",accuracy(fit.splinefc4mk4,c4ck4$away_team_odd)[2,2])


forecast(c4fk4$away_team_odd,model = c4k4arfima)
# Cluster 1 k = 2 ----
sprintf("RMSE of Spline Model (cluster 1 of k=2) on Far Point: %f",forecast::accuracy(fit.splinefc1mk2,c1fk2$home_team_odd)[2,2])
sprintf("RMSE of Spline Model (cluster 1 of k=2) on Close Point: %f",forecast::accuracy(fit.splinefc1mk2,c1ck2$home_team_odd)[2,2])
rmse(c1fk2$home_team_odd,forecast(c1fk2$home_team_odd,h = count(c1fk2)$n,model = c1mk2ts011)$fitted)
rmse(c1ck2$home_team_odd,forecast(c1ck2$home_team_odd,h = count(c1ck2)$n,model = c1mk2ts011)$fitted)
accuracy(forecast::Arima(c1fk2$home_team_odd,model = c1mk2ts012))
accuracy(forecast::Arima(c1ck2$home_team_odd,model = c1mk2ts012))
# Cluster 2 k = 2 ----
accuracy(forecast::Arima(c2fk2$away_team_odd,model = c2mk2tsL010))
accuracy(forecast::Arima(c2ck2$away_team_odd,model = c2mk2tsL010))



####visuals----------------
require(ggplot2)
require(reshape2)

o <- melt(data.frame(cbind(c(1:nrow(moreWidened)),moreWidened)),
          id.vars = 'c.1.nrow.moreWidened..',variable.name = 'series')
ggplot(o,aes(x = c.1.nrow.moreWidened..,y = value, color = series)) + geom_line()

for (i in c(4)){
  r <- t(as.data.frame(moreWidenedT)[row.names(as.data.frame(moreWidenedT)) %in% row.names(as.data.frame(data.frame(k4$clustering) %>% filter(k4.clustering == i))),])
  r2 <- melt(data.frame(cbind(c(1:nrow(r)),r)),id.vars = 'V1',variable.name = 'series') %>% 
    filter(series != 'away_team_odd667658')
  ggplot(r2,aes(x = V1,y = value,color = series)) + geom_line()
}


library(TSstudio)
# Plotting Predictions + ts for close,far of c1k2 ARIMA 012 Best----

c <- cbind(c1fk2$home_team_odd,forecast::Arima(c1fk2$home_team_odd,model = c1mk2ts012)$fitted)
colnames(c) <- c("Far Actual","Far Predicted")
ts_plot(c,Ytitle = "Odds",Xtitle = "Time")

c <- cbind(c1ck2$home_team_odd,forecast::Arima(c1ck2$home_team_odd,model = c1mk2ts012)$fitted)
colnames(c) <- c("Close Actual","Close Predicted")
ts_plot(c,Ytitle = "Odds",Xtitle = "Time")

# Plotting Predictions + ts for close,far of c2k2 ARFIMA Best----

c <- cbind(c2fk2$away_team_odd,arfima(c2fk2$away_team_odd,model = arfima.c2mk2)$fitted)
colnames(c) <- c("Far Actual","Far Predicted")
ts_plot(c,Ytitle = "Odds",Xtitle = "Time")

c <- cbind(c2ck2$away_team_odd,arfima(c2ck2$away_team_odd,model = arfima.c2mk2)$fitted)
colnames(c) <- c("Close Actual","Close Predicted")
ts_plot(c,Ytitle = "Odds",Xtitle = "Time")

# Plotting Predictions + ts for close,far of c1k4 ARIMA Best----
c <- cbind(c1fk4$home_team_odd,forecast::Arima(c1fk4$home_team_odd,model = c1mk4tsL011)$fitted)
colnames(c) <- c("Far Actual","Far Predicted")
ts_plot(c,Ytitle = "Odds",Xtitle = "Time")

c <- cbind(c1ck4$away_team_odd,forecast::Arima(c1ck4$away_team_odd,model = c1mk4tsL011)$fitted)
colnames(c) <- c("Close Actual","Close Predicted")
ts_plot(c,Ytitle = "Odds",Xtitle = "Time")

# Plotting Predictions + ts for close,far of c2k4 ARIMA Best----
c <- cbind(c2fk4$home_team_odd,forecast::Arima(c2fk4$home_team_odd,model = c2mk4tsL111)$fitted)
colnames(c) <- c("Far Actual","Far Predicted")
ts_plot(c,Ytitle = "Odds",Xtitle = "Time")

c <- cbind(c2ck4$home_team_odd,forecast::Arima(c2ck4$home_team_odd,model = c2mk4tsL111)$fitted)
colnames(c) <- c("Close Actual","Close Predicted")
ts_plot(c,Ytitle = "Odds",Xtitle = "Time")

# Plotting Predictions + ts for close,far of c3k4 ARIMA Best ----
c <- cbind(c3fk4$home_team_odd,forecast::Arima(c3fk4$home_team_odd,model = c3mk4tsL011)$fitted)
colnames(c) <- c("Far Actual","Far Predicted")
ts_plot(c,Ytitle = "Odds",Xtitle = "Time")

c <- cbind(c3ck4$away_team_odd,forecast::Arima(c3ck4$away_team_odd,model = c3mk4tsL011)$fitted)
colnames(c) <- c("Close Actual","Close Predicted")
ts_plot(c,Ytitle = "Odds",Xtitle = "Time")

# Plotting Predictions + ts for close,far of c4k4 ARFIMA Best ----

c <- cbind(c4fk4$away_team_odd,arfima(c4fk4$away_team_odd,model = arfima.c4mk4)$fitted)
colnames(c) <- c("Far Actual","Far Predicted")
ts_plot(c,Ytitle = "Odds",Xtitle = "Time")

c <- cbind(c4ck4$away_team_odd,arfima(c4ck4$away_team_odd,model = arfima.c4mk4)$fitted)
colnames(c) <- c("Close Actual","Close Predicted")
ts_plot(c,Ytitle = "Odds",Xtitle = "Time")