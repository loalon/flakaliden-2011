library(dplyr)


data$treatment <- substr(colnames(data), 4,4) #ifelse(grepl('C',(colnames(data))), 'C','F')

data.t <- as.data.frame(t(data))
data.t$treatment <- substr(colnames(data), 4,4)
data.t$time <- substr(colnames(data), 2,3)

#minidata <- data.t [1:12, c(1:50, 311102, 311103)]
# dim(data.t) 
# data.t$treatment <- substr(colnames(data), 4,4)[1:30]
# data.t$time <- substr(colnames(data), 2,3)[1:30]
#ifelse(grepl('C',(colnames(data))), 'C','F')
byTreatment <- data.t %>%
  group_by(treatment)  %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)


byTime <- data.t %>%
  group_by(time)  %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

treatNames <- byTreatment$treatment
timeNames <- byTime$time

tempTime <- byTime[,-1]
tempTreatment <- byTreatment[,-1]

byTime.t <-data.frame( t(tempTime))
colnames(byTime.t) <- timeNames

byTreatment.t <- data.frame(t(tempTreatment))
colnames(byTreatment.t) <- treatNames

data$treatment <- ifelse(byTreatment.t$C>byTreatment.t$F, "C", "F")
data$time <- colnames(byTime.t)[apply(byTime.t, 1, which.max)]

dataPrime <- data
temp <- byTime.t[1:100,]
colnames(temp)[apply(temp, 1, which.max)]
