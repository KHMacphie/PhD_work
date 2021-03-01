#using NCtempsSpread
rm(list=ls())
library(dplyr)
library(readr)
setwd('/Users/s1205615/')
NCtempsSpread <- read_csv("~/Downloads/NCtempsSpread.csv")

tempData <- NCtempsSpread %>% 
  filter(year < 2021) #removing 2021

#creating the slide function
slideFunct <- function(data, window, step){   #data= daily temps, window=duration, step=change in start
  total <- length(data) #number of columns in data
  spots <- seq(from=1, to=(total-window), by=step)  #sequence the different starting columns
  result <- vector(length = length(spots))  #mean temp for each window
  for(i in 1:length(spots)){
    result[i] <- mean(data[spots[i]:(spots[i]+window)])  #this makes the windows 1 day longer than the 'window' term (can add -1 if you dont want that to be the case)
  }
  return(result)
}

  #using the slide function for 1992
data <- as.vector (t (tempData[2, 3:368])) # all temps for 1992
window <- 7 #duration (1 less than the number of days it is using)
step <- 1 #change start point by

slideFunct(data, window, step) #vector of mean temp for each window in year 1992


store <- data.frame(NA_col = rep(NA, 52)) #empty df, nrow=number of windows
i <- NA
for (i in 1:length(tempData$year)){ #changed from 52 to the number of years- sorry that might have been my confusion when we chatted
  window <- 7 
  step <- 7
  data <- as.vector(t(tempData[i, 3:368])) # i goes through each year in turn
  store[,i] <- slideFunct(data, window, step) #store mean temps as column
  #rownames(store)[i] <- paste(tempData[i, 2]) #is this meant to be naming rows as start date? this is trying to name them as year but each row now is a window and columns are years
  #names(store) <- c(1:ncol(store)) #is this trying to name the columns? if so needs [i] somewhere or can be done outside of the loop afterwards
}

##It wasnt really working for me with the lines of code for naming things but without them it worked fine
#it didn't reproduce the problem you were having before but I think you changed a few bits while we chatted
#I think with the naming rows and columns you may have it the wrong way around?
#This loop is going through years as columns not windows so the rows should be the window ID 
#(sorry I may have suggested changing the start of the loop when we spoke, I thought it was looping through windows not years)

#I think to name the rows you need spots[i] from the function so doing this without it as a function initially may be easier (or adding more of the whole process to the function if you are feeling fancy).

#Currently each window is ending up as a row, and each column a year- keeping the format where years are rows might be easier to work with
#Though I guess you can just flip the df afterwards!

#some adjustments to this set up to solve row and column name problems. 
store <- data.frame(NA_col = rep(NA, 52)) 
i <- NA
window <- 7 #these 2 rows dont really need to be in the loop if they arent changing throughout (maybe it will for different durations though)
step <- 7
total <- length(tempData[, 3:368]) #these two have been taken from the function because they are relevant to the identity of the window
spots <- seq(from=1, to=(total-window), by=step)

for (i in 1:length(tempData$year)){ #changed from 52 to the number of years
  data <- as.vector(t(tempData[i, 3:368])) #each year in turn
  store[,i] <- slideFunct(data, window, step) #store mean temps as column
}

rownames(store) <- colnames(tempData)[spots+2] #now each row is named the start date for each window- will need to change when different durations are added
colnames(store) <- tempData$year #and the columns are years


##Other thoughts: (we spoke about these the other day but just a reminder)     
#needs to bridge the new year with windows (like marsham set up so -x days to x days)
#ordinal date to resolve leap year problem 
#I think this'd be easier to set up on the marsham data set because you'll be changing this df to be the same format marsham already has
#The duration of the window will need to vary as well 

#Window ranges suggestion 
#Duration: 2 weeks-8weeks by 2 weekly increments 
#Start date: shifting by a week at a time (maybe 2 weeks if that's making a lot of models!)
