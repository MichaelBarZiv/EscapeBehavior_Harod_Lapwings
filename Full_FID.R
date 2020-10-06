
######## 1) Preparing the FID data ####
# This script prepare my data of FID by arranging the predator CSV to the right mesourments and combine it with the human 
# FID to create one file for analyzing.
####
#### 1.1) Open file and library ####

setwd("C:/Users/OrrS4/Dropbox (1)/Orr lab - miki's shared workspace/PhD/FID/predator")


#Open file:
FIDPredator <- read.csv("FIDpredatorR13.csv")
# Labraries:

library(leaflet)
library(tidyverse)
library(ggplot2)
library(ggplot2)
library(lattice)
library(lme4)
library(gplots)
library(multcomp)
library(dplyr)
library(AICcmodavg)#for aictab
library(lubridate) # to work with time
library(rptR)
library(qpcR)
library(sp)
library(AICcmodavg)





#### 1.2) Change column names to the right variables####
# I need to creat new column with the variables i want:
colnames (FIDPredator)

# Rename some of the column names:

names(FIDPredator)[names(FIDPredator) == "EsDme"] <- "EsDMe"    
names(FIDPredator)[names(FIDPredator) == "RuD"] <- "RuDMe"    
names(FIDPredator)[names(FIDPredator) == "SecRuD"] <- "SecRuDMe"    
names(FIDPredator)[names(FIDPredator) == "SecEsD"] <- "SecEsDMe"    

colnames (FIDPredator)


#### 1.3) FID and RuD (will later become both as FirstFlee) and a column for running or flying ####

#First - The actual FID:

FIDPredator$FID <- FIDPredator$StD - FIDPredator$FIDcar

str(FIDPredator)
hist(FIDPredator$FID)
summary(FIDPredator$FID)

minusFID <- subset(FIDPredator, FID < 0) # One row with minus
#The actual Running distance:

FIDPredator$RuD <- FIDPredator$StD - FIDPredator$RuDMe

colnames (FIDPredator)
str(FIDPredator)
hist(FIDPredator$RuD)
summary(FIDPredator$RuD)
minusRuD <- subset(FIDPredator, RuD < 0) # One row with minus

# I need to create a new column "FirstFlee" - that take into considuration the FID, but if there is
# RuD it takes the RuD

FIDPredator$noimportant2 <- ifelse( FIDPredator$RuD >= 0, FIDPredator$RuD, NA)
FIDPredator$notimportant <- ifelse( is.na(FIDPredator$RuD), FIDPredator$FID, NA)

FIDPredator$FirstFlee <- rowSums(FIDPredator[,c("noimportant2", "notimportant")], na.rm=TRUE)*NA^!rowSums(!is.na(FIDPredator[,c("noimportant2", "notimportant")]))
hist(FIDPredator$FirstFlee)
summary(FIDPredator$FirstFlee)
#and removing the column i dont need:

FIDPredator = subset(FIDPredator, select = -c(noimportant2,notimportant) )

# Remove NA from "first flee"
FIDPredator <- FIDPredator[complete.cases(FIDPredator$FirstFlee),]
summary(FIDPredator$FirstFlee)

# Now I need a value of if run by foot or fly:

FIDPredator$RunFly1 <- rep("Fly",nrow(FIDPredator))
FIDPredator$RunFly1[FIDPredator$RuD > 0] <- "Run"
FIDPredator$RunFly1[is.na(FIDPredator$FirstFlee)] <- NA
summary(FIDPredator$RunFly1)

table(FIDPredator$RunFly1)


#### 1.4) Triangular calculation for distance from bird to predator####

#Now i would like to know the escape distance from the car. I will need to know
#the engle that the lapwing flow from the car:

list(FIDPredator$SecDirection)
list(FIDPredator$FirstDirection)

str(FIDPredator$FirstDirection)
str(FIDPredator$SecDirection)

FIDPredator$FirstDirection <- as.numeric(as.character(FIDPredator$FirstDirection))

# Choosing the direction of the lapwing:

#direct1 <- FIDPredator$FirstDirection - 180
#direct2 <- FIDPredator$SecDirection - 180


#directions <- as.data.frame(cbind(direct1,direct2))
#directions <- na.omit(directions) 
#directions$ID <- seq.int(nrow(directions))
#directions$sum <- directions$direct1 - directions$direct2


FIDPredator$FirstDirection - FIDPredator$SecDirection
FIDPredator$SecDirection - FIDPredator$FirstDirection

FIDPredator$FirstAngle <- abs(FIDPredator$FirstDirection - FIDPredator$SecDirection)
colnames (FIDPredator)
summary(FIDPredator$FirstAngle)

FIDPredator$PredaSecDirec <- abs(FIDPredator$PredatorDirection - FIDPredator$SecDirection)
colnames (FIDPredator)
summary(FIDPredator$PredaSecDirec)


#Escape distance = ROOT(EsDMe^2 + FIDcar^2 - 2*FIDcar*EsDMe·cos(FirstAngle))##
# Need to fix it to the new angles
root1 <- FIDPredator$EsDMe^2
root2 <- FIDPredator$FIDcar^2
cos1 <- 2*FIDPredator$EsDMe*FIDPredator$FIDcar*cos(FIDPredator$FirstAngle*pi/180)
FIDPredator$EsD0.1 <- sqrt(root1 + root2 - (cos1))
FIDPredator$SecStD <- FIDPredator$EsD

##FIRST SECFID: 0.1
#Now I have a new tirangle, with EsD, FIDcar, Firstangle. I want to calculate the second FID. 
# The car approach the bird, and I have now SecFIDcar.

#Option 1: calculate with angle of predator:

#ROOT(b2 + c2 - 2bc·cos(A))

root1
root3 <- FIDPredator$SecFIDcar^2
cos2 <- 2*FIDPredator$EsDMe*FIDPredator$SecFIDcar*cos(FIDPredator$PredaSecDirec*pi/180)
FIDPredator$SecFID0.1 <- sqrt(root1 + root3 - (cos2))

##SECOND SECFID: 0.2
# to first I will have to calculate the angle from bird to car and me (i will call it B)
# I will need a = EsD, b = EsDMe, c = FIDcar
# arccos((EsD^2 + EsDMe^2 - FIDcar^2)/(2*EsD*EsDMe))

root4 <- FIDPredator$EsD0.1^2
root1
root2
base <- 2*FIDPredator$EsD0.1*FIDPredator$EsDMe
first1 <- root4 + root1 - root2
arcos <- acos(first1/base)

FIDPredator$SecondAngle <- arcos*57.2957795

# I want to create now the SecFID - will be calculated from 
#SecondAngle, c = SecFIDcar, b = EsDMe.
#a = (c·sin(A))/sin(C)
# C = SecondAngle, c = SecFIDcar, A = ThirdAngle.

# I will need first to know the third angle:
# A = 180 - SecondAngle - forthAngle.... need the forth angle then...
# B = arcsin((b*sin(C)/c) - I need b = EsdMe, C = SecondAngle, c = secFIDcar

sinus <- sin(FIDPredator$SecondAngle*pi/180)

secondes <- FIDPredator$EsDMe*sinus

secondb <- (secondes/FIDPredator$SecFIDcar)

FIDPredator$ForthAngle <- 180 - asin(secondb)*57.2957795

FIDPredator$ThirdAngle <- abs(180 - FIDPredator$ForthAngle - FIDPredator$SecondAngle)

#sROOT(EsDme2 + SecFIDcar2 - 2SecFIDcar*EsDme·cos(A))

root1
root3
cos3 <- 2*FIDPredator$EsDMe*FIDPredator$SecFIDcar*cos(FIDPredator$ThirdAngle*pi/180)
FIDPredator$SecFID0.2 <- sqrt(root1 + root3 - (cos3))

# Now to combine SecFID0.1 and 0.2 to one column, with FID0.1 the dominance:

FIDPredator$noimportant2 <- ifelse( FIDPredator$SecFID0.1 >= 0, FIDPredator$SecFID0.1, NA)
FIDPredator$notimportant <- ifelse( is.na(FIDPredator$SecFID0.1), FIDPredator$SecFID0.2, NA)

FIDPredator$SecFID <- rowSums(FIDPredator[,c("noimportant2", "notimportant")], na.rm=TRUE)*NA^!rowSums(!is.na(FIDPredator[,c("noimportant2", "notimportant")]))

#and removing the column i dont need:

FIDPredator = subset(FIDPredator, select = -c(noimportant2,notimportant, SecFID0.1,SecFID0.2,EsD) )

names(FIDPredator)[names(FIDPredator) == "EsD0.1"] <- "EsD"    

# Now after I have the second FID, I need to take into considuration the second RuD:

cos5 <- 2*FIDPredator$EsDMe*FIDPredator$SecRuDMe*cos(FIDPredator$PredaSecDirec*pi/180)
FIDPredator$SecRuD0.2 <- sqrt(root1 + root3 - (cos5))

cos4 <- 2*FIDPredator$EsDMe*FIDPredator$SecRuDMe*cos(FIDPredator$ThirdAngle*pi/180)
FIDPredator$SecRuD0.1 <- sqrt(root1 + root3 - (cos4))

try <- subset(FIDPredator, select = c("SecRuDMe","SecRuD0.2","SecRuD0.1","SecFID"))

# Only the SecRud0.1 work here - So I will make hime the only column:

FIDPredator$SecRuD <- FIDPredator$SecRuD0.1

FIDPredator = subset(FIDPredator, select = -c(SecRuD0.1,SecRuD0.2 ) )

# Now to combine the two column run 2 and fid 2 to a single column:


FIDPredator$noimportant2 <- ifelse( FIDPredator$SecRuD >= 0, FIDPredator$SecRuD, NA)
FIDPredator$notimportant <- ifelse( is.na(FIDPredator$SecRuD), FIDPredator$SecFID, NA)

FIDPredator$SecFlee <- rowSums(FIDPredator[,c("noimportant2", "notimportant")], na.rm=TRUE)*NA^!rowSums(!is.na(FIDPredator[,c("noimportant2", "notimportant")]))

FIDPredator = subset(FIDPredator, select = -c(noimportant2,notimportant ))

#Fly2:

FIDPredator$RunFly2 <- rep("Fly",nrow(FIDPredator))
FIDPredator$RunFly2[FIDPredator$SecRuD >= 0] <- "Run"
FIDPredator$RunFly2[is.na(FIDPredator$SecFlee)] <- NA

#  Calculating the second distance flee (SecEsD):

FIDPredator$SecFID
FIDPredator$SecEsDMe
FIDPredator$SecFIDcar
FIDPredator$ThirDirection

FIDPredator$ThirDirection <- as.numeric(as.character(FIDPredator$ThirDirection))
esddirection <- abs(FIDPredator$PredatorDirection - FIDPredator$ThirDirection)

# SecEsD = sqr(SecEsDMe^2 + SecFIDcar^2 - 2SecEsDMe*SecFIDcar·cos(esddirection)) 

root5 <- FIDPredator$SecEsDMe^2
root6 <- FIDPredator$SecFIDcar^2

cos6 <- 2*FIDPredator$SecEsDMe*FIDPredator$SecFIDcar*cos(esddirection*pi/180)

FIDPredator$SecEsD <- sqrt(root5 + root6 - (cos6))


#### 1.5) Method column####

# I will add a new column with the name of the approacher

FIDPredator$FIDmethod <- rep("Jackale",nrow(FIDPredator))

# I need to make a commend that if Notes == "Control", FIDmethod == "Control"

FIDPredator$FIDmethod[FIDPredator$Notes=='Control'] <- 'Control'

#Some graph to understand more my data:

table(FIDPredator$Habitat)

plot(FIDPredator$FID~FIDPredator$SecFID)

plot(FIDPredator$FirstFlee~FIDPredator$SecFlee)

ggplot(FIDPredator, aes(x = SecFID, y = FID, group = Habitat, colour = Habitat)) +
  geom_point() + geom_line() + theme_bw()

#Remove "Test" column from dataframe

FIDPredator <- FIDPredator[FIDPredator$Notes != "Test",]


#### 1.6) Map location as repeated measures ####

#Create MapLoc for FIDpredator:

#libraries:

library(leaflet)
library(sp)
library(rgdal)
library(lubridate)

#functions:

#this function round numbers so the digit i would like to use (not mine):
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  (df)
}

# I will start by making a grid of the map with ovelaping FID test - around 100 meter radius

coor2d <- data.frame(FIDPredator)
coor2d <- FIDPredator[complete.cases(coor2d$long),]

# make the data as real coordinats on the map

itm<-"+init=epsg:23031 +proj=tmerc +lat_0=31.73439361111111 +lon_0=35.20451694444445 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-48,55,52,0,0,0,0 +units=m +no_defs"
coordinates(coor2d) <- c("lat", "long")
proj4string(coor2d) <- CRS(itm)
res2 <- spTransform(coor2d, CRS(itm))

#cheking how the coordinate look:
plot(res2)

#next, I will transform my coordinate into the right UTM:

cord.UTM2 <- spTransform(res2, CRS("+init=epsg:23031"))
plot(cord.UTM2, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)

#now that i have the coordinates as UTM, i can transform it again into data frame so i can work with the data:

coord.UTM.2.2 <- data.frame(cord.UTM2)

# now i will round the coordinates so i can make a gride of every 100 meter radius, and consider the FIDs that are in this radius
# as repeated measures:

rounddf1.1 <- round_df(coord.UTM.2.2,1) #round 1000 meters
rounddf.11 <- round_df(coord.UTM.2.2,2) #round 100 meters
rounddf2.1 <- round_df(coord.UTM.2.2,3) #round 10 meters

# I will see how the plot look to choose the best plot

plot(rounddf1.1$lat~rounddf1.1$long, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95,col = "black", cex=1)
plot(rounddf.11$lat~rounddf.11$long, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95,col = "black", cex=1)
plot(rounddf2.1$lat~rounddf2.1$long, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95,col = "black", cex=1)

# I will make a new column that connect the two new coordinates to the first data frame (coor) and rename it:

rounddf.11$MapLoc <- paste(rounddf.11$lat,rounddf.11$long,sep = "_")

try <- subset(rounddf.11,select =  c("MapLoc","lat","long"))
table(try$MapLoc,try$long)

#After i saw it match the rounding numbers, I will add it to the main dataframe:

rounddf.11$MapLoc <- paste(rounddf.11$lat,rounddf.11$long,sep = "_")
roundfortest <- rounddf.11[c("MapLoc","num")]

Test <- merge(roundfortest,FIDPredator,by="num")



#### 1.7) Working on human FID####


#Adding the human FID data:


setwd("C:/Users/OrrS4/Dropbox (1)/Orr lab - miki's shared workspace/PhD/FID")

FID <-  read.csv("Coordination for statistics.csv")

colnames (FIDPredator)
colnames (FID)

# Changing the names to the column names of FIDpredator:

names(FID)[names(FID) == "X.2"] <- "num"    
names(FID)[names(FID) == "Y"] <- "lat"    
names(FID)[names(FID) == "X"] <- "long"    
names(FID)[names(FID) == "Number_of_lapwings"] <- "GroupSize"    
names(FID)[names(FID) == "Starting_distance"] <- "StD"    
names(FID)[names(FID) == "Running_distance"] <- "RuD"    
names(FID)[names(FID) == "Attack.escape"] <- "AttackEscape"    
names(FID)[names(FID) == "Where_escaped"] <- "HideOrNot"    
names(FID)[names(FID) == "Escape_distance"] <- "EsD"    
names(FID)[names(FID) == "Starting_distance2"] <- "SecStD"    
names(FID)[names(FID) == "Running_distance_2"] <- "SecRuD"    
names(FID)[names(FID) == "FID_2"] <- "SecFID"    
names(FID)[names(FID) == "Attack.escape2"] <- "SecAttackEscape"    
names(FID)[names(FID) == "Where_escaped2"] <- "SecHideOrNot"    
names(FID)[names(FID) == "rounddf.X_Y"] <- "MapLoc"    
names(FID)[names(FID) == "Escape_distance.1"] <- "SecEsD"    

FID$FIDmethod <- rep("Human",nrow(FID))

#make column for FID also with fly and escape:

# I need to create a new column "FirstFlee" - that take into considuration the FID, but if there is
# RuD it takes the RuD

FID$noimportant2 <- ifelse( FID$RuD > 0, FID$RuD, NA)
FID$notimportant <- ifelse( is.na(FID$noimportant2), FID$FID, NA)

FID$FirstFlee <- rowSums(FID[,c("noimportant2", "notimportant")], na.rm=TRUE)*NA^!rowSums(!is.na(FID[,c("noimportant2", "notimportant")]))

#and removing the column i dont need:

FID = subset(FID, select = -c(noimportant2,notimportant) )

# Now I need a value of if run by foot or fly:

FID$RunFly1 <- rep("Fly",nrow(FID))
FID$RunFly1[FID$RuD > 0] <- "Run"
FID$RunFly1[is.na(FID$FirstFlee)] <- NA

# Make the second FID with RuD as SecFlee:

FID$noimportant2 <- ifelse( FID$SecRuD > 0, FID$SecRuD, NA)
FID$notimportant <- ifelse( FID$SecRuD <=0, FID$SecFID, NA)

FID$SecFlee <- rowSums(FID[,c("noimportant2", "notimportant")], na.rm=TRUE) *NA^!rowSums(!is.na(FID[,c("noimportant2", "notimportant")]))

FID = subset(FID, select = -c(noimportant2,notimportant ))

#Fly2:

FID$RunFly2 <- rep("Fly",nrow(FID))
FID$RunFly2[FID$SecRuD > 0] <- "Run"

FID$RunFly2[is.na(FID$SecFlee)] <- NA


##########  Making a grid of the map with ovelaping FID test - around 100 meter radius

coor2d <- data.frame(FID)

# make the data as real coordinats on the map

itm<-"+init=epsg:23031 +proj=tmerc +lat_0=31.73439361111111 +lon_0=35.20451694444445 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-48,55,52,0,0,0,0 +units=m +no_defs"
coordinates(coor2d) <- c("lat", "long")
proj4string(coor2d) <- CRS(itm)
res2 <- spTransform(coor2d, CRS(itm))

#cheking how the coordinate look:
plot(res2)

#next, I will transform my coordinate into the right UTM:

cord.UTM2 <- spTransform(res2, CRS("+init=epsg:23031"))
plot(cord.UTM2, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)

#now that i have the coordinates as UTM, i can transform it again into data frame so i can work with the data:

coord.UTM.2.2 <- data.frame(cord.UTM2)

# now i will round the coordinates so i can make a gride of every 100 meter radius, and consider the FIDs that are in this radius
# as repeated measures:

rounddf1.1 <- round_df(coord.UTM.2.2,1) #round 1000 meters
rounddf.1 <- round_df(coord.UTM.2.2,2) #round 100 meters
rounddf2.1 <- round_df(coord.UTM.2.2,3) #round 10 meters

# I will see how the plot look to choose the best plot

plot(rounddf1.1$lat~rounddf1.1$long, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95,col = "black", cex=1)
plot(rounddf.1$lat~rounddf.1$long, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95,col = "black", cex=1)
plot(rounddf2.1$lat~rounddf2.1$long, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95,col = "black", cex=1)

# I will make a new column that connect the two new coordinates to the first data frame (coor) and rename it:

rounddf.1$MapLoc <- paste(rounddf.1$lat,rounddf.1$long,sep = "_")

try <- subset(rounddf.1,select =  c("MapLoc","lat","long"))
table(try$MapLoc,try$long)

#After i saw it match the rounding numbers, I will add it to the main dataframe:

FID$MapLoc2 <- paste(rounddf.1$lat,rounddf.1$long,sep = "_")

try <- as.data.frame(table(FID$MapLoc2,FID$MapLoc))

#The old MapLoc was not arranged in the right order! the new MapLoc is the new locations:

FID$MapLoc <-FID$MapLoc2




#### 1.8) Combine the two data into one dataframe and arranging it####

colnames (FID)
colnames (FIDPredator)

Date2 <- subset(FID, select =c("Date","Name","Time", "lat", "long","Habitat","GroupSize","StD", 
                               "RuD", "FID","AttackEscape","HideOrNot","FIDmethod", "EsD","SecStD",
                               "SecRuD","SecFID","SecAttackEscape","SecHideOrNot", "SecEsD","Nesting",
                               "MapLoc","FirstFlee", "RunFly1","SecFlee","RunFly2","Notes"))
Date3 <- subset(Test, select =c("Date","Name","Time", "lat", "long","Habitat","GroupSize","StD", 
                                "RuD", "FID","AttackEscape","HideOrNot","FIDmethod", "EsD","SecStD",
                                "SecRuD","SecFID","SecAttackEscape","SecHideOrNot", "SecEsD","Nesting",
                                "MapLoc","FirstFlee", "RunFly1","SecFlee","RunFly2","Notes"))

FIDbin <- rbind(Date2,Date3)
summary(FIDbin)

#Changing the names of habitats so it will match values of predator:

FIDbin$Habitat[FIDbin$Habitat == "Agriculture"] = "Open field"
FIDbin$Habitat[FIDbin$Habitat == "Cow house"] = "Urban"

#Changing the hideornot of the human:

FIDbin$HideOrNot[FIDbin$HideOrNot == "Across a lake"] = "Hide"
FIDbin$HideOrNot[FIDbin$HideOrNot == "Bushes"] = "Hide"
FIDbin$HideOrNot[FIDbin$HideOrNot == "Other (in notes)"] = "Hide"
FIDbin$HideOrNot[FIDbin$HideOrNot == "Few meters"] = "Didn't hide"
FIDbin$HideOrNot[FIDbin$HideOrNot==""]  <- NA
table(FIDbin$HideOrNot)


FIDbin$SecHideOrNot[FIDbin$SecHideOrNot == "Across a lake"] = "Hide"
FIDbin$SecHideOrNot[FIDbin$SecHideOrNot == "Bushes"] = "Hide"
FIDbin$SecHideOrNot[FIDbin$SecHideOrNot == "Other (in notes)"] = "Hide"
FIDbin$SecHideOrNot[FIDbin$SecHideOrNot == "Few meters"] = "Didn't hide"
FIDbin$SecHideOrNot[FIDbin$SecHideOrNot==""]  <- NA
#FIDbin$SecHideOrNot[FIDbin$SecHideOrNot=-c("","Few meters","Other (in notes)","Bushes","Across a lake")]

table(FIDbin$SecHideOrNot)

table(FIDbin$SecRuD)
table(FIDbin$RuD)

FIDbin$RuD[FIDbin$RuD <= 0 ]  <- NA
FIDbin$SecRuD[FIDbin$SecRuD <= 0 ]  <- NA

table(FIDbin$FIDmethod)

# Rounding all the numbers:
FIDbin <- FIDbin %>% mutate_at(vars(SecFlee, SecEsD, SecFID, SecRuD, SecStD, EsD), funs(round(., 0)))


#The data is ready to be used for calculations!

str(FIDbin)

#Date - date format V
#Name - factor V
#Time - POXIS V
#lat - number - degrees V
#lon - number - degree V
#Habitat - factor V
#GroupSize - number V
#StD - number V
#RuD - Number V
#FID - Number V
#AttackEscape - Factor - see if its significant - if not, remove V
#HideNot - Factor V
#FIDmethod - change to factors V
#EsD - Num V
#SecStD - Num V
# SecRuD - Num V
# SecFID - Num V
# SecAttackEscape - Factor - see if its significant - if not, remove V V
# SecHideOrNot - Factor V
# SecEsD - Num V
# Nesting - See how much nest i have - not significant - renmove them V
#MapLoc - Factor V
# FirstFlee - num V 
# RunFly1 - Factor - change V
# SecFlee - Num V
# RunFly2 - Factor V
# Notes - not important

#Need to add:
# MWD - Walking distance - StD - FirstFlee

# Now to fix the data:

# Date:
FIDbin$Date <- as.POSIXct(as.character(FIDbin$Date), format = "%d/%m/%Y" , tz="UTC")

# Time:
FIDbin$Time <- as.POSIXct(as.character(FIDbin$Time),  format = "%H:%M:%S" , tz="UTC")

str(FIDbin)

FIDbin$month <- month(FIDbin$Date)
plot(table(FIDbin$Date, FIDbin$FIDmethod))

# Plot that show in when jackal and human were tested:
ggplot(FIDbin, aes(x = factor(month), y = FIDmethod)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("month")

# Change all the int into numbers:

FIDbin[8:11] <- lapply(FIDbin[8:11], as.numeric)

names(FIDbin)
# General data:

table(FIDbin$Habitat)

# How much data I have from each test:
table(FIDbin$FIDmethod)
control <- subset(FIDbin, FIDmethod == "Control" )
table(control$Habitat)
Jackal <- subset(FIDbin, FIDmethod == "Jackale" )
table(Jackal$Habitat)
Human <- subset(FIDbin, FIDmethod == "Human" )
table(Human$Habitat)

# Month tested:

table(FIDbin$month, FIDbin$FIDmethod)
table(FIDbin$FIDmethod, FIDbin$Nesting)
table(FIDbin$FIDmethod, FIDbin$HideOrNot) # Good data
table(FIDbin$month, FIDbin$HideOrNot)
table(FIDbin$RunFly1, FIDbin$HideOrNot)
table(FIDbin$Nesting, FIDbin$HideOrNot)



#### 1.9) Mapping the locations####

### I will start by making a plot with the general areas where I tested the birds:

FIDmaploc <- FIDbin[c("FIDmethod","MapLoc", "lat", "long")]

plot(FIDmaploc$lat~FIDmaploc$long)
#this function from earlier if I dont have it again:
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  (df)
}

# I will start by making a grid of the map with ovelaping FID test - around 100 meter radius

coor2d <- FIDmaploc
coor2d <- FIDmaploc[complete.cases(coor2d$long),]

# make the data as real coordinats on the map

itm<-"+init=epsg:23031 +proj=tmerc +lat_0=31.73439361111111 +lon_0=35.20451694444445 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-48,55,52,0,0,0,0 +units=m +no_defs"

coordinates(coor2d) <- c("long", "lat")
proj4string(coor2d) <- CRS(itm)
res2 <- spTransform(coor2d, CRS("+init=epsg:23031"))

#cheking how the coordinate look:
plot(res2, axes = TRUE)
ll <- leaflet(data=coor2d) %>% addTiles %>% addCircleMarkers(weight = 6 ,fillOpacity = 0.9, radius = 4, color = "Red") %>%
  addProviderTiles("Stamen.Terrain", group="Positron") %>%
  addProviderTiles("Stamen.Toner", group = "Toner") %>%
  addProviderTiles('Esri.WorldImagery', group = "Satellite") %>%
  addProviderTiles('Esri.WorldShadedRelief', group = "Relief") %>%
  addLayersControl(coor2d,baseGroups = c("Positron", "Toner", "Satellite", "Relief"),
                   options = layersControlOptions(collapsed = FALSE, weight = 1,fillOpacity = 0.5))%>%
  addScaleBar(position = c("bottomleft"), options = scaleBarOptions(maxWidth = 100, metric = TRUE, imperial = FALSE,
                                                                    updateWhenIdle = TRUE))

ll #open map

#next, I will transform my coordinate into the right UTM:

res3 <- data.frame(res2)

rounddf.11 <- round_df(res3,2) #round 100 meters

rounddf.11$MapLoc2 <- paste(rounddf.11$lat,rounddf.11$long,sep = "_")

# have to each coordinate a ID number:

listmap <- unique(rounddf.11$MapLoc2)

for (i in 1:length(listmap)) {
  rounddf.11$MapNumber[rounddf.11$MapLoc2 == listmap[i]] <- i
}

coordinates(rounddf.11) <- c("long", "lat") #coordinate of the grid
proj4string(rounddf.11) <- CRS("+init=epsg:23031")
llpd <- spTransform(rounddf.11, CRS(itm))
plot(rounddf.11, axes = TRUE)
plot(llpd, axes = TRUE)

# Icon of jackal and human for the map:
IconsFID <- iconList(
  Jackale = makeIcon("wild-animals.png", iconWidth = 25, iconHeight = 25),
  Human = makeIcon("walk.png", iconWidth = 25, iconHeight = 25))

ll <- leaflet(data=llpd) %>% addTiles%>% addProviderTiles("OpenStreetMap.DE")  %>% addMarkers(data=llpd, icon = ~IconsFID) %>% 
  addScaleBar(position = c("bottomleft"), options = scaleBarOptions(maxWidth = 100, metric = TRUE, imperial = FALSE,
                                                                    updateWhenIdle = TRUE))  %>%
  addLabelOnlyMarkers(
    label=~MapNumber, 
    labelOptions = labelOptions(noHide = TRUE, textOnly = TRUE, direction = "bottom",style = list(
      "color" = "Blue")))
ll #open map

FIDbin$Maploc2 <- rounddf.11$MapLoc2
FIDbin$MapNumber <-rounddf.11$MapNumber

plot(table(FIDbin$MapNumber), bty="n" , main="FID tested locations",
     xlab="Locations",
     ylab="Freq")


#### 1.10) Removing unuable data####

# Next to look into data I might remove:

# I dont want nesting data - lets see how it looks:

table(FIDbin$Nesting,FIDbin$FIDmethod)

# Too little yes - remove the lines where yes and maybe are.

FIDbin2 <- FIDbin[ !(FIDbin$Nesting %in% c("Yes","Maybe")), ]
table(FIDbin2$Nesting,FIDbin2$FIDmethod)


FIDbin2 <- subset(FIDbin2, select = -c(Nesting))

# Now lets see how the data look now

names(FIDbin2)

table(FIDbin2$AttackEscape,FIDbin2$FIDmethod) 
table(FIDbin2$Habitat,FIDbin2$FIDmethod)  # Change the "Other" to a normal habitat

Other <- subset(FIDbin2, Habitat== "Other (in notes)" ) # 32.52954 35.43082 - water pond

FIDbin2$Habitat[FIDbin2$Habitat == "Other (in notes)"] <- "Water pond"

FIDbin3 <- FIDbin2[ !(FIDbin2$FIDmethod %in% c("Control")), ] # Remove control for now

table(FIDbin3$Habitat,FIDbin3$FIDmethod)


ggplot(data = FIDbin3, aes(x = FIDmethod, y = FirstFlee, color = Habitat)) + 
  geom_boxplot()

names(FIDbin3)

# Remove the other column that are not intresting:

FIDbin4 <- subset(FIDbin3, select = -c(Name,lat,long,AttackEscape,SecAttackEscape, Notes,RuD,FID,SecRuD,SecFID)) # Running distance and FID are BOTH FID! no need for them to be separated

names(FIDbin4)

ggplot(data = FIDbin4, aes(x = StD, y = FirstFlee, color = Habitat)) + 
  geom_point() + geom_smooth(method=lm, fullrange=TRUE)+ facet_grid(~FIDmethod)

ggplot(FIDbin4, aes(y = Habitat, x = FirstFlee,color=FIDmethod)) + 
  geom_jitter() +
  theme_minimal()
str(FIDbin4)

table(FIDbin4$MapLoc,FIDbin4$Habitat )

names(FIDbin4)


#### 1.11)Creating the WD####

# Buffer walking:

FIDbin4$WD <- FIDbin4$StD - FIDbin4$FirstFlee

plot(FIDbin4$StD~FIDbin4$WD)


#### 1.12) Save the data####

write.csv(FIDbin4, file = "FID after cleaning.csv")

######## 2) Preparing the FID data########
#### 2.1) Open file and library ####
library(tidyverse)
library(modelr) 
library(rgl)
library(shiny)
library(knitr)
library(qpcR)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggplot2)
library(lattice)
library(lme4)
library(gplots)
library(multcomp)
library(dplyr)
library(AICcmodavg)#for aictab
## library(lmerTest)#without this package i cannot tell the p-value in lme4
library(lubridate) # to work with time
library(rptR)
library(qpcR)
library(sp)
library(AICcmodavg)
library(ggpubr)

# some more variables:

FIDbin4$Minuteday <- minute(FIDbin4$Time) +60*hour(FIDbin4$Time) # Time in minutes

# Groups as factors
FIDbin4$Groups <- NA
FIDbin4$Groups[FIDbin4$GroupSize == 1] <- "Single"
FIDbin4$Groups[FIDbin4$GroupSize == 2] <- "Pair"
FIDbin4$Groups[FIDbin4$GroupSize > 2] <- "group"

FIDbin4$secWD <- FIDbin4$SecStD -FIDbin4$SecFlee 

#setwd("C:/Users/OrrS4/Dropbox (1)/Orr lab - miki's shared workspace/PhD/FID")
#FIDbin <- read.csv("FID after cleaning.csv")

str(FIDbin)
str(FIDbin4)
names(FIDbin4)
FIDbin4$Habitat <- factor(FIDbin4$Habitat)
FIDbin4$HideOrNot <- factor(FIDbin4$HideOrNot)
FIDbin4$SecHideOrNot <- factor(FIDbin4$SecHideOrNot)
FIDbin4$SecHideOrNot <- factor(FIDbin4$SecHideOrNot)

##################V 2.2) Is the FID affected by habitat? # FirstFlee, Habitat ####

plot(FIDbin4$FirstFlee~FIDbin4$Habitat)
plot(FIDbin4$FirstFlee,FIDbin4$FIDmethod)

table(FIDbin4$Habitat, FIDbin4$FIDmethod)

hist(FIDbin4$FirstFlee)
ggboxplot(FIDbin4, x = "Habitat", y = "FirstFlee", color = "FIDmethod")
ggboxplot(FIDbin4, x = "FIDmethod", y = "FirstFlee", color = "Habitat")

ggline(FIDbin4, x = "Habitat", y = "FirstFlee", color = "FIDmethod",
       add = c("mean_se", "dotplot"))

ggline(FIDbin4, x = "FIDmethod", y = "FirstFlee", color = "Habitat",
       add = c("mean_se", "dotplot"))

res.aov2 <- aov(FirstFlee ~ FIDmethod + Habitat, data = FIDbin4)
res.aov3 <- aov(FirstFlee ~ FIDmethod + Habitat + FIDmethod:Habitat, data = FIDbin4)

summary(res.aov2)
summary(res.aov3) # Not interact # Dont use this model

TukeyHSD(res.aov2, which = "Habitat")
#                          diff        lwr       upr     p adj
#Urban-Open field      -10.690984 -15.125972 -6.255995 0.0000001
#Water pond-Open field   2.326362  -1.697962  6.350686 0.3625177
#Water pond-Urban       13.017345   9.207344 16.827347 0.0000000
# Urban is different between the two
TukeyHSD(res.aov2, which = "FIDmethod")

plot(res.aov2, 1) # 96, 92 and 37 are outliners
plot(res.aov2, 2)

# Extract the residuals
aov_residuals <- residuals(object = res.aov2)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals ) # Data is normal

ggplot(data=FIDbin4, aes (x =  FIDmethod, y=FirstFlee , fill = factor(Habitat, level = c('Urban', 'Open field', 'Water pond')))) + 
  geom_boxplot(notch = TRUE) +scale_fill_manual(values=c("Grey", "Green", "Blue")) +
  theme_minimal() + ggtitle("FID between habitat") +
  xlab("Habitat") + ylab("FID")  + labs(fill= "Habitat")

ggplot(data=FIDbin4, aes (x =  reorder(Habitat,FirstFlee), y=FirstFlee)) + 
  geom_boxplot(fill='#A4A4A4', color="black") +
  theme_minimal() + ggtitle("FID between habitat") +
  xlab("Habitat") + ylab("FID")  + labs(fill= "Habitat") +
  theme_classic()

res.aovtry <- aov(FirstFlee ~ Habitat, data = FIDbin4)

summary(res.aovtry)

# The method as well as the habitat affect FID - with urban having a strong significant over the others

##################V 2.3) Is the EsD affected by the habitat? Method? EsD, Habitat, FIDmethod ####

plot(FIDbin4$EsD~FIDbin4$Habitat)
plot(FIDbin4$EsD~FIDbin4$FIDmethod)

hist(FIDbin4$EsD) # Not normal distributade

FIDbin4$LogEsD <- log(FIDbin4$EsD)
hist(FIDbin4$LogEsD)
shapiro.test(FIDbin4$LogEsD) # Data is normal

ggboxplot(FIDbin4, x = "Habitat", y = "LogEsD", color = "FIDmethod")
ggboxplot(FIDbin4, x = "FIDmethod", y = "LogEsD", color = "Habitat")

ggline(FIDbin4, x = "Habitat", y = "LogEsD", color = "FIDmethod",
       add = c("mean_se", "dotplot"))

ggline(FIDbin4, x = "FIDmethod", y = "LogEsD", color = "Habitat",
       add = c("mean_se", "dotplot"))

res.aov4 <- aov(LogEsD ~ FIDmethod + Habitat, data = FIDbin4)
res.aov5 <- aov(LogEsD ~ FIDmethod + Habitat + FIDmethod:Habitat, data = FIDbin4)

summary(res.aov4)
summary(res.aov5) # Not interact - almost significant # Dont use this model

TukeyHSD(res.aov4, which = "Habitat")
#                          diff        lwr       upr     p adj
#Urban-Open field      -0.4591689 -0.855452167 -0.0628856 0.0185841
#Water pond-Open field -0.1164277 -0.485857901  0.2530025 0.7369160
#Water pond-Urban       0.3427412  0.006453105  0.6790293 0.0446658
# Small differences between urban and such

TukeyHSD(res.aov4, which = "FIDmethod") # Not significant

plot(res.aov4, 1) # 26, 32 and 241 are outliners
plot(res.aov4, 2)

# Extract the residuals
aov_residuals2 <- residuals(object = res.aov4)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals2 ) # Data is normal

ggplot(data=FIDbin4, aes (x =  factor(Habitat, level = c('Urban', 'Open field', 'Water pond')), y=LogEsD , fill = factor(Habitat, level = c('Urban', 'Open field', 'Water pond')))) + 
  geom_boxplot(notch = TRUE) +scale_fill_manual(values=c("Grey", "Green", "Blue")) +
  theme_minimal() + ggtitle("LogDF between habitat") +
  xlab("Habitat") + ylab("DF (log)")  + labs(fill= "Habitat")
(0.0185841 + 0.0446658)/2

modeDF <- glm(log(EsD) ~ FIDmethod + Habitat , data=FIDbin4, family = "gaussian")
summary(modeDF)

# sEsD

FIDbinDF <- subset(FIDbin4, SecEsD >= 0 )

hist(FIDbinDF$SecEsD)

FIDbinDF$SeclogEsD <- log(FIDbinDF$SecEsD)
hist(FIDbinDF$SeclogEsD)
shapiro.test(FIDbinDF$SeclogEsD) # Data is normalish

ggboxplot(FIDbinDF, x = "Habitat", y = "SeclogEsD", color = "FIDmethod")
ggboxplot(FIDbinDF, x = "FIDmethod", y = "SeclogEsD", color = "Habitat")

ggline(FIDbinDF, x = "Habitat", y = "SeclogEsD", color = "FIDmethod",
       add = c("mean_se", "dotplot"))

ggline(FIDbinDF, x = "FIDmethod", y = "SeclogEsD", color = "Habitat",
       add = c("mean_se", "dotplot"))

res.aov25 <- aov(SeclogEsD ~ FIDmethod + Habitat, data = FIDbinDF)
res.aov26 <- aov(SeclogEsD ~ FIDmethod + Habitat + FIDmethod:Habitat, data = FIDbinDF)

summary(res.aov25)
summary(res.aov26) # Not interact - almost significant # Dont use this model

TukeyHSD(res.aov25, which = "Habitat")
#                          diff        lwr       upr     p adj
#Urban-Open field      -0.4591689 -0.855452167 -0.0628856 0.0185841
#Water pond-Open field -0.1164277 -0.485857901  0.2530025 0.7369160
#Water pond-Urban       0.3427412  0.006453105  0.6790293 0.0446658
# Small differences between urban and such

TukeyHSD(res.aov25, which = "FIDmethod") # Not significant

plot(res.aov25, 1) # 26, 32 and 241 are outliners
plot(res.aov25, 2)

# Extract the residuals
aov_residuals2 <- residuals(object = res.aov25)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals2 ) # Data is normal

ggplot(data=FIDbinDF, aes (x =  factor(Habitat, level = c('Urban', 'Open field', 'Water pond')), y=SeclogEsD , fill = factor(Habitat, level = c('Urban', 'Open field', 'Water pond')))) + 
  geom_boxplot() +scale_fill_manual(values=c("Grey", "Green", "Blue")) +
  theme_minimal() + ggtitle("sLogDF between habitat") +
  xlab("Habitat") + ylab("sDF (log)")  + labs(fill= "Habitat")

modesDF <- glm(log(SecEsD) ~ FIDmethod + Habitat , data=FIDbin4, family = "gaussian")

summary(modesDF)

# There is a small differences between habitat, the method does not affect it. Urban have a small significant over the others

##################V 2.4) Is the StD is affected by habitat or method? StD, FIDmethod, Habitat ####

plot(FIDbin4$StD~FIDbin4$Habitat)
plot(FIDbin4$StD~FIDbin4$FIDmethod)

hist(FIDbin4$StD) # Not normal distributade
shapiro.test(FIDbin4$StD) # Normal

ggboxplot(FIDbin4, x = "Habitat", y = "StD", color = "FIDmethod")
ggboxplot(FIDbin4, x = "FIDmethod", y = "StD", color = "Habitat")

ggline(FIDbin4, x = "Habitat", y = "StD", color = "FIDmethod",
       add = c("mean_se", "dotplot"))

ggline(FIDbin4, x = "FIDmethod", y = "StD", color = "Habitat",
       add = c("mean_se", "dotplot"))

res.aov6 <- aov(StD ~ FIDmethod + Habitat, data = FIDbin4)
res.aov7 <- aov(StD ~ FIDmethod + Habitat + FIDmethod:Habitat, data = FIDbin4)

summary(res.aov6)
summary(res.aov7) # Not interact - almost significant # Dont use this model

TukeyHSD(res.aov4, which = "Habitat") # No differences between water pond and 
# open field - small differences between urban and the other habitats

#                             diff        lwr       upr     p adj
#Urban-Open field      -11.5268409 -19.937021 -3.116661 0.0039387
#Water pond-Open field  -0.7600851  -8.424516  6.904346 0.9703797
#Water pond-Urban       10.7667558   3.528634 18.004877 0.0015271

TukeyHSD(res.aov6, which = "FIDmethod") # No differences between the methods and the habitats

plot(res.aov6, 1) # 263, 92 and 53 are outliners
plot(res.aov6, 2)

model.tables(res.aov6, type="means", se = TRUE)

# Extract the residuals
aov_residuals3 <- residuals(object = res.aov6)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals3) # Data is normal

# There is a difference abetween the stating distance - the urban starting distance is different...

plot(Jackal$StD~Jackal$Habitat)
plot(Human$StD~Human$Habitat)

# Urban StD is always smaller - as a result of the habitat structure

urban <- subset(FIDbin4, Habitat == "Urban")
Noturban <- subset(FIDbin4, Habitat != "Urban")

str(Noturban)
table(Noturban$Habitat)

summary(Noturban$StD)

table(Noturban$StD) # There is a starting point of 155 meter - maybe remove
par(mfrow= c(1,2))

hist(Noturban$StD)
table(urban$StD)
hist(urban$StD)
par(mfrow= c(1,1))

ggplot(data=FIDbin4, aes (x =  FIDmethod, y=StD , fill = factor(Habitat, level = c('Urban', 'Open field', 'Water pond')))) + 
  geom_boxplot(notch = TRUE) +scale_fill_manual(values=c("Grey", "Green", "Blue")) +
  theme_minimal() + ggtitle("SD between habitat and method of approach") +
  xlab("Approach methods") + ylab("SD")  + labs(fill= "Habitat")

modesd <- glm(StD ~ FIDmethod + Habitat , data=FIDbin4, family = "gaussian")

summary(modesd)

modeWD <- glm(WD ~ FIDmethod + Habitat , data=FIDbin4, family = "gaussian")
summary(modeWD)

modesd$rank

#### sSD

FIDbinsSD <- subset(FIDbin4,SecStD >= 0 )
par(mfrow = c(1, 1))

plot(FIDbinsSD$SecStD~FIDbinsSD$Habitat)
plot(FIDbinsSD$secWD~FIDbinsSD$Habitat)
#plot(FIDbinsSD$SecStD~FIDbinsSD$FIDmethod)

hist(FIDbinsSD$SecStD) 
shapiro.test(FIDbinsSD$SecStD) # Normal

ggboxplot(FIDbinsSD, x = "Habitat", y = "SecStD", color = "FIDmethod")
ggboxplot(FIDbinsSD, x = "FIDmethod", y = "SecStD", color = "Habitat")

ggline(FIDbinsSD, x = "Habitat", y = "SecStD", color = "FIDmethod",
       add = c("mean_se", "dotplot"))

ggline(FIDbinsSD, x = "FIDmethod", y = "SecStD", color = "Habitat",
       add = c("mean_se", "dotplot"))

res.aov21 <- aov(SecStD ~ FIDmethod + Habitat, data = FIDbinsSD)
res.aov22 <- aov(SecStD ~ FIDmethod + Habitat + FIDmethod:Habitat, data = FIDbinsSD)

summary(res.aov21)
summary(res.aov22) # Not interact - almost significant # Dont use this model

TukeyHSD(res.aov21, which = "Habitat") # No differences between water pond and 
# open field - small differences between urban and the other habitats

#                             diff        lwr       upr     p adj
#Urban-Open field      -11.5268409 -19.937021 -3.116661 0.0039387
#Water pond-Open field  -0.7600851  -8.424516  6.904346 0.9703797
#Water pond-Urban       10.7667558   3.528634 18.004877 0.0015271

TukeyHSD(res.aov21, which = "FIDmethod") # No differences between the methods and the habitats

plot(res.aov21, 1) # 90, 70 and 104 are outliners
plot(res.aov21, 2)

model.tables(res.aov21, type="means", se = TRUE)

# Extract the residuals
res.aov21 <- residuals(object = res.aov6)
# Run Shapiro-Wilk test
shapiro.test(x = res.aov21) # Data is normal

ggplot(data=FIDbinsSD, aes (x =  FIDmethod, y=SecStD , fill = factor(Habitat, level = c('Urban', 'Open field', 'Water pond')))) + 
  geom_boxplot() +scale_fill_manual(values=c("Grey", "Green", "Blue")) +
  theme_minimal() + ggtitle("sSD between habitat and method of approach") +
  xlab("Approach methods") + ylab("sSD (AD)")  + labs(fill= "Habitat")

ggplot(data=FIDbinsSD, aes (x =  factor(Habitat, level = c('Urban', 'Open field', 'Water pond')), y=SecStD , fill = FIDmethod)) + 
  geom_boxplot() +scale_fill_manual(values=c("Grey", "Green", "Blue")) +
  theme_minimal() + ggtitle("SD between habitat and method of approach") +
  xlab("Approach methods") + ylab("SD")  + labs(fill= "Habitat")

glmssd <- glm(formula = SecStD ~ FIDmethod + Habitat, family = "gaussian", 
              data = FIDbin4)

summary(glmssd)

#### sWD

FIDbinsWD <- subset(FIDbin4,secWD >= 0 )
par(mfrow = c(1, 1))

plot(FIDbinsWD$SecStD~FIDbinsWD$Habitat)
plot(FIDbinsWD$secWD~FIDbinsWD$Habitat)
#plot(FIDbinsSD$SecStD~FIDbinsSD$FIDmethod)

hist(FIDbinsWD$secWD) # Will log it
FIDbinsWD$LogsWD <- log(FIDbinsWD$secWD)
hist(FIDbinsWD$LogsWD)

shapiro.test(FIDbinsWD$LogsWD) # Normal

ggboxplot(FIDbinsWD, x = "Habitat", y = "LogsWD", color = "FIDmethod")
ggboxplot(FIDbinsWD, x = "FIDmethod", y = "LogsWD", color = "Habitat")

ggline(FIDbinsWD, x = "Habitat", y = "LogsWD", color = "FIDmethod",
       add = c("mean_se", "dotplot"))

ggline(FIDbinsWD, x = "FIDmethod", y = "LogsWD", color = "Habitat",
       add = c("mean_se", "dotplot"))

res.aov23 <- aov(LogsWD ~ FIDmethod + Habitat, data = FIDbinsWD)
res.aov24 <- aov(LogsWD ~ FIDmethod + Habitat + FIDmethod:Habitat, data = FIDbinsWD)

summary(res.aov23)
summary(res.aov24) # Not interact

TukeyHSD(res.aov23, which = "Habitat") # No differences between water pond and 
# open field - small differences between urban and the other habitats

#                             diff        lwr       upr     p adj
#Urban-Open field      -11.5268409 -19.937021 -3.116661 0.0039387
#Water pond-Open field  -0.7600851  -8.424516  6.904346 0.9703797
#Water pond-Urban       10.7667558   3.528634 18.004877 0.0015271

TukeyHSD(res.aov23, which = "FIDmethod") # No differences between the methods and the habitats

plot(res.aov23, 1) # 232, 191 and 116 are outliners
plot(res.aov23, 2)

model.tables(res.aov23, type="means", se = TRUE)

# Extract the residuals
res.aov23 <- residuals(object = res.aov23)
# Run Shapiro-Wilk test
shapiro.test(x = res.aov21) # Data is normal

ggplot(data=FIDbinsWD, aes (x =  FIDmethod, y=LogsWD , fill = factor(Habitat, level = c('Urban', 'Open field', 'Water pond')))) + 
  geom_boxplot() +scale_fill_manual(values=c("Grey", "Green", "Blue")) +
  theme_minimal() + ggtitle("sWD between habitat and method of approach") +
  xlab("Approach methods") + ylab("sWD")  + labs(fill= "Habitat")

ggplot(data=FIDbinsWD, aes (x =  factor(Habitat, level = c('Urban', 'Open field', 'Water pond')), y=LogsWD , fill = FIDmethod)) + 
  geom_boxplot() +scale_fill_manual(values=c("Grey", "Green", "Blue")) +
  theme_minimal() + ggtitle("SD between habitat and method of approach") +
  xlab("Approach methods") + ylab("SD")  + labs(fill= "Habitat")

##################V 2.5) hidden or flying vs. habitats vs. method # Habitat, HideOrNot, FIDmethod # Return to it later ####

# Is the hidden or not affected by the habitat, and does it change by methods?                                                                                      ., scales = "free")

# Logistic Mixed Effects Model with Interaction Term did not work (the groups gave singularities)#

# Logistic regression:
# Each category is dependent on the previous one. Habitat -> Run/fly -> Hide/not.
# Depended variable - hide or not as binary. 

# For testing what affect hiding:
FIDbin4$HideZero[FIDbin4$HideOrNot == "Hide"] <- 0
FIDbin4$HideZero[FIDbin4$HideOrNot == "Didn't hide"] <- 1

# For testing what affect running:
FIDbin4$RunZero[FIDbin4$RunFly1 == "Run"] <- 0
FIDbin4$RunZero[FIDbin4$RunFly1 == "Fly"] <- 1

FIDbin7 <- subset(FIDbin4,HideZero >=0 )
FIDbin8 <- subset(FIDbin4,RunZero >=0 ) # To avoid NAs

table(FIDbin8$RunZero,FIDbin8$Habitat,FIDbin8$FIDmethod) #seems like more running in urban
table(FIDbin7$HideZero,FIDbin7$Habitat,FIDbin7$FIDmethod) #seems no diffrences

names(FIDbin7)
names(FIDbin8)
# creating models:
# With random effect
h0 <- glmer(
  HideZero ~ (1 | MapLoc), 
  data = FIDbin7, 
  family = binomial) # boundary fit - cant use random effect

# Only the running
h1 <- glm(HideZero ~ RunFly1, 
          data = FIDbin7, 
          family = binomial) # verry effective (of course)
summary(h1)
#All that can effect the hiding
h2 <- glm(HideZero ~ RunFly1*scale(EsD) + RunFly1  + scale(StD) + scale(EsD) + scale(FirstFlee) + scale(month) + scale(WD) + scale(Minuteday) + Groups + Habitat*FIDmethod + Habitat*RunFly1 ,  # Running and method might interact with habitat
          data = FIDbin7, 
          family = binomial)

anova(h2, test="Chisq") #Removing interaction with running, WD, FID, month

h3 <- glm(HideZero ~ RunFly1 + scale(StD) + scale(EsD) + scale(Minuteday) + Groups + Habitat*FIDmethod,  # Running and method might interact with habitat
          data = FIDbin7, 
          family = binomial)

anova(h3, test="Chisq") #Removing groups as well

h4 <- glm(HideZero ~ RunFly1 + scale(StD) + scale(EsD) + scale(Minuteday) + Habitat*FIDmethod,  # Running and method might interact with habitat
          data = FIDbin7, 
          family = binomial)

anova(h4, test="Chisq") # Remove Habitat and method

h5 <- glm(HideZero ~ RunFly1 + scale(StD) + scale(EsD) ,  # Running and method might interact with habitat
          data = FIDbin7, 
          family = binomial)


anoh5 <- anova(h5, test="Chisq") 
summary(h5) 

model = glm(HideZero ~ RunFly1 * scale(EsD), data = FIDbin7, family="binomial")
model3 = glm(HideZero ~ RunFly1 + EsD, data = FIDbin7, family="binomial")
model4 = glm(HideZero ~ EsD, data = FIDbin7, family="binomial")
model4 = glm(HideZero ~ RunFly1, data = FIDbin7, family="binomial")

summary(model4)
anova(h5, test="Chisq") 
anova(model2, test="Chisq") 
anova(model3, test="Chisq") 

# The Desicion to hide was not effected by habitat or method - but by the EsD (DF = 145, P < 0.05), if it run or flow 
# (DF = 147, p = 0.018) and in a small amount by the SD (DF = 146, p = 0.035)

FIDbin9 <- subset(FIDbin7, EsD >= 0 )

#model2
summary(model)
summary(model3)

# new.data for predicting new observations

# Creating binary model:

SEsE <- scale(FIDbin9$EsD)
summary(SEsE)
summary(FIDbin9$EsD)
summary(FIDbin9$RunFly1)

new.data <- expand.grid(EsD    = seq(5, 270, length = n), 
                        RunFly1 = c("Run", "Fly"))

preds <- predict(model3, newdata = new.data, type = 'response',se = TRUE)

new.data$pred.full <- preds$fit

new.data$ymin <- new.data$pred.full - 2*preds$se.fit
new.data$ymax <- new.data$pred.full + 2*preds$se.fit  

# Plotting
ggplot(FIDbin9,aes(x = EsD, y = HideZero, col = RunFly1, fill = RunFly1)) + 
  geom_point(size = 2.5) + 
  geom_ribbon(data = new.data,aes(y = pred.full, ymin = ymin, ymax = ymax),alpha = 0.7) +
  xlab("DF") +
  ylab('Hide/Did not hide') + labs(col = "Run or fly", fill = "Run or fly")

ggplot(FIDbin9,aes(x = EsD, y = HideZero, col = RunFly1)) + 
  geom_point(size = 2.5) + 
  geom_smooth(data = new.data,aes(y = pred.full, ymin = ymin, ymax = ymax),alpha = 0.7) +
  xlab("DF") +
  ylab('Hide/Did not hide') + labs(col = "Run or fly")

### Running effect:

names(FIDbin8)

r1 <- glm(
  RunZero ~ Habitat*FIDmethod + month + FIDmethod + scale(StD) + scale(FirstFlee) + scale(Minuteday) + Groups, 
  data = FIDbin8, 
  family = binomial) #boundary fit

summary(r1)
anova(r1, test="Chisq") #Removing everything that is not significant

r2 <- glm(
  RunZero ~ Habitat, 
  data = FIDbin8, 
  family = binomial) #boundary fit

anova(r2, test="Chisq") # The only thing that affect running is the habitat

r3 <- glm(
  RunZero ~ month, 
  data = FIDbin8, 
  family = binomial) #boundary fit
anova(r3, test="Chisq") # The only thing that affect running is the habitat

summary(r2)
cq <- chisq.test(FIDbin8$Habitat,FIDbin8$RunFly1)
cq2 <- chisq.test(FIDbin8$HideOrNot,FIDbin8$RunFly1)

# The decision to fly or not was more affected in urban habitats (Z = -2.113 +- 0.4497, p = 0.0346), but was 

RunNot <- as.data.frame(table(FIDbin8$Habitat,FIDbin8$RunFly1))
HideOrNot <-(as.data.frame(table(FIDbin8$HideOrNot,FIDbin8$RunFly1)))

# Plotting:
ggplot(HideOrNot, aes(x=  factor(Var1, level = c('Urban', 'Open field', 'Water pond')), y=Freq, fill = Var2))+
  geom_bar(stat='identity', position=position_dodge())+
  ylab("Freq") + xlab("Habitat") + labs( fill = "Run or fly")+
  geom_text(aes(label=Freq), size=3.5, vjust=-0.4,
            position = position_dodge (0.8))+
  theme_minimal()

# How much percentage did lapwings run:

21/68*100
8/67*100
5/140*100

# Lapwing decided to flee by running 30% in urban habitat, while in open field it was only 11% and in water ponds 3.5%

###sFly,sHide:

# Is the hidden or not affected by the habitat, and does it change by methods?                                                                                      ., scales = "free")

# Logistic Mixed Effects Model with Interaction Term did not work (the groups gave singularities)#

# Logistic regression:
# Each category is dependent on the previous one. Habitat -> Run/fly -> Hide/not.
# Depended variable - hide or not as binary. 

# For testing what affect hiding:
FIDbin4$HideZero2[FIDbin4$SecHideOrNot == "Hide"] <- 0
FIDbin4$HideZero2[FIDbin4$SecHideOrNot == "Didn't hide"] <- 1

# For testing what affect running:
FIDbin4$RunZero2[FIDbin4$RunFly2 == "Run"] <- 0
FIDbin4$RunZero2[FIDbin4$RunFly2 == "Fly"] <- 1

FIDbin10 <- subset(FIDbin4,HideZero2 >=0 )
FIDbin11 <- subset(FIDbin4,RunFly2 >=0 ) # To avoid NAs

table(FIDbin10$RunZero2,FIDbin10$Habitat,FIDbin10$FIDmethod) #seems like more running in urban
table(FIDbin11$HideZero2,FIDbin11$Habitat,FIDbin11$FIDmethod) #seems no diffrences

names(FIDbin10)
names(FIDbin11)
# creating models:
# With random effect
h01 <- glmer(
  HideZero2 ~ (1 | MapLoc), 
  data = FIDbin10, 
  family = binomial) # boundary fit - cant use random effect

# Only the running
h11 <- glm(HideZero2 ~ RunFly2, 
           data = FIDbin10, 
           family = binomial) # not significant
summary(h11)

anova(h22, test="Chisq") #Removing interaction with running, WD, FID, month

h33 <- glm(HideZero2 ~ RunFly2 + scale(SecStD) + scale(SecFlee) + scale(Minuteday) + Habitat + FIDmethod,  # Running and method might interact with habitat
           data = FIDbin10, 
           family = binomial) # Method and habitat
summary(h33)
anova(h33, test="Chisq") # Minutes of day and sFID are least significant

h44 <- glm(HideZero2 ~ RunFly2 + scale(SecStD) + Habitat + FIDmethod,  # Running and method might interact with habitat
           data = FIDbin10, 
           family = binomial) # Only significant is with habitat, removing interactions
summary(h44)
anova(h44, test="Chisq") # Minutes of day and sFID are least significant

h55 <- glm(HideZero2 ~ RunFly2 + Habitat + FIDmethod,  # Running and method might interact with habitat
           data = FIDbin10, 
           family = binomial) 
summary(h55)
anova(h55, test="Chisq")


h66 <- glm(HideZero2 ~ RunFly2 + Habitat ,  # Running and method might interact with habitat
           data = FIDbin10, 
           family = binomial) 
summary(h66)
anova(h66, test="Chisq")

h77 <- glm(HideZero2 ~ Habitat ,  # Running and method might interact with habitat
           data = FIDbin10, 
           family = binomial) # Habitat the only thing significant. chi-square test

h77 <- anova(h77, test="Chisq") 
summary(h77) 

#

summary(h77)
cq2 <- chisq.test(FIDbin10$Habitat,FIDbin10$HideZero2)


# 

hidenot2 <- as.data.frame(table(FIDbin10$Habitat,FIDbin10$HideZero2))

# Plotting:
ggplot(hidenot2, aes(x=  factor(Var1, level = c('Urban', 'Open field', 'Water pond')), y=Freq, fill = Var2))+
  geom_bar(stat='identity', position=position_dodge())+
  ylab("Freq") + xlab("Habitat") + labs( fill = "Hide or not 2")+
  geom_text(aes(label=Freq), size=3.5, vjust=-0.4,
            position = position_dodge (0.8))+
  theme_minimal()

# How much percentage did lapwings hide:

12/21*100 #<- # 57.14% in urban did not hide
6/12*100 #<- # 50% in open field did not hide
12/25*100 #<- # 48% of lapwings from water pond hide the second time

# Second running:

table(FIDbin11$RunZero2,FIDbin11$Habitat,FIDbin11$FIDmethod) #seems no diffrences

names(FIDbin11)
# creating models:
# With random effect
h001 <- glmer(
  RunZero2 ~ (1 | MapLoc), 
  data = FIDbin11, 
  family = binomial) # boundary fit - cant use random effect


anova(h22, test="Chisq") #Removing interaction with running, WD, FID, month

h333 <- glm(RunZero2 ~ SecHideOrNot+ scale(SecStD) + scale(SecFlee) + scale(Minuteday) + Habitat + FIDmethod,  # Running and method might interact with habitat
            data = FIDbin11, 
            family = binomial) # Method and habitat
summary(h333)
anova(h333, test="Chisq") # Minutes of day and sFID are least significant

h444 <- glm(RunZero2 ~ SecHideOrNot+ scale(SecStD) + scale(SecFlee) + Habitat + FIDmethod,  # Running and method might interact with habitat
            data = FIDbin11, 
            family = binomial) # Method and habitat
summary(h444)
anova(h444, test="Chisq") # Minutes of day and sFID are least significant

h555 <- glm(RunZero2 ~ SecHideOrNot + scale(SecFlee) + Habitat + FIDmethod,  # Running and method might interact with habitat
            data = FIDbin11, 
            family = binomial) # Method and habitat
summary(h555)
anova(h555, test="Chisq") # Minutes of day and sFID are least significant

h555 <- glm(RunZero2 ~ SecHideOrNot + scale(SecFlee) + Habitat + FIDmethod,  # Running and method might interact with habitat
            data = FIDbin11, 
            family = binomial) # Method and habitat
table(FIDbin11$RunFly2, FIDbin11$Habitat)

h66 <- glm(HideZero2 ~ RunFly2 + Habitat ,  # Running and method might interact with habitat
           data = FIDbin10, 
           family = binomial) 
summary(h66)
anova(h66, test="Chisq")

h77 <- glm(HideZero2 ~ Habitat ,  # Running and method might interact with habitat
           data = FIDbin10, 
           family = binomial) # Habitat the only thing significant. chi-square test

h77 <- anova(h77, test="Chisq") 
summary(h77) 

11/89*100 - 100
#

summary(h77)
cq2 <- chisq.test(FIDbin10$Habitat,FIDbin10$HideZero2)


# 

hidenot2 <- as.data.frame(table(FIDbin10$Habitat,FIDbin10$HideZero2))

# Plotting:
ggplot(hidenot2, aes(x=  factor(Var1, level = c('Urban', 'Open field', 'Water pond')), y=Freq, fill = Var2))+
  geom_bar(stat='identity', position=position_dodge())+
  ylab("Freq") + xlab("Habitat") + labs( fill = "Hide or not 2")+
  geom_text(aes(label=Freq), size=3.5, vjust=-0.4,
            position = position_dodge (0.8))+
  theme_minimal()

# How much percentage did lapwings hide:

12/21*100 #<- # 57.14% in urban did not hide
6/12*100 #<- # 50% in open field did not hide
12/25*100 #<- # 48% of lapwings from water pond hide the second time




model = glm(HideZero ~ RunFly1 * scale(EsD), data = FIDbin7, family="binomial")
model3 = glm(HideZero ~ RunFly1 + EsD, data = FIDbin7, family="binomial")

anova(h5, test="Chisq") 
anova(model2, test="Chisq") 
anova(model3, test="Chisq") 

# The Desicion to hide was not effected by habitat or method - but by the EsD (DF = 145, P < 0.05), if it run or flow 
# (DF = 147, p = 0.018) and in a small amount by the SD (DF = 146, p = 0.035)

FIDbin9 <- subset(FIDbin7, EsD >= 0 )

#model2

summary(model2)

# new.data for predicting new observations

# Creating binary model:

SEsE <- scale(FIDbin9$EsD)
summary(SEsE)
summary(FIDbin9$EsD)
summary(FIDbin9$RunFly1)

new.data <- expand.grid(EsD    = seq(5, 270, length = n), 
                        RunFly1 = c("Run", "Fly"))

preds <- predict(model3, newdata = new.data, type = 'response',se = TRUE)

new.data$pred.full <- preds$fit

new.data$ymin <- new.data$pred.full - 2*preds$se.fit
new.data$ymax <- new.data$pred.full + 2*preds$se.fit  

# Plotting
ggplot(FIDbin9,aes(x = EsD, y = HideZero, col = RunFly1, fill = RunFly1)) + 
  geom_point(size = 2.5) + 
  geom_ribbon(data = new.data,aes(y = pred.full, ymin = ymin, ymax = ymax),alpha = 0.7) +
  xlab("DF") +
  ylab('Hide/Did not hide') + labs(col = "Run or fly", fill = "Run or fly")

ggplot(FIDbin9,aes(x = EsD, y = HideZero, col = RunFly1)) + 
  geom_point(size = 2.5) + 
  geom_smooth(data = new.data,aes(y = pred.full, ymin = ymin, ymax = ymax),alpha = 0.7) +
  xlab("DF") +
  ylab('Hide/Did not hide') + labs(col = "Run or fly")

### Running effect:

names(FIDbin8)

r1 <- glm(
  RunZero ~ Habitat*FIDmethod + month + FIDmethod + scale(StD) + scale(FirstFlee) + scale(Minuteday) + Groups, 
  data = FIDbin8, 
  family = binomial) #boundary fit

summary(r1)
anova(r1, test="Chisq") #Removing everything that is not significant

r2 <- glm(
  RunZero ~ Habitat, 
  data = FIDbin8, 
  family = binomial) #boundary fit

anova(r2, test="Chisq") # The only thing that affect running is the habitat

r3 <- glm(
  RunZero ~ month, 
  data = FIDbin8, 
  family = binomial) #boundary fit
anova(r3, test="Chisq") # The only thing that affect running is the habitat

summary(r2)
cq <- chisq.test(FIDbin8$Habitat,FIDbin8$RunFly1)

# The decision to fly or not was more affected in urban habitats (Z = -2.113 +- 0.4497, p = 0.0346), but was 

RunNot <- as.data.frame(table(FIDbin8$Habitat,FIDbin8$RunFly1))

# Plotting:
ggplot(RunNot, aes(x=  factor(Var1, level = c('Urban', 'Open field', 'Water pond')), y=Freq, fill = Var2))+
  geom_bar(stat='identity', position=position_dodge(), colour="black")+
  ylab("Freq") + xlab("Habitat") + labs( fill = "Run or fly")+
  geom_text(aes(label=Freq), size=3.5, vjust=-0.4,
            position = position_dodge (0.8))+
  theme_minimal() + scale_fill_grey(start = 0.3, end = 0.9)

# How much percentage did lapwings run:

21/68*100
8/67*100
5/140*100
##################V 2.6) Is the WD affected by the habitat? Method? WD, Habitat, FIDmethod ####


plot(FIDbin4$WD~FIDbin4$Habitat)
plot(FIDbin4$WD~FIDbin4$FIDmethod)

hist(FIDbin4$WD) # Not normal distributade

shapiro.test(FIDbin4$WD) # Data is normal

ggboxplot(FIDbin4, x = "Habitat", y = "WD", color = "FIDmethod")
ggboxplot(FIDbin4, x = "FIDmethod", y = "WD", color = "Habitat")

ggline(FIDbin4, x = "Habitat", y = "WD", color = "FIDmethod",
       add = c("mean_se", "dotplot"))

ggline(FIDbin4, x = "FIDmethod", y = "WD", color = "Habitat",
       add = c("mean_se", "dotplot"))

res.aov8 <- aov(WD ~ FIDmethod + Habitat, data = FIDbin4)
res.aov9 <- aov(WD ~ FIDmethod + Habitat + FIDmethod:Habitat, data = FIDbin4)

summary(res.aov8)
summary(res.aov9) # Not interact - almost significant # Dont use this model

TukeyHSD(res.aov8, which = "FIDmethod") # No differences in approach in habitats, but between methods
#                  diff      lwr      upr    p adj
#Jackale-Human 11.05824 5.840458 16.27602 3.97e-05
# Small differences between urban and such

plot(res.aov8, 1) # 53, 123 and 292 are outliners
plot(res.aov8, 2)

# Extract the residuals
aov_residuals3 <- residuals(object = res.aov8)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals3 ) # Data is normal


ggplot(data=FIDbin4, aes (x =  FIDmethod, y=WD , fill = factor(Habitat, level = c('Urban', 'Open field', 'Water pond')))) + 
  geom_boxplot(notch = TRUE) +scale_fill_manual(values=c("Grey", "Green", "Blue")) +
  theme_minimal() + ggtitle("WD between habitat and method of approach") +
  xlab("Approach methods") + ylab("WD")  + labs(fill= "Habitat")
##################V 2.7) Is the escape distance affected by flying? or getting to cover? EsD, RunFly1, HideOrNot ####

FIDbin4$HideOrNot[FIDbin4$HideOrNot == "NA"] <- NA
FIDbin4$RunFly1[FIDbin4$RunFly1 == "NA"] <- NA

plot(FIDbin4$EsD~FIDbin4$RunFly1)

table(FIDbin4$HideOrNot)
table(FIDbin4$RunFly1,FIDbin4$HideOrNot)

plot(FIDbin4$EsD~FIDbin4$HideOrNot)

hist(FIDbin4$EsD) # Not normal distributade

shapiro.test(FIDbin4$LogEsD) # Data is normal

ggboxplot(FIDbin4, x = "RunFly1", y = "LogEsD", color = "HideOrNot")
ggboxplot(FIDbin4, x = "HideOrNot", y = "LogEsD", color = "RunFly1")

ggline(FIDbin4, x = "RunFly1", y = "LogEsD", color = "HideOrNot",
       add = c("mean_se", "dotplot"))

ggline(FIDbin4, x = "HideOrNot", y = "LogEsD", color = "RunFly1",
       add = c("mean_se", "dotplot"))


res.aov10 <- aov(LogEsD ~ HideOrNot + RunFly1, data = FIDbin4)
res.aov11 <- aov(LogEsD ~ HideOrNot + RunFly1 + RunFly1:HideOrNot, data = FIDbin4)

summary(res.aov10) # Both are significant for the escape distance
summary(res.aov11) # Not interact - almost significant # Dont use this model

TukeyHSD(res.aov10, which = "HideOrNot") # No differences in approach in habitats, but between methods
#                  diff       lwr       upr   p adj
#Hide-Didn't hide 0.5692001 0.3366129 0.8017873 3.3e-06

TukeyHSD(res.aov10, which = "RunFly1")

#              diff        lwr        upr     p adj
#Run-Fly -0.5145157 -0.8439939 -0.1850376 0.0024256
plot(res.aov10, 1) # 194, 241 and 227 are outliners
plot(res.aov10, 2)

# Extract the residuals
aov_residuals4 <- residuals(object = res.aov10)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals4 ) # Data is normal

FIDbin5 <- na.omit(subset(FIDbin4, select = c(LogEsD,RunFly1, HideOrNot)))
ggplot(data = FIDbin5, aes(x = HideOrNot, y = LogEsD, color =RunFly1)) + 
  geom_boxplot()

# When the lapwings flew they had longer EsD, as well as when they hide

##################V 2.8) Is FID 2 correlate with FID 1? between habitats and methods? FirstFlee, SecFlee, FIDmethod, Habitat ####

hist(FIDbin4$SecFlee)
names(FIDbin4)

FIDbin6 <- subset(FIDbin4,SecFlee >=0 )

# Map locations for FID2:
plot(table(FIDbin6$MapLoc), main="sFID tested locations",
     xlab="Locations after rounded",
     ylab="Time tested")

str(as.factor(FIDbin4$MapLoc))


names(FIDbin6)

modFb0 <- lmer(SecFlee ~ (1 | MapLoc), data = FIDbin6, REML = FALSE)
# The group is not consider this time, as after the first FID they might not be next to him. I added the first FID and if it run or fly the first time. hiding is not considered as if it hide there is not second FID
modFb1 <- lmer(SecFlee ~ (1 | MapLoc) + Habitat + scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD), data = FIDbin6, REML = FALSE)
# Too large.
modFb1 <- lmer(SecFlee ~ (1 | MapLoc) + Habitat + scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee), data = FIDbin6, REML = FALSE)
# WD is problematic. removed it. still singularities... will make a lm instead.

# Added interactions between the habitat, method and first FID
modFb1 <- glm(SecFlee ~ Habitat + scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian",  data = FIDbin6)
modFb2 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb3 <- glm(SecFlee ~ Habitat                 + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb4 <- glm(SecFlee ~ Habitat + scale(SecStD)             + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb5 <- glm(SecFlee ~ Habitat + scale(SecStD) + FIDmethod                + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb6 <- glm(SecFlee ~ Habitat + scale(SecStD) + FIDmethod + scale(month)                    + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb7 <- glm(SecFlee ~ Habitat + scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday)           + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb8 <- glm(SecFlee ~ Habitat + scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1                    + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb9 <- glm(SecFlee ~ Habitat + scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee)                + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb10 <-glm(SecFlee ~ Habitat + scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD)                     + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb11 <-glm(SecFlee ~ Habitat + scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat                     + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb12 <-glm(SecFlee ~ Habitat + scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee                       ,family = "gaussian", data = FIDbin6)

AICall <- AIC(modFb1, modFb2,modFb3,modFb4,modFb5,modFb6,modFb7,modFb8,modFb9,modFb10,
              modFb11,modFb12)
AICall <- AICall[order(AICall$AIC),] # Best model is 3 
AICall # modFb2- habitat alone have no significant for secnond FID

modFb21 <- glm(SecFlee ~                           FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb22 <- glm(SecFlee ~           scale(SecStD)             + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb23 <- glm(SecFlee ~           scale(SecStD) + FIDmethod                + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb24 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month)                    + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb25 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday)           + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb26 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1                    + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb27 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee)                + FIDmethod*Habitat + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb28 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD)                     + Habitat*FirstFlee + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb29 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat                     + FIDmethod*FirstFlee ,family = "gaussian", data = FIDbin6)
modFb20 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD) + FIDmethod*Habitat + Habitat*FirstFlee                       ,family = "gaussian", data = FIDbin6)

AICall <- AIC(modFb1, modFb2,modFb3,modFb4,modFb5,modFb6,modFb7,modFb8,modFb9,modFb10,
              modFb11,modFb12,modFb21,modFb22,modFb23,modFb24,modFb25,modFb26,modFb27,modFb28,
              modFb29,modFb20)
AICall <- AICall[order(AICall$AIC),] # Best model is 3 
AICall # modFb29;modFb28;modFb20- Habitat*FID, Method*Habitat, Fid*Method - all interactions need to be removed

modFb200 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD)                        ,family = "gaussian", data = FIDbin6)
modFb201 <- glm(SecFlee ~                           FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD)                        ,family = "gaussian", data = FIDbin6)
modFb202 <- glm(SecFlee ~           scale(SecStD)             + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD)                        ,family = "gaussian", data = FIDbin6)
modFb203 <- glm(SecFlee ~           scale(SecStD) + FIDmethod                + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD)                        ,family = "gaussian", data = FIDbin6)
modFb204 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month)                    + RunFly1 + scale(FirstFlee) + scale(secWD)                        ,family = "gaussian", data = FIDbin6)
modFb205 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday)           + scale(FirstFlee) + scale(secWD)                        ,family = "gaussian", data = FIDbin6)
modFb206 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1                    + scale(secWD)                        ,family = "gaussian", data = FIDbin6)
modFb207 <- glm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee)                                       ,family = "gaussian", data = FIDbin6)

AICall <- AIC(modFb1, modFb2,modFb3,modFb4,modFb5,modFb6,modFb7,modFb8,modFb9,modFb10,
              modFb11,modFb12,modFb21,modFb22,modFb23,modFb24,modFb25,modFb26,modFb27,modFb28,
              modFb29,modFb20,modFb200,modFb201,modFb202,modFb203,modFb204,modFb205,modFb206,modFb207)
AICall <- AICall[order(AICall$AIC),] # Best model is 205 
AICall # Everything but the habitat affect the secFID... so I will take the next model to test after

# This will be my base model: modFb200 <- lm(SecFlee ~           scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD)                        , data = FIDbin6)
modFb2000 <- glm(SecFlee ~  scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD),family = "gaussian", data = FIDbin6)
modFb2001 <- glm(SecFlee ~  scale(SecStD)              + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD),family = "gaussian", data = FIDbin6)
modFb2002 <- glm(SecFlee ~  scale(SecStD) + FIDmethod                 + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD),family = "gaussian", data = FIDbin6)
modFb2003 <- glm(SecFlee ~  scale(SecStD) + FIDmethod + scale(month)                     + RunFly1 + scale(FirstFlee) + scale(secWD),family = "gaussian", data = FIDbin6)
modFb2004 <- glm(SecFlee ~  scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday)            + scale(FirstFlee) + scale(secWD),family = "gaussian", data = FIDbin6)
modFb2005 <- glm(SecFlee ~  scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1                     + scale(secWD),family = "gaussian", data = FIDbin6)
modFb2006 <- glm(SecFlee ~  scale(SecStD) + FIDmethod + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee)                ,family = "gaussian", data = FIDbin6)

AICall2 <- AIC(modFb2000,modFb2001,modFb2003,modFb2004,modFb2005,modFb2006)
AICall2 <- AICall2[order(AICall2$AIC),] 
AICall2

modFb20011 <- glm(SecFlee ~                               scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD),family = "gaussian", data = FIDbin6)
modFb20012 <- glm(SecFlee ~  scale(SecStD)                             + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD),family = "gaussian", data = FIDbin6)
modFb20013 <- glm(SecFlee ~  scale(SecStD)              + scale(month)                    + RunFly1 + scale(FirstFlee) + scale(secWD),family = "gaussian", data = FIDbin6)
modFb20014 <- glm(SecFlee ~  scale(SecStD)              + scale(month) + scale(Minuteday)           + scale(FirstFlee) + scale(secWD),family = "gaussian", data = FIDbin6)
modFb20015 <- glm(SecFlee ~  scale(SecStD)              + scale(month) + scale(Minuteday) + RunFly1                    + scale(secWD),family = "gaussian", data = FIDbin6)
modFb20016 <- glm(SecFlee ~  scale(SecStD)              + scale(month) + scale(Minuteday) + RunFly1 + scale(FirstFlee)               ,family = "gaussian", data = FIDbin6)

AICall2 <- AIC(modFb2000,modFb2001,modFb2003,modFb2004,modFb2005,modFb2006,
               modFb20011,modFb20012,modFb20013,modFb20014,modFb20015,modFb20016)
AICall2 <- AICall2[order(AICall2$AIC),] 
AICall2

modFb200121 <- glm(SecFlee ~                                           + scale(Minuteday) + RunFly1 + scale(FirstFlee) + scale(secWD),family = "gaussian", data = FIDbin6)
modFb200122 <- glm(SecFlee ~  scale(SecStD)                                                + RunFly1 + scale(FirstFlee) + scale(secWD),family = "gaussian", data = FIDbin6)
modFb200123 <- glm(SecFlee ~  scale(SecStD)                             + scale(Minuteday)           + scale(FirstFlee) + scale(secWD),family = "gaussian", data = FIDbin6)
modFb200124 <- glm(SecFlee ~  scale(SecStD)                             + scale(Minuteday) + RunFly1                    + scale(secWD),family = "gaussian", data = FIDbin6)
modFb200125 <- glm(SecFlee ~  scale(SecStD)                             + scale(Minuteday) + RunFly1 + scale(FirstFlee)               ,family = "gaussian", data = FIDbin6)


AICall2 <- AIC(modFb2000,modFb2001,modFb2003,modFb2004,modFb2005,modFb2006,
               modFb200121,modFb200122,modFb200123,modFb200124,modFb200125,modFb20011,modFb20012,modFb20013,modFb20014,modFb20015,modFb20016)
AICall2 <- AICall2[order(AICall2$AIC),] 
AICall2

modFb2001241 <- glm(SecFlee ~                                             scale(Minuteday) + RunFly1                    + scale(secWD),family = "gaussian", data = FIDbin6)
modFb2001242 <- glm(SecFlee ~  scale(SecStD)                                                + RunFly1                    + scale(secWD),family = "gaussian", data = FIDbin6)
modFb2001243 <- glm(SecFlee ~  scale(SecStD)                             + scale(Minuteday)                              + scale(secWD),family = "gaussian", data = FIDbin6)
modFb2001244 <- glm(SecFlee ~  scale(SecStD)                             + scale(Minuteday) + RunFly1                                  ,family = "gaussian", data = FIDbin6)

AICall2 <- AIC(modFb2000,modFb2001,modFb2003,modFb2004,modFb2005,modFb2006,
               modFb2001241,modFb2001242,modFb2001243,modFb2001244,modFb200121,modFb200122,modFb200123,modFb200124,modFb200125,modFb20011,modFb20012,modFb20013,modFb20014,modFb20015,modFb20016)
AICall2 <- AICall2[order(AICall2$AIC),] 
AICall2

modFb20012431 <- glm(SecFlee ~                                            + scale(Minuteday)                              + scale(secWD),family = "gaussian", data = FIDbin6)
modFb20012432 <- glm(SecFlee ~  scale(SecStD)                                                                             + scale(secWD),family = "gaussian", data = FIDbin6)
modFb20012433 <- glm(SecFlee ~  scale(SecStD)                             + scale(Minuteday)                                            ,family = "gaussian", data = FIDbin6)

AICall2 <- AIC(tryonemore,modFb1, modFb2,modFb3,modFb4,modFb5,modFb6,modFb7,modFb8,modFb9,modFb10,
               modFb11,modFb12,modFb22,modFb23,modFb24,modFb25,modFb26,modFb28,
               modFb29,modFb20,modFb201,modFb203,modFb20012431,modFb20012432,modFb20012433,modFb2000,modFb2001,modFb2003,modFb2004,modFb2005,modFb2006,
               modFb2001241,modFb2001242,modFb2001243,modFb2001244,modFb200121,modFb200122,modFb200123,modFb200124,modFb200125,modFb20011,modFb20012,modFb20013,modFb20014,modFb20015,modFb20016)
AICall2 <- AICall2[order(AICall2$AIC),] 
AICall2 # Best model: WD,SD and day time (?)

modFb200124321 <- glm(SecFlee ~                                                                                            + scale(secWD),family = "gaussian", data = FIDbin6)
modFb200124322 <- glm(SecFlee ~  scale(SecStD)                                                                                           ,family = "gaussian", data = FIDbin6)
modFb200124323 <- glm(SecFlee ~  scale(SecStD)*scale(secWD),family = "gaussian", data = FIDbin6)

AICall2 <- AIC(modFb200124321,modFb200124322, modFb1, modFb2,modFb3,modFb4,modFb5,modFb6,modFb7,modFb8,modFb9,modFb10,
               modFb11,modFb12,modFb22,modFb23,modFb24,modFb25,modFb26,modFb28,
               modFb29,modFb20,modFb201,modFb203,modFb20012431,modFb20012432,modFb20012433,modFb2000,modFb2001,modFb2003,modFb2004,modFb2005,modFb2006,
               modFb2001241,modFb2001242,modFb2001243,modFb2001244,modFb200121,modFb200122,modFb200123,modFb200124,modFb200125,modFb20011,modFb20012,modFb20013,modFb20014,modFb20015,modFb20016)
AICall2 <- AICall2[order(AICall2$AIC),] 
AICall2 # Best model: WD,SD and day time (?)

AICall2$delta <-  AICall2$AIC - (-4505.3661)
#install.packages("rgl")
AICall2$weights<- akaike.weights(AICall2$AIC)$weights

AICall2$weights <- format(round(AICall2$weights, 3), nsmall = 3)
AICall2$delta <- format(round(AICall2$delta, 3), nsmall = 3)

AICall2$weights <- as.numeric(AICall2$weights)
AICall2$weights <- 100*AICall2$weights
AICall2 # The best model to describe is SD with WD, and the time have a small affect to.
summary(modFb20012432)
summary(modFb2001243)
tryonemore <- glm(SecFlee ~ scale(SecStD) + scale(FirstFlee) + scale(secWD) + scale(Minuteday), family = "gaussian", data = FIDbin6) # ok so no
2.486e+16 -1.046e+16
modFb2001243m <- lm(SecFlee ~  scale(SecStD)                             + scale(Minuteday)                             , data = FIDbin6)
2e-16
2.699e+01
2.755e+01
-2.284e+01
summary(modFb2001243m)
# FID 2 is mostly affected by SecSD, Day time and SecWD
#SecFlee, SecStD, Minuteday, secWD
options(scipen = 5)
(1.086e-15)
par(mfrow = c(2, 2))

ggplot(FIDbin6,aes(y = SecFlee, x = SecStD)) + 
  #facet_grid(~Minuteday) +
  geom_point(alpha = 0.8) +
  stat_smooth(method = "lm") +
  #geom_hline(yintercept=0, linetype="dashed") +
  theme_bw()

ggplot(FIDbin4,aes(y = SecFlee, x = FirstFlee)) + 
  #facet_grid(~Minuteday) +
  geom_point(alpha = 0.8) +
  stat_smooth(method = "lm") +
  #geom_hline(yintercept=0, linetype="dashed") +
  theme_bw()

#just correlation:
FID12 <- lm(SecFlee ~ FirstFlee, data = FIDbin4 )
summary(FID12) # It is correlated but not the best fitted model

ggplot(FIDbin6,aes(y = SecFlee, x = secWD, col = Habitat)) + geom_boxplot() +
  theme_bw()
################## 2.9) Repeat the same tests on FID 2  ####

# I will test repeatability in FID, EsD (If possible), WD, fly and hide between habitat and method of approach

#1) FID~FID2:
# I need to create a new data.frame with Test1 + Test2 and a new column name test number

# Create a repeat of 1 for the first fid:

FIDbin4$repeats <- 1
FIDbin4$Ind <- seq.int(nrow(FIDbin4))
FIDbin4$AllFID <- FIDbin4$FirstFlee
FIDbin4$AllStD <- FIDbin4$StD
FIDbin4$AllHide <- FIDbin4$HideOrNot
FIDbin4$AllEsD <- FIDbin4$EsD
FIDbin4$AllFly <- FIDbin4$RunFly1
FIDbin4$AllWD <- FIDbin4$WD

FID2 <- subset(FIDbin4, SecFlee >= 0)

FID2$repeats <- 2
FID2$AllFID <- FID2$SecFlee
FID2$AllStD <- FID2$SecStD
FID2$AllHide <- FID2$SecHideOrNot
FID2$AllEsD <- FID2$SecEsD
FID2$AllFly <- FID2$RunFly2
FID2$AllWD <- FID2$secWD
FID2$Ind

individuals <- unique(FID2$Ind)
newFID <- list()
for (i in individuals) {
  re <- subset(FIDbin4,Ind == i )
  newFID[[i]] <- re
}
newFID <- do.call(rbind.data.frame, newFID)

FIDrepeats <- rbind(FID2,newFID)

names(FIDrepeats)

drops <- c("secWD","RunFly2","SecEsD","SecHideOrNot","SecStD","SecFlee","WD","RunFly1",
           "EsD", "HideOrNot", "StD", "FirstFlee")
FIDrepeats <- FIDrepeats[ , !(names(FIDrepeats) %in% drops)]
# 

library(lme4)
library(rptR)

# FID
names(FIDrepeats)
Rep1 <- rpt(AllFID ~ (1 | Ind), grname = c("Ind"), data = FIDrepeats, datatype = "Gaussian", nboot = 1000, 
            npermut = 0)
plot(Rep1, cex.main = 1)

summary(Rep1)

Rep2 <- rpt(AllFID ~ (1 | Ind) + (1|MapLoc), grname = c("Ind", "MapLoc"), data = FIDrepeats, datatype = "Gaussian", nboot = 1000, 
            npermut = 0)

summary(Rep2)

Rep3 <- rpt(AllFID ~ (1 | Ind) + Habitat + FIDmethod, grname = c("Ind"), data = FIDrepeats, datatype = "Gaussian", nboot = 1000, 
            npermut = 0)

Rep4 <- rpt(AllFID ~ (1 | Ind) + Habitat*FIDmethod, grname = c("Ind"), data = FIDrepeats, datatype = "Gaussian", nboot = 1000, 
            npermut = 0)

Rep5 <- rpt(AllFID ~ (1 | Ind) + Habitat, grname = c("Ind"), data = FIDrepeats, datatype = "Gaussian", nboot = 1000, 
            npermut = 0)

Rep6 <- rpt(AllFID ~ (1 | Ind) + FIDmethod, grname = c("Ind"), data = FIDrepeats, datatype = "Gaussian", nboot = 1000, 
            npermut = 0)

# Only repeatable measure of individuals vs. FID is very repeatable (r = 0.49, se = 0.077, boot = 1000)
ggplot(data = FIDrepeats, aes(x = Ind, y = AllFID, group = FIDmethod)) + geom_boxplot(notch=TRUE)
ggplot(data = FIDrepeats, aes(x = Ind, y = Time, group = Habitat)) + geom_boxplot(notch=FALSE)
ggplot(data = FIDrepeats, aes(x = Ind, y = Time, group = Lapwing)) + geom_boxplot(notch=FALSE)
ggplot(data = FIDrepeats, aes(x = Ind, y = Time, group = Lapwing, color = Lapwing)) + geom_point(size=2,position=position_jitter(w=.2)) +theme_minimal()

names(FIDrepeats)
ggplot(FIDrepeats, aes(x = repeats, y = AllFID, color = FIDmethod)) +
  facet_grid(~Habitat) +
  geom_point(size = 3) + 
  geom_line(aes(group=Ind)) +
  
  theme_minimal()
names(FIDrepeats)
#AllFID ~ Habitat + FIDmethod
Rep2 <- rpt(AllFID ~ (1 | Ind)  + Habitat + FIDmethod , grname = c("Ind"), data = FIDrepeats, datatype = "Gaussian", nboot = 1000, 
            npermut = 0)
summary(Rep2)
Rep3 <- rpt(AllFID ~ (1 | MapLoc)  + Habitat + FIDmethod , grname = c("Lapwing"), data = FIDrepeats, datatype = "Gaussian", nboot = 1000, 
            npermut = 0)

ggplot(FIDrepeats, aes(x = repeats, y = AllFID, color = FIDmethod)) +
  facet_grid(~Habitat) +
  geom_point(size = 3) + 
  geom_line(aes(group=MapLoc)) +
  theme_minimal()

#StD

FIDrepStd <- subset(FIDrepeats,AllStD >= 0 )
RepAllStD <- rpt(AllStD ~ (1 | Ind)  , grname = c("Ind"), data = FIDrepStd, datatype = "Gaussian", nboot = 1000, 
                 npermut = 0)

plot(RepAllStD) # 0.343

ggplot(FIDrepStd, aes(x = repeats, y = AllStD, color = FIDmethod)) +
  facet_grid(~Habitat) +
  geom_point(size = 3) + 
  geom_line(aes(group=MapLoc)) +
  theme_minimal()

#WD

FIDrepWD <- subset(FIDrepeats,AllWD >= 0 )
RepAllWD <- rpt(AllWD ~ (1 | Ind)  , grname = c("Ind"), data = FIDrepWD, datatype = "Gaussian", nboot = 1000, 
                npermut = 0)

plot(RepAllWD) # Nope

# EsD
names(FIDrepeats)
FIDrepEsD <- subset(FIDrepeats,AllEsD >= 0 )

RepAllEsD <- rpt(AllEsD ~ (1 | Ind)  , grname = c("Ind"), data = FIDrepEsD, datatype = "Gaussian", nboot = 1000, 
                 npermut = 0)
plot(RepAllEsD) # Yes 0.371

ggplot(FIDrepEsD, aes(x = repeats, y = AllEsD, color = FIDmethod)) +
  facet_grid(~Habitat) +
  geom_point(size = 3) + 
  geom_line(aes(group=MapLoc)) +
  theme_minimal()

# EsD
names(FIDrepeats)

# For testing what affect hiding:
FIDrepeats$HideZero2[FIDrepeats$AllHide == "Hide"] <- 0
FIDrepeats$HideZero2[FIDrepeats$AllHide == "Didn't hide"] <- 1

# For testing what affect running:
FIDrepeats$RunZero2[FIDrepeats$AllFly == "Run"] <- 0
FIDrepeats$RunZero2[FIDrepeats$AllFly == "Fly"] <- 1

FIDrepHide <- subset(FIDrepeats,HideZero2 >= 0 )

RepAllHide <- rpt(HideZero2 ~ (1 | Ind)  , grname = c("Ind"), data = FIDrepHide, datatype = "Binary", nboot = 1000, 
                  npermut = 0)

plot(RepAllHide) # Of course it is not, if one hide it wont have a second test...

FIDrepRun <- subset(FIDrepeats,RunZero2 >= 0 )

RepAllHide <- rpt(RunZero2 ~ (1 | Ind)  , grname = c("Ind"), data = FIDrepRun, datatype = "Binary", nboot = 1000, 
                  npermut = 0)
plot(RepAllHide) # Not repeatable

################## 2.10) What affect the FID the most? FirstFlee~ Date, Time, Habitat, GroupSize, StD, FIDmethod, random(MapLoc), month,Groups,Minuteday,WD  ####
# starting with GLMM for first flee:
names(FIDbin4)
scale(FIDbin4$month)
hist(FIDbin4$WD)
mod1 <- lmer(FirstFlee ~ (1 | MapLoc)                                                                                   , data = FIDbin4, REML = FALSE)
mod2 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD) + FIDmethod + scale(month) + Groups + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod3 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(WD) + FIDmethod + scale(month) + Groups + scale(Minuteday) , data = FIDbin4, REML = FALSE)

AICall <- AIC(mod1,mod2,mod3)
AICall <- AICall[order(AICall$AIC),] 
AICall
summary(mod2)

mod21 <- lmer(FirstFlee ~ (1 | MapLoc) +           scale(StD) + FIDmethod + scale(month) + Groups + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod22 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat              + FIDmethod + scale(month) + Groups + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod23 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD)             + scale(month) + Groups + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod24 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD) + FIDmethod                + Groups + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod25 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD) + FIDmethod + scale(month)          + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod26 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD) + FIDmethod + scale(month) + Groups                    , data = FIDbin4, REML = FALSE)

AICall <- AIC(mod1,mod2,mod3, mod21, mod22, mod23, mod24, mod25, mod26)
AICall <- AICall[order(AICall$AIC),] 
AICall
summary(mod25)

mod251 <- lmer(FirstFlee ~ (1 | MapLoc) +           scale(StD) + FIDmethod + scale(month)          + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod252 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat              + FIDmethod + scale(month)          + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod253 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD)             + scale(month)          + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod254 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD) + FIDmethod                         + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod255 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD) + FIDmethod                                            , data = FIDbin4, REML = FALSE)
mod241 <- lmer(FirstFlee ~ (1 | MapLoc)           + scale(StD) + FIDmethod                + Groups + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod242 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat              + FIDmethod                + Groups + scale(Minuteday) , data = FIDbin4, REML = FALSE)
#mod243 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD)                            + Groups + scale(Minuteday) , data = FIDbin4, REML = FALSE)
#mod244 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD)                            + Groups + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod261 <- lmer(FirstFlee ~ (1 | MapLoc)           + scale(StD) + FIDmethod + scale(month) + Groups                    , data = FIDbin4, REML = FALSE)
mod262 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat              + FIDmethod + scale(month) + Groups                    , data = FIDbin4, REML = FALSE)
#mod263 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD)             + scale(month) + Groups                    , data = FIDbin4, REML = FALSE)
mod264 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD) + FIDmethod                + Groups                    , data = FIDbin4, REML = FALSE)

# Group together with StD when there is no time dont work
AICall <- AIC(mod264,mod262,mod261,mod242,mod241,
              mod255,mod254,mod253,mod252,mod251,
              mod1,mod2,mod3, mod21, mod22, mod23, mod24, mod25, mod26)

AICall <- AICall[order(AICall$AIC),] 
AICall #255,254

mod2551 <- lmer(FirstFlee ~ (1 | MapLoc)           + scale(StD) + FIDmethod                                            , data = FIDbin4, REML = FALSE)
mod2552 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat              + FIDmethod                                            , data = FIDbin4, REML = FALSE)
#mod2553 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD)                                                        , data = FIDbin4, REML = FALSE)
mod2541 <- lmer(FirstFlee ~ (1 | MapLoc)           + scale(StD) + FIDmethod                         + scale(Minuteday) , data = FIDbin4, REML = FALSE)
mod2542 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat              + FIDmethod                         + scale(Minuteday) , data = FIDbin4, REML = FALSE)
#mod2543 <- lmer(FirstFlee ~ (1 | MapLoc) + Habitat + scale(StD)                                     + scale(Minuteday) , data = FIDbin4, REML = FALSE)

AICall <- AIC(mod2551,mod2552,mod2541,mod2542,mod264,mod262,mod261,mod242,mod241,
              mod255,mod254,mod253,mod252,mod251,
              mod1,mod2,mod3, mod21, mod22, mod23, mod24, mod25, mod26)

AICall <- AICall[order(AICall$AIC),] 
AICall #255,254

summary(mod255)

AICall$delta <-  AICall$AIC - 2374.065
#install.packages("rgl")
AICall$weights<- akaike.weights(AICall$AIC)$weights

AICall$weights <- format(round(AICall$weights, 3), nsmall = 3)
AICall$delta <- format(round(AICall$delta, 3), nsmall = 3)

AICall$weights <- as.numeric(AICall$weights)
AICall$weights <- 100*AICall$weights
AICall
summary(AICall)
summary(mod255)

library(lme4)
library(effects)

ggplot(FIDbin4,aes(FirstFlee, StD, group=MapLoc, col=Habitat, shape=FIDmethod )) + 
  #facet_grid(~N) +
  geom_line(aes(y=FirstFlee, lty=StD), size=0.8) +
  geom_point(alpha = 0.3) + 
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw()

ggplot(FIDbin4,aes(y = FirstFlee, x = StD)) + 
  facet_grid(~Habitat) +
  geom_point(alpha = 0.8) +
  stat_smooth(method = "lm",aes(fill = FIDmethod,  col = FIDmethod)) +
  #geom_hline(yintercept=0, linetype="dashed") +
  theme_bw()

ggplot(FIDbin4,aes(y = FirstFlee, x = StD, group = Habitat, fill = Habitat)) + 
  facet_grid(~FIDmethod) +
  geom_point(alpha = 0.8) +
  stat_smooth(method = "lm", colour="black") + 
  scale_fill_manual(values = c("green", "grey", "blue")) +
  theme_bw()




