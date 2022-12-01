##HMM Models##

## In this file, we provide the R code used for our manuscript "Behaviorally Mediated Coexistence of Ocelots, Bobcats, and Coyotes Using Hidden Markov Models" 
## published in Oikos. The following code prepares data for a hidden Markov analysis, fits the HMM, performs a resource selection function to examine 
## differences in habitat selection across behaviors, and examines temporal overlap in behavior within and across species. Code to create the plots in the manuscript
## is included as well. The code provided is an example for one species (ocelots), we repeated this same process for bobcats and coyotes. 



#### Load Packages####
library(moveHMM)
library(amt)
library(tibble)
library(activity)
library(mapview)
library(sp)
library(adehabitatHR)
library(sf)
library(raster)
library(SDraw)
library(maptools)
library(lme4)
library(MuMIn)
library(activity)
library(card)
library(terra)
library(overlap)
library(circular)
library(scales)




ocelots <- read.csv("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Ocelot & Bobcat Data/OcelotsAnalysisData.csv")




ocelots <- ocelots[!is.na(ocelots$Latitude),]

#Convert from Lat/Long to UTM
cord.dec2 <- SpatialPoints(cbind(ocelots$Longitude, ocelots$Latitude),
                           proj4string = CRS("+proj=longlat +datum=WGS84"))
cord.UTM2 <- spTransform(cord.dec2, CRS("+proj=utm +zone=14N +datum=WGS84"))
cord.UTM2 <- as.data.frame(cord.UTM2)
#Add UTM Columns
ocelots$Easting <- cord.UTM2$coords.x1
ocelots$Northing <- cord.UTM2$coords.x2
head(ocelots)



ocelots$Easting <- ocelots$Easting/1000
ocelots$Northing <- ocelots$Northing/1000


head(ocelots)
table(ocelots$ID)

##Convert Date & Time into One Proper Format
ocelots <- ocelots[order(ocelots$ID,ocelots$Date,ocelots$Time),]

ocelots$Time <- format(strptime(ocelots$Time, "%H:%M:%S"), "%H:%M:%S")

ocelots$Date <- as.Date(ocelots$Date, "%m/%d/%Y")


#Combine Date/Time into one column
date_info <- with(ocelots, paste((ocelots$Date), ocelots$Time))
date_time <- strptime(date_info, "%Y-%m-%d %H:%M:%S")
#as.data.frame(date_time)
ocelots$Date2 <- date_time
ocelots$Date3 <- as.POSIXct(ocelots$Date2)
ocelots$Date_numeric <- as.numeric(ocelots$Date3)


ocelots <- ocelots[order(ocelots$ID,ocelots$Date_numeric),]

#Populate Time lag column
ocelots$Time.lag=(ave(ocelots$Date_numeric, ocelots$ID, FUN=function(x) c(0, diff(x)))/60)

ocelots$DT <-as.POSIXct(ocelots$Date2, tz = "America/Chicago")
ocelots <- ocelots[order(ocelots$ID, ocelots$DT),]
ocelots$Date=as.Date(ocelots$Date)


ocelots$Time.lag <- round(ocelots$Time.lag)

ocelots <- ocelots[ocelots$Time.lag < 123,]
#ocelots <- ocelots[ocelots$Step < 1,]


#ocelots <- ocelots[!(is.na(ocelots$Steps)),]
tail(ocelots)
head(ocelots)
ocelots <- ocelots[!(is.na(ocelots$Time.lag)),]


ochmm <- prepData(ocelots, type = "UTM", coordNames = c("Easting", "Northing"))

plot(ochmm)




##Examine variable means/distributions
mean(ochmm$step, na.rm = T)
sd <- sd(ochmm$step, na.rm = T)
n <- nrow(ochmm)
se <- sd/sqrt(n)

mean(ochmm$angle, na.rm = T)
hist(ochmm$step)
hist(ochmm$angle)

max(ochmm$step, na.rm = T)
min(ochmm$step, na.rm = T)
##Standardize Temp as a test covariate





####Candidate Models####
##2 State##
mu2a <- c(0.1,0.9) # step mean (two parameters: one for each state)
sigma2a <- c(0.1,0.5) # step SD
zeromass2a <- c(0.1,0.05) # step zero-mass
stepPar2a <- c(mu2a,sigma2a, zeromass2a)

angleMean2a <- c(pi,0) # angle mean
kappa2a <- c(1,1) # angle concentration
anglePar2a <- c(angleMean2a,kappa2a)

m2a <- fitHMM(data=ochmm,nbStates=2,stepPar0=stepPar2a,
              anglePar0=anglePar2a)
##########################################
mu2b <- c(0.01,0.8) # step mean (two parameters: one for each state)
sigma2b <- c(0.1,0.4) # step SD
zeromass2b <- c(0.1,0.05) # step zero-mass
stepPar2b <- c(mu2b,sigma2b, zeromass2b)

angleMean2b <- c(pi,0) # angle mean
kappa2b <- c(1,1) # angle concentration
anglePar2b <- c(angleMean2b,kappa2b)

m2b <- fitHMM(data=ochmm,nbStates=2,stepPar0=stepPar2b,
              anglePar0=anglePar2b)
##########################################
mu2c <- c(0.05,1) # step mean (two parameters: one for each state)
sigma2c <- c(0.05,0.4) # step SD
zeromass2c <- c(0.1,0.05) # step zero-mass
stepPar2c <- c(mu2c,sigma2c, zeromass2c)

angleMean2c <- c(pi,0) # angle mean
kappa2c <- c(1,1) # angle concentration
anglePar2c <- c(angleMean2c,kappa2c)

m2c <- fitHMM(data=ochmm,nbStates=2,stepPar0=stepPar2c,
              anglePar0=anglePar2c)
##########################################

mu2d <- c(0.1,0.7) # step mean (two parameters: one for each state)
sigma2d <- c(0.05,0.4) # step SD
zeromass2d <- c(0.1,0.05) # step zero-mass
stepPar2d <- c(mu2d,sigma2d, zeromass2d)

angleMean2d <- c(pi,0) # angle mean
kappa2d <- c(1,1) # angle concentration
anglePar2d <- c(angleMean2d,kappa2d)

m2d <- fitHMM(data=ochmm,nbStates=2,stepPar0=stepPar2d,
              anglePar0=anglePar2d)
##########################################

mu2e <- c(0.05,0.75) # step mean (two parameters: one for each state)
sigma2e <- c(0.02,0.3) # step SD
zeromass2e <- c(0.1,0.05) # step zero-mass
stepPar2e <- c(mu2e,sigma2e, zeromass2e)

angleMean2e <- c(pi,0) # angle mean
kappa2e <- c(1,1) # angle concentration
anglePar2e <- c(angleMean2e,kappa2e)

m2e <- fitHMM(data=ochmm,nbStates=2,stepPar0=stepPar2e,
              anglePar0=anglePar2e)

##########################################

mu2f <- c(0.15,0.9) # step mean (two parameters: one for each state)
sigma2f <- c(0.1,0.5) # step SD
zeromass2f <- c(0.1,0.05) # step zero-mass
stepPar2f <- c(mu2f,sigma2f, zeromass2f)

angleMean2f <- c(pi,0) # angle mean
kappa2f <- c(1,1) # angle concentration
anglePar2f <- c(angleMean2f,kappa2f)

m2f <- fitHMM(data=ochmm,nbStates=2,stepPar0=stepPar2f,
              anglePar0=anglePar2f)

###########################################


##Residuals plot##
pr <- pseudoRes(m2a)
# plot pseudo-residuals
plotPR(m2a)


##3 State##
##Package defaults

###########################################################

mu3b <- c(0.01,0.3,.7)
sigma3b <- c(0.01,0.2,0.5)
zeromass3b <- c(0.05,0.0001,0.0001)
stepPar3b <- c(mu3b,sigma3b,zeromass3b)
angleMean3b <- c(pi,pi,0)
kappa3b <- c(1,1,1)
anglePar3b <- c(angleMean3b,kappa3b)

m3b <- fitHMM(data=ochmm,nbStates=3,stepPar0=stepPar3b,
              anglePar0=anglePar3b)
###########################################################

mu3c <- c(0.01,0.5,1)
sigma3c <- c(0.01,0.3,0.5)
zeromass3c <- c(0.05,0.0001,0.0001)
stepPar3c <- c(mu3c,sigma3c,zeromass3c)
angleMean3c <- c(pi,pi,0)
kappa3c <- c(1,1,1)
anglePar3c <- c(angleMean3c,kappa3c)

m3c <- fitHMM(data=ochmm,nbStates=3,stepPar0=stepPar3c,
              anglePar0=anglePar3c)

###########################################################

mu3d <- c(0.01,0.25,0.8)
sigma3d <- c(0.01,0.25,0.6)
zeromass3d <- c(0.05,0.0001,0.0001)
stepPar3d <- c(mu3d,sigma3d,zeromass3d)
angleMean3d <- c(pi,pi,0)
kappa3d <- c(1,1,1)
anglePar3d <- c(angleMean3d,kappa3d)

m3d <- fitHMM(data=ochmm,nbStates=3,stepPar0=stepPar3d,
              anglePar0=anglePar3d)

########################################################


mu3e <- c(0.05,0.5,1.0)
sigma3e <- c(0.03,0.3,0.6)
zeromass3e <- c(0.05,0.0001,0.0001)
stepPar3e <- c(mu3e,sigma3e,zeromass3e)
angleMean3e <- c(pi,pi,0)
kappa3e <- c(1,1,1)
anglePar3e <- c(angleMean3e,kappa3e)

m3e <- fitHMM(data=ochmm,nbStates=3,stepPar0=stepPar3e,
              anglePar0=anglePar3e)

############################################################

mu3f <- c(0.05,0.35,0.9)
sigma3f <- c(0.01,0.25,0.6)
zeromass3f <- c(0.05,0.0001,0.0001)
stepPar3f <- c(mu3f,sigma3f,zeromass3f)
angleMean3f <- c(pi,pi,0)
kappa3f <- c(1,1,1)
anglePar3f <- c(angleMean3f,kappa3f)

m3f <- fitHMM(data=ochmm,nbStates=3,stepPar0=stepPar3f,
              anglePar0=anglePar3f)

############################################################

AIC(m2a, m2b, m2c, m2d, m2e, m2f)

AIC(m2a, m2b, m2c, m2d, m2e, m2f, m3a, m3b, m3c, m3d, m3e, m3f)
AIC(m2a, m3a)
model.sel(m2a, m3a)


##Residuals plot##
pr <- pseudoRes(m2a)
# plot pseudo-residuals
plotPR(m2a)


pr <- pseudoRes(m3a)
# plot pseudo-residuals
plotPR(m3a)

states <- viterbi(m3a)

#assign states to GPS locs
ochmm_States<-ochmm
ochmm_States$State<-states


summary(states)
ochmm_States<-as.data.frame(ochmm_States)
mapview(ochmm_States, xcol="x", ycol="y", zcol="State")

ochmm_States$Used <- "1"



##############################################
############## Creating MCPs/Random Points #################
####Ocelots####
ocelotsmcp <- ocelots
ocelotsmcp <-data.frame(splitData(ocelotsmcp, ocelotsmcp$ID))
ocelots.sp <-SpatialPointsDataFrame(ocelotsmcp[,c('Longitude','Latitude')],
                                    proj4string=CRS('+proj=longlat'),
                                    data=ocelotsmcp)

ocelots.proj<-spTransform(ocelots.sp, CRS("+proj=utm +zone=14"))


#E10F
e10f.proj <- ocelots.proj[ocelots.proj$ID == "E10F",]
e10fmcp95 <- mcp(e10f.proj, percent = 95, unin = "m", unout = "km")
crs(e10fmcp95) <- CRS("+proj=utm +zone=14")
plot(e10fmcp95)
writeOGR(e10fmcp95, "C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots", "e10fmcp95", "ESRI Shapefile")

id_e10f <- e10f.proj$ID[1]

k_e10f <- nrow(e10f.proj)
e10fmcp<- readOGR("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots/e10fmcp95.shp")

randoms_e10f <- spsample(e10fmcp, 5*k_e10f, type = "random")
randoms_e10f <- as.data.frame(randoms_e10f)
head(randoms_e10f)
randoms_e10f$ID <- id_e10f
randoms_e10f$x <- randoms_e10f$x/1000
randoms_e10f$y <- randoms_e10f$y/1000
randoms_e10f$State <- floor(runif(5*k_e10f, min = 1, max = 4))
randoms_e10f$Used <- "0"


#E15M
e15m.proj <- ocelots.proj[ocelots.proj$ID == "E15M",]
e15mmcp95 <- mcp(e15m.proj, percent = 95, unin = "m", unout = "km")
crs(e15mmcp95) <- CRS("+proj=utm +zone=14")
plot(e15mmcp95)
writeOGR(e15mmcp95, "C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots", "e15mmcp95", "ESRI Shapefile")

id_e15m <- e15m.proj$ID[1]

k_e15m <- nrow(e15m.proj)
e15mmcp<- readOGR("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots/e15mmcp95.shp")

randoms_e15m <- spsample(e15mmcp, 5*k_e15m, type = "random")
randoms_e15m <- as.data.frame(randoms_e15m)
head(randoms_e15m)
randoms_e15m$ID <- id_e15m
randoms_e15m$x <- randoms_e15m$x/1000
randoms_e15m$y <- randoms_e15m$y/1000
randoms_e15m$State <- floor(runif(5*k_e15m, min = 1, max = 4))
randoms_e15m$Used <- "0"


#E17F
e17f.proj <- ocelots.proj[ocelots.proj$ID == "E17F",]
e17fmcp95 <- mcp(e17f.proj, percent = 95, unin = "m", unout = "km")
crs(e17fmcp95) <- CRS("+proj=utm +zone=14")
plot(e17fmcp95)
writeOGR(e17fmcp95, "C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots", "e17fmcp95", "ESRI Shapefile")

id_e17f <- e17f.proj$ID[1]

k_e17f <- nrow(e17f.proj)
e17fmcp<- readOGR("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots/e17fmcp95.shp")

randoms_e17f <- spsample(e17fmcp, 5*k_e17f, type = "random")
randoms_e17f <- as.data.frame(randoms_e17f)
head(randoms_e17f)
randoms_e17f$ID <- id_e17f
randoms_e17f$x <- randoms_e17f$x/1000
randoms_e17f$y <- randoms_e17f$y/1000
randoms_e17f$State <- floor(runif(5*k_e17f, min = 1, max = 4))
randoms_e17f$Used <- "0"


#E18F
e18f.proj <- ocelots.proj[ocelots.proj$ID == "E18F",]
e18fmcp95 <- mcp(e18f.proj, percent = 95, unin = "m", unout = "km")
crs(e18fmcp95) <- CRS("+proj=utm +zone=14")
plot(e18fmcp95)
writeOGR(e18fmcp95, "C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots", "e18fmcp95", "ESRI Shapefile")

id_e18f <- e18f.proj$ID[1]

k_e18f <- nrow(e18f.proj)
e18fmcp<- readOGR("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots/e18fmcp95.shp")

randoms_e18f <- spsample(e18fmcp, 5*k_e18f, type = "random")
randoms_e18f <- as.data.frame(randoms_e18f)
head(randoms_e18f)
randoms_e18f$ID <- id_e18f
randoms_e18f$x <- randoms_e18f$x/1000
randoms_e18f$y <- randoms_e18f$y/1000
randoms_e18f$State <- floor(runif(5*k_e18f, min = 1, max = 4))
randoms_e18f$Used <- "0"

write.csv(randoms_e18f, file = "E18FRandomPoints.csv")

#E19M
e19m.proj <- ocelots.proj[ocelots.proj$ID == "E19M",]
e19mmcp95 <- mcp(e19m.proj, percent = 95, unin = "m", unout = "km")
crs(e19mmcp95) <- CRS("+proj=utm +zone=14")
plot(e19mmcp95)
writeOGR(e19mmcp95, "C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots", "e19mmcp95", "ESRI Shapefile")

id_e19m <- e19m.proj$ID[1]

k_e19m <- nrow(e19m.proj)
e19mmcp<- readOGR("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots/e19mmcp95.shp")

randoms_e19m <- spsample(e19mmcp, 5*k_e19m, type = "random")
randoms_e19m <- as.data.frame(randoms_e19m)
head(randoms_e19m)
randoms_e19m$ID <- id_e19m
randoms_e19m$x <- randoms_e19m$x/1000
randoms_e19m$y <- randoms_e19m$y/1000
randoms_e19m$State <- floor(runif(5*k_e19m, min = 1, max = 4))
randoms_e19m$Used <- "0"


#E20F
e20f.proj <- ocelots.proj[ocelots.proj$ID == "E20F",]
e20fmcp95 <- mcp(e20f.proj, percent = 95, unin = "m", unout = "km")
crs(e20fmcp95) <- CRS("+proj=utm +zone=14")
plot(e20fmcp95)
writeOGR(e20fmcp95, "C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots", "e20fmcp95", "ESRI Shapefile")

id_e20f <- e20f.proj$ID[1]

k_e20f <- nrow(e20f.proj)
e20fmcp<- readOGR("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots/e20fmcp95.shp")

randoms_e20f <- spsample(e20fmcp, 5*k_e20f, type = "random")
randoms_e20f <- as.data.frame(randoms_e20f)
head(randoms_e20f)
randoms_e20f$ID <- id_e20f
randoms_e20f$x <- randoms_e20f$x/1000
randoms_e20f$y <- randoms_e20f$y/1000
randoms_e20f$State <- floor(runif(5*k_e20f, min = 1, max = 4))
randoms_e20f$Used <- "0"



#E21F
e21f.proj <- ocelots.proj[ocelots.proj$ID == "E21F",]
e21fmcp95 <- mcp(e21f.proj, percent = 95, unin = "m", unout = "km")
crs(e21fmcp95) <- CRS("+proj=utm +zone=14")
plot(e21fmcp95)
writeOGR(e21fmcp95, "C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots", "e21fmcp95", "ESRI Shapefile")

id_e21f <- e21f.proj$ID[1]

k_e21f <- nrow(e21f.proj)
e21fmcp<- readOGR("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots/e21fmcp95.shp")

randoms_e21f <- spsample(e21fmcp, 5*k_e21f, type = "random")
randoms_e21f <- as.data.frame(randoms_e21f)
head(randoms_e21f)
randoms_e21f$ID <- id_e21f
randoms_e21f$x <- randoms_e21f$x/1000
randoms_e21f$y <- randoms_e21f$y/1000
randoms_e21f$State <- floor(runif(5*k_e21f, min = 1, max = 4))
randoms_e21f$Used <- "0"



#Y22F
y22f.proj <- ocelots.proj[ocelots.proj$ID == "Y22F",]
y22fmcp95 <- mcp(y22f.proj, percent = 95, unin = "m", unout = "km")
crs(y22fmcp95) <- CRS("+proj=utm +zone=14")
plot(y22fmcp95)
writeOGR(y22fmcp95, "C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots", "y22fmcp95", "ESRI Shapefile")

id_y22f <- y22f.proj$ID[1]

k_y22f <- nrow(y22f.proj)
y22fmcp<- readOGR("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots/y22fmcp95.shp")

randoms_y22f <- spsample(y22fmcp, 5*k_y22f, type = "random")
randoms_y22f <- as.data.frame(randoms_y22f)
head(randoms_y22f)
randoms_y22f$ID <- id_y22f
randoms_y22f$x <- randoms_y22f$x/1000
randoms_y22f$y <- randoms_y22f$y/1000
randoms_y22f$State <- floor(runif(5*k_y22f, min = 1, max = 4))
randoms_y22f$Used <- "0"



#Y27M
y27m.proj <- ocelots.proj[ocelots.proj$ID == "Y27M",]
y27mmcp95 <- mcp(y27m.proj, percent = 95, unin = "m", unout = "km")
crs(y27mmcp95) <- CRS("+proj=utm +zone=14")
plot(y27mmcp95)
writeOGR(y27mmcp95, "C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots", "y27mmcp95", "ESRI Shapefile")

id_y27m <- y27m.proj$ID[1]

k_y27m <- nrow(y27m.proj)
y27mmcp<- readOGR("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots/y27mmcp95.shp")

randoms_y27m <- spsample(y27mmcp, 5*k_y27m, type = "random")
randoms_y27m <- as.data.frame(randoms_y27m)
head(randoms_y27m)
randoms_y27m$ID <- id_y27m
randoms_y27m$x <- randoms_y27m$x/1000
randoms_y27m$y <- randoms_y27m$y/1000
randoms_y27m$State <- floor(runif(5*k_y27m, min = 1, max = 4))
randoms_y27m$Used <- "0"


#Y29F
y29f.proj <- ocelots.proj[ocelots.proj$ID == "Y29F",]
y29fmcp95 <- mcp(y29f.proj, percent = 95, unin = "m", unout = "km")
crs(y29fmcp95) <- CRS("+proj=utm +zone=14")
plot(y29fmcp95)
writeOGR(y29fmcp95, "C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots", "y29fmcp95", "ESRI Shapefile")

id_y29f <- y29f.proj$ID[1]

k_y29f <- nrow(y29f.proj)
y29fmcp<- readOGR("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots/y29fmcp95.shp")

randoms_y29f <- spsample(y29fmcp, 5*k_y29f, type = "random")
randoms_y29f <- as.data.frame(randoms_y29f)
head(randoms_y29f)
randoms_y29f$ID <- id_y29f
randoms_y29f$x <- randoms_y29f$x/1000
randoms_y29f$y <- randoms_y29f$y/1000
randoms_y29f$State <- floor(runif(5*k_y29f, min = 1, max = 4))
randoms_y29f$Used <- "0"


#Y30F
y30f.proj <- ocelots.proj[ocelots.proj$ID == "Y30F",]
y30fmcp95 <- mcp(y30f.proj, percent = 95, unin = "m", unout = "km")
crs(y30fmcp95) <- CRS("+proj=utm +zone=14")
plot(y30fmcp95)
writeOGR(y30fmcp95, "C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots", "y30fmcp95", "ESRI Shapefile")

id_y30f <- y30f.proj$ID[1]

k_y30f <- nrow(y30f.proj)
y30fmcp<- readOGR("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots/y30fmcp95.shp")

randoms_y30f <- spsample(y30fmcp, 5*k_y30f, type = "random")
randoms_y30f <- as.data.frame(randoms_y30f)
head(randoms_y30f)
randoms_y30f$ID <- id_y30f
randoms_y30f$x <- randoms_y30f$x/1000
randoms_y30f$y <- randoms_y30f$y/1000
randoms_y30f$State <- floor(runif(5*k_y30f, min = 1, max = 4))
randoms_y30f$Used <- "0"


#Y5M
y5m.proj <- ocelots.proj[ocelots.proj$ID == "Y5M",]
y5mmcp95 <- mcp(y5m.proj, percent = 95, unin = "m", unout = "km")
crs(y5mmcp95) <- CRS("+proj=utm +zone=14")
plot(y5mmcp95)
writeOGR(y5mmcp95, "C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots", "y5mmcp95", "ESRI Shapefile")

id_y5m <- y5m.proj$ID[1]

k_y5m <- nrow(y5m.proj)
y5mmcp<- readOGR("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/MCPs_RSF/Ocelots/y5mmcp95.shp")

randoms_y5m <- spsample(y5mmcp, 5*k_y5m, type = "random")
randoms_y5m <- as.data.frame(randoms_y5m)
head(randoms_y5m)
randoms_y5m$ID <- id_y5m
randoms_y5m$x <- randoms_y5m$x/1000
randoms_y5m$y <- randoms_y5m$y/1000
randoms_y5m$State <- floor(runif(5*k_y5m, min = 1, max = 4))
randoms_y5m$Used <- "0"


ocelotrandoms <- rbind(randoms_e10f, randoms_e17f, randoms_e18f, randoms_e19m, randoms_e20f,
                       randoms_y27m, randoms_y30f, randoms_y5m)
table(ocelotrandoms$ID)

write.csv(ocelotrandoms, file = "RandomPoints_Ocelots.csv")

ocelotrandoms <- read.csv("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/RandomPoints_Ocelots.csv")
table(ocelots$ID)
table(ocelotrandoms$ID)


#### With Covars/Models ####
##Ocelots##

ocelotsrsf <- read.table("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/OcelotsCovars.txt", sep = ",", header = T)
head(ocelotsrsf)
tail(ocelotsrsf)
dim(ocelotsrsf)
ocelotsrsf <- ocelotsrsf[, c(2:12)]
ocelotsrsf$State <- as.factor(ocelotsrsf$State)

colnames(ocelotsrsf)[9] <- "CanopyCover"
colnames(ocelotsrsf)[10] <- "DistHeavy"
colnames(ocelotsrsf)[11] <- "DistLow"

ocelotsrsf$DistLowScale <- scale(ocelotsrsf$DistLow)
ocelotsrsf$DistHeavyScale <- scale(ocelotsrsf$DistHeavy)

####Density Rasters####

ocelotsrsf <-SpatialPointsDataFrame(ocelotsrsf[,c('X_UTM','Y_UTM')],
                                    proj4string=CRS("+proj=utm +zone=14"),
                                    data=ocelotsrsf)


###############################
#Combine Used and Available

ocelotused <- read.csv("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/HMMData_Ocelots2.csv")
head(ocelotused)


ocelotrandom <- read.csv("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/RandomPoints_Ocelots.csv")
head(ocelotrandom)

ocelotused <- ocelotused[, c(2, 5:8, 12, 18, 19)]

ocelotsrsf <- rbind(ocelotrandom, ocelotused)
table(ocelotsrsf$Used)
head(ocelotsrsf)
ocelotsrsf$X_UTM <- ocelotsrsf$x*1000
ocelotsrsf$Y_UTM <- ocelotsrsf$y*1000
ocelotsrsf <-SpatialPointsDataFrame(ocelotsrsf[,c('X_UTM','Y_UTM')],
                                    proj4string=CRS("+proj=utm +zone=14"),
                                    data=ocelotsrsf)






Density_1m <- raster("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/10m_Rasters/Density1m_10m.tif")
Density_2m <- raster("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/10m_Rasters/Density2m_10m.tif")
Density_3m <- raster("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/10m_Rasters/Density3m_10m.tif")
Density_4m <- raster("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/10m_Rasters/Density4m_10m.tif")
Density_5m <- raster("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/10m_Rasters/Density5m_10m.tif")
Density_Over5m <- raster("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/10m_Rasters/DensityOver5m_10m.tif")

DistHeavy <- raster("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/10m_Rasters/DistHeavy_10m.tif")
DistLow <- raster("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/10m_Rasters/DistLow_10m.tif")
CanopyCover <- raster("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/10m_Rasters/CanCover_10m.tif")



ocelotsrsf$Dens1m <- terra::extract(Density_1m, ocelotsrsf)
ocelotsrsf$Dens2m <- terra::extract(Density_2m, ocelotsrsf)
ocelotsrsf$Dens3m <- terra::extract(Density_3m, ocelotsrsf)
ocelotsrsf$Dens4m <- terra::extract(Density_4m, ocelotsrsf)
ocelotsrsf$Dens5m <- terra::extract(Density_5m, ocelotsrsf)
ocelotsrsf$DensOver5 <- terra::extract(Density_Over5m, ocelotsrsf)

ocelotsrsf$DistHeavy <- terra::extract(DistHeavy, ocelotsrsf)
ocelotsrsf$DistLow <- terra::extract(DistLow, ocelotsrsf)
ocelotsrsf$CanopyCover <- terra::extract(CanopyCover, ocelotsrsf)

ocelotsrsf <- as.data.frame(ocelotsrsf)
ocelotsrsf <- ocelotsrsf[, c(1:19)]

ocelotsrsf[, 11:19][is.na(ocelotsrsf[, 11:19])] <- 0
ocelotsrsf <- ocelotsrsf[,1:19]


ocelotsrsf$DensMedium <- ocelotsrsf$Dens2m + ocelotsrsf$Dens3m
ocelotsrsf$DensHigh <- ocelotsrsf$Dens4m + ocelotsrsf$Dens5m + ocelotsrsf$DensOver5

ocelotsrsf$DensHighScale <- scale(ocelotsrsf$DensHigh, center = T)
ocelotsrsf$DensMediumScale <- scale(ocelotsrsf$DensMedium, center = T)
ocelotsrsf$Dens1mScale <- scale(ocelotsrsf$Dens1m, center = T)
ocelotsrsf$CanopyCoverScale <- scale(ocelotsrsf$CanopyCover, center = T)

ocelotsrsf$DistLowScale <- scale(ocelotsrsf$DistLow, center = T)
ocelotsrsf$DistHeavyScale <- scale(ocelotsrsf$DistHeavy, center = T)

ocelotsrsf$State <- as.factor(ocelotsrsf$State)

write.csv(ocelotsrsf, file = "HM_RSFData_Ocelots.csv")

ocelotsrsf <- read.csv("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/HM_RSFData_Ocelots.csv")

## Run RSF Model
x18 <- glmer(Used ~ (1|ID) + Dens1mScale*State + DistHeavyScale*State + CanopyCoverScale*State + DistLowScale*State, family = "binomial", data = ocelotsrsf)


summary(x18)


#Check correlation
cor(ocelotsrsf$Dens1m, ocelotsrsf$DensMedium)
cor(ocelotsrsf$DensMedium, ocelotsrsf$CanopyCover)
cor(ocelotsrsf$Dens1m, ocelotsrsf$CanopyCover)
cor(ocelotsrsf$DistHeavy, ocelotsrsf$CanopyCover)
cor(ocelotsrsf$DistLow, ocelotsrsf$CanopyCover)
cor(ocelotsrsf$DensHigh, ocelotsrsf$DensMedium)

##Predictive Plots


#State1
LogState1 <- (-3.331529 + 0.872291*ocelotsrsf$Dens1mScale -0.048751*0 + 0.974756*0 + 0.470521*ocelotsrsf$DensMediumScale
              +2.016112*ocelotsrsf$CanopyCover -0.358124*0 -0.702409 *0 -0.063880*0 -0.145897*0 -0.219852*0 -0.973426*0)
OddsState1 <- exp(LogState1)
ProbState1 <- (OddsState1/(1 + OddsState1))
ocelotsrsf$ProbState1 <- ProbState1

#State2
LogState2 <- (-3.331529 + 0.872291*ocelotsrsf$Dens1mScale -0.048751*1 + 0.974756*0 + 0.470521*ocelotsrsf$DensMediumScale
              +2.016112*ocelotsrsf$CanopyCover -0.358124*1 -0.702409 *0 -0.063880*1 -0.145897*0 -0.219852*1 -0.973426*0)
OddsState2 <- exp(LogState2)
ProbState2 <- (OddsState2/(1 + OddsState2))
ocelotsrsf$ProbState2 <- ProbState2


#State3
LogState3 <- (-3.331529 + 0.872291*ocelotsrsf$Dens1mScale -0.048751*0 + 0.974756*1 + 0.470521*ocelotsrsf$DensMediumScale
              +2.016112*ocelotsrsf$CanopyCover -0.358124*0 -0.702409 *1 -0.063880*0 -0.145897*1 -0.219852*0 -0.973426*1)
OddsState3 <- exp(LogState3)
ProbState3 <- (OddsState3/(1 + OddsState3))
ocelotsrsf$ProbState3 <- ProbState3




ocused <- ocelotsrsf[ocelotsrsf$Used == 1,]
head(ocused)

ocused <- ocused[ocused$CanopyCover < 1.00000000001,]

ggplot(ocused, aes(x=CanopyCover, y=ProbState1)) +  geom_point(color = "black") +
  stat_smooth(method="glm", method.args=list(family="binomial"), se=F, size = 1.2) +
  #stat_smooth(method="glm", aes(fill = Year),method.args=list(family="binomial"), se=T, size = 1.2) +
  labs(x = "Canopy Cover", y = "Probability of Use")+
  theme(axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold"),axis.title.x = element_text(face = "bold", size = 12),axis.title.y = element_text(face = "bold", size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



##MultiLine Plot
#Canopy Cover
ggplot(ocused, aes(x=CanopyCover, y = ProbState1, color = "Resting")) +   geom_smooth() + geom_smooth(aes(y = ProbState2, color="Hunting")) +
  geom_smooth() + geom_smooth(aes(y = ProbState3, color="Exploring")) +
  labs(x = "Canopy Cover", y = "Probability of Use") +
  labs(color="Ocelot Behavior") +
  xlim (0,1) + theme(axis.title.x = element_text( size=12), axis.title.y = element_text( size = 12) ) +
  theme(axis.text = element_text(size = 12), panel.grid.major = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white"), axis.line = element_line(color = "black")) 

#Low Cover
ggplot(ocused, aes(x=Dens1m, y = ProbState1, color = "Resting")) +   geom_smooth() + geom_smooth(aes(y = ProbState2, color="Hunting")) +
  geom_smooth() + geom_smooth(aes(y = ProbState3, color="Exploring")) +
  labs(x = "0 - 1m Cover Density", y = "Probability of Use") +
  labs(color="Ocelot Behavior") +
  xlim (0,50) + theme(axis.title.x = element_text( size=12), axis.title.y = element_text( size = 12) ) +
  theme(axis.text = element_text(size = 12), panel.grid.major = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white"), axis.line = element_line(color = "black")) 


#Medium Cover
ggplot(ocused, aes(x=DensMedium, y = ProbState1, color = "Resting")) +   geom_smooth() + geom_smooth(aes(y = ProbState2, color="Hunting")) +
  geom_smooth() + geom_smooth(aes(y = ProbState3, color="Exploring")) +
  labs(x = "1 - 3m Cover Density", y = "Probability of Use") +
  labs(color="Ocelot Behavior") +
  xlim (0,100) + theme(axis.title.x = element_text( size=12), axis.title.y = element_text( size = 12) ) +
  theme(axis.text = element_text(size = 12), panel.grid.major = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white"), axis.line = element_line(color = "black")) 





####Temporal Overlap####
##Ocelots
ocelotsrsf <- read.table("C:/Users/KUMS2113/Desktop/A&M Stuff/Research/Hidden-Markov/OcelotsCovars.txt", sep = ",", header = T)
ocelotstemp <- ocelotsrsf[ocelotsrsf$Used == "1",]
head(ocelotstemp)
table(ocelotstemp$Used)

timetest <- ocelotstemp$Time
Times2 <- format(strptime(timetest, "%H:%M:%S"), "%H:%M")
Times3 <- gettime(Times2, format = "%H:%M", scale = "hour")

ocelotstemp <- cbind(ocelotstemp, Times3)


ocelotstemp$State <- as.factor(ocelotstemp$State)

colnames(ocelotstemp)[7] <- "Ocelot_Behavior"
ocelotstemp$Ocelot_Behavior <- as.character(ocelotstemp$Ocelot_Behavior)
ocelotstemp$Ocelot_Behavior[ocelotstemp$Ocelot_Behavior == "1"] <- "Resting"
ocelotstemp$Ocelot_Behavior[ocelotstemp$Ocelot_Behavior == "2"] <- "Hunting"
ocelotstemp$Ocelot_Behavior[ocelotstemp$Ocelot_Behavior == "3"] <- "Exploring"

head(ocelotstemp)
ocelotstemp$TimeScale <- rescale(ocelotstemp$Times3, to = c(0,1))
ocelotstemp$TimeRad <- ocelotstemp$TimeScale*2*pi



####Coefficients of Overlap

oc1 <- ocelotstemp[ocelotstemp$Ocelot_Behavior == "Resting",]
oc1times <- oc1[, 30]
oc1times <- as.circular(oc1times)

oc2 <- ocelotstemp[ocelotstemp$Ocelot_Behavior == "Hunting",]
oc2times <- oc2[, 30]
oc2times <- as.circular(oc2times)

oc3 <- ocelotstemp[ocelotstemp$Ocelot_Behavior == "Exploring",]
oc3times <- oc3[, 30]
oc3times <- as.circular(oc3times)

overlapEst(oc1$TimeRad, oc2$TimeRad)
overlapEst(oc1$TimeRad, oc3$TimeRad)
overlapEst(oc2$TimeRad, oc3$TimeRad)





##WATSON U TEST##
watson.two.test(oc1times, oc2times)
watson.two.test(oc1times, oc3times)
watson.two.test(oc2times, oc3times)

watson.two.test(bob1times, bob2times)
watson.two.test(bob1times, bob3times)
watson.two.test(bob2times, bob3times)

watson.two.test(coy1times, coy2times)
watson.two.test(coy1times, coy3times)
watson.two.test(coy2times, coy3times)

###
watson.two.test(oc1times, bob1times)
watson.two.test(oc2times, bob2times)
watson.two.test(oc3times, bob3times)


watson.two.test(oc1times, coy1times)
watson.two.test(oc2times, coy2times)
watson.two.test(oc3times, coy3times)

watson.two.test(coy1times, bob1times)
watson.two.test(coy2times, bob2times)
watson.two.test(coy3times, bob3times)



ocelotstemp$TimeRad <- as.circular(ocelotstemp$TimeRad)

#Make Transparent Colors
mycol2 <- rgb(255, 0, 0, max = 255, alpha = 100)
mycol3 <- rgb(255, 255, 0, max = 255, alpha = 100)
mycol4 <- rgb(0, 0, 255, max = 255, alpha = 100)

rose.diag(oc3times, col = mycol4, bins = 12, prop = 2.1, units = "hours") #FIRST
rose.diag(oc1times, col = mycol3, bins = 12,  add = TRUE, prop = 2.1, units = "hours") #THIRD
rose.diag(oc2times, col = mycol2, bins = 12, add = TRUE, prop = 2.1, units = "hours") #SECOND



hist(ocelotstemp$TimeRad)




######################################################################
################# PREDICTIVE PLOTS  #####################


####Ocelots####


####DistHeavy####
beta = matrix(fixef(x18))

myCanopyCoverScale     = as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.50))
myDistLowScale = as.numeric(quantile(ocelotsrsf$DistLowScale,0.50))
myDens1mScale = as.numeric(quantile(ocelotsrsf$Dens1mScale,0.50))



# For a plot:

# for State 1
MyInt             = c(rep(1,100))
MyCanopyCoverScale     = c(rep(myCanopyCoverScale,100))
MyState2          = c(rep(0,100))
MyState3          = c(rep(0,100))
MyDistLowScale = c(rep(myDistLowScale,100))
MyDens1mScale = c(rep(myDens1mScale,100))
MyDistHeavyScale     = seq(from = min(ocelotsrsf$DistHeavyScale), 
                           to = max(ocelotsrsf$DistHeavyScale), length.out = 100)

MyCanopyCoverScaleState2       = c(rep((myCanopyCoverScale*0),100))
MyCanopyCoverScaleState3       = c(rep((myCanopyCoverScale*0),100))

MyDistLowScaleState2   = c(rep((myDistLowScale*0),100))
MyDistLowScaleState3   = c(rep((myDistLowScale*0),100))

MyDens1mScaleState2   = c(rep((myDens1mScale*0),100))
MyDens1mScaleState3   = c(rep((myDens1mScale*0),100))


MyDistHeavyScaleState2       = c((MyDistHeavyScale*0))
MyDistHeavyScaleState3       = c((MyDistHeavyScale*0))

MyX1 = cbind(
  MyInt,MyDens1mScale,MyState2,MyState3, MyDistHeavyScale,MyCanopyCoverScale, 
  MyDistLowScale,MyDens1mScaleState2, MyDens1mScaleState3, MyDistHeavyScaleState2,
  MyDistHeavyScaleState3,MyCanopyCoverScaleState2,MyCanopyCoverScaleState3,MyDistLowScaleState2,MyDistLowScaleState3 
  
)

yhat1 = MyX1%*%beta
phat1 = 1/(1 + exp(-1*yhat1))

plot(MyDistHeavyScale,phat1,type = 'l',
     mgp=c(2.5,1,0),
     cex.lab = 1.3,
     ylab = "Probability of use",
     xlab = "Distance to Heavy Cover", xlim = c(-0.5, 5.1), ylim = c(0, 0.15)
)



text(4.9,0.029,"State 1\n(resting)")

# For state 2:
MyInt             = c(rep(1,100))
MyCanopyCoverScale     = c(rep(myCanopyCoverScale,100))
MyState2          = c(rep(1,100))
MyState3          = c(rep(0,100))
MyDistLowScale = c(rep(myDistLowScale,100))
MyDens1mScale = c(rep(myDens1mScale,100))

MyDistHeavyScale     = seq(from = min(ocelotsrsf$DistHeavyScale), 
                           to = max(ocelotsrsf$DistHeavyScale), length.out = 100)

MyCanopyCoverScaleState2       = c(rep((myCanopyCoverScale*1),100))
MyCanopyCoverScaleState3       = c(rep((myCanopyCoverScale*0),100))

MyDistLowScaleState2   = c(rep((myDistLowScale*1),100))
MyDistLowScaleState3   = c(rep((myDistLowScale*0),100))

MyDens1mScaleState2   = c(rep((myDens1mScale*1),100))
MyDens1mScaleState3   = c(rep((myDens1mScale*0),100))

MyDistHeavyScaleState2       = c((MyDistHeavyScale*1))
MyDistHeavyScaleState3       = c((MyDistHeavyScale*0))

MyX2 = cbind(
  MyInt,MyDens1mScale,MyState2,MyState3, MyDistHeavyScale,MyCanopyCoverScale, 
  MyDistLowScale,MyDens1mScaleState2, MyDens1mScaleState3, MyDistHeavyScaleState2,
  MyDistHeavyScaleState3,MyCanopyCoverScaleState2,MyCanopyCoverScaleState3,MyDistLowScaleState2,MyDistLowScaleState3
  
)

yhat2 = MyX2%*%beta
phat2 = 1/(1 + exp(-1*yhat2))
points(MyDistHeavyScale,phat2,type = 'l',lty=2)
text(4.9,0.065,"State 2\n(hunting)")


# For state 3:
MyInt             = c(rep(1,100))
MyCanopyCoverScale     = c(rep(myCanopyCoverScale,100))
MyState2          = c(rep(0,100))
MyState3          = c(rep(1,100))
MyDistLowScale = c(rep(myDistLowScale,100))
MyDens1mScale = c(rep(myDens1mScale,100))

MyDistHeavyScale     = seq(from = min(ocelotsrsf$DistHeavyScale), 
                           to = max(ocelotsrsf$DistHeavyScale), length.out = 100)

MyCanopyCoverScaleState2       = c(rep((myCanopyCoverScale*0),100))
MyCanopyCoverScaleState3       = c(rep((myCanopyCoverScale*1),100))

MyDistLowScaleState2   = c(rep((myDistLowScale*0),100))
MyDistLowScaleState3   = c(rep((myDistLowScale*1),100))

MyDens1mScaleState2   = c(rep((myDens1mScale*0),100))
MyDens1mScaleState3   = c(rep((myDens1mScale*1),100))


MyDistHeavyScaleState2       = c((MyDistHeavyScale*0))
MyDistHeavyScaleState3       = c((MyDistHeavyScale*1))

MyX3 = cbind(
  MyInt,MyDens1mScale,MyState2,MyState3, MyDistHeavyScale,MyCanopyCoverScale, 
  MyDistLowScale,MyDens1mScaleState2, MyDens1mScaleState3, MyDistHeavyScaleState2,
  MyDistHeavyScaleState3,MyCanopyCoverScaleState2,MyCanopyCoverScaleState3,MyDistLowScaleState2,MyDistLowScaleState3
  
)

yhat3 = MyX3%*%beta
phat3 = 1/(1 + exp(-1*yhat3))
points(MyDistHeavyScale,phat3,type = 'l',lty=3)
text(4.9,0.007,"State 3\n(exploring)")

arrows(as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.50)),0.117,
       as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.50)),-0.015,length = 0.1, lty = 5, col = "dark grey")
text(as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.50)),0.12,
     adj=c(0,0),srt=90,cex = 0.7,
     bquote("Median :" ==
              .(round(as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.50)),3))
     ))

arrows(as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.25)),0.117,
       as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.25)),-0.015,length = 0.1, lty = 5, col = "dark grey")
text(as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.25)),0.12,
     adj=c(0,0),srt=90,cex = 0.7,
     bquote("Lower quartile:" ==
              .(round(as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.25)),3))
     ))

arrows(as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.75)),0.117,
       as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.75)),-0.015, length = 0.1,lty = 5, col = "dark grey")
text(as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.75)),0.12,
     adj=c(0,0),srt=90,cex = 0.7,
     bquote("Upper quartile:" ==
              .(round(as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.75)),3))
     ))


####Canopy Cover ####

#Wester Plot Code - Extended for 4 Variables#
beta = matrix(fixef(x18))

myDistHeavyScale     = as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.50))
myDens1mScale = as.numeric(quantile(ocelotsrsf$Dens1mScale,0.50))
myDistLowScale = as.numeric(quantile(ocelotsrsf$DistLowScale,0.50))



# For a plot:

# for State 1
MyInt             = c(rep(1,100))
MyDistHeavyScale     = c(rep(myDistHeavyScale,100))
MyState2          = c(rep(0,100))
MyState3          = c(rep(0,100))
MyDens1mScale = c(rep(myDens1mScale,100))
MyDistLowScale = c(rep(myDistLowScale,100))
MyCanopyCoverScale     = seq(from = min(ocelotsrsf$CanopyCoverScale), 
                             to = max(ocelotsrsf$CanopyCoverScale), length.out = 100)

MyDistHeavyScaleState2       = c(rep((myDistHeavyScale*0),100))
MyDistHeavyScaleState3       = c(rep((myDistHeavyScale*0),100))

MyDens1mScaleState2   = c(rep((myDens1mScale*0),100))
MyDens1mScaleState3   = c(rep((myDens1mScale*0),100))

MyDistLowScaleState2   = c(rep((myDistLowScale*0),100))
MyDistLowScaleState3   = c(rep((myDistLowScale*0),100))


MyCanopyCoverScaleState2       = c((MyCanopyCoverScale*0))
MyCanopyCoverScaleState3       = c((MyCanopyCoverScale*0))

MyX1 = cbind(
  MyInt,MyDens1mScale,MyState2,MyState3, MyDistHeavyScale,MyCanopyCoverScale,MyDistLowScale,
  MyDens1mScaleState2,MyDens1mScaleState3,MyDistHeavyScaleState2,MyDistHeavyScaleState3,
  MyCanopyCoverScaleState2,MyCanopyCoverScaleState3,MyDistLowScaleState2, MyDistLowScaleState3
)

yhat1 = MyX1%*%beta
phat1 = 1/(1 + exp(-1*yhat1))

plot(MyCanopyCoverScale,phat1,type = 'l',
     mgp=c(2.5,1,0),
     cex.lab = 1.3,
     ylab = "Probability of use",
     xlab = "Canopy Cover"
)


text(1.45,0.34,"State 1\n(resting)")

# For state 2:
MyInt             = c(rep(1,100))
MyDistHeavyScale     = c(rep(myDistHeavyScale,100))
MyState2          = c(rep(1,100))
MyState3          = c(rep(0,100))
MyDens1mScale = c(rep(myDens1mScale,100))
MyDistLowScale = c(rep(myDistLowScale,100))

MyCanopyCoverScale     = seq(from = min(ocelotsrsf$CanopyCoverScale), 
                             to = max(ocelotsrsf$CanopyCoverScale), length.out = 100)

MyDistHeavyScaleState2       = c(rep((myDistHeavyScale*1),100))
MyDistHeavyScaleState3       = c(rep((myDistHeavyScale*0),100))

MyDens1mScaleState2   = c(rep((myDens1mScale*1),100))
MyDens1mScaleState3   = c(rep((myDens1mScale*0),100))

MyDistLowScaleState2   = c(rep((myDistLowScale*1),100))
MyDistLowScaleState3   = c(rep((myDistLowScale*0),100))

MyCanopyCoverScaleState2       = c((MyCanopyCoverScale*1))
MyCanopyCoverScaleState3       = c((MyCanopyCoverScale*0))

MyX2 = cbind(
  MyInt,MyDens1mScale,MyState2,MyState3, MyDistHeavyScale,MyCanopyCoverScale,MyDistLowScale,
  MyDens1mScaleState2,MyDens1mScaleState3,MyDistHeavyScaleState2,MyDistHeavyScaleState3,
  MyCanopyCoverScaleState2,MyCanopyCoverScaleState3,MyDistLowScaleState2, MyDistLowScaleState3
)

yhat2 = MyX2%*%beta
phat2 = 1/(1 + exp(-1*yhat2))
points(MyCanopyCoverScale,phat2,type = 'l',lty=2)
text(1.42,0.24,"State 2\n(hunting)")


# For state 3:
MyInt             = c(rep(1,100))
MyDistHeavyScale     = c(rep(myDistHeavyScale,100))
MyState2          = c(rep(0,100))
MyState3          = c(rep(1,100))
MyDens1mScale = c(rep(myDens1mScale,100))
MyDistLowScale = c(rep(myDistLowScale,100))

MyCanopyCoverScale     = seq(from = min(ocelotsrsf$CanopyCoverScale), 
                             to = max(ocelotsrsf$CanopyCoverScale), length.out = 100)

MyDistHeavyScaleState2       = c(rep((myDistHeavyScale*0),100))
MyDistHeavyScaleState3       = c(rep((myDistHeavyScale*1),100))

MyDens1mScaleState2   = c(rep((myDens1mScale*0),100))
MyDens1mScaleState3   = c(rep((myDens1mScale*1),100))

MyDistLowScaleState2   = c(rep((myDistLowScale*0),100))
MyDistLowScaleState3   = c(rep((myDistLowScale*1),100))


MyCanopyCoverScaleState2       = c((MyCanopyCoverScale*0))
MyCanopyCoverScaleState3       = c((MyCanopyCoverScale*1))

MyX3 = cbind(
  MyInt,MyDens1mScale,MyState2,MyState3, MyDistHeavyScale,MyCanopyCoverScale,MyDistLowScale,
  MyDens1mScaleState2,MyDens1mScaleState3,MyDistHeavyScaleState2,MyDistHeavyScaleState3,
  MyCanopyCoverScaleState2,MyCanopyCoverScaleState3,MyDistLowScaleState2, MyDistLowScaleState3
)

yhat3 = MyX3%*%beta
phat3 = 1/(1 + exp(-1*yhat3))
points(MyCanopyCoverScale,phat3,type = 'l',lty=5)
text(1.42,0.17,"State 3\n(exploring)")

arrows(as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.50)),0.31,
       as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.50)),-0.2,length = 0.1, lty = 5, col = "dark grey")
text(as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.50)),0.32,
     adj=c(0,0),srt=90,cex = 0.7,
     bquote("Medium:" ==
              .(round(as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.50)),3))
     ))

arrows(as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.25)),0.31,
       as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.25)),-0.2,length = 0.1, lty = 5,col = "dark grey")
text(as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.25)),0.32,
     adj=c(0,0),srt=90,cex = 0.7,
     bquote("Lower quartile:" ==
              .(round(as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.25)),3))
     ))

arrows(as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.75)),0.31,
       as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.75)),-0.2,length = 0.1, lty = 5,col = "dark grey")
text(as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.75)),0.32,
     adj=c(0,0),srt=90,cex = 0.7,
     bquote("Upper quartile:" ==
              .(round(as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.75)),3))
     ))








####DistLow ####

#Wester Plot Code - Extended for 4 Variables#
beta = matrix(fixef(x18))

myDistHeavyScale     = as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.50))
myDens1mScale = as.numeric(quantile(ocelotsrsf$Dens1mScale,0.50))
myCanopyCoverScale = as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.50))



# For a plot:

# for State 1
MyInt             = c(rep(1,100))
MyDistHeavyScale     = c(rep(myDistHeavyScale,100))
MyState2          = c(rep(0,100))
MyState3          = c(rep(0,100))
MyDens1mScale = c(rep(myDens1mScale,100))
MyCanopyCoverScale = c(rep(myCanopyCoverScale,100))
MyDistLowScale     = seq(from = min(ocelotsrsf$DistLowScale), 
                         to = max(ocelotsrsf$DistLowScale), length.out = 100)

MyDistHeavyScaleState2       = c(rep((myDistHeavyScale*0),100))
MyDistHeavyScaleState3       = c(rep((myDistHeavyScale*0),100))

MyDens1mScaleState2   = c(rep((myDens1mScale*0),100))
MyDens1mScaleState3   = c(rep((myDens1mScale*0),100))

MyCanopyCoverScaleState2   = c(rep((myCanopyCoverScale*0),100))
MyCanopyCoverScaleState3   = c(rep((myCanopyCoverScale*0),100))


MyDistLowScaleState2       = c((MyDistLowScale*0))
MyDistLowScaleState3       = c((MyDistLowScale*0))

MyX1 = cbind(
  MyInt,MyDens1mScale,MyState2,MyState3, MyDistHeavyScale,MyCanopyCoverScale,MyDistLowScale,
  MyDens1mScaleState2,MyDens1mScaleState3,MyDistHeavyScaleState2,MyDistHeavyScaleState3,
  MyCanopyCoverScaleState2,MyCanopyCoverScaleState3,MyDistLowScaleState2, MyDistLowScaleState3
  
)

yhat1 = MyX1%*%beta
phat1 = 1/(1 + exp(-1*yhat1))

plot(MyDistLowScale,phat1,type = 'l',
     mgp=c(2.5,1,0),
     cex.lab = 1.3,
     ylab = "Probability of use",
     xlab = "Distance to Open Areas", ylim = c(0, 0.45)
)



text(8,0.33,"State 1\n(resting)")

# For state 2:
MyInt             = c(rep(1,100))
MyDistHeavyScale     = c(rep(myDistHeavyScale,100))
MyState2          = c(rep(1,100))
MyState3          = c(rep(0,100))
MyDens1mScale = c(rep(myDens1mScale,100))
MyCanopyCoverScale = c(rep(myCanopyCoverScale,100))

MyDistLowScale     = seq(from = min(ocelotsrsf$DistLowScale), 
                         to = max(ocelotsrsf$DistLowScale), length.out = 100)

MyDistHeavyScaleState2       = c(rep((myDistHeavyScale*1),100))
MyDistHeavyScaleState3       = c(rep((myDistHeavyScale*0),100))

MyDens1mScaleState2   = c(rep((myDens1mScale*1),100))
MyDens1mScaleState3   = c(rep((myDens1mScale*0),100))

MyCanopyCoverScaleState2   = c(rep((myCanopyCoverScale*1),100))
MyCanopyCoverScaleState3   = c(rep((myCanopyCoverScale*0),100))

MyDistLowScaleState2       = c((MyDistLowScale*1))
MyDistLowScaleState3       = c((MyDistLowScale*0))

MyX2 = cbind(
  MyInt,MyDens1mScale,MyState2,MyState3, MyDistHeavyScale,MyCanopyCoverScale,MyDistLowScale,
  MyDens1mScaleState2,MyDens1mScaleState3,MyDistHeavyScaleState2,MyDistHeavyScaleState3,
  MyCanopyCoverScaleState2,MyCanopyCoverScaleState3,MyDistLowScaleState2, MyDistLowScaleState3
  
)

yhat2 = MyX2%*%beta
phat2 = 1/(1 + exp(-1*yhat2))
points(MyDistLowScale,phat2,type = 'l',lty=2)
text(7.4,0.43,"State 2\n(hunting)")


# For state 3:
MyInt             = c(rep(1,100))
MyDistHeavyScale     = c(rep(myDistHeavyScale,100))
MyState2          = c(rep(0,100))
MyState3          = c(rep(1,100))
MyDens1mScale = c(rep(myDens1mScale,100))
MyCanopyCoverScale = c(rep(myCanopyCoverScale,100))

MyDistLowScale     = seq(from = min(ocelotsrsf$DistLowScale), 
                         to = max(ocelotsrsf$DistLowScale), length.out = 100)

MyDistHeavyScaleState2       = c(rep((myDistHeavyScale*0),100))
MyDistHeavyScaleState3       = c(rep((myDistHeavyScale*1),100))

MyDens1mScaleState2   = c(rep((myDens1mScale*0),100))
MyDens1mScaleState3   = c(rep((myDens1mScale*1),100))

MyCanopyCoverScaleState2   = c(rep((myCanopyCoverScale*0),100))
MyCanopyCoverScaleState3   = c(rep((myCanopyCoverScale*1),100))


MyDistLowScaleState2       = c((MyDistLowScale*0))
MyDistLowScaleState3       = c((MyDistLowScale*1))

MyX3 = cbind(
  MyInt,MyDens1mScale,MyState2,MyState3, MyDistHeavyScale,MyCanopyCoverScale,MyDistLowScale,
  MyDens1mScaleState2,MyDens1mScaleState3,MyDistHeavyScaleState2,MyDistHeavyScaleState3,
  MyCanopyCoverScaleState2,MyCanopyCoverScaleState3,MyDistLowScaleState2, MyDistLowScaleState3
  
)

yhat3 = MyX3%*%beta
phat3 = 1/(1 + exp(-1*yhat3))
points(MyDistLowScale,phat3,type = 'l',lty=3)
text(8,0.09,"State 3\n(exploring)")

arrows(as.numeric(quantile(ocelotsrsf$DistLowScale,0.50)),0.3,
       as.numeric(quantile(ocelotsrsf$DistLowScale,0.50)),-0.2,length = 0.1, lty = 5,col = "dark grey")
text(as.numeric(quantile(ocelotsrsf$DistLowScale,0.50)),0.31,
     adj=c(0,0),srt=90,cex = 0.7,
     bquote("Median:" ==
              .(round(as.numeric(quantile(ocelotsrsf$DistLowScale,0.50)),3))
     ))

arrows(as.numeric(quantile(ocelotsrsf$DistLowScale,0.25)),0.3,
       as.numeric(quantile(ocelotsrsf$DistLowScale,0.25)),-0.2,length = 0.1, lty = 5,col = "dark grey")
text(as.numeric(quantile(ocelotsrsf$DistLowScale,0.25)),0.31,
     adj=c(0,0),srt=90,cex = 0.7,
     bquote("Lower quartile:" ==
              .(round(as.numeric(quantile(ocelotsrsf$DistLowScale,0.25)),3))
     ))

arrows(as.numeric(quantile(ocelotsrsf$DistLowScale,0.75)),0.3,
       as.numeric(quantile(ocelotsrsf$DistLowScale,0.75)),-0.2,length = 0.1, lty = 5,col = "dark grey")
text(as.numeric(quantile(ocelotsrsf$DistLowScale,0.75)),0.31,
     adj=c(0,0),srt=90,cex = 0.7,
     bquote("Upper quartile:" ==
              .(round(as.numeric(quantile(ocelotsrsf$DistLowScale,0.75)),3))
     ))





####Dens1m Veg ####

#Wester Plot Code - Extended for 4 Variables#
beta = matrix(fixef(x18))

myDistLowScale     = as.numeric(quantile(ocelotsrsf$DistLowScale,0.50))
myDistHeavyScale = as.numeric(quantile(ocelotsrsf$DistHeavyScale,0.50))
myCanopyCoverScale = as.numeric(quantile(ocelotsrsf$CanopyCoverScale,0.50))



# For a plot:

# for State 1
MyInt             = c(rep(1,100))
MyDistLowScale     = c(rep(myDistLowScale,100))
MyState2          = c(rep(0,100))
MyState3          = c(rep(0,100))
MyDistHeavyScale = c(rep(myDistHeavyScale,100))
MyCanopyCoverScale = c(rep(myCanopyCoverScale,100))
MyDens1mScale     = seq(from = min(ocelotsrsf$Dens1mScale), 
                        to = max(ocelotsrsf$Dens1mScale), length.out = 100)

MyDistLowScaleState2       = c(rep((myDistLowScale*0),100))
MyDistLowScaleState3       = c(rep((myDistLowScale*0),100))

MyDistHeavyScaleState2   = c(rep((myDistHeavyScale*0),100))
MyDistHeavyScaleState3   = c(rep((myDistHeavyScale*0),100))

MyCanopyCoverScaleState2   = c(rep((myCanopyCoverScale*0),100))
MyCanopyCoverScaleState3   = c(rep((myCanopyCoverScale*0),100))


MyDens1mScaleState2       = c((MyDens1mScale*0))
MyDens1mScaleState3       = c((MyDens1mScale*0))

MyX1 = cbind(
  MyInt,MyDens1mScale,MyState2,MyState3, MyDistHeavyScale,MyCanopyCoverScale,MyDistLowScale,
  MyDens1mScaleState2,MyDens1mScaleState3,MyDistHeavyScaleState2,MyDistHeavyScaleState3,
  MyCanopyCoverScaleState2,MyCanopyCoverScaleState3,MyDistLowScaleState2, MyDistLowScaleState3
  
)

yhat1 = MyX1%*%beta
phat1 = 1/(1 + exp(-1*yhat1))

plot(MyDens1mScale,phat1,type = 'l',
     mgp=c(2.5,1,0),
     cex.lab = 1.3,
     ylab = "Probability of use",
     xlab = "Density 0 - 1 m Veg"
)



text(10,0.96,"State 1\n(resting)")

# For state 2:
MyInt             = c(rep(1,100))
MyDistLowScale     = c(rep(myDistLowScale,100))
MyState2          = c(rep(1,100))
MyState3          = c(rep(0,100))
MyDistHeavyScale = c(rep(myDistHeavyScale,100))
MyCanopyCoverScale = c(rep(myCanopyCoverScale,100))

MyDens1mScale     = seq(from = min(ocelotsrsf$Dens1mScale), 
                        to = max(ocelotsrsf$Dens1mScale), length.out = 100)

MyDistLowScaleState2       = c(rep((myDistLowScale*1),100))
MyDistLowScaleState3       = c(rep((myDistLowScale*0),100))

MyDistHeavyScaleState2   = c(rep((myDistHeavyScale*1),100))
MyDistHeavyScaleState3   = c(rep((myDistHeavyScale*0),100))

MyCanopyCoverScaleState2   = c(rep((myCanopyCoverScale*1),100))
MyCanopyCoverScaleState3   = c(rep((myCanopyCoverScale*0),100))

MyDens1mScaleState2       = c((MyDens1mScale*1))
MyDens1mScaleState3       = c((MyDens1mScale*0))

MyX2 = cbind(
  MyInt,MyDens1mScale,MyState2,MyState3, MyDistHeavyScale,MyCanopyCoverScale,MyDistLowScale,
  MyDens1mScaleState2,MyDens1mScaleState3,MyDistHeavyScaleState2,MyDistHeavyScaleState3,
  MyCanopyCoverScaleState2,MyCanopyCoverScaleState3,MyDistLowScaleState2, MyDistLowScaleState3
  
)

yhat2 = MyX2%*%beta
phat2 = 1/(1 + exp(-1*yhat2))
points(MyDens1mScale,phat2,type = 'l',lty=2)
text(11.4,0.89,"State 2\n(hunting)")


# For state 3:
MyInt             = c(rep(1,100))
MyDistLowScale     = c(rep(myDistLowScale,100))
MyState2          = c(rep(0,100))
MyState3          = c(rep(1,100))
MyDistHeavyScale = c(rep(myDistHeavyScale,100))
MyCanopyCoverScale = c(rep(myCanopyCoverScale,100))

MyDens1mScale     = seq(from = min(ocelotsrsf$Dens1mScale), 
                        to = max(ocelotsrsf$Dens1mScale), length.out = 100)

MyDistLowScaleState2       = c(rep((myDistLowScale*0),100))
MyDistLowScaleState3       = c(rep((myDistLowScale*1),100))

MyDistHeavyScaleState2   = c(rep((myDistHeavyScale*0),100))
MyDistHeavyScaleState3   = c(rep((myDistHeavyScale*1),100))

MyCanopyCoverScaleState2   = c(rep((myCanopyCoverScale*0),100))
MyCanopyCoverScaleState3   = c(rep((myCanopyCoverScale*1),100))


MyDens1mScaleState2       = c((MyDens1mScale*0))
MyDens1mScaleState3       = c((MyDens1mScale*1))

MyX3 = cbind(
  MyInt,MyDens1mScale,MyState2,MyState3, MyDistHeavyScale,MyCanopyCoverScale,MyDistLowScale,
  MyDens1mScaleState2,MyDens1mScaleState3,MyDistHeavyScaleState2,MyDistHeavyScaleState3,
  MyCanopyCoverScaleState2,MyCanopyCoverScaleState3,MyDistLowScaleState2, MyDistLowScaleState3
  
)

yhat3 = MyX3%*%beta
phat3 = 1/(1 + exp(-1*yhat3))
points(MyDens1mScale,phat3,type = 'l',lty=3)
text(11.3,0.24,"State 3\n(exploring)")

arrows(as.numeric(quantile(ocelotsrsf$Dens1mScale,0.50)),0.7,
       as.numeric(quantile(ocelotsrsf$Dens1mScale,0.50)),-0.015,length = 0.1, lty = 5,col = "dark grey")
text(as.numeric(quantile(ocelotsrsf$Dens1mScale,0.50)),0.71,
     adj=c(0,0),srt=90,cex = 0.7,
     bquote("Median:" ==
              .(round(as.numeric(quantile(ocelotsrsf$Dens1mScale,0.50)),3))
     ))

arrows(as.numeric(quantile(ocelotsrsf$Dens1mScale,0.25)),0.7,
       as.numeric(quantile(ocelotsrsf$Dens1mScale,0.25)),-0.015,length = 0.1, lty = 5,col = "dark grey")
text(as.numeric(quantile(ocelotsrsf$Dens1mScale,0.25)),0.71,
     adj=c(0,0),srt=90,cex = 0.7,
     bquote("Lower quartile:" ==
              .(round(as.numeric(quantile(ocelotsrsf$Dens1mScale,0.25)),3))
     ))

arrows(as.numeric(quantile(ocelotsrsf$Dens1mScale,0.75)),0.7,
       as.numeric(quantile(ocelotsrsf$Dens1mScale,0.75)),-0.015,length = 0.1, lty = 5,col = "dark grey")
text(as.numeric(quantile(ocelotsrsf$Dens1mScale,0.75)),0.71,
     adj=c(0,0),srt=90,cex = 0.7,
     bquote("Upper quartile:" ==
              .(round(as.numeric(quantile(ocelotsrsf$Dens1mScale,0.75)),3))
     ))







