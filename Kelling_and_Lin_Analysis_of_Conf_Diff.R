###
### Replication Code:
### Analysis of Conflict Diffusion over Continuous Space
### Claire Kelling and YiJyun Lin
###
### Last Updated: 5/30/19
###

#Clear the workspace
rm(list = ls())

#Libraries
library(sp)
library(gstat)
library(fields)
library(classInt)
library(maps)
library(acs)
library(tigris)
library(spdep)
library(ggplot2)
library(dplyr)
library(ade4)
library(ggmap)
library(rgdal)
library(spatstat)
library(stpp)
library(xtable)
library(rgeos)
library(GiNA)
library(readxl)
library(pracma)
library(raster)
library(elevatr)
library(plot3D)

#####
## Load Data and Initial Setup
#####
# Set working directory to location of the repository
setwd("")
# Load conflict data
conf_data <- read_excel("data/ACLED (South Sudan).xlsx")

# save Google api key to get map of South Sudan for plotting
register_google(key = "YOUR_API_KEY")
South_Sudan <- ggmap(get_map(location=c(29.950669,7.9621765), zoom=6))

#Only the event types we are interested in for this analysis
conf_data <- conf_data[which(conf_data$EVENT_TYPE %in% c("Battle-No change of territory", "Battle-Non-state actor overtakes territory",
                                                         "Battle-Government regains territory", "Violence against civilians", 
                                                         "Riots/Protests")),]

#shape file for South Sudan
south_sud <- raster::getData('GADM', country='SSD', level=0)

# Adding additional columns to make into Spatial Points
conf_data$LONGITUDE1 <- conf_data$LONGITUDE
conf_data$LATITUDE1 <- conf_data$LATITUDE

# Convert the files to spatial points and polygons with the same projection
coordinates(conf_data) <- ~LONGITUDE1+LATITUDE1
proj4string(conf_data) <- proj4string(south_sud)


#####
## Preliminary Plots
#####
#First, I will convert the time variable to something we can work with in R
conf_data$date <- as.POSIXct(conf_data$EVENT_DATE, format="%Y-%m-%d",tz=Sys.timezone())

# FIGURE 1: Histogram of event frequency over time by event category
ggplot(conf_data@data, aes(x=date, fill = EVENT_TYPE)) + geom_histogram(bins=50)+
  labs(title = "Event Types over Time", x= "Date", y = "Count")+ scale_fill_brewer("Event Type", palette="Spectral")#, guide = FALSE)

## Aggregate by block group for plotting, used in legend of Figure 1
agg_dat <- plyr::count(conf_data@data, c('EVENT_TYPE'))


# FIGURE 2a: Distribution of conflict events over space, location only
point_proc <- South_Sudan  + geom_point(aes(x = LONGITUDE, y = LATITUDE), size = 1,
                                        data = conf_data@data, col = "blue", alpha =0.1) + coord_equal() +
  ggtitle("Conflict Events, South Sudan")
point_proc


#####
## Duration setup
#####
event_data <- conf_data@data

event_data <- event_data[,c("date", "LATITUDE", "LONGITUDE", "FATALITIES", "NOTES", "EVENT_TYPE", "ACTOR1", "ACTOR2")]
event_data <- unique(event_data) #take out a handful of repeated rows

### Find events that last longer than one day by location and description
dup_event <- which(duplicated(event_data[,c("LATITUDE", "LONGITUDE", "NOTES")])) 
duration_data <- event_data[,c("LATITUDE", "LONGITUDE", "NOTES")] %>%
  group_by(LATITUDE, LONGITUDE, NOTES) %>%
  summarize(duration = n())

duration_data <- left_join(event_data[-dup_event,], duration_data)
nrow(duration_data) #how many events in our final dataset
range(duration_data$duration) #the range of durations
mean(duration_data$duration) #mean duration
1-nrow(duration_data[which(duration_data$duration > 1),])/nrow(duration_data) #percent that only last one day

# FIGURE 2b: Distribution of conflict events over space, with duration as size of point
duration_data <- duration_data %>% dplyr::rename("Duration (Days)" = "duration")
dur_point_proc <- South_Sudan  + geom_point(aes(x = LONGITUDE, y = LATITUDE, size = `Duration (Days)`),
                                        data = duration_data, col = "blue", alpha =0.1) + coord_equal() +
  ggtitle("Duration of Conflict Events")
dur_point_proc


#####
## Actors
#####
actors <-unique(c(conf_data@data$ACTOR1, conf_data@data$ACTOR2))

#TABLE 1: Top actors by frequency of conflict event involvement
agg_act <- plyr::count(c(conf_data@data$ACTOR1, conf_data@data$ACTOR2))
agg_act <- head(agg_act[order(-agg_act$freq),], n =11)
xtable(agg_act)

#####
## Test for complete spatial randomness
#####
# Transform our data into ppm object
sud_owin <- as.owin(south_sud)
xyzt <- as.matrix(duration_data[,c("LONGITUDE", "LATITUDE")])
conf_ppp <- as.ppp(xyzt, sud_owin) #22 outside of the specified window

# simultaneous
#   sig level = (nrank/(1+nsim))
f_aff_sim <- envelope(conf_ppp, fun = Fest, global = TRUE, nrank = 20, nsim = 800)

# pointwise
#   sig level = 2*(nrank/(1+nsim))
f_aff_pnt <- envelope(conf_ppp, fun = Fest, global = FALSE, nrank = 20, nsim = 800)

# FIGURE 3: F function Envelope, Simultaneous and Pointwise
par(mfrow = c(1,2))
plot(f_aff_sim, main = "F function Envelope, Simultaneous", ylab = "F function")
plot(f_aff_pnt, main = "F function Envelope, Pointwise", ylab = "F function") 

#####
## Fit a Point Process Model - Inhomogenous Models with varying intensity functions
#####

## Incorporate temporal element
# Create the dataset
xyzt <- duration_data[,c("LONGITUDE", "LATITUDE",
                         "date")]
range(duration_data$date) #this data occurs over the span of 7.5 years
colnames(xyzt) <- c("x", "y", "t")
xyzt$t <- as.numeric(as.Date(xyzt$t))
min(xyzt$t) #convert to numeric

#converting to months
xyzt$month <- (xyzt$t - min(xyzt$t))/30
xyzt <- xyzt[,-c(3)] #take out t so that month is the time variable
xyzt$month <- round(xyzt$month)
xyzt <- as.data.frame(xyzt)
colnames(xyzt) <- c('x', 'y', 't')
keep <- xyzt

#Need to create the boundary of Sudan with points
sud_bound <- fortify(south_sud)
sud_bound <- sud_bound[,c("long", "lat")]


#####
## Estimate and Plot Intensity Estimate, both static and over months
#####
#vignette: https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/inst/doc/stpp.pdf?revision=61&root=stpp&pathrev=61
#page 12
sud_bound <- as.matrix(sud_bound)
keep <- as.matrix(keep)
h <- mse2d(as.points(keep[, 1:2]), sud_bound, nsmse = 30, range = 3000)
h <- h$h[which.min(h$mse)]

#make h = 0.01 smaller for more precise estimation
nx <- 1000
ny <- 1000
nt <- 91 #length(unique(keep[,3]))

Ls <- kernel2d(as.points(keep[, 1:2]), sud_bound, 0.7, nx = nx, ny = ny)
Lt <- dim(keep)[1] * density(keep[, 3], n = nt)$y
Lst <- array(0, dim = c(nx, ny, nt))
for(k in 1:nt) Lst[,,k] <- Ls$z * Lt[k] / dim(keep)[1]

#FIGURE 4a: spatial intensity function estimate
par(mfrow = c(1,1))
raster::image(Ls$x, Ls$y, Ls$z, col=brewer.pal(11,"RdBu"), xlab = "", ylab ="", xaxt = "n", yaxt = "n")#,
colkey(side = 3, clim = c(0, 420), col=brewer.pal(11,"RdBu"), clab = "Density Estimate", 
       length = 0.5, width = 0.5, add = T) #save 5x5 inches (pdf)
plot(south_sud, add = T)

#FIGURE 4b: Temporal Intensity Function Estimate
plot(1:nt, Lt, col="white", xlab = "month", ylab = "temporal trend")
lines(Lt)

#creating estimates of the spatial intensity over time
#https://cran.r-project.org/web/packages/splancs/splancs.pdf
t <- as.data.frame(as.numeric(keep[,3]))
colnames(t) <- "t"

#keep needs to be a matrix
nx <- 0.1

b3d <- kernel3d(keep[,1:2], t$t, seq(24.151930, 35.8699, nx), seq(3.480999, 12.21900, nx),
                seq(1,90,12), 1, 0.05)

t_length <- length(seq(1,90,12))

#Find points outside of South Sudan
x_point <- c(t(meshgrid(seq(24.151930, 35.8699, nx), seq(3.480999, 12.21900, nx))$X))
y_point <- c(t(meshgrid(seq(24.151930, 35.8699, nx), seq(3.480999, 12.21900, nx))$Y))
points <- as.data.frame(cbind(x_point, y_point))
coordinates(points) <- ~x_point+y_point
proj4string(points) <- proj4string(south_sud)
in_points <- over(points, south_sud)
in_points[which(!is.na(in_points[,1])),] <- 1
in_vec <- in_points[,1]
in_mat <- matrix(in_vec, nrow = length(seq(24.151930, 35.8699, nx)), ncol = length(seq(3.480999, 12.21900, nx)))


for(i in 1:t_length){
  new_df <- b3d$v[,,i]
  new_df[which(is.na(in_mat))] <- NA
  b3d$v[,,i] <- new_df
}

oldpar <- par(mfrow=c(1,8), mai = c(0.3,0,0.3,0))
fig_caption <- c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018")

all_z_val <- c(b3d$v[,,1], b3d$v[,,2], b3d$v[,,3] ,b3d$v[,,4], b3d$v[,,5], b3d$v[,,6], b3d$v[,,7], b3d$v[,,8])
max(all_z_val, na.rm = T)
min(all_z_val, na.rm = T)

#FIGURE 5a: common scale over the complete time scale
for(i in 1:t_length){
  graphics::image(seq(24.151930, 35.8699, nx), seq(3.480999, 12.21900, nx), b3d$v[,,i], zlim=c(0,max(all_z_val, na.rm = T)),
                  asp=1, xlab="", ylab="",  col=brewer.pal(11,"RdBu"), xaxt = "n", yaxt = "n", 
                  main=fig_caption[i])
  colkey(side = 1, clim = c(0, max(all_z_val, na.rm = T)), col=brewer.pal(11,"RdBu"), clab = "Density", 
         length = 0.5, width = 0.5, add = T, dist = -0.3) #save 5x5 inches (pdf)
  plot(south_sud, add = T)
} 
par(oldpar) #3x15

oldpar <- par(mfrow=c(1,8), mai = c(0.3,0,0.3,0))

#FIGURE 5b: scale determined by year
for(i in 1:t_length){
  #i <- 1
  graphics::image(seq(24.151930, 35.8699, nx), seq(3.480999, 12.21900, nx), b3d$v[,,i], 
                  asp=1, xlab="", ylab="",  col=brewer.pal(11,"RdBu"), xaxt = "n", yaxt = "n", 
                  main=fig_caption[i])
  colkey(side = 1, clim = c(0, max(b3d$v[,,i], na.rm = T)), col=brewer.pal(11,"RdBu"), clab = "Density", 
         length = 0.5, width = 0.5, add = T, dist = -0.3) #save 5x5 inches (pdf)
  #plot(south_sud, add = T)
} 


#######################
## Covariance Function Modeling
#######################

#####
## BY YEAR
#####
### PPP over 8 years
xyzt_year1 <- as.matrix(duration_data[which(grepl("2011", duration_data$date)),
                                      c("LONGITUDE", "LATITUDE")])
xyzt_year2 <- as.matrix(duration_data[which(grepl("2012", duration_data$date)),
                                      c("LONGITUDE", "LATITUDE")])
xyzt_year3 <- as.matrix(duration_data[which(grepl("2013", duration_data$date)),
                                      c("LONGITUDE", "LATITUDE")])
xyzt_year4 <- as.matrix(duration_data[which(grepl("2014", duration_data$date)),
                                      c("LONGITUDE", "LATITUDE")])
xyzt_year5 <- as.matrix(duration_data[which(grepl("2015", duration_data$date)),
                                      c("LONGITUDE", "LATITUDE")])
xyzt_year6 <- as.matrix(duration_data[which(grepl("2016", duration_data$date)),
                                      c("LONGITUDE", "LATITUDE")])
xyzt_year7 <- as.matrix(duration_data[which(grepl("2017", duration_data$date)),
                                      c("LONGITUDE", "LATITUDE")])
xyzt_year8 <- as.matrix(duration_data[which(grepl("2018", duration_data$date)),
                                      c("LONGITUDE", "LATITUDE")])
year1_ppp <- as.ppp(xyzt_year1, sud_owin) 
year2_ppp <- as.ppp(xyzt_year2, sud_owin) 
year3_ppp <- as.ppp(xyzt_year3, sud_owin) 
year4_ppp <- as.ppp(xyzt_year4, sud_owin) 
year5_ppp <- as.ppp(xyzt_year5, sud_owin) 
year6_ppp <- as.ppp(xyzt_year6, sud_owin) 
year7_ppp <- as.ppp(xyzt_year7, sud_owin) 
year8_ppp <- as.ppp(xyzt_year8, sud_owin) 


lgcp_est1 <- lgcp.estpcf(year1_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est2 <- lgcp.estpcf(year2_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est3 <- lgcp.estpcf(year3_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est4 <- lgcp.estpcf(year4_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est5 <- lgcp.estpcf(year5_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est6 <- lgcp.estpcf(year6_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est7 <- lgcp.estpcf(year7_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est8 <- lgcp.estpcf(year8_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))

#https://rdrr.io/cran/spatstat/man/lgcp.estpcf.html
cov_f <- function(r, lgcp_est, year){
  cov_function <- as.numeric(lgcp_est$par[1]) * exp(-r/as.numeric(lgcp_est$par[2]))
  dat <- as.data.frame(cbind(r, cov_function))
  dat$year <- year
  return(dat)
}

r <- seq(0,0.45, 0.005)
full_dat <- rbind(cov_f(r,lgcp_est1, "2011"),
                  cov_f(r,lgcp_est2, "2012"),
                  cov_f(r,lgcp_est3, "2013"),
                  cov_f(r,lgcp_est4, "2014"),
                  cov_f(r,lgcp_est5, "2015"),
                  cov_f(r,lgcp_est6, "2016"),
                  cov_f(r,lgcp_est7, "2017"),
                  cov_f(r,lgcp_est8, "2018"))#,

full_dat$Year <- full_dat$year

# FIGURE 6: Covariance function by year
ggplot(data=full_dat, aes(x = r, y = cov_function, col = Year)) +geom_line(size =1.0001) + #, aes(linetype= Year)
  labs(title = "Covariance Function Estimate", x = "r", y = "C(r)")+ 
  geom_hline(yintercept=0.1, linetype = "dashed", size = 1)+ scale_color_brewer(palette="Accent")#+ theme_bw() #4x7

eff_range_dat <- NULL
for(i in unique(full_dat$year)){
  #i <- 2011
  eff_range_ind <- min(which(full_dat$cov_function[which(full_dat$year == i)] < 0.1))
  eff_range <- full_dat$r[which(full_dat$year == i)][eff_range_ind]
  eff_range_dat <- rbind(eff_range_dat, c(i, eff_range))
}

# TABLE 3: Effective Range by Year
eff_range_dat <- as.data.frame(eff_range_dat)
colnames(eff_range_dat) <- c("Year", "Effective Range")
xtable(eff_range_dat)

#####
## By Actor Type
#####

#read in the data that classifies actors
type_act <- read.csv(file = "data/actors(re-classify).csv")

#there are some events that only have one actor involved
length(which(is.na(duration_data$ACTOR2)))

# Cleaning actor type data
type_act <- type_act[,c(2,3)]
colnames(type_act) <- c("actor", "type")
type_act$type[which(type_act$actor == "Militia (Amanya)")] <- 2
type_act$type[which(type_act$actor == "Militia (David Kuburin)")] <- 2

#state, non-state, combination
actor1 <- as.data.frame(duration_data$ACTOR1)
colnames(actor1) <- c("actor1")
actor1 <- left_join(actor1, type_act, by = c(actor1 = "actor"))
colnames(actor1) <- c("actor1", "type_act1")

actor2 <- as.data.frame(duration_data$ACTOR2)
colnames(actor2) <- c("actor2")
actor2 <- left_join(actor2, type_act, by = c(actor2 = "actor"))
colnames(actor2) <- c("actor2", "type_act2")

duration_data <- cbind(duration_data, actor1$type_act1, actor2$type_act2)
colnames(duration_data)[(ncol(duration_data)-1):ncol(duration_data)] <- c("type_act1", "type_act2")
duration_data$type_act1 <- as.character(duration_data$type_act1)
duration_data$type_act2 <- as.character(duration_data$type_act2)

#1 = state actor
#2 = non-state actor
duration_data$`Actors Involved` <- NA
for(i in 1:nrow(duration_data)){
  print(i)
  if(!is.na(duration_data$type_act1[i]) & !is.na(duration_data$type_act2[i])){
    if(duration_data$type_act1[i] == "1" & duration_data$type_act2[i] == "1"){
      duration_data$`Actors Involved`[i] <- "Only State"
    }else if(duration_data$type_act1[i] == "1" & duration_data$type_act2[i] == "2"){
      duration_data$`Actors Involved`[i] <- "State and Non-State"
    }else if(duration_data$type_act1[i] == "2" & duration_data$type_act2[i] == "1"){
      duration_data$`Actors Involved`[i] <- "State and Non-State"
    }else if(duration_data$type_act1[i] == "2" & duration_data$type_act2[i] == "2"){
      duration_data$`Actors Involved`[i] <- "Only Non-State"
    }
  }else if(is.na(duration_data$type_act1[i]) & !is.na(duration_data$type_act2[i])){
    if(duration_data$type_act2[i] == "1" & is.na(duration_data$type_act1[i])){
      duration_data$`Actors Involved`[i] <- "Only State"
    }else if(duration_data$type_act2[i] == "2" & is.na(duration_data$type_act1[i])){
      duration_data$`Actors Involved`[i] <- "Only Non-State"
    }
  }else if(!is.na(duration_data$type_act1[i]) & is.na(duration_data$type_act2[i])){
    if(duration_data$type_act1[i] == "2" & is.na(duration_data$type_act2[i])){
      duration_data$`Actors Involved`[i] <- "Only Non-State"
    }else if(duration_data$type_act1[i] == "1" & is.na(duration_data$type_act2[i])){
      duration_data$`Actors Involved`[i] <- "Only State"
    }
  }
}

length(unique(duration_data$`Actors Involved`))

#reclassify when civilians involved
unique(duration_data$ACTOR1[which(is.na(duration_data$type_act1))]) #civilians, rioters, protesters
unique(duration_data$ACTOR2[which(is.na(duration_data$type_act2))]) #civilians, protesters, rioters, and prisoners

#need to make when no second actor present, doesn't think this is against civilians
duration_data$type_act2[which(is.na(duration_data$ACTOR2))] <- "None"

duration_data$`Actors Involved`[which(is.na(duration_data$type_act1) & is.na(duration_data$type_act2))] <- "Only Civilian"
duration_data$`Actors Involved`[which(is.na(duration_data$type_act1) & duration_data$type_act2 == "1")] <- "State and Civilian"
duration_data$`Actors Involved`[which(is.na(duration_data$type_act2) & duration_data$type_act1 == "1")] <- "State and Civilian"
duration_data$`Actors Involved`[which(is.na(duration_data$type_act1) & duration_data$type_act2 == "2")] <- "Non-State and Civilian"
duration_data$`Actors Involved`[which(is.na(duration_data$type_act2) & duration_data$type_act1 == "2")] <- "Non-State and Civilian"

duration_data$`Actors Involved`[which(is.na(duration_data$`Actors Involved`))] <- "Only Civilian"

#percent with only one actor involved
length(which(duration_data$type_act2 == "None"))/nrow(duration_data)
length(which(is.na(duration_data$ACTOR2)))/nrow(duration_data)

#TABLE 2: Counts of conflict events by actor dyads
agg <- duration_data %>% count(`Actors Involved`)
xtable(agg)

#finding effective range for each
xyzt_state <- as.matrix(duration_data[which(duration_data$`Actors Involved` == "Only State"),
                                      c("LONGITUDE", "LATITUDE")])
xyzt_nonstate <- as.matrix(duration_data[which(duration_data$`Actors Involved` == "Only Non-State"),
                                         c("LONGITUDE", "LATITUDE")])
xyzt_civ <- as.matrix(duration_data[which(duration_data$`Actors Involved` == "Only Civilian"),
                                    c("LONGITUDE", "LATITUDE")])
xyzt_state_non <- as.matrix(duration_data[which(duration_data$`Actors Involved` == "State and Non-State"),
                                          c("LONGITUDE", "LATITUDE")])
xyzt_state_civ <- as.matrix(duration_data[which(duration_data$`Actors Involved` == "State and Civilian"),
                                          c("LONGITUDE", "LATITUDE")])
xyzt_non_civ <- as.matrix(duration_data[which(duration_data$`Actors Involved` == "Non-State and Civilian"),
                                        c("LONGITUDE", "LATITUDE")])

state_ppp <- as.ppp(xyzt_state, sud_owin) 
nonstate_ppp <- as.ppp(xyzt_nonstate, sud_owin) 
civ_ppp <- as.ppp(xyzt_civ, sud_owin) 
state_non_ppp <- as.ppp(xyzt_state_non, sud_owin) 
state_civ_ppp <- as.ppp(xyzt_state_civ, sud_owin) 
non_civ_ppp <- as.ppp(xyzt_non_civ, sud_owin) 


lgcp_est1 <- lgcp.estpcf(state_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est2 <- lgcp.estpcf(nonstate_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est3 <- lgcp.estpcf(civ_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est4 <- lgcp.estpcf(state_non_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est5 <- lgcp.estpcf(state_civ_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est6 <- lgcp.estpcf(non_civ_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))



#https://rdrr.io/cran/spatstat/man/lgcp.estpcf.html
cov_f <- function(r, lgcp_est, year){
  cov_function <- as.numeric(lgcp_est$par[1]) * exp(-r/as.numeric(lgcp_est$par[2]))
  dat <- as.data.frame(cbind(r, cov_function))
  dat$`Actors Involved` <- year
  return(dat)
}

r <- seq(0,0.45, 0.005) #0.35 for graph
full_dat <- rbind(cov_f(r,lgcp_est1, "Only State"),
                  cov_f(r,lgcp_est2, "Only Non-State"),
                  cov_f(r,lgcp_est3, "Only Civilian"),
                  cov_f(r,lgcp_est4, "State and Non-State"),
                  cov_f(r,lgcp_est5, "State and Civilian"),
                  cov_f(r,lgcp_est6, "Non-State and Civilian"))#,

# FIGURE 7: Covariance function by Actors Involved
ggplot(data=full_dat, aes(x = r, y = cov_function, col = `Actors Involved`)) +geom_line(size =1.0001) + 
  labs(title = "Covariance Function Estimate by Actors Involved", x = "r", y = "C(r)")+ 
  geom_hline(yintercept=0.1, linetype = "dashed", size = 1) + scale_color_brewer(palette="Accent")

eff_range_dat <- NULL
for(i in unique(full_dat$`Actors Involved`)){
  eff_range_ind <- min(which(full_dat$cov_function[which(full_dat$`Actors Involved` == i)] < 0.1))
  eff_range <- full_dat$r[which(full_dat$`Actors Involved` == i)][eff_range_ind]
  eff_range_dat <- rbind(eff_range_dat, c(i, eff_range))
}

# TABLE 4: Effective range estimate by actors involved
eff_range_dat <- as.data.frame(eff_range_dat)
colnames(eff_range_dat) <- c("Actors Involved", "Effective Range")
xtable(eff_range_dat)



#####
## By Conflict Type
#####
#finding effective range for each
xyzt_civil <- as.matrix(duration_data[which(duration_data$EVENT_TYPE == "Violence against civilians"),
                                      c("LONGITUDE", "LATITUDE")])
xyzt_riot <- as.matrix(duration_data[which(duration_data$EVENT_TYPE == "Riots/Protests"),
                                     c("LONGITUDE", "LATITUDE")])
xyzt_battg <- as.matrix(duration_data[which(duration_data$EVENT_TYPE == "Battle-Government regains territory"),
                                      c("LONGITUDE", "LATITUDE")])
xyzt_battnc <- as.matrix(duration_data[which(duration_data$EVENT_TYPE == "Battle-No change of territory"),
                                       c("LONGITUDE", "LATITUDE")])
xyzt_battn <- as.matrix(duration_data[which(duration_data$EVENT_TYPE == "Battle-Non-state actor overtakes territory"),
                                      c("LONGITUDE", "LATITUDE")])


civil_ppp <- as.ppp(xyzt_civil, sud_owin) 
riot_ppp <- as.ppp(xyzt_riot, sud_owin) 
battg_ppp <- as.ppp(xyzt_battg, sud_owin) 
battnc_ppp <- as.ppp(xyzt_battnc, sud_owin) 
battn_ppp <- as.ppp(xyzt_battn, sud_owin) 



lgcp_est1 <- lgcp.estpcf(civil_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est2 <- lgcp.estpcf(riot_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est3 <- lgcp.estpcf(battg_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est4 <- lgcp.estpcf(battnc_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est5 <- lgcp.estpcf(battn_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))



#https://rdrr.io/cran/spatstat/man/lgcp.estpcf.html
cov_f <- function(r, lgcp_est, year){
  cov_function <- as.numeric(lgcp_est$par[1]) * exp(-r/as.numeric(lgcp_est$par[2]))
  dat <- as.data.frame(cbind(r, cov_function))
  dat$`Conflict Type- Specific` <- year
  return(dat)
}

r <- seq(0,0.65, 0.005) #0.55 for plotting
full_dat <- rbind(cov_f(r,lgcp_est1, "Violence against civilians"),
                  cov_f(r,lgcp_est2, "Riots/Protests"),
                  cov_f(r,lgcp_est3, "Battle-Government regains territory"),
                  cov_f(r,lgcp_est4, "Battle-No change of territory"),
                  cov_f(r,lgcp_est5, "Battle-Non-state actor overtakes territory"))#,

# FIGURE 8: Covariance function by conflict type
ggplot(data=full_dat, aes(x = r, y = cov_function, col = `Conflict Type- Specific`)) +geom_line(size =1.0001) + 
  labs(title = "Covariance Function Estimate by Conflict Type", x = "r", y = "C(r)")+ 
  geom_hline(yintercept=0.1, linetype = "dashed", size = 1) + scale_color_brewer(palette="Accent")

eff_range_dat <- NULL
for(i in unique(full_dat$`Conflict Type- Specific`)){
  eff_range_ind <- min(which(full_dat$cov_function[which(full_dat$`Conflict Type- Specific` == i)] < 0.1))
  eff_range <- full_dat$r[which(full_dat$`Conflict Type- Specific` == i)][eff_range_ind]
  eff_range_dat <- rbind(eff_range_dat, c(i, eff_range))
}
# TABLE 5: Effective range by conflict type
eff_range_dat <- as.data.frame(eff_range_dat)
colnames(eff_range_dat) <- c("Conflict Type- Specific", "Effective Range")
xtable(eff_range_dat)


#####
### By Duration
#####
xyzt_sh <- as.matrix(duration_data[which(duration_data$duration == 1),
                                   c("LONGITUDE", "LATITUDE")])
xyzt_long <- as.matrix(duration_data[which(duration_data$duration > 1),
                                     c("LONGITUDE", "LATITUDE")])

nrow(xyzt_sh) + nrow(xyzt_long) == nrow(duration_data)

sh_ppp <- as.ppp(xyzt_sh, sud_owin) 
long_ppp <- as.ppp(xyzt_long, sud_owin) 


lgcp_est1 <- lgcp.estpcf(sh_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))
lgcp_est2 <- lgcp.estpcf(long_ppp, c(var = 1, scale =0.1),
                         covmodel = list(model = "exponential"))


#https://rdrr.io/cran/spatstat/man/lgcp.estpcf.html
cov_f <- function(r, lgcp_est, year){
  cov_function <- as.numeric(lgcp_est$par[1]) * exp(-r/as.numeric(lgcp_est$par[2]))
  dat <- as.data.frame(cbind(r, cov_function))
  dat$`Duration Length` <- year
  return(dat)
}

r <- seq(0,0.35, 0.005)
full_dat <- rbind(cov_f(r,lgcp_est1, "One Day"),
                  cov_f(r,lgcp_est2, "More than One Day"))#,

# FIGURE 9: Covariance function by duration length
ggplot(data=full_dat, aes(x = r, y = cov_function, col = `Duration Length`)) +geom_line(size =1.0001) + 
  labs(title = "Covariance Function Estimate by Duration Length", x = "r", y = "C(r)")+ 
  geom_hline(yintercept=0.1, linetype = "dashed", size = 1) + scale_color_brewer(palette="Accent")#+ theme_bw() #4x7

eff_range_dat <- NULL
for(i in unique(full_dat$`Duration Length`)){
  eff_range_ind <- min(which(full_dat$cov_function[which(full_dat$`Duration Length` == i)] < 0.1))
  eff_range <- full_dat$r[which(full_dat$`Duration Length` == i)][eff_range_ind]
  eff_range_dat <- rbind(eff_range_dat, c(i, eff_range))
}
# TABLE 6: Effective range by duration length
eff_range_dat <- as.data.frame(eff_range_dat)
colnames(eff_range_dat) <- c("Duration Length", "Effective Range")
xtable(eff_range_dat)