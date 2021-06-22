#################################################################################
# This R script describes a PVA model that incorporates autocorrelation into    #
# the analysis of local demographic information for the Sonoran desert tortoise #
# to produce spatially explicit estimates of demography and viability at        #
# relatively fine spatial scales across the species range in Arizona.           #
#                                                                               #
# The PVA model consists of:                                                    #
# 1. A spatial hierarchical model for estimating survival and transition rates  #
# 2. Range-wide estimate of recruitment (data were insufficient to obtain       # 
#    spatially explicit estimates of recruitment)                               #
# 3. A stage-structured population model that integrates local demographic      #
#    rates into estimates of rate of population change                          #
# 4. A simulation model that forecasts viability using the estimated rate of    #
#    population change and its uncertainty                                      #
#################################################################################


#--- Loads libraries ---#
#install.packages(c("Matrix", "spam", "tripack", "spdep", "gpclib", "deldir", "spatstat", "maptools", "boot", "gplots", "lattice"))
library(spdep)     
library(spatstat)      
library(maptools)
library(R2OpenBUGS)
library(boot)
library(gplots)
library(lattice)
library(coda)


#===========================================================================================#
# 1. Spatial hierarchical model for generating spatially explicit estimates of adult (phi2) #
#    and juvenile survival (phi1) and juvenile-to-adult transitions (psi12)                 #                                #
#===========================================================================================#

#-------------------------------------#
# 1.1 Reads in encounter history data #
#-------------------------------------#

#set working directory
setwd("C:/SDT Spatial PVA/")

y.st <- read.table("Sonoran_desert_tortoise_enc_hist_data_1977-2008.csv", header=TRUE, sep=",", na.strings=c(NA,"NA","N","."))

# limits data to capture histories
CH <- subset(y.st, select = c(y1987:y2008))
rCH <- as.matrix(CH)
rCH[rCH == 0] <- 3 # changes code for unobserved tortoises from 0 to 3 

# creates vectors containing the occasion of first and last capture
n_ind <- dim(rCH)[1]
n_yrs <- dim(rCH)[2]
first.cap <- rep(NA, n_ind)
for(i in 1:n_ind){
   g <- rCH[i,]
   first.cap[i] <- min((1:n_yrs)[!is.na(g) & g!=3])
   }

#-------------------------------------------#
# 1.2 Creates spatial information for model #
#-------------------------------------------#

#--- Initializes information used in the formation of the spatial grid ---#
# latitude and longitude of grid corners for Sonoran desert tortoises
long_w <- -115.5  #eastern California and Nevada
long_e <- -108.5  #eastern New Mexico
lat_n <- 37.5     #southern Utah
lat_s <- 30.5     #northern Mexico

# size of grid cells in degrees
cell_size <- 0.25

# number of e/w cells, n/s cells, and total number of cells in spatial grid
x_cells <- abs(long_w - long_e)/cell_size
y_cells <- abs(lat_n - lat_s)/cell_size
grid_cells <- x_cells*y_cells


#--- Creates spatial grid that encompasses the plot locations ---#
# creates spatial grid that spans from long_w to long_e and from lat_n to lat_s; the grid size is x_cells
# from west to east and y_cells from south to north
xgrid <- GridTopology(c(long_w,lat_s), c(cell_size, cell_size), c(x_cells, y_cells))
xdata <- as.data.frame(c(1:grid_cells))
SDTgrid <- SpatialGridDataFrame(xgrid, xdata)
#image(SDTgrid, col = "lightgrey")
writeAsciiGrid(SDTgrid, "SDTgrid.txt")
modelgrid=readAsciiGrid("SDTgrid.txt")


# define plot locations boundary with clickpoly (make sure to select a 
# large enough area so that it extends to the southern and western 
# boundaries of Arizona)
par(mfrow=c(1,1))
dlocs <- subset(y.st, select = c(long, lat))
head(dlocs)
long=as.vector(dlocs$long)
lat=as.vector(dlocs$lat)
plot(coordinates(modelgrid), pch = 20, cex = 0.5, col = "grey70")
points(dlocs,pch=20, cex=1)
bound <-clickpoly(add=TRUE) # just left click around the perimeter then press the escape key to connect the last 2 pts
#plot(bound)


# defines spatial grid by overlaying convex hull polygon from point locations on model grid
# identifies which points from modelgrid fall in the polygon bound
tmp <- point.in.polygon(coordinates(modelgrid)[,1],coordinates(modelgrid)[,2], bound$bdry[[1]]$x, bound$bdry[[1]]$y)
polygon.mask <- tmp
modelgrid$SDTgrid.txt<-polygon.mask
modelgrid$SDTgrid.txt[modelgrid$SDTgrid.txt == 0] <- NA

# coerces modelgrid to a SpatialPixelsDataFrame class
xp <- as(modelgrid, "SpatialPixelsDataFrame")
#image(xp)

# converts a grid of SpatialPixels into a SpatialPolygons object, which can be
# transformed to a different projection or datum with spTransform in package rgdal
pp <- as(xp, "SpatialPolygons")
#plot(pp)

# exports a SpatialPolygons object into a S-Plus map format to be imported by OpenBUGS
sp2WB(pp, filename="modelgrid.txt")

# create neighbor list for each grid cell, default = 'queen's' neighborhood
neighbors <- poly2nb(pp) 


# subsets plot and grid data
# identifies which monitoring plots and points from modelgrid fall in the polygon boundaries
data.subset<-point.in.polygon(dlocs[,1], dlocs[,2], bound$bdry[[1]]$x, bound$bdry[[1]]$y)
grid.subset<- coordinates(modelgrid)[!is.na(modelgrid$SDTgrid.txt),]

# subsets tortoise data based on plots in spatial grid
data<-y.st[data.subset==1,]


#--- assign tortoises to grid cells ---#

# gridid creates a vector to store the identifier for the grid point that is closest to the plot
# where the tortoise was captured

# spDistsN1 (from package "sp") creates a vector (dvec) for each individual that
# contains the distance between the plot at which the individual was captured
# (dlocs[i,]) to each grid point (grid.subset) 
 
# gridid is the integer grid id to which each plot belongs; ids are assigned to the plot based
# on the index of dvec that corresponds to the grid point for which the distance is smallest
# (i.e., grid point that is closest to the capture plot)  

gridid <- rep(NA,nrow(dlocs))
for(i in 1:nrow(dlocs)){
   #creates a vector of distances between plot i and all grid cells
   dvec <- spDistsN1(grid.subset,as.matrix(dlocs[i,],ncol=2))
   #assigns the grid cell number with the smallest distance to each plot 
   gridid[i] <- (1:length(dvec))[dvec==min(dvec)] 
   }


#-- Objects needed for CAR model ---#
# vector containing the no. of neighbors of each grid cell (see plot(pp) for grid cells and neighbors)
num <- sapply(neighbors, length) 
# describes adjacency structure of grid by converting grid ids from list form into a simple vector
adj <- unlist(neighbors) 
# length of weights = no. adjacent pairs of grid cells 
sumNumNeigh <- length(unlist(neighbors)) 
# total number of grid cells in the polygon
ngrid <- length(num)



#-----------------------------------------------#
# 1.3  OpenBUGS Program Code for Spatial Model  #
#-----------------------------------------------# 

# write file with OpenBUGS model specification
modelFilename = 'spatial multistate model.txt'
cat('
model{

     #--- Establishes priors ---#
 
     # precision parameters for survival and transition probabilities

     tau.phi1 ~ dgamma(0.1, 0.1)
     sd.phi1 <- 1/pow(tau.phi1, 0.5)
     
     tau.phi2 ~ dgamma(0.1, 0.1)
     sd.phi2 <- 1/pow(tau.phi2, 0.5)

     tau.psi12 ~ dgamma(0.1, 0.1)
     sd.psi12 <- 1/pow(tau.psi12, 0.5)


     # creates a vector the same length as adj[] that contains unnormalized weights associated
     # with each pair of areas; sets all weights equal to 1, giving the standard Besag, York,
     # and Mollie (1991) model (see GeoBUGS manual)
     for(j in 1:sumNumNeigh) {weights[j] <- 1}

     # spatially autocorrelated zero-mean random effects
     u.phi1[1:ngrid] ~ car.normal(adj[], weights[], num[], tau.phi1)
     u.phi2[1:ngrid] ~ car.normal(adj[], weights[], num[], tau.phi2)
     u.psi12[1:ngrid] ~ car.normal(adj[], weights[], num[], tau.psi12)

     
     # Priors on p, phi & psi in each year and state
  	 p1 ~ dunif(0,1)  # logit transformation is not necessary because
 	   p2 ~ dunif(0,1)  # there is no interest in spatial variation in p

     phi10 ~ dunif(0,1)   # juvenile survival
     lphi10 <- log(phi10/(1-phi10))
  	 phi20 ~ dunif(0,1)   # adult survival
     lphi20 <- log(phi20/(1-phi20))

     psi120 ~ dunif(0,1)          # juvenile becomes adult
     lpsi120 <- log(psi120/(1-psi120))


     for (i in 1:nind) {
         
         # spatial models for survival and growth probabilities
         logit(phi1[i]) <- lphi10 + u.phi1[gridid[i]]
         logit(phi2[i]) <- lphi20 + u.phi2[gridid[i]]
         logit(psi12[i]) <- lpsi120 + u.psi12[gridid[i]]

         # Define parameters
         # Define probabilities of state S(t+1) given S(t)
         # First index = states at time t-1, last index = states at time t
         ps[1,i,1] <- phi1[i] * (1-psi12[i])
         ps[1,i,2] <- phi1[i] * psi12[i]
         ps[1,i,3] <- 1-phi1[i]
         ps[2,i,1] <- 0
         ps[2,i,2] <- phi2[i]
         ps[2,i,3] <- 1-phi2[i]
         ps[3,i,1] <- 0
         ps[3,i,2] <- 0
         ps[3,i,3] <- 1
         } #closes i loop


     # Define probabilities of O(t) given S(t)
     # First index = states at time t, last index = observations at time t
     po[1,1] <- p1
     po[1,2] <- 0
     po[1,3] <- 1-p1
     po[2,1] <- 0
     po[2,2] <- p2
     po[2,3] <- 1-p2
     po[3,1] <- 0
     po[3,2] <- 0
     po[3,3] <- 1

     # State-space model likelihood
     for (i in 1:nind){
         z[i,f[i]] <- Y[i,f[i]]
         for (t in (f[i]+1):n.occasions){  # loop over time
             # State equation: draw S(t) given S(t-1)
             z[i,t] ~ dcat(ps[z[i,t-1], i,])

             # Observation equation: draw O(t) given S(t)
             Y[i,t] ~ dcat(po[z[i,t],])
             } #closes t loop
     } # #closes i loop
        
        
} #ends model statement
', fill=TRUE, file=modelFilename)

# Bundle data
bugs.data <- list(Y = rCH, f = first.cap, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], gridid = gridid, num = num, adj=adj, sumNumNeigh = sumNumNeigh, ngrid = ngrid)

# Function, which creates a matrix of initial values for latent states z
ch.init <- function(ch, f){
   for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
   return(ch)
   }

# initial values for WinBUGS
inits <- function(){
list(tau.phi1 = 1, tau.phi2 = 1, tau.psi12 = 1, phi10 = runif(1), phi20 = runif(1), psi120 = runif(1), p1 = runif(1), p2 = runif(1), z = ch.init(rCH, first.cap),
     u.phi1 = rnorm(ngrid,0,0), u.phi2 = rnorm(ngrid,0,0), u.psi12 = rnorm(ngrid,0,0))
}

#parameters to trace
parameters <- c("phi10", "phi20", "p1", "p2", "psi120", "tau.phi1", "tau.phi2", "sd.phi1", "sd.phi2", "tau.psi12", "sd.psi12", "u.phi1", "u.phi2", "u.psi12")

# MCMC settings
nchains <- 3
nthin <- 3
niter <- 10000
nburn <- 9000

# Do the MCMC stuff calling OpenBUGS from R and store results in 'out'
out <- bugs (data = bugs.data, inits = inits, parameters.to.save = parameters, model.file = 'spatial multistate model.txt',
             n.chains = nchains, n.thin = nthin, n.iter = niter, n.burnin = nburn, debug = TRUE, OpenBUGS.pgm = "C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe")

# save coda files (Inference --> Sample --> Coda) before closing OpenBUGS [e.g., coda_index.txt, coda_chain1.txt]
# to calculate convergence diagnostics (e.g., BGR)

# examines output
print(out)
plot(out)

# looks at values of gelman diagnostocs
outmcmc <- read.coda.interactive() # at prompt add file names without quotes
gelman.diag(outmcmc)


#--------------------------------------------------#
# 1.4  Makes ASCII grids for mapping phi and psi   #
#--------------------------------------------------# 

grid.data <- function(dem0, u, stat, out) {

     #--- Reconstructs phi and psi estimates (e.g., logit(phi1) <- lphi10 + u.phi1[])---#
     # logit transforms range-wide mean of phi and psi from each iteration
     ldem0<- logit(dem0)
     # adds zero-mean, site-specific, random effect determined by CAR model value from each iteration to to give logit(phi)
     # for each grid cell in each iteration
     ldem.grid <- ldem0 + u
     # converts phi from logit to linear scale
     dem.grid <- inv.logit(ldem.grid)

     #--- Calculates average or sd of phi and psi across iterations for each grid cell ---#
     st.dem.grid <- apply(dem.grid, 2, stat)
     # combines estimate of phi or psi for each grid cell with the coordinates of each grid cell
     st.dem.grid <- cbind(grid.subset, st.dem.grid)
     st.dem.grid <- as.data.frame(st.dem.grid)
     coordinates(st.dem.grid)=c("s1","s2")
     gridded(st.dem.grid)<-TRUE
     st.dem.grid <- as(st.dem.grid, "SpatialGridDataFrame")
     
     
     # write ascii grid file for import into arcgis
     write.asciigrid(st.dem.grid, fname = out)
     
     return(st.dem.grid)
} #end grid.data

mn.phi1.grid <- grid.data(out$sims.list$phi10, out$sims.list$u.phi1, mean, "mn.phi1.grid.asc")
mn.phi2.grid <- grid.data(out$sims.list$phi20, out$sims.list$u.phi2, mean, "mn.phi2.grid.asc")
mn.psi12.grid <- grid.data(out$sims.list$psi120, out$sims.list$u.psi12, mean, "mn.psi12.grid.asc")

sd.phi1.grid <- grid.data(out$sims.list$phi10, out$sims.list$u.phi1, sd, "sd.phi1.grid.asc")
sd.phi2.grid <- grid.data(out$sims.list$phi20, out$sims.list$u.phi2, sd, "sd.phi2.grid.asc")
sd.psi12.grid <- grid.data(out$sims.list$psi120, out$sims.list$u.psi12, sd, "sd.psi12.grid.asc")




#====================================#
# 2. Generates recruitment estimates #
#====================================#

# Data were insufficient to obtain spatially-explicit estimates of recruitment, so we applied a single estimate
# across the entire range. Because data for estimating some components of recruitment were sparse, we also explored
# the effects of multiple range-wide recruitment estimates on viability. Each recruitment estimate is represented 
# as 10000 MCMC samples from posterior distributions. The analysis for obtaining the posterior distributions of 
# recruitment can be found in Campbell et al. (2015).

recruit <- read.table("Sonoran_desert_tortoise_recruitment_data.csv", header=TRUE, sep=",", na.strings=c(NA,"NA","N","."))
# selects recruitment estimate of interest (e.g., R_0.32 selects the posterior distribution with a mean of
# 0.32 one-year-old females per female per year.
recruit <- recruit[, c("R_0.32")]
recruit <- as.matrix(recruit)

#========================================================================================#
# 3. Stage-structure population model that combines samples from posterior distributions # 
#    of survival (phi1 and phi2), transition (psi12), and recruitment into posterior     #
#    distributions of rate of population change (lambda)                                 #
#========================================================================================#

#-----------------------------------------------------------------------------------------#
# 3.1 Takes samples from the posterior distributions of phi1, phi2, psi12 and recruitment #
#-----------------------------------------------------------------------------------------#

#--- Creates the posterior distributions for phi1, phi2, and psi12 ---#

dem.rate <- function(dem0, u.dem) { 
  
  # logit transforms range-wide mean of phi and psi from each iteration
  ldem0<- logit(dem0)
  # adds zero-mean, site-specific, random effect determined by CAR model value from each iteration to give logit(phi)
  # for each grid cell in each iteration
  ldem <- ldem0 + u.dem
  # converts phi from logit to linear scale
  dem <- inv.logit(ldem)
  
  return(dem)
} #end dem.rate


phi1.pd <- dem.rate(out$sims.list$phi10, out$sims.list$u.phi1)
phi2.pd <- dem.rate(out$sims.list$phi20, out$sims.list$u.phi2)
psi12.pd <- dem.rate(out$sims.list$psi120, out$sims.list$u.psi12)

#--- Creates a sample of phi1, phi2, and psi12 from their posterior distribution ---#

# number of grid cells
n.cells <- dim(phi1.pd)[2]
# number of samples taken from the posterior distribution
n.samp <- 3000

# number of elements in sample matrix
n <- n.cells * n.samp

phi1 <- matrix (data = rep(NA, n), nrow = n.samp, ncol = n.cells)
phi2 <- matrix (data = rep(NA, n), nrow = n.samp, ncol = n.cells)
psi12 <- matrix (data = rep(NA, n), nrow = n.samp, ncol = n.cells)

# takes a sample from the posterior distribution of each demographic rate of each grid cell
for (i in 1:n.cells){
  phi1[,i] <- sample (phi1.pd[,i], n.samp, replace = TRUE)
  phi2[,i] <- sample (phi2.pd[,i], n.samp, replace = TRUE)
  psi12[,i] <- sample (psi12.pd[,i], n.samp, replace = TRUE)
}


#--- Creates a sample of recruitment values from their posterior distribution ---#
# data were insuffucient to allow for spatial variation in recruitment 
# so the same sample is used for each grid cell
recruit.s <- sample(recruit, n.samp, replace = TRUE)
recruit <- matrix(data = recruit.s, nrow = n.samp, ncol = n.cells)


#----------------------------------------------------------------#
# 3.2 Estimates lambda using a stage-structured population model #
#----------------------------------------------------------------#

#--- Calculates lambda for each iteration for each grid cell ---#

lambda <- array (rep(NA, n), c(n.samp, n.cells))

for (i in 1:n.samp){
  for (j in 1:n.cells){
    # combines demographic rates into a population projection matrix
    proj.mat <- matrix(c(phi1[i,j]*(1-psi12[i,j]), recruit[i,j],
                         phi1[i,j]* psi12[i,j]   , phi2[i,j])
                       , nrow = 2, ncol = 2, byrow = TRUE)
    # derives eigen values and vectors for the matrix
    eigen.out <- eigen(proj.mat)
    # puts lambda into its proper grid cell
    lambda[i,j] <- eigen.out$values[1]
  }  # end j loop
} # end i loop


#---------------------------------------------#
# 3.3  Makes ASCII grids for mapping lambda   #
#---------------------------------------------# 

#--- Creates an ASCII grid for mean and standard deviation of lambda ---#
grid.data <- function(dat, stat, out) {
  
  # calculate mean and sd of lambda for each grid cell
  st.lambda <- apply(dat, 2, stat)
  # combines mean or sd of lambda for each grid cell with the coordinates of each grid cell
  st.lambda <- cbind(grid.subset, st.lambda)
  st.lambda <- as.data.frame(st.lambda)
  coordinates(st.lambda)=c("s1","s2")
  gridded(st.lambda)<-TRUE
  st.lambda <- as(st.lambda, "SpatialGridDataFrame")
  
  
  # write ascii grid file for import into arcgis
  write.asciigrid(st.lambda, fname = out)
  
  return(st.lambda)
} #end grid.data

mn.lambda <- grid.data(lambda, mean, "mn.lambda.grid.asc")
sd.lambda <- grid.data(lambda, sd, "sd.lambda.grid.asc")


#=====================================================================================#
# 4. Simulation model that forecasts viability using the estimated rate of population #
#    change and its uncertainty.                                                      #
#=====================================================================================#

#--------------------------------------------------------------------#
# 4.1 Estimates probability of extirpation (p.ex) for each grid cell # 
#     for a given starting population size (N0)                      #
#--------------------------------------------------------------------#

#--- Sets model parameters ---#
n.cell <- dim(lambda)[2]

N0 <- 500 # initial population size of each grid cell
ex.th <- 10 # extinction threshold
n.yr <- 100 # number of years the population is projected
n.iter <- 1000 # number of iterations

n1 <- n.yr * n.iter
n2 <- n.iter * n.cell

pr.ex <- rep(NA, n.cell)

#--- Runs simulation model ---#
for (j in 1:n.cell){ 
  
  # sample of lambdas for cell j
  temp <- matrix(sample(lambda[,j], n1, replace = TRUE), nrow = n.iter, ncol = n.yr)

  #cumulative product of lambdas
  cum.pr <- t(apply(temp, 1, cumprod)) 
  
  # population size at time t
  Nt <- N0 * cum.pr
  
  # identifies occasions when population fall below extinction threshold
  Nt[Nt < ex.th] <- 0
  
  # vector containing the first year in which the species dropped below the extinction threshold
  first0 <- rep(NA, n.iter)
  for(i in 1:n.iter){
    g <- Nt[i,]
    first0[i] <- min((1:n.yr)[g == 0])
  }
  
  # iterations without an extinction appear as Inf and are changed to NA
  is.na(first0) <- !is.finite(first0) 
  
  # ensures that once a population falls below the extinction threshold it goes extinct
  # changes all occasions after first drop below the threshold to zero
  for(i in 1:n.iter){
    for (y in 1:n.yr){
      if (!is.na(first0[i]) & y >= first0[i]) Nt[i,y] <- 0
    }
  }
  
  # final population size
  Nf <- Nt[,n.yr]
  
  # proportion of iterations that resulted in extinction
  pr.ex[j] <- length(Nf[Nf==0])/n.iter
}


#------------------------------------------------------------------#
# 4.2  Makes ASCII grids for mapping probabilities of extirpation  #
#------------------------------------------------------------------# 

#--- Creates an ASCII grid ---#
grid.data <- function(dat, out) {
  
  # combines data for each grid cell with the coordinates of each grid cell
  x <- cbind(grid.subset, dat)
  x <- as.data.frame(x)
  coordinates(x)=c("s1","s2")
  gridded(x)<-TRUE
  x <- as(x, "SpatialGridDataFrame")
  
  # write ascii grid file for import into arcgis
  write.asciigrid(x, fname = out)
  
  return(x)
} #end grid.data


grid.data(pr.ex, paste("pr.extinct", N0, ".grid.asc", sep = ""))


