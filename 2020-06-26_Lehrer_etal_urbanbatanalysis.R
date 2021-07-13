# source utility functions
source("2020-06-26_Lehrer_etal_bats2018_utility_funcs.R")

## Read in observation and sampling data

# Observation data - y array
yarray <- readRDS("2020-06-26_Lehrer_etal_yarray.rds")

# Sampling data - j matrix
jmat <- readRDS("2020-06-26_Lehrer_etal_jmatrix.rds")

# Create species list to use throughout
species_list <- rownames(yarray)


## Site-level covariates

# Read in data set of site level covariates
sitecovs <- read.csv("2020-06-26_Lehrer_etal_HabitatCoV_2017.csv", 
                     header = TRUE, stringsAsFactors = FALSE)

# Order by sites to match y-array and j-matrix
sitecovs <- sitecovs[order(sitecovs$site),]

# Calculating building density (number of buildings in each buffer)
building_df <- readRDS("2020-06-26_Lehrer_etal_building_xyz_allsites.rds")
nbuildings <- do.call(rbind, lapply(building_df, function(x){nrow(x)}))
nbuildings <- nbuildings[order(row.names(nbuildings))]
nbuildings <- log(nbuildings)

# Calculate median building height
med_height <- do.call(rbind, lapply(building_df, function(x){median(x$Z)}))
med_height <- med_height[order(row.names(med_height))]

# Calculate structural complexity
sci_tmp <- lapply(building_df, function(x) {
  delaunayn(x[,2:4], output.options = "Fa")})
sci <- do.call(rbind, lapply(sci_tmp, function(x){
  sum(x$areas)/(pi*1000^2)}))
sci <- sci[order(row.names(sci))]

# Add building data to dataframe
sitecovs$nbuildings <- nbuildings
sitecovs$med_height <- med_height
sitecovs$sci <- sci

# Calculate Median Sound Pressure Level for each site

# source PAMGuide and customized Meta function
#source("./Meta_auto.R")
#source("./PAMGuide.R")

# # create directory vector to indicate where .wav files our located
# dir <- "/Volumes/Fidino_backup/_Ambient Calls/"
# # temporary file to write output
# dir_write <- "./tmp"

# Create vector of sites to indicate where .wav files were recorded
sites <- as.character(sitecovs$site)[-28]

# # loop through and caluclate Median SPL for each site sampled
# med_SPL <- rep(NA, length(sites))
# for(i in 1:length(sites)){
#   med_SPL[i] <- Meta_auto(directory = dir, dir_write = dir_write, site = sites[i], 
#    atype = "Broadband", calib=1, envi="Air", 
#    ctype="TS", Mh=-40, G=48, vADC=2.5)
#   }

#saveRDS(med_SPL, "./Data/2018-12-19_median_soundlevelpressure.rds")

# Read in SPL data calculated from above
med_SPL <- readRDS("2020-06-26_Lehrer_etal_median_soundlevelpressure.rds")

# Combine with sites
spl <- data.frame(site=as.character(sites), spl = med_SPL)

spl <- spl[order(spl$site), ]

# Add to covariate dataframe
# We are missing a site so we will use the average after scaling
sitecovs$spl <- c(spl[,2], NA)

# Scale the continuous variables to have mean 0 & sd 1
sitecovs$imperv_1000 <- with(sitecovs, scale(imperv_1000))
sitecovs$roads_1000 <- with(sitecovs, scale(roads_1000))
sitecovs$edge_1000 <- with(sitecovs, scale(edge_1000))
sitecovs$dist2water <- with(sitecovs, scale(dist2water))
sitecovs$canopy_1000 <- with(sitecovs, scale(canopy_1000))
sitecovs$sci <- with(sitecovs, scale(sci))
sitecovs$med_height <- with(sitecovs, scale(med_height))
sitecovs$nbuildings <- with(sitecovs, scale(nbuildings))
sitecovs$spl <- with(sitecovs, scale(spl))

# Give WHEA average SPL of urban sites
sitecovs$spl[28] <- 0

# Convert trt to 0-1 (wants to be 1-2, so -1)
sitecovs$trt<-as.numeric(as.factor(sitecovs$trt))-1


## Observation-level covariates

# Read in observation-level covariates for detecion model
obscovs <- read.csv("2020-06-26_Lehrer_etal_DetCovs2013_2017.csv", 
                    header = TRUE, stringsAsFactors = FALSE)
# Order by site to match sitecovs, y-array, and j-matrix
obscovs <- obscovs[order(obscovs$site),]

# Organize average precipitation covariate
precip_mat <- as.matrix(obscovs[,grep("AvgPrec", colnames(obscovs))])

# Calculate the mean and sd
mpcpm <- mean(precip_mat, na.rm = TRUE)
sdpcpm <- sd(precip_mat, na.rm = TRUE)

# Scale mean 0 sd 1
precip <- (precip_mat - mpcpm) / sdpcpm

# Give misssing data the mean value
precip[is.na(precip)] <- 0

# Organize detector model covariate
detectors <- as.matrix(obscovs[,grep("Det", colnames(obscovs))])

# 0's and NA are unknown detectors. Change 0 to NA to combine
detectors[detectors==0] <- NA

# Give all NA's some label. Here we say 11111 to represent unkown or other.
detectors[is.na(detectors)] <- 11111

# Organize detectors as unique values
# Get unique values
lvl <- unique(as.numeric(detectors))

# Apply a unique level to all the detectors
tmp <- data.frame(apply(detectors, 2, factor, levels = lvl))

# Make them numeric
for(i in 1:ncol(tmp)){
  tmp[,i] <- as.numeric(tmp[,i])
}

# Turn into a matrix to use in JAGS
detector <- as.matrix(tmp)


## Set up data for JAGS

# Create some dimensions
nsite <- ncol(yarray)         # number of sites
nsession <- dim(yarray)[3]    # number of sampling sessions
nspec <- length(species_list) # number of species

# Make site covariates a matrix for jags
scovs <- as.matrix(sitecovs[,c("canopy_1000",
                               "edge_1000",
                               "dist2water",
                               "spl",
                               "sci",
                               "med_height",
                               "nbuildings")])

# Take a look at colinearity among covariates
cor(scovs)

# Data list for JAGS
data.list <- list(y = yarray, jmat=as.matrix(jmat),
                  nspec = nspec, nsite = nsite, nsession = nsession,
                  precip = precip, ndetector = length(lvl),
                  sitecovs = scovs, detector = detector,
                  ncof = ncol(scovs))


## JAGS settings

# Create initial values
# Here we use naive occupancy for initial starting values for z in model
zinit <- yarray
zinit[zinit>1] <- 1

# Function that generates starting values for parameters
ncof <- ncol(scovs)
inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = zinit,
      b0 = rnorm(nspec),
      a0 = rnorm(nspec),
      a_pcp = rnorm(nspec),
      b = matrix(rnorm(nspec*ncof), ncol=ncof, nrow = nspec),
      a_det = rnorm(length(lvl)),
      sigma = rgamma(nspec, 1, 1),
      u = matrix(rnorm(nspec*nsession), ncol=nsession, nrow=nspec),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

# Parameters to monitor
params <- c("b0", "b", "a0", "a_pcp", "a_det", "psi", "p", "u", "sigma", 
            "lambda", "lambda_det", "rich")

# Settings
# Adjust chains according to computer specs
n_chains <- parallel::detectCores()-8
adapt_steps <- 50
burn_in <- 150
sample_steps <- 1500
thin_steps <- 5

# Run the JAGS model.
mod_mcmc_lasso <- as.mcmc.list(run.jags(model = "2020-06-26_Lehrer_etal_JAGS_occu_LASSO.R",
                                        monitor = params, 
                                        data = data.list,  
                                        inits = inits, 
                                        n.chains = n_chains,
                                        adapt = adapt_steps,
                                        burnin = burn_in, 
                                        sample = ceiling(sample_steps / n_chains),
                                        thin = thin_steps,
                                        summarise = FALSE,
                                        plots = FALSE,
                                        method = "parallel"))


## Results

# Create a matrix of the beta values
mcmc_matrix <- as.matrix(mod_mcmc_lasso, chains = TRUE)
beta_matrix <- mcmc_matrix[,grep("^b",colnames(mcmc_matrix))]
# Remove model results and mcmc matrix to save memory
rm(mcmc_matrix); rm(mod_mcmc_lasso)

# Explore quantiles for each beta value
quants_beta <- apply(beta_matrix, 2, function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))})

# Print species name to check orders
species_list
colnames(scovs)

# 90% credible intervals 
quants90 <- apply(beta_matrix, 2, function(x) {
  quantile(x, probs = c(0.05, 0.5, 0.95))})

# Proportion of posterior below or above 0
prop_below <- apply(beta_matrix, 2, function(x){
  sum(x < 0)/length(x)})

# Function that plots posterior distributions
plotBetas <- function(x,parameter){
  # Create empty plot
  par(mar=c(5,2,4.5,2), xpd = NA)
  plot(1:10, ylim=c(0.75,nspec+0.25), yaxt="n", xaxt="n", xlim=c(-3,3), xlab="", ylab="",
       pch=16, col='white', bty="n")
  
  # Each plot is a model parameter
  title(parameter, line = 2.5, font = 2, cex.main = 2.5)
  
  # y-axis
  axis(1, at=seq(-3,3,1), cex.axis=1.5)
  
  # y-label
  if(i == 6 | i == 7){
    text(x = 0, y = -0.75, labels = "Model coefficient", font = 2, cex = 2)
  }
  
  # Species lables down the middle
  if(i == 1 |i == 3 |i == 5 | i == 7){
    text(x = rep(3.5, 7) + 0.04, y = (1:7 + 0.4) , labels = species_list, pos = 1, 
         font = 2, cex = 2)
  }
  
  # Create underlines for each species
  for(i in 1:nspec){
    lines(x=c(-3,3),y=c(i,i), lwd = 1, col="lightgrey")
  }
  
  # Intercept line at 0
  lines(x=c(0,0), y=c(0.5,8),lty = 2, lwd = 1.5, col="black")
  
  # Plot transperent background of entire posterior distribution
  for(i in 1:nspec){
    # Create a vector for the particular species_list and parameter
    b_vec <- beta_matrix[,paste0("b[",i,",",x,"]")]
    # Plot distribution
    my_vioplot(b_vec, horizontal = TRUE, at = i, side = "right", add=TRUE, h=0.25,
               wex = 2, colMed = "#FF842A", lwdMed=4, col="#E6F0F5")
  }
  
  # Plot 95% propbability distribution
  for(i in 1:nspec){
    # Create a vector for the particular species_list and parameter
    b_vec <- beta_matrix[,paste0("b[",i,",",x,"]")]
    # Plot only the 95% BCI
    my_vioplot(b_vec[between(b_vec,quantile(b_vec, probs = 0.025), quantile(b_vec, probs = 0.975))],
               horizontal = TRUE, at = i, side = "right", add=TRUE, h=0.25, wex = 2, 
               border = "#0A5EBD", colMed = "#FF842A", lwdMed=5, 
               col="#0A5EBD")
  }
}

# Set up panel layout
m_layout <- matrix(1:8, nrow=4, ncol=2, byrow = TRUE)
# Plot with device
tiff("figure2.tiff", width = 12, height = 16, units = "in",
     res = 300, compression = "lzw")
# Call layout
layout(m_layout, widths=rep(6, ncol(m_layout)), heights = rep(4,nrow(m_layout)))
# Main title names
titles <- c("Canopy cover", "Total forest edge", "Distance to water", 
            bquote(bold("Sound pressure level (dB re 20 "~mu*"Pa)")),
            "Building structural complexity", "Median building height",
            "Log number of buildings")

# Loop through plotBetas function to plot each parameter
for(i in 1:dim(scovs)[2]){
  plotBetas(i,parameter = titles[i])
}

# Turn off device
dev.off()


## Plotting figure that compares sound and water

# Predicting across sound when closest distance to water

# Function to back transform covariates of interest
transformValue <- function(x = NULL, cov_vec = NULL){
  # Mean of covariate vector
  mu <- as.numeric(attributes(cov_vec)[2])
  # Standard deveation of covariate vector
  stdev <- as.numeric(attributes(cov_vec)[3])
  # Transform the real value of interest to scaled value 
  toreturn <- (x-mu)/stdev
  return(toreturn)
}

# Create the scaled value of distance to water = 0
water_best <- transformValue(0, sitecovs$dist2water)

# Create matrix of beta values for covariates of interest
predictImpact <- function(sequence, anthro_cov, cov_vec, spec_num, cov_num){
  
  # Create matrix of values for linear predictor
  # A scaled sequence of the anthropogenic covariate
  new_cov <- transformValue(sequence, cov_vec)
  # New matrix of covariates
  pred_covs <- cbind(1, water_best, new_cov)
  
  # Extract our MCMC samples for the covariates of interests (intercept, water, anthro)
  pred_betas <- beta_matrix[,c(paste0("b0[",spec_num,"]"),
                               paste0("b[",spec_num,",3]"),
                               paste0("b[",spec_num,",",cov_num,"]"))]
  
  # Matrix math to predict across sequence of values
  predict <- pred_betas %*% t(pred_covs)
  
  # Calculate quantile values from MCMC samples for plotting
  quants95 <- apply(predict, 2, quantile, probs=c(0.025, 0.5, 0.975))
  quants90 <- apply(predict, 2, quantile, probs=c(0.05, 0.5, 0.95))
  
  return(list(beta_values = pred_betas, 
              predicted_values = predict,
              quants95 = quants95,
              quants90 = quants90))
}

# Species to plot
spec_to_plot <- c(1,4,6)

# Create list of predictions for the species of interests
prediction_list_spl <- vector("list", length(spec_to_plot))
for(i in 1:length(spec_to_plot)){
  prediction_list_spl[[i]] <-predictImpact(sequence = seq(50, 100, 0.01),
                                           anthro_cov = spl$spl,
                                           cov_vec = sitecovs$spl,
                                           spec_num = spec_to_plot[i],
                                           cov_num = 4)
}

## Plotting
# Set up layout and device
m_layout2 <- matrix(1:2, nrow = 1, ncol = 2, byrow = TRUE)
tiff("figure3.tiff", width = 6, height = 3, units = "in",
     res = 300, compression = "lzw")
# Call layout
layout(m_layout2)
# Set margins
par(mar=c(3.5,3,2.5,0.5))

# Plot the spl level where water no longer has a positive benefit
plot(1:10, ylim=c(30,130), yaxt="n", xaxt="n", xlim=c(0.75,3.25), 
     xlab="", ylab="", pch=16, col='white')
# x-axis
axis(1, at=1:length(spec_to_plot), labels = c("Epfu","Lano","Nyhu"), 
     cex.axis = 0.75, mgp = c(3,0.3,0))
# y-axis
axis(2, at=seq(30,130,10), cex.axis = 0.75, mgp = c(3,0.6,0), las = TRUE)
mtext(text = bquote(bold("Sound pressure level (dB re 20 "~mu*"Pa)")), side = 2, line = 1.75,
      cex = 0.75)

# Set colors for each species
species_col = c("#CA3074", "#F6C667", "#332851")

# Plot median and 95% BCI 
for(i in 1:ncol(spl_quants)){
  points(i,spl_quants[2,i], pch = 16, col = species_col[i], cex = 1)
  arrows(x0 = i, y0 = spl_quants[1,i], y1 = spl_quants[3,i], 
         code= 3, angle = 90, length = 0.03, lwd = 2, col = species_col[i])
}

# Plot sound at water 0m
plot(plogis(prediction_list_spl[[1]]$quants95[2,])~new_seq, type = "l", lwd = 1.5, 
     ylim = c(0,1), xlim = c(50,85), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col="#CA3074")
# x-axis
axis(1, at = seq(50,90,5), cex.axis = 0.75, mgp = c(3,0.3,0))
mtext(text = bquote(bold("Sound pressure level (dB re 20 "~mu*"Pa)")), side = 1, line = 1.5,
      cex = 0.75)
#  y-axis
axis(2, at = seq(0,1,0.2), cex.axis = 0.75, las = TRUE, mgp = c(3,0.6,0))
mtext(text = bquote(bold(psi~"at 0 m from water source")), side = 2, line = 1.75, cex = 0.75)
# 95% credible intervals
polygonsFnc(prediction_list_spl[[1]]$quants95[1,], prediction_list_spl[[1]]$quants95[3,], "#CA3074")
lines(plogis(prediction_list_spl[[2]]$quants95[2,])~new_seq, lwd = 1.5, col="#F6C667")
polygonsFnc(prediction_list_spl[[2]]$quants95[1,], prediction_list_spl[[2]]$quants95[3,], "#F6C667")
lines(plogis(prediction_list_spl[[3]]$quants95[2,])~new_seq, lwd = 1.5, col="#332851")
polygonsFnc(prediction_list_spl[[3]]$quants95[1,], prediction_list_spl[[3]]$quants95[3,], "#332851")

# Indicate 50% probability of occupancy
abline(h=0.5, lty = 2, lwd = 0.75)

# Turn off device
dev.off()


# Gettting point values for graphs
value_fig3 <- rbind(spl_quants[1,], spl_quants[2,], spl_quants[3,])

## Getting values where occupancy crosses 0.5
df_50 <- matrix(NA, ncol= 3, nrow=3)
for(i in 1:length(spec_to_plot)){
  df <- data.frame(new_seq = new_seq, 
                   lo_psi = plogis(prediction_list_spl[[i]]$quants95[1,]),
                   med_psi = plogis(prediction_list_spl[[i]]$quants95[2,]),
                   hi_psi = plogis(prediction_list_spl[[i]]$quants95[3,]))
  df_50[i,] <- c(lo = df[which(df$lo_psi < 0.5)[1],1], 
                 med = df[which(df$med_psi < 0.5)[1],1],
                 hi = df[which(df$hi_psi < 0.5)[1],1])
}
