# author: Robert Ladwig
# date: 10/07/2020
# title: GLM Workshop 

#### Workshop setup ####
cat("\f")
rm(list = ls())

# if you're using Rstudio:
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
#setwd('./example')

# overview of files for this workshop in the example folder
list.files()


# we will need these packages
library(glmtools)
library(GLM3r)
library(rLakeAnalyzer)
library(tidyverse)
library(lubridate)
library(ncdf4)


# overview of glmtools functions
#   | Function       | Title           |
#   | ------------- |:-------------|
#   | `calibrate_sim` | Calibrates GLM-AED2 variables to improve fit between observed and simulated data |
#   | `compare_to_field` | compare metric for GLM vs field observations |
#   | `get_evaporation`  | get evaporation from GLM simulation |
#   | `get_hypsography` | retrieve hypsography information |
#   | `get_ice` | get ice depth from GLM simulation |
#   | `get_nml_value` | gets a nml value according to an arg_name |
#   | `get_surface_height` | get surface height from GLM simulation |
#   | `get_var` | get a variable from a GLM simulation |
#   | `get_wind` | get wind speed from GLM simulation |
#   | `model_diagnostics` | run diagnostics on model results |
#   | `plot_var_compare` | Plot matching heatmaps for modeled and observed variables |
#   | `plot_var_nc` | plot variables from a GLM simulation |
#   | `plot_var_df` | plot variables from a data.frame |
#   | `read_field_obs` | read in field data into a data.frame |
#   | `read_nml` | read in a GLM simulation `*.nml` file |
#   | `resample_sim` | get subset of time from a generic timeseries data.frame |
#   | `resample_to_field` | match GLM water temperatures with field observations |
#   | `set_nml` | sets values in nml object |
#   | `sim_metrics` | get possible metrics for comparing GLM outputs to field |
#   | `summarize_sim` | creates GLM simulation summary outputs |
#   | `validate_sim` | run diagnostics on model results vs observations |
#   | `write_nml` | write GLM `*.nml` for a GLM simulation |

# check out which R version we're currently using
glm_version()

#### Example 1: reading the namelist file into R  ####
glm_template = 'glm3_rb_Bathy.nml' 
sim_folder <- getwd()
setwd(sim_folder)
sim_folder
out_file <- file.path(sim_folder, "output","output.nc")
out_file
#sim_vars(out_file) #Check the list of simulated variables and units

#field_data <- file.path(sim_folder,"bcs","tt.csv")
field_data <- file.path(sim_folder,"bcs","Dam_2020_Hr.csv")
field_data
file.copy(glm_template, 'glm3.nml', overwrite = TRUE)
nml_file <- file.path(sim_folder, 'glm3.nml')
nml_file <- paste0(sim_folder,"/glm3.nml") # This step sets the nml_file for your 




#### Example 2: first visualisations ####
# run GLM
GLM3r::run_glm(sim_folder, verbose = T)

temp_rmse <- compare_to_field(nc_file = out_file, 
                              field_file = field_data,
                              metric = 'water.temperature', 
                              as_value = FALSE, 
                              precision= 'hours')

print(paste('Total time period (uncalibrated):',round(temp_rmse,2),'deg C RMSE'))

var = 'temp'         # variable to which we apply the calibration procedure
path = getwd()       # simulation path/folder
nml_file = nml_file  # path of the nml configuration file that you want to calibrate on
glm_file = nml_file # # path of the gml configuration file
# which parameter do you want to calibrate? a sensitivity analysis helps
calib_setup <- data.frame('pars' = as.character(c('wind_factor','Kw')),
                          'lb' = c(0.7,0.1),
                          'ub' = c(2,0.8),
                          'x0' = c(1,0.3))
print(calib_setup)
glmcmd = NULL        # command to be used, default applies the GLM3r function
# glmcmd = '/Users/robertladwig/Documents/AquaticEcoDynamics_gfort/GLM/glm'        # custom path to executable
# Optional variables
first.attempt = TRUE # if TRUE, deletes all local csv-files that stores the 
#outcome of previous calibration runs
(period = get_calib_periods(nml_file, ratio = 7)) # define a period for the calibration, 
# this supports a split-sample calibration (e.g. calibration and validation period)
# the ratio value is the ratio of calibration period to validation period, ratio 2 means two parts to calibrationand 1 part to calibration
print(period)
scaling = TRUE       # scaling of the variables in a space of [0,10]; TRUE for CMA-ES
verbose = TRUE
method = 'CMA-ES'    # optimization method, choose either `CMA-ES` or `Nelder-Mead`
metric = 'RMSE'      # objective function to be minimized, here the root-mean square error
target.fit = 0.1     # refers to a target fit of 2.0 degrees Celsius (stops when RMSE is below that)
target.iter = 100    # refers to a maximum run of 20 calibration iterations (stops after that many runs)
plotting = TRUE      # if TRUE, script will automatically save the contour plots
output = out_file    # path of the output file
field_file = field_data # path of the field data
conversion.factor = 1 # conversion factor for the output, e.g. 1 for water temp.
additional.var = var

calibrate_sim(var = 'temp', path = getwd(), 
              field_file = field_file, 
              nml_file = nml_file, 
              glm_file = glm_file, 
              calib_setup = calib_setup, 
              glmcmd = NULL, first.attempt = TRUE, 
              period = period, 
              scaling = TRUE, method = 'CMA-ES', metric = 'RMSE', 
              target.fit = 0.1, target.iter = 100, 
              plotting = TRUE, 
              output = output, 
              verbose = TRUE,
              conversion.factor = 1, 
              additional.var = var)


plot_contour <- function(mod_nc, reference = "surface", h, var, unit, tlevels){
  ncin <- nc_open(mod_nc)
  watdep <- ncvar_get(ncin, "z")
  wattemp <- ncvar_get(ncin, var)
  time <- ncvar_get(ncin, "time")
  time.units <- ncatt_get(ncin, "time", "units")
  #sub("^.*\\s2","",time.units$value)
  time.start <- as.POSIXct(strptime(sub("hours since ","",time.units$value), 
                                    format = "%Y-%m-%d %H:%M:%S"))
  datetime <- time.start + time*3600
  layer <- ncvar_get(ncin, "NS")
  nc_close(ncin)
  watdep[which(watdep == max(watdep))] <- NaN
  wattemp[which(wattemp == max(wattemp))] <- NaN
  sim.watdep <- 0.0*watdep - 999
  for (i in 1:length(datetime)){
    max_depth <- watdep[layer[i],i]
    sim.watdep[1,i] <- max_depth - watdep[1,i]/2 
    for (j in 2:layer[i]){
      sim.watdep[j, i] <- max_depth - (watdep[j,i] + watdep[j-1, i ])/2 
    }
  }
  int.dep <- rev(seq(0.25,round(max(watdep,na.rm = TRUE)),0.1))
  sim.wattemp <- matrix(0, nrow = length(int.dep), ncol= length(datetime))
  for (i in 1:length(datetime)){
    sim.approx <- approx(na.omit(sim.watdep[,i]), na.omit(wattemp[,i]), int.dep)
    sim.wattemp[1:length(int.dep),i] <- sim.approx$y
  }
  if ((median(apply(sim.watdep,2,max,na.rm=TRUE))+median(apply(sim.watdep,2,sd,na.rm=TRUE)))<
      max(sim.watdep[,1],na.rm=TRUE)){
    max.plot.dep <- ceiling(median(apply(sim.watdep,2,max,na.rm=TRUE)))
  } else {
    max.plot.dep <-ceiling(max(sim.watdep[,1],na.rm=TRUE))
  }
  if (max.plot.dep<=20){iter.plot.dep = 1} else if (max.plot.dep<=50){iter.plot.dep = 5} else (iter.plot.dep = 10)
  spectral.colors <-  colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral') )
  inMeter <- function(x) {paste0(x, " m")}
  inCelsius <- function(x) {paste0(x, paste(" ",unit,sep=''))}
  filled.contour(x=as.numeric(datetime), y=(int.dep)*(-1), z = t(sim.wattemp),
                 levels=tlevels,
                 col=rev(spectral.colors(length(tlevels))), main=h, cex = 1.5, cex.main = 3., 
                 plot.axes = {axis(1, labels=format(pretty(datetime,20), "%Y-%b"), at = as.numeric(pretty(datetime,20)), cex.axis=2,las=1.5);
                   axis(2, labels=inMeter(rev(seq(0,max.plot.dep,iter.plot.dep))*(-1)), at = rev(seq(0,max(max.plot.dep),iter.plot.dep))*(-1), cex.axis=1.5)},
                 key.axes = {axis(4,at=unique(tlevels),labels=(inCelsius(unique(tlevels))))})
}


plot_contour(mod_nc=out_file, reference = "surface", h="Test", var="temp", unit="celsius", tlevels=seq(0, 30, 1))

ncin <- nc_open(out_file)
wattemp <- ncvar_get(ncin, var)
plot(wattemp[1,], ylim=c(0,30))
points(wattemp[20,], col = 'red')
