
##################################
##################################
### ISOC Map Ribbon Separation ###
##################################
##################################


#########################
### Read in Libraries ###
#########################

library(dplyr)
library(lubridate)
library(ggplot2)
library(circular)
library(pracma)
library(pbs)
library(reshape2)
library(truncnorm)
library(sn)
library(parallel)
library(foreach)
library(doParallel)
library(Rcpp)
library(optimx)
library(bootstrap)
library(fitConic)

###############################
### Source Helper Functions ###
###############################
current_directory = dirname(rstudioapi::getSourceEditorContext()$path)
library(parallel)
library(foreach)
library(doParallel)
source(paste0(current_directory,'/Internal_Functions/Helper_Functions.R'))
source(paste0(current_directory,'/Internal_Functions/Coordinate_Transform_Functions.R'))
source(paste0(current_directory,'/Internal_Functions/Ribbon_Center_Functions.R'))
source(paste0(current_directory,'/Internal_Functions/Ribbon_Separation_Functions.R'))


#############################
### Set Up File Structure ###
#############################

if(!dir.exists(paste0(current_directory,'/ISOC'))){
  dir.create(paste0(current_directory,'/ISOC'))
}
if(!dir.exists(paste0(current_directory,'/ISOC/Rotations'))){
  dir.create(paste0(current_directory,'/ISOC/Rotations'))
}
if(!dir.exists(paste0(current_directory,'/ISOC/Final_Separations'))){
  dir.create(paste0(current_directory,'/ISOC/Final_Separations'))
}

#########################################
### Calculate Reverse Rotation Matrix ###
#########################################
model_ecliptic = expand.grid(lat = seq(-87,87,6), lon = seq(3,357,6))
WEIGHT_MAT = Coord_Transform_Matrix(MAP =model_ecliptic, RIBBON_CENTER = c(221, 39), REVERSE = T, nbreaks = 10)
save(list = c('WEIGHT_MAT'),file =paste0(current_directory,'/ISOC/Rotations/RotationMatrix.RData'))


###########################
### Perform Separations ###
###########################

DATA_TO_SEP = read.csv(paste0(current_directory, '/input_isoc.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                          ESA = unique(DATA_TO_SEP$ESA))

myncores <- detectCores()-4
mysockettype <- "SOCK"
cl <- parallel::makeCluster(myncores, type = mysockettype,methods = F)
setDefaultCluster(cl)
registerDoParallel(cl)
candidate_fits <- foreach(
  j=1:nrow(MAPS_TO_RUN), .packages=c("dplyr","data.table","coda", "pracma", "reshape2", "splines", "pbs",
                                     "truncnorm", "optimx", "circular", "extraDistr", "sn" , "truncnorm",
                                     "mgcv","imagine", "splines2", "zoo", "Matrix"),
  .combine='rbind')%dopar%{
    DATA_TO_SEP = read.csv(paste0(current_directory, '/input_isoc.csv'))
    MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                              ESA = unique(DATA_TO_SEP$ESA))
    Time_Group = MAPS_TO_RUN$Time_Group[j]
    ESA = MAPS_TO_RUN$ESA[j]
    MAP_KEY = paste0(Time_Group, '_', ESA)
    model_ecliptic = DATA_TO_SEP[DATA_TO_SEP$Time_Group == Time_Group &  DATA_TO_SEP$ESA == ESA,]
    model_ecliptic_esa4 = DATA_TO_SEP[DATA_TO_SEP$Time_Group == Time_Group &  DATA_TO_SEP$ESA == 4,]
    
    
    library(Matrix)
    library(R.utils)
    library(readr)
    library(data.table)
    load(paste0(current_directory,'/ISOC/Rotations/RotationMatrix.RData'))
    dim(WEIGHT_MAT)
    
    
    ### Working Ribbon Center for Separation
    RIBBON_CENTER = c(221,39)
    
    ### Rotate Maps into Working Ribbon-Centered Frame
    AZIMUTHAL = unique(model_ecliptic$lon)
    POLAR = unique(model_ecliptic$lat)+90
    model_input = Coord_Transform_Map(model_ecliptic, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', nbreaks = 10, STAND_DEV = TRUE, REVERSE = FALSE)
    model_input$azimuthal_discrete = model_input$lon
    model_input$polar_discrete = model_input$lat + 90
    model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
    model_input$ESA = ESA
    
    model_input_esa4 = Coord_Transform_Map(model_ecliptic_esa4, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', nbreaks = 10, STAND_DEV = TRUE, REVERSE = FALSE)
    model_input_esa4$azimuthal_discrete = model_input_esa4$lon
    model_input_esa4$polar_discrete = model_input_esa4$lat + 90
    model_input_esa4 = model_input_esa4[order(model_input_esa4$azimuthal_discrete, model_input_esa4$polar_discrete),]

    ### Estimate ESA 4 map peaks, used to stabilize peak estimation
    model_input_esa4$eligible_angles = as.numeric(model_input_esa4$polar_discrete >= 80 & model_input_esa4$polar_discrete <= 120 )
    SUBSET = model_input_esa4
    SUBSET = SUBSET %>% dplyr::group_by(azimuthal_discrete) %>% 
      dplyr::mutate(phi_max = fitsplines_center_restricted(polar_discrete, ena_rate, polar_discrete[eligible_angles==1])[1])
    SUBSET = SUBSET[!duplicated(SUBSET$azimuthal_discrete),]
    SUBSET = data.frame(SUBSET[order(as.numeric(as.character(SUBSET$azimuthal_discrete))),])
    SUBSET$phi_max[c(0,abs(diff(SUBSET$phi_max)))> 12 | is.na(SUBSET$phi_max)]=NA#don't allow a greater than 12 degree jump between neighboring pixels. 
    SUBSET$phi_max = as.numeric(forecast::tsclean(ts(SUBSET$phi_max)))
    PREDS_ANCHOR = Get_Loess_PHI_ANCHOR(SUBSET, span = 0.1)
    PHI_ANCHOR <<- as.vector(PREDS_ANCHOR)
 
    #######################
    ### Get Separations ###
    #######################
    RESULTS_SAVE = NULL
    VALS = expand.grid(THRESH_VALS = c(0.6, 0.7, 0.8, 0.9),
                       BOUNDS_VALS = c(20, 25, 30, 35, 40, 45, 50))
    for(v in 1:length(VALS[,1])){
      
      ### Perform Separation (No bias)
      RESULTS_LIST = GAM_Separation(RIBBON_CENTER, 
                                    THRESH = VALS$THRESH_VALS[v], 
                                    BUFFER = 1, 
                                    RIBBON_BOUNDS = VALS$BOUNDS_VALS[v], 
                                    PEAK_ALL = data.frame(lon_new = AZIMUTHAL, lat_new = PHI_ANCHOR-90), 
                                    THESEUS_DELINE = F, 
                                    GAMK = 15)
      RESULTS = RESULTS_LIST[[1]][,c('polar_discrete','azimuthal_discrete','gdf_final','ribbon_final','mask')]
      RESULTS = RESULTS[order(RESULTS$azimuthal_discrete, RESULTS$polar_discrete),]
      RESULTS$sobel = get_sobel_isoc(RESULTS$azimuthal_discrete, RESULTS$polar_discrete, RESULTS$gdf_final)
      RESULTS_SAVE = rbind(RESULTS_SAVE, data.frame(RESULTS[,c('polar_discrete','azimuthal_discrete','gdf_final','ribbon_final','mask', 'sobel')], THRESH = VALS$THRESH_VALS[v], BOUNDS = VALS$BOUNDS_VALS[v]))
      
      rm(list=c('RESULTS_LIST','RESULTS'))
    }
    
    ##############################
    ### Goodness of Separation ###
    ##############################
    
    ### Sobel Metric ###
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(azimuthal_discrete, polar_discrete) %>% 
      dplyr::mutate(sobel_med = median(sobel,na.rm=T))#,
    RESULTS_SAVE$sobel_stand = RESULTS_SAVE$sobel - RESULTS_SAVE$sobel_med
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(THRESH, BOUNDS) %>% 
      dplyr::mutate(sobel_stand_sq = sum(sobel_stand^2, na.rm=T))
    
    ### GDF Flux Metric
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(azimuthal_discrete, polar_discrete) %>% 
      dplyr::mutate(gdf_med = median(gdf_final,na.rm=T))
    RESULTS_SAVE$gdf_stand = RESULTS_SAVE$gdf_final - RESULTS_SAVE$gdf_med
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(THRESH,BOUNDS) %>% 
      dplyr::mutate(gdf_stand_sq = sum(gdf_stand^2, na.rm=T))
    COR_SUB = RESULTS_SAVE[!duplicated(paste0(RESULTS_SAVE$THRESH, RESULTS_SAVE$BOUNDS)),]
    
    ### Distance Metric
    COR_SUB = COR_SUB %>% dplyr::group_by(1) %>% 
      dplyr::mutate(distance_metric = sqrt(0.5*((gdf_stand_sq/max(gdf_stand_sq))^2) + 0.5*((sobel_stand_sq/max(sobel_stand_sq))^2)))
    COR_SUB$eligible = as.numeric(COR_SUB$BOUNDS >= 25 & COR_SUB$BOUNDS <= 40)
    
    ##################################
    ### Subset to Best Separations ###
    ##################################
    RESULTS_SAVE = merge(RESULTS_SAVE, COR_SUB[,c('distance_metric', 'eligible', 'BOUNDS', 'THRESH')],
                         by = c('BOUNDS','THRESH'), all.x = T, all.y = F)
    
    RESULTS_SAVE$distance_diff = RESULTS_SAVE$distance_metric - min(RESULTS_SAVE$distance_metric)
    RESULTS_SAVE = RESULTS_SAVE[RESULTS_SAVE$distance_metric <= quantile(RESULTS_SAVE$distance_metric, 0.25),]
    
    RESULTS_SAVE$wt = 1-RESULTS_SAVE$distance_diff
    RESULTS_SAVE = RESULTS_SAVE[RESULTS_SAVE$eligible == 1 & !is.na(RESULTS_SAVE$eligible),] #DIFFERENCE NOTED
    RESULTS_SAVE = RESULTS_SAVE %>% dplyr::group_by(azimuthal_discrete, polar_discrete) %>%
      dplyr::mutate(wt_stand = wt/sum(wt)) 
    
    ###########################
    ### Convert to Ecliptic ###
    ###########################
    
    RESULTS_SAVE$lon = RESULTS_SAVE$azimuthal_discrete
    RESULTS_SAVE$lat = RESULTS_SAVE$polar_discrete-90
    
    model_ecliptic = model_ecliptic[order(model_ecliptic$lon, model_ecliptic$lat),]
    VALS_New = data.frame(RESULTS_SAVE[!duplicated(paste0(RESULTS_SAVE$THRESH, '_', RESULTS_SAVE$BOUNDS)),c('THRESH','BOUNDS')])
    RESULTS_ECLIPTIC = NULL
    for(v in 1:length(VALS_New[,1])){
      RESULTS = RESULTS_SAVE[RESULTS_SAVE$THRESH == VALS_New$THRESH[v] & RESULTS_SAVE$BOUNDS == VALS_New$BOUNDS[v],]
      RESULTS = RESULTS[order(RESULTS$azimuthal_discrete, RESULTS$polar_discrete),]
      
      TEMP = model_ecliptic[,c('lon','lat')]
      TEMP$data_ecliptic = model_ecliptic$ena_rate
      TEMP$data_ecliptic_se = model_ecliptic$sd_ena_rate
      TEMP$prop_ribbon = as.numeric(WEIGHT_MAT %*% (RESULTS$ribbon_final/(RESULTS$ribbon_final+RESULTS$gdf_final)))
      TEMP$est_gdf = TEMP$data_ecliptic*(1-TEMP$prop_ribbon)
      TEMP$est_ribbon = TEMP$data_ecliptic*TEMP$prop_ribbon
      
      RESULTS_ECLIPTIC = rbind(RESULTS_ECLIPTIC, data.frame(TEMP, THRESH = VALS_New$THRESH[v], BOUNDS = VALS_New$BOUNDS[v], wt_stand = unique(RESULTS$wt_stand)))
    }    
    
    ################
    ### Ensemble ###
    ################
    RESULTS_ECLIPTIC = RESULTS_ECLIPTIC %>% dplyr::group_by(lon,lat) %>%
      dplyr::mutate(est_gdf_ensemble = stats::weighted.mean(est_gdf, wt_stand),
                    est_ribbon_ensemble = stats::weighted.mean(est_ribbon, wt_stand))
    RESULTS_ENSEMBLE = RESULTS_ECLIPTIC[!duplicated(paste0(RESULTS_ECLIPTIC$lat, '_', RESULTS_ECLIPTIC$lon)),]
    RESULTS_ENSEMBLE$est_gdf = RESULTS_ENSEMBLE$est_gdf_ensemble
    RESULTS_ENSEMBLE$est_ribbon = RESULTS_ENSEMBLE$est_ribbon_ensemble
    RESULTS_ENSEMBLE = RESULTS_ENSEMBLE[order(RESULTS_ENSEMBLE$lon, RESULTS_ENSEMBLE$lat),]
    
    OUTPUTNAME = paste0(current_directory,'/ISOC/Final_Separations/SeparationEcliptic_',MAP_KEY,'.csv')
    write.csv(x=RESULTS_ENSEMBLE[,c('lon','lat','est_gdf','est_ribbon','data_ecliptic','data_ecliptic_se')],file =OUTPUTNAME,quote = FALSE, row.names = FALSE)
    gc()

}

stopCluster(cl)




##################################################
### Rotate Separations to Ribbon-Centric Frame ###
##################################################

### Rotate Separations into Updated Ribbon-Centric Frame
DATA_TO_SEP = read.csv(paste0(current_directory, '/output_isoc_ecliptic.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                          ESA = unique(DATA_TO_SEP$ESA))
for(j in c(1:length(MAPS_TO_RUN[,1]))){
  Time_Group = MAPS_TO_RUN$Time_Group[j]
  ESA = MAPS_TO_RUN$ESA[j]
  MAP_KEY = paste0(Time_Group, '_', ESA)
  MAP_KEY_esa4 = paste0(Time_Group, '_', 4)
  
  ### Output Final Ecliptic Separation
  RESULTS = DATA_TO_SEP[DATA_TO_SEP$Time_Group == Time_Group &  DATA_TO_SEP$ESA == ESA,]
  RIBBON_CENTER = c(221,39)
  
  RESULTS$ena_rate = RESULTS$est_gdf
  TEMP = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = F, nbreaks = 10)
  TEMP$est_gdf = TEMP$ena_rate
  TEMP = TEMP[order(TEMP$lon, TEMP$lat),]
  
  RESULTS$ena_rate = RESULTS$est_ribbon
  TEMP2 = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = F, nbreaks = 10)
  TEMP2$est_ribbon = TEMP2$ena_rate
  TEMP2 = TEMP2[order(TEMP2$lon, TEMP2$lat),]
  
  RESULTS$ena_rate = (RESULTS$se_gdf^2)
  TEMP4 = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = T, nbreaks = 10)
  TEMP4$se_gdf = sqrt(TEMP4$ena_rate)
  TEMP4 = TEMP4[order(TEMP4$lon, TEMP4$lat),]
  
  RESULTS$ena_rate = (RESULTS$se_ribbon^2)
  TEMP5 = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = T, nbreaks = 10)
  TEMP5$se_ribbon = sqrt(TEMP5$ena_rate)
  TEMP5 = TEMP5[order(TEMP5$lon, TEMP5$lat),]
  
  RESULTS2 = data.frame(TEMP[,c('lon', 'lat', 'est_gdf')], est_ribbon=TEMP2[,c('est_ribbon')],
                        se_ribbon = TEMP5$se_ribbon, se_gdf = TEMP4$se_gdf)
  
  ### Output Final Ribbon-Centered Separation
  OUTPUTNAME_RC = paste0(current_directory,'/ISOC/Final_Separations/SeparationRC_', MAP_KEY, '.csv')
  write.csv(x=RESULTS2[,c('lon', 'lat', 'est_gdf','se_gdf', 'est_ribbon','se_ribbon')],file =OUTPUTNAME_RC,quote = FALSE, row.names = FALSE)
  print(j)
}


#########################
### Merge Separations ###
#########################

DATA_TO_SEP = read.csv(paste0(current_directory, '/input_isoc.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                          ESA = unique(DATA_TO_SEP$ESA))
RESULTS_ECLIPTIC = RESULTS_RC = NULL
for(j in c(1:length(MAPS_TO_RUN[,1]))){
  Time_Group = MAPS_TO_RUN$Time_Group[j]
  ESA = MAPS_TO_RUN$ESA[j]
  MAP_KEY = paste0(Time_Group, '_', ESA)
  
  FILENAME = paste0(current_directory,'/ISOC/Final_Separations/SeparationEcliptic_',MAP_KEY, '.csv')
  if(file.exists(FILENAME)){
    RESULTS = read.csv(FILENAME)
    RESULTS$Time_Group = Time_Group
    RESULTS$ESA = ESA
    RESULTS_ECLIPTIC = rbind(RESULTS_ECLIPTIC, RESULTS)
  }
  
  FILENAME = paste0(current_directory,'/ISOC/Final_Separations/SeparationRC_',MAP_KEY, '.csv')
  if(file.exists(FILENAME)){
    RESULTS2 = read.csv(FILENAME)
    RESULTS2$Time_Group = Time_Group
    RESULTS2$ESA = ESA
    RESULTS_RC = rbind(RESULTS_RC, RESULTS2)
  }
  print(j)
}
write.csv(x=RESULTS_ECLIPTIC,file =paste0(current_directory,'/output_isoc_ecliptic.csv'),quote = FALSE, row.names = FALSE)
write.csv(x=RESULTS_RC,file =paste0(current_directory,'/output_isoc_rc.csv'),quote = FALSE, row.names = FALSE)

