
#####################################
#####################################
### Theseus Map Ribbon Separation ###
#####################################
#####################################


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

if(!dir.exists(paste0(current_directory,'/Theseus'))){
  dir.create(paste0(current_directory,'/Theseus'))
}
if(!dir.exists(paste0(current_directory,'/Theseus/Final_Separations'))){
  dir.create(paste0(current_directory,'/Theseus/Final_Separations'))
}
if(!dir.exists(paste0(current_directory,'/Theseus/Centers'))){
  dir.create(paste0(current_directory,'/Theseus/Centers'))
}
if(!dir.exists(paste0(current_directory,'/Theseus/Rotations'))){
  dir.create(paste0(current_directory,'/Theseus/Rotations'))
}
if(!dir.exists(paste0(current_directory,'/Theseus/Reisenfeld_Separations'))){
  dir.create(paste0(current_directory,'/Theseus/Reisenfeld_Separations'))
}
if(!dir.exists(paste0(current_directory,'/Theseus/Centers_Separated'))){
  dir.create(paste0(current_directory,'/Theseus/Centers_Separated'))
}


#############################
### Unzip files as needed ###
#############################
library(GEOquery)
if(!file.exists(paste0(current_directory, '/input_theseus.csv'))){
  gunzip(filename =paste0(current_directory, '/input_theseus.gz'), destname = paste0(current_directory, '/input_theseus.csv'))
}
if(!file.exists(paste0(current_directory, '/output_theseus_ecliptic_reisenfeld.csv'))){
  gunzip(filename =paste0(current_directory, '/output_theseus_ecliptic_reisenfeld.gz'), destname = paste0(current_directory, '/output_theseus_ecliptic_reisenfeld.csv'))
}



#####################################
### Estimate ESA 4 Ribbon Centers ### 
#####################################

DATA_TO_SEP = read.csv(paste0(current_directory, '/input_theseus.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                          ESA = unique(DATA_TO_SEP$ESA))
MAPS_TO_RUN = MAPS_TO_RUN[MAPS_TO_RUN$ESA == 4,] #only run for ESA 4 maps

myncores <- detectCores()-4
mysockettype <- "SOCK"
cl <- parallel::makeCluster(myncores, type = mysockettype,methods = F)
setDefaultCluster(cl)
registerDoParallel(cl)
candidate_fits <- foreach(
  j=c(1:length(MAPS_TO_RUN[,1])), .packages=c("dplyr","data.table","coda", "pracma", "reshape2", "splines", "pbs",
                                              "truncnorm", "optimx", "circular", "extraDistr", "sn" , "truncnorm",
                                              "mgcv","imagine", "splines2", "zoo", "Matrix"),
  .combine='rbind')%dopar%{
    
    DATA_TO_SEP = read.csv(paste0(current_directory, '/input_theseus.csv'))
    MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                              ESA = unique(DATA_TO_SEP$ESA))
    MAPS_TO_RUN = MAPS_TO_RUN[MAPS_TO_RUN$ESA == 4,] #only run for ESA 4 maps
    Time_Group = MAPS_TO_RUN$Time_Group[j]
    ESA = MAPS_TO_RUN$ESA[j]
    MAP_KEY = paste0(Time_Group, '_', ESA)
    model_ecliptic = DATA_TO_SEP[DATA_TO_SEP$Time_Group == Time_Group &  DATA_TO_SEP$ESA == ESA,]
    
    iter_num = 5
    vals_center = c(221, 39) #initial working ribbon center
    for(iter in c(1:iter_num)){
      ### Rotate Map into Working Ribbon-Centered Frame
      AZIMUTHAL = unique(model_ecliptic$lon)
      POLAR = unique(model_ecliptic$lat)+90
      model_input = Coord_Transform_Map(model_ecliptic, RIBBON_CENTER = vals_center, METHOD = 'Inversion', nbreaks = 10, STAND_DEV = TRUE, REVERSE = FALSE)
      model_input$azimuthal_discrete = model_input$lon
      model_input$polar_discrete = model_input$lat + 90
      model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
      model_input$ESA = ESA
      if(iter == iter_num){
        RESULTSTEMP = Estimate_Ribbon_Center(model_input, WORKING_CENTER = vals_center, DRAW = F)
        write.csv(x=RESULTSTEMP[[1]],file = paste0(current_directory,'/Theseus/Centers/Center_',MAP_KEY,'.csv'),quote = FALSE, row.names = FALSE)
      }else{
        RESULTSTEMP = Estimate_Ribbon_Center(model_input, WORKING_CENTER = vals_center, DRAW = F)
      }
      vals_center = as.vector(unlist(apply(RESULTSTEMP[[1]],2,mean)))
    }
    
  }
stopCluster(cl)



###########################################
### Calculate Reverse Rotation Matrices ###
###########################################
DATA_TO_SEP = read.csv(paste0(current_directory, '/input_theseus.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                          ESA = unique(DATA_TO_SEP$ESA))
MAPS_TO_RUN = MAPS_TO_RUN[MAPS_TO_RUN$ESA == 4,] #only run for ESA 4 maps
for(j in 1:nrow(MAPS_TO_RUN)){
  Time_Group = MAPS_TO_RUN$Time_Group[j]
  ESA = MAPS_TO_RUN$ESA[j]
  MAP_KEY = paste0(Time_Group, '_', ESA)
  RIBBON_CENTER = unlist(read.csv(paste0(current_directory,'/Theseus/Centers/Center_',MAP_KEY,'.csv')))
  model_ecliptic = expand.grid(lat = seq(-89,89,2), lon = seq(1,359,2))
  WEIGHT_MAT = Coord_Transform_Matrix(MAP =model_ecliptic, RIBBON_CENTER = RIBBON_CENTER, REVERSE = T, nbreaks = 10)
  save(list = c('WEIGHT_MAT'),file = paste0(current_directory,'/Theseus/Rotations/RotationMatrix_',MAP_KEY,'.RData'))
  print(j)
}

###########################
### Perform Separations ###
###########################

DATA_TO_SEP = read.csv(paste0(current_directory, '/input_theseus.csv'))
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
                                     "mgcv","imagine", "splines2", "zoo", "Matrix", "readr","R.utils"),
  .combine='rbind')%dopar%{
    
    DATA_TO_SEP = read.csv(paste0(current_directory, '/input_theseus.csv'))
    MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                              ESA = unique(DATA_TO_SEP$ESA))
    Time_Group = MAPS_TO_RUN$Time_Group[j]
    ESA = MAPS_TO_RUN$ESA[j]
    MAP_KEY = paste0(Time_Group, '_', ESA)
    model_ecliptic = DATA_TO_SEP[DATA_TO_SEP$Time_Group == Time_Group &  DATA_TO_SEP$ESA == ESA,]

    ### Rotate Maps into Working Ribbon-Centered Frame
    AZIMUTHAL = unique(model_ecliptic$lon)
    POLAR = unique(model_ecliptic$lat)+90
    ROOTVAL_4 = paste0(Time_Group,'_4')
    RIBBON_CENTER = unlist(read.csv(paste0(current_directory,'/Theseus/Centers/Center_',ROOTVAL_4,'.csv')))
    model_input = Coord_Transform_Map(model_ecliptic, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', nbreaks = 10, STAND_DEV = TRUE, REVERSE = FALSE)
    model_input$azimuthal_discrete = model_input$lon
    model_input$polar_discrete = model_input$lat + 90
    model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
    model_input$ESA = ESA

    ### Estimate ESA 4 map peaks, used to stabilize peak estimation
    model_input$eligible_angles = as.numeric(model_input$polar_discrete >= 80 & model_input$polar_discrete <= 120 )
    SUBSET = model_input
    SUBSET = SUBSET %>% dplyr::group_by(azimuthal_discrete) %>% 
      dplyr::mutate(phi_max = fitsplines_center_restricted(polar_discrete, ena_rate, polar_discrete[eligible_angles==1])[1])
    SUBSET = SUBSET[!duplicated(SUBSET$azimuthal_discrete),]
    SUBSET = data.frame(SUBSET[order(as.numeric(as.character(SUBSET$azimuthal_discrete))),])
    SUBSET$phi_max[c(0,abs(diff(SUBSET$phi_max)))> 4 | is.na(SUBSET$phi_max)]=NA
    SUBSET$phi_max = as.numeric(forecast::tsclean(ts(SUBSET$phi_max)))
    PREDS_ANCHOR = Get_Loess_PHI_ANCHOR(SUBSET, span = 0.1)
    PHI_ANCHOR <<- as.vector(PREDS_ANCHOR)
    
    
    library(Matrix)
    library(R.utils)
    library(readr)
    library(data.table)
    load(paste0(current_directory,'/Theseus/Rotations/RotationMatrix_',paste0(substr(MAP_KEY, 1, nchar(MAP_KEY)-1),'4'),'.RData'))
    dim(WEIGHT_MAT)

    #######################
    ### Get Separations ###
    #######################
    VALS = expand.grid(THRESH_VALS = c(0.35, 0.55, 0.75, 0.9), 
                       BOUNDS_VALS = c(30, 35, 40, 45))
    RESULTS_SAVE = NULL
    for(v in 1:length(VALS[,1])){
      
      ### Perform Separation 
      RESULTS_LIST = GAM_Separation(RIBBON_CENTER, 
                                    THRESH = VALS$THRESH_VALS[v], 
                                    BUFFER = 2, 
                                    RIBBON_BOUNDS = VALS$BOUNDS_VALS[v], 
                                    PEAK_ALL = data.frame(lon_new = AZIMUTHAL, lat_new = PHI_ANCHOR-90), 
                                    AZI_KNOTS = 24,
                                    THESEUS_DELINE = T, 
                                    GAMK = 25)
      
      RESULTS = RESULTS_LIST[[1]][,c('polar_discrete','azimuthal_discrete','gdf_final','ribbon_final','mask')]
      RESULTS = RESULTS[order(RESULTS$azimuthal_discrete, RESULTS$polar_discrete),]
      RESULTS$sobel = get_sobel(RESULTS$azimuthal_discrete, RESULTS$polar_discrete, RESULTS$gdf_final)
      
      RESULTS_SAVE = rbind(RESULTS_SAVE, data.frame(RESULTS[,c('polar_discrete','azimuthal_discrete','gdf_final','ribbon_final','mask','sobel')], THRESH = VALS$THRESH_VALS[v], BOUNDS = VALS$BOUNDS_VALS[v]))
      
      rm(list=c('RESULTS_LIST','RESULTS'))
    }    
    
    gc()
    
    
    ##############################
    ### Goodness of Separation ###
    ##############################
    
    ### Sobel Metric ###
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(azimuthal_discrete, polar_discrete) %>% 
      dplyr::mutate(sobel_med = median(sobel,na.rm=T))#,
    RESULTS_SAVE$sobel_stand = RESULTS_SAVE$sobel - RESULTS_SAVE$sobel_med
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(THRESH, BOUNDS) %>% 
      dplyr::mutate(sobel_stand_sq = sum(sobel_stand^2, na.rm=T))
    
    ### Correlation Metric
    COR_SUB = RESULTS_SAVE[RESULTS_SAVE$mask == 1 & RESULTS_SAVE$ribbon_final > 0 & !is.na(RESULTS_SAVE$ribbon_final),]
    COR_SUB = COR_SUB %>% dplyr::group_by(THRESH,BOUNDS) %>% 
      dplyr::mutate(cor_proposed = cor(ribbon_final, gdf_final, use = 'pairwise.complete.obs'))
    COR_SUB = COR_SUB[!duplicated(paste0(COR_SUB$THRESH, COR_SUB$BOUNDS)),]
    
    ### Distance Metric
    COR_SUB = COR_SUB %>% dplyr::group_by(1) %>% 
      dplyr::mutate(distance_metric = sqrt(0.5*(cor_proposed^2) + 0.5*((sobel_stand_sq/max(sobel_stand_sq,na.rm=T))^2)))
    COR_SUB$eligible = as.numeric(COR_SUB$BOUNDS >= 30 & COR_SUB$BOUNDS <= 45)
    
    ##################################
    ### Subset to Best Separations ###
    ##################################
    RESULTS_SAVE = merge(RESULTS_SAVE, COR_SUB[,c('distance_metric', 'eligible', 'BOUNDS', 'THRESH')],
                         by = c('BOUNDS','THRESH'), all.x = T, all.y = F)
    
    RESULTS_SAVE$distance_diff = RESULTS_SAVE$distance_metric - min(RESULTS_SAVE$distance_metric,na.rm=T)
    RESULTS_SAVE = RESULTS_SAVE[RESULTS_SAVE$distance_metric <= quantile(RESULTS_SAVE$distance_metric, 0.25),]
    
    RESULTS_SAVE$wt = 1-RESULTS_SAVE$distance_diff
    RESULTS_SAVE = RESULTS_SAVE[RESULTS_SAVE$eligible == 1 & !is.na(RESULTS_SAVE$eligible),]
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
      TEMP$est_ribbon = as.numeric(WEIGHT_MAT %*% RESULTS$ribbon_final)
      TEMP$est_gdf = as.numeric(WEIGHT_MAT %*% RESULTS$gdf_final)
      TEMP$prop_ribbon = as.numeric(WEIGHT_MAT %*% (RESULTS$ribbon_final/(RESULTS$ribbon_final+RESULTS$gdf_final)))
      TEMP = TEMP[order(TEMP$lon, TEMP$lat),]
      TEMP$data = TEMP$est_gdf + TEMP$est_ribbon
      TEMP$error = TEMP$data_ecliptic - TEMP$data
      TEMP$mask = as.numeric(TEMP$est_ribbon > 0)
      
      ### Identify exterior and interior edge-of-mask pixels
      TEMP = TEMP %>% dplyr::group_by(lon) %>%
        dplyr::mutate(min_mask = min(lat[mask == 1], na.rm=T), max_mask = max(lat[mask == 1], na.rm=T))
      TEMP = TEMP %>% dplyr::group_by(lon) %>%
        dplyr::mutate(edge_mask_exterior = as.numeric(lat %in% (min_mask - c(2,4,6)) | lat %in% (max_mask + c(2,4,6))),
                      edge_mask_interior = as.numeric(lat %in% (min_mask + c(0,2,4)) | lat %in% (max_mask - c(0,2,4))))
      
      ### Deal with rotational blurring, putting remainder in ribbon
      TEMP$est_gdf[TEMP$est_gdf>TEMP$data_ecliptic & !is.na(TEMP$data_ecliptic)] = TEMP$data_ecliptic[TEMP$est_gdf>TEMP$data_ecliptic & !is.na(TEMP$data_ecliptic)]
      TEMP$est_gdf[TEMP$mask==0 & TEMP$edge_mask_exterior==0] = TEMP$data_ecliptic[TEMP$mask==0 & TEMP$edge_mask_exterior==0]
      TEMP$est_ribbon = TEMP$data_ecliptic - TEMP$est_gdf
      
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
    
    OUTPUTNAME = paste0(current_directory,'/Theseus/Final_Separations/SeparationEcliptic_',MAP_KEY,'.csv')
    write.csv(x=RESULTS_ENSEMBLE[,c('lon','lat','est_gdf','est_ribbon','data_ecliptic','data_ecliptic_se')],file =OUTPUTNAME,quote = FALSE, row.names = FALSE)
    gc()
    
  }
stopCluster(cl)






##############################
### Uncertainty Estimation ###
##############################

### At this point in the data analysis pipeline, we estimate the 
### uncertainties of the separated maps and output the results
### into a large dataset called output_theseus_ecliptic.csv.
### Due to the large computational burden associated with this step,
### we use cluster computing for this step. We provide results 
### without code for this step.


########################################################
### Reestimate ESA 4 Ribbon Centers after Separation ###
########################################################

DATA_TO_SEP = read.csv(paste0(current_directory, '/output_theseus_ecliptic.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                          ESA = unique(DATA_TO_SEP$ESA))
MAPS_TO_RUN = MAPS_TO_RUN[MAPS_TO_RUN$ESA == 4,] #only run for ESA 4 maps

myncores <- detectCores()-4
mysockettype <- "SOCK"
cl <- parallel::makeCluster(myncores, type = mysockettype,methods = F)
setDefaultCluster(cl)
registerDoParallel(cl)
candidate_fits <- foreach(
  j=c(1:length(MAPS_TO_RUN[,1])), .packages=c("dplyr","data.table","coda", "pracma", "reshape2", "splines", "pbs",
                                              "truncnorm", "optimx", "circular", "extraDistr", "sn" , "truncnorm",
                                              "mgcv","imagine", "splines2", "zoo", "Matrix"),
  .combine='rbind')%dopar%{
    
    DATA_TO_SEP = read.csv(paste0(current_directory, '/output_theseus_ecliptic.csv'))
    MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                              ESA = unique(DATA_TO_SEP$ESA))
    MAPS_TO_RUN = MAPS_TO_RUN[MAPS_TO_RUN$ESA == 4,] #only run for ESA 4 maps
    Time_Group = MAPS_TO_RUN$Time_Group[j]
    ESA = MAPS_TO_RUN$ESA[j]
    MAP_KEY = paste0(Time_Group, '_', ESA)
    MAP_KEY_esa4 = paste0(Time_Group, '_', 4)
    
    model_ecliptic = DATA_TO_SEP[DATA_TO_SEP$Time_Group == Time_Group &  DATA_TO_SEP$ESA == ESA,]
    model_ecliptic$ena_rate = model_ecliptic$est_ribbon ### variable names used in Estimate_Ribbon_Center
 
    ### Use Estimated Ribbon Center for ESA 4 as Initial Working Ribbon Center
    vals_center = unlist(read.csv(paste0(current_directory,'/Theseus/Centers/Center_',MAP_KEY_esa4,'.csv')))

    iter_num = 5
    for(iter in c(1:iter_num)){
      ### Rotate Map into Working Ribbon-Centered Frame
      AZIMUTHAL = unique(model_ecliptic$lon)
      POLAR = unique(model_ecliptic$lat)+90
      model_input = Coord_Transform_Map(model_ecliptic, RIBBON_CENTER = vals_center, METHOD = 'Inversion', nbreaks = 10, STAND_DEV = FALSE, REVERSE = FALSE)
      model_input$azimuthal_discrete = model_input$lon
      model_input$polar_discrete = model_input$lat + 90
      model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
      model_input$ESA = ESA
      if(iter == iter_num){
        RESULTSTEMP = Estimate_Ribbon_Center(model_input, WORKING_CENTER = vals_center, DRAW = F)
        write.csv(x=RESULTSTEMP[[1]],file = paste0(current_directory,'/Theseus/Centers_Separated/Center_',MAP_KEY,'.csv'),quote = FALSE, row.names = FALSE)
      }else{
        RESULTSTEMP = Estimate_Ribbon_Center(model_input, WORKING_CENTER = vals_center, DRAW = F)
      }
      vals_center = as.vector(unlist(apply(RESULTSTEMP[[1]],2,mean)))
    }
  }
stopCluster(cl)



##################################################
### Rotate Separations into RC Frame and Merge ###
##################################################
DATA_TO_SEP = read.csv(paste0(current_directory, '/output_theseus_ecliptic.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                          ESA = unique(DATA_TO_SEP$ESA))
for(j in c(1:length(MAPS_TO_RUN[,1]))){
    Time_Group = MAPS_TO_RUN$Time_Group[j]
    ESA = MAPS_TO_RUN$ESA[j]
    MAP_KEY = paste0(Time_Group, '_', ESA)
    MAP_KEY_esa4 = paste0(Time_Group, '_', 4)
    
    ### Output Final Ecliptic Separation

    RESULTS = DATA_TO_SEP[DATA_TO_SEP$Time_Group == Time_Group &  DATA_TO_SEP$ESA == ESA,]
    CENTERNAME = paste0(current_directory,'/Theseus/Centers_Separated/Center_',MAP_KEY_esa4,'.csv')
    RIBBON_CENTER = unlist(read.csv(file = CENTERNAME))

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
    OUTPUTNAME_RC = paste0(current_directory,'/Theseus/Final_Separations/SeparationRC_', MAP_KEY, '.csv')
    write.csv(x=RESULTS2[,c('lon', 'lat', 'est_gdf','se_gdf', 'est_ribbon','se_ribbon')],file =OUTPUTNAME_RC,quote = FALSE, row.names = FALSE)
    print(j)
  }


RESULTS_RC = NULL
for(j in c(1:length(MAPS_TO_RUN[,1]))){
  Time_Group = MAPS_TO_RUN$Time_Group[j]
  ESA = MAPS_TO_RUN$ESA[j]
  MAP_KEY = paste0(Time_Group, '_', ESA)
  
  FILENAME = paste0(current_directory,'/Theseus/Final_Separations/SeparationRC_',MAP_KEY, '.csv')
  if(file.exists(FILENAME)){
    RESULTS2 = read.csv(FILENAME)
    RESULTS2$Time_Group = Time_Group
    RESULTS2$ESA = ESA
    RESULTS_RC = rbind(RESULTS_RC, RESULTS2)
  }
  print(j)
}
write.csv(x=RESULTS_RC,file =paste0(current_directory,'/output_theseus_rc.csv'),quote = FALSE, row.names = FALSE)









#############################################################
### Rotate Reisenfeld Separations into RC Frame and Merge ### 
#############################################################

DATA_TO_ROTATE = read.csv(paste0(current_directory, '/output_theseus_ecliptic_reisenfeld.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_ROTATE$Time_Group),
                          ESA = unique(DATA_TO_ROTATE$ESA))

myncores <- detectCores()-4
mysockettype <- "SOCK"
cl <- parallel::makeCluster(myncores, type = mysockettype,methods = F)
setDefaultCluster(cl)
registerDoParallel(cl)
candidate_fits <- foreach(
  j=c(1:length(MAPS_TO_RUN[,1])), .packages=c("dplyr","data.table","coda", "pracma", "reshape2", "splines", "pbs",
                                              "truncnorm", "optimx", "circular", "extraDistr", "sn" , "truncnorm",
                                              "mgcv","imagine", "splines2", "zoo", "Matrix"),
  .combine='rbind')%dopar%{
    
    DATA_TO_ROTATE = read.csv(paste0(current_directory, '/output_theseus_ecliptic_reisenfeld.csv'))
    MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_ROTATE$Time_Group),
                              ESA = unique(DATA_TO_ROTATE$ESA))
    
    Time_Group = MAPS_TO_RUN$Time_Group[j]
    ESA = MAPS_TO_RUN$ESA[j]
    MAP_KEY = paste0(Time_Group, '_', ESA)
    MAP_KEY_esa4 = paste0(Time_Group, '_', 4)
    RESULTS = DATA_TO_ROTATE[DATA_TO_ROTATE$Time_Group == Time_Group &  DATA_TO_ROTATE$ESA == ESA,]
    
    ### Use Estimated Ribbon Center for ESA 4
    CENTERNAME = paste0(current_directory,'/Theseus/Centers_Separated/Center_',MAP_KEY_esa4,'.csv')
    RIBBON_CENTER = unlist(read.csv(file = CENTERNAME))  

    RESULTS$ena_rate = RESULTS$est_gdf_reisenfeld
    TEMP = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = F, nbreaks = 10)
    TEMP$est_gdf_reisenfeld = TEMP$ena_rate
    TEMP = TEMP[order(TEMP$lon, TEMP$lat),]
    
    RESULTS$ena_rate = RESULTS$est_ribbon_reisenfeld
    TEMP2 = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = F, nbreaks = 10)
    TEMP2$est_ribbon_reisenfeld = TEMP2$ena_rate
    TEMP2 = TEMP2[order(TEMP2$lon, TEMP2$lat),]
    
    RESULTS$ena_rate = (RESULTS$se_gdf_reisenfeld^2)
    TEMP4 = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = T, nbreaks = 10)
    TEMP4$se_gdf_reisenfeld = sqrt(TEMP4$ena_rate)
    TEMP4 = TEMP4[order(TEMP4$lon, TEMP4$lat),]
    
    RESULTS$ena_rate = (RESULTS$se_ribbon_reisenfeld^2)
    TEMP5 = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = T, nbreaks = 10)
    TEMP5$se_ribbon_reisenfeld = sqrt(TEMP5$ena_rate)
    TEMP5 = TEMP5[order(TEMP5$lon, TEMP5$lat),]
    
    RESULTS2 = data.frame(TEMP[,c('lon', 'lat', 'est_gdf_reisenfeld')], est_ribbon_reisenfeld=TEMP2[,c('est_ribbon_reisenfeld')],
                          se_ribbon_reisenfeld = TEMP5$se_ribbon_reisenfeld, se_gdf_reisenfeld = TEMP4$se_gdf_reisenfeld)
    
    ### Output Final Ribbon-Centered Separation
    OUTPUTNAME_RC = paste0(current_directory,'/Theseus/Reisenfeld_Separations/ReisenfeldSeparationRC_', MAP_KEY, '.csv')
    write.csv(x=RESULTS2[,c('lon', 'lat', 'est_gdf_reisenfeld','se_gdf_reisenfeld', 'est_ribbon_reisenfeld','se_ribbon_reisenfeld')],file =OUTPUTNAME_RC,quote = FALSE, row.names = FALSE)
    
  }
stopCluster(cl)



RESULTS_RC = NULL
for(j in c(1:length(MAPS_TO_RUN[,1]))){
  Time_Group = MAPS_TO_RUN$Time_Group[j]
  ESA = MAPS_TO_RUN$ESA[j]
  MAP_KEY = paste0(Time_Group, '_', ESA)
  
  FILENAME = paste0(current_directory,'/Theseus/Reisenfeld_Separations/ReisenfeldSeparationRC_',MAP_KEY, '.csv')
  if(file.exists(FILENAME)){
    RESULTS2 = read.csv(FILENAME)
    RESULTS2$Time_Group = Time_Group
    RESULTS2$ESA = ESA
    RESULTS_RC = rbind(RESULTS_RC, RESULTS2)
  }
  print(j)
}
write.csv(x=RESULTS_RC,file =paste0(current_directory,'/output_theseus_rc_reisenfeld.csv'),quote = FALSE, row.names = FALSE)



