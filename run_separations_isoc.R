
#######################################
#######################################
### Ribbon Separation for ISOC Maps ###
#######################################
#######################################

current_directory = '/Users/lbeesley/Desktop/Github_IBEX_RibbonSep'
### Written by Dr. Lauren J Beesley, PhD
### Contact: lvandervort@lanl.gov

#########################
### Read in Libraries ###
#########################

library(parallel)
library(foreach)
library(doParallel)
library(dplyr)

###############################
### Source Helper Functions ###
###############################

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
if(!dir.exists(paste0(current_directory,'/ISOC/Sensitivity'))){
  dir.create(paste0(current_directory,'/ISOC/Sensitivity'))
}
if(!dir.exists(paste0(current_directory,'/ISOC/Final_Separations'))){
  dir.create(paste0(current_directory,'/ISOC/Final_Separations'))
}


###########################
### Perform Separations ###  (may take several hours to run)
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
  j=c(1:length(MAPS_TO_RUN[,1])), .packages=c("dplyr","data.table","coda", "pracma", "reshape2", "splines", "pbs",
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
    
    
    VALS = expand.grid(THRESH_VALS = c(0.6, 0.7, 0.8, 0.9),
                       BOUNDS_VALS = c(20, 25, 30, 35, 40, 45, 50))
    for(k in 1:length(VALS[,1])){
      OUTPUTNAME = paste0(current_directory,'/ISOC/Sensitivity/SeparationRC_',MAP_KEY,'_', VALS$THRESH_VALS[k],'_',VALS$BOUNDS_VALS[k],'.csv')
      OUTPUTNAME_ECLIPTIC = paste0(current_directory,'/ISOC/Sensitivity/SeparationEcliptic_',MAP_KEY,'_', VALS$THRESH_VALS[k],'_',VALS$BOUNDS_VALS[k],'.csv')
      if(!file.exists(OUTPUTNAME)){
        
        ### Perform Separation
        RESULTS_LIST = GAM_Separation(RIBBON_CENTER, 
                                      THRESH = VALS$THRESH_VALS[k], 
                                      BUFFER = 1, 
                                      RIBBON_BOUNDS = VALS$BOUNDS_VALS[k], 
                                      PEAK_ALL = data.frame(lon_new = AZIMUTHAL, lat_new = PHI_ANCHOR-90), 
                                      THESEUS_DELINE = F, 
                                      GAMK = 15)
        
        RESULTS = RESULTS_LIST[[1]]
        RESULTS$sobel = get_sobel_isoc(RESULTS$azimuthal_discrete, RESULTS$polar_discrete, RESULTS$gdf_final)
        VARS_TO_SAVE = c('azimuthal_discrete', 'polar_discrete', 'gdf_final', 'ribbon_final','se.fit', 'se.jack', 'mask','ena_rate', 'se_ena_rate','sobel')
        write.csv(x=RESULTS,file =OUTPUTNAME,quote = FALSE, row.names = FALSE)
        
        ### Convert to Ecliptic
        RESULTS$lon = RESULTS$azimuthal_discrete
        RESULTS$lat = RESULTS$polar_discrete-90
        
        RESULTS$ena_rate = RESULTS$gdf_final
        TEMP = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = T, nbreaks = 10)
        TEMP$est_gdf = TEMP$ena_rate
        TEMP = TEMP[order(TEMP$lon, TEMP$lat),]
        
        RESULTS$ena_rate = RESULTS$ribbon_final
        TEMP2 = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = T, nbreaks = 10)
        TEMP2$est_ribbon = TEMP2$ena_rate
        TEMP2 = TEMP2[order(TEMP2$lon, TEMP2$lat),]
        
        RESULTS$ena_rate = RESULTS$ribbon_final/(RESULTS$ribbon_final+RESULTS$gdf_final)
        TEMP3 = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = T, nbreaks = 10)
        TEMP3$prop_ribbon = TEMP3$ena_rate
        TEMP3 = TEMP3[order(TEMP3$lon, TEMP3$lat),]
        
        RESULTS$ena_rate = (RESULTS$se.fit^2)
        TEMP4 = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = T, nbreaks = 10)
        TEMP4$se.fit = sqrt(TEMP4$ena_rate)
        TEMP4 = TEMP4[order(TEMP4$lon, TEMP4$lat),]
        
        RESULTS$ena_rate = (RESULTS$se.jack^2)
        TEMP7 = Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = T, nbreaks = 10)
        TEMP7$se.jack = sqrt(TEMP7$ena_rate)
        TEMP7 = TEMP7[order(TEMP7$lon, TEMP7$lat),]
        
        RESULTS_ECLIPTIC = data.frame(TEMP[,c('lon', 'lat', 'est_gdf' )], est_ribbon=TEMP2[,c('est_ribbon')], 
                                      prop_ribbon=TEMP3[,c('prop_ribbon')],
                                      input_data = model_ecliptic$ena_rate,
                                      input_data_se = model_ecliptic$sd_ena_rate,
                                      se.fit = TEMP4$se.fit,
                                      se.jack = TEMP7$se.jack)
        
        write.csv(x=RESULTS_ECLIPTIC,file =OUTPUTNAME_ECLIPTIC,quote = FALSE, row.names = FALSE)
      }
    }
  }
stopCluster(cl)


################################################
### Calculate Goodness-of-Separation Metrics ###
################################################


DATA_TO_SEP = read.csv(paste0(current_directory, '/input_isoc.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                          ESA = unique(DATA_TO_SEP$ESA))
GOODNESS_METRIC = NULL
for(j in c(1:length(MAPS_TO_RUN[,1]))){
    Time_Group = MAPS_TO_RUN$Time_Group[j]
    ESA = MAPS_TO_RUN$ESA[j]
    MAP_KEY = paste0(Time_Group, '_', ESA)
    
    VALS = expand.grid(THRESH_VALS = c(0.6, 0.7, 0.8, 0.9),
                       BOUNDS_VALS = c(20, 25, 30, 35, 40, 45, 50))
    RESULTS_SAVE = NULL
    for(k in 1:length(VALS[,1])){
      OUTPUTNAME = paste0(current_directory,'/ISOC/Sensitivity/SeparationRC_',MAP_KEY,'_', VALS$THRESH_VALS[k],'_',VALS$BOUNDS_VALS[k],'.csv')
      if(file.exists(OUTPUTNAME)){
        RESULTS_TEMP = read.csv(file =OUTPUTNAME)
        RESULTS_SAVE = rbind(RESULTS_SAVE, data.frame(RESULTS_TEMP, THRESH = VALS$THRESH_VALS[k], BOUNDS = VALS$BOUNDS_VALS[k]))
      }
    }
    
    RESULTS_SAVE$eligible = as.numeric(RESULTS_SAVE$BOUNDS >= 25 & RESULTS_SAVE$BOUNDS <= 40)
    
    ### Sobel Metric ###
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(azimuthal_discrete, polar_discrete) %>% 
      dplyr::mutate(sobel_med = median(sobel[eligible==1],na.rm=T))
    RESULTS_SAVE$sobel_stand = RESULTS_SAVE$sobel - RESULTS_SAVE$sobel_med
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(THRESH, BOUNDS) %>% 
      dplyr::mutate(sobel_stand_sq = sum(sobel_stand^2, na.rm=T))
    
    ### GDF Flux Metric
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(azimuthal_discrete, polar_discrete) %>% 
      dplyr::mutate(gdf_med = median(gdf_final[eligible==1],na.rm=T))
    RESULTS_SAVE$gdf_stand = RESULTS_SAVE$gdf_final - RESULTS_SAVE$gdf_med
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(THRESH, BOUNDS) %>% 
      dplyr::mutate(gdf_stand_sq = sum(gdf_stand^2, na.rm=T))
    
    ### Distance Metric
    COR_SUB = RESULTS_SAVE
    COR_SUB = COR_SUB[!duplicated(paste0(COR_SUB$THRESH, COR_SUB$BOUNDS)),]
    COR_SUB$distance_metric = sqrt(((COR_SUB$sobel_stand_sq/max(COR_SUB$sobel_stand_sq, na.rm=T))^2)*0.5+ ((COR_SUB$gdf_stand_sq/max(COR_SUB$gdf_stand_sq, na.rm=T))^2)*0.5)
    GOODNESS_METRIC = rbind(GOODNESS_METRIC, data.frame(COR_SUB[,c('eligible', 'distance_metric', 'BOUNDS', 'THRESH')], Time_Group = Time_Group, ESA = ESA))
}
write.csv(x=GOODNESS_METRIC,file =paste0(current_directory,'/ISOC/Sensitivity/Goodness_Metrics.csv'),quote = FALSE, row.names = FALSE)




#############################
### Get Final Separations ###
#############################


DATA_TO_SEP = read.csv(paste0(current_directory, '/input_isoc.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                          ESA = unique(DATA_TO_SEP$ESA))

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

    DATA_TO_SEP = read.csv(paste0(current_directory, '/input_isoc.csv'))
    MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                              ESA = unique(DATA_TO_SEP$ESA))
    Time_Group = MAPS_TO_RUN$Time_Group[j]
    ESA = MAPS_TO_RUN$ESA[j]
    MAP_KEY = paste0(Time_Group, '_', ESA)
  
    RESULTS_ECLIPTIC_TEMP = NULL
    VALS = expand.grid(THRESH_VALS = c(0.6, 0.7, 0.8, 0.9),
                       BOUNDS_VALS = c(20, 25, 30, 35, 40, 45, 50))
    for(k in 1:length(VALS[,1])){
      OUTPUTNAME_ECLIPTIC = paste0(current_directory,'/ISOC/Sensitivity/SeparationEcliptic_',MAP_KEY, '_',VALS$THRESH[k],'_',VALS$BOUNDS[k],'.csv')
      if(file.exists(OUTPUTNAME_ECLIPTIC)){
        RESULTS = read.csv(file =OUTPUTNAME_ECLIPTIC)
        RESULTS = data.frame(RESULTS, Time_Group = Time_Group, ESA = ESA)

        #Rotation then Partition Method (note: differs from approach for 2-degree maps)
        RESULTS$est_gdf = RESULTS$input_data*(1-RESULTS$prop_ribbon)
        RESULTS$est_ribbon = RESULTS$input_data*RESULTS$prop_ribbon

        RESULTS$se.fit[RESULTS$est_ribbon == 0 & !is.na(RESULTS$est_ribbon)] = 0
        RESULTS$se.jack[RESULTS$est_ribbon == 0 & !is.na(RESULTS$est_ribbon)] = 0
        
        ### Don't allow additional uncertainty to be greater than uniform
        RESULTS$var_attr_sep = (RESULTS$se.fit^2)+(RESULTS$se.jack^2)
        overage = RESULTS$var_attr_sep - ((RESULTS$input_data^2)/12)
        overage[overage<0 & !is.na(overage)]=0
        overage_jack = ifelse(overage == 0 | is.na(overage), overage, overage - (RESULTS$se.fit^2))
        overage_jack[overage_jack>(RESULTS$se.jack^2) & !is.na(overage_jack)]=(RESULTS$se.jack^2)[overage_jack>(RESULTS$se.jack^2) & !is.na(overage_jack)]
        RESULTS$se.jack.trunc = sqrt((RESULTS$se.jack^2) - overage_jack)
        
        RESULTS$BOUNDS = VALS$BOUNDS[k]
        RESULTS$THRESH = VALS$THRESH[k]
        RESULTS_ECLIPTIC_TEMP = rbind(RESULTS_ECLIPTIC_TEMP, RESULTS)
      }
    }

    
    ############################
    ### Get Ensembled Values ###
    ############################
    GOODNESS_METRIC = read.csv(file =paste0(current_directory,'/ISOC/Sensitivity/Goodness_Metrics.csv'))
    GOODNESS_METRIC = GOODNESS_METRIC[GOODNESS_METRIC$Time_Group == Time_Group & GOODNESS_METRIC$ESA == ESA,]
    RESULTS_ECLIPTIC_TEMP = merge(RESULTS_ECLIPTIC_TEMP, GOODNESS_METRIC[,c('distance_metric', 'eligible', 'BOUNDS', 'THRESH')],
                                  by = c('BOUNDS','THRESH'), all.x = T, all.y = F)
    RESULTS_ECLIPTIC_TEMP$distance_metric[RESULTS_ECLIPTIC_TEMP$distance_metric>1]=0
    RESULTS_ECLIPTIC_TEMP$distance_metric[RESULTS_ECLIPTIC_TEMP$distance_metric<0]=0
    RESULTS_ECLIPTIC_TEMP$distance_diff = RESULTS_ECLIPTIC_TEMP$distance_metric - min(RESULTS_ECLIPTIC_TEMP$distance_metric, na.rm=T)
    RESULTS_ECLIPTIC_TEMP = RESULTS_ECLIPTIC_TEMP[RESULTS_ECLIPTIC_TEMP$distance_metric <= quantile(RESULTS_ECLIPTIC_TEMP$distance_metric, 0.25),]
    
    RESULTS_ECLIPTIC_TEMP$wt = 1-RESULTS_ECLIPTIC_TEMP$distance_diff
    RESULTS_ECLIPTIC_TEMP = RESULTS_ECLIPTIC_TEMP[RESULTS_ECLIPTIC_TEMP$eligible == 1,]
    RESULTS_ECLIPTIC_TEMP = RESULTS_ECLIPTIC_TEMP %>% dplyr::group_by(lon,lat) %>%
      dplyr::mutate(wt_stand = wt/sum(wt))
    
    RESULTS_ECLIPTIC_TEMP = RESULTS_ECLIPTIC_TEMP %>% dplyr::group_by(lon,lat) %>%
      dplyr::mutate(est_gdf_ensemble = stats::weighted.mean(est_gdf, wt_stand),
                    est_ribbon_ensemble = stats::weighted.mean(est_ribbon, wt_stand),
                    se.fit_ensemble = sqrt(stats::weighted.mean(se.fit^2, wt_stand)),
                    se.jack.trunc_ensemble = sqrt(stats::weighted.mean(se.jack.trunc^2, wt_stand)))
    
    RESULTS_ENSEMBLE = RESULTS_ECLIPTIC_TEMP[!duplicated(paste0(RESULTS_ECLIPTIC_TEMP$lat, '_', RESULTS_ECLIPTIC_TEMP$lon)),]
    RESULTS_ENSEMBLE$est_gdf = RESULTS_ENSEMBLE$est_gdf_ensemble
    RESULTS_ENSEMBLE$est_ribbon = RESULTS_ENSEMBLE$est_ribbon_ensemble
    RESULTS_ENSEMBLE$se.fit = RESULTS_ENSEMBLE$se.fit_ensemble
    RESULTS_ENSEMBLE$se.jack.trunc = RESULTS_ENSEMBLE$se.jack.trunc_ensemble
    RESULTS_ENSEMBLE = RESULTS_ENSEMBLE[order(RESULTS_ENSEMBLE$lon, RESULTS_ENSEMBLE$lat),]
    RESULTS_ENSEMBLE$se_ribbon = sqrt(((RESULTS_ENSEMBLE$input_data_se^2)*(RESULTS_ENSEMBLE$est_ribbon/RESULTS_ENSEMBLE$input_data))+(RESULTS_ENSEMBLE$se.fit^2)+(RESULTS_ENSEMBLE$se.jack.trunc^2))
    RESULTS_ENSEMBLE$se_gdf = sqrt(((RESULTS_ENSEMBLE$input_data_se^2)*(RESULTS_ENSEMBLE$est_gdf/RESULTS_ENSEMBLE$input_data))+(RESULTS_ENSEMBLE$se.fit^2)+(RESULTS_ENSEMBLE$se.jack.trunc^2))
    
    ### Output Final Ecliptic Separation
    OUTPUTNAME_ECLIPTIC = paste0(current_directory,'/ISOC/Final_Separations/SeparationEcliptic_', MAP_KEY, '.csv')
    write.csv(x=RESULTS_ENSEMBLE[,c('lon', 'lat', 'est_gdf','se_gdf', 'est_ribbon','se_ribbon', 'input_data','input_data_se')],file =OUTPUTNAME_ECLIPTIC,quote = FALSE, row.names = FALSE)
    
    RESULTS = RESULTS_ENSEMBLE
    
    RIBBON_CENTER = c(221, 39)
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
    
    RESULTS$ena_rate = (RESULTS$input_data_se^2)
    TEMP6= Coord_Transform_Map(MAP =RESULTS, RIBBON_CENTER = RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = FALSE, REVERSE = T, nbreaks = 10)
    TEMP6$se_data = sqrt(TEMP6$ena_rate)
    TEMP6 = TEMP6[order(TEMP6$lon, TEMP6$lat),]
    
    
    RESULTS2 = data.frame(TEMP[,c('lon', 'lat', 'est_gdf')], est_ribbon=TEMP2[,c('est_ribbon')],
                          input_data = TEMP$est_gdf + TEMP2$est_ribbon,
                          input_data_se = TEMP6$se_data,
                          se_ribbon = TEMP5$se_ribbon, se_gdf = TEMP4$se_gdf)
    
    ### Output Final Ribbon-Centered Separation
    OUTPUTNAME_RC = paste0(current_directory,'/ISOC/Final_Separations/SeparationRC_', MAP_KEY, '.csv')
    write.csv(x=RESULTS2[,c('lon', 'lat', 'est_gdf','se_gdf', 'est_ribbon','se_ribbon','input_data','input_data_se')],file =OUTPUTNAME_RC,quote = FALSE, row.names = FALSE)
    
}
stopCluster(cl)


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



