
################################################################
################################################################
### Ribbon Separation and Center Estimation for Theseus Maps ###
################################################################
################################################################

### Written by Dr. Lauren J Beesley, PhD
### Contact: lvandervort@lanl.gov

current_directory = '/Users/lbeesley/Desktop/Github_IBEX_RibbonSep'


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

if(!dir.exists(paste0(current_directory,'/Theseus'))){
  dir.create(paste0(current_directory,'/Theseus'))
}
if(!dir.exists(paste0(current_directory,'/Theseus/Sensitivity'))){
  dir.create(paste0(current_directory,'/Theseus/Sensitivity'))
}
if(!dir.exists(paste0(current_directory,'/Theseus/Final_Separations'))){
  dir.create(paste0(current_directory,'/Theseus/Final_Separations'))
}
if(!dir.exists(paste0(current_directory,'/Theseus/Centers'))){
  dir.create(paste0(current_directory,'/Theseus/Centers'))
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
        RESULTSTEMP = Estimate_Ribbon_Center(model_input, WORKING_CENTER = vals_center, ndraws = 100, JACKDRAW = TRUE)
        write.csv(x=RESULTSTEMP[[1]],file = paste0(current_directory,'/Theseus/Centers/Center_',MAP_KEY,'.csv'),quote = FALSE, row.names = FALSE)
      }else{
        RESULTSTEMP = Estimate_Ribbon_Center(model_input, WORKING_CENTER = vals_center, ndraws = 20, JACKDRAW = FALSE)
      }
      vals_center = as.vector(unlist(apply(RESULTSTEMP[[1]],2,mean)))
    }

  }
stopCluster(cl)



###########################
### Perform Separations ### (recommend running overnight)
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
  j=c(1:length(MAPS_TO_RUN[,1])), .packages=c("dplyr","data.table","coda", "pracma", "reshape2", "splines", "pbs",
                                              "truncnorm", "optimx", "circular", "extraDistr", "sn" , "truncnorm",
                                              "mgcv","imagine", "splines2", "zoo", "Matrix"),
  .combine='rbind')%dopar%{
    
    DATA_TO_SEP = read.csv(paste0(current_directory, '/input_theseus.csv'))
    MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                              ESA = unique(DATA_TO_SEP$ESA))
    Time_Group = MAPS_TO_RUN$Time_Group[j]
    ESA = MAPS_TO_RUN$ESA[j]
    MAP_KEY = paste0(Time_Group, '_', ESA)
    MAP_KEY_esa4 = paste0(Time_Group, '_', 4)
    model_ecliptic = DATA_TO_SEP[DATA_TO_SEP$Time_Group == Time_Group &  DATA_TO_SEP$ESA == ESA,]
    model_ecliptic_esa4 = DATA_TO_SEP[DATA_TO_SEP$Time_Group == Time_Group &  DATA_TO_SEP$ESA == 4,]
    
    ### Use Estimated Ribbon Center for ESA 4
    CENTERNAME = paste0(current_directory,'/Theseus/Centers/Center_',MAP_KEY_esa4,'.csv')
    RIBBON_CENTER = read.csv(file = CENTERNAME)  
    RIBBON_CENTER = as.vector(unlist(apply(RIBBON_CENTER,2,mean)))
    
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
    
    ### Estimate ESA 4 Peaks
    model_input_esa4$eligible_angles = as.numeric(model_input_esa4$polar_discrete >= 80 & model_input_esa4$polar_discrete <= 120 )
    SUBSET = model_input_esa4
    SUBSET = SUBSET %>% dplyr::group_by(azimuthal_discrete) %>% 
      dplyr::mutate(phi_max = fitsplines_center_restricted(polar_discrete, ena_rate, polar_discrete[eligible_angles==1])[1])
    SUBSET = SUBSET[!duplicated(SUBSET$azimuthal_discrete),]
    SUBSET = data.frame(SUBSET[order(as.numeric(as.character(SUBSET$azimuthal_discrete))),])
    SUBSET$phi_max[c(0,abs(diff(SUBSET$phi_max)))> 4 | is.na(SUBSET$phi_max)]=NA
    SUBSET$phi_max = as.numeric(forecast::tsclean(ts(SUBSET$phi_max)))
    PREDS_ANCHOR = Get_Loess_PHI_ANCHOR(SUBSET, span = 0.1)
    PHI_ANCHOR <<- as.vector(PREDS_ANCHOR)

    
    VALS = expand.grid(THRESH_VALS = c(0.35, 0.55, 0.75, 0.9), 
                       BOUNDS_VALS = c(30, 35, 40, 45))
    for(k in 1:length(VALS[,1])){
      OUTPUTNAME = paste0(current_directory,'/Theseus/Sensitivity/SeparationRC_',MAP_KEY,'_', VALS$THRESH_VALS[k],'_',VALS$BOUNDS_VALS[k],'.csv')
      OUTPUTNAME_ECLIPTIC = paste0(current_directory,'/Theseus/Sensitivity/SeparationEcliptic_',MAP_KEY,'_', VALS$THRESH_VALS[k],'_',VALS$BOUNDS_VALS[k],'.csv')
      #if(!file.exists(OUTPUTNAME)){
        
        ### Perform Separation
        RESULTS_LIST = GAM_Separation(RIBBON_CENTER, 
                                      THRESH = VALS$THRESH_VALS[k], 
                                      BUFFER = 2, 
                                      RIBBON_BOUNDS = VALS$BOUNDS_VALS[k], 
                                      PEAK_ALL = data.frame(lon_new = AZIMUTHAL, lat_new = PHI_ANCHOR-90), 
                                      AZI_KNOTS = 24,
                                      THESEUS_DELINE = T, 
                                      GAMK = 25, 
                                      PTERMS = 0)
        
        RESULTS = RESULTS_LIST[[1]]
        RESULTS$sobel = get_sobel(RESULTS$azimuthal_discrete, RESULTS$polar_discrete, RESULTS$gdf_final)
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
      #}
    }
  }
stopCluster(cl)


################################################
### Calculate Goodness-of-Separation Metrics ###
################################################

DATA_TO_SEP = read.csv(paste0(current_directory, '/input_theseus.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                          ESA = unique(DATA_TO_SEP$ESA))
GOODNESS_METRIC = NULL
for(j in c(1:length(MAPS_TO_RUN[,1]))){
    Time_Group = MAPS_TO_RUN$Time_Group[j]
    ESA = MAPS_TO_RUN$ESA[j]
    MAP_KEY = paste0(Time_Group, '_', ESA)
    
    VALS = expand.grid(THRESH_VALS = c(0.35, 0.55, 0.75, 0.9), 
                       BOUNDS_VALS = c(30, 35, 40, 45))
    RESULTS_SAVE = NULL
    for(k in 1:length(VALS[,1])){
      OUTPUTNAME = paste0(current_directory,'/Theseus/Sensitivity/SeparationRC_',MAP_KEY,'_', VALS$THRESH_VALS[k],'_',VALS$BOUNDS_VALS[k],'.csv')
      if(file.exists(OUTPUTNAME)){
        RESULTS_TEMP = read.csv(file =OUTPUTNAME)
        RESULTS_SAVE = rbind(RESULTS_SAVE, data.frame(RESULTS_TEMP, THRESH = VALS$THRESH_VALS[k], BOUNDS = VALS$BOUNDS_VALS[k]))
      }
    }
    
    RESULTS_SAVE$eligible = as.numeric(RESULTS_SAVE$BOUNDS >= 30 & RESULTS_SAVE$BOUNDS <= 45)
    
    ### Defining the outer mask region (artificial edge due to edge-of-mask effects)
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(azimuthal_discrete) %>% 
      dplyr::mutate(mask_upper = max(polar_discrete[mask == 1]))
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(azimuthal_discrete) %>% 
      dplyr::mutate(mask_lower = min(polar_discrete[mask == 1]))
    RESULTS_SAVE$mask_upper1 = RESULTS_SAVE$mask_upper + 2
    RESULTS_SAVE$mask_upper2 = RESULTS_SAVE$mask_upper + 4
    RESULTS_SAVE$mask_upper3 = RESULTS_SAVE$mask_upper + 6
    RESULTS_SAVE$mask_lower1 = RESULTS_SAVE$mask_lower - 2
    RESULTS_SAVE$mask_lower2 = RESULTS_SAVE$mask_lower - 4
    RESULTS_SAVE$mask_lower3 = RESULTS_SAVE$mask_lower - 6
    RESULTS_SAVE$edge_outer = as.numeric(RESULTS_SAVE$polar_discrete == RESULTS_SAVE$mask_upper | 
                                           RESULTS_SAVE$polar_discrete == RESULTS_SAVE$mask_upper1 | 
                                           RESULTS_SAVE$polar_discrete == RESULTS_SAVE$mask_upper2 | 
                                           RESULTS_SAVE$polar_discrete == RESULTS_SAVE$mask_upper3 | 
                                           RESULTS_SAVE$polar_discrete == RESULTS_SAVE$mask_lower | 
                                           RESULTS_SAVE$polar_discrete == RESULTS_SAVE$mask_lower1 | 
                                           RESULTS_SAVE$polar_discrete == RESULTS_SAVE$mask_lower2 | 
                                           RESULTS_SAVE$polar_discrete == RESULTS_SAVE$mask_lower3 )
    
    ### Sobel Metric ### (excludes edge_outer to avoid bias due to edge-of-mask effects)
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(azimuthal_discrete, polar_discrete) %>% 
      dplyr::mutate(sobel_med = median(sobel[eligible == 1],na.rm=T))
    RESULTS_SAVE$sobel_stand = RESULTS_SAVE$sobel - RESULTS_SAVE$sobel_med
    RESULTS_SAVE =RESULTS_SAVE %>% dplyr::group_by(THRESH, BOUNDS) %>% 
      dplyr::mutate(sobel_stand_sq = sum(sobel_stand[edge_outer==0]^2, na.rm=T))

    ### Correlation Metric
    COR_SUB = RESULTS_SAVE[RESULTS_SAVE$mask == 1 & RESULTS_SAVE$ribbon_final > 0,]
    COR_SUB = COR_SUB %>% dplyr::group_by(THRESH,BOUNDS) %>% 
      dplyr::mutate(cor_proposed = cor(ribbon_final, gdf_final))
    COR_SUB = COR_SUB[!duplicated(paste0(COR_SUB$THRESH, COR_SUB$BOUNDS)),]
    
    ### Distance Metric
    COR_SUB$distance_metric = sqrt(((COR_SUB$sobel_stand_sq/max(COR_SUB$sobel_stand_sq, na.rm=T))^2)*0.5+ (COR_SUB$cor_proposed^2)*0.5)
    GOODNESS_METRIC = rbind(GOODNESS_METRIC, data.frame(COR_SUB[,c('eligible', 'distance_metric', 'BOUNDS', 'THRESH')], Time_Group = Time_Group, ESA = ESA))
}
write.csv(x=GOODNESS_METRIC,file =paste0(current_directory,'/Theseus/Sensitivity/Goodness_Metrics.csv'),quote = FALSE, row.names = FALSE)



# ggplot(RESULTS_SAVE)+
#   geom_raster(aes(x=azimuthal_discrete, y = polar_discrete, fill = gdf_final, color = gdf_final))+
#   facet_grid(THRESH~BOUNDS)+
#   scale_color_gradientn('ENAs/sec',colors = ibex_palette$hex ,limits = c(0,0.16),na.value = 'black' )+
#   scale_fill_gradientn('ENAs/sec',colors = ibex_palette$hex ,limits = c(0,0.16),na.value = 'black' )+
#   theme_classic()+
#   geom_text(aes(x=200,y=90,label = round(distance_metric,2)), color = 'white', data = COR_SUB)
# ggplot(COR_SUB)+
#   geom_raster(aes(x=THRESH, y=BOUNDS, fill = cor_proposed, color = cor_proposed))
# ggplot(COR_SUB)+
#   geom_raster(aes(x=THRESH, y=BOUNDS, fill = sobel_stand_sq, color = sobel_stand_sq))



#############################
### Get Final Separations ### 
#############################

DATA_TO_SEP = read.csv(paste0(current_directory, '/input_theseus.csv'))
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

    DATA_TO_SEP = read.csv(paste0(current_directory, '/input_theseus.csv'))
    MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                              ESA = unique(DATA_TO_SEP$ESA))
    Time_Group = MAPS_TO_RUN$Time_Group[j]
    ESA = MAPS_TO_RUN$ESA[j]
    MAP_KEY = paste0(Time_Group, '_', ESA)
    MAP_KEY_esa4 = paste0(Time_Group, '_', 4)
    
    RESULTS_ECLIPTIC_TEMP = NULL
    VALS = expand.grid(THRESH_VALS = c(0.6, 0.7, 0.8, 0.9),
                       BOUNDS_VALS = c(20, 25, 30, 35, 40, 45, 50))
    for(k in 1:length(VALS[,1])){
      OUTPUTNAME_ECLIPTIC = paste0(current_directory,'/Theseus/Sensitivity/SeparationEcliptic_',MAP_KEY, '_',VALS$THRESH[k],'_',VALS$BOUNDS[k],'.csv')
      if(file.exists(OUTPUTNAME_ECLIPTIC)){
        RESULTS = read.csv(file =OUTPUTNAME_ECLIPTIC)
        RESULTS = data.frame(RESULTS, Time_Group = Time_Group, ESA = ESA)

        #Custom method to ensure gdf + ribbon = total (i.e., deal with rotational blurring)
        #Approach puts all blurring error within ribbon region + buffer into ribbon
        RESULTS$mask = as.numeric(RESULTS$est_ribbon > 0)
        RESULTS = RESULTS %>% dplyr::group_by(lon) %>%
          dplyr::mutate(min_mask = min(lat[mask == 1], na.rm=T), max_mask = max(lat[mask == 1], na.rm=T))
        RESULTS = RESULTS %>% dplyr::group_by(lon) %>%
          dplyr::mutate(edge_mask_exterior = as.numeric(lat %in% (min_mask - c(2,4,6)) | lat %in% (max_mask + c(2,4,6))),
                        edge_mask_interior = as.numeric(lat %in% (min_mask + c(0,2,4)) | lat %in% (max_mask - c(0,2,4))))
        RESULTS$est_gdf[RESULTS$est_gdf>RESULTS$input_data] = RESULTS$input_data[RESULTS$est_gdf>RESULTS$input_data]
        RESULTS$est_gdf[RESULTS$mask==0 & RESULTS$edge_mask_exterior==0] = RESULTS$input_data[RESULTS$mask==0 & RESULTS$edge_mask_exterior==0]
        RESULTS$est_ribbon = RESULTS$input_data - RESULTS$est_gdf
        
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
    GOODNESS_METRIC = read.csv(file =paste0(current_directory,'/Theseus/Sensitivity/Goodness_Metrics.csv'))
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
    OUTPUTNAME_ECLIPTIC = paste0(current_directory,'/Theseus/Final_Separations/SeparationEcliptic_', MAP_KEY, '.csv')
    write.csv(x=RESULTS_ENSEMBLE[,c('lon', 'lat', 'est_gdf','se_gdf', 'est_ribbon','se_ribbon', 'input_data','input_data_se')],file =OUTPUTNAME_ECLIPTIC,quote = FALSE, row.names = FALSE)
    
    RESULTS = RESULTS_ENSEMBLE
    CENTERNAME = paste0(current_directory,'/Theseus/Centers/Center_',MAP_KEY_esa4,'.csv')
    RIBBON_CENTER = read.csv(file = CENTERNAME)  
    RIBBON_CENTER = as.vector(unlist(apply(RIBBON_CENTER,2,mean)))
    
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
    OUTPUTNAME_RC = paste0(current_directory,'/Theseus/Final_Separations/SeparationRC_', MAP_KEY, '.csv')
    write.csv(x=RESULTS2[,c('lon', 'lat', 'est_gdf','se_gdf', 'est_ribbon','se_ribbon','input_data','input_data_se')],file =OUTPUTNAME_RC,quote = FALSE, row.names = FALSE)
    
}
stopCluster(cl)


#########################
### Merge Separations ###
#########################

DATA_TO_SEP = read.csv(paste0(current_directory, '/input_theseus.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                          ESA = unique(DATA_TO_SEP$ESA))
RESULTS_ECLIPTIC = RESULTS_RC = NULL
for(j in c(1:length(MAPS_TO_RUN[,1]))){
  Time_Group = MAPS_TO_RUN$Time_Group[j]
  ESA = MAPS_TO_RUN$ESA[j]
  MAP_KEY = paste0(Time_Group, '_', ESA)
  
  FILENAME = paste0(current_directory,'/Theseus/Final_Separations/SeparationEcliptic_',MAP_KEY, '.csv')
  if(file.exists(FILENAME)){
    RESULTS = read.csv(FILENAME)
    RESULTS$Time_Group = Time_Group
    RESULTS$ESA = ESA
    RESULTS_ECLIPTIC = rbind(RESULTS_ECLIPTIC, RESULTS)
  }
  
  FILENAME = paste0(current_directory,'/Theseus/Final_Separations/SeparationRC_',MAP_KEY, '.csv')
  if(file.exists(FILENAME)){
    RESULTS2 = read.csv(FILENAME)
    RESULTS2$Time_Group = Time_Group
    RESULTS2$ESA = ESA
    RESULTS_RC = rbind(RESULTS_RC, RESULTS2)
  }
  print(j)
}
### Remove input_data and input_data_se columns for storage reasons
write.csv(x=subset(RESULTS_ECLIPTIC, select = -c(input_data, input_data_se)),file =paste0(current_directory,'/output_theseus_ecliptic.csv'),quote = FALSE, row.names = FALSE)
write.csv(x=subset(RESULTS_RC, select = -c(input_data, input_data_se)),file =paste0(current_directory,'/output_theseus_rc.csv'),quote = FALSE, row.names = FALSE)




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
    CENTERNAME = paste0(current_directory,'/Theseus/Centers/Center_',MAP_KEY_esa4,'.csv')
    RIBBON_CENTER = read.csv(file = CENTERNAME)  
    RIBBON_CENTER = as.vector(unlist(apply(RIBBON_CENTER,2,mean)))

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
    model_ecliptic$sd_ena_rate = model_ecliptic$se_ribbon ### variable names used in Estimate_Ribbon_Center
    ### Code can't handle exactly zero uncertainty
    model_ecliptic$sd_ena_rate[model_ecliptic$sd_ena_rate == 0 | is.na(model_ecliptic$sd_ena_rate)] = min(model_ecliptic$sd_ena_rate[!is.na(model_ecliptic$sd_ena_rate) & model_ecliptic$sd_ena_rate > 0])
    
    
    ### Use Estimated Ribbon Center for ESA 4 as Initial Working Ribbon Center
    CENTERNAME = paste0(current_directory,'/Theseus/Centers/Center_',MAP_KEY_esa4,'.csv')
    vals_center = read.csv(file = CENTERNAME)  
    vals_center = as.vector(unlist(apply(vals_center,2,mean)))
    
    iter_num = 5
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
        RESULTSTEMP = Estimate_Ribbon_Center(model_input, WORKING_CENTER = vals_center, ndraws = 100, JACKDRAW = TRUE)
        write.csv(x=RESULTSTEMP[[1]],file = paste0(current_directory,'/Theseus/Centers_Separated/Center_',MAP_KEY,'.csv'),quote = FALSE, row.names = FALSE)
      }else{
        RESULTSTEMP = Estimate_Ribbon_Center(model_input, WORKING_CENTER = vals_center, ndraws = 20, JACKDRAW = FALSE)
      }
      vals_center = as.vector(unlist(apply(RESULTSTEMP[[1]],2,mean)))
    }
  }
stopCluster(cl)




### Rotate Separations into Updated Ribbon-Centric Frame
DATA_TO_SEP = read.csv(paste0(current_directory, '/output_theseus_ecliptic.csv'))
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
    
    DATA_TO_SEP = read.csv(paste0(current_directory, '/output_theseus_ecliptic.csv'))
    MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                              ESA = unique(DATA_TO_SEP$ESA))
    Time_Group = MAPS_TO_RUN$Time_Group[j]
    ESA = MAPS_TO_RUN$ESA[j]
    MAP_KEY = paste0(Time_Group, '_', ESA)
    MAP_KEY_esa4 = paste0(Time_Group, '_', 4)
    
    ### Output Final Ecliptic Separation
    OUTPUTNAME_ECLIPTIC = paste0(current_directory,'/Theseus/Final_Separations/SeparationEcliptic_', MAP_KEY, '.csv')
    RESULTS_ENSEMBLE = read.csv(OUTPUTNAME_ECLIPTIC)

    RESULTS = RESULTS_ENSEMBLE
    CENTERNAME = paste0(current_directory,'/Theseus/Centers_Separated/Center_',MAP_KEY_esa4,'.csv')
    RIBBON_CENTER = read.csv(file = CENTERNAME)  
    RIBBON_CENTER = as.vector(unlist(apply(RIBBON_CENTER,2,mean)))
    
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
    OUTPUTNAME_RC = paste0(current_directory,'/Theseus/Final_Separations/SeparationRCUpdated_', MAP_KEY, '.csv')
    write.csv(x=RESULTS2[,c('lon', 'lat', 'est_gdf','se_gdf', 'est_ribbon','se_ribbon','input_data','input_data_se')],file =OUTPUTNAME_RC,quote = FALSE, row.names = FALSE)
    
  }
stopCluster(cl)


