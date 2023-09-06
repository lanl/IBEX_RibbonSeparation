############################################
############################################
### Plotting Selected Simulation Results ###
############################################
############################################

current_directory = '/Users/lbeesley/Desktop/Github_IBEX_RibbonSep'

#########################
### Read in Libraries ###
#########################

library(ggplot2)
library(scales)
library(gridtext)
library(cowplot)
library(MetBrewer)
### IBEX color palette
ibex_palette = read.csv(paste0(current_directory,'/ibex_palette.csv'))

#############################
### Set Up File Structure ###
#############################

if(!dir.exists(paste0(current_directory,'/Plots'))){
  dir.create(paste0(current_directory,'/Plots'))
}

#############################
### Unzip files as needed ###
#############################

if(!file.exists(paste0(current_directory, '/output_simulations_ecliptic.csv'))){
  gunzip(filename =paste0(current_directory, '/output_simulations_ecliptic.gz'), destname = paste0(current_directory, '/output_simulations_ecliptic.csv'))
}
if(!file.exists(paste0(current_directory, '/output_simulations_ecliptic_reisenfeld.csv'))){
  gunzip(filename =paste0(current_directory, '/output_simulations_ecliptic_reisenfeld.gz'), destname = paste0(current_directory, '/output_simulations_ecliptic_reisenfeld.csv'))
}
if(!file.exists(paste0(current_directory, '/output_simulations_rc.csv'))){
  gunzip(filename =paste0(current_directory, '/output_simulations_rc.gz'), destname = paste0(current_directory, '/output_simulations_rc.csv'))
}
if(!file.exists(paste0(current_directory, '/output_simulations_rc_reisenfeld.csv'))){
  gunzip(filename =paste0(current_directory, '/output_simulations_rc_reisenfeld.gz'), destname = paste0(current_directory, '/output_simulations_rc_reisenfeld.csv'))
}

###########################################
###########################################
### Plots in Ecliptic Coordinate System ###
###########################################
###########################################

RESULTS = data.frame(data.table::fread(file = paste0(current_directory,'/output_simulations_ecliptic.csv')))

### Merge Reisenfeld et al. (2021) Results
REISENFELD_ECLIPTIC = data.frame(data.table::fread(file = paste0(current_directory,'/output_simulations_ecliptic_reisenfeld.csv')))
RESULTS = merge(RESULTS, subset(REISENFELD_ECLIPTIC,select = -c(truth_gdf,truth_ribbon,input_data,input_data_se)), by = c('lon', 'lat', 'Map_Stat','Map_Class', 'ESA'), all.x = T, all.y = T)

### Standard IBEX plots show the heliospheric nose as the approx. central longitude, 
### with longitude increasing from right to left
center = 180-(360 - 255)
orig_360 = seq(0,359,85)[-1]
new_360 = orig_360-center
new_360[new_360<0]=new_360[new_360<0]+360
RESULTS$ecliptic_lon_center = RESULTS$lon - center + 0.01
RESULTS$ecliptic_lon_center[RESULTS$ecliptic_lon_center<0.01]=RESULTS$ecliptic_lon_center[RESULTS$ecliptic_lon_center<0.01]+360
RESULTS$ecliptic_lon_center = RESULTS$ecliptic_lon_center + 3

##################################################
### Figure 6: Proposed Separations in Ecliptic ###
##################################################
SUBSET = RESULTS[RESULTS$Map_Stat == 'Ideal Input' & RESULTS$ESA %in% c('esa4', 'noesa'),]
p0 = ggplot(SUBSET)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = input_data, color = input_data))+
  theme_classic()+
  facet_grid(.~Map_Class)+
  xlab('Longitude')+ylab('Latitude')+
  theme(legend.position = 'top')+
  labs(title = 'Input Maps')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))
legend <- cowplot::get_legend( p0 )
p0 = p0+theme(legend.position = "none")
p1 = ggplot(SUBSET)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = est_ribbon, color = est_ribbon))+
  theme_classic()+
  facet_grid(.~Map_Class)+
  xlab('Longitude')+ylab('Latitude')+
  theme(legend.position = 'top')+
  labs(title = 'Ribbon Maps')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(legend.position = "none")
p2 = ggplot(SUBSET)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = est_gdf, color = est_gdf))+
  theme_classic()+
  facet_grid(.~Map_Class)+
  xlab('Longitude')+ylab('Latitude')+
  theme(legend.position = 'top')+
  labs(title = 'GDF Maps')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(legend.position = "none")
p4 = ggplot(SUBSET)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = truth_gdf, color = truth_gdf))+
  theme_classic()+
  facet_grid(.~Map_Class)+
  xlab('Longitude')+ylab('Latitude')+
  theme(legend.position = 'top')+
  labs(title = 'True GDF Maps')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(legend.position = "none")
p3 = cowplot::plot_grid(legend,p0,p1,p2,p4,ncol = 1, nrow = 5, rel_heights = c(0.25,1,1,1,1))
pdf(file = paste0(current_directory, '/Plots/Proposed_Separations_Ecliptic.pdf'), width = 11, height = 11)
p3
dev.off()




########################################################
### Reisenfeld et al. (2021) Separations in Ecliptic ###
########################################################
SUBSET = RESULTS[RESULTS$Map_Stat == 'Ideal Input' & RESULTS$ESA %in% c('esa4', 'noesa'),]
p0 = ggplot(SUBSET)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = input_data, color = input_data))+
  theme_classic()+
  facet_grid(.~Map_Class)+
  xlab('Longitude')+ylab('Latitude')+
  theme(legend.position = 'top')+
  labs(title = 'Input Maps')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))
legend <- cowplot::get_legend( p0 )
p0 = p0+theme(legend.position = "none")
p1 = ggplot(SUBSET)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = est_ribbon_reisenfeld, color = est_ribbon_reisenfeld))+
  theme_classic()+
  facet_grid(.~Map_Class)+
  xlab('Longitude')+ylab('Latitude')+
  theme(legend.position = 'top')+
  labs(title = 'Ribbon Maps')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(legend.position = "none")
p2 = ggplot(SUBSET)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = est_gdf_reisenfeld, color = est_gdf_reisenfeld))+
  theme_classic()+
  facet_grid(.~Map_Class)+
  xlab('Longitude')+ylab('Latitude')+
  theme(legend.position = 'top')+
  labs(title = 'GDF Maps')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(legend.position = "none")
p3 = cowplot::plot_grid(legend,p0,p1,p2,ncol = 1, nrow = 4, rel_heights = c(0.25,1,1,1))
pdf(file = paste0(current_directory, '/Plots/Reisenfeld_Separations_Ecliptic.pdf'), width = 9, height = 8)
p3
dev.off()




#################################################
#################################################
### Plots in Ribbon-Centric Coordinate System ###
#################################################
#################################################

RESULTS_RC = data.frame(data.table::fread(file = paste0(current_directory,'/output_simulations_rc.csv')))
RESULTS_RC$azimuthal_discrete = RESULTS_RC$lon
RESULTS_RC$polar_discrete = RESULTS_RC$lat + 90

### Merge Reisenfeld et al. (2021) Results
REISENFELD_RC = data.frame(data.table::fread(file = paste0(current_directory,'/output_simulations_rc_reisenfeld.csv')))
RESULTS_RC = merge(RESULTS_RC, subset(REISENFELD_RC,select = -c(truth_gdf,truth_ribbon,input_data,input_data_se)), by = c('lon', 'lat', 'Map_Stat','Map_Class', 'ESA'), all.x = T, all.y = T)


### Ribbon limits used for summarizing performance
RIBBON_LOWER = 80
RIBBON_UPPER = 130


###########################################################
### Supp. Figure E1a: Proposed Separations in RC Frame ####
###########################################################
SUBSET = RESULTS_RC[RESULTS_RC$Map_Stat == 'Ideal Input' & RESULTS_RC$ESA %in% c('esa4', 'noesa'),]
p0 = ggplot(SUBSET)+
  geom_tile(aes(x=azimuthal_discrete, y=polar_discrete, fill = input_data, color = input_data))+
  theme_classic()+
  facet_grid(.~Map_Class)+
  xlab('Azimuthal Angle')+ylab('Polar Angle')+
  theme(legend.position = 'top')+
  labs(title = 'Input Maps')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))
legend <- cowplot::get_legend( p0 )
p0 = p0+theme(legend.position = "none")
p1 = ggplot(SUBSET)+
  geom_tile(aes(x=azimuthal_discrete, y=polar_discrete, fill = est_ribbon, color = est_ribbon))+
  theme_classic()+
  facet_grid(.~Map_Class)+
  xlab('Azimuthal Angle')+ylab('Polar Angle')+
  theme(legend.position = 'top')+
  labs(title = 'Ribbon Maps')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(legend.position = "none")
p2 = ggplot(SUBSET)+
  geom_tile(aes(x=azimuthal_discrete, y=polar_discrete, fill = est_gdf, color = est_gdf))+
  theme_classic()+
  facet_grid(.~Map_Class)+
  xlab('Azimuthal Angle')+ylab('Polar Angle')+
  theme(legend.position = 'top')+
  labs(title = 'GDF Maps')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(SUBSET$input_data,na.rm=T)), na.value = 'black' )+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(legend.position = "none")
p3 = cowplot::plot_grid(legend,p0,p1,p2,ncol = 1, nrow = 4, rel_heights = c(0.25,1,1,1))
pdf(file = paste0(current_directory, '/Plots/Proposed_Separations_RC.pdf'), width = 9, height = 8)
p3
dev.off()


########################################
### Figure 7a: Average Percent Error ###
########################################
RESULTS_RC$bias = abs(RESULTS_RC$est_gdf - RESULTS_RC$truth_gdf)*100
RESULTS_RC$percent_error = 100*(RESULTS_RC$est_gdf - RESULTS_RC$truth_gdf)/RESULTS_RC$truth_gdf
RESULTS_RC$abs_percent_error = abs(RESULTS_RC$percent_error)
ERROR_OVERALL = aggregate(cbind(abs_percent_error, bias)~Map_Class + Map_Stat, FUN = function(x){mean(x,na.rm=T)}, na.action = na.pass, data = RESULTS_RC[RESULTS_RC$polar_discrete >= RIBBON_LOWER &RESULTS_RC$polar_discrete <= RIBBON_UPPER, ])
RESULTS_RC$bias = abs(RESULTS_RC$est_gdf_reisenfeld - RESULTS_RC$truth_gdf)*100
RESULTS_RC$percent_error_reisenfeld = 100*(RESULTS_RC$est_gdf_reisenfeld - RESULTS_RC$truth_gdf)/RESULTS_RC$truth_gdf
RESULTS_RC$abs_percent_error = abs(RESULTS_RC$percent_error_reisenfeld)
ERROR_OVERALL_reisenfeld = aggregate(cbind(abs_percent_error, bias)~Map_Class + Map_Stat, FUN = function(x){mean(x,na.rm=T)}, na.action = na.pass, data = RESULTS_RC[RESULTS_RC$polar_discrete >= RIBBON_LOWER &RESULTS_RC$polar_discrete <= RIBBON_UPPER, ])
TO_PLOT = rbind(data.frame(Type = 'Proposed', ERROR_OVERALL),
                data.frame(Type = 'Reisenfeld\net al. (2021)',ERROR_OVERALL_reisenfeld ))
TO_PLOT = TO_PLOT[TO_PLOT$Map_Stat %in% c('Estimated Input', 'Estimated Input, 3x', 'Ideal Input'),]
TO_PLOT = TO_PLOT %>% dplyr::group_by(Map_Class, Map_Stat) %>% dplyr::mutate(best_val = Type[which.min(abs_percent_error)])
p1 = ggplot(TO_PLOT)+
  geom_tile(aes(x=factor(Type), y = factor(Map_Class), fill =factor(as.numeric(Type == best_val))), color = 'black', alpha = 0.5)+
  geom_text(aes(x=factor(Type), y = factor(Map_Class), label =round(abs_percent_error,2)), color = 'black', size = 3)+
  geom_text(aes(x=factor(Type), y = factor(Map_Class), label =round(abs_percent_error,2)), color = 'black', size = 3, fontface = 'bold', data = TO_PLOT[TO_PLOT$Type == TO_PLOT$best_val,])+
  xlab('')+ylab('')+
  labs(title = 'Average Absolute Percent Error')+
  scale_fill_manual(name = '', values = c('tomato','olivedrab2'), limits = factor(c(0,1)))+
  guides(color = 'none', fill = 'none')+
  theme(panel.background = element_rect(fill = NA, color = "black"),
        strip.background =element_rect(fill="white", color = "black"))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), limits = c('S3: Varying Profile', 'S2: Spatial Retention', 'S1: Weak Scattering'))+
  facet_grid(.~Map_Stat)+
  theme(axis.ticks = element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())
pdf(file = paste0(current_directory, '/Plots/PercentError.pdf'), width = 7, height = 2)
p1
dev.off()
# 


########################################
### Figure 7b: Spearman Correlations ###
########################################
ROC_SUB = RESULTS_RC[RESULTS_RC$polar_discrete >= RIBBON_LOWER & RESULTS_RC$polar_discrete <= RIBBON_UPPER,]
ROC_SUB = ROC_SUB %>% dplyr::group_by(Map_Stat, Map_Class) %>% 
  dplyr::mutate(cor_proposed = cor(truth_gdf, est_gdf, method = 'spearman'),
                cor_reisenfeld = cor(truth_gdf, est_gdf_reisenfeld, method = 'spearman'),
                lin_proposed = as.numeric(unlist(DescTools::CCC(truth_gdf, est_gdf)$rho.c[1])),
                lin_reisenfeld = as.numeric(unlist(DescTools::CCC(truth_gdf, est_gdf_reisenfeld)$rho.c[1])))
RESULTS = ROC_SUB[!duplicated(paste0(ROC_SUB$Map_Class, ROC_SUB$Map_Stat)),]
TO_PLOT = rbind(data.frame(Type = 'Proposed', data.frame(RESULTS[,c('Map_Class', 'Map_Stat')], est = RESULTS$cor_proposed, lin = RESULTS$lin_proposed)),
                data.frame(Type = 'Reisenfeld\net al. (2021)',data.frame(RESULTS[,c('Map_Class', 'Map_Stat')], est = RESULTS$cor_reisenfeld, lin = RESULTS$lin_reisenfeld)))
TO_PLOT = TO_PLOT[TO_PLOT$Map_Stat %in% c('Estimated Input', 'Estimated Input, 3x', 'Ideal Input'),]
TO_PLOT = TO_PLOT %>% dplyr::group_by(Map_Class, Map_Stat) %>% dplyr::mutate(best_val = Type[which.max(est)])
p1 = ggplot(TO_PLOT)+
  geom_tile(aes(x=factor(Type), y = factor(Map_Class),  fill =factor(as.numeric(Type == best_val))), color = 'black', alpha = 0.5)+
  geom_text(aes(x=factor(Type), y = factor(Map_Class), label =round(est,3)), color = 'black', size = 3)+
  geom_text(aes(x=factor(Type), y = factor(Map_Class), label =round(est,3)), color = 'black', size = 3, fontface = 'bold', data = TO_PLOT[TO_PLOT$Type == TO_PLOT$best_val,])+
  xlab('')+ylab('')+
  labs(title = 'Correlation between estimated and true GDF under ribbon')+
  scale_fill_manual(name = '', values = c('tomato','olivedrab2'), limits = factor(c(0,1)))+
  guides(color = 'none', fill = 'none')+
  theme(panel.background = element_rect(fill = NA, color = "black"),
        strip.background =element_rect(fill="white", color = "black"))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), limits = c('S3: Varying Profile', 'S2: Spatial Retention', 'S1: Weak Scattering'))+
  facet_grid(.~Map_Stat)+
  theme(axis.ticks = element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())
pdf(file = paste0(current_directory, '/Plots/Correlations.pdf'), width = 7, height = 2)
p1
dev.off()




###########################################
### Figure 7c: Weighted Interval Scores ###
###########################################
WIS_normal = function(x){
  y = x[1]
  mu = x[2]
  sigma2 = x[3]
  alpha_grid = c(0.02, 0.05, seq(0.1,0.9,0.1))
  lowers = mu - qnorm(alpha_grid, mean = 0, sd = 1, lower.tail = FALSE)*sqrt(sigma2)
  uppers = mu + qnorm(alpha_grid, mean = 0, sd = 1, lower.tail = FALSE)*sqrt(sigma2)
  interval_score = (uppers-lowers)+(2/alpha_grid)*(lowers - y)*as.numeric(y<lowers)+
    (2/alpha_grid)*(y - uppers)*as.numeric(y>uppers)
  weighted_interval_score = (1/(length(alpha_grid)+0.5))*0.5*abs(y-mu) +
    (1/(length(alpha_grid)+0.5))*sum(0.5*alpha_grid*interval_score)
  return(weighted_interval_score)
}
RESULTS_RC$wis = apply(cbind(RESULTS_RC$truth_gdf,RESULTS_RC$est_gdf, RESULTS_RC$se_gdf^2),1,FUN = WIS_normal)
RESULTS_RC$wis_reisenfeld = apply(cbind(RESULTS_RC$truth_gdf,RESULTS_RC$est_gdf_reisenfeld, RESULTS_RC$se_gdf_reisenfeld^2),1,FUN = WIS_normal)
WIS_OVERALL = aggregate(cbind(wis, wis_reisenfeld)~Map_Class + Map_Stat, FUN = function(x){mean(x,na.rm=T)}, na.action = na.pass, data = RESULTS_RC[RESULTS_RC$polar_discrete >= RIBBON_LOWER &RESULTS_RC$polar_discrete <= RIBBON_UPPER, ])
TO_PLOT = rbind(data.frame(Type = 'Proposed', data.frame(WIS_OVERALL[,c('Map_Class', 'Map_Stat')], wis = WIS_OVERALL$wis)),
                data.frame(Type = 'Reisenfeld\net al. (2021)',data.frame(WIS_OVERALL[,c('Map_Class', 'Map_Stat')], wis = WIS_OVERALL$wis_reisenfeld)))
TO_PLOT = TO_PLOT[TO_PLOT$Map_Stat %in% c('Estimated Input', 'Estimated Input, 3x', 'Ideal Input'),]
TO_PLOT = TO_PLOT %>% dplyr::group_by(Map_Class, Map_Stat) %>% dplyr::mutate(best_val = Type[which.min(wis)])
p1 = ggplot(TO_PLOT)+
  geom_tile(aes(x=factor(Type), y = factor(Map_Class),  fill =factor(as.numeric(Type == best_val))), color = 'black', alpha = 0.5)+
  geom_text(aes(x=factor(Type), y = factor(Map_Class), label =round(wis*100,2)), color = 'black', size = 3)+
  geom_text(aes(x=factor(Type), y = factor(Map_Class), label =round(wis*100,2)), color = 'black', size = 3, fontface = 'bold', data = TO_PLOT[TO_PLOT$Type == TO_PLOT$best_val,])+
  xlab('')+ylab('')+
  labs(title = '100 x Weighted interval score for estimated GDF')+
  scale_fill_manual(name = '', values = c('tomato','olivedrab2'), limits = factor(c(0,1)))+
  guides(color = 'none', fill = 'none')+
  theme(panel.background = element_rect(fill = NA, color = "black"),
        strip.background =element_rect(fill="white", color = "black"))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), limits = c('S3: Varying Profile', 'S2: Spatial Retention', 'S1: Weak Scattering'))+
  facet_grid(.~Map_Stat)+
  theme(axis.ticks = element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())
pdf(file = paste0(current_directory, '/Plots/WIS.pdf'), width = 7, height = 2)
p1
dev.off()




############################################
### Figure 7d: Coverage of 95% Intervals ###
############################################
RESULTS_RC$lowers = RESULTS_RC$est_gdf - 1.96*RESULTS_RC$se_gdf
RESULTS_RC$uppers = RESULTS_RC$est_gdf + 1.96*RESULTS_RC$se_gdf
RESULTS_RC$covers = as.numeric(RESULTS_RC$truth_gdf <= RESULTS_RC$uppers & RESULTS_RC$truth_gdf >= RESULTS_RC$lowers)
ERROR_OVERALL = aggregate(cbind(covers)~Map_Class + Map_Stat, FUN = function(x){mean(x,na.rm=T)}, na.action = na.pass, data = RESULTS_RC[RESULTS_RC$polar_discrete >= RIBBON_LOWER &RESULTS_RC$polar_discrete <= RIBBON_UPPER, ])
RESULTS_RC$lowers = RESULTS_RC$est_gdf_reisenfeld - 1.96*RESULTS_RC$se_gdf_reisenfeld
RESULTS_RC$uppers = RESULTS_RC$est_gdf_reisenfeld + 1.96*RESULTS_RC$se_gdf_reisenfeld
RESULTS_RC$covers = as.numeric(RESULTS_RC$truth_gdf <= RESULTS_RC$uppers & RESULTS_RC$truth_gdf >= RESULTS_RC$lowers)
ERROR_OVERALL_reisenfeld = aggregate(cbind(covers)~Map_Class + Map_Stat, FUN = function(x){mean(x,na.rm=T)}, na.action = na.pass, data = RESULTS_RC[RESULTS_RC$polar_discrete >= RIBBON_LOWER &RESULTS_RC$polar_discrete <= RIBBON_UPPER, ])
TO_PLOT = rbind(data.frame(Type = 'Proposed', ERROR_OVERALL),
                data.frame(Type = 'Reisenfeld\net al. (2021)',ERROR_OVERALL_reisenfeld ))
TO_PLOT = TO_PLOT[TO_PLOT$Map_Stat %in% c('Estimated Input', 'Estimated Input, 3x', 'Ideal Input'),]
TO_PLOT = TO_PLOT %>% dplyr::group_by(Map_Class, Map_Stat) %>% dplyr::mutate(best_val = Type[which.min(0.95 - covers)])
p1 = ggplot(TO_PLOT)+
  geom_tile(aes(x=factor(Type), y = factor(Map_Class),  fill =factor(as.numeric(Type == best_val))), color = 'black', alpha = 0.5)+
  geom_text(aes(x=factor(Type), y = factor(Map_Class), label =round(covers,2)), color = 'black', size = 3)+
  geom_text(aes(x=factor(Type), y = factor(Map_Class), label =round(covers,2)), color = 'black', size = 3, fontface = 'bold', data = TO_PLOT[TO_PLOT$Type == TO_PLOT$best_val,])+
  xlab('')+ylab('')+
  labs(title = 'Coverage of 95% Confidence Intervals')+
  scale_fill_manual(name = '', values = c('tomato','olivedrab2'), limits = factor(c(0,1)))+
  guides(color = 'none', fill = 'none')+
  theme(panel.background = element_rect(fill = NA, color = "black"),
        strip.background =element_rect(fill="white", color = "black"))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), limits = c('S3: Varying Profile', 'S2: Spatial Retention', 'S1: Weak Scattering'))+
  facet_grid(.~Map_Stat)+
  theme(axis.ticks = element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())
pdf(file = paste0(current_directory, '/Plots/Coverage.pdf'), width = 7, height = 2)
p1
dev.off()



#################################
### Supp. Figure E5: Skewness ###
#################################
estSkew = function(x,y){
  spline_func = splines::periodicSpline(x,y,period = 360, ord = 4L)
  x2= seq(min(x),max(x))
  y2 = predict(spline_func,x2)$y
  y2[y2<0]=0
  k = Weighted.Desc.Stat::w.skewness(x2, y2/sum(y2,na.rm=T))
  return(k)
}
RESULTS_RC = RESULTS_RC[order(RESULTS_RC$azimuthal_discrete, RESULTS_RC$polar_discrete),]
RESULTS_RC = RESULTS_RC %>% dplyr::group_by(Map_Stat, Map_Class, azimuthal_discrete, ESA) %>% 
  dplyr::mutate(sum_val = sum(est_ribbon,na.rm=T),
                sum_val_reisenfeld = sum(est_ribbon_reisenfeld,na.rm=T),
                sum_val_truth = sum(truth_ribbon,na.rm=T))
RESULTS_RC = RESULTS_RC %>% dplyr::group_by(Map_Stat, Map_Class, azimuthal_discrete, ESA) %>% 
  dplyr::mutate(peak_loc = polar_discrete[which.max(est_ribbon)],
                peak_loc_reisenfeld = polar_discrete[which.max(est_ribbon_reisenfeld)],
                peak_loc_truth = polar_discrete[which.max(truth_ribbon)])
DAT_SUB = RESULTS_RC
DAT_SUB = DAT_SUB[DAT_SUB$ESA %in% c('esa4', 'noesa'),]
DAT_SUB =  DAT_SUB %>% dplyr::group_by(Map_Stat, Map_Class, azimuthal_discrete, ESA) %>%
  dplyr::mutate(est_ribbon0= c(rev(cummin(rev(est_ribbon[polar_discrete <= peak_loc]))),cummin(est_ribbon[polar_discrete > peak_loc])),
                est_ribbon_reisenfeld0= c(rev(cummin(rev(est_ribbon_reisenfeld[polar_discrete <= peak_loc]))),cummin(est_ribbon_reisenfeld[polar_discrete > peak_loc])),
                truth_ribbon0= c(rev(cummin(rev(truth_ribbon[polar_discrete <= peak_loc]))),cummin(truth_ribbon[polar_discrete > peak_loc])))
DAT_SUB$est_ribbon = DAT_SUB$est_ribbon0
DAT_SUB$est_ribbon_reisenfeld = DAT_SUB$est_ribbon_reisenfeld0
DAT_SUB$truth_ribbon = DAT_SUB$truth_ribbon0
DAT_SUB = DAT_SUB %>% dplyr::group_by(Map_Stat, Map_Class, azimuthal_discrete, ESA) %>% 
  dplyr::mutate(k_skew = estSkew(polar_discrete, est_ribbon),
                k_skew_reisenfeld = estSkew(polar_discrete, est_ribbon_reisenfeld),
                k_skew_truth = estSkew(polar_discrete, truth_ribbon))
DAT_SUB = DAT_SUB[!duplicated(paste0(DAT_SUB$azimuthal_discrete, DAT_SUB$Map_Class, DAT_SUB$Map_Stat, DAT_SUB$ESA)),]
TO_PLOT = rbind(data.frame(Type = 'Proposed', data.frame(DAT_SUB[,c('Map_Class', 'Map_Stat', 'azimuthal_discrete')], skew = DAT_SUB$k_skew, bias = abs(DAT_SUB$k_skew - DAT_SUB$k_skew_truth))),
                data.frame(Type = 'Reisenfeld\net al. (2021)',data.frame(DAT_SUB[,c('Map_Class', 'Map_Stat',  'azimuthal_discrete')], skew = DAT_SUB$k_skew_reisenfeld, bias = abs(DAT_SUB$k_skew_reisenfeld - DAT_SUB$k_skew_truth))),
                data.frame(Type = 'Truth',data.frame(DAT_SUB[,c('Map_Class', 'Map_Stat',  'azimuthal_discrete')], skew = DAT_SUB$k_skew_truth, bias = abs(DAT_SUB$k_skew_truth - DAT_SUB$k_skew_truth))))
p0=ggplot(TO_PLOT[TO_PLOT$Map_Stat %in% c('Estimated Input','Estimated Input, 3x', 'Ideal Input') & !is.na(TO_PLOT$skew),])+
  geom_hline(yintercept = 0, linetype = 2, color = 'gray', size = 1.5)+
  geom_boxplot(aes(x=factor(Type), y =skew, fill = factor(Type)),color = 'black',width = 0.8, position = position_dodge2(0.85, preserve = "single"), outlier.shape = NA)+
  geom_point(aes(x=factor(Type), y =skew, fill = factor(Type)),shape = 21,color = 'black', position = position_jitterdodge(jitter.width = 0.9, jitter.height = 0), alpha = 0.3)+
  facet_grid(Map_Class~Map_Stat)+
  xlab('Simulation Scenario')+ylab('Ribbon Skewness')+
  labs(title = 'Estimated Skewness per Azimuthal Angle')+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  guides(color = 'none', fill = 'none')+
  theme(legend.position = 'top')+
  coord_cartesian(ylim = c(-1.5,0.5))
pdf(file = paste0(current_directory, '/Plots/Skewness.pdf'), width = 8, height = 5)
p0
dev.off()



###############################################
### Supp. Figure E6: Full Width at Half Max ###
###############################################
estFWHM = function(x,y){
  spline_func = splines::periodicSpline(x,y,period = 360, ord = 4L)
  x2= seq(min(x),max(x),0.5)
  y2 = predict(spline_func,x2)$y
  y2[y2<0]=0
  x2 = x2[y2 >= 0.5]
  return(max(x2)-min(x2))
}
RESULTS_RC = RESULTS_RC[order(RESULTS_RC$azimuthal_discrete, RESULTS_RC$polar_discrete),]
DAT_SUB = RESULTS_RC
DAT_SUB = DAT_SUB[DAT_SUB$ESA %in% c('esa4', 'noesa'),]
DAT_SUB =  DAT_SUB %>% dplyr::group_by(Map_Stat, Map_Class, azimuthal_discrete, ESA) %>%
  dplyr::mutate(est_ribbon0= c(rev(cummin(rev(est_ribbon[polar_discrete <= peak_loc]))),cummin(est_ribbon[polar_discrete > peak_loc])),
                est_ribbon_reisenfeld0= c(rev(cummin(rev(est_ribbon_reisenfeld[polar_discrete <= peak_loc]))),cummin(est_ribbon_reisenfeld[polar_discrete > peak_loc])),
                truth_ribbon0= c(rev(cummin(rev(truth_ribbon[polar_discrete <= peak_loc]))),cummin(truth_ribbon[polar_discrete > peak_loc])))
DAT_SUB$est_ribbon = DAT_SUB$est_ribbon0
DAT_SUB$est_ribbon_reisenfeld = DAT_SUB$est_ribbon_reisenfeld0
DAT_SUB$truth_ribbon = DAT_SUB$truth_ribbon0
DAT_SUB = DAT_SUB %>% dplyr::group_by(Map_Stat, Map_Class, azimuthal_discrete, ESA) %>% 
  dplyr::mutate(max_val = max(est_ribbon,na.rm=T),
                max_val_reisenfeld = max(est_ribbon_reisenfeld,na.rm=T),
                max_val_truth = max(truth_ribbon,na.rm=T))
DAT_SUB$k = DAT_SUB$est_ribbon/DAT_SUB$max_val
DAT_SUB$k_reisenfeld = DAT_SUB$est_ribbon_reisenfeld/DAT_SUB$max_val_reisenfeld
DAT_SUB$k_truth = DAT_SUB$truth_ribbon/DAT_SUB$max_val_truth
DAT_SUB = DAT_SUB %>% dplyr::group_by(Map_Stat, Map_Class, azimuthal_discrete, ESA) %>% 
  dplyr::mutate(k_skew = estFWHM(polar_discrete, k),
                k_skew_reisenfeld = estFWHM(polar_discrete, k_reisenfeld),
                k_skew_truth = estFWHM(polar_discrete, k_truth))
DAT_SUB = DAT_SUB[!duplicated(paste0(DAT_SUB$azimuthal_discrete, DAT_SUB$Map_Class, DAT_SUB$Map_Stat, DAT_SUB$ESA)),]
TO_PLOT = rbind(data.frame(Type = 'Proposed', data.frame(DAT_SUB[,c('Map_Class', 'Map_Stat', 'azimuthal_discrete')], skew = DAT_SUB$k_skew, bias = abs(DAT_SUB$k_skew - DAT_SUB$k_skew_truth))),
                data.frame(Type = 'Reisenfeld\net al. (2021)',data.frame(DAT_SUB[,c('Map_Class', 'Map_Stat', 'azimuthal_discrete')], skew = DAT_SUB$k_skew_reisenfeld, bias = abs(DAT_SUB$k_skew_reisenfeld - DAT_SUB$k_skew_truth))),
                data.frame(Type = 'Truth',data.frame(DAT_SUB[,c('Map_Class', 'Map_Stat', 'azimuthal_discrete')], skew = DAT_SUB$k_skew_truth, bias = abs(DAT_SUB$k_skew_truth - DAT_SUB$k_skew_truth))))
p0=ggplot(TO_PLOT[TO_PLOT$Map_Stat %in% c('Estimated Input','Estimated Input, 3x', 'Ideal Input') & !is.na(TO_PLOT$skew),])+
  geom_boxplot(aes(x=factor(Type), y =skew, fill = factor(Type)),color = 'black',width = 0.8, position = position_dodge2(0.85, preserve = "single"), outlier.shape = NA)+
  geom_point(aes(x=factor(Type), y =skew, fill = factor(Type)),shape = 21,color = 'black', position = position_jitterdodge(0.9), alpha = 0.3)+
  facet_grid(Map_Class~Map_Stat)+
  xlab('Simulation Scenario')+ylab('Ribbon Width (FWHM)')+
  labs(title = 'Estimated Width per Azimuthal Angle')+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  guides(color = 'none', fill = 'none')+
  theme(legend.position = 'top')+
  coord_cartesian(ylim = c(0,25))
pdf(file = paste0(current_directory, '/Plots/FWHM.pdf'), width = 8, height = 5)
p0
dev.off()







########################################
### Figure 8a: Ribbon Cross-Sections ###
########################################
RESULTS_RC = RESULTS_RC[order(RESULTS_RC$azimuthal_discrete, RESULTS_RC$polar_discrete),]
TO_PLOT = rbind(data.frame(Type = 'Proposed', data.frame(RESULTS_RC[,c('Map_Class','Map_Stat', 'azimuthal_discrete','polar_discrete')], ribbon = RESULTS_RC$est_ribbon, residual = RESULTS_RC$est_gdf )),
                data.frame(Type = 'Reisenfeld et al. (2021)', data.frame(RESULTS_RC[,c('Map_Class','Map_Stat', 'azimuthal_discrete','polar_discrete')], ribbon = RESULTS_RC$est_ribbon_reisenfeld, residual = RESULTS_RC$est_gdf_reisenfeld )),
                data.frame(Type = 'Truth', data.frame(RESULTS_RC[,c('Map_Class','Map_Stat', 'azimuthal_discrete','polar_discrete')], ribbon = RESULTS_RC$truth_ribbon,  residual = RESULTS_RC$truth_gdf)))
p2 = ggplot(TO_PLOT[TO_PLOT$Map_Stat == 'Ideal Input',])+
  geom_hline(yintercept = seq(0.05,0.1,0.05), color = 'gray', size = 0.2, linetype = 2)+
  geom_line(aes(x=polar_discrete, y=ribbon, group = factor(azimuthal_discrete), color = azimuthal_discrete ), size = 0.7, alpha = 0.7)+
  theme_classic()+
  xlab('Polar angle')+ylab('ENA rate')+
  theme(legend.position = 'top')+
  facet_grid(Map_Class~Type, scales = 'free_y')+
  coord_cartesian(xlim = c(75,140))+
  scale_color_gradientn('Azimuthal angle',colors = met.brewer("Johnson"))
ggsave(p2, file = paste0(current_directory, '/Plots/Ribbon_CrossSections.pdf'), width = 10, height = 6)


#####################################
### Figure 8b: GDF Cross-Sections ###
#####################################
TO_PLOT = rbind(data.frame(Type = 'Proposed', data.frame(RESULTS_RC[,c('Map_Class','Map_Stat', 'azimuthal_discrete','polar_discrete')], ribbon = RESULTS_RC$est_ribbon, residual = RESULTS_RC$est_gdf )),
                data.frame(Type = 'Reisenfeld et al. (2021)', data.frame(RESULTS_RC[,c('Map_Class','Map_Stat', 'azimuthal_discrete','polar_discrete')], ribbon = RESULTS_RC$est_ribbon_reisenfeld, residual = RESULTS_RC$est_gdf_reisenfeld )),
                data.frame(Type = 'Truth', data.frame(RESULTS_RC[,c('Map_Class','Map_Stat', 'azimuthal_discrete','polar_discrete')], ribbon = RESULTS_RC$truth_ribbon,  residual = RESULTS_RC$truth_gdf)))
p3 = ggplot(TO_PLOT[TO_PLOT$Map_Stat == 'Ideal Input',])+
  geom_hline(yintercept = seq(0.05,0.1,0.05), color = 'gray', size = 0.2, linetype = 2)+
  geom_line(aes(x=polar_discrete, y=residual, group = factor(azimuthal_discrete), color = azimuthal_discrete ), size = 0.7, alpha = 0.7)+
  theme_classic()+
  xlab('Polar angle')+ylab('ENA rate')+
  theme(legend.position = 'top')+
  facet_grid(Map_Class~Type, scales = 'free_y')+
  coord_cartesian(xlim = c(75,140))+
  scale_color_gradientn('Azimuthal angle',colors = met.brewer("Johnson"))
ggsave(p3, file = paste0(current_directory, '/Plots/GDF_CrossSections.pdf'), width = 10, height = 6)

