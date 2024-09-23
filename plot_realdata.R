###########################################
###########################################
### Plotting Selected Real Data Results ###
###########################################
###########################################

current_directory = dirname(rstudioapi::getSourceEditorContext()$path)

#########################
### Read in Libraries ###
#########################

library(ggplot2)
library(scales)
library(gridtext)
library(cowplot)
library(MetBrewer)
library(ggpubr)
library(ggforce)
library(colorspace)
library(pracma)
source(paste0(current_directory,'/Internal_Functions/Coordinate_Transform_Functions.R'))
source(paste0(current_directory,'/Internal_Functions/Ribbon_Center_Functions.R'))
### IBEX color palette
ibex_palette = read.csv(paste0(current_directory,'/ibex_palette.csv'))

#############################
### Set Up File Structure ###
#############################

if(!dir.exists(paste0(current_directory,'/Plots'))){
  dir.create(paste0(current_directory,'/Plots'))
}
if(!dir.exists(paste0(current_directory,'/Plots/RealData'))){
  dir.create(paste0(current_directory,'/Plots/RealData'))
}


#############################
### Unzip files as needed ###
#############################
library(GEOquery)
if(!file.exists(paste0(current_directory, '/output_theseus_ecliptic_part1.csv'))){
  gunzip(filename =paste0(current_directory, '/output_theseus_ecliptic_part1.gz'), destname = paste0(current_directory, '/output_theseus_ecliptic_part1.csv'))
}
if(!file.exists(paste0(current_directory, '/output_theseus_ecliptic_part2.csv'))){
  gunzip(filename =paste0(current_directory, '/output_theseus_ecliptic_part2.gz'), destname = paste0(current_directory, '/output_theseus_ecliptic_part2.csv'))
}
if(!file.exists(paste0(current_directory, '/output_theseus_ecliptic_reisenfeld.csv'))){
  gunzip(filename =paste0(current_directory, '/output_theseus_ecliptic_reisenfeld.gz'), destname = paste0(current_directory, '/output_theseus_ecliptic_reisenfeld.csv'))
}
if(!file.exists(paste0(current_directory, '/output_theseus_rc.csv'))){
  gunzip(filename =paste0(current_directory, '/output_theseus_rc.gz'), destname = paste0(current_directory, '/output_theseus_rc.csv'))
}
if(!file.exists(paste0(current_directory, '/output_theseus_rc_reisenfeld.csv'))){
  gunzip(filename =paste0(current_directory, '/output_theseus_rc_reisenfeld.gz'), destname = paste0(current_directory, '/output_theseus_rc_reisenfeld.csv'))
}


#################################
### Read in and Merge Results ###
#################################

THESEUS_ECLIPTIC1 = read.csv(paste0(current_directory, '/output_theseus_ecliptic_part1.csv'))
THESEUS_ECLIPTIC2 = read.csv(paste0(current_directory, '/output_theseus_ecliptic_part2.csv'))
THESEUS_ECLIPTIC = rbind(THESEUS_ECLIPTIC1, THESEUS_ECLIPTIC2)

THESEUS_RC = data.frame(data.table::fread(file = paste0(current_directory,'/output_theseus_rc.csv')))
THESEUS_ECLIPTIC$input_data = THESEUS_ECLIPTIC$est_ribbon + THESEUS_ECLIPTIC$est_gdf
THESEUS_RC$input_data = THESEUS_RC$est_ribbon + THESEUS_RC$est_gdf

ISOC_ECLIPTIC = data.frame(data.table::fread(file = paste0(current_directory,'/output_isoc_ecliptic.csv')))
ISOC_RC = data.frame(data.table::fread(file = paste0(current_directory,'/output_isoc_rc.csv')))

### Merge Reisenfeld et al. (2021) Results
REISENFELD_ECLIPTIC = data.frame(data.table::fread(file = paste0(current_directory,'/output_theseus_ecliptic_reisenfeld.csv')))
THESEUS_ECLIPTIC = merge(THESEUS_ECLIPTIC, REISENFELD_ECLIPTIC, by = c('lon', 'lat', 'Time_Group', 'ESA'), all.x = T, all.y = T)
REISENFELD_RC = data.frame(data.table::fread(file = paste0(current_directory,'/output_theseus_rc_reisenfeld.csv')))
THESEUS_RC = merge(THESEUS_RC, REISENFELD_RC, by = c('lon', 'lat', 'Time_Group', 'ESA'), all.x = T, all.y = T)

REISENFELD_ECLIPTIC = data.frame(data.table::fread(file = paste0(current_directory,'/output_isoc_ecliptic_reisenfeld.csv')))
ISOC_ECLIPTIC = merge(ISOC_ECLIPTIC, REISENFELD_ECLIPTIC, by = c('lon', 'lat', 'Time_Group', 'ESA'), all.x = T, all.y = T)
REISENFELD_RC = data.frame(data.table::fread(file = paste0(current_directory,'/output_isoc_rc_reisenfeld.csv')))
ISOC_RC = merge(ISOC_RC, REISENFELD_RC, by = c('lon', 'lat', 'Time_Group', 'ESA'), all.x = T, all.y = T)

THESEUS_RC$azimuthal_discrete = THESEUS_RC$lon
THESEUS_RC$polar_discrete = THESEUS_RC$lat+90
ISOC_RC$azimuthal_discrete = ISOC_RC$lon
ISOC_RC$polar_discrete = ISOC_RC$lat+90


### Convert ISOC Flux Maps to ENA Rates
Gj_triples = c(0.00013, 0.00037, 0.00073, 0.0014, 0.0025, 0.0042)
Gj_doubles = c(0.00031, 0.00085, 0.0015, 0.0029, 0.0045, 0.0071)
Ej = c(0.45, 0.71, 1.10, 1.74, 2.73, 4.29)
ISOC_ECLIPTIC$est_gdf = ISOC_ECLIPTIC$est_gdf*Gj_triples[ISOC_ECLIPTIC$ESA]*Ej[ISOC_ECLIPTIC$ESA]
ISOC_ECLIPTIC$est_ribbon = ISOC_ECLIPTIC$est_ribbon*Gj_triples[ISOC_ECLIPTIC$ESA]*Ej[ISOC_ECLIPTIC$ESA]
ISOC_ECLIPTIC$se_gdf = ISOC_ECLIPTIC$se_gdf*Gj_triples[ISOC_ECLIPTIC$ESA]*Ej[ISOC_ECLIPTIC$ESA]
ISOC_ECLIPTIC$se_ribbon = ISOC_ECLIPTIC$se_ribbon*Gj_triples[ISOC_ECLIPTIC$ESA]*Ej[ISOC_ECLIPTIC$ESA]
ISOC_ECLIPTIC$input_data = ISOC_ECLIPTIC$est_gdf + ISOC_ECLIPTIC$est_ribbon


ISOC_RC$est_gdf = ISOC_RC$est_gdf*Gj_triples[ISOC_RC$ESA]*Ej[ISOC_RC$ESA]
ISOC_RC$est_ribbon = ISOC_RC$est_ribbon*Gj_triples[ISOC_RC$ESA]*Ej[ISOC_RC$ESA]
ISOC_RC$se_gdf = ISOC_RC$se_gdf*Gj_triples[ISOC_RC$ESA]*Ej[ISOC_RC$ESA]
ISOC_RC$se_ribbon = ISOC_RC$se_ribbon*Gj_triples[ISOC_RC$ESA]*Ej[ISOC_RC$ESA]


###########################################
###########################################
### Plots in Ecliptic Coordinate System ###
###########################################
###########################################

### Standard IBEX plots in ecliptic coordinates show the heliospheric 
### nose as the approx. central longitude, 
### with longitude increasing from right to left
center = 180-(360 - 255)
orig_360 = seq(0,359,85)#[-1]
new_360 = orig_360-center
new_360[new_360<0]=new_360[new_360<0]+360
THESEUS_ECLIPTIC$ecliptic_lon_center = THESEUS_ECLIPTIC$lon - center + 0.01
THESEUS_ECLIPTIC$ecliptic_lon_center[THESEUS_ECLIPTIC$ecliptic_lon_center<0.01]=THESEUS_ECLIPTIC$ecliptic_lon_center[THESEUS_ECLIPTIC$ecliptic_lon_center<0.01]+360
THESEUS_ECLIPTIC$ecliptic_lon_center = THESEUS_ECLIPTIC$ecliptic_lon_center + 1
ISOC_ECLIPTIC$ecliptic_lon_center = ISOC_ECLIPTIC$lon - center + 0.01
ISOC_ECLIPTIC$ecliptic_lon_center[ISOC_ECLIPTIC$ecliptic_lon_center<0.01]=ISOC_ECLIPTIC$ecliptic_lon_center[ISOC_ECLIPTIC$ecliptic_lon_center<0.01]+360
ISOC_ECLIPTIC$ecliptic_lon_center = ISOC_ECLIPTIC$ecliptic_lon_center + 3


############################
### Annoted ENA Rate Map ### 
############################
THESEUS_SUB = THESEUS_ECLIPTIC[THESEUS_ECLIPTIC$Time_Group == '2010A' & THESEUS_ECLIPTIC$ESA == 5,c('ecliptic_lon_center', 'lat', 'input_data')]
MAP_REV = Coord_Transform_Points(MAP =expand.grid(lon = seq(0,360,1), lat = c(-10,35)), RIBBON_CENTER = c(221,39),REVERSE = T)
MAP_REV$ecliptic_lon_center = MAP_REV$lon_new - center + 0.01
MAP_REV$ecliptic_lon_center[MAP_REV$ecliptic_lon_center<0.01]=MAP_REV$ecliptic_lon_center[MAP_REV$ecliptic_lon_center<0.01]+360
MAP_REV$ecliptic_lon_center = MAP_REV$ecliptic_lon_center + 3
MAP_REV = MAP_REV[order(MAP_REV$ecliptic_lon_center),]
p1 = ggplot(THESEUS_SUB)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = input_data, color = input_data))+
  theme_classic()+
  xlab('Longitude')+ylab('Latitude')+
  theme(legend.position = 'bottom', legend.key.size = unit(1.2,'cm'), axis.line.x = element_blank())+
  labs(title = '')+
  scale_color_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black' )+
  scale_fill_gradientn ('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  geom_circle(aes(x0=180, y0=-12, r=65), inherit.aes=FALSE, color = 'gray', linetype = 1, size = 6, n = 40, fill = 'lightgray', alpha = 0.002)+
  geom_ellipse(aes(x0=359, y0=0, a=50,b=90, angle = 0), inherit.aes=FALSE, color = 'gray', linetype = 1, size = 6, n = 40, fill = 'lightgray', alpha = 0.002)+
  geom_ellipse(aes(x0=1, y0=0, a=50,b=90, angle = 0), inherit.aes=FALSE, color = 'gray', linetype = 1, size = 6, n = 40, fill = 'lightgray', alpha = 0.002)+
  annotate("text",x=170, y = 20, label = 'Nose\nRegion', color = 'darkgray', size = 4, hjust = 0.5)+
  annotate("text",x=25, y = -10, label = 'Tail\nRegion', color = 'darkgray', size = 4, hjust = 0.5)+
  annotate("text",x=335, y = -10, label = 'Tail\nRegion', color = 'darkgray', size = 4, hjust = 0.5)+
  annotate("text",x=305, y = 70, label = 'Ribbon', color = 'black', size = 4, hjust = 0.5, angle = 0)+ 
  annotate("text",x=70, y = 10, label = 'Ribbon', color = 'black', size = 4, hjust = 0.5, angle = 45)+ 
  annotate("text",x=180, y = -30, label = 'Ribbon', color = 'black', size = 4, hjust = 0.5, angle = -20)+ 
  geom_point(aes(x=ecliptic_lon_center, y =lat_new), data = MAP_REV, color = 'black', size = 0.5)+
  coord_fixed(ratio = 1, clip = 'on', xlim=c(360,0), ylim = c(-89,89))+
  guides(color = 'none', fill = 'none')
pdf(file = paste0(current_directory, '/Plots/RealData/Annotated_2010.pdf'), width = 6, height = 4)
p1
dev.off()





###################################
### Comparison of Two 2010 Maps ### (without ribbon annotations)
###################################
ISOC_SUB = ISOC_ECLIPTIC[ISOC_ECLIPTIC$Time_Group == '2010' & ISOC_ECLIPTIC$ESA == 5,c('ecliptic_lon_center', 'lat', 'input_data')]
THESEUS_SUB = THESEUS_ECLIPTIC[THESEUS_ECLIPTIC$Time_Group == '2010A' & THESEUS_ECLIPTIC$ESA == 5,c('ecliptic_lon_center', 'lat', 'input_data')]
# LIMITS = c(0,max(c(THESEUS_SUB$input_data, ISOC_SUB$input_data), na.rm=T))
LIMITS = c(0,0.35)
p0 = ggplot(ISOC_SUB)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = input_data, color = input_data))+
  theme_classic()+
  xlab('Longitude')+ylab('Latitude')+
  labs(title = 'ISOC')+
  #theme(legend.position = 'bottom', legend.key.size = unit(1.2,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black',limits = LIMITS,oob = scales::oob_squish )+
  scale_fill_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black' ,limits = LIMITS,oob = scales::oob_squish )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  coord_fixed(ratio = 1, clip = 'on', xlim=c(360,1), ylim = c(-89,89))+
  theme(legend.position = 'none')
p1 = ggplot(THESEUS_SUB)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = input_data, color = input_data))+
  theme_classic()+
  xlab('Longitude')+ylab('Latitude')+
  #theme(legend.position = 'bottom', legend.key.size = unit(1.2,'cm'), axis.line.x = element_blank())+
  labs(title = 'Theseus')+
  scale_color_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black',limits = LIMITS ,oob = scales::oob_squish )+
  scale_fill_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black' ,limits = LIMITS,oob = scales::oob_squish )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  coord_fixed(ratio = 1, clip = 'on', xlim=c(360,0), ylim = c(-89,89))+
  theme(legend.position = 'right')
p3 = ggpubr::ggarrange(p0,p1, ncol = 2, common.legend=F, heights = c(1,1), widths = c(1,1.22))
pdf(file = paste0(current_directory, '/Plots/RealData/Compare_2010.pdf'), width = 10, height = 5)
p3
dev.off()




#################################
### Theseus ESA 4 Separations ### 
#################################
THESEUS_SUB = THESEUS_ECLIPTIC[THESEUS_ECLIPTIC$ESA == 4,c('ecliptic_lon_center', 'lat', 'est_gdf','est_ribbon','Time_Group')]
p1 = ggplot(THESEUS_SUB)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = est_gdf, color = est_gdf))+
  theme_classic()+
  xlab('Longitude')+ylab('Latitude')+
  labs(title = 'GDF Map')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENAs/sec',colors = ibex_palette$hex , na.value = 'black' )+
  scale_fill_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  facet_grid(gsub('A','',Time_Group)~.)+
  coord_fixed(ratio = 1, clip = 'on', xlim=c(360,1), ylim = c(-89,89))
p2 = ggplot(THESEUS_SUB)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = est_ribbon, color = est_ribbon))+
  theme_classic()+
  xlab('Longitude')+ylab('Latitude')+
  labs(title = 'Ribbon Map')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENAs/sec',colors = ibex_palette$hex , na.value = 'black' )+
  scale_fill_gradientn('ENAs/sec',colors = ibex_palette$hex , na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  facet_grid(gsub('A','',Time_Group)~.)+
  coord_fixed(ratio = 1, clip = 'on', xlim=c(360,1), ylim = c(-89,89))
p0 = ggplot(THESEUS_SUB)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = est_gdf+est_ribbon, color = est_gdf+est_ribbon))+
  theme_classic()+
  xlab('Longitude')+ylab('Latitude')+
  labs(title = 'Total Map')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENAs/sec',colors = ibex_palette$hex , na.value = 'black' )+
  scale_fill_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  facet_grid(gsub('A','',Time_Group)~.)+
  coord_fixed(ratio = 1, clip = 'on', xlim=c(360,1), ylim = c(-89,89))
p3 = ggpubr::ggarrange(p0,p1,p2, ncol = 3, common.legend=F)
ggsave(p3,filename = paste0(current_directory, '/Plots/RealData/Theseus_ESA4.pdf'), width = 11, height = 16)





##############################
### ISOC ESA 4 Separations ### 
##############################
ISOC_SUB = ISOC_ECLIPTIC[ISOC_ECLIPTIC$ESA == 4,c('ecliptic_lon_center', 'lat', 'est_gdf','est_ribbon','Time_Group')]
p1 = ggplot(ISOC_SUB)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = est_gdf, color = est_gdf))+
  theme_classic()+
  xlab('Longitude')+ylab('Latitude')+
  labs(title = 'GDF Map')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENAs/sec',colors = ibex_palette$hex , na.value = 'black' )+
  scale_fill_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  facet_grid(gsub('A','',Time_Group)~.)+
  coord_fixed(ratio = 1, clip = 'on', xlim=c(360,1), ylim = c(-89,89))
p2 = ggplot(ISOC_SUB)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = est_ribbon, color = est_ribbon))+
  theme_classic()+
  xlab('Longitude')+ylab('Latitude')+
  labs(title = 'Ribbon Map')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENAs/sec',colors = ibex_palette$hex , na.value = 'black' )+
  scale_fill_gradientn('ENAs/sec',colors = ibex_palette$hex , na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  facet_grid(gsub('A','',Time_Group)~.)+
  coord_fixed(ratio = 1, clip = 'on', xlim=c(360,1), ylim = c(-89,89))
p0 = ggplot(ISOC_SUB)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = est_gdf+est_ribbon, color = est_gdf+est_ribbon))+
  theme_classic()+
  xlab('Longitude')+ylab('Latitude')+
  labs(title = 'Total Map')+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'), axis.line.x = element_blank())+
  scale_color_gradientn('ENAs/sec',colors = ibex_palette$hex , na.value = 'black' )+
  scale_fill_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black' )+
  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  facet_grid(gsub('A','',Time_Group)~.)+
  coord_fixed(ratio = 1, clip = 'on', xlim=c(360,1), ylim = c(-89,89))
p3 = ggpubr::ggarrange(p0,p1,p2, ncol = 3, common.legend=F)
ggsave(p3,filename = paste0(current_directory, '/Plots/RealData/ISOC_ESA4.pdf'), width = 11, height = 16)








###############################
### ISOC Yearly Separations ### 
###############################
for(time in unique(ISOC_ECLIPTIC$Time_Group)){
  ISOC_SUB = ISOC_ECLIPTIC[ISOC_ECLIPTIC$Time_Group == time,c('ecliptic_lon_center', 'lat', 'est_gdf','est_ribbon','ESA')]
  TO_PLOT = rbind(data.frame(ISOC_SUB[,c('lat', 'ecliptic_lon_center', 'ESA')], ena_rate = ISOC_SUB$est_ribbon+ISOC_SUB$est_gdf, Type = 'Total'),
                  data.frame(ISOC_SUB[,c('lat', 'ecliptic_lon_center', 'ESA')], ena_rate = ISOC_SUB$est_ribbon, Type = 'Ribbon'),
                  data.frame(ISOC_SUB[,c('lat', 'ecliptic_lon_center', 'ESA')], ena_rate = ISOC_SUB$est_gdf, Type = 'GDF'))
  TO_PLOT$Type = factor(TO_PLOT$Type, levels = c('Total', 'Ribbon', 'GDF'))
  p0 = replicate(5,list(NULL))
  for(i in 1:5){
    p0[[i]] = ggplot(TO_PLOT[TO_PLOT$ESA == i+1,])+
      geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = ena_rate, color = ena_rate))+
      theme_classic()+
      xlab('Longitude')+ylab('Latitude')+
      theme(legend.position = 'right', legend.key.size = unit(0.5,'cm'), axis.line.x = element_blank())+
      scale_color_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black' )+
      scale_fill_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black' )+
      scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
      scale_y_continuous(expand = c(0,0))+
      theme(panel.background = element_rect(fill = NA, color = "black"))+
      facet_grid(paste0('ESA ',ESA)~Type)+
      coord_fixed(ratio = 1, clip = 'on', xlim=c(360,1), ylim = c(-89,89))
  }
  p1 = ggpubr::ggarrange(p0[[1]],p0[[2]],p0[[3]],p0[[4]],p0[[5]], ncol=1, nrow =5, common.legend = F)
  ggsave(p1,filename = paste0(current_directory, '/Plots/RealData/ISOC_',time,'.pdf'), width = 10, height = 10)
}







##################################
### Theseus Yearly Separations ###  (data for first 6 months of each year)
##################################
for(time in unique(THESEUS_ECLIPTIC$Time_Group)){
  THESEUS_SUB = THESEUS_ECLIPTIC[THESEUS_ECLIPTIC$Time_Group == time,c('ecliptic_lon_center', 'lat', 'est_gdf','est_ribbon','input_data','ESA')]
  TO_PLOT = rbind(data.frame(THESEUS_SUB[,c('lat', 'ecliptic_lon_center', 'ESA')], ena_rate = THESEUS_SUB$input_data, Type = 'Total'),
                  data.frame(THESEUS_SUB[,c('lat', 'ecliptic_lon_center', 'ESA')], ena_rate = THESEUS_SUB$est_ribbon, Type = 'Ribbon'),
                  data.frame(THESEUS_SUB[,c('lat', 'ecliptic_lon_center', 'ESA')], ena_rate = THESEUS_SUB$est_gdf, Type = 'GDF'))
  TO_PLOT$Type = factor(TO_PLOT$Type, levels = c('Total', 'Ribbon', 'GDF'))
  p0 = replicate(5,list(NULL))
  for(i in 1:5){
    p0[[i]] = ggplot(TO_PLOT[TO_PLOT$ESA == i+1,])+
      geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = ena_rate, color = ena_rate))+
      theme_classic()+
      xlab('Longitude')+ylab('Latitude')+
      theme(legend.position = 'right', legend.key.size = unit(0.5,'cm'), axis.line.x = element_blank())+
      scale_color_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black' )+
      scale_fill_gradientn('ENAs/sec',colors = ibex_palette$hex ,na.value = 'black' )+
      scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
      scale_y_continuous(expand = c(0,0))+
      theme(panel.background = element_rect(fill = NA, color = "black"))+
      facet_grid(paste0('ESA ',ESA)~Type)+
      coord_fixed(ratio = 1, clip = 'on', xlim=c(360,1), ylim = c(-89,89))
  }
  p1 = ggpubr::ggarrange(p0[[1]],p0[[2]],p0[[3]],p0[[4]],p0[[5]], ncol=1, nrow =5, common.legend = F)
  ggsave(p1,filename = paste0(current_directory, '/Plots/RealData/Theseus_',time,'.pdf'), width = 10, height = 10)
}





##############################
### Height of Ribbon Peaks ### 
##############################
PEAKS_THESEUS = THESEUS_ECLIPTIC %>% dplyr::group_by(Time_Group, ESA) %>% 
  dplyr::mutate(maxval = max(est_ribbon, na.rm=T))
PEAKS_THESEUS = PEAKS_THESEUS[!duplicated(paste0(PEAKS_THESEUS$Time_Group, PEAKS_THESEUS$ESA)),]
PEAKS_ISOC = ISOC_ECLIPTIC %>% dplyr::group_by(Time_Group, ESA) %>% 
  dplyr::mutate(maxval = max(est_ribbon, na.rm=T))
PEAKS_ISOC = PEAKS_ISOC[!duplicated(paste0(PEAKS_ISOC$Time_Group, PEAKS_ISOC$ESA)),]
COLORS = colorspace::diverging_hcl(5,'Tofino')
p2 = ggplot(PEAKS_THESEUS)+
  geom_line(aes(x=Time_Group, y = maxval, color = factor(ESA), group = factor(ESA)), size = 2)+
  theme_classic()+
  xlab('Map')+ylab('Peak Ribbon ENA Rate')+
  theme(legend.position = c(0.65,0.8), legend.direction = "horizontal")+
  scale_fill_manual('ESA', values = COLORS)+
  scale_color_manual('ESA', values = COLORS)+
  labs(title = 'Comparison of estimated ribbon peak heights')+
  theme(panel.background = element_rect(fill = NA, color = "black"))
ggsave(p2, file = paste0(current_directory, '/Plots/RealData/Theseus_Peaks.pdf'), width = 6, height = 3)
p2 = ggplot(PEAKS_ISOC)+
  geom_line(aes(x=Time_Group, y = maxval, color = factor(ESA), group = factor(ESA)), size = 2)+
  theme_classic()+
  xlab('Map')+ylab('Peak Ribbon ENA Rate')+
  theme(legend.position = c(0.65,0.8), legend.direction = "horizontal")+
  scale_fill_manual('ESA', values = COLORS)+
  scale_color_manual('ESA', values = COLORS)+
  labs(title = 'Comparison of estimated ribbon peak heights')+
  theme(panel.background = element_rect(fill = NA, color = "black"))
ggsave(p2, file = paste0(current_directory, '/Plots/RealData/ISOC_Peaks.pdf'), width = 6, height = 3)








#################################################
#################################################
### Plots in Ribbon-Centric Coordinate System ###
#################################################
#################################################

estFWHM = function(x,y){
  spline_func = splines::periodicSpline(x,y,period = 360, ord = 4L)
  x2= seq(min(x),max(x),0.5)
  y2 = predict(spline_func,x2)$y
  y2[y2<0]=0
  x2 = x2[y2 >= 0.5]
  return(max(x2)-min(x2))
}

estSkew = function(x,y){
  spline_func = splines::periodicSpline(x,y,period = 360, ord = 4L)
  x2= seq(min(x),max(x))
  y2 = predict(spline_func,x2)$y
  y2[y2<0]=0
  k = Weighted.Desc.Stat::w.skewness(x2, y2/sum(y2,na.rm=T))
  return(k)
}

#############################
### Theseus Ribbon Widths ### 
#############################

THESEUS_WIDTHS = THESEUS_RC[order(THESEUS_RC$azimuthal_discrete, THESEUS_RC$polar_discrete),]
THESEUS_WIDTHS = THESEUS_WIDTHS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>% 
  dplyr::mutate(peak_loc = polar_discrete[which.max(est_ribbon)])
THESEUS_WIDTHS =  THESEUS_WIDTHS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>%
  dplyr::mutate(est_ribbon0= c(rev(cummin(rev(est_ribbon[polar_discrete <= peak_loc]))),cummin(est_ribbon[polar_discrete > peak_loc])))
THESEUS_WIDTHS$est_ribbon = THESEUS_WIDTHS$est_ribbon0
THESEUS_WIDTHS = THESEUS_WIDTHS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>% 
  dplyr::mutate(max_val = max(est_ribbon,na.rm=T))
THESEUS_WIDTHS$k = THESEUS_WIDTHS$est_ribbon/THESEUS_WIDTHS$max_val
THESEUS_WIDTHS = THESEUS_WIDTHS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>% 
  dplyr::mutate(k_width = estFWHM(polar_discrete, k))
THESEUS_WIDTHS = THESEUS_WIDTHS[!duplicated(paste0(THESEUS_WIDTHS$azimuthal_discrete, THESEUS_WIDTHS$Time_Group, THESEUS_WIDTHS$ESA)),]

THESEUS_SKEWS = THESEUS_RC[order(THESEUS_RC$azimuthal_discrete, THESEUS_RC$polar_discrete),]
THESEUS_SKEWS = THESEUS_SKEWS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>% 
  dplyr::mutate(peak_loc = polar_discrete[which.max(est_ribbon)])
THESEUS_SKEWS =  THESEUS_SKEWS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>%
  dplyr::mutate(est_ribbon0= c(rev(cummin(rev(est_ribbon[polar_discrete <= peak_loc]))),cummin(est_ribbon[polar_discrete > peak_loc])))
THESEUS_SKEWS$est_ribbon = THESEUS_SKEWS$est_ribbon0
THESEUS_SKEWS = THESEUS_SKEWS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>% 
  dplyr::mutate(k_skew = estSkew(polar_discrete, est_ribbon))
THESEUS_SKEWS = THESEUS_SKEWS[!duplicated(paste0(THESEUS_SKEWS$azimuthal_discrete, THESEUS_SKEWS$Time_Group, THESEUS_SKEWS$ESA)),]

p1 = ggplot(THESEUS_SKEWS)+
  geom_hline(yintercept = 0, linetype = 1, size = 2, color = 'gray')+
  geom_point(aes(x=factor(ESA), y = k_skew, color = azimuthal_discrete, group = azimuthal_discrete), shape = 16,
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.5), alpha = 1)+
  geom_boxplot(aes(x=factor(ESA), y = k_skew),size = 1.5 , fill = NA, outlier.shape = NA, color = 'gray20')+
  geom_boxplot(aes(x=factor(ESA), y = k_skew),size = 0.7 , fill = NA, outlier.shape = NA, color = 'white')+
  xlab('Energy Step (ESA)')+ylab('Skewness')+
  theme_classic()+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'))+
  labs(title = 'Skewness')+
  scale_fill_gradientn('Azimuthal angle',colors = met.brewer("Johnson"))+
  scale_color_gradientn('Azimuthal angle',colors = met.brewer("Johnson"))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(legend.position = 'top')+
  coord_cartesian(ylim = c(-0.9,0.9))
p2 = ggplot(THESEUS_WIDTHS)+
  geom_point(aes(x=factor(ESA), y = k_width, color = azimuthal_discrete, group = azimuthal_discrete), shape = 16,
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.5), alpha = 1)+
  geom_boxplot(aes(x=factor(ESA), y = k_width),size = 1.5 , fill = NA, outlier.shape = NA, color = 'gray20')+
  geom_boxplot(aes(x=factor(ESA), y = k_width),size = 0.7 , fill = NA, outlier.shape = NA, color = 'white')+
  xlab('Energy Step (ESA)')+ylab('Width (Full Width at Half Max)')+
  theme_classic()+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'))+
  labs(title = 'Width')+
  scale_fill_gradientn('Azimuthal angle',colors = met.brewer("Johnson"))+
  scale_color_gradientn('Azimuthal angle',colors = met.brewer("Johnson"))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(legend.position = 'top')
p3 = ggpubr::ggarrange(p1,p2, ncol = 2, nrow = 1, common.legend = T)
ggsave(p3,filename = paste0(current_directory, '/Plots/RealData/Theseus_WidthSkew.pdf'), width = 10, height = 4)

aggregate(k_width~ESA,FUN = median, data = THESEUS_WIDTHS)
# ESA k_width
# 1   2    25.0
# 2   3    19.0
# 3   4    15.5
# 4   5    23.0
# 5   6    31.5



##########################
### ISOC Ribbon Widths ### 
##########################

ISOC_WIDTHS = ISOC_RC[order(ISOC_RC$azimuthal_discrete, ISOC_RC$polar_discrete),]
ISOC_WIDTHS = ISOC_WIDTHS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>% 
  dplyr::mutate(peak_loc = polar_discrete[which.max(est_ribbon)])
ISOC_WIDTHS =  ISOC_WIDTHS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>%
  dplyr::mutate(est_ribbon0= c(rev(cummin(rev(est_ribbon[polar_discrete <= peak_loc]))),cummin(est_ribbon[polar_discrete > peak_loc])))
ISOC_WIDTHS$est_ribbon = ISOC_WIDTHS$est_ribbon0
ISOC_WIDTHS = ISOC_WIDTHS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>% 
  dplyr::mutate(max_val = max(est_ribbon,na.rm=T))
ISOC_WIDTHS$k = ISOC_WIDTHS$est_ribbon/ISOC_WIDTHS$max_val
ISOC_WIDTHS = ISOC_WIDTHS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>% 
  dplyr::mutate(k_width = estFWHM(polar_discrete, k))
ISOC_WIDTHS = ISOC_WIDTHS[!duplicated(paste0(ISOC_WIDTHS$azimuthal_discrete, ISOC_WIDTHS$Time_Group, ISOC_WIDTHS$ESA)),]

ISOC_SKEWS = ISOC_RC[order(ISOC_RC$azimuthal_discrete, ISOC_RC$polar_discrete),]
ISOC_SKEWS = ISOC_SKEWS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>% 
  dplyr::mutate(peak_loc = polar_discrete[which.max(est_ribbon)])
ISOC_SKEWS =  ISOC_SKEWS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>%
  dplyr::mutate(est_ribbon0= c(rev(cummin(rev(est_ribbon[polar_discrete <= peak_loc]))),cummin(est_ribbon[polar_discrete > peak_loc])))
ISOC_SKEWS$est_ribbon = ISOC_SKEWS$est_ribbon0
ISOC_SKEWS = ISOC_SKEWS %>% dplyr::group_by(Time_Group, ESA, azimuthal_discrete) %>% 
  dplyr::mutate(k_skew = estSkew(polar_discrete, est_ribbon))
ISOC_SKEWS = ISOC_SKEWS[!duplicated(paste0(ISOC_SKEWS$azimuthal_discrete, ISOC_SKEWS$Time_Group, ISOC_SKEWS$ESA)),]

p1 = ggplot(ISOC_SKEWS)+
  geom_hline(yintercept = 0, linetype = 1, size = 2, color = 'gray')+
  geom_point(aes(x=factor(ESA), y = k_skew, color = azimuthal_discrete, group = azimuthal_discrete), shape = 16,
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.5), alpha = 1)+
  geom_boxplot(aes(x=factor(ESA), y = k_skew),size = 1.5 , fill = NA, outlier.shape = NA, color = 'gray20')+
  geom_boxplot(aes(x=factor(ESA), y = k_skew),size = 0.7 , fill = NA, outlier.shape = NA, color = 'white')+
  xlab('Energy Step (ESA)')+ylab('Skewness')+
  theme_classic()+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'))+
  labs(title = 'Skewness')+
  scale_fill_gradientn('Azimuthal angle',colors = met.brewer("Johnson"))+
  scale_color_gradientn('Azimuthal angle',colors = met.brewer("Johnson"))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(legend.position = 'top')+
  coord_cartesian(ylim = c(-0.9,0.9))
p2 = ggplot(ISOC_WIDTHS)+
  geom_point(aes(x=factor(ESA), y = k_width, color = azimuthal_discrete, group = azimuthal_discrete), shape = 16,
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.5), alpha = 1)+
  geom_boxplot(aes(x=factor(ESA), y = k_width),size = 1.5 , fill = NA, outlier.shape = NA, color = 'gray20')+
  geom_boxplot(aes(x=factor(ESA), y = k_width),size = 0.7 , fill = NA, outlier.shape = NA, color = 'white')+
  xlab('Energy Step (ESA)')+ylab('Width (Full Width at Half Max)')+
  theme_classic()+
  theme(legend.position = 'top', legend.key.size = unit(1,'cm'))+
  labs(title = 'Width')+
  scale_fill_gradientn('Azimuthal angle',colors = met.brewer("Johnson"))+
  scale_color_gradientn('Azimuthal angle',colors = met.brewer("Johnson"))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  theme(legend.position = 'top')
p3 = ggpubr::ggarrange(p1,p2, ncol = 2, nrow = 1, common.legend = T)
ggsave(p3,filename = paste0(current_directory, '/Plots/RealData/ISOC_WidthSkew.pdf'), width = 10, height = 4)




aggregate(k_width~ESA,FUN = median, data = ISOC_WIDTHS)
# ESA k_width
# 1   2    20.0
# 2   3    20.5
# 3   4    20.5
# 4   5    24.0
# 5   6    27.0



############################
### Theseus 2019 Example ### 
############################

TEMP = THESEUS_ECLIPTIC[THESEUS_ECLIPTIC$Time_Group == '2019A' & THESEUS_ECLIPTIC$ESA == 4,]
horizontallinesdf <- expand.grid(ecliptic_lat = seq(-60,60,30),
                                 ecliptic_lon = seq(0,360,30))
LEVELS = round(seq(0,0.25,length.out = 5),2)
p0 = ggplot(TEMP)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = as.numeric(input_data), color = as.numeric(input_data)))+
  theme_classic()+
  xlab('')+ylab('')+
  #labs(title = 'Input Map')+
  coord_map("mollweide", orientation = c(90,0,180), clip = "off")+
  #scale_fill_gradientn('ENA Rate',breaks = LEVELS,colors = c('black', 'blue', 'gold', 'orange', 'red') ,limits = c(0,max(LEVELS)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), na.value = 'black' )+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), na.value = 'black' )+
  theme(plot.margin = margin(t=5,5,5,5, "lines")) +
  theme(legend.direction = "horizontal") +
  theme(legend.position = c(0.5, 1.5))+
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())+
  guides(color = 'none', fill = 'none')+
  scale_x_reverse(expand = c(0,0))+
  geom_line(aes(x=ecliptic_lon, y=ecliptic_lat, group=ecliptic_lat), color = I('white'), data = horizontallinesdf, alpha = 0.7, size = 0.2)+
  geom_vline(xintercept = seq(0.0001,360-0.0001,length.out = 13), color = 'white', alpha = 0.7, size = 0.2)
p1 = ggplot(TEMP)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = as.numeric(est_ribbon), color = as.numeric(est_ribbon)))+
  theme_classic()+
  theme(plot.margin = margin(t=5,5,5,5, "lines")) +
  theme(legend.direction = "horizontal") +
  theme(legend.position = c(0.5, 1.5))+
  xlab('')+ylab('')+
  #labs(title = 'Ribbon Map')+
  coord_map("mollweide", orientation = c(90,0,180), clip = "off")+
  #scale_fill_gradientn('ENA Rate',breaks = LEVELS,colors = c('black', 'blue', 'gold', 'orange', 'red') ,limits = c(0,max(LEVELS)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), na.value = 'black' )+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), na.value = 'black' )+
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())+
  guides(color = 'none', fill = 'none')+
  scale_x_reverse(expand = c(0,0))+
  geom_line(aes(x=ecliptic_lon, y=ecliptic_lat, group=ecliptic_lat), color = I('white'), data = horizontallinesdf, alpha = 0.7, size = 0.2)+
  geom_vline(xintercept = seq(0.0001,360-0.0001,length.out = 13), color = 'white', alpha = 0.7, size = 0.2)
p2 = ggplot(TEMP)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = as.numeric(est_gdf), color = as.numeric(est_gdf)))+
  theme_classic()+
  theme(plot.margin = margin(t=5,5,5,5, "lines")) +
  theme(legend.direction = "horizontal") +
  theme(legend.position = c(0.5, 1.5))+
  xlab('')+ylab('')+
  #labs(title = 'GDF Map')+
  coord_map("mollweide", orientation = c(90,0,180),clip = "off" )+
  #scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(RESULTS_ECLIPTIC$data[RESULTS_ECLIPTIC$ESA == ESA])), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), na.value = 'black' )+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), na.value = 'black' )+
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())+
  guides(color = 'none', fill = 'none')+
  scale_x_reverse(expand = c(0,0))+
  geom_line(aes(x=ecliptic_lon, y=ecliptic_lat, group=ecliptic_lat), color = I('white'), data = horizontallinesdf, alpha = 0.7, size = 0.2)+
  geom_vline(xintercept = seq(0.0001,360-0.0001,length.out = 13), color = 'white', alpha = 0.7, size = 0.2)
p3 = get_legend(ggplot(TEMP)+
                  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = as.numeric(est_gdf)))+
                  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), na.value = 'black' )+
                  theme(legend.position = 'top', legend.key.width = unit(1.5, "cm")))
ggsave(p0,file = paste0(current_directory, '/Plots/RealData/Theseus_ExTotal.png'), width = 4, height = 3)
ggsave(p1,file = paste0(current_directory, '/Plots/RealData/Theseus_ExRibbon.png'), width = 4, height = 3)
ggsave(p2,file = paste0(current_directory, '/Plots/RealData/Theseus_ExResidual.png'), width = 4, height = 3)
ggsave(p3,file = paste0(current_directory, '/Plots/RealData/Theseus_ExLegend.png'), width = 4, height = 2)



TEMP = THESEUS_RC[THESEUS_RC$Time_Group == '2019A' & THESEUS_RC$ESA == 4,]
horizontallinesdf <- expand.grid(ecliptic_lat = seq(-60,60,30),
                                 ecliptic_lon = seq(0,360,30))
LEVELS = round(seq(0,0.25,length.out = 5),2)
p0 = ggplot(TEMP)+
  geom_tile(aes(x=lon, y=lat, fill = as.numeric(input_data), color = as.numeric(input_data)))+
  theme_classic()+
  xlab('')+ylab('')+
  #labs(title = 'Input Map')+
  coord_map("mollweide", orientation = c(90,0,180), clip = "off")+
  #scale_fill_gradientn('ENA Rate',breaks = LEVELS,colors = c('black', 'blue', 'gold', 'orange', 'red') ,limits = c(0,max(LEVELS)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), na.value = 'black' )+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), na.value = 'black' )+
  theme(plot.margin = margin(t=5,5,5,5, "lines")) +
  theme(legend.direction = "horizontal") +
  theme(legend.position = c(0.5, 1.5))+
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())+
  guides(color = 'none', fill = 'none')+
  scale_x_continuous(expand = c(0,0))+
  geom_line(aes(x=ecliptic_lon, y=ecliptic_lat, group=ecliptic_lat), color = I('white'), data = horizontallinesdf, alpha = 0.7, size = 0.2)+
  geom_vline(xintercept = seq(0.0001,360-0.0001,length.out = 13), color = 'white', alpha = 0.7, size = 0.2)
ggsave(p0,file = paste0(current_directory, '/Plots/RealData/Theseus_ExRotated.png'), width = 4, height = 3)


# TEMP$test = TEMP$est_ribbon_reisenfeld+TEMP$est_gdf_reisenfeld
# TEMP = THESEUS_RC[THESEUS_RC$Time_Group == '2019A' & THESEUS_RC$ESA == 4,]
# TEMP$test = TEMP$est_ribbon_reisenfeld+TEMP$est_gdf_reisenfeld
# TEMP = ISOC_RC[ISOC_RC$Time_Group == '2019' & ISOC_RC$ESA == 4,]
# TEMP$test = TEMP$est_ribbon_reisenfeld+TEMP$est_gdf_reisenfeld
# TEMP = ISOC_ECLIPTIC[ISOC_ECLIPTIC$Time_Group == '2019' & ISOC_ECLIPTIC$ESA == 4,]
# TEMP$test = TEMP$est_ribbon_reisenfeld+TEMP$est_gdf_reisenfeld

TEMP = THESEUS_ECLIPTIC[THESEUS_ECLIPTIC$Time_Group == '2019A' & THESEUS_ECLIPTIC$ESA == 4,]
horizontallinesdf <- expand.grid(ecliptic_lat = seq(-60,60,30),
                                 ecliptic_lon = seq(0,360,30))
LEVELS = round(seq(0,0.25,length.out = 5),2)
p1 = ggplot(TEMP)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = as.numeric(est_ribbon_reisenfeld), color = as.numeric(est_ribbon_reisenfeld)))+
  theme_classic()+
  theme(plot.margin = margin(t=5,5,5,5, "lines")) +
  theme(legend.direction = "horizontal") +
  theme(legend.position = c(0.5, 1.5))+
  xlab('')+ylab('')+
  #labs(title = 'Ribbon Map')+
  coord_map("mollweide", orientation = c(90,0,180), clip = "off")+
  #scale_fill_gradientn('ENA Rate',breaks = LEVELS,colors = c('black', 'blue', 'gold', 'orange', 'red') ,limits = c(0,max(LEVELS)), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), na.value = 'black' )+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), na.value = 'black' )+
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())+
  guides(color = 'none', fill = 'none')+
  scale_x_reverse(expand = c(0,0))+
  geom_line(aes(x=ecliptic_lon, y=ecliptic_lat, group=ecliptic_lat), color = I('white'), data = horizontallinesdf, alpha = 0.7, size = 0.2)+
  geom_vline(xintercept = seq(0.0001,360-0.0001,length.out = 13), color = 'white', alpha = 0.7, size = 0.2)
p2 = ggplot(TEMP)+
  geom_tile(aes(x=ecliptic_lon_center, y=lat, fill = as.numeric(est_gdf_reisenfeld), color = as.numeric(est_gdf_reisenfeld)))+
  theme_classic()+
  theme(plot.margin = margin(t=5,5,5,5, "lines")) +
  theme(legend.direction = "horizontal") +
  theme(legend.position = c(0.5, 1.5))+
  xlab('')+ylab('')+
  #labs(title = 'GDF Map')+
  coord_map("mollweide", orientation = c(90,0,180),clip = "off" )+
  #scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(RESULTS_ECLIPTIC$data[RESULTS_ECLIPTIC$ESA == ESA])), na.value = 'black' )+
  scale_fill_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), oob = scales::oob_squish )+
  scale_color_gradientn('ENA Rate',colors = ibex_palette$hex ,limits = c(0,max(TEMP$input_data)), oob = scales::oob_squish  )+
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())+
  guides(color = 'none', fill = 'none')+
  scale_x_reverse(expand = c(0,0))+
  geom_line(aes(x=ecliptic_lon, y=ecliptic_lat, group=ecliptic_lat), color = I('white'), data = horizontallinesdf, alpha = 0.7, size = 0.2)+
  geom_vline(xintercept = seq(0.0001,360-0.0001,length.out = 13), color = 'white', alpha = 0.7, size = 0.2)

ggsave(p1,file = paste0(current_directory, '/Plots/RealData/Theseus_ExRibbon_Reisenfeld.png'), width = 4, height = 3)
ggsave(p2,file = paste0(current_directory, '/Plots/RealData/Theseus_ExResidual_Reisenfeld.png'), width = 4, height = 3)




################
### RC Frame ###
################
DATA_TO_SEP = read.csv(paste0(current_directory, '/output_theseus_ecliptic.csv'))
MAPS_TO_RUN = expand.grid(Time_Group = unique(DATA_TO_SEP$Time_Group),
                          ESA = unique(DATA_TO_SEP$ESA))
MAPS_TO_RUN = MAPS_TO_RUN[MAPS_TO_RUN$ESA == 4,] #only run for ESA 4 maps
RESULTS_RC = NULL
for(j in 1:length(MAPS_TO_RUN[,1])){
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


DAT_SUB = RESULTS_RC[order(RESULTS_RC$lon, RESULTS_RC$lat),]
p2 = ggplot(DAT_SUB[DAT_SUB$ESA == 4,])+
  geom_line(aes(x=lat+90, y=est_ribbon, group = factor(lon), color = lon ), size = 0.3, alpha = 0.7)+
  theme_classic()+
  xlab('Polar Angle')+ylab('ENAs/sec')+
  facet_wrap(Time_Group~.)+
  theme(legend.position = 'top')+
  coord_cartesian(xlim = c(75,140))+
  scale_color_gradientn('Azimuthal Angle',colors = met.brewer("Johnson"))
ggsave(p2,file = paste0(current_directory, '/Plots/RealData/Theseus_UpdatedRCMapsRam4.png'), width = 10, height = 8)












