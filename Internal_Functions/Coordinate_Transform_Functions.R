

get_mapval = function(x,MAP_TEMP){
  temp = MAP_TEMP[as.numeric(unlist(MAP_TEMP$lon_lower)) <= as.numeric(x[1]) & as.numeric(unlist(MAP_TEMP$lon_upper))>as.numeric(x[1]),]
  temp = temp[as.numeric(unlist(temp$lat_lower)) <= as.numeric(x[2]) & as.numeric(unlist(temp$lat_upper))>as.numeric(x[2]),]
  return(rbind(temp[1,c('ena_rate')]))
}  
get_mapval_sd = function(x,MAP_TEMP){
  temp = MAP_TEMP[as.numeric(unlist(MAP_TEMP$lon_lower)) <= as.numeric(x[1]) & as.numeric(unlist(MAP_TEMP$lon_upper))>as.numeric(x[1]),]
  temp = temp[as.numeric(unlist(temp$lat_lower)) <= as.numeric(x[2]) & as.numeric(unlist(temp$lat_upper))>as.numeric(x[2]),]
  return(rbind(temp[1,c('ena_rate', 'sd_ena_rate')]))
} 

### RIBBON_CENTER = c(lon, lat) , e.g. 221, 39.9
### Assumes latitude in [-90,90]
### Assumes longitude in [0,360]
### MAP must have columns lon, lat, ena_rate
### Column sd_ena_rate optional
### Reverse indicates whether we want to go from ribbon-centric to ecliptic but only know the working ribbon center for the current map
Coord_Transform_Map = function(MAP, RIBBON_CENTER, METHOD = 'Inversion', STAND_DEV = TRUE, nbreaks, REVERSE = FALSE){
  K = unique(diff(sort(unique(MAP$lon))))
  if(is.null(nbreaks)){stop('Specify number of breaks per pixel')}
  if(REVERSE){
    TRANSFORM_MAT_REVERSE = davenport_rotation(yaw = RIBBON_CENTER[1]*pi/180, pitch = -(90-RIBBON_CENTER[2])*pi/180)
    TRANSFORM_MAT = solve(TRANSFORM_MAT_REVERSE)
  }else{
    TRANSFORM_MAT_REVERSE = solve(davenport_rotation(yaw = RIBBON_CENTER[1]*pi/180, pitch = -(90-RIBBON_CENTER[2])*pi/180))
    TRANSFORM_MAT = solve(TRANSFORM_MAT_REVERSE)
  }
  
  if(METHOD == 'Inversion'){
    ### Define pixel breaks
    MAP_PARTITION = expand.grid(lon_new= seq(0, 360, K/nbreaks), lat_new = seq(-90, 90, K/nbreaks))
    MAP_PARTITION$input1_new = as.numeric(MAP_PARTITION$lon_new - 180)*pi/180
    MAP_PARTITION$input2_new = as.numeric(MAP_PARTITION$lat_new)*pi/180
    ### Invert to current coordinate system
    MAP_CART = data.frame(pracma::sph2cart(as.matrix(data.frame(MAP_PARTITION[,c('input1_new', 'input2_new')],1)) ))  
    names(MAP_CART)=c('x_new', 'y_new', 'z_new')
    temp = as.matrix(MAP_CART[,c('x_new','y_new','z_new')]) %*% TRANSFORM_MAT_REVERSE
    MAP_PARTITION$x_old = temp[,1]
    MAP_PARTITION$y_old = temp[,2]
    MAP_PARTITION$z_old = temp[,3]
    ### Transform to lat/lon in current coordinate system
    MAP_TRANS_SPHERICAL = data.frame(cart2sph(as.matrix(MAP_PARTITION[,c('x_old', 'y_old', 'z_old')])))
    MAP_PARTITION$lon_old = MAP_TRANS_SPHERICAL$theta*180/pi + 180
    MAP_PARTITION$lat_old= (MAP_TRANS_SPHERICAL$phi)*180/pi
    MAP_PARTITION$lon_binned_old = cut(MAP_PARTITION$lon_old, breaks = seq(0, 360, K), labels = sort(unique(MAP$lon)), include.lowest = T)
    MAP_PARTITION$lat_binned_old = cut(MAP_PARTITION$lat_old, breaks = seq(-90, 90, K), labels = sort(unique(MAP$lat)), include.lowest = T)
    ### Calculate the ena rate at each corner
    if(!STAND_DEV){
      MAP_PARTITION = merge(MAP_PARTITION, data.frame(lon_binned_old = MAP$lon, lat_binned_old = MAP$lat, ena_rate = MAP$ena_rate), all.x = T, all.y = F)
    }else{
      MAP_PARTITION = merge(MAP_PARTITION, data.frame(lon_binned_old = MAP$lon, lat_binned_old = MAP$lat, ena_rate = MAP$ena_rate, sd_ena_rate = MAP$sd_ena_rate), all.x = T, all.y = F)
    }
  }else if(METHOD == 'Partition'){
    ### Define pixel breaks
    MAP_PARTITION = expand.grid(lon_old= seq(0, 360, K/nbreaks), lat_old = seq(-90, 90, K/nbreaks))
    MAP_PARTITION$lon_binned_old = cut(MAP_PARTITION$lon_old, breaks = seq(0, 360, K), labels = sort(unique(MAP$lon)), include.lowest = T)
    MAP_PARTITION$lat_binned_old = cut(MAP_PARTITION$lat_old, breaks = seq(-90, 90, K), labels = sort(unique(MAP$lat)), include.lowest = T)
    if(!STAND_DEV){
      MAP_PARTITION = merge(MAP_PARTITION, data.frame(lon_binned_old = MAP$lon, lat_binned_old = MAP$lat, ena_rate = MAP$ena_rate), all.x = T, all.y = F)
    }else{
      MAP_PARTITION = merge(MAP_PARTITION, data.frame(lon_binned_old = MAP$lon, lat_binned_old = MAP$lat, ena_rate = MAP$ena_rate, sd_ena_rate = MAP$sd_ena_rate), all.x = T, all.y = F)
    }
    ### Convert to new coordinate system
    MAP_PARTITION$input1_old = as.numeric(MAP_PARTITION$lon_old - 180)*pi/180
    MAP_PARTITION$input2_old = as.numeric(MAP_PARTITION$lat_old)*pi/180
    MAP_CART = data.frame(pracma::sph2cart(as.matrix(data.frame(MAP_PARTITION[,c('input1_old', 'input2_old')],1)) ))  
    names(MAP_CART)=c('x_old', 'y_old', 'z_old')
    ### Transform to new lat/lon in current coordinate system
    temp = as.matrix(MAP_CART[,c('x_old','y_old','z_old')]) %*% TRANSFORM_MAT
    MAP_PARTITION$x_new = temp[,1]
    MAP_PARTITION$y_new = temp[,2]
    MAP_PARTITION$z_new = temp[,3]
    MAP_PARTITION_SPHERICAL = cart2sph(as.matrix(MAP_PARTITION[,c('x_new', 'y_new', 'z_new')]))
    MAP_PARTITION$lon_new = MAP_PARTITION_SPHERICAL[,1]*180/pi+180
    MAP_PARTITION$lat_new= MAP_PARTITION_SPHERICAL[,2]*180/pi 
  }

  ### Define new pixel groupings 
  MAP_PARTITION$lon_binned_new = cut(MAP_PARTITION$lon_new, breaks = seq(0, 360, K), labels = sort(unique(MAP$lon)), include.lowest = T)
  MAP_PARTITION$lat_binned_new = cut(MAP_PARTITION$lat_new, breaks = seq(-90, 90, K), labels = sort(unique(MAP$lat)), include.lowest = T)
  ### Aggregate into pixels of equal size to input maps
  if(!STAND_DEV){
    MAP_PARTITION_AGG = aggregate(ena_rate~lon_binned_new+lat_binned_new, FUN = function(x){mean(x,na.rm=T)}, drop = F, na.action = na.pass, data = MAP_PARTITION)
  }else{
    MAP_PARTITION$var_ena_rate = (MAP_PARTITION$sd_ena_rate^2)
    MAP_PARTITION_AGG = aggregate(cbind(ena_rate,var_ena_rate)~lon_binned_new+lat_binned_new, FUN = function(x){mean(x,na.rm=T)}, drop = F, na.action = na.pass, data = MAP_PARTITION)
    MAP_PARTITION_AGG$sd_ena_rate = sqrt(MAP_PARTITION_AGG$var_ena_rate)
  }
  names(MAP_PARTITION_AGG)[names(MAP_PARTITION_AGG)=='lon_binned_new']='lon'
  names(MAP_PARTITION_AGG)[names(MAP_PARTITION_AGG)=='lat_binned_new']='lat'
  MAP_PARTITION_AGG$lon = as.numeric(as.character(MAP_PARTITION_AGG$lon))
  MAP_PARTITION_AGG$lat = as.numeric(as.character(MAP_PARTITION_AGG$lat))
  MAP_PARTITION_AGG = MAP_PARTITION_AGG[order(MAP_PARTITION_AGG$lon, MAP_PARTITION_AGG$lat),]
  if(!STAND_DEV){
    return(MAP_PARTITION_AGG[,c('lon', 'lat', 'ena_rate')])
  }else{
    return(MAP_PARTITION_AGG[,c('lon', 'lat', 'ena_rate','sd_ena_rate')])
    
  }
}


### RIBBON_CENTER = c(lon, lat) , e.g. 221, 39.9
### Assumes latitude in [-90,90]
### Assumes longitude in [0,360]
### MAP must have columns lon, lat
### Reverse indicates whether we want to go from ribbon-centric to ecliptic but only know the working ribbon center for the current map
Coord_Transform_Points = function(MAP, RIBBON_CENTER, REVERSE = FALSE){
  K = unique(diff(sort(unique(MAP$lon))))
  MAP$lon_old = MAP$lon
  MAP$lat_old = MAP$lat
  MAP_CART = data.frame(t(apply(MAP[,c('lon_old', 'lat_old')], 1, FUN = function(x){pracma::sph2cart(as.vector(c(as.numeric(x[1]-180)*pi/180,as.numeric(x[2])*pi/180, 1)))})))
  names(MAP_CART)=c('x_old', 'y_old', 'z_old')
  if(REVERSE){
    TRANSFORM_MAT = solve(davenport_rotation(yaw = RIBBON_CENTER[1]*pi/180, pitch = -(90-RIBBON_CENTER[2])*pi/180))
  }else{
    TRANSFORM_MAT = davenport_rotation(yaw = RIBBON_CENTER[1]*pi/180, pitch = -(90-RIBBON_CENTER[2])*pi/180)
  }
  MAP_TRANS = as.matrix(MAP_CART[,c('x_old','y_old','z_old')]) %*% TRANSFORM_MAT
  MAP$x_new = MAP_TRANS[,1]
  MAP$y_new = MAP_TRANS[,2]
  MAP$z_new = MAP_TRANS[,3]
  
  MAP_TRANS_SPHERICAL = data.frame(cart2sph(as.matrix(MAP[,c('x_new', 'y_new', 'z_new')])))
  MAP$lon_new = MAP_TRANS_SPHERICAL$theta*180/pi + 180
  MAP$lat_new= MAP_TRANS_SPHERICAL$phi*180/pi
  return(MAP[,c('lon_old', 'lat_old', 'lon_new', 'lat_new')])
}






Coord_Transform_Matrix = function(MAP, RIBBON_CENTER, nbreaks, REVERSE = FALSE){
  K = unique(diff(sort(unique(MAP$lon))))
  if(is.null(nbreaks)){stop('Specify number of breaks per pixel')}
  if(REVERSE){
    TRANSFORM_MAT_REVERSE = davenport_rotation(yaw = RIBBON_CENTER[1]*pi/180, pitch = -(90-RIBBON_CENTER[2])*pi/180)
    TRANSFORM_MAT = solve(TRANSFORM_MAT_REVERSE)
  }else{
    TRANSFORM_MAT_REVERSE = solve(davenport_rotation(yaw = RIBBON_CENTER[1]*pi/180, pitch = -(90-RIBBON_CENTER[2])*pi/180))
    TRANSFORM_MAT = solve(TRANSFORM_MAT_REVERSE)
  }
  
  ### Define pixel breaks
  MAP_PARTITION = expand.grid(lon_new= seq(0, 360, K/nbreaks), lat_new = seq(-90, 90, K/nbreaks))
  MAP_PARTITION$input1_new = as.numeric(MAP_PARTITION$lon_new - 180)*pi/180
  MAP_PARTITION$input2_new = as.numeric(MAP_PARTITION$lat_new)*pi/180
  ### Invert to current coordinate system
  MAP_CART = data.frame(pracma::sph2cart(as.matrix(data.frame(MAP_PARTITION[,c('input1_new', 'input2_new')],1)) ))  
  names(MAP_CART)=c('x_new', 'y_new', 'z_new')
  temp = as.matrix(MAP_CART[,c('x_new','y_new','z_new')]) %*% TRANSFORM_MAT_REVERSE
  MAP_PARTITION$x_old = temp[,1]
  MAP_PARTITION$y_old = temp[,2]
  MAP_PARTITION$z_old = temp[,3]
  ### Transform to lat/lon in current coordinate system
  MAP_TRANS_SPHERICAL = data.frame(cart2sph(as.matrix(MAP_PARTITION[,c('x_old', 'y_old', 'z_old')])))
  MAP_PARTITION$lon_old = MAP_TRANS_SPHERICAL$theta*180/pi + 180
  MAP_PARTITION$lat_old= (MAP_TRANS_SPHERICAL$phi)*180/pi
  MAP_PARTITION$lon_binned_old = cut(MAP_PARTITION$lon_old, breaks = seq(0, 360, K), labels = sort(unique(MAP$lon)), include.lowest = T)
  MAP_PARTITION$lat_binned_old = cut(MAP_PARTITION$lat_old, breaks = seq(-90, 90, K), labels = sort(unique(MAP$lat)), include.lowest = T)
  
  ### Define new pixel groupings 
  MAP_PARTITION$lon_binned_new = cut(MAP_PARTITION$lon_new, breaks = seq(0, 360, K), labels = sort(unique(MAP$lon)), include.lowest = T)
  MAP_PARTITION$lat_binned_new = cut(MAP_PARTITION$lat_new, breaks = seq(-90, 90, K), labels = sort(unique(MAP$lat)), include.lowest = T)
  
  
  MAP_PARTITION = MAP_PARTITION[,c('lon_binned_old', 'lat_binned_old','lon_binned_new', 'lat_binned_new')]
  MAP_PARTITION$key_old = paste0(MAP_PARTITION$lon_binned_old, '_', MAP_PARTITION$lat_binned_old)
  MAP_PARTITION$key_new = paste0(MAP_PARTITION$lon_binned_new, '_', MAP_PARTITION$lat_binned_new)
  MAP_PARTITION$key_all = paste0(MAP_PARTITION$key_old, '_', MAP_PARTITION$key_new)
  
  MAP_PARTITION = MAP_PARTITION %>% dplyr::group_by(key_all) %>% 
    dplyr::mutate(count = length(key_all))
  MAP_PARTITION = MAP_PARTITION[!duplicated(MAP_PARTITION$key_all),]
  MAP$key_old = paste0(MAP$lon, '_', MAP$lat)
  MAP_PARTITION$key_old_num = as.numeric(factor(MAP_PARTITION$key_old, levels = MAP$key_old))
  MAP_PARTITION$key_new_num = as.numeric(factor(MAP_PARTITION$key_new, levels = MAP$key_old))
  

  N = length(unlist(MAP[,1]))
  TEMP = MAP_PARTITION[,c('key_new_num', 'key_old_num', 'count')]
  names(TEMP) = c('x','y','z')
  TEST = with(TEMP, Matrix::sparseMatrix(i = as.numeric(x), j=as.numeric(y), x=z, dims = c(N,N)))
  SUMS = Matrix::rowSums(TEST)
  sweep_sparse <- function(x, margin, stats, fun = "*") {
    f <- match.fun(fun)
    if (margin == 1) {
      idx <- x@i + 1
    } else {
      idx <- x@j + 1
    }
    x@x <- f(x@x, stats[idx])
    return(x)
  }
  OUT = sweep_sparse(TEST,1,1/as.vector(unlist(SUMS)), '*')
  rm(list=c('TEMP','TEST','SUMS','MAP_PARTITION'))
  return(OUT)
}

