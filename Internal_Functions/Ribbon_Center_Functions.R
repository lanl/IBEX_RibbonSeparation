##########################################
##########################################
### Ribbon Center Estimation Functions ###
##########################################
##########################################

### Developed by Dr. Lauren J. Beesley, PhD
### Contact: lvandervort@lanl.gov
### Last Updated: 08/24/2023

### This script provides functions for implementing ribbon center estimation

###################################
### Estimate Best Fitting Plane ###
###################################

### Code adapted from https://meshlogic.github.io/posts/jupyter/curve-fitting/fitting-a-circle-to-cluster-of-3d-points/
GET_PLANE = function(MAT){
  MAT$x_rotated = MAT$x_rotated - MAT$x_rotated_mean
  MAT$y_rotated = MAT$y_rotated - MAT$y_rotated_mean
  MAT$z_rotated = MAT$z_rotated - MAT$z_rotated_mean
  SVD = svd(MAT[,c('x_rotated','y_rotated','z_rotated')])
  normal_vec = SVD$v[,3]
  if(normal_vec[3]<0){
    normal_vec = -normal_vec
  }
  return(normal_vec)
}

##############################
### Estimate ribbon center ###
##############################

### Code adapted from https://meshlogic.github.io/posts/jupyter/curve-fitting/fitting-a-circle-to-cluster-of-3d-points/
GET_CENTER = function(MAT){
  normal_vector = MAT[,c('x_plane', 'y_plane', 'z_plane')][1,]
  ### Perform rodrigues rotation
  P = MAT[,c('x_rotated', 'y_rotated', 'z_rotated')]-MAT[,c('x_rotated_mean', 'y_rotated_mean', 'z_rotated_mean')]
  #Following two runs should give same P_rot
  R = rodrigues_matrix(P, n0=normal_vector, n1 = c(0,0,1))
  P_rot = as.matrix(P) %*% t(R)
  ### Estimate ellipse center
  test = fitConic::fitConic(X = P_rot[,1], Y = P_rot[,2], conicType = 'e')$parA
  # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2288654/
  A <- matrix(c(2*test[1], test[2], test[2], 2*test[3]), nrow=2, ncol=2, byrow=T )
  b <- matrix(c(-test[4], -test[5]), nrow=2, ncol=1, byrow=T)
  soln <- solve(A) %*% b
  center_estimated = soln[1:2]
  zvals_test = seq(0,1,0.001)
  center_3D = as.matrix(data.frame(a=as.vector(unlist(center_estimated)[1]),
                                   b=as.vector(unlist(center_estimated)[2]),
                                   zvals_test)) %*% solve(t(R))
  center_3D[,1] = center_3D[,1]+MAT$x_rotated_mean[1]
  center_3D[,2] = center_3D[,2]+MAT$y_rotated_mean[1]
  center_3D[,3] = center_3D[,3]+MAT$z_rotated_mean[1]
  spherical = cart2sph(center_3D)
  return(center_3D[which.min(abs(spherical[,3]-1)),])
}



GET_CENTER_CIRCULAR = function(MAT){
  normal_vector = MAT[,c('x_plane', 'y_plane', 'z_plane')][1,]
  ### Perform rodrigues rotation
  P = MAT[,c('x_rotated', 'y_rotated', 'z_rotated')]-MAT[,c('x_rotated_mean', 'y_rotated_mean', 'z_rotated_mean')]
  #Following two runs should give same P_rot
  R = rodrigues_matrix(P, n0=normal_vector, n1 = c(0,0,1))
  P_rot = as.matrix(P) %*% t(R)
  
  fit = pracma::circlefit(xp = P_rot[,1], yp = P_rot[,2])

  center_estimated = fit[1:2]
  
  zvals_test = seq(0,1,0.001)
  center_3D = as.matrix(data.frame(a=as.vector(unlist(center_estimated)[1]),
                                   b=as.vector(unlist(center_estimated)[2]),
                                   zvals_test)) %*% solve(t(R))
  center_3D[,1] = center_3D[,1]+MAT$x_rotated_mean[1]
  center_3D[,2] = center_3D[,2]+MAT$y_rotated_mean[1]
  center_3D[,3] = center_3D[,3]+MAT$z_rotated_mean[1]
  spherical = cart2sph(center_3D)
  return(center_3D[which.min(abs(spherical[,3]-1)),])
}


GET_CENTER_CIRCULAR_LACKOFFIT = function(MAT){
  normal_vector = MAT[,c('x_plane', 'y_plane', 'z_plane')][1,]
  ### Perform rodrigues rotation
  P = MAT[,c('x_rotated', 'y_rotated', 'z_rotated')]-MAT[,c('x_rotated_mean', 'y_rotated_mean', 'z_rotated_mean')]
  #Following two runs should give same P_rot
  R = rodrigues_matrix(P, n0=normal_vector, n1 = c(0,0,1))
  P_rot = as.matrix(P) %*% t(R)
  
  fit = pracma::circlefit(xp = P_rot[,1], yp = P_rot[,2])
  
  center_estimated = fit[1:2]
  
  ### Get predictions
  w<- seq(0, 2*pi, len=360)
  xx <- fit[3] * cos(w) + fit[1]
  yy <- fit[3] * sin(w) + fit[2]
  preds_estimated = data.frame(x = xx, y = yy)
  get_dist = function(z,preds_estimated){
    calc_dist = sqrt(((preds_estimated$x - z[1])^2) + ((preds_estimated$y - z[2])^2))
    return(min(calc_dist,na.rm=T))
  }
  DISTS = apply(P_rot[,1:2], 1, FUN = get_dist, preds_estimated)
  Sig = rbind(c(var(DISTS, na.rm=T),0),c(0,var(DISTS, na.rm=T)))
  center_estimated = MASS::mvrnorm(n=1, mu = center_estimated, Sigma = Sig)
  
  zvals_test = seq(0,1,0.001)
  center_3D = as.matrix(data.frame(a=as.vector(unlist(center_estimated)[1]),
                                   b=as.vector(unlist(center_estimated)[2]),
                                   zvals_test)) %*% solve(t(R))
  center_3D[,1] = center_3D[,1]+MAT$x_rotated_mean[1]
  center_3D[,2] = center_3D[,2]+MAT$y_rotated_mean[1]
  center_3D[,3] = center_3D[,3]+MAT$z_rotated_mean[1]
  spherical = cart2sph(center_3D)
  return(center_3D[which.min(abs(spherical[,3]-1)),])
}






GET_CENTER_LACKOFFIT = function(MAT, N = 1){
  normal_vector = MAT[,c('x_plane', 'y_plane', 'z_plane')][1,]
  ### Perform rodrigues rotation
  P = MAT[,c('x_rotated', 'y_rotated', 'z_rotated')]-MAT[,c('x_rotated_mean', 'y_rotated_mean', 'z_rotated_mean')]
  #Following two runs should give same P_rot
  R = rodrigues_matrix(P, n0=as.vector(unlist(normal_vector)), n1 = c(0,0,1))
  P_rot = as.matrix(P) %*% t(R)
  ### Estimate ellipse center
  fit = fitConic::fitConic(X = P_rot[,1], Y = P_rot[,2], conicType = 'e')
  test = fit$parA
  # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2288654/
  A <- matrix(c(2*test[1], test[2], test[2], 2*test[3]), nrow=2, ncol=2, byrow=T )
  b <- matrix(c(-test[4], -test[5]), nrow=2, ncol=1, byrow=T)
  soln <- solve(A) %*% b
  center_estimated = soln[1:2]
  
  ### Get predictions
  preds = fitConic::createConic(x= P_rot[,1], param = fit$parA, conicType = 'e')
  preds_estimated = data.frame(x = preds[,1], y = preds[,2])
  get_dist = function(z,preds_estimated){
    calc_dist = sqrt(((preds_estimated$x - z[1])^2) + ((preds_estimated$y - z[2])^2))
    return(min(calc_dist,na.rm=T))
  }
  DISTS = apply(P_rot[,1:2], 1, FUN = get_dist, preds_estimated)
  Sig = rbind(c(var(DISTS, na.rm=T),0),c(0,var(DISTS, na.rm=T)))
  center_estimated = MASS::mvrnorm(n=1, mu = center_estimated, Sigma = Sig)

  if(is.matrix(center_estimated)){
    center_3D_save = matrix(NA, ncol = 3, nrow = nrow(center_estimated))
  }else{
    center_estimated = rbind(center_estimated)
    center_3D_save = matrix(NA, ncol = 3, nrow = 1)
  }
  zvals_test = seq(0,1,0.001)
  for(i in 1:nrow(center_estimated)){
    center_3D = as.matrix(data.frame(a=as.vector(unlist(center_estimated[i,1])),
                                     b=as.vector(unlist(center_estimated[i,2])),
                                     zvals_test)) %*% solve(t(R))
    center_3D[,1] = center_3D[,1]+MAT$x_rotated_mean[1]
    center_3D[,2] = center_3D[,2]+MAT$y_rotated_mean[1]
    center_3D[,3] = center_3D[,3]+MAT$z_rotated_mean[1]
    spherical = cart2sph(center_3D)
    center_3D_save[i,] = center_3D[which.min(abs(spherical[,3]-1)),]
  }
  return(center_3D_save)
}



##################################
### Perform Rodrigues Rotation ###
##################################
rodrigues_rotation = function(P,n0,n1){
  # Get vector of rotation k and angle theta
  n0 = as.vector(as.numeric(n0/sqrt(sum(n0^2))))
  n1 = as.vector(as.numeric(n1/sqrt(sum(n1^2))))
  k = pracma::cross(n0,n1)
  k = k/sqrt(sum(k^2))
  theta = acos(dot(n0,n1))
  K = matrix(rep(k, length(P[,1])), nrow = length(P[,1]), byrow = T)
  P_rot = P*cos(theta) + cross(K,as.matrix(P))*sin(theta) + K*apply(K*P,1,sum)*(1-cos(theta))
  return(P_rot)
}


rodrigues_matrix = function(P,n0,n1){
  # Get vector of rotation k and angle theta
  n0 = as.vector(as.numeric(n0/sqrt(sum(n0^2))))
  n1 = as.vector(as.numeric(n1/sqrt(sum(n1^2))))
  k = pracma::cross(n0,n1)
  k = k/sqrt(sum(k^2))
  theta = acos(dot(n0,n1))
  K = matrix(0,3,3)
  K[1,2] = -k[3]
  K[1,3] = k[2]
  K[2,1] = k[3]
  K[2,3] = -k[1]
  K[3,1] = -k[2]
  K[3,2] = k[1]
  R = diag(3)+sin(theta)*K + (1-cos(theta))*K %*% K
  return(R)
}




##################################
### Perform Davenport Rotation ###
##################################
#https://en.wikipedia.org/wiki/Davenport_chained_rotations#Tait–Bryan_chained_rotations
davenport_rotation=function(roll=0, pitch, yaw){
  ROLL_MAT = rbind(c(1,0,0), c(0,cos(roll),-sin(roll)), c(0,sin(roll), cos(roll)))
  PITCH_MAT = rbind(c(cos(pitch),0,sin(pitch)),c(0,1,0), c(-sin(pitch), 0, cos(pitch)))
  YAW_MAT = rbind(c(cos(yaw), -sin(yaw), 0),c(sin(yaw), cos(yaw),0), c(0,0,1))
  TRANSFORM_MAT = YAW_MAT %*% PITCH_MAT %*% ROLL_MAT
  return(TRANSFORM_MAT)
}





#####################################
#####################################
### Using Ribbon Separation Peaks ###
#####################################
#####################################
fitsplines_center = function(x,y){
  spline_func = splines::periodicSpline(x,y,period = 360, ord = 4L)
  x2= seq(min(x),max(x),0.01)
  y2 = predict(spline_func,x2)$y
  result = x2[which.max(y2)]
  if(length(result)==0){result = NA}
  return(result)
}

fitsplines_fwqm_center = function(x,y){
  spline_func = splines::periodicSpline(x,obj2=y/max(y,na.rm=T),period = 360, ord = 4L)
  x2= seq(min(x),max(x),0.01)
  y2 = predict(spline_func,x2)$y
  y2[y2<0]=0
  peak_loc = x2[which.max(y2)]
  y_mono = c(rev(cummin(rev(y2[x2 <= peak_loc]))),cummin(y2[x2 > peak_loc]))
  x2 = x2[y_mono >= 0.25]
  result = (max(x2)+min(x2))/2
  if(length(result)==0){result = NA}
  return(result)
}


loesssmooth_center = function(temp, SPAN = 0.3){
  temp = rbind(temp[,c('azimuthal_discrete', 'fits')],
               data.frame(azimuthal_discrete =temp$azimuthal_discrete - 360, fits = temp$fits),
               data.frame(azimuthal_discrete =temp$azimuthal_discrete + 360, fits = temp$fits))
  fit = loess(fits~azimuthal_discrete, data = temp, degree = 2, span = SPAN)
  preds = predict(fit, newdata = data.frame(azimuthal_discrete = AZIMUTHAL))
  return(preds)
}
Estimate_Ribbon_Center = function(model_input, WORKING_CENTER, ndraws = 1000, JACKDRAW = FALSE, EXCLUDE = TRUE, SUMMARY_TYPE = 'Peak', ELLIPSE = TRUE){
  
  if(SUMMARY_TYPE == 'Peak'){
    my_peak_fun = fitsplines_center
  }else if(SUMMARY_TYPE == 'FWQM'){
    my_peak_fun = fitsplines_fwqm_center
  }
  
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  ### Require peaks between these polar angles in working ribbon center frame
  MIN_THRESH = 80
  MAX_THRESH = 130
  K = unique(diff(AZIMUTHAL))
  model_input$include = as.numeric(model_input$polar_discrete >= MIN_THRESH & model_input$polar_discrete <= MAX_THRESH)
  SUBSET = model_input[model_input$include == 1,]
  ### Exclude specified azimuthal angles, following Funsten et al. (2013)
  SUBSET$FunstenExclude1 = as.numeric(SUBSET$azimuthal_discrete >= 126 & SUBSET$azimuthal_discrete <= 192)
  if(EXCLUDE){
    SUBSET = SUBSET[SUBSET$FunstenExclude1 == 0,]
  }
  SUBSET = SUBSET %>% dplyr::group_by(azimuthal_discrete) %>% dplyr::mutate(fits_anchor = my_peak_fun(polar_discrete,ena_rate))
  SUBSET_GLOBAL = SUBSET[!duplicated(SUBSET$azimuthal_discrete),]
  SUBSET_GLOBAL = SUBSET_GLOBAL[order(as.numeric(as.character(SUBSET_GLOBAL$azimuthal_discrete))),]
  SUBSET_GLOBAL$fits_anchor = as.numeric(forecast::tsclean(stats::ts(SUBSET_GLOBAL$fits_anchor)))
  PEAKS_GLOBAL = as.vector(loesssmooth_center(data.frame(azimuthal_discrete=SUBSET_GLOBAL$azimuthal_discrete,
                                                         fits = SUBSET_GLOBAL$fits_anchor), SPAN = 0.1))
  SUBSET = SUBSET[order(SUBSET$azimuthal_discrete, SUBSET$polar_discrete),]

  ### Approximate the covariance
  DATA_DRAW=reshape(data = data.frame(SUBSET[,c('azimuthal_discrete', 'polar_discrete', 'ena_rate')]), idvar = 'polar_discrete', timevar = 'azimuthal_discrete', v.names = c('ena_rate'), direction = 'wide')
  DATA_DRAW = DATA_DRAW[order(DATA_DRAW[,1]),]
  DATA_DRAW = DATA_DRAW[,-1]
  CORMAT_lat = cor(t(DATA_DRAW), method = 'spearman', use = 'pairwise.complete.obs')
  CORMAT_lat[is.na(CORMAT_lat)]=0
  CORMAT_lat = as.matrix(Matrix::nearPD(CORMAT_lat, corr = T)$mat, maxit = 1000)
  CORMAT_lon = cor(DATA_DRAW, method = 'spearman', use = 'pairwise.complete.obs')
  CORMAT_lon[is.na(CORMAT_lon)]=0
  CORMAT_lon = as.matrix(Matrix::nearPD(CORMAT_lon, corr = T)$mat, maxit = 1000)
  Q = Matrix::kronecker(CORMAT_lon, CORMAT_lat)
  SUBSET$sd_ena_rate[is.na(SUBSET$sd_ena_rate)]= min(SUBSET$sd_ena_rate[!is.na(SUBSET$sd_ena_rate)])
  SIGMA = as.matrix(as.vector(unlist(SUBSET$sd_ena_rate))) %*% t(as.matrix(as.vector(unlist(SUBSET$sd_ena_rate))))
  COVAR = Q*SIGMA
  rm('SIGMA')
  COVAR_CHOL = chol(COVAR)
 
  
  
  ######################
  ### Estimate Peaks ###
  ######################
  PEAKS = matrix(NA, nrow = length(AZIMUTHAL), ncol = ndraws)
  for(iter in 1:ndraws){
    #Draw map using provided standard errors
    SUBSET$ENA_RATES = as.vector(mvnfast::rmvn(1, mu = as.vector(unlist(SUBSET$ena_rate)), sigma = COVAR_CHOL, ncores = 1, isChol = TRUE, A = NULL, kpnames = FALSE))
    #Estimate initial peaks using cubic interpolation
    SUBSET = SUBSET %>% dplyr::group_by(azimuthal_discrete) %>% dplyr::mutate(fits = my_peak_fun(polar_discrete,ENA_RATES))
    

    SUBSET2 = SUBSET[!duplicated(SUBSET$azimuthal_discrete),]
    SUBSET2 = SUBSET2[order(as.numeric(as.character(SUBSET2$azimuthal_discrete))),]
    SUBSET2$fits = as.numeric(forecast::tsclean(stats::ts(SUBSET2$fits)))
    PEAKS[,iter] = as.vector(loesssmooth_center(data.frame(SUBSET2[,c('fits', 'azimuthal_discrete')]), SPAN = 0.1))
  }
  
  rm(list = c('COVAR_CHOL', 'COVAR'))

  gc()
  
  PHI_MAX_LONG = data.frame(azimuthal_discrete = rep(AZIMUTHAL, ndraws), 
                            iter = rep(c(1:ndraws), each = length(AZIMUTHAL)), 
                            polar_discrete = as.vector(PEAKS), 
                            id = rep(c(1:length(AZIMUTHAL)), each = ndraws))
  PHI_MAX_LONG = rbind(PHI_MAX_LONG, data.frame(azimuthal_discrete = AZIMUTHAL, 
                                                iter = 0, 
                                                polar_discrete = as.vector(PEAKS_GLOBAL), 
                                                id = c(1:length(AZIMUTHAL))))
  
  #################################
  ### Convert to Ecliptic Frame ###
  #################################
  ### Convert to Cartesian
  # azimuthal_discrete 0-360
  # phi_max 0-180
  PHI_MAX_CART = data.frame(t(apply(PHI_MAX_LONG[,c('azimuthal_discrete', 'polar_discrete')], 1, FUN = function(x){pracma::sph2cart(as.vector(c(as.numeric(x[1]-180)*pi/180,as.numeric(x[2]-90)*pi/180, 1)))})))
  names(PHI_MAX_CART)=c('x', 'y', 'z')
  PHI_MAX_CART$iter = PHI_MAX_LONG$iter

  ### Convert to ecliptic frame
  TRANSFORM_MAT_REVERSE = solve(davenport_rotation(yaw = WORKING_CENTER[1]*pi/180, pitch = -(90-WORKING_CENTER[2])*pi/180))
  temp = as.matrix(PHI_MAX_CART[,c('x','y','z')]) %*% TRANSFORM_MAT_REVERSE
  PHI_MAX_CART$x_rotated = temp[,1]
  PHI_MAX_CART$y_rotated = temp[,2]
  PHI_MAX_CART$z_rotated = temp[,3]
  
  #############################################
  ### Estimate Best Fitting Plane to Points ###
  #############################################
  #https://meshlogic.github.io/posts/jupyter/curve-fitting/fitting-a-circle-to-cluster-of-3d-points/
  
  ### Calculate the mean-centering for each iteration
  CENTER_ITER = aggregate(cbind(x_rotated,y_rotated,z_rotated)~iter, data = PHI_MAX_CART, FUN = mean)
  PHI_MAX_CART = merge(PHI_MAX_CART, data.frame(iter = CENTER_ITER$iter, x_rotated_mean = CENTER_ITER$x_rotated, y_rotated_mean = CENTER_ITER$y_rotated, z_rotated_mean = CENTER_ITER$z_rotated), by = c('iter'), all.x = T)
  
  ### Estimate plane normal vector (limit to positive-z normal)
  PLANE_CART = matrix(NA, ncol = 3, nrow = ndraws)
  for(i in c(1:ndraws)){
    PLANE_CART[i,] =as.vector(unlist(GET_PLANE(MAT = PHI_MAX_CART[PHI_MAX_CART$iter == i,])))
    #print(i)
  }
  PLANE_CART = data.frame(PLANE_CART)
  names(PLANE_CART) = c('x_plane', 'y_plane', 'z_plane')
  PLANE_CART$iter =c(1:ndraws)
  
  PLANE_CART_OVERALL = data.frame(t(unlist(GET_PLANE(MAT = PHI_MAX_CART[PHI_MAX_CART$iter == 0,]))))
  names(PLANE_CART_OVERALL) = c('x_plane', 'y_plane', 'z_plane')
  PLANE_CART_OVERALL$iter =0
  
  ### Merge data
  PHI_MAX_CART = merge(PHI_MAX_CART, rbind(PLANE_CART[,c('iter','x_plane', 'y_plane', 'z_plane')],
                                           PLANE_CART_OVERALL[,c('iter','x_plane', 'y_plane', 'z_plane')]), by = c('iter'))
  
  ######################################################################
  ### Project Data to Best Fitting Plane and Estimate Ellipse Center ###
  ######################################################################
  
  ELLIPSE_CENTERS = matrix(NA, ncol = 3, nrow = ndraws)
  if(ELLIPSE){
  for(i in c(1:ndraws)){
    if(JACKDRAW){
     ELLIPSE_CENTERS[i,] = as.vector(unlist(GET_CENTER_LACKOFFIT(MAT = PHI_MAX_CART[PHI_MAX_CART$iter == i,])))
    }else{
     ELLIPSE_CENTERS[i,] =as.vector(unlist(GET_CENTER(MAT = PHI_MAX_CART[PHI_MAX_CART$iter == i,])))
    }
  }
  ELLIPSE_CENTER_OVERALL = as.vector(unlist(GET_CENTER(MAT = PHI_MAX_CART[PHI_MAX_CART$iter == 0,])))
  }else{
    for(i in c(1:ndraws)){
      if(JACKDRAW){
        ELLIPSE_CENTERS[i,] = as.vector(unlist(GET_CENTER_CIRCULAR_LACKOFFIT(MAT = PHI_MAX_CART[PHI_MAX_CART$iter == i,])))
      }else{
        ELLIPSE_CENTERS[i,] =as.vector(unlist(GET_CENTER_CIRCULAR(MAT = PHI_MAX_CART[PHI_MAX_CART$iter == i,])))
      }
    }
    ELLIPSE_CENTER_OVERALL = as.vector(unlist(GET_CENTER_CIRCULAR(MAT = PHI_MAX_CART[PHI_MAX_CART$iter == 0,])))
  }
  
  #test = apply(as.matrix(c(1:ndraws)), 1, FUN = function(x,MAT){as.vector(unlist(GET_CENTER(MAT = MAT[MAT$iter == x,])))}, MAT = PHI_MAX_CART)
  
  ELLIPSE_CENTERS = data.frame(ELLIPSE_CENTERS)
  names(ELLIPSE_CENTERS) = c('x_center', 'y_center', 'z_center')
  
  ELLIPSE_CENTER_OVERALL = data.frame(t(ELLIPSE_CENTER_OVERALL))
  names(ELLIPSE_CENTER_OVERALL) = c('x_center', 'y_center', 'z_center')
  
  
  ##########################
  ### Convert to Lat/Lon ###
  ##########################
  ELLIPSE_CENTERS_SPHERICAL = data.frame(cart2sph(as.matrix(ELLIPSE_CENTERS)))
  ELLIPSE_CENTERS_SPHERICAL$lon = ELLIPSE_CENTERS_SPHERICAL[,1]*180/pi + 180
  ELLIPSE_CENTERS_SPHERICAL$ecliptic_lat= (ELLIPSE_CENTERS_SPHERICAL[,2])*180/pi
  
  ELLIPSE_CENTERS_SPHERICAL_OVERALL = data.frame(t(cart2sph(as.matrix(ELLIPSE_CENTER_OVERALL))))
  ELLIPSE_CENTERS_SPHERICAL_OVERALL$lon = ELLIPSE_CENTERS_SPHERICAL_OVERALL[,1]*180/pi + 180
  ELLIPSE_CENTERS_SPHERICAL_OVERALL$ecliptic_lat= (ELLIPSE_CENTERS_SPHERICAL_OVERALL[,2])*180/pi
  
  PEAKS_SPHERICAL = data.frame(cart2sph(as.matrix(PHI_MAX_CART[,c('x_rotated', 'y_rotated', 'z_rotated')])))
  PEAKS_SPHERICAL$lon = PEAKS_SPHERICAL[,1]*180/pi + 180
  PEAKS_SPHERICAL$ecliptic_lat= (PEAKS_SPHERICAL[,2])*180/pi
  PEAKS_SPHERICAL$iter = PHI_MAX_CART$iter
  
  return(list(ELLIPSE = ELLIPSE_CENTERS_SPHERICAL[,c('lon', 'ecliptic_lat')], OVERALL = ELLIPSE_CENTERS_SPHERICAL_OVERALL[,c('lon', 'ecliptic_lat')],
              PEAKS = PEAKS_SPHERICAL[,c('lon', 'ecliptic_lat', 'iter')]))
}


 

########################################
### Funsten Ribbon Center Estimation ###
########################################

LeastSquares = function(param, DAT){
  A = exp(param[1])
  B = exp(param[2])
  phi_max = param[3]
  s = param[4]
  Observed = DAT$ENA_RATES
  Model = A + B*exp(-(DAT$polar_discrete-phi_max)^2/(2*exp(s)))
  return(sum((Observed-Model)^2)*1000)
} 

Funsten_RibbonEst = function(x, AZIMUTHAL, TEMP){
  cur_index = which(AZIMUTHAL == x)
  if(cur_index == 1){
    include_azi = c(AZIMUTHAL[1], AZIMUTHAL[2], AZIMUTHAL[length(AZIMUTHAL)])
  }else if(cur_index == length(AZIMUTHAL)){
    include_azi = c(AZIMUTHAL[length(AZIMUTHAL)], AZIMUTHAL[length(AZIMUTHAL)-1], AZIMUTHAL[1])
  }else{
    include_azi = c(AZIMUTHAL[cur_index],AZIMUTHAL[cur_index+1], AZIMUTHAL[cur_index-1])
  }
  DAT = TEMP[TEMP$azimuthal_discrete %in% include_azi,]
  param = c(log(min(DAT$ENA_RATES)), log(max(DAT$ENA_RATES)-min(DAT$ENA_RATES)), mean(DAT$phi_max_discrete),3 )
  if(!is.null(DAT[1,'par1'])){
    if(!is.na(DAT[1,'par1'])){
      param = as.vector(unlist(DAT[1,c('par1', 'par2', 'par3', 'par4')]))
    }
  }
  fit = optimx(par = param, fn = LeastSquares, DAT = DAT, method = 'newuoa')
  hess_val = try(solve(gHgen(par = c(fit$p1, fit$p2, fit$p3, fit$p4), LeastSquares, gr=NULL, hess=NULL,
                             control=list(ktrace=0), DAT = DAT)$Hn)[3,3])
  return(c(fit$p1,fit$p2,ifelse(is.na(fit$p3), DAT$phi_max_discrete[1],fit$p3),fit$p4,ifelse(class(hess_val) == 'try-error', 100, hess_val)))
}

Estimate_Funsten = function(model_input, ESA, WORKING_CENTER, ndraws, ELLIPSE = TRUE, EXCLUDE = T){
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  PHI_MAX_LONG = data.frame(azimuthal_discrete =AZIMUTHAL, iter = 1, polar_discrete = NA,  id = c(1:length(AZIMUTHAL)), BOVERA = NA, HESS = NA)
  MIN_THRESH = ifelse(ESA == 6, 78,90) 
  MAX_THRESH = ifelse(ESA == 6, 120,120)
  K = diff(AZIMUTHAL)[1]
  
  ### Fit the Gaussian profiles
  TEMP = data.frame(model_input, ENA_RATES=model_input$ena_rate)[model_input$polar_discrete >= MIN_THRESH & model_input$polar_discrete <= MAX_THRESH,]    
  temp = aggregate(ENA_RATES~azimuthal_discrete,data = TEMP, FUN = function(x,POLAR){POLAR[which.max(x)]},POLAR = POLAR[POLAR>=MIN_THRESH & POLAR <= MAX_THRESH])
  TEMP = merge(TEMP, data.frame(azimuthal_discrete = temp$azimuthal_discrete, phi_max_discrete = temp$ENA_RATES), by = 'azimuthal_discrete', all.x = T)

  VALUES = apply(as.matrix(unique(TEMP$azimuthal_discrete)),1,FUN = Funsten_RibbonEst, AZIMUTHAL, TEMP=TEMP[,c('polar_discrete', 'azimuthal_discrete', 'ENA_RATES', 'phi_max_discrete')])
  PHI_MAX_LONG$polar_discrete = VALUES[3,]
  PHI_MAX_LONG$BOVERA = exp(VALUES[2,])/exp(VALUES[1,])
  PHI_MAX_LONG$HESS = VALUES[5,]
  
  ### Additional Exclusions from Funsten
  if(EXCLUDE){
    PHI_MAX_LONG$FunstenExclude1 = as.numeric(PHI_MAX_LONG$azimuthal_discrete >= 126 & PHI_MAX_LONG$azimuthal_discrete <= 192)
  }else(
    PHI_MAX_LONG$FunstenExclude1 = 0 
  )
  PHI_MAX_LONG$FunstenExclude2 = as.numeric(PHI_MAX_LONG$BOVERA <= 0.05)
  PHI_MAX_LONG$FunstenExclude2[is.na(PHI_MAX_LONG$FunstenExclude2)]=1
  PHI_MAX_LONG$FunstenExclude2[is.na(PHI_MAX_LONG$BOVERA)]=1
  
  PHI_MAX_LONG = PHI_MAX_LONG[!is.na(PHI_MAX_LONG$HESS) & PHI_MAX_LONG$HESS > 0,]
  PHI_MAX_LONG$SE = sqrt(PHI_MAX_LONG$HESS) #note: already inverted!!
  ERROR_THRESH = ifelse(ESA == 6, 10,10) ### Modified from Funsten et al.!!!!

  PHI_MAX_LONG$FunstenExclude3 = as.numeric(PHI_MAX_LONG$SE>=ERROR_THRESH)
  PHI_MAX_LONG$FunstenExclude3[is.na(PHI_MAX_LONG$FunstenExclude3) | is.na(PHI_MAX_LONG$SE)]=1
  PHI_MAX_LONG$FunstenExclude = as.numeric(PHI_MAX_LONG$FunstenExclude1|PHI_MAX_LONG$FunstenExclude2|PHI_MAX_LONG$FunstenExclude3)
  PHI_MAX_LONG = PHI_MAX_LONG[PHI_MAX_LONG$FunstenExclude == 0,]

  ############################
  ### Convert to Cartesian ###
  ############################
  PHI_PLANE = as.matrix(PHI_MAX_LONG[,c('azimuthal_discrete', 'polar_discrete')])
  PHI_PLANE[,1] = (PHI_PLANE[,1]-180)*pi/180
  PHI_PLANE = pracma::pol2cart(PHI_PLANE)

  ################################
  ### Fit Ellipse in Cartesian ###
  ################################
  # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2288654/
  
  PHI_PLANE = data.frame(PHI_PLANE)
  SCALING = max(c(sd(PHI_PLANE$x), sd(PHI_PLANE$y)))
  PHI_PLANE$x_stand = PHI_PLANE$x/SCALING
  PHI_PLANE$y_stand = PHI_PLANE$y/SCALING

  if(ELLIPSE){
  ### Get center
  fit = fitConic::fitConic(X = PHI_PLANE$x_stand, Y = PHI_PLANE$y_stand, conicType = 'e')
  test = fit$parA
  A <- matrix(c(2*test[1], test[2], test[2], 2*test[3]), nrow=2, ncol=2, byrow=T )
  b <- matrix(c(-test[4], -test[5]), nrow=2, ncol=1, byrow=T)
  soln <- solve(A) %*% b
  center_estimated = soln[1:2]
  }else{
    fit = pracma::circlefit(xp = PHI_PLANE$x_stand, yp = PHI_PLANE$y_stand)
    center_estimated = fit[1:2]
  }
  
  ### Convert center back to lat/lon
  center_estimated = data.frame(x = center_estimated[1], y = center_estimated[2])
  center_estimated$x = center_estimated$x*SCALING
  center_estimated$y = center_estimated$y*SCALING
  center_estimated = data.frame(t(cart2pol(as.matrix(center_estimated[,c('x','y')]) )))
  names(center_estimated) = c('phi','r')
  center_estimated$polar_discrete = center_estimated$r
  center_estimated$azimuthal_discrete = center_estimated$phi*180/pi + 180
  if(ELLIPSE){
    preds = fitConic::createConic(x= PHI_PLANE$x_stand, param = fit$parA, conicType = 'e')
    preds_estimated = data.frame(x = preds[,1], y = preds[,2])
  }else{
    w<- seq(0, 2*pi, len=360)
    xx <- fit[3] * cos(w) + fit[1]
    yy <- fit[3] * sin(w) + fit[2]
    preds_estimated = data.frame(x = xx, y = yy)
  }
  preds_estimated$x = preds_estimated$x*SCALING
  preds_estimated$y = preds_estimated$y*SCALING
  preds_estimated = data.frame(cart2pol(as.matrix(preds_estimated[,c('x','y')]) ))
  names(preds_estimated) = c('phi','r')
  preds_estimated$polar_discrete = preds_estimated$r
  preds_estimated$azimuthal_discrete = preds_estimated$phi*180/pi + 180
  
  ### Match up fitted values
  PHI_MAX_LONG = PHI_MAX_LONG %>% dplyr::group_by(azimuthal_discrete) %>%
    dplyr::mutate(fitted_preds = preds_estimated$polar_discrete[which.min(abs(preds_estimated$azimuthal_discrete -azimuthal_discrete))])
  PHI_MAX_LONG$error = PHI_MAX_LONG$polar_discrete - PHI_MAX_LONG$fitted_preds

  ## Draw centers
  var_est =  var(PHI_MAX_LONG$error, na.rm=T) 
  center_drawn  = MASS::mvrnorm(n = ndraws, mu = c(center_estimated$polar_discrete, center_estimated$azimuthal_discrete), 
                                Sigma = rbind(c(var_est,0), c(0,var_est)))
  center_drawn = data.frame(center_drawn)
  names(center_drawn) = c('polar_discrete', 'azimuthal_discrete')

  #################################
  ### Convert to Ecliptic Frame ###
  #################################
  ### Convert to Cartesian
  ### + 90 deals with issues with polar coordinate systems
  TO_ROTATE = rbind(center_estimated[,c('azimuthal_discrete', 'polar_discrete')], center_drawn[,c('azimuthal_discrete', 'polar_discrete')])
  PHI_MAX_CART = data.frame(t(apply(TO_ROTATE, 1, FUN = function(x){pracma::sph2cart(as.vector(c(as.numeric(x[1]-180)*pi/180,as.numeric(x[2] + 90)*pi/180, 1)))})))
  names(PHI_MAX_CART)=c('x', 'y', 'z')
  
  ### Convert to ecliptic frame
  TRANSFORM_MAT_REVERSE = solve(davenport_rotation(yaw = WORKING_CENTER[1]*pi/180, pitch = -(90-WORKING_CENTER[2])*pi/180))
  temp = as.matrix(PHI_MAX_CART[,c('x','y','z')]) %*% TRANSFORM_MAT_REVERSE
  PHI_MAX_CART$x_rotated = temp[,1]
  PHI_MAX_CART$y_rotated = temp[,2]
  PHI_MAX_CART$z_rotated = temp[,3]
  
  
  ##########################
  ### Convert to Lat/Lon ###
  ##########################
  ELLIPSE_CENTERS_SPHERICAL = data.frame(cart2sph(as.matrix(PHI_MAX_CART[,c('x_rotated', 'y_rotated', 'z_rotated')])))
  ELLIPSE_CENTERS_SPHERICAL$lon = ELLIPSE_CENTERS_SPHERICAL[,1]*180/pi + 180
  ELLIPSE_CENTERS_SPHERICAL$ecliptic_lat= (ELLIPSE_CENTERS_SPHERICAL[,2])*180/pi
  
  
  return(list(ELLIPSE = ELLIPSE_CENTERS_SPHERICAL[,c('lon', 'ecliptic_lat')]))
}



