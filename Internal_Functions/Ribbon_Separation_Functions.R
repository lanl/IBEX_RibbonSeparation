###################################
###################################
### Ribbon Separation Functions ###
###################################
###################################
### Developed by Dr. Lauren J. Beesley, PhD
### Contact: lvandervort@lanl.gov
### Last Updated: 08/11/2023

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

Get_Loess_PHI_ANCHOR = function(temp, span = 0.1){
  temp = rbind(temp[,c('azimuthal_discrete', 'phi_max')],
               data.frame(azimuthal_discrete =temp$azimuthal_discrete - 360, phi_max = temp$phi_max),
               data.frame(azimuthal_discrete =temp$azimuthal_discrete + 360, phi_max = temp$phi_max))
  fit = loess(phi_max~azimuthal_discrete, data = temp, degree = 2, span = span) 
  preds = predict(fit, newdata = data.frame(azimuthal_discrete = AZIMUTHAL))
  return(preds)
}

fitsplines_center_restricted = function(x,y, possible){
  x = x[!is.na(y) & !is.nan(y) & !is.infinite(y)]
  y = y[!is.na(y) & !is.nan(y) & !is.infinite(y)]
  spline_func = splines::periodicSpline(x,y,period = 360, ord = 4L)
  xnew= seq(min(x),max(x),0.01)
  xnew = xnew[xnew >= min(possible, na.rm=T) & xnew <= max(possible, na.rm=T)]
  z = predict(spline_func,xnew)$y
  peak_est = pracma::findpeaks(z, npeaks = 1)[,1:2]
  result = c(xnew[peak_est[2]],peak_est[1],max(y,na.rm=T))
  if(length(peak_est)==0){result = rep(NA,3)}
  result[is.infinite(result)]=NA
  return(result)
}

evalsplines_center_restricted = function(x,y, xnew){
  x = x[!is.na(y) & !is.nan(y) & !is.infinite(y)]
  y = y[!is.na(y) & !is.nan(y) & !is.infinite(y)]
  spline_func = splines::periodicSpline(x,y,period = 360, ord = 4L)
  z = predict(spline_func,xnew)$y
  return(z)
}

Get_Distance_Features = function(model_input, RIBBON_CENTER){
  model_input = data.frame(model_input)
  
  ### Convert to ecliptic
  TEMP= Coord_Transform_Points(model_input[,c('lon','lat')], RIBBON_CENTER = RIBBON_CENTER, REVERSE = TRUE)
  model_input$lon_ecliptic = TEMP$lon_new
  model_input$lat_ecliptic = TEMP$lat_new
  
  ### Get cartesian coordinates
  tail = c(-104.2+180,5.16)
  nose = c(75.8+180,-5.16)
  PHI_MAX_CART = data.frame(t(apply(rbind(tail, nose), 1, FUN = function(x){pracma::sph2cart(as.vector(c(as.numeric(x[1]-180)*pi/180,as.numeric(x[2])*pi/180, 1)))})))
  names(PHI_MAX_CART)=c('x_rotated', 'y_rotated', 'z_rotated')
  
  ### Calculate a point with latitude = 0 that is midway between the nose and tail. Technically, there are two,
  ### but will not matter which one we choose. Other possible value comes from negative sqrt. 
  mid_point = apply(PHI_MAX_CART,2,mean)
  mid_point = mid_point / sqrt(sum(mid_point^2))
  PHI_MAX_CART = rbind(PHI_MAX_CART, mid_point)
  PHI_MAX_CART$iter = 1
  
  ### Calculate best fitting plane between points
  CENTER_ITER = aggregate(cbind(x_rotated,y_rotated,z_rotated)~iter, data = PHI_MAX_CART, FUN = mean)
  PHI_MAX_CART = merge(PHI_MAX_CART, data.frame(iter = CENTER_ITER$iter, x_rotated_mean = CENTER_ITER$x_rotated, y_rotated_mean = CENTER_ITER$y_rotated, z_rotated_mean = CENTER_ITER$z_rotated), by = c('iter'), all.x = T)
  PLANE_CART = matrix(NA, ncol = 3, nrow = 1)
  PLANE_CART =as.vector(unlist(GET_PLANE(MAT = PHI_MAX_CART)))
  PLANE_CART = data.frame(rbind(PLANE_CART))
  names(PLANE_CART) = c('x_plane', 'y_plane', 'z_plane')
  PLANE_CART$iter = 1
  PHI_MAX_CART = merge(PHI_MAX_CART, PLANE_CART[,c('iter','x_plane', 'y_plane', 'z_plane')], by = c('iter'))
  
  ### Calculate location of new north pole
  ELLIPSE_CENTER=as.vector(unlist(GET_CENTER_CIRCULAR(MAT = PHI_MAX_CART)))
  ELLIPSE_CENTERS_SPHERICAL = data.frame(rbind(cart2sph(as.matrix(PLANE_CART[1:3]))))
  ELLIPSE_CENTERS_SPHERICAL$lon = ELLIPSE_CENTERS_SPHERICAL[,1]*180/pi + 180
  ELLIPSE_CENTERS_SPHERICAL$ecliptic_lat= (ELLIPSE_CENTERS_SPHERICAL[,2])*180/pi

  ### Rotate into new coordinate system, where midpoint doesn't move
  TEST = data.frame(model_input[,c('lon_ecliptic', 'lat_ecliptic')])
  names(TEST) = c('lon','lat')
  TESTNEW = Coord_Transform_Points(TEST, RIBBON_CENTER = c(ELLIPSE_CENTERS_SPHERICAL$lon, ELLIPSE_CENTERS_SPHERICAL$ecliptic_lat), REVERSE = FALSE)
  model_input$lon_nosetail = TESTNEW$lon_new
  model_input$lat_nosetail = TESTNEW$lat_new

  return(model_input)
}








GAM_Separation = function(RIBBON_CENTER = NULL, 
                          THRESH = 0.75,  BUFFER = 2, RIBBON_BOUNDS = 30, KNOTS = 10, 
                          PEAK_ALL = NULL, AZI_KNOTS = 24, GAMK = 50, THESEUS_DELINE = F,
                          PTERMS = 10){
  #######################################
  ### Calculate Cartesian Coordinates ###
  #######################################
  if('z' %in% names(model_input)){model_input = subset(model_input, select = -c(x,y,z))}
  MODEL_CART = data.frame(t(apply(model_input[,c('azimuthal_discrete', 'polar_discrete')], 1, FUN = function(x){pracma::sph2cart(as.vector(c(as.numeric(x[1]-180)*pi/180,as.numeric(x[2]-90)*pi/180, 1)))})))
  names(MODEL_CART)=c('x', 'y', 'z')
  model_input = data.frame(model_input, MODEL_CART)
  K = unique(diff(AZIMUTHAL))
  model_input$ena_rate[is.nan(model_input$ena_rate)]=NA
  
  ###########################
  ### Get Initial Phi_Max ###
  ###########################
  if('phi_max_esa4' %in% names(model_input)){model_input = subset(model_input, select = -c(phi_max_esa4))}
  model_input = merge(model_input, data.frame(azimuthal_discrete = AZIMUTHAL, phi_max_esa4 = PEAK_ALL$lat_new + 90), by = 'azimuthal_discrete', all.x = T, all.y = T)
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  if(K==6){
    model_input$eligible_angles = as.numeric(model_input$polar_discrete <= (model_input$phi_max_esa4+24) & model_input$polar_discrete >=  (model_input$phi_max_esa4-24) & model_input$polar_discrete >= 90  )
  }else{
    model_input$eligible_angles = as.numeric(model_input$polar_discrete <= (model_input$phi_max_esa4+10) & model_input$polar_discrete >=  (model_input$phi_max_esa4-10) & model_input$polar_discrete >= 90 )
  }
  SUBSET = model_input
  SUBSET = SUBSET %>% dplyr::group_by(azimuthal_discrete) %>% 
    dplyr::mutate(phi_max = fitsplines_center_restricted(polar_discrete, ena_rate, polar_discrete[eligible_angles==1])[1])
  SUBSET = SUBSET[!duplicated(SUBSET$azimuthal_discrete),]
  SUBSET = data.frame(SUBSET[order(as.numeric(as.character(SUBSET$azimuthal_discrete))),])
  SUBSET$phi_max[c(0,abs(diff(SUBSET$phi_max)))> 12 | is.na(SUBSET$phi_max)]=NA#don't allow a greater than 12 degree jump between neighboring pixels. 
  SUBSET$phi_max = as.numeric(forecast::tsclean(ts(SUBSET$phi_max)))
  PREDS_ANCHOR = Get_Loess_PHI_ANCHOR(SUBSET, span = 0.1)
  PHI_ANCHOR <<- as.vector(PREDS_ANCHOR)
  if('phi_max' %in% names(model_input)){
    model_input = subset(model_input, select = -c(phi_max))
  }
  model_input = merge(model_input,data.frame(azimuthal_discrete = AZIMUTHAL, 
                                             phi_max=  PHI_ANCHOR), by = 'azimuthal_discrete',
                      all.x = T, all.y = F)
  model_input$dist_from_peak =model_input$polar_discrete- model_input$phi_max
  model_input$phi_max_orig = model_input$phi_max
  
  ######################################
  ### Define the Primary Ribbon Mask ###
  ######################################
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  
  ### Edge detection and inclusion/exclusions
  if(K==2){
    BREAKS =seq(-181, 181, 2)
  }else{
    BREAKS =seq(-183, 183, 6)
  }
  model_input$dist_from_peak_int = as.numeric(as.character(cut(model_input$dist_from_peak, breaks=  BREAKS, right = FALSE, labels = BREAKS[-length(BREAKS)]+1)))
  BREAKS_LABELS = sort(unique(model_input$dist_from_peak_int ))
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$dist_from_peak_int),]
  model_input = model_input %>% dplyr::group_by(polar_discrete) %>% dplyr::mutate(ena_smoothed = predict(loess(ena_rate~azimuthal_discrete, span = 0.5), newdata = data.frame(azimuthal_discrete = AZIMUTHAL)))
  model_input$dist_from_peak_int = as.numeric(as.character(model_input$dist_from_peak_int))
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$dist_from_peak_int),]
  
  DATA_DRAW=reshape(data = data.frame(model_input[,c('azimuthal_discrete', 'dist_from_peak_int', 'ena_smoothed')]), idvar = 'dist_from_peak_int', timevar = 'azimuthal_discrete', v.names = c('ena_smoothed'), direction = 'wide')
  DATA_DRAW = DATA_DRAW[order(DATA_DRAW[,1]),]
  DATA_DRAW = as.matrix(DATA_DRAW[,-1])
  DATA_DRAW[is.na(DATA_DRAW)] = 0
  kernel1 <- matrix(c(1,2,1, 0,  0,  0,-1,-2,-1), nrow = 3, byrow = T)
  convolutionExample1  <- imagine::convolution2D(X = t(DATA_DRAW), kernel = kernel1)
  convolutionExample2  <- imagine::convolution2D(X = t(DATA_DRAW), kernel = t(kernel1))
  im_laplace = sqrt(convolutionExample1^2 + convolutionExample2^2)
  im_laplace_sub = im_laplace
  im_laplace_sub[1,] = im_laplace_sub[2,]
  im_laplace_sub[length(im_laplace_sub[,1]),] = im_laplace_sub[length(im_laplace_sub[,1])-1,]
  im_laplace_sub[,abs(BREAKS_LABELS)>RIBBON_BOUNDS] = 0
  im_laplace_save = im_laplace_sub
  im_laplace_sub[im_laplace_sub<quantile(im_laplace_sub[im_laplace_sub>0], THRESH)]=0
  im_laplace_save2 = im_laplace_sub
  if(sum(im_laplace_sub)==0){
    im_laplace_sub[,abs(BREAKS_LABELS)<=RIBBON_BOUNDS]=1
  }
  
  ### Define conservative mask band
  UPPERS = apply(t(im_laplace_sub),2,FUN = function(x,BREAKS_LABELS){max(BREAKS_LABELS[x>0],na.rm=T)},BREAKS_LABELS)
  UPPERS[is.infinite(UPPERS)] = max(UPPERS[!is.infinite(UPPERS)])
  LOWERS = apply(t(im_laplace_sub),2,FUN = function(x,BREAKS_LABELS){min(BREAKS_LABELS[x>0],na.rm=T)},BREAKS_LABELS)
  LOWERS[is.infinite(LOWERS)] = min(LOWERS[!is.infinite(LOWERS)])
  if('upper_cut' %in% names(model_input)){model_input = subset(model_input, select = -c(upper_cut,lower_cut))}
  model_input = merge(model_input, data.frame(azimuthal_discrete = AZIMUTHAL,
                                              upper_cut = max(UPPERS,na.rm=T)+BUFFER*unique(diff(AZIMUTHAL)), lower_cut = min(LOWERS,na.rm=T)-BUFFER*unique(diff(AZIMUTHAL))),
                      by = 'azimuthal_discrete', all.x = T, all.y = T)
  model_input$mask = as.numeric(model_input$dist_from_peak_int <= model_input$upper_cut &
                                  model_input$dist_from_peak_int >= model_input$lower_cut)
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  
  
  ###############################
  ### Initial Estimate of GDF ### 
  ###############################
  
  ### Get features related to nose/tail BV plane
  mybasis <- "sos" 
  model_input = Get_Distance_Features(model_input, RIBBON_CENTER = RIBBON_CENTER)
  model_input$abs_lat_nosetail = abs(model_input$lat_nosetail)
  model_input$ind_above = as.numeric(model_input$lat_nosetail>=0)
  model_input = model_input %>% dplyr::group_by(azimuthal_discrete) %>% dplyr::mutate(minval = min(ena_rate[mask == 1 & !is.na(ena_rate)], na.rm=T))
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  library(mgcv) #guards against dpnorm issue for cnorm function---function name conflict issue
  model_input$lower = ifelse(model_input$mask == 0, model_input$ena_rate, model_input$minval)/100
  model_input$upper = model_input$ena_rate/100
  model_input$GAMK = GAMK
  #Divide by 100 due to cnorm bug that dies when y is too large
  
  ### Model 1: Symmetry
  TEMP = rbind(data.frame(model_input,ind = 1),data.frame(model_input,ind = 0))
  TEMP$abs_lat_nosetail[TEMP$ind==0] = -TEMP$abs_lat_nosetail[TEMP$ind==0]
  TEMP$yi = as.matrix(cbind(TEMP$lower,TEMP$upper))
  TEMP = TEMP[!is.na(TEMP$ena_rate),]
  gam_mod2 <- mgcv::bam(yi ~ s(abs_lat_nosetail, lon_nosetail, bs=mybasis, m= 4, k = GAMK),
                        data=TEMP,
                        family = mgcv::cnorm, discrete= T) #fit using mask, too
  PRED2 <- predict(gam_mod2, newdata = model_input, type="response", se.fit = T)
  model_input$init2 <- PRED2$fit*100
  model_input$init2_se <- PRED2$se.fit*100
  
  ### Model 2: No Symmetry
  gam_mod3b <- mgcv::bam(yi ~ s(lat_nosetail, lon_nosetail, bs=mybasis, m= 4, k = GAMK),
                         data=TEMP,
                         family = mgcv::cnorm, discrete= T) #fit using mask, too
  PRED3 <- predict(gam_mod3b, newdata = model_input, type="response", se.fit = T)
  model_input$init3 <- PRED3$fit*100
  model_input$init3_se <- PRED3$se.fit*100
  
  ### Combine: weighting as function of abs(lon_nosetail - 180)
  model_input$dist = as.numeric(abs(model_input$lon_nosetail-180))
  model_input$dist = model_input$dist / max(model_input$dist)
  
  model_input$gdf_init = as.numeric(model_input$init2)*(1-model_input$dist)+
    as.numeric(model_input$init3)*(model_input$dist)
  model_input$se.fit = sqrt( as.numeric(model_input$init2_se^2)*((1-model_input$dist)^2)+
                               as.numeric(model_input$init3_se^2)*(model_input$dist^2)+
                               2*model_input$dist*(1-model_input$dist)*as.numeric(model_input$init2_se)*as.numeric(model_input$init3_se))
  #Var(aX+bY)=a2Var(X)+b2Var(Y)+2abCov(X,Y) <= a2Var(X)+b2Var(Y)+2ab*SE(X)*SE(Y), using correlation = 1
  Xp3 <- model.matrix(gam_mod3b)
  vcovpred_lower = Xp3[TEMP$ind == 1,] %*% t(chol(gam_mod3b$Vp))
  rm(list = c('Xp3'))
  
  
  
  ################################
  ### Re-estimate Ribbon Peaks ###
  ################################
  
  model_input$ribbon_init = model_input$ena_rate - model_input$gdf_init
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  SUBSET = model_input
  if(K==6){
    SUBSET$eligible_angles = as.numeric(model_input$polar_discrete <= (model_input$phi_max_esa4+24) & model_input$polar_discrete >=  (model_input$phi_max_esa4-24) & model_input$polar_discrete >= 90  )
  }else{
    SUBSET$eligible_angles = as.numeric(model_input$polar_discrete <= (model_input$phi_max_esa4+10) & model_input$polar_discrete >=  (model_input$phi_max_esa4-10) & model_input$polar_discrete >= 90 )
  }
  SUBSET = SUBSET %>% dplyr::group_by(azimuthal_discrete) %>% 
    dplyr::mutate(phi_max = fitsplines_center_restricted(polar_discrete, ribbon_init, polar_discrete[eligible_angles==1])[1])
  SUBSET = SUBSET %>% dplyr::group_by(azimuthal_discrete) %>%
    dplyr::mutate(min_eligible = min(polar_discrete[eligible_angles==1], na.rm=T))
  SUBSET_SAVE = SUBSET
  SUBSET = SUBSET[!duplicated(SUBSET$azimuthal_discrete),]
  
  ### Clean and Smooth
  if(K == 6){
    SUBSET$phi_max[c(0,abs(diff(SUBSET$phi_max)))> 12 | is.na(SUBSET$phi_max)]=NA#don't allow a greater than 12 degree jump between neighboring pixels. 
    SUBSET$phi_max[abs(SUBSET$phi_max - SUBSET$phi_max_orig) > 12 | is.na(SUBSET$phi_max) | is.na(SUBSET$phi_max_orig)]=NA
  }else{
    SUBSET$phi_max[c(0,abs(diff(SUBSET$phi_max)))> 4 | is.na(SUBSET$phi_max)]=NA#don't allow a greater than 12 degree jump between neighboring pixels. 
    SUBSET$phi_max[abs(SUBSET$phi_max - SUBSET$phi_max_orig) > 4 | is.na(SUBSET$phi_max) | is.na(SUBSET$phi_max_orig)]=NA
  }
  SUBSET = data.frame(SUBSET[order(as.numeric(as.character(SUBSET$azimuthal_discrete))),])
  SUBSET$phi_max = as.numeric(forecast::tsclean(ts(SUBSET$phi_max)))
  PREDS_ANCHOR = Get_Loess_PHI_ANCHOR(SUBSET, span = 0.1)
  PHI_ANCHOR <<- as.vector(PREDS_ANCHOR)
  SUBSET$phi_max = PHI_ANCHOR
  
  ### Merge
  if('phi_max' %in% names(model_input)){model_input = subset(model_input, select = -c(phi_max))}
  model_input = merge(model_input,data.frame(azimuthal_discrete = AZIMUTHAL, 
                                             phi_max=  PHI_ANCHOR), by = 'azimuthal_discrete',
                      all.x = T, all.y = F)
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  model_input$dist_from_peak = model_input$polar_discrete - model_input$phi_max
  
  
  ##############################################
  ### Re-masking, Only if THESEUS_DELINE = T ###
  ##############################################
  if(THESEUS_DELINE){
    model_input$ribbon_init = model_input$ena_rate - model_input$gdf_init
    model_input =  model_input %>% dplyr::group_by(azimuthal_discrete) %>%
      dplyr::mutate(ribbon_init0= c(rev(cummin(rev(ribbon_init[dist_from_peak <= 0]))),cummin(ribbon_init[dist_from_peak > 0])))
    #fit Gaussian to each vertical stripe
    fit_Gaussian = function(x, model_input){
      DAT = model_input[model_input$azimuthal_discrete == x,]
      LeastSquares = function(param, DAT){
        B = exp(param[1])
        s1 = param[2]
        s2 = param[3]
        Model = B*exp(-(DAT$polar_discrete-DAT$phi_max)^2/(2*exp(s1)))*as.numeric(DAT$polar_discrete >= DAT$phi_max)
        Model = Model + B*exp(-(DAT$polar_discrete-DAT$phi_max)^2/(2*exp(s2)))*as.numeric(DAT$polar_discrete < DAT$phi_max)
        return(sum((DAT$ribbon_init0-Model)^2,na.rm=T)*1000)
      }
      param = c(min(DAT$ribbon_init0,na.rm=T),3,3 )
      fit = try(optimx(par = param, fn = LeastSquares, DAT = DAT, method = 'newuoa'))
      if(class(fit)[1]!='try-error'){
        return(c(exp(fit$p1), fit$p2, fit$p3))
      }else{
        return(c(NA,NA))
      }
    }
    library(optimx)
    GAUSS = apply(cbind(AZIMUTHAL),1,fit_Gaussian,model_input) 
    SD_CUTOFF_UPPER = sqrt(median(exp(GAUSS[2,]),na.rm=T))
    SD_CUTOFF_LOWER = sqrt(median(exp(GAUSS[3,]),na.rm=T))
    
    model_input$mask_temp = as.numeric(model_input$dist_from_peak_int <= (2*SD_CUTOFF_UPPER + BUFFER*unique(diff(AZIMUTHAL))) & 
                                         model_input$dist_from_peak_int >= (-2*SD_CUTOFF_LOWER -BUFFER*unique(diff(AZIMUTHAL))) ) 
    model_input$mask_temp[model_input$mask_temp == 1 & model_input$mask == 0]=0
    model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  }else{
    model_input$mask_temp = model_input$mask
  }
  
  
  ##################################
  ### Implement Delining via PPR ### 
  ##################################
  if(THESEUS_DELINE){
    model_input$ribbon_init = model_input$ena_rate - model_input$gdf_init
    if(PTERMS > 0){
      ppr_resid_mod <- ppr(ribbon_init ~ x+y+z,
                         nterms = PTERMS,sm.method = 'gcvspline',
                         #weights = ppr_wt,
                         data = model_input[model_input$mask_temp == 0,])
      model_input$ppr_gdf_bias <- predict(ppr_resid_mod, newdata=model_input, type="response")
    }else{
      model_input$ppr_gdf_bias = 0 
    }
    ### don't allow prediction to be larger than ribbon residual
    model_input$ppr_gdf_bias[model_input$ppr_gdf_bias>model_input$ribbon_init & !is.na(model_input$ribbon_init)] = model_input$ribbon_init[model_input$ppr_gdf_bias>model_input$ribbon_init & !is.na(model_input$ribbon_init)] #ADDED 8/7/23
    
    model_input$ribbon_init_temp = model_input$ribbon_init-model_input$ppr_gdf_bias
    model_input$ribbon_init_temp[model_input$ribbon_init_temp<0] = 0
    
    model_input =  model_input %>% dplyr::group_by(azimuthal_discrete) %>%
      dplyr::mutate(ribbon_init0= c(rev(cummin(rev(ribbon_init_temp[dist_from_peak <= 0]))),cummin(ribbon_init_temp[dist_from_peak > 0])))
    model_input$ribbon_init_temp = model_input$ribbon_init0
    
    ### Choose whether to keep PPR modification
    if(sum(model_input$ribbon_init_temp, na.rm=T)>=0.5*sum(model_input$ribbon_init)){
      model_input$ribbon_init = model_input$ribbon_init_temp
    }
  }
  model_input$ribbon_init[model_input$mask == 0]=0 
  model_input$ribbon_init[model_input$ribbon_init<0 & !is.na(model_input$ribbon_init)]=0 
  model_input$ribbon_init[model_input$ribbon_init>model_input$ena_rate & !is.na(model_input$ena_rate)]=model_input$ena_rate[model_input$ribbon_init>model_input$ena_rate & !is.na(model_input$ena_rate)]=0
  
  ####################################################
  ### Estimate Additional Fine-Scale Uncertainties ###
  ####################################################
  
  ### Jackknife variance estimation, without implementing the PPR correction
  my_jack = function(x, data, THESEUS_DELINE){
    data$ribbon_init_temp = data$ena_rate - data$gdf_init
    data_sub = data[data$group != x,]
    ppr_resid_mod <- ppr(ribbon_init_temp ~ x+y+z,
                         nterms = 10,sm.method = 'gcvspline',
                         data = data_sub[data_sub$mask_temp == 0,])
    
    preds = predict(ppr_resid_mod, newdata=data, type="response")
    preds[preds>data$ribbon_init_temp & !is.na(data$ribbon_init_temp)] = data$ribbon_init_temp[preds>data$ribbon_init_temp & !is.na(data$ribbon_init_temp)] #ADDED 8/7/23
    preds[preds<0] = 0
    
    #don't allow negative gdf
    gdf_temp = data$gdf_init + preds
    gdf_temp[gdf_temp < 0 & !is.na(data$ena_rate)] = 0 #should never actually happen....
    gdf_temp[gdf_temp > data$ena_rate & !is.na(data$ena_rate)] = data$ena_rate[gdf_temp > data$ena_rate & !is.na(data$ena_rate)]
    preds_mod = gdf_temp - data$gdf_init
    return(preds_mod)
  }
  set.seed(1026)
  model_input$group = sample(x=c(1:10),size = length(unlist(model_input[,1])), replace = T)
  gc()
  JACK = matrix(NA, ncol = 10, nrow = length(unlist(model_input[,1])))
  for(x in 1:10){
    jack_temp = try(my_jack(x, model_input, THESEUS_DELINE))
    if(class(jack_temp)[1] != 'try-error'){
      JACK[,x] = jack_temp
    }
  }
  JACK = JACK[,apply(JACK,2,FUN = function(x){sum(!complete.cases(x))}) == 0]
  gc()
  #using Kott 2001 delete a group jackknife paper, pg 522-3.
  t_est = apply(JACK,1,FUN = function(x,w){weighted.mean(x,w)}, w = as.vector(table(model_input$group))/length(unlist(model_input$group)))
  R = length(JACK[1,])
  JACK_DIFF = sweep(JACK, 1, t_est, FUN = '-')
  JACKSE = sqrt(((R-1)/R)*apply(JACK_DIFF^2,1,sum))
  model_input$se.jack = JACKSE
  JACKCOV_L = sqrt((R-1)/R)*JACK_DIFF
  
  ### Don't allow additional uncertainty to be greater than uniform
  var_max = ((model_input$ena_rate - model_input$gdf_init)^2)/12
  model_input$se.jack = ifelse((model_input$se.jack^2)>var_max & !is.na(var_max), sqrt(var_max), model_input$se.jack)
  
  
  #################################################
  ### Reestimating Ribbon Peaks and Intensities ###
  #################################################
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  SUBSET = model_input
  if(K==6){
    SUBSET$eligible_angles = as.numeric(model_input$polar_discrete <= (model_input$phi_max_esa4+24) & model_input$polar_discrete >=  (model_input$phi_max_esa4-24) & model_input$polar_discrete >= 90  )
  }else{
    SUBSET$eligible_angles = as.numeric(model_input$polar_discrete <= (model_input$phi_max_esa4+10) & model_input$polar_discrete >=  (model_input$phi_max_esa4-10) & model_input$polar_discrete >= 90 )
  }
  SUBSET = SUBSET %>% dplyr::group_by(azimuthal_discrete) %>% 
    dplyr::mutate(phi_max = fitsplines_center_restricted(polar_discrete, ribbon_init, polar_discrete[eligible_angles==1])[1],
                  flux_max = fitsplines_center_restricted(polar_discrete, ribbon_init, polar_discrete[eligible_angles==1])[2],
                  obs_max = fitsplines_center_restricted(polar_discrete, ribbon_init, polar_discrete[eligible_angles==1])[3])
  SUBSET = SUBSET %>% dplyr::group_by(azimuthal_discrete) %>%
    dplyr::mutate(min_eligible = min(polar_discrete[eligible_angles==1], na.rm=T))
  SUBSET_SAVE = SUBSET
  SUBSET = SUBSET[!duplicated(SUBSET$azimuthal_discrete),]
  
  ### Clean and Smooth
  if(K == 6){
    SUBSET$phi_max[c(0,abs(diff(SUBSET$phi_max)))> 12 | is.na(SUBSET$phi_max)]=NA#don't allow a greater than 12 degree jump between neighboring pixels. 
    SUBSET$phi_max[abs(SUBSET$phi_max - SUBSET$phi_max_orig) > 12 | is.na(SUBSET$phi_max) | is.na(SUBSET$phi_max_orig)]=NA
  }else{
    SUBSET$phi_max[c(0,abs(diff(SUBSET$phi_max)))> 4 | is.na(SUBSET$phi_max)]=NA#don't allow a greater than 12 degree jump between neighboring pixels. 
    if(sum(abs(SUBSET$phi_max - SUBSET$phi_max_orig) > 4 | is.na(SUBSET$phi_max) | is.na(SUBSET$phi_max_orig),na.rm=T)<180){
      SUBSET$phi_max[abs(SUBSET$phi_max - SUBSET$phi_max_orig) > 4 | is.na(SUBSET$phi_max) | is.na(SUBSET$phi_max_orig)]=NA
    }else{
      SUBSET$phi_max[abs(SUBSET$phi_max - SUBSET$phi_max_orig) > 12 | is.na(SUBSET$phi_max) | is.na(SUBSET$phi_max_orig)]=NA
    }
  }
  SUBSET = data.frame(SUBSET[order(as.numeric(as.character(SUBSET$azimuthal_discrete))),])
  SUBSET$phi_max = as.numeric(forecast::tsclean(ts(SUBSET$phi_max)))
  PREDS_ANCHOR = Get_Loess_PHI_ANCHOR(SUBSET, span = 0.1)
  PHI_ANCHOR <<- as.vector(PREDS_ANCHOR)
  SUBSET$phi_max = PHI_ANCHOR
  
  ### Fix flux_max estimates
  if('phi_max_fixed' %in% names(SUBSET_SAVE)){SUBSET_SAVE = subset(SUBSET_SAVE, select = -c(phi_max_fixed))}
  SUBSET_SAVE = merge(SUBSET_SAVE, data.frame(azimuthal_discrete = SUBSET$azimuthal_discrete, 
                                              phi_max_fixed = SUBSET$phi_max), by = 'azimuthal_discrete')
  SUBSET_SAVE = SUBSET_SAVE %>% dplyr::group_by(azimuthal_discrete) %>% 
    dplyr::mutate(flux_max = evalsplines_center_restricted(polar_discrete, ribbon_init, phi_max_fixed))
  SUBSET_MAX = SUBSET_SAVE
  SUBSET_MAX = SUBSET_MAX[!is.na(SUBSET_MAX$ena_rate) & !is.nan(SUBSET_MAX$ena_rate),]
  SUBSET_MAX = SUBSET_MAX[!duplicated(SUBSET_MAX$azimuthal_discrete),]
  
  ### Merge
  if('phi_max' %in% names(model_input)){model_input = subset(model_input, select = -c(phi_max))}
  model_input = merge(model_input,data.frame(azimuthal_discrete = AZIMUTHAL, 
                                             phi_max=  PHI_ANCHOR), by = 'azimuthal_discrete',
                      all.x = T, all.y = F)
  if('flux_max' %in% names(model_input)){model_input = subset(model_input, select = -c(flux_max))}
  model_input = merge(model_input,data.frame(azimuthal_discrete = SUBSET_MAX$azimuthal_discrete, 
                                             flux_max=  SUBSET_MAX$flux_max), by = 'azimuthal_discrete',
                      all.x = T, all.y = F)
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  
  model_input$dist_from_peak = model_input$polar_discrete - model_input$phi_max
  model_input$above_peak = as.numeric(model_input$dist_from_peak>0)
  model_input$flux_max[model_input$flux_max<0.001 & !is.na(model_input$flux_max)]=0
  model_input$ribbon_init_scaled = model_input$ribbon_init/model_input$flux_max
  model_input$ribbon_init_scaled[is.na(model_input$ribbon_init_scaled)|is.infinite(model_input$ribbon_init_scaled)]=0
  model_input$ribbon_init_scaled[model_input$ribbon_init_scaled<0]=0
  
  ########################
  ### Fit Ribbon Model ###
  ########################
  
  model_input$dist1 = abs(model_input$dist_from_peak*as.numeric(model_input$above_peak==1)) # iSpline monotone increasing
  model_input$dist2 = abs(-model_input$dist_from_peak*as.numeric(model_input$above_peak==0))
  model_input$UPPER_THRESH = max(model_input$dist1,na.rm=T)
  model_input$LOWER_THRESH = max(model_input$dist2,na.rm=T)
  model_input$KNOTS = KNOTS
  model_input$AZI_KNOTS = AZI_KNOTS
  model_input$ribbon_init_scaled_inv = 1-model_input$ribbon_init_scaled
  if(K==6){
    ### The following line is needed for ISOC. Helps a lot with stability
    model_input$ribbon_init_scaled_inv[model_input$ribbon_init_scaled_inv < (-0.25) & (model_input$dist1 > 4 | model_input$dist2 > 4)] = NA
  }
  fit = mgcv::gam(ribbon_init_scaled_inv~splines2::iSpline(dist1, knots = seq(0.01,unique(UPPER_THRESH)-1, length.out = unique(KNOTS)),
                                                           Boundary.knots = c(0,unique(UPPER_THRESH)), intercept = FALSE)*pbs(x=azimuthal_discrete, df = unique(AZI_KNOTS), periodic = TRUE, Boundary.knots = c(0,360),intercept = TRUE)+
                    splines2::iSpline(dist2, knots = seq(0.01,unique(LOWER_THRESH)-1, length.out = unique(KNOTS)),
                                      Boundary.knots = c(0,unique(LOWER_THRESH)), intercept = FALSE)*pbs(x=azimuthal_discrete, df = unique(AZI_KNOTS), periodic = TRUE, Boundary.knots = c(0,360),intercept = TRUE),
                  data = model_input, select = T)
  model_input$ribbon_profile = predict(fit,newdata = model_input[,c('azimuthal_discrete','dist1','dist2', 'UPPER_THRESH','LOWER_THRESH','KNOTS', 'AZI_KNOTS')],  se.fit = F, type = 'response' )
  model_input$ribbon_profile = 1-model_input$ribbon_profile
  model_input$ribbon_profile[model_input$ribbon_profile<0]=0
  if(K==6){
    ### The following line is needed for ISOC. Helps a lot with stability
    model_input$ribbon_profile[ model_input$ribbon_profile>1.25] = 1 
  }
  model_input$ribbon_profile[abs(model_input$dist_from_peak)<1 & model_input$ribbon_profile<model_input$max_val]=model_input$max_val[abs(model_input$dist_from_peak)<1 & model_input$ribbon_profile<model_input$max_val]
  model_input$ribbon_profile[model_input$mask==0]=0
  model_input$ribbon_fitted = model_input$ribbon_profile*model_input$flux_max
  model_input$ribbon_fitted[model_input$mask == 0]=0 #
  model_input$ribbon_fitted[model_input$ribbon_fitted<0 & !is.na(model_input$ribbon_fitted)]=0 
  
  
  #######################
  ### Refit GDF Model ###
  #######################
  model_input = model_input %>% dplyr::group_by(azimuthal_discrete) %>% dplyr::mutate(mask_lower = min(polar_discrete[mask == 1],na.rm=T),
                                                                                      mask_upper = max(polar_discrete[mask == 1],na.rm=T))
  model_input$gdf_est2 = model_input$ena_rate - model_input$ribbon_fitted
  model_input$gdf_est2[model_input$gdf_est2<0]=0
  model_input$gdf_est2[model_input$gdf_est2>model_input$ena_rate & !is.na(model_input$ena_rate)] = model_input$ena_rate[model_input$gdf_est2>model_input$ena_rate & !is.na(model_input$ena_rate)]
  model_input$gdf_est2[model_input$mask == 0] = model_input$ena_rate[model_input$mask == 0]
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  
  ### Define edges to handle edge-of-mask effects
  if(K == 6){
    model_input = model_input %>% dplyr::group_by(azimuthal_discrete) %>% 
      dplyr::mutate(edge_inner = as.numeric(polar_discrete %in%  c(unique(mask_upper)-0:(1*K)) | 
                                              polar_discrete %in%  c(unique(mask_lower)+0:(1*K))))    
    model_input = model_input %>% dplyr::group_by(azimuthal_discrete) %>% 
      dplyr::mutate(edge_outer = as.numeric(polar_discrete %in%  c(unique(mask_upper)+0:(1*K)) | 
                                              polar_discrete %in%  c(unique(mask_lower)-0:(1*K))))    
  }else{
    model_input = model_input %>% dplyr::group_by(azimuthal_discrete) %>% 
      dplyr::mutate(edge_inner = as.numeric(polar_discrete %in%  c(unique(mask_upper)-0:(3*K)) | 
                                              polar_discrete %in%  c(unique(mask_lower)+0:(3*K))))     
    model_input = model_input %>% dplyr::group_by(azimuthal_discrete) %>% 
      dplyr::mutate(edge_outer = as.numeric(polar_discrete %in%  c(unique(mask_upper)+0:(3*K)) | 
                                              polar_discrete %in%  c(unique(mask_lower)-0:(3*K))))  
    model_input = model_input %>% dplyr::group_by(azimuthal_discrete) %>% 
      dplyr::mutate(edge_actual = as.numeric(polar_discrete %in%  c(unique(mask_upper)+c(-K:K)) | 
                                               polar_discrete %in%  c(unique(mask_lower)+c(-K:K)))) 
  }
  
  model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
  
  ### Re-estimate GDF with GAM
  mybasis <- "sos"
  gam_mod2 <- mgcv::gam(gdf_est2 ~ s(lat, lon,bs=mybasis) ,
                        data=model_input, family = gaussian())
  model_input$gdf_init2 <- predict(gam_mod2, newdata = model_input, type="response")
  gc()
  
  ### Refine GDF with PPR
  model_input$resid_gdf = model_input$gdf_est2 - model_input$gdf_init2
  ppr_resid_mod <- ppr(resid_gdf ~ x+y+z,
                       nterms = 100, sm.method = 'gcvspline',
                       data = model_input[model_input$edge_inner == 0,])
  model_input$ppr_gdf_bias <- predict(ppr_resid_mod, newdata=model_input, type="response")
  
  model_input$gdf_est3 = model_input$gdf_init2+model_input$ppr_gdf_bias
  model_input$gdf_est3[model_input$gdf_est3>model_input$ena_rate & !is.na(model_input$ena_rate)] = model_input$ena_rate[model_input$gdf_est3>model_input$ena_rate & !is.na(model_input$ena_rate)]
  model_input$gdf_est3[is.na(model_input$ena_rate)]=NA
  model_input$gdf_est3[model_input$mask == 0]= model_input$ena_rate[model_input$mask == 0]
  model_input$gdf_est3[model_input$gdf_est3<0 & !is.na(model_input$gdf_est3)]=0
  model_input$gdf_final = model_input$gdf_est3
  model_input$ribbon_final = model_input$ena_rate - model_input$gdf_final
  
  #########################################
  ### Final Edge and Monotonicity Fixes ###
  #########################################
  if(K == 2){
    
    Get_InterpEdge_GDF = function(lat, val, edge_exclude){
      temp = data.frame(lat = lat, val = val, edge_exclude = edge_exclude)
      spline_func = stats::splinefun(temp$val[edge_exclude == 0]~temp$lat[edge_exclude == 0], method = 'natural')
      preds = spline_func(temp$lat)
      return(preds)
    }
    model_input = model_input[order(model_input$azimuthal_discrete, model_input$polar_discrete),]
    
    model_input = model_input %>% dplyr::group_by(azimuthal_discrete) %>%
      dplyr::mutate(gdf_final_smoothed = Get_InterpEdge_GDF(polar_discrete, gdf_final, as.numeric(edge_actual == 1)))
    
    model_input$gdf_final_smoothed[is.na(model_input$gdf_final_smoothed)]=NA
    model_input$gdf_final_smoothed[model_input$gdf_final_smoothed<0 & !is.na(model_input$gdf_final_smoothed)]=0
    model_input$gdf_final_smoothed[model_input$gdf_final_smoothed>model_input$ena_rate & !is.na(model_input$ena_rate)] = model_input$ena_rate[model_input$gdf_final_smoothed>model_input$ena_rate & !is.na(model_input$ena_rate)]
    model_input$gdf_final = model_input$gdf_final_smoothed
    model_input$ribbon_final = model_input$ena_rate - model_input$gdf_final
  }
  
  #############################
  ### Write Results to File ###
  #############################
  RESULTS = model_input[,c('polar_discrete', 'azimuthal_discrete','lon', 'lat', 
                           'ena_rate', 'sd_ena_rate',
                           'gdf_final', 
                           'ribbon_profile', 'ribbon_fitted','ribbon_final', 
                           'mask', 'se.fit', 'se.jack', 'mask_temp')]
  RESULTS = RESULTS[order(RESULTS$azimuthal_discrete, RESULTS$polar_discrete),]
  gc()
  return(list(RESULTS = RESULTS, L_G_GAM = vcovpred_lower, L_G_JACK = JACKCOV_L))
}

