########################
########################
### Helper Functions ###
########################
########################


### Performs Sobel Decomposition
get_sobel = function(azimuthal_discrete, polar_discrete, ena_rate){
  MAP = data.frame(azimuthal_discrete = azimuthal_discrete, 
                   polar_discrete = polar_discrete, 
                   ena_rate = ena_rate)
  DATA_DRAW=reshape(data = data.frame(MAP[,c('azimuthal_discrete', 'polar_discrete', 'ena_rate')]), idvar = 'polar_discrete', timevar = 'azimuthal_discrete', v.names = c('ena_rate'), direction = 'wide')
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
  im_laplace_sub[,1:2]=min(im_laplace_sub)
  im_laplace_sub[,89:90]=min(im_laplace_sub)
  return(as.vector(t(im_laplace_sub)))
}

### Performs Sobel Decomposition for ISOC
get_sobel_isoc = function(azimuthal_discrete, polar_discrete, ena_rate){
  MAP = data.frame(azimuthal_discrete = azimuthal_discrete, 
                   polar_discrete = polar_discrete, 
                   ena_rate = ena_rate)
  DATA_DRAW=reshape(data = data.frame(MAP[,c('azimuthal_discrete', 'polar_discrete', 'ena_rate')]), idvar = 'polar_discrete', timevar = 'azimuthal_discrete', v.names = c('ena_rate'), direction = 'wide')
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
  im_laplace_sub[,1:2]=min(im_laplace_sub,na.rm=T)
  im_laplace_sub[,29:30]=min(im_laplace_sub,na.rm=T)
  return(as.vector(t(im_laplace_sub)))
}
