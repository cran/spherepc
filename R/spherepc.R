

############################################################################
### Spherepc R package: dimension reduction methods in a sphere
############################################################################

### install rgl, sphereplot, and geosphere packages--------------------------

#install.packages("rgl")
#install.packages("sphereplot")
#install.packages("geosphere")
library(rgl)
library(sphereplot)
library(geosphere)

##########################################################################
# Define functions
##########################################################################

  Trans.Euclid <- function(vec){                # Input: 'vec' is 2-dim spherical coordinate (longitude, latitude). Output: 3-dim Euclidean coordinate.
    vec <- as.numeric(vec)
    if (length(vec) != 2) stop("length of vector is not 2")
    return( c(cospi(vec[1]/180) * cospi(vec[2]/180), cospi(vec[2]/180) * sinpi(vec[1]/180), sinpi(vec[2]/180)) )
  }

  
  Trans.sph <- function(vec){                  # Input: 'vec' is 3-dim Euclidean coordinate. Output: 2-dim spherical coordinate (longitude, latitude).
    vec <- as.numeric(vec)
    if (sum(vec^2)  < 1e-15) stop("Input should not be (0, 0, 0)")
    if (length(vec) != 3) stop("length of vector is not 3")
    vec <- vec/(sum(vec^2))^(1/2)               # Normalize so that 'vec' is in the unit sphere.
    if ((vec[2] >= 0) && (vec[1] >= 0)){
      if (vec[1]^2 + vec[2]^2 == 0){
        if (vec[3] > 0){
          return(c(0, 90))
        } else {
          return(c(0, -90))
        }
       } else {
         return (c(atan(vec[2]/ vec[1]), asin(vec[3])) * 180/pi)
       }
    }
    if ((vec[2] >= 0) && (vec[1] < 0)){
      return (c(pi + atan(vec[2]/vec[1]), asin(vec[3])) * 180/pi)
    }
    if ((vec[2] < 0) && (vec[1] >= 0)){
      return (c(atan(vec[2]/vec[1]), asin(vec[3])) * 180/pi)
    }
    if ((vec[2] < 0) && (vec[1] < 0)){
      return (c(atan(vec[2]/vec[1]) - pi, asin(vec[3])) * 180/pi)
    }
  }

  Logmap <- function(vec){                                                                                 # Input: 3-dim Euclidean coordinate. Output: 2-dim Euclidean coordinate.
    if (sum(vec^2) < 1e-15) stop("Input is not the zero.")
    vec <- as.numeric(vec)
    vec <- vec/(sum(vec^2))^(1/2)                                                                           # Normalize so that 'vec' is in the unit sphere.
    if (length(vec) != 3) stop("length of vector is not 3")
    if (sum(vec == c(0, 0, -1)) == 3) stop("input is the cut locus of (0, 0, 1)")                           # Input should not be (0, 0, -1).
    z <- vec[3]
    x <- vec[1]
    y <- vec[2]
    result <- c(x * asin((1 - z^2)^(1/2))/(1 - z^2)^(1/2), y * asin((1 - z^2)^(1/2))/(1 - z^2)^(1/2))
    if (is.nan(result[1]) == TRUE){
        return (c(0, 0))
      }
        return (result)
  }

  Expmap <- function(vec){                                                                                 # Input: 2 dim-Euclidean vector. Output: 3 dim-Euclidean vector.
    vec <- as.numeric(vec)
    if (length(vec) != 2) stop("number of vector component is not 2")
    if(sum(vec^2) < 1e-15){
      return(c(0, 0, 1))
    } else {
      norm <- (sum(vec^2))^(1/2)
      return( c(vec[1] * sin(norm)/norm, vec[2] * sin(norm)/norm, cos(norm)) )
    }
  }

  Crossprod <- function(vec1, vec2){                                                                      # Cross product of 3-dim two vectors.
    if (length(vec1) != 3 | length(vec2) != 3) stop('number of vector component is not 3')
    x <- vec1[2] * vec2[3] - vec1[3] * vec2[2]
    y <- vec1[3] * vec2[1] - vec1[1] * vec2[3]
    z <- vec1[1] * vec2[2] - vec1[2] * vec2[1]
    result <- c(x, y, z)
    return (result)
  }

  Rotate <- function(pt1, pt2){                                                                           # Input 'pt1', 'pt2': spatial locations (longitude, latitude). Output: 3 dim-Euclidean coordinate.
    pt1 <- as.numeric(pt1)                                                                                 # Rotation function: Rotate 'pt2' by the extent from which rotate 'pt1' to spherical coordinate (0, 90) (= (0, 0, 1) Euclidean coordinate.)
    pt2 <- as.numeric(pt2)
    a <- Trans.Euclid(pt1)
    v <- Trans.Euclid(pt2)
    axis <- Crossprod(a, c(0, 0, 1))                                                                       # Rotation axis.
    if (sum(axis^2) < 1e-15){
      return(Trans.Euclid(pt2))
    }
    axis <- axis/(sum(axis^2))^(1/2)                                                                       # Normalized unit axis.
    theta <- abs(pi/2 - pt1[2] * pi/180)                                                                   # Rotation angle.
    v.rot <- v * cos(theta) + Crossprod(axis, v) * sin(theta) + axis * (sum(axis * v)) * (1 - cos(theta))  # Rodrigues' rotation formula.
    return (v.rot)
  }
 
  Rotate.inv <- function(pt1, pt2){                                                                       # Input 'pt1', 'pt2': measurements of angular form of (longitude, latitude). Output: 3-dim Euclidean vector.
    pt1 <- as.numeric(pt1)                                                                                 # Inverse Rotation function: inverse rotation of the 'Rotate' function.
    pt2 <- as.numeric(pt2)
    a <- Trans.Euclid(pt1)
    v <- Trans.Euclid(pt2)
    axis <- Crossprod(c(0, 0, 1), a)
    if (sum(axis^2) < 1e-15){
      return(Trans.Euclid(pt2))
    }
    axis <- axis/(sum(axis^2))^(1/2)                                                                       # normalized unit axis
    theta <- abs(pi/2 - pt1[2] * pi/180)                                                                   # rotation angle.
    v.rot <- v * cos(theta) + Crossprod(axis, v) * sin(theta) + axis * (sum(axis * v)) * (1 - cos(theta))  # Rodrigues' rotation formula
    return (v.rot)
  }

  Kernel.quartic <- function(vec){                                     # kernel function; if 'vec' is a vector, it returns a vector
    result <- c()
    for (i in 1:length(vec)){
      if (abs(vec[i]) <= 1){
        result[i] <- (1 - vec[i]^2)^2
      } else {
        result[i] <- 0
      }
    }
      return(result)
  }

  Kernel.indicator <- function(vec){                                   # Kernel function. If 'vec' is vector, it returns vector.
    result <- c()
    for (i in 1:length(vec)){
     if (abs(vec[i]) <= 1){
       result[i] <- 1
     } else {
       result[i] <- 0
     }
    }
    return(result)
  }

  Kernel.Gaussian <- function(vec){                                    # Kernel function. If 'vec' is a vector, it returns a vector.
    result <- c()
    for (i in 1:length(vec)){
      if (abs(vec[i]) <= 1){
        result[i] <- exp(-vec[i]^2)
      } else {
        result[i] <- 0
      }
    }
    return(result)
  }

  Cal.recon <- function(data, line){                                  # calculating reconstruction error from data to line.
    r <- 6378137                                                      # earth radius (m)
    rownames(data) <- NULL                                            # 'data' and 'line' should be n * 2 matrix or data frame.
    if (nrow(line) == 0){
       return (0)
    } else if (nrow(line) == 1){
       return (sum((geosphere::distGeo(data, line)/r)^2))
    } else {
       proj <- geosphere::dist2Line(data, line, distfun = geosphere::distGeo)
       return (sum((proj[, 1]/r)^2))
    }
  }

  ExtrinsicMean <- function(data, weights = rep(1, nrow(data))){                        # 'data': longitude/latitude matrix or dataframe with two column
    if (sum(weights^2)^(1/2) < 1e-15) stop("weight should not be the zero vector.")     # 'weights': a weight vector whose length is the rownumber of points
    weights <- weights/sum(weights)                                                     #  normalize 'weights' so that its components add up to 1
    mean.Euc <- c(0, 0, 0)
    for (i in 1:nrow(data)){
      mean.Euc <- mean.Euc + Trans.Euclid(data[i, ]) * weights[i]
    }
    if (sum(mean.Euc^2)^(1/2) < 1e-15){
      stop("extrinsic mean is not well-defined. Data should be more localized.")
    }
    mean.Euc <- mean.Euc/sum(mean.Euc^2)^(1/2)
    mean <- Trans.sph(mean.Euc)
    return(as.numeric(mean))
  }

  IntrinsicMean <- function(data, weights = rep(1, nrow(data)), thres = 1e-5){         # 'data': longitude/latitude matrix or dataframe with two column.
    if (sum(weights^2)^(1/2) < 1e-15) stop("weight should not be the zero vector.")     # 'weights': n-dimensional vector.
    mu <- data[1, , drop = F]                                                          # Initialize mean as a point.
    delta.mu <- c(1, 0)
    while (sum((delta.mu)^2)^(1/2) > thres){
      weights.normal <- weights/sum(weights)                                           # Normalize of 'weights' so that its components add up to 1.
      summation <- c(0, 0)
        for (m in 1:length(weights)){
          rot <- Rotate(mu, data[m, ])
          summation <- summation + weights.normal[m] * Logmap(rot)
        }
      delta.mu <- summation
      exp <- Expmap(delta.mu)
      mu.Euc <- Rotate.inv(mu, Trans.sph(exp))
      mu <- Trans.sph(mu.Euc)
    }
    return(as.numeric(mu))
  }

### PrincipalCircle
  PrincipalCircle <- function(data, step.size = 1e-3, thres = 1e-5, maxit = 10000){
    # 'data': a matrix or data.frame of data points represented by angular form of (longitude/latitude)
    # To convergence of this algorithm, step.size is recommended below 0.01.
    circle <- IntrinsicMean(data) * pi/180                                          # Initialized by center of data to avoid being trapped at local minima. (radian measurement)
    circle[3] <- pi/2
    data <- data * pi/180                                                           # Transformation from angle into radian.
    delta <- 1
    iter <- 0
    while ((delta > thres) && (iter < maxit)){
      iter <- iter + 1
      grad <- c(0, 0, 0)
      for(i in 1:nrow(data)){
        temp <- sin(circle[2]) * sin(data[i, 2]) + cos(circle[2]) * cos(data[i, 2]) * cos(circle[1] - data[i, 1])
        names(temp) <- NULL
        if (temp == 1){
          next
        }
        grad[1] <- grad[1] + (acos(temp) - circle[3]) * 1/sqrt(1-temp^2) * cos(data[i, 2]) * cos(circle[2]) * sin(circle[1] - data[i, 1])
        grad[2] <- grad[2] - (acos(temp) - circle[3]) * 1/sqrt(1-temp^2) * (cos(circle[2]) * sin(data[i, 2]) - sin(circle[2]) * cos(data[i, 2]) * cos(circle[1] - data[i, 1]))
        grad[3] <- grad[3] + circle[3] - acos(temp)
      }
      grad[1:2] <- grad[1:2] * pi/180
      circle <- circle - step.size * 2 * grad
      if (circle[3] > pi) stop("step.size should be lowered.")
      delta <- sqrt(sum((step.size * 2 * grad)^2))
    }
    circle[1:2] <- circle[1:2] * 180/pi
    return(circle)                                                                 # Output: angular form of (longitude, latitude)
  }

  GenerateCircle <- function(center, radius, T = 1000){
    # center should be an angular form of (longitude, latitude).
    # the radius r is in [0, pi]
    # the number of points in circle is T
    lon <- seq(0, 360, length.out = T + 1)
    generate.circle <- matrix(nrow = T, ncol = 2)
    for (i in 1:T){
      generate.circle[i, ] <- Trans.sph(Rotate.inv(center, c(lon[i], (pi/2 - radius) * 180/pi)))
    }
    return(generate.circle)
  }

  PGA <- function(data, col1 = "blue", col2 = "red"){
    # 'data' should be a matrix or data frame with 2 column (longitude/latitude)
    center <- IntrinsicMean(data)
    ant <- Trans.sph(-Trans.Euclid(center))
    cov <- matrix(c(0, 0, 0, 0), nrow = 2)
    for (i in 1:nrow(data)){
       z <- Rotate(center, data[i, ])
       y <- Logmap(z)
       cov <- cov + y %*% t(y)
    }
    eivec <- eigen(cov)$vectors[, 1]
    direc <- eivec/(sum(eivec^2))^(1/2)
    pt <- Expmap(direc)
    pt.sph <- Trans.sph(pt)
    Euc <- Rotate.inv(center, pt.sph)
    mid <- Trans.sph(Euc)
    mid.ant <- Trans.sph(-Euc)
    line <- geosphere::makeLine(rbind(center, mid, ant, mid.ant, center))
    # plot
    sphereplot::rgl.sphgrid(col.long = 'black', col.lat = 'black', radaxis = TRUE)
    sphereplot::rgl.sphpoints(data, radius = 1, col = col1, size = 12)
    sphereplot::rgl.sphpoints(line, radius = 1, col = col2, size = 9)
    fit <- list(line = line)
    return(fit)
  }


  Proj.Hauberg <- function(data, line){
     # this function performs the approximated projections
     # each data point is projected onto the nearest points of line (not exactly)
     proj <- matrix(nrow = nrow(data), ncol = 2)
     for (i in 1:nrow(data)){
       proj[i, ] <- line[which.min(geosphere::distGeo(data[i, ], line)), ]
     }
     return(proj)
  }

  Dist.pt <- function(data){
     # this function calculates the number of distinct points
     # 'data': matrix or data.frame of data points represented by angular form of (longitude, latitude).
     num <- nrow(data)
     if (num == 1){
       return(1)
     }
     dist.num <- 1
     for (i in 2:num){
       temp <- c()
       temp <- colSums(t(data[1:(i - 1), , drop = F]) == data[i, ]) < 2
       prod <- 1
       for (j in 1:(i - 1)){
         prod <- prod * temp[j]
       }
       dist.num <- dist.num + prod
     }
   return(dist.num)
  }

  ########################################
  ### Spherical Principal Curves (global)
  #######################################
  SPC <- function(data, q = 0.1, T = nrow(data), step.size = 1e-3, maxit = 10, type = "Intrinsic", thres = 1e-2,
                 deletePoints = FALSE, plot.proj = FALSE, kernel = "quartic", 
                 col1 = "blue", col2 = "green", col3 = "red") {
    # 'data': matrix or data.frame of data points represented by angular form of (longitude, latitude).
    # 'q': Smoothing parameter in Expectation step.
    # 'T': The number of points making up circle.
    # 'step.size': It is recommended below 0.01 owing to convergence of this algorithm.
    # 'maxit': The maximum number of iterations
    # 'deletePoints': TRUE or FALSE. If it is TRUE, then for each expectation step delete the points which do not have adjacent data.
    # If 'deletePoints' is FALSE, leave the points which do not have adjacent data for each expectation step.
    # 'col1' is color of 'data' and 'col2' is that of points making up the principal curves; in addition, 'col3' is the color of the principal curves.

    r <- 6378137                                                # earth radius (m)
    PC <- PrincipalCircle(data, step.size = step.size)
    prin.curves <- GenerateCircle(PC[1:2], PC[3], T = T)        # initialize principal curves as principal circle
    if (nrow(prin.curves) == 1){
      rss <- sum((geosphere::distGeo(data, prin.curves)/r)^2)   # residual sum of squares (rss)
    } else {
      proj <- geosphere::dist2Line(data, prin.curves, distfun = geosphere::distGeo)
      rownames(data) <- NULL
      rss <- Cal.recon(data, prin.curves)
    }
    rss.new <- rss.old <- rss
    relative.change <- 1
    iter <- 0                                                   # iteration of principal curve algorithm
    
    while ((relative.change >= thres) && (rss.old != 0) && (iter < maxit + 1)){
      ## projection step
      T.temp <- nrow(prin.curves)
      if (T.temp > 1){
        proj.temp <- geosphere::dist2Line(data, prin.curves, distfun = geosphere::distGeo)
      } else {
        proj.temp <- cbind(cbind(geosphere::distGeo(data, prin.curves), rep(prin.curves[, 1], nrow(data))), rep(prin.curves[, 2], nrow(data)))
      }
      
      curve.length <- 0
      for (t in 1:(T.temp - 1)){
        curve.length <- curve.length + geosphere::distGeo(prin.curves[t, ], prin.curves[t + 1, ])/r
      }
      if (deletePoints == FALSE){
        curve.length <- curve.length + geosphere::distGeo(prin.curves[T.temp, ], prin.curves[1, ])/r         # length of principal curves
      }
      
      ## Expectation step
      for (t in 1:T.temp){
        choice <- weights <- c()
        choice <- geosphere::distGeo(proj.temp[, 2:3], prin.curves[t, ])/r < q * curve.length
        adjacent <- as.data.frame(data[choice, , drop = F])                                                  # adjacent points
        if (is.na(sum(choice)) == TRUE) {warning("Errors, 'choice' have not a number or non available.")}
        if (sum(choice) == 0){                                                                               # If a point dosen't have adjacent data, there are two options: deleting and leaving the point (deletePoints)
          if (deletePoints == FALSE){
            prin.curves[t, ] <- prin.curves[t, ]                                                             # delete the points which do not have adjacent data each expectation step
          } else {        # deletePoints == TRUE
            prin.curves[t, ] <- c(NA, NA)                                                                    # leave the points which do not have adjacent data for each expectation step
          }
          
        } else {
          
          if ((kernel == "quartic") | (kernel == "Quartic")){
            weights <- Kernel.quartic(geosphere::distGeo(proj.temp[choice, 2:3], prin.curves[t, ])/(r * q * curve.length))
          } else if ((kernel == "Gaussian") | (kernel == "gaussian")){
            weights <- Kernel.Gaussian(geosphere::distGeo(proj.temp[choice, 2:3], prin.curves[t, ])/(r * q * curve.length))
          } else if ((kernel == "indicator") | (kernel == "Indicator")){
            weights <- Kernel.indicator(geosphere::distGeo(proj.temp[choice, 2:3], prin.curves[t, ])/(r * q * curve.length))
          } else {
            stop("Errors, kernel should be quartic, Gaussian, and indicator.")
          }
          
          if ((type == "Intrinsic") | (type == "intrinsic")){
            mean <- IntrinsicMean(adjacent, weights)                                                         # calculating intrinsic mean of adjacent pts
            prin.curves[t, ] <- mean
          } else if ((type == "Extrinsic") | (type == "extrinsic")){
            mean <- ExtrinsicMean(adjacent, weights)                                                         # calculating extrinsic mean of adjacent pts
            prin.curves[t, ] <- mean
          } else {
            stop("Errors, type should be Intrinsic or Extrinsic")
          }
          
        }
      }
      
      prin.curves <- prin.curves[!is.na(prin.curves[, 1]), , drop = F]
      rss.new <- Cal.recon(data, prin.curves)
      relative.change <- abs((rss.old - rss.new)/rss.old)
      rss.old <- rss.new
      iter <- iter + 1
    }
   
    proj <- geosphere::dist2Line(data, prin.curves)
    distnum <- Dist.pt(proj[, 2:3, drop = F])                                              # The number of distinct projections
    if (nrow(prin.curves) > 1){
      if (deletePoints == FALSE){ 
        line <- geosphere::makePoly(prin.curves)
      } else {     # deletePoints == TRUE
        line <- geosphere::makeLine(prin.curves)     
      } 
    } else {
      line <- prin.curves
    }
    
    if (iter < maxit){
      converged <- TRUE
    } else {
      converged <- FALSE 
    }
    
    sphereplot::rgl.sphgrid(col.long = "black", col.lat = "black", radaxis = FALSE)                                                             # plot the data and principal curves.
    sphereplot::rgl.sphpoints(data[, 1], data[, 2], radius = 1, col = col1, size = 12)
    sphereplot::rgl.sphpoints(prin.curves[, 1], prin.curves[, 2], radius = 1, col = col2, size = 4)
    sphereplot::rgl.sphpoints(line[, 1], line[, 2], radius = 1, col = col3, size = 6)
    
    if (plot.proj == TRUE){
      line.proj <- matrix(ncol = 2)
      for (i in 1:nrow(proj)){                                                            # plot the projection line
        if (abs(sum((data[i, ] - proj[i, 2:3])^2)) < 1e-10){
          line.proj <- data[i, , drop = F]
        } else {
          line.proj <- geosphere::makeLine(rbind(data[i, ], data[i, ], proj[i, 2:3]))
        }
        sphereplot::rgl.sphpoints(line.proj[, 1], line.proj[, 2], radius = 1, col = "black", size = 1)
      }
    }
    
    fit <- list(prin.curves = prin.curves, line = line, converged = converged, iteration = iter,
                recon.error = rss.new, num.dist.pt = distnum)
    return(fit)
 }

 
 SPC.Hauberg <- function(data, q = 0.1, T = nrow(data), step.size = 1e-3, maxit = 10, type = "Intrinsic",
                         thres = 1e-2, deletePoints = FALSE, plot.proj = FALSE, kernel = "quartic",
                         col1 = "blue", col2 = "green", col3 = "red") {
   # 'data': matrix or data.frame of data points represented by angular form of (longitude, latitude)
   # 'q': smoothing parameter in expectation step
   # 'T': the number of points in initial circle
   # 'step.size': step size of the 'PrincipalCircle' function. It is recommended to be below 0.01 owing to convergence of the algorithm of 'PrincipalCircle' function
   # 'maxit': the number of iterations
   # 'deletePoints': TRUE or FALSE. If it is TRUE, then for each expecatation step delete the points which do not have adjacent data
   #  if 'deletePoints' is FALSE, leave the points which do not have adjacent data for each expectation step
   # 'col1' is color of data and 'col2' is that of points making up the resulting principal curves; in addition, 'col3' is that of principal curves.

   r <- 6378137                                                        # earth radius(m)
   PC <- PrincipalCircle(data, step.size = step.size)
   prin.curves <- GenerateCircle(PC[1:2], PC[3], T = T)                # principal circle is set to be an initialization of principal curves.
   if (nrow(prin.curves) > 1){
     proj <- Proj.Hauberg(data, prin.curves)
     rss <- sum((geosphere::distGeo(data, proj)/r)^2)                  # residual sum of squares (approximated rss)
   } else {
     rss <- sum((geosphere::distGeo(data, prin.curves)/r)^2)
   }
   rss.new <- rss.old <- rss
   relative.change <- 1
   iter <- 0                                                           # iteration of principal curve algorithm.
   
   while ((relative.change >= thres) && (rss.old != 0) && (iter < maxit + 1)){
     ## projection step
     T.temp <- nrow(prin.curves)
     if (T.temp > 1){
       proj.temp <- Proj.Hauberg(data, prin.curves)
     } else {
       proj.temp <- cbind(geosphere::distGeo(data, prin.curves), rep(prin.curves[, 1], nrow(data)), rep(prin.curves[, 2], nrow(data)))
     }
     
     curve.length <- 0
     for (t in 1:(T.temp - 1)){
       curve.length <- curve.length + geosphere::distGeo(prin.curves[t, ], prin.curves[t + 1, ])/r
     }
     if (deletePoints == FALSE){
       curve.length <- curve.length + geosphere::distGeo(prin.curves[T.temp, ], prin.curves[1, ])/r         # length of principal curves
     }
     
     ## Expectation step
     for (t in 1:T.temp){
       choice <- weights <- c()
       choice <- geosphere::distGeo(proj.temp, prin.curves[t, ])/r < q * curve.length
       adjacent <- as.data.frame(data[choice, , drop = F])                                                  # adjacent points
       if (is.na(sum(choice)) == TRUE) {warning("Errors, object 'choice' have not a number or non available.")}
       if (sum(choice) == 0){                                                                               # If a point dosen't have adjacent data, there are two options: deleting and leaving the point. (deletePoints)
         if (deletePoints == FALSE){
           prin.curves[t, ] <- prin.curves[t, ]                                                             # leave the points which do not have adjacent data for each expectation step
         } else {
           prin.curves[t, ] <- c(NA, NA)                                                                    # delete the points which do not have adjacent data for each expectation step
         }
         
       } else {
         if ((kernel == "quartic") | (kernel == "Quartic")){
           weights <- Kernel.quartic(geosphere::distGeo(proj.temp[choice, , drop = F], prin.curves[t, ]) / (r * q * curve.length))
         } else if ((kernel == "Gaussian") | (kernel == "gaussian")){
           weights <- Kernel.Gaussian(geosphere::distGeo(proj.temp[choice, , drop = F], prin.curves[t, ]) / (r * q * curve.length))
         } else if ((kernel == "indicator") | (kernel == "Indicator")){
           weights <- Kernel.indicator(geosphere::distGeo(proj.temp[choice, , drop = F], prin.curves[t, ]) / (r * q * curve.length))
         } else {
           stop("Errors, kernel should be quartic, Gaussian, and indicator.")
         }
         
         if ((type == "Intrinsic") | (type == "intrinsic")){
           mean <- IntrinsicMean(adjacent, weights)                                                         # calculating intrinsic mean of adjacent pts
           prin.curves[t, ] <- mean
         } else if ((type == "Extrinsic") | (type == "extrinsic")){
           mean <- ExtrinsicMean(adjacent, weights)                                                         # calculating extrinsic mean of adjacent pts
           prin.curves[t, ] <- mean
         } else {
           stop("Errors, type should be Intrinsic or Extrinsic")
         }
         
       }
     }
     
     prin.curves <- prin.curves[!is.na(prin.curves[, 1]), , drop = F]                        # deleting NA rows
     proj.temp <- Proj.Hauberg(data, prin.curves)
     rss.new <- sum((geosphere::distGeo(data, proj.temp)/r)^2)
     relative.change <- abs((rss.old - rss.new)/rss.old)                                     # relative change (stop when it is below 'thres')
     rss.old <- rss.new
     iter <- iter + 1
   }
   
   proj <- Proj.Hauberg(data, prin.curves)
   distnum <- Dist.pt(proj)                                                                  # the number of distinct projections
   if (nrow(prin.curves) > 1){
     if (deletePoints == FALSE){ 
       line <- geosphere::makePoly(prin.curves)     
     } else {            # deletePoints == TRUE
       line <- geosphere::makeLine(prin.curves)   
     }
   } else {
     line <- prin.curves
   }
   
   if (iter < maxit){
     converged <- TRUE
   } else {
     converged <- FALSE
   }
   
   sphereplot::rgl.sphgrid(col.long = "black", col.lat = "black")                                                                                # Plot the data and principal curves.
   sphereplot::rgl.sphpoints(data[, 1], data[, 2], radius = 1, col = col1, size = 12)
   sphereplot::rgl.sphpoints(prin.curves[, 1], prin.curves[, 2], radius = 1, col = col2, size = 4)
   sphereplot::rgl.sphpoints(line[, 1], line[, 2], radius = 1, col = col3, size = 6)
   
   if (plot.proj == TRUE){
     line.proj <- matrix(ncol = 2)
     for (i in 1:nrow(proj)){                                                                               # Plot the projection line
       if (abs(sum((data[i, ] - proj[i, ])^2)) < 1e-15){
         line.proj <- data[i, , drop = F]
       } else {
         line.proj <- geosphere::makeLine(rbind(data[i, ], data[i, ], proj[i, ]))
       }
       sphereplot::rgl.sphpoints(line.proj[, 1], line.proj[, 2], radius = 1, col = "black", size = 1)
     }
   }
   fit <- list(prin.curves = prin.curves, line = line, converged = converged, iteration = iter,
               recon.error = rss.new, num.dist.point = distnum)
   return(fit)
 }



  ######################################
  ### Local principal geodesics (Local)
  ######################################

  LPG <- function(data, scale = 0.04, tau = scale/3, nu = 0, maxpt = 500, seed = NULL,
         kernel = "indicator", thres = 1e-4, col1 = "blue", col2 = "green", col3 = "red") {
     # 'Data': a matrix or data.frame of data points represented by angular form of (longitude, latitude).
     # 'scale' is a scale parameter (LPG scaleing region).
     #  For each LPG step, the curve proceed by 'tau'. (step size of the algorithm)
     # 'nu' represents viscosity of the given data.
     # 'thres' is a threshold in the 'IntrinsicMean' function.
     set.seed(seed)                                                        # Set seed number.
     nu <- nu
     h <- scale
     tau <- tau                                                            # Step size of the algorithm.
     r <- 6378137                                                          # Earth radius(m).
     H <- r * h                                                            # scaleing distance on earth(m).
     num.curves <- 0                                                       # The number of the resulting curves.
     # create blank matrix and list
     dim.redu <- line <- line2 <- matrix(NA, nrow = 10000, ncol = 2)
     dim.redu.LPG <- list()
     data.raw <- data
    
     ### start LPG
     repeat{
       if (nrow(data) == 0){
         break
       }

       # first step-------------------------------------------------------------------------------
       num.data <- nrow(data)                                              # The number of data.
       X <- Y <- matrix(nrow = 10000, ncol = 2)                                # Limit nrow set 10000.
       ran.num <- sample(1:num.data, 1)                                    # Choose initial number randomly.
       Y[1, ] <- X[1, ] <- as.matrix(data[ran.num, , drop = F])            # Random sampling X1 from data.
       cov.local <- S <- matrix(c(0, 0, 0, 0), nrow = 2)
       num.related <- ker <- ker.sum <- 0
       direc <- direc2 <- matrix(nrow = 10000, ncol = 2)                       # Limit of repeat is 10000.
       cohesive <- matrix(0, nrow = 2, ncol = 1)
       mean.initial <- c(0, 0, 0)
       adjacent <- data[geosphere::distGeo(X[1, ], data) < H, , drop = F]  # Find the center (intrinsic mean) of adjacent points.
       X[1, ] <- Y[1, ] <- IntrinsicMean(adjacent, rep(1, nrow(adjacent)), thres)

       for (i in 1:num.data){
         dist <- geosphere::distGeo(X[1, ], data[i, ])
         if (dist < H){
           adj <- Rotate(c(X[1, 1], X[1, 2]), c(data[i, 1], data[i, 2]))
           adj.plane <- Logmap(adj)
           if ((kernel == "indicator") | (kernel == "Indicator")){
              ker <- Kernel.indicator(dist/H)
           } else if ((kernel == "quartic") | (kernel == "Quartic")){
              ker <- Kernel.quartic(dist/H)
           } else if ((kernel == "Gaussian") | (kernel == "gaussian")){
              ker <- Kernel.Gaussian(dist/H)
           } else {
              stop("kernal should be either quartic, indicator or Gaussian")
           }
           cohesive <- cohesive + adj.plane * ker                          # Cohesive force.
           S <- S + adj.plane %*% t(adj.plane) * ker
           num.related <- num.related + 1
           ker.sum <- ker.sum + ker
         }
      }
      if (num.related == 0){
        dim.redu.temp <- t(as.matrix(X[1, ]))
        dim.redu  <- rbind(as.matrix(na.omit(dim.redu)), dim.redu.temp)
        line <- rbind(as.matrix(na.omit(line)), dim.redu.temp)
        data <- data[-ran.num, , drop = F]
      } else {
         cov.local <- S/ker.sum                                           # Locally kernel-weighted tangent covariance at scale h.
         eivec <- eigen(cov.local)$vectors[, 1]                           # Leading eigen vector of local tangent covariance.
         if (sum(cohesive^2) > 0){
           cohesive <- cohesive/(sum(cohesive^2))^(1/2)                    # Making cohesive unit vector.
         } else {
           cohesive <- c(0, 0)
         }
         force <- (cohesive * nu + eivec)                                  # Summation of cohesive and maximal variation and 'nu' is viscous parameter.
         direc[1, ] <- force/(sum(force^2))^(1/2)                          # Normalized forward first direction.
         direc2[1, ] <- -direc[1, ]                                        # Normalized backward first direction.
         forward <- tau * direc[1, ]
         backward <- tau * direc2[1, ]
         pt <- Expmap(forward)
         pt.sph <- Trans.sph(pt)
         pt.2 <- Expmap(backward)
         pt.2.sph <- Trans.sph(pt.2)
         Euc <- Rotate.inv(X[1, ], c(pt.sph[1], pt.sph[2]))                # Euclidean points(R^3) of X[2, ].
         Euc2 <- Rotate.inv(X[1, ], c(pt.2.sph[1], pt.2.sph[2]))
         X[2, ] <- Trans.sph(Euc)                                          # Get spherical coordinate (longitude / latitude).
         Y[2, ] <- Trans.sph(Euc2)

      # forward diriection step---------------------------------------------------------------------------
      j <- 2                                                            # Iteration of forward direction step.
      repeat{
        num.related <- ker <- ker.sum <- 0
        cov.local <- S <- matrix(c(0, 0, 0, 0), nrow = 2)
        cohesive <- sum.adj.plane <- c(0, 0)
        for (i in 1:num.data){
          dist <- geosphere::distGeo(X[j, ], data[i, ])
          adj.plane <- c(0, 0)
          if (dist < H){
            adj <- Rotate(c(X[j, 1], X[j, 2]), c(data[i, 1], data[i, 2]))                 # Adjacent data of X[j, ].
            adj.plane <- Logmap(adj)                                                      # Tangent plane coordinate of the adjacent data.
            if ((kernel == "indicator") | (kernel == "Indicator")){
               ker <- Kernel.indicator(dist/H)
            } else if ((kernel == "quartic") | (kernel == "Quartic")){
               ker <- Kernel.quartic(dist/H)
            } else if ((kernel == "Gaussian") | (kernel == "gaussian")){
               ker <- Kernel.Gaussian(dist/H)
            } else{
               stop("kernal should be either quartic, indicator or Gaussian")
            }
            cohesive <- cohesive + adj.plane * ker
            S <- S + adj.plane %*% t(adj.plane) * ker
            ker.sum <- ker.sum + ker
            num.related <- num.related + 1
          }
          sum.adj.plane <- sum.adj.plane + adj.plane
        }
        if ((j > maxpt) | (num.related == 0)){
          break
        } else {
          loc.center <- sum.adj.plane/num.related                                        # Center of adjacent data. (local centering)
          cov.local <- S/ker.sum - (loc.center) %*% t(loc.center)                        # Locally kernel-weighted tangent covariance at scale h.
          if (sum(cohesive^2) > 0){
            cohesive <- cohesive/(sum(cohesive^2))^(1/2)                                  # Making cohesive unit vector.
          } else {
            cohesive <- c(0, 0)
          }
          eivec <- eigen(cov.local)$vectors[, 1]                                          # Leading eigenvector of local tangent covariance.
          force1 <- cohesive * nu + eivec                                                 # Summation of cohesive and maximal variation.
          force2 <- cohesive * nu - eivec                                               
          if (sum(direc[j - 1, ] * eivec) >= 0){
            direc[j, ] <- force1/(sum(force1^2))^(1/2)                                    # Making unit vector.
          } else {
            direc[j, ] <- force2/(sum(force2^2))^(1/2)
          }
          forward <- tau * direc[j, ]
          pt <- Expmap(forward)
          pt.sph <- Trans.sph(pt)
          Euc <- Rotate.inv(X[j, ], c(pt.sph[1], pt.sph[2]))                              # Euclidean point of X[j, ].
          X[j + 1, ] <- Trans.sph(Euc)                                                    # Get a spatial coordinate.
          j <- j + 1
          }
       }
       X <- na.omit(X)

       # backward direction step-----------------------------------------------------------------
       j <- 2                                                                              # Iteration of backward direction step.
       repeat{
          num.related <- ker <- ker.sum <- 0
          cov.local <- S <- matrix(c(0, 0, 0, 0), nrow = 2)
          cohesive <- matrix(0, nrow = 2, ncol = 1)
          sum.adj.plane <- c(0, 0)
          for (i in 1:num.data){
             dist <- geosphere::distGeo(Y[j, ], data[i, ])
             adj.plane <- c(0, 0)
             if (geosphere::distGeo(Y[j, ], data[i, ]) < H){
                adj <- Rotate(c(Y[j, 1], Y[j, 2]), c(data[i, 1], data[i, 2]))                 # Adjacent data of Y[j, ].
                adj.plane <- Logmap(adj)                                                      # Tangent plane coordinate of adjacent data.
                if ((kernel == "indicator") | (kernel == "Indicator")){
                  ker <- Kernel.indicator(dist/H)
                } else if ((kernel == "quartic") | (kernel == "Quartic")){
                  ker <- Kernel.quartic(dist/H)
                } else if ((kernel == "Gaussian") | (kernel == "gaussian")){
                  ker <- Kernel.Gaussian(dist/H)
                } else {
                  stop("kernal should be either quartic, indicator or Gaussian")
                }
                cohesive <- cohesive + adj.plane * ker
                S <- S + adj.plane %*% t(adj.plane) * ker
                ker.sum <- ker.sum + ker
                num.related <- num.related + 1
             }
             sum.adj.plane <- sum.adj.plane + adj.plane 
           }
            if ((j > maxpt) | (num.related == 0)){                                           # Stop the procedure when there is no neighborhood points.
               break
            } else {
              loc.center <- sum.adj.plane/num.related                                       # Center of adjacent data. (local centering)
              cov.local <- S/ker.sum - (loc.center) %*% t(loc.center)                       # Local tangent covariance at scale h.
              eivec <- eigen(cov.local)$vectors[, 1]                                        # Leading eigenvector of local tangent covariance.
              if (sum(cohesive^2) > 0){
                cohesive <- cohesive/(sum(cohesive^2))^(1/2)                                # Making cohesive unit vector.
              } else {
                cohesive <- c(0, 0)
              }
              force1 <- cohesive * nu + eivec                                               # Summation of cohesive and maximal variation.
              force2 <- cohesive * nu - eivec                                     
              if (sum(direc2[j - 1, ] * eivec) >= 0){
                direc2[j, ] <- force1/(sum(force1^2))^(1/2)
              } else {
                direc2[j, ] <- force2/(sum(force2^2))^(1/2)
              }
              backward <- tau * direc2[j, ]
              pt <- Expmap(backward)
              pt.sph <- Trans.sph(pt)
              Euc <- Rotate.inv(Y[j, ], c(pt.sph[1], pt.sph[2]))                             # Euclidean point of Y[j, ].
              Y[j + 1, ] <- Trans.sph(Euc)                                                   # Transform three-dimensional Euclidean coordinate into spherical coordinate (longitude,latitude).
              j <- j + 1
            }
        }
        Y <- na.omit(Y)                                                                    # Deleting NA.
        Y.rev <- Y                                                                         # Y row reverse.
        for (i in 1:nrow(Y.rev)){
           temp <- nrow(Y.rev)
           Y.rev[i, ] <- Y[temp + 1 - i, ]
        }
        dim.redu.temp <- rbind(Y.rev, X)                                                   # Connected componemt of dimension reduction.
        dim.redu <- rbind(as.matrix(na.omit(dim.redu)), dim.redu.temp)                     # Sum of connected component of dimension reduction.
        line <- rbind(as.matrix(na.omit(line)), geosphere::makeLine(dim.redu.temp))        # Line of dim.redu.
        proj <- geosphere::dist2Line(data, dim.redu, distfun = geosphere::distGeo)         # proj[, 2:3] are projection point.
        data <- data[proj[, 1] > H, , drop = F]                                            # Remaining data set.
      }
       dim.redu.LPG[[num.curves + 1]] <- dim.redu.temp                                      # Storage of dim.redu.temp.
       num.curves <- num.curves + 1
     }
     # plot the result
     sphereplot::rgl.sphgrid(col.long = 'black', col.lat = 'black', radaxis = FALSE)                         
     sphereplot::rgl.sphpoints(data.raw, radius = 1, col = col1, size = 12)
     sphereplot::rgl.sphpoints(dim.redu, radius = 1, col = col2, size = 4)
     sphereplot::rgl.sphpoints(line, radius = 1, col = col3, size = 6)
     
     fit <- list(num.curves = num.curves, prin.curves = dim.redu.LPG, line = line)
     return(fit)
     
  }



 