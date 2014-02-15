#' convert a matrix to a named vector in which the names are taken from the matrix row and column names
#' @export
matrix2vector = function(matrix) {unlist(lapply(1:ncol(x), function(i) lapply(1:nrow(x), function(j) {
  result = x[j,i]
  names(result) = paste(rownames(x)[j], colnames(x)[i], sep="_")
  return(result)
})))
}

#' create a timestamp of the current time in format that can be used in filename 
#' @export
timestampNow <- function() paste(unlist(strsplit(unlist(strsplit(as.character(Sys.time()), " ")), ":")), collapse="-")

# vincenty bearing based on: http://www.r-bloggers.com/great-circle-distance-calculations-in-r/ and http://www.movable-type.co.uk/scripts/latlong-vincenty.html
# Calculates the bearing at starting point in order to reach destination point following geodesic line over WGS-84 ellipsoid using Vincenty method
#' calculate Vincenty bearing
#' @export
vincB <- function(lon1, lat1, lon2, lat2) {
  toRad = pi/180
  toDeg = 180/pi
  # WGS-84 ellipsoid parameters
  a <- 6378137         # length of major axis of the ellipsoid (radius at equator)
  b <- 6356752.314245  # length of minor axis of the ellipsoid (radius at the poles)
  f <- 1/298.257223563 # flattening of the ellipsoid
  
  L <- (lon2-lon1)*toRad # difference in longitude
  U1 <- atan((1-f) * tan(lat1*toRad)) # reduced latitude
  U2 <- atan((1-f) * tan(lat2*toRad)) # reduced latitude
  sinU1 <- sin(U1)
  cosU1 <- cos(U1)
  sinU2 <- sin(U2)
  cosU2 <- cos(U2)
  
  cosSqAlpha <- NULL
  sinSigma <- NULL
  cosSigma <- NULL
  cos2SigmaM <- NULL
  sigma <- NULL
  
  lambda <- L
  lambdaP <- 0
  iterLimit <- 100
  while (abs(lambda-lambdaP) > 1e-12 & iterLimit>0) {
    sinLambda <- sin(lambda)
    cosLambda <- cos(lambda)
    sinSigma <- sqrt( (cosU2*sinLambda) * (cosU2*sinLambda) +
                        (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda) )
    if (sinSigma==0) return(0)  # Co-incident points
    cosSigma <- sinU1*sinU2 + cosU1*cosU2*cosLambda
    sigma <- atan2(sinSigma, cosSigma)
    sinAlpha <- cosU1 * cosU2 * sinLambda / sinSigma
    cosSqAlpha <- 1 - sinAlpha*sinAlpha
    cos2SigmaM <- cosSigma - 2*sinU1*sinU2/cosSqAlpha
    if (is.na(cos2SigmaM)) cos2SigmaM <- 0  # Equatorial line: cosSqAlpha=0
    C <- f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
    lambdaP <- lambda
    lambda <- L + (1-C) * f * sinAlpha *
      (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))
    iterLimit <- iterLimit - 1
  }
  if (iterLimit==0) return(NA)  # formula failed to converge
  
  fwdAz = atan2(cosU2*sinLambda,  cosU1*sinU2-sinU1*cosU2*cosLambda)
  
  return(fwdAz*toDeg) # bearing
}


#' convert seconds into formatted hour:minute:second
#' @export
formatSecs = function(secs){
 hr <- floor(secs/(60*60))
 min <-  floor((secs - (hr*60*60))/60)
 sec <- round(secs - ((hr*60*60)+(min*60)),digits=2)
return(paste(hr,min,sec,sep=':'))
}

