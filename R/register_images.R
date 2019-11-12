#' Perform image registration for fly-view, window, and fluo-view
#'
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' register_images()
#'

register_images <- function(frid, fvimgl, flimgrt, fvimgbwbrfh, angles, output, zoom, center, cores=1, saveRDS=F, reuse=F){
  message("Performing image registration...")
  if(file.exists(paste0(output, "_regimgi.RDS"))==T &
     file.exists(paste0(output, "_regresi.RDS"))==T & reuse==T){
    message("Loading RDS file")
    regimgi <- readRDS(paste0(output, "_regimgi.RDS"))
    regresi <- readRDS(paste0(output, "_regresi.RDS"))
  }else{
    # Prepare fvimg for registration
    fvimgli <- resize(255 - fvimgl, dim(fvimgl)[1]*zoom)
    centermask <- drawCircle(matrix(0,dim(fvimgli)[1],dim(fvimgli)[2]), dim(fvimgli)[1]/2,
                             dim(fvimgli)[2]/2, dim(fvimgli)[1]/2-1, col=1, fill=1)
    # Create first image, which will be the target in registration
    fvimgrt1sti <- EBImage::rotate(fvimgli[,,1], angles[frid[1]]*180/pi, output.dim=dim(fvimgli)[1:2])
    # Prepare cores
    if(cores!=1){
      library(doParallel)
      registerDoParallel(cores=cores)
      message(paste0("Using ", getDoParWorkers(), " cores..."))
    }

    # Rotate flyview
    rot <- fvimgli
    for (r in 1:dim(rot)[3]){
      rot[,,r] <- RNiftyReg::rotate(fvimgli[,,r], -angles[frid[r]], anchor = c("center"))
    }

    # Run image registration
    regresi <- list()
    if(cores==1){
      for(rg in 1:dim(fvimgli)[3]){
        regresi[[rg]] <- niftyreg(rot[,,rg], fvimgrt1sti,
                                  scope="rigid", symmetric=F)
      }
    }else{
      regresi <- foreach(rg = 1:dim(fvimgli)[3]) %dopar% niftyreg(fvimgli[,,rg], fvimgrt1sti,  scope="rigid", symmetric=F, internal=FALSE)
    }
    regimgi <- array(sapply(regresi, function(x) x$image), dim=dim(fvimgli))
    regimgi[which(is.na(regimgi)==T)] <- 0
    if(saveRDS==T){
      saveRDS(regimgi, paste0(output, "_regimgi.RDS"))
      saveRDS(regresi, paste0(output, "_regresi.RDS"))
    }
    writeImage((255-regimgi)/255, file=paste0(output, "_regimgi.tif"))
    rm(fvimgli)
  }

  # Apply affine transformation to fluo-view
  message("Transforming fluo-view...")
  if(file.exists(paste0(output, "_flimgreg.RDS"))==T & reuse==T){
    message("Loading RDS file")
    flimgreg <- readRDS(paste0(output, "_flimgreg.RDS"))

  }else{

    # Pad flimgrt to match the size of fvimg
    if(dim(regimgi)[1] > dim(flimgrt)[1]){
      flimgpad <- regimgi*0
      flimgpad[round((dim(regimgi)[1]-dim(flimgrt)[1])/2):
                 (round((dim(regimgi)[1]-dim(flimgrt)[1])/2)+dim(flimgrt)[1]-1),
               round((dim(regimgi)[2]-dim(flimgrt)[2])/2):
                 (round((dim(regimgi)[2]-dim(flimgrt)[2])/2)+dim(flimgrt)[2]-1),
               1:dim(flimgrt)[3]] <- flimgrt

    }else{
      flimgpad <- flimgrt[round((dim(flimgrt)[1]-dim(regimgi)[1])/2):
                            (round((dim(flimgrt)[1]-dim(regimgi)[1])/2)+dim(regimgi)[1]-1),
                          round((dim(flimgrt)[2]-dim(regimgi)[2])/2):
                            (round((dim(flimgrt)[2]-dim(regimgi)[2])/2)+dim(regimgi)[2]-1),]
    }
    flimgpadmv <- translate(flimgpad, center)

    # Rotate fluoview
    flimgrot <- flimgpadmv
    for (r in 1:dim(flimgrot)[3]){
      flimgrot[,,r] <- RNiftyReg::rotate(flimgpadmv[,,r], -angles[frid[r]], anchor = c("center"))
    }

    flimgrgres <- list()
    if(cores==1){
      for(app in 1:dim(fvimgl)[3]){
        flimgrgres[[app]] <- applyTransform(forward(regresi[[app]]), flimgrot[,,app])
      }
    }else{
      flimgrgres <- foreach(app = 1:dim(fvimgl)[3]) %dopar%  applyTransform(forward(regresi[[app]]), flimgrot[,,app])
    }
    flimgreg <- array(unlist(flimgrgres), dim=dim(flimgrot))
    writeImage(normalize(flimgreg), file=paste0(output, "_flimgreg.tif"))
    rm(flimgrgres)
    rm(flimgpad)
    rm(flimgpadmv)
    if(saveRDS==T){
      saveRDS(flimgreg, paste0(output, "_flimgreg.RDS"))
    }
  }

  # Apply affine transformation to window mask
  message("Transforming window images...")
  if(file.exists(paste0(output, "_fvimgbwbrfhregimg.RDS"))==T & reuse==T){
    message("Loading RDS file")
    fvimgbwbrfhregimg <- readRDS(paste0(output, "_fvimgbwbrfhregimg.RDS"))
  }else{

    fvimgbwbrfhrs <- resize(fvimgbwbrfh, dim(fvimgl)[1]*zoom)
    rm(fvimgbwbrfh)

    rotwin <- fvimgbwbrfhrs
    for (rwin in 1:dim(rotwin)[3]){
      rotwin[,,rwin] <- RNiftyReg::rotate(fvimgbwbrfhrs[,,rwin], -angles[frid[rwin]], anchor = c("center"))
    }


    fvimgbwbrfhreg <- list()
    if(cores==1){
      for(app in 1:dim(fvimgl)[3]){
        fvimgbwbrfhreg[[app]] <- applyTransform(forward(regresi[[app]]), rotwin[,,app])
      }
    }else{
      fvimgbwbrfhreg <- foreach(app = 1:dim(fvimgl)[3]) %dopar% applyTransform(forward(regresi[[app]]), rotwin[,,app])
    }
    fvimgbwbrfhregimg <- array(unlist(fvimgbwbrfhreg), dim=dim(fvimgbwbrfhrs))
    fvimgbwbrfhregimg <- fvimgbwbrfhregimg >= 0.5
    writeImage((255-regimgi)/255*(1-fvimgbwbrfhregimg) + normalize(flimgreg),
               file=paste0(output, "_fvfvwindflregimg.tif"))
    rm(fvimgbwbrfhrs)
    rm(fvimgbwbrfhreg)
    rm(regresi)
    if(saveRDS==T){
      saveRDS(fvimgbwbrfhregimg, paste0(output, "_fvimgbwbrfhregimg.RDS"))
    }
  }
  return(list("regimgi"=regimgi, "flimgreg"=flimgreg, "fvimgbwbrfhregimg"=fvimgbwbrfhregimg))
}
