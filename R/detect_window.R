#' Perform image processing to segment the imaging window
#'
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' detect_window()

detect_window <- function(fvimgl, output, reuse=F){

  message("Performing window detection...")
  if(file.exists(paste0(output, "_fvimgbwbrfh.RDS"))==T &
     file.exists(paste0(output, "_ftrs.RDS"))==T & reuse==T){
    message("Loading RDS file")
    fvimgbwbrfh <- readRDS(paste0(output, "_fvimgbwbrfh.RDS"))
    ftrs  <- readRDS(paste0(output, "_ftrs.RDS"))

  }else{
    # fvimgl: raw fly-view image stack
    # fvimgbw: Find "bright" pixels (px that are 0.1 greater than the 30x30 boxcar filtered version of the image)
    boxcar = 30
    fvimgbw <- thresh(fvimgl, boxcar, boxcar, 0.1)
    writeImage(fvimgl[,,1]/255, file=paste0(output, "_fvimgl.png"))
    writeImage(fvimgbw[,,1], file=paste0(output, "_fvimgbw.png"))

    # fvimgbwc: focus on central circular aperture (100 px)
    centermask <- drawCircle(matrix(0,dim(fvimgl)[1],dim(fvimgl)[2]), dim(fvimgl)[1]/2, dim(fvimgl)[2]/2, 100, col=1, fill=1)
    fvimgbwc <- dipr::ssweep(fvimgbw, centermask, op="*")
    writeImage(fvimgbwc[,,1], file=paste0(output, "_fvimgbwc.png"))
    rm(fvimgbw)

    # fvimgbd: mask of body (use fixed threshold of 200)
    fvimgbd <- (255-fvimgl) > 180
    writeImage(fvimgbd[,,1], file=paste0(output, "_fvimgbd.png"))

    # fvimgbwhd: intersection of fvimgbwc and fvimgbd
    fvimgbwhd <- fvimgbwc * fvimgbd
    writeImage(fvimgbwhd[,,1], file=paste0(output, "_fvimgbwhd.png"))
    rm(fvimgbd)
    rm(fvimgbwc)


    kern1 <- makeBrush(3, shape="diamond")
    fvimgbwhdo <- opening(fvimgbwhd, kern1)
    rm(fvimgbwhd)
    writeImage(fvimgbwhdo[,,1], file=paste0(output, "_fvimgbwhdo.png"))
    fvimgbwhdlb <- bwlabel(fvimgbwhdo) # assign an integer label to each cluster of connected pixels
    fvimgbwbrfh <- fillHull(fvimgbwhdlb) # fill in empty pixels withing each cluster
    rm(fvimgbwhdo)

    # Calculate object size
    message("Calculating window size")
    # fvimgbwbrfh: mask for the largest pixel cluster
    ftrs <- dipr::sfeatures(fvimgbwbrfh, fvimgbwbrfh)
    maxobj <- lapply(ftrs, function(x) x[which(x[, 'm.pxs'] == max(x[, 'm.pxs'])),])
    nonmaxobjid <- lapply(ftrs, function(x) which(x[, 'm.pxs'] != max(x[, 'm.pxs'])))
    fvimgbwbrfh <- rmObjects(fvimgbwbrfh, nonmaxobjid) > 0
    kern1 <- makeBrush(3, shape="diamond")
    fvimgbwbrfh <- dilate(fvimgbwbrfh, kern1)
    saveRDS(fvimgbwbrfh, paste0(output, "_fvimgbwbrfh.RDS"))
    saveRDS(ftrs, paste0(output, "_ftrs.RDS"))
    writeImage(fvimgbwbrfh[,,1], file=paste0(output, "_fvimgbwbrfh.png"))
    rm(fvimgbwhdlb)
  }
  return(fvimgbwbrfh)
}
