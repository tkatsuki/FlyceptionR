#' Detect flash
#'
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' detect_flash()
#'

detect_flash <- function(input, output, type=c("fluo", "fly", "arena"), flash_thresh, reuse=F){
  ## TODO
  ## Adaptive ROI

  if(type=="fluo"){
    if(file.exists(paste0(output, "_flimgint.RDS"))==T & reuse==T){
      message("Loading from RDS file")
      flimgint <- readRDS(paste0(output, "_flimgint.RDS"))
    } else{
    message(sprintf("Reading %s", input))
    flimgint <- dipr::readTIFF2(input, intensity=T)
    png(file=paste0(output, "_flflash.png"), width=400, height=400)
    plot(flimgint)
    dev.off()
    nframesfl <- dipr::readTIFF2(input, getFrames=T)
    message(paste0("Number of frames in fluo-view: ", nframesfl))
    }

    if (flash_thresh==0) {
      # Auto-detect threshold.
      # (Find the largest intensity gap among the brightest 10% of frames.)
      brightestFrames = sort(flimgint, decreasing = TRUE)[1:(length(flimgint)/10)]
      frameDiff = -diff(brightestFrames)
      biggestDiff = which(frameDiff == max(frameDiff))[1]
      flash_thresh = mean(brightestFrames[biggestDiff:(biggestDiff+1)])
      message(sprintf("Flash thresh autodetected as %f", flash_thresh))
    } else {
      message(sprintf("Flash thresh is %f", flash_thresh))
    }


    flflashes <- which(flimgint > flash_thresh)
    # correct for double-flashes from rolling shutter (pick second from each pair)
    flflashes = c(flflashes[which(diff(flflashes) != 1)], tail(flflashes,1))
    flimgflash <- min(flflashes)
    if(flimgflash==Inf) stop("Flash was not detected in fluo-view.")
    message(sprintf("Flash was detected in fluo-view frames: %s", paste(flflashes, collapse=" ")))
    return(list("flflashes"=flflashes, "nframesfl"=nframesfl))
  }

  # Detect flash in fly-view
  if(type=="fly"){
    if(file.exists(paste0(output, "_fvimgsubint.RDS"))==T & reuse==T){
      message("Loading from RDS file")
      fvimgsubint <- readRDS(paste0(output, "_fvimgsubint.RDS"))
      nframesfv <- readRDS(paste0(output, "_nframesfv.RDS"))
    } else{
      message(sprintf("Reading %s", input))
      # Load only diagonal ROIs
      nframesfv <- dipr::readFMF(input, getFrames=T)
      fvimgsub1 <- dipr::readFMF(input, crop=c(5,10,5,10))
      fvimgsub2 <- dipr::readFMF(input, crop=c(220,225,220,225))
      fvimgsubint1 <- colMeans(fvimgsub1, dim=2)
      fvimgsubint2 <- colMeans(fvimgsub2, dim=2)
      rm(fvimgsub1)
      rm(fvimgsub2)
      fvimgsubint <- ifelse(fvimgsubint1 > fvimgsubint2, fvimgsubint1, fvimgsubint2)
      png(file=paste0(output, "_fvflash.png"), width=400, height=400)
      plot(fvimgsubint)
      dev.off()
      saveRDS(fvimgsubint, file=paste0(output, "_fvimgsubint.RDS"))
      saveRDS(nframesfv, file=paste0(output, "_nframesfv.RDS"))
    }

    if (flash_thresh==0) {
      # Auto-detect threshold.
      # (Find the largest intensity gap among the brightest 10% of frames.)
      brightestFrames = sort(fvimgsubint, decreasing = TRUE)[1:(length(fvimgsubint)/10)]
      frameDiff = -diff(brightestFrames)
      biggestDiff = which(frameDiff == max(frameDiff))[1]
      flash_thresh = mean(brightestFrames[biggestDiff:(biggestDiff+1)])
      message(sprintf("Flash thresh autodetected as %f", flash_thresh))
    } else {
      message(sprintf("Flash thresh is %f", flash_thresh))
    }

    fvflashes <- which(fvimgsubint > flash_thresh)
    fvimgflash <- min(fvflashes)
    if(fvimgflash==Inf) stop("Flash was not detected in fly-view.")
    message(sprintf("Flash was detected in fly-view frames: %s", paste(fvflashes, collapse=" ")))
    return(list("fvflashes"=fvflashes, "nframesfv"=nframesfv))
  }

  # Detect flash in arena-view
  if(type=="arena"){
    if(file.exists(paste0(output, "_avimgsubint.RDS"))==T & reuse==T){
      message("Loading from RDS file")
      avimgsubint <- readRDS(paste0(output, "_avimgsubint.RDS"))
      nframesav <- readRDS(paste0(output, "_nframesav.RDS"))
     }else{
      message(sprintf("Reading %s", input))
      nframesav <- dipr::readFMF(input, getFrames=T)
      avimgsub <- dipr::readFMF(input, crop=c(5,10,5,10))
      avimgsubint <- colMeans(avimgsub, dim=2)
      rm(avimgsub)
      png(file=paste0(output, "_avflash.png"), width=400, height=400)
      plot(avimgsubint)
      dev.off()
      saveRDS(avimgsubint, file=paste0(output, "_avimgsubint.RDS"))
      saveRDS(nframesav, file=paste0(output, "_nframesav.RDS"))
    }

    if (flash_thresh==0) {
      # Auto-detect threshold.
      # (Find the largest intensity gap among the brightest 10% of frames.)
      brightestFrames = sort(avimgsubint, decreasing = TRUE)[1:(length(avimgsubint)/10)]
      frameDiff = -diff(brightestFrames)
      biggestDiff = which(frameDiff == max(frameDiff))[1]
      flash_thresh = mean(brightestFrames[biggestDiff:(biggestDiff+1)])
      message(sprintf("Flash thresh autodetected as %f", flash_thresh))
    } else {
      message(sprintf("Flash thresh is %f", flash_thresh))
    }


    avflashes <- which(avimgsubint > flash_thresh)
    avimgflash <- min(avflashes)
    if(avimgflash==Inf) stop("Flash was not detected in arena-view.")
    message(sprintf("Flash was detected in arena-view frames: %s", paste(avflashes, collapse=" ")))
    return(list("avflashes"=avflashes, "nframesav"=nframesav))
  }
}
