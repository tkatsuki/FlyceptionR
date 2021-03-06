#' Analyze trajectory of flies
#'
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' analyze_trajectories()
#'

analyze_trajectories <- function(dir, output, fpsfv, interaction=F){
  fvtrj <- read.table(paste0(dir, list.files(dir, pattern="fv-traj-")))
  trj <- fvtrj[,c(2,3)]
  trja <- read.table(paste0(dir, list.files(dir, pattern="av-traj-")), colClasses = "character")
  trjancol <- ncol(trja)

  # Read trajactory coordinates from av-trj file
  if(trjancol==3|trjancol==4){
    trja <- trja[,2:3]
  }else if(trjancol==6|trjancol==7){
    trja <- trja[,c(2,3,5,6)]
  }else if(trjancol==5) {
    trja <- trja[,c(2,3,4,5)]
  }else {
    stop("Unknown av-trj format")
  }

  trja <- as.data.frame(sapply(trja,gsub,pattern="\\[",replacement=""), stringsAsFactors=F)
  trja <- as.data.frame(sapply(trja,gsub,pattern="\\]",replacement=""), stringsAsFactors=F)
  trja <- sapply(trja, as.numeric)
  headpos <- fvtrj[,c(4,5)]
  edgepos <- fvtrj[,c(6,7)]
  angles <- atan2((headpos - edgepos)[,1], (headpos - edgepos)[,2])
  distance <- dipr::trackDistance(trj)
  distance <- zoo::rollmedian(distance, k=5)
  speed <- zoo::rollsum(distance, k=200)/200*fpsfv
  error <- sqrt(diff(headpos[,1])^2 + diff(headpos[,2])^2)
  png(file=paste0(output, "_error.png"), width=400, height=400)
  plot(error)
  dev.off()

  if(interaction==T){
    flydist <- sqrt((trja[,1] - trja[,3])^2 + (trja[,2] - trja[,4])^2)
    png(file=paste0(output, "_flydist.png"), width=400, height=400)
    plot(flydist, ylab="Distance between flies (mm)")
    dev.off()
  } else {
    flydist <- NA
  }
  return(list("speed"=speed, "error"=error, "trja"=trja, "flydist"=flydist, "angles"=angles))
}
