#' Create delta F over F plot for ROI
#'
#' @param obj A target image of Image object or an array.
#' @param ref A reference image of Image object or an array.
#' @export
#' @examples
#' measureROI()
#'

measureROI <- function(img, mask, output, goodfr=F){
  ROI <- dipr::ssweep(img, mask, op="*")
  len <- length(which(mask==1))
  meanint <- colSums(ROI, dim=2)/len
  ROIF0 <- mean(meanint[1:5])
  deltaROIF <- meanint - ROIF0
  ROIdFF0 <- deltaROIF/ROIF0 * 100
  dat <- data.frame(x=1:dim(img)[3], y=ROIdFF0)
  p <- ggplot(data=dat, aes(x=x, y=y)) +
    geom_smooth(method="loess", span = 0.4, level=0.95) +
    ylim(-5, 10)
  ggsave(filename = paste0(output, "_dFF0int_ROI.pdf"), width = 8, height = 8)

  if(goodfr!=F){
    ROIdFF0[!goodfr] <- NA
    dat_goodfr <- data.frame(x=1:dim(img)[3], y=ROIdFF0)
    p <- ggplot(data=dat, aes(x=x, y=y)) +
      geom_smooth(method="loess", span = 0.4, level=0.95) +
      ylim(-5, 10)
    ggsave(filename = paste0(output, "_dFF0int_ROI_goodfr.pdf"), width = 8, height = 8)
  }
  return(ROIdFF0)
}
