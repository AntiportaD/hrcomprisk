#' Bootstrap for estimated CIF
#'
#' This function is based on the CRCumInc from this XX package.
#' Returns an object type dataframe with the following variables and order:
#' col1   col2    col3  col4
#' ttF1.lower  ttF1.upper  ttF2.lower  ttF2.upper
#' @param df A number.
#' @param exit A number.
#' @param event A number.
#' @param exposure A number.
#' @param entry A number.
#' @param weights A number.
#' @param ipwvars A number.
#' @param rep A number.
#' @param seed A number.
#' @param print.attr A number.
#' @importFrom "stats" "approx" "as.formula" "glm" "predict" "quantile" "sd" "stepfun"
#' @return Estimating CIF per event and exposure level
#' @export
bootCRCumInc <- function (df, exit, event, exposure, entry = NULL, weights = NULL, ipwvars=NULL, rep = 0, seed=54321, print.attr=T)
{
  df$exit <- df[[deparse(substitute(exit))]]
  df$event <- df[[deparse(substitute(event))]]
  if(is.factor(df$event)) df$event <- as.numeric(df$event)-1
  df$exposure <- df[[deparse(substitute(exposure))]]
  if(is.factor(df$exposure)) df$exposure <- as.numeric(df$exposure)-1
  df$entry <- df[[deparse(substitute(entry))]]
  df$weights <- df[[deparse(substitute(weights))]]
  set.seed(seed)
  e3 <- any(df$event==3)
  e4 <- any(df$event==4)
  if (rep == 0) {
    stop("\n", "`n` boostrapping repetitions missing", call. = FALSE)
  }
  else {
    nboot <- rep
    bI1o.all <- NULL
    bI1x.all <- NULL
    bI2o.all <- NULL
    bI2x.all <- NULL
    R1.all <- NULL
    R2.all <- NULL
    bI3o.all <- NULL
    bI3x.all <- NULL
    R3.all <- NULL
    bI4o.all <- NULL
    bI4x.all <- NULL
    R4.all <- NULL
    ttx.all <- NULL
    datx <- df[which(df$exposure == 1), ]
    dato <- df[which(df$exposure == 0), ]
    for (b in 1:nboot) {
      b.data <- rbind(datx[sample(1:nrow(datx),replace=T),],dato[sample(1:nrow(dato),replace=T),])
      b.data$exit <- jitter(b.data$exit) #needed???


      #dat_boot <- CRCumInc(b.data, exit, event, exposure, entry, weights, ipwvars, print.attr=F)
      dat_boot <- do.call(CRCumInc,list(b.data,substitute(exit), substitute(event), substitute(exposure), substitute(entry), substitute(weights), substitute(ipwvars), F))
      ttx.all[[b]] <- dat_boot$time
      bI1o.all[[b]] <- dat_boot$CIoinc_1
      bI2o.all[[b]] <- dat_boot$CIoinc_2
      if(e3) bI3o.all[[b]] <- dat_boot$CIoinc_3
      else bI3o.all[[b]] <- 0*bI2o.all[[b]]
      if(e4) bI4o.all[[b]] <- dat_boot$CIoinc_4
      else bI4o.all[[b]] <- 0*bI2o.all[[b]]
      bI1x.all[[b]] <- dat_boot$CIxinc_1
      bI2x.all[[b]] <- dat_boot$CIxinc_2
      if(e3) bI3x.all[[b]] <- dat_boot$CIxinc_3
      else bI3x.all[[b]] <- 0*bI2x.all[[b]]
      if(e4) bI4x.all[[b]] <- dat_boot$CIxinc_4
      else bI4x.all[[b]] <- 0*bI2x.all[[b]]
      R1.all[[b]] <- (1 - (bI2x.all[[b]]+bI3x.all[[b]]+bI4x.all[[b]])/(1 - bI1x.all[[b]]))/(1 - (bI2o.all[[b]]+bI3o.all[[b]]+bI4o.all[[b]])/(1 - bI1o.all[[b]]))
      R1.all[[b]][!is.finite(R1.all[[b]])] <- 99 #needed?
      R2.all[[b]] <- (1 - (bI1x.all[[b]]+bI3x.all[[b]]+bI4x.all[[b]])/(1 - bI2x.all[[b]]))/(1 - (bI1o.all[[b]]+bI3o.all[[b]]+bI4o.all[[b]])/(1 - bI2o.all[[b]]))
      R2.all[[b]][!is.finite(R2.all[[b]])] <- 99
      if(e3){
        R3.all[[b]] <- (1 - (bI2x.all[[b]]+bI1x.all[[b]]+bI4x.all[[b]])/(1 - bI3x.all[[b]]))/(1 - (bI2o.all[[b]]+bI1o.all[[b]]+bI4o.all[[b]])/(1 - bI3o.all[[b]]))
        R3.all[[b]][!is.finite(R3.all[[b]])] <- 99
      }
      if(e4){
        R4.all[[b]] <- (1 - (bI2x.all[[b]]+bI1x.all[[b]]+bI3x.all[[b]])/(1 - bI4x.all[[b]]))/(1 - (bI2o.all[[b]]+bI1o.all[[b]]+bI3o.all[[b]])/(1 - bI4o.all[[b]]))
        R4.all[[b]][!is.finite(R4.all[[b]])] <- 99
      }
      rm(dat_boot)
    }
    ttx <- sort(c(0, df$exit[df$event > 0]))
    R1.lower <- NULL
    R1.upper <- NULL
    R2.lower <- NULL
    R2.upper <- NULL
    if(e3){
      R3.lower <- NULL
      R3.upper <- NULL
    }
    if(e4){
      R4.lower <- NULL
      R4.upper <- NULL
    }
    for (i in 1:length(ttx)) {
      tty1 <- NULL
      tty2 <- NULL
      tty3 <- NULL
      tty4 <- NULL
      for (b in 1:nboot) {
        tty1[b] <- approx(ttx.all[[b]], R1.all[[b]], xout = ttx[i], method = "constant", f = 0)$y
        tty2[b] <- approx(ttx.all[[b]], R2.all[[b]], xout = ttx[i], method = "constant", f = 0)$y
        if(e3) tty3[b] <- approx(ttx.all[[b]], R3.all[[b]], xout = ttx[i], method = "constant", f = 0)$y
        if(e4) tty4[b] <- approx(ttx.all[[b]], R4.all[[b]], xout = ttx[i], method = "constant", f = 0)$y
      }
      ttq <- quantile(tty1, probs = c(0.025, 0.975), na.rm = T)
      R1.lower[i] <- ttq[1]
      R1.upper[i] <- ttq[2]
      ttq <- quantile(tty2, probs = c(0.025, 0.975), na.rm = T)
      R2.lower[i] <- ttq[1]
      R2.upper[i] <- ttq[2]
      if(e3){
        ttq <- quantile(tty3, probs = c(0.025, 0.975), na.rm = T)
        R3.lower[i] <- ttq[1]
        R3.upper[i] <- ttq[2]
      }
      if(e4){
        ttq <- quantile(tty4, probs = c(0.025, 0.975), na.rm = T)
        R4.lower[i] <- ttq[1]
        R4.upper[i] <- ttq[2]
      }
    }
    CRCumInc_boot <- as.data.frame(cbind(R1.lower, R1.upper, R2.lower, R2.upper))
    if(e3) CRCumInc_boot <- as.data.frame(cbind(CRCumInc_boot,R3.lower,R3.upper))
    if(e4) CRCumInc_boot <- as.data.frame(cbind(CRCumInc_boot,R4.lower,R4.upper))
    if(print.attr) print(attributes(CRCumInc_boot))
    invisible(CRCumInc_boot)
  }
}
