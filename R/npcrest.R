#' Wrapper Function
#'
#' Main function
#' See help for other functions: check.packages, checking_data, data_CIF, fig_cif
#' @param df A number.
#' @param exit A number.
#' @param event A number.
#' @param exposure A number.
#' @param entry A number.
#' @param weights A number.
#' @param ipwvars A number.
#' @param maxtime A number.
#' @param rep A number.
#' @param eoi A number.
#' @importFrom "grDevices"  "gray"
#' @importFrom "graphics" "abline" "axis" "box" "lines" "mtext" "par" "plot" "text"
#' @importFrom "stats" "approx" "as.formula" "glm" "predict" "quantile" "sd" "stepfun"
#' @return npcrest
#' @export
npcrest <- function (df, exit, event, exposure, entry = NULL, weights = NULL, ipwvars=NULL, maxtime = Inf, rep = NULL, eoi = -1)
{
  datcheck(df, deparse(substitute(exit)), deparse(substitute(event)), deparse(substitute(exposure)), deparse(substitute(entry)), deparse(substitute(weights)), ipwvars, eoi = eoi)
  #Ensure that EOI matches its factor level
  evt <- df[[deparse(substitute(event))]]
  if(is.factor(evt) & eoi>0) eoi <- (as.numeric(evt[evt == eoi])-1)[1]
  myCIF <- do.call(CRCumInc,list(df,substitute(exit), substitute(event), substitute(exposure), substitute(entry), substitute(weights), substitute(ipwvars)))

  if (!is.null(rep)) {
    boot_myCIF <- do.call(bootCRCumInc,list(df, substitute(exit), substitute(event), substitute(exposure), substitute(entry), substitute(weights), substitute(ipwvars), rep))
    plotCIF(myCIF, maxtime = maxtime, ci = boot_myCIF, eoi = eoi)
  }
  if (is.null(rep)) {
    plotCIF(myCIF, maxtime = maxtime, ci = NULL, eoi = eoi)
  }
  invisible(myCIF)
}
