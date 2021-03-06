#' CKID dataset
#'
#' A dataset containing time, socieconomic and outcome variables of 626 subjects from the
#' Chronic Kidney Disease in Children (CKiD) Study.
#'
#' @format A data frame with 626 rows and 13 variables:
#' \describe{
#'   \item{b1nb0}{Binary indicator for race: black=1, non-black=0}
#'   \item{entry}{Years since onset of chronic kidney disease at entry into study}
#'   \item{event}{Renal replacement therapy indicator: 0=none, 1=dialysis, 2=transplant}
#'   \item{exit}{Years since onset of chronic kidney disease at event/censoring time}
#'   \item{foodassist}{Binary indicator for use of food assistance}
#'   \item{inckd}{Years in study (=exit-entry)}
#'   \item{incomegt75}{Household income > $75,000 per year}
#'   \item{incomelt30}{Household income < $30,000 per year}
#'   \item{lps}{Binary indicator of low birth weight, premature birth, or small for gestational age}
#'   \item{male1fe0}{Binary indicator for sex: male=1, female=0}
#'   \item{matedultcoll}{Maternal education less than college}
#'   \item{privatemd}{Binary indicator for private doctor}
#'   \item{public}{Binary indicator for public insurance}s
#' }
#' @source \url{https://statepi.jhsph.edu/ckid/ckid.html}
"dat_ckid"
