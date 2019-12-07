#' Data Check for Wrapper Function
#'
#' Data needs to be of class: dataframe
#' Data in Long format (1 row per ID)
#' Package needs at LEAST colnames without any specific order:
#' 		 "exposure" : binary (0, 1), numeric
#' 		 "entry"	: continuous, numeric
#' 		 "exit"		: continuous, numeric
#' 		 "event"	: binary (>0), numeric
#' @param df A number.
#' @param qexit A number.
#' @param qevent A number.
#' @param qexposure A number.
#' @param qentry A number.
#' @param qweights A number.
#' @param qipwvars A number.
#' @param eoi A number.
#' @return Check dataset

datcheck<- function(df, qexit, qevent, qexposure, qentry, qweights, qipwvars, eoi = -1){
  errors  <-character()
  changes <-character()
  if(!is.data.frame(df)) {
    msg <- paste("Data is not a data frame class object.\n")
    errors <- c(errors, msg)
  }
  if (!(length(errors) == 0)) {
    stop("\n", errors)
  }
  if(!(qexposure %in% colnames(df))) {
    msg <- paste(qexposure, "not in dataset\n")
    errors <- c(errors, msg)
  }
  if(qexposure %in% colnames(df)) {
    if(!is.numeric(df[[qexposure]]) &  !is.factor(df[[qexposure]])) {
      msg <- paste("Exposure needs to be numeric or factor\n")
      errors <- c(errors, msg)
    }
    if (length(unique(df[[qexposure]]))!=2) {
      msg <- paste("Exposure does not have 2 levels\n")
      errors <- c(errors, msg)
    }
  }
  if(!(qevent %in% colnames(df))) {
    msg <- paste(qevent, "not in dataset\n")
    errors <- c(errors, msg)
  }
  if(qevent %in% colnames(df)) {
    if (!(is.numeric(df[[qevent]]) | is.factor(df[[qevent]]))) {
      msg <- paste("Event needs to be numeric or factor\n")
      errors <- c(errors, msg)
    }
    if(length(unique(df[[qevent]]))>=5) {
      msg <- paste("Event has more than 4 levels\n")
      errors <- c(errors, msg)
    }
  }
  # Recode factor and tell people about each level?
  if(!(qexit %in% colnames(df))) {
    msg <- paste(qexit, "not in dataset\n")
    errors <- c(errors, msg)
  }
  if (qentry!="NULL"){
    if(!(qentry %in% colnames(df))) {
      msg <- paste(qentry, "not in dataset\n")
      errors <- c(errors, msg)
    }
  }
  if (qweights!="NULL"){
    if(!(qweights %in% colnames(df))) {
      msg <- paste(qweights, "not in dataset\n")
      errors <- c(errors, msg)
    }
  }
  if (length(qipwvars)>0){
    wtvars<-unique(qipwvars)
    for (w in wtvars){
      if(!(w %in% colnames(df))) {
        msg <- paste(deparse(substitute(w)), "not in dataset\n")
        errors <- c(errors, msg)
      }
    }
  }
  if (eoi!=-1){
    if(!(eoi %in% df[[qevent]])) {
      msg <- paste("eoi", eoi, "not an event type in", qevent, "\n")
      errors <- c(errors, msg)
    }
  }
  #Stop Message Output
  if (!(length(errors) == 0)) {
    stop("\n", errors, call. = FALSE)
  }
  #Modifications
  if(is.factor(df[[qexposure]])) {
    expt<-unique(sort(as.numeric(df[[qexposure]])-1))
    msg.txt <- paste(qexposure,":\t",
                     levels(df[[qexposure]])[1], "\t", "\t",expt[1], "\n",
                     "\t","\t",levels(df[[qexposure]])[2], "\t", "\t", expt[2], "\n")
    changes <- c(changes, msg.txt)
  }
  if(is.factor(df[[qevent]])){
    evt<-unique(sort(as.numeric(df[[qevent]])-1))
    msg.evt	<- paste(qevent,":\t",
                     levels(df[[qevent]])[1], "\t", "\t",evt[1], "\n",
                     "\t", "\t",levels(df[[qevent]])[2], "\t", "\t", evt[2], "\n",
                     "\t", "\t",levels(df[[qevent]])[3], "\t", "\t", evt[3], "\n")
    if(length(unique(df[[qevent]]))>3){
      evt.2<-unique(evt[evt>2])
      msg.evt1<-character()
      for(i in evt.2){
        msg.evtx<- paste("\t", "\t",
                         levels(df[[qevent]])[i], "\t", "\t",evt[i], "\n",
                         msg.evt1<-c(msg.evt1, msg.evtx))
      }
      changes <- c(changes, msg.evt, msg.evt1)
    }
    else{
      changes <- c(changes, msg.evt)
    }
  }
  #Modifications message output
  if (!(length(changes) == 0)) {
    message("\n", "*Modifications\n",
            "Var\t", "Factor Val\t",
            "Numeric Val\n", changes, "***", "\n")
  }
}
