#' Plot Raw Chip for GenePix
#'
#' This function plot GenePix Arrays
#'
#' @param dataset GenePix data.
#' @param coltype color refelcts labeling, green corresponds to IgG; red corresponds to IgM.
#' @param col choose the column you want to look at.
#' @param block choose the block you want to look at closer.
#' @param ttl specific title.
#'
#' @export
#'
QuickPlot <- function(dataset,coltype,col,ttl=""){

    # coltype=ifelse(grepl("g", deparse(substitute(dataset))), "green", "red")
    if(coltype=="green") colfunc<-colorRamp(c("lightgreen","darkgreen"))
    if(coltype=="red") colfunc<-colorRamp(c("red","darkred"))

    rawS = dataset
    x= rawS[,col]
    rx=rank(x)/(length(x)+1)
    rgb_cols=colfunc(rx)
    cols=rgb(rgb_cols,maxColorValue = 256)
    mainplot <- plot(rawS[,6],rawS[,7],
                     main=c(ttl,col),
                     xlab="x-coord",ylab="y-coord",
                     pch=20,cex=0.5,col=cols)

}

#' @export
quickPlot <- function(dataset,coltype,col,block){

    # coltype=ifelse(grepl("g", deparse(substitute(dataset))), "green", "red")
    if(coltype=="green") colfunc<-colorRamp(c("lightgreen","darkgreen"))
    if(coltype=="red") colfunc<-colorRamp(c("red","darkred"))

    rawS = dataset
    x= rawS[,col]
    rx=rank(x)/(length(x)+1)
    rgb_cols=colfunc(rx)
    cols=rgb(rgb_cols,maxColorValue = 256)

    blockplot<- plot(subset(rawS,Block==block)[,6],subset(rawS,Block==block)[,7],
                     main=paste(col," | Block: ",block),
                     xlab="x-coord",ylab="y-coord",
                     pch=15,cex=1.2,col=cols)

}


#' @importFrom shiny runApp
#' @export
quickVis <- function() {
    path <- system.file(package = "AutoAntigen", "QuickPlot")
    runApp(path)
}


#' @importFrom shiny runApp
#' @export
quickVis2 <- function() {
    path <- system.file(package = "AutoAntigen", "QuickPlot2")
    runApp(path)
}



