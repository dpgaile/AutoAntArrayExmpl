#' Extract Data From GenePix List Object
#'
#' @param GPlist GenePix data list.
#' @param SignalName Name of the GenePix column that will be used to quantify the signal.
#' @param SampleIndices The indices of the samples to be selected.
#'
#' @export
#'
ExtractFromGenePix <- function(GPlist,SignalName,SampleIndices){
   # GPlist=S1g; SignalName="F532.Median"; SampleIndices=1:10

    # load data into convenient containers.
    # only select the features with two replicates
    nonSpcl=which(rowSums(!is.na(GPlist[[1]][,1,]))==2)

    out=list(
        X_raw_1=GPlist[[which(names(GPlist)==SignalName)]][nonSpcl,SampleIndices,1],
        X_raw_2=GPlist[[which(names(GPlist)==SignalName)]][nonSpcl,SampleIndices,2]
    )

    class(out)="AutoAnt"
    return(out)
}

#' Log Transform the Data
#'
#' @param AAdat An Autoantigen Array object
#' @param const The consant to include in the log transfom. cnst=NA means no constant
#' @param delta option to specify delta exactly
#' @param offsets option to supply sample specific offsets
#'
#' @export
#'
LogTransform <- function(AAdat,const=1,delta=NA,offsets=NA){

    X_raw_1=AAdat$X_raw_1
    X_raw_2=AAdat$X_raw_2
    if(!is.na(offsets[1])){
        for(j in 1:(length(offsets))){
            X_raw_1[,j]=X_raw_1[,j]-offsets[j]
            X_raw_2[,j]=X_raw_2[,j]-offsets[j]
        }
        AAdat$offsets=offsets
        AAdat$X_raw_adj_1=X_raw_1
        AAdat$X_raw_adj_2=X_raw_2
        }

    if(is.na(delta)){
        delta=0
        if(!is.na(const)){
            delta=-min(cbind(X_raw_1,X_raw_2))+const
        }
    }

    AAdat$W_1=log(X_raw_1+delta)
    AAdat$W_2=log(X_raw_2+delta)
    AAdat$delta=delta

    return(AAdat)
}

# some of my own quick functions..

#'  Calculate Tukey Tri-mean
#'
#' @param x a vector of values
#' @param type the method of quantile estimation that will be employed
#'
#' @export
#'
TriMean=function(x,type=8) sum(c(1/4,1/2,1/4)*quantile(x,probs=c(0.25,0.5,0.75),type=type))

#'  Calculate Tukey Tri-means for the Rows of a matrix
#'
#' @param X a matrix of values
#' @param type the method of quantile estimation that will be employed
#'
#' @export
#'
RowTriMeans=function(X,type=8) apply(X,1,TriMean)


#' Plot AA data by mean
#'
#' @param AAdat An Autoantigen Array object
#' @param SampleIndices The indices of the samples to be selected.
#' @param clrs colors for the samples. If length 1 then they are all that color
#' @param pchs point type for the samples. If length 1 then they are all that point type
#' @param MeanSubset subset index of samples by which to order the means of the features
#' @param useMat which data matrix to use for the data points. Currently, "log", "raw", and "resid" are available
#' @param cex the size of the plotted points
#'
#' @export
#'
AAplotByMean <- function(AAdat,SampleIndices,clrs=1,pchs=3,MeanSubset=NA, useMat="log", cex=0.5,...){
    # SampleIndices=1:10; clrs=1; pchs=3; MeanSubset=6:10

    nsmpl=length(SampleIndices)
    if(length(clrs)==1) clrs=rep(clrs,nsmpl)
    if(length(pchs)==1) pchs=rep(pchs,nsmpl)

    if(is.na(MeanSubset[1])) MeanSubset=SampleIndices

    if(useMat=="log"){
        W_1=AAdat$W_1
        W_2=AAdat$W_2
        RW=RowTriMeans(cbind(W_1[,MeanSubset],W_2[,MeanSubset]))
        xlbl="TriMean log(Signal)"
        ylbl="log(Signal)"
    }

    if(useMat=="raw"){
        W_1=AAdat$X_raw_1
        W_2=AAdat$X_raw_2
        RW=RowTriMeans(cbind(W_1[,MeanSubset],W_2[,MeanSubset]))
        xlbl="TriMean Signal"
        ylbl="Signal"
    }

    if(useMat=="resid"){
        W_1=AAdat$W_1
        W_2=AAdat$W_2
        RW=RowTriMeans(cbind(W_1[,MeanSubset],W_2[,MeanSubset]))
        # get residual matrices
        RW_MAT=matrix(rep(RW,dim(W_1)[2]),ncol=dim(W_1)[2])
        W_1=W_1-RW_MAT
        W_2=W_2-RW_MAT
        xlbl="TriMean log(Signal)"
        ylbl="log(Signal) - (TriMean log(Signal))"

    }

    plot(rep(RW,2*nsmpl),cbind(W_1[,SampleIndices],W_2[,SampleIndices]),
         type="n",xlab=xlbl,ylab=ylbl,...)

    if(useMat=="resid") abline(h=0,col="steelblue",lwd=3)

    cnt=0
    for(j in SampleIndices){
        cnt=cnt+1
        points(RW,W_1[,j],pch=pchs[cnt],col=clrs[cnt],cex=cex)
        points(RW,W_2[,j],pch=pchs[cnt],col=clrs[cnt],cex=cex)
    }

}


#' Plot Densities of AA data
#'
#' @param AAdat An Autoantigen Array object
#' @param SampleIndices The indices of the samples to be selected.
#' @param clrs colors for the samples. If length 1 then they are all that color
#' @param useMat which data matrix to use for the data points. Currently, "log", "raw", and "resid" are available
#' @param xlims the plotting range for the figure. A value of NA will ustilize the range of the data
#' @param returnModes flag to return estimates of the density modes
#'
#' @export
#'
AAplotDens <- function(AAdat,SampleIndices,clrs=1,useMat="log",xlims=NA, returnModes=T,...){
    # SampleIndices=1:10; clrs=1

    nsmpl=length(SampleIndices)
    if(length(clrs)==1) clrs=rep(clrs,nsmpl)

    if(useMat=="log"){
        W_1=AAdat$W_1
        W_2=AAdat$W_2
        xlbl="TriMean log(Signal)"
    }

    if(useMat=="raw"){
        W_1=AAdat$X_raw_1
        W_2=AAdat$X_raw_2
        xlbl="TriMean Signal"
     }

    if(useMat=="resid"){
        W_1=AAdat$W_1
        W_2=AAdat$W_2
        RW=RowTriMeans(cbind(W_1,W_2))
        # get residual matrices
        RW_MAT=matrix(rep(RW,dim(W_1)[2]),ncol=dim(W_1)[2])
        W_1=W_1-RW_MAT
        W_2=W_2-RW_MAT
        xlbl="Residuals wrt TriMean log(Signal)"
    }

    ylbl="Estimated Density"

    DensRay=array(NA,dim=c(1024,dim(W_1)[2]))
    Modes=rep(NA,dim(W_1)[2])
    Wrnge=range(as.vector(cbind(W_1[,SampleIndices],W_2[,SampleIndices])))
    for(j in SampleIndices){
        tdens=density(c(W_1[,j],W_2[,j]),n=1024,from=Wrnge[1],to=Wrnge[2])
        DensRay[,j]=tdens$y
        Modes[j]=tdens$x[which.max(tdens$y)]
    }

    if(is.na(xlims[1])) xlims=range(tdens$x)

    plot(rep(tdens$x,nsmpl),as.vector(DensRay[,SampleIndices]),
         type="n",xlab=xlbl,ylab=ylbl,xlim=xlims,...)

    if(useMat=="resid") abline(h=0,col="steelblue",lwd=3)

    cnt=0
    for(j in SampleIndices){
        cnt=cnt+1
        lines(tdens$x,DensRay[,j],col=clrs[cnt])
    }

    if(returnModes) return(Modes)

}




#' Center AA data
#'
#' @param AAdat An Autoantigen Array object
#' @param SampleIndices The indices of the samples to be centered.
#' @param AAIndices The indices of the autoantigens used to center the data. If NA, then it is al features
#'
#' @export
#'
AAcenter <- function(AAdat,SampleIndices,AAIndices=NA){
    # SampleIndices=1:5; AAIndices=NA

    nsmpl=length(SampleIndices)
    nAA=dim(AAdat$W_1)[1]

    if(is.na(AAIndices[1])) AAIndices=1:nAA

    W_1=AAdat$W_1[,SampleIndices]
    W_2=AAdat$W_2[,SampleIndices]

    RW=RowTriMeans(cbind(W_1,W_2))
    # get residual matrices
    RW_MAT=matrix(rep(RW,nsmpl),ncol=nsmpl)
    res_1=W_1-RW_MAT
    res_2=W_2-RW_MAT

    # get sample specific centers
    Cnt=RowTriMeans(t(rbind(res_1,res_2)))

    # get corrected matrices..
    Cnt_MAT=matrix(rep(Cnt,nAA),ncol=nsmpl,byrow=T)
    AAdat$W_1[,SampleIndices]=W_1-Cnt_MAT
    AAdat$W_2[,SampleIndices]=W_2-Cnt_MAT

    # now, calculate the new residuals and store them as well
    RW=RowTriMeans(cbind(W_1,W_2))
    # get residual matrices
    RW_MAT=matrix(rep(RW,nsmpl),ncol=nsmpl)

    # if res_1 and res_2 do not exist, then create them
    if(!exists("res_1",where=AAdat)){
        AAdat$res_1=array(NA,dim=dim(AAdat$X_raw_1))
        AAdat$res_2=array(NA,dim=dim(AAdat$X_raw_2))
    }

    AAdat$res_1[,SampleIndices]=AAdat$W_1-RW_MAT
    AAdat$res_2[,SampleIndices]=AAdat$W_2-RW_MAT

    return(AAdat)
}



