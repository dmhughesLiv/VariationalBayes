########## R-function: summMCMC ##########

# Does a compact graphical summary of Markov chain Monte Carlo output.

# Last changed: 04 MAY 2018 by M.P.Wand.

summMCMC <- function(xList,EPSfileName,PDFfileName,plotInd=1,
                     parNames,columnHeadings,colourVersion=TRUE,
                     credLevel=0.95,columnHeadCex=NULL,
                     paletteNum=1,numerSummCex=NULL,BGRsttPos=10,
                     BGRyRange=c(0.95,1.25),BGRtickPos=1.2,
                     BGRlogTransf,BGRlogitTransf,KDExlim,KDEvertLine=TRUE,
                     KDEvertLineCol="black",addTruthToKDE=NULL)
{
   # Define required functions:

   empty.panel <- function()
   {
      plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xaxt="n",
           yaxt="n",xlab="",ylab="",bty="o")
      invisible()
   }

   logit <- function(x)
   {
      if (any((x<=0)|(x>=1))) stop("argument out of range")
      return(log(x/(1-x)))
   }
  
   # Set dimension and other variables:

   options(warn=-1)
   x <- lapply(xList,as.matrix)
   numPar <- ncol(x[[plotInd]])
   sampSize <- nrow(x[[plotInd]])
   numChains <- length(x)

   # Divert figure to EPS or PDF file if filename specified:

   if (!missing(EPSfileName))
   {
      if (is.null(columnHeadCex)) columnHeadCex <- 2.4
      if (is.null(numerSummCex)) numerSummCex <- 1.0
      postscript(EPSfileName,horizontal=FALSE,width=11,height=(numPar+1))
   }
       
   if (!missing(PDFfileName))
   {
      if (is.null(columnHeadCex)) columnHeadCex <- 2.9
      if (is.null(numerSummCex)) numerSummCex <- 1.3
      pdf(PDFfileName,width=9,height=(numPar+1))
   }

   if (missing(EPSfileName)&missing(PDFfileName))
   {
      if (is.null(columnHeadCex)) columnHeadCex <- 2.9
      if (is.null(numerSummCex)) numerSummCex <- 1.3
   }

   if (numChains>1)
   {
      columnHeadCex <- 6*columnHeadCex/7
      numerSummCex <- 6*numerSummCex/7
   }    
   
   if (missing(columnHeadings))
      columnHeadings <- c("parameter","trace","lag 1","acf","BGR",
                           "density","summary")   

   if (colourVersion)
   {
      if (paletteNum==1)
      {
         columnCols <- c("purple4","tomato","mediumblue",
                         "olivedrab4","DarkCyan","darkred","DarkGreen")
         if (KDEvertLine&missing(KDEvertLineCol))
            KDEvertLineCol <- "darkgoldenrod"
      }
      if (paletteNum==2)
      {
         columnCols <- c("darkmagenta", "green4","darkorange","dodgerblue",
                         "darkgoldenrod1","red","navy")
         if (KDEvertLine&missing(KDEvertLineCol))
            KDEvertLineCol <- "DarkGreen"
      }
   }

   panels.per.par <- length(columnHeadings)
   if (numChains==1) panels.per.par <- panels.per.par - 1

   if (!colourVersion)
   {
      columnCols <- rep(NA,numPar)
      for (i in 1:length(columnHeadings))
         columnCols[i] <- "black"
   }

   op <- par()
   par(mfrow=c((numPar+1),panels.per.par))

   headInds <- 1:7
   if (numChains==1) headInds <- c(1:4,6,7)
   for (i in 1:panels.per.par)
   {
      par(ann=F,mar=rep(0,4),xaxt="n",yaxt="n",xpd=TRUE)
      empty.panel()
      text(0.5,0.5,columnHeadings[headInds[i]],cex=columnHeadCex,
           col=columnCols[headInds[i]])
   }

   for (j in 1:numPar)
   {
      # Write the variable name:

      par(ann=F,mar=c(0,0,0,0),xaxt="n",yaxt="n")

      empty.panel()
      if (length(parNames[[j]])==1)
         text(0.5,0.5,parNames[[j]][1],cex=3.0,col=columnCols[1])
      if (length(parNames[[j]])==2)
      {
         text(0.5,0.7,parNames[[j]][1],cex=1.9,col=columnCols[1])
         text(0.5,0.3,parNames[[j]][2],cex=1.9,col=columnCols[1])
      }
      if (length(parNames[[j]])==3)
      {
         text(0.5,0.8,parNames[[j]][1],cex=1.55,col=columnCols[1])
         text(0.5,0.5,parNames[[j]][2],cex=1.55,col=columnCols[1])
         text(0.5,0.2,parNames[[j]][3],cex=1.55,col=columnCols[1])
      }

      # Do the trace plot:

      plot(x[[plotInd]][,j],xlab="",ylab="",type="l",col=columnCols[2])
      
      # Do the lag 1 plot:

      plot(x[[plotInd]][1:(sampSize-1),j],
           x[[plotInd]][2:sampSize,j],xlab="",ylab="",type="n")

      points(x[[plotInd]][1:(sampSize-1),j],
             x[[plotInd]][2:sampSize,j],pch=1,cex=0.5,col=columnCols[3])
      
      # Do the autocorrelation function plot:

      ci.col.val <- "black"
      if (colourVersion) ci.col.val <- "blue"
      acf(x[[plotInd]][,j],lag.max=20,col=columnCols[4],
          lwd=2,ci.col=ci.col.val)
      
      if (numChains>1) # Multiple chains.
      {

         # Do the Gelman-Rubin plot

         x.list <- list(x[[1]][,j])
            for (k in 2:numChains)
               x.list <- c(x.list,list(x[[k]][,j]))

         if ((!missing(BGRlogTransf)))
         {
            if (any(j==BGRlogTransf))
            {
               x.list <- list(log(x[[1]][,j]))
                  for (k in 2:numChains)
                     x.list <- c(x.list,list(log(x[[k]][,j])))
            }
         }
         if ((!missing(BGRlogitTransf)))
         {
            if (any(j==BGRlogitTransf))
            {
               x.list <- list(log(x[[1]][,j]))
                  for (k in 2:numChains)
                     x.list <- c(x.list,list(logit(x[[k]][,j])))
            }
         }

         BGR.outp <- BGRinterval(x.list,BGRsttPos)

         # Extract x-axis and y-axis components of `BGR.outp'.
         
         BGR.x <- BGR.outp$x
         BGR.y <- BGR.outp$numer/BGR.outp$denom
 
         lb.wt <- 0.20*length(BGR.x)
         plot(0,0,type="n",xlim=c(-lb.wt,length(BGR.x)),ylim=BGRyRange,
              xlab="",ylab="", bty="l")
         lines(BGR.x,BGR.y,lwd=1,col=columnCols[5],err=-1)
         lines(rep(0,2),BGRyRange,err=-1)

         text(-0.75*lb.wt,BGRtickPos,as.character(BGRtickPos),srt=90,cex=0.8)
         text(-0.75*lb.wt,1,"1.00",srt=90,cex=0.8)
         abline(h=1,lwd=1)
         lines(c(-lb.wt/2,0),rep(BGRtickPos,2))
      }

      # Do the density plot

      h <- dpik(x[[plotInd]][,j])
      est <- bkde(x[[plotInd]][,j],bandwidth=h)

      if (!missing(KDExlim))
         xlim.val <-  KDExlim[[j]]   

      if (missing(KDExlim))
         xlim.val <-  range(est$x)

      lb.ht <- 0.2*(max(est$y)-min(est$y))
      ylim.val <- c(-lb.ht,1.1*max(est$y))
      plot(est,type="l",xlab="",ylab="",ylim= ylim.val,
           col=columnCols[6],xlim=xlim.val,lwd=2)
      lines(c(min(est$x),max(est$x)),rep(0,2))

      if (KDEvertLine)
         lines(rep(0,2),c(0,max(est$y)),lwd=2,err=-1,col=KDEvertLineCol)
      
      if (!is.null(addTruthToKDE))
         lines(rep(addTruthToKDE[j],2),c(0,max(est$y)),
               lwd=1,lty=2,,err=-1,col=KDEvertLineCol)
     
      
      x.labels <- pretty(x[[plotInd]][,j],n=3)
      for (i in 1:length(x.labels))
      {
         text(x.labels[i],-0.6*lb.ht,as.character(x.labels[i]),cex=0.8)
         lines(rep(x.labels[i],2),c(0,-0.1*lb.ht))
      }   

      empty.panel()
      text(0.5,0.75,
      paste("posterior mean: ",
             as.character(signif(mean(x[[plotInd]][,j]),3)),sep=""),
             col=columnCols[7],cex=numerSummCex)
      text(0.5,0.5,paste(100*credLevel,"% credible interval: ",
           sep=""),col=columnCols[7],cex=numerSummCex)

      text(0.5,0.25,
      paste("(",
      as.character(signif(quantile(x[[plotInd]][,j],(1-credLevel)/2),3)),",",
      as.character(signif(quantile(x[[plotInd]][,j],(1+credLevel)/2),3)),")",
      sep=""),col=columnCols[7],cex=numerSummCex)
   }

   par(op)

   invisible()
}

########## End of summMCMC ##########
