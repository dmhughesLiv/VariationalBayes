########## R script: ISDR Model ##########

# For conducting a Bayesian Joint mixed
# model analysis for PBC Data

# Last changed: 22 April 2021

# Set tasks:
doMCMC <- TRUE
doMFVB <- TRUE
doMCMCvsMFVB <- TRUE
doSummMCMC<-TRUE

# Set MCMC sample size paramaters:

nWarm <- 5000
nKept <- 10000
nThin <- 10

# Load required packages:

library(lattice)    ; library(MASS)
# Generate data:
library(mixAK)
data(PBCseq)
PBCseq$lbili<-scale(PBCseq$lbili)
PBCseq$lalbumin<-scale(PBCseq$lalbumin)
PBCseq$lalk.phos<-scale(PBCseq$lalk.phos)
PBCseq$lchol<-scale(PBCseq$lchol)
PBCseq$lsgot<-scale(PBCseq$lsgot)
PBCseq$lprotime<-scale(PBCseq$lprotime)
PBCseq$lplatelet<-scale(PBCseq$lplatelet)
PBCseq$age<-scale(PBCseq$age)
PBCseq$month<-scale(PBCseq$month)

set.seed(34567)
# Set hyperparameters:
sigmaBeta <- 1e5 ; AU <- 1e5 ; Aeps1 <- 1e5;

if (doMCMC)
{
   # Perform MCMC:
   MCMCStart<-Sys.time()
   MCMCmod<-GLMM_MCMC(y = PBCseq[,c("lbili","lalbumin","lalk.phos","lchol","lsgot","lprotime","lplatelet","ascites","hepatom","spiders")],                             
                      dist=c("gaussian","gaussian","gaussian","gaussian","gaussian","gaussian","gaussian","binomial(logit)","binomial(logit)","binomial(logit)"),                                               
                      id=PBCseq[,"id"],					                                        
                      x=list(lbili="empty",lalbumin="empty",lalk.phos="empty",lchol="empty",lsgot="empty",lprotime="empty",lplatelet="empty",ascites=PBCseq[,"month"],
                             hepatom=PBCseq[,"month"],spiders=PBCseq[,"month"]),    
                      z=list(lbili=PBCseq[,"month"],lalbumin=PBCseq[,"month"],lalk.phos=PBCseq[,"month"],lchol=PBCseq[,"month"],lsgot=PBCseq[,"month"],
                             lprotime=PBCseq[,"month"],lplatelet=PBCseq[,"month"],ascites="empty",hepatom="empty",spiders="empty"),
                      random.intercept = c(lbili = TRUE,lalbumin=TRUE,lalk.phos=TRUE,lchol = TRUE,lsgot=TRUE,lprotime = TRUE,lplatelet=TRUE,ascites=TRUE,hepatom=TRUE,
                                           spiders=TRUE),                          
                      prior.b = list(Kmax = 1, xi = rep(0, 17), D = diag(rep(sigmaBeta, 17)), zeta = 18, gD = rep(0.5, 17), hD = rep(4/AU^2, 17)),
                      prior.eps = list(zeta = rep(1,7), g = rep(1,7), h = rep(1/Aeps1^2,7)),
                      nMCMC = c(burn = nWarm, keep = nKept, thin = nThin, info = 1000),
                      PED=FALSE)
   MCMCEnd<-Sys.time()
   MCMCTime<-MCMCEnd-MCMCStart
   
   # Save and extract relevant MCMC samples:
P<-MCMCmod$p+MCMCmod$q+1
N<-sum(MCMCmod$R)                   ##Number of markers
Eps<-MCMCmod$R[1]                   ##Number of Gaussian Error Terms.
Q<-MCMCmod$dimb
Q1<-MCMCmod$q+1
CumQ<-c(1,cumsum(Q1)+1)
CumP<-c(1,cumsum(MCMCmod$p)+1)
parms<-list()

##Assign Fixed Effects
MCMCOut<-matrix()
if(length(MCMCmod$alpha)>0){for(i in 1:N){
                                          MCMCOut<-cbind(MCMCOut,MCMCmod$mixture_b[,CumQ[i]:(CumQ[(i+1)]-1)])
                                          if(MCMCmod$p[i]>0){MCMCOut<-cbind(MCMCOut,MCMCmod$alpha[,CumP[i]:(CumP[(i+1)]-1)])}
                                         }}
if(length(MCMCmod$alpha)==0){for(i in 1:N){MCMCOut<-cbind(MCMCOut,MCMCmod$mixture_b[,CumQ[i]:(CumQ[(i+1)]-1)])}}
MCMCOut<-MCMCOut[,-1]

FP1<-rep(1:N,times=P)
FP2<-sequence(P)
for(i in 1:length(FP1)){parms[[paste("beta",FP1[i],".",FP2[i],"MCMC",sep="")]]<-MCMCOut[,i]}

##Assign Error Terms 
for(k in 1:Eps){parms[[paste("sigmaEpsMCMC",k,sep="")]]<-MCMCmod$sigma_eps[,k]}

##Assign Random Effects Terms
##Add Diagonals to the parms list
t1<-vector(length=Q)
for(s in 1:Q){
t1[s]<-(s*((Q+1)-.5*(s-1)))
parms[[paste("Sigma",s,".",s,"MCMC",sep="")]]<-MCMCmod$mixture_b[,t1[s]]^2
}

##Add off-diagonals to parms list.
Pos1<-rep(1:((Q)-1),c(((Q)-1):1))
Pos2<-NULL
for(i in 1:((Q)-1)){Pos2<-c(Pos2,(i+1):(Q))}
for(s1 in 1:length(Pos1)){
parms[[paste("Sigma",Pos1[s1],".",Pos2[s1],"MCMC",sep="")]]<-MCMCmod$mixture_b[,-c(1:(Q),t1)][,s1]*MCMCmod$mixture_b[,t1[Pos1[s1]]]*MCMCmod$mixture_b[,t1[Pos2[s1]]]
}

##Now create list of names
parNamesVal<-list()
for(i in 1:length(FP1)){parNamesVal[[i]]<-parse(text=paste("beta[list(",FP1[i],",",FP2[i],")]",sep=""))}
for(i in (length(FP1)+1):(length(FP1)+Eps)){parNamesVal[[i]]<-parse(text=paste("sigma[epsilon[",i-length(FP1),"]]^2",sep=""))}
for(i in (length(FP1)+Eps+1):(length(FP1)+Eps+(Q))){parNamesVal[[i]]<-parse(text=paste("Sigma[list(",i-(length(FP1)+Eps),",",i-(length(FP1)+Eps),")]",sep=""))}
for(i in 1:length(Pos1)){parNamesVal[[i+(length(FP1)+Eps+(Q))]]<-parse(text=paste("Sigma[list(",Pos1[i],",",Pos2[i],")]",sep=""))}

save(MCMCmod,MCMCTime,parms,parNamesVal,file=paste("PBC",N,"MCMCmod.RData",sep=""))
  }

if (doSummMCMC)
{
   if(!doMCMC){
   load(paste("PBC",10,"MCMCmod.RData",sep=""))
   P<-MCMCmod$p+MCMCmod$q+1
   N<-sum(MCMCmod$R)                   ##Number of markers
   Eps<-MCMCmod$R[1]                   ##Number of Gaussian Error Terms.
   Q<-MCMCmod$dimb
   Q1<-MCMCmod$q+1
   CumQ<-c(1,cumsum(Q1)+1)
   CumP<-c(1,cumsum(MCMCmod$p)+1)
   FP1<-rep(1:N,times=P)
   FP2<-sequence(P)
   Pos1<-rep(1:((Q)-1),c(((Q)-1):1))
   Pos2<-NULL}
    
   # Do parameters plot:
   source("summMCMC.r")
   library(KernSmooth)
   parms1<-list(do.call(cbind,parms))

   ##Fixed Effects Summaries
   for(i in 1:ceiling(sum(P)/10)){
   summMCMC(list(parms1[[1]][,(1+10*(i-1)):min((10*i),sum(P))]),parNames = parNamesVal[(1+10*(i-1)):min((10*i),sum(P))],EPSfileName=paste("PBCMod",N,"FixedEffects",i,".eps",sep=""))
   dev.off()}

   ##ErrorSummaries
   summMCMC(list(parms1[[1]][,(length(FP1)+1):(length(FP1)+Eps)]),parNames = parNamesVal[(length(FP1)+1):(length(FP1)+Eps)],EPSfileName=paste("PBCMod",N,"ErrorTerms",".eps",sep=""))

   ##Random Effects Summaries
   D<-Q*(Q+1)/2   ##Number of random effects parameters.
   StPt1<-length(parms)-D
   for(i in 1:ceiling(D/10)){
   summMCMC(list(parms1[[1]][,StPt1+(1+10*(i-1)):min((10*i),D)]),parNames = parNamesVal[StPt1+(1+10*(i-1)):min((10*i),D)],EPSfileName=paste("PBCMod",N,"RandomEffects",i,".eps",sep=""))
   dev.off()}
}

if (doMFVB)
{
library(magic)
source("mglmm_new.R")

modVB<-mglmmVB(y = PBCseq[,c("lbili","lalbumin","lalk.phos","lchol","lsgot","lprotime","lplatelet","ascites","hepatom","spiders")],
               x=list(lbili=cbind(rep(1,nrow(PBCseq)),PBCseq[,c("month")]),lalbumin=cbind(rep(1,nrow(PBCseq)),PBCseq[,c("month")]),
                      lalk.phos=cbind(rep(1,nrow(PBCseq)),PBCseq[,c("month")]),lchol=cbind(rep(1,nrow(PBCseq)),PBCseq[,c("month")]),    
                      lsgot=cbind(rep(1,nrow(PBCseq)),PBCseq[,c("month")]),lprotime=cbind(rep(1,nrow(PBCseq)),PBCseq[,c("month")]),    
                      lplatelet=cbind(rep(1,nrow(PBCseq)),PBCseq[,c("month")]),ascites=cbind(rep(1,nrow(PBCseq)),PBCseq[,c("month")]),    
                      hepatom=cbind(rep(1,nrow(PBCseq)),PBCseq[,c("month")]),spiders=cbind(rep(1,nrow(PBCseq)),PBCseq[,c("month")])),    
               z=list(lbili=cbind(rep(1,nrow(PBCseq)),PBCseq[,"month"]),lalbumin=cbind(rep(1,nrow(PBCseq)),PBCseq[,"month"]),lalk.phos=cbind(rep(1,nrow(PBCseq)),PBCseq[,"month"]),
                      lchol=cbind(rep(1,nrow(PBCseq)),PBCseq[,"month"]),lsgot=cbind(rep(1,nrow(PBCseq)),PBCseq[,"month"]),lprotime=cbind(rep(1,nrow(PBCseq)),PBCseq[,"month"]),
                      lplatelet=cbind(rep(1,nrow(PBCseq)),PBCseq[,"month"]),ascites=cbind(rep(1,nrow(PBCseq))),hepatom=cbind(rep(1,nrow(PBCseq))),
                      spiders=cbind(rep(1,nrow(PBCseq)))),
               id=PBCseq[,"id"],responseType=c("Gaussian","Gaussian","Gaussian","Gaussian","Gaussian","Gaussian","Gaussian","Bernoulli","Bernoulli","Bernoulli"),
               method="KMW",doStreamlined=TRUE,maxIter=500)

save(modVB,file=paste("PBC",N,"MFVBmod.RData",sep=""))
}
          
if (doMCMCvsMFVB)
{
        #if(!doMCMC){load(paste("PBC",N,"mod.RData",sep=""))}
    source("accVarApp.r")
    source("trapint.r")
    mu.q.betauG <-modVB$mu.q.betauG
    Sigma.q.betauG <- modVB$Sigma.q.betauG
    B.q.SigmauR <- modVB$B.q.SigmaR
    A.q.sigsq.eps <- modVB$A.q.sigsq.eps
    B.q.sigsq.eps <- modVB$B.q.sigsq.eps

###Accuracy Assessment for ISDR data
libname <- c("MASS", "magic", "Matrix", "rstan", "gdata", "lattice", "splines","KernSmooth", "stringr")
vapply(libname, require, character.only = TRUE, 0)

ng <- 101      # Grid size.
lwdVal <- 2    # Line width.
cexVal <- 1  # Font size.
colVAval <- "darkorange" ; colMCval <- "blue"
  colTruthVal <- "green4" ; colAccVal <- "purple4"
  colBoxplots <- "lightblue"
  cex.mainVal <- 2.5 ; cex.accurVal <- 2.5; cex.labVal <- 2.5
  cex.axisVal <- 2.5 ; cexVal <- 2.5


##Fixed effects Accuracy Assessment
parms<-list(do.call(cbind,parms))
accBetax<-rep(NA,sum(P))


for(i in 1:sum(P)){
accBetax[i] <- accVarApp(c(mu.q.betauG[i], Sigma.q.betauG[i,i]),
                        parNamesVal[i],                        
                        parms[[1]][,i], "Normal",
                        plotDensities=FALSE,colVA = colVAval, colMC = colMCval,
                        colTruth = colTruthVal, colAcc = colAccVal,
                        cex.mainVal = cex.mainVal, cex.accur = cex.accurVal,
                        cex.labVal = cex.labVal)
}

##Error Term Accuracy
accSigsqEps<-rep(NA,Eps)
for(i in (1+sum(P)):(Eps+sum(P))){
accSigsqEps[i-sum(P)] <- accVarApp(c(A.q.sigsq.eps[i-sum(P)], B.q.sigsq.eps[i-sum(P)]),
                             parNamesVal[i],
                             parms[[1]][,i]^2,"InvGamma",
                             plotDensities=FALSE,colVA = colVAval, colMC = colMCval,
                             colTruth = colTruthVal, colAcc = colAccVal,
                             cex.mainVal = cex.mainVal, cex.accur = cex.accurVal,
                             cex.labVal = cex.labVal)
}

##Covariance Matrix Accuracy Assessment.
##No longer demonstrate the traditional Accuracy plots. Instead calculate Accuracy and display on a heat map with covariance means alongside!
Pos1<-rep(1:((Q)-1),c(((Q)-1):1))
Pos2<-NULL
for(i in 1:((Q)-1)){Pos2<-c(Pos2,(i+1):(Q))}


##Covariance Matrix Off Diagonal Accuracy...
m<-length(modVB$id)
StartPoint<-length(parNamesVal)-length(Pos1)
A.q.SigmaR12 <- (2 + Q - 1) + m
A.q.SigmauR11 <- ((2 + Q - 1 + m) - Q + 1) / 2
accSigmaR<-NULL

for(i in 1:length(Pos1)){
 accSigmaR[i] <- accVarApp(list(A.q.SigmaR12,B.q.SigmauR,Pos1[i],Pos2[i]),
                           parNamesVal[(StartPoint+i)],
                           parms[[1]][,(StartPoint+i)], "OffDiagWishart",
                           plotDensities=FALSE,colVA = colVAval, colMC = colMCval,
                           colTruth = colTruthVal, colAcc = colAccVal,
                           cex.mainVal = cex.mainVal, cex.accur = cex.accurVal,
                           cex.labVal = cex.labVal)
}

##Covariance matrix Diagonal Accuracy...
StartPoint<-length(parNamesVal)-length(Pos1)-Q
accSigmaRDiag<-NULL
for(i in 1:Q){
   accSigmaRDiag[i] <- accVarApp(c(A.q.SigmauR11, B.q.SigmauR[i,i] / 2),
                           parNamesVal[(StartPoint+i)],
                           parms[[1]][,(StartPoint+i)], "InvGamma",
                           plotDensities=FALSE,colVA = colVAval, colMC = colMCval,
                           colTruth = colTruthVal, colAcc = colAccVal,
                           cex.mainVal = cex.mainVal, cex.accur = cex.accurVal,
                           cex.labVal = cex.labVal)
}

#Heat Maps
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
namesU<-c("log(Bilirubin):Int","log(Bilirubin):Slope","log(Albumin):Int","log(Albumin):Slope","log(alk.phos):Int","log(alk.phos):Slope",
"log(Cholesterol):Int","log(Cholesterol):Slope","log(ALT):Int","log(ALT):Slope","log(Clotting):Int","log(Clotting):Slope","Platelets:Int",
"Platelets:Slope","Ascites:Int","Hepatom:Int","Spiders:Int")

m1 <- matrix(, Q, Q)
m1[lower.tri(m1, diag=FALSE)] <- accSigmaR
m2 <- t(m1)
m2[lower.tri(m2, diag=FALSE)] <- accSigmaR
diag(m2)<-accSigmaRDiag
colnames(m2)<-namesU[1:Q]
rownames(m2)<-namesU[1:Q]

#Get upper triangle of the correlation matrix
  get_upper_tri <- function(m2){
    m2[lower.tri(m2)]<- NA
    return(m2)
  }

upper_tri <- get_upper_tri(m2)
melted_m2 <- melt(upper_tri, na.rm = TRUE)

AccPlot <- ggplot(melted_m2, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white",show.legend=TRUE)+
 scale_fill_gradient2(low = "red", high = "green", mid = "white", 
   midpoint = 50, limit = c(0,100), space = "Lab", 
    name="Accuracy (%)") +
  theme_minimal()+ # minimal theme
 ggtitle("(b) Covariance Matrix Accuracy (%)") +
 theme(plot.title = element_text(size=12, face="bold",hjust = 0.5)) +
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 10, hjust = 1))+
 theme(axis.text.y = element_text(angle = 0, vjust = 1, 
    size = 10, hjust = 1))+
 coord_fixed() + 
geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.65, 0.7),
  legend.direction = "horizontal",
        plot.margin=unit(c(-0.1,-0.9,0.5,-0.3), "cm"))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))


##Similar Idea for Covariance...
CovMat<-B.q.SigmauR/m
colnames(CovMat)<-namesU[1:Q]
rownames(CovMat)<-namesU[1:Q]

##Convert to correlation matrix
CorrMat<-round(cov2cor(CovMat),digits=2)
upper_tri <- get_upper_tri(CorrMat)
melted_CorrMat <- melt(upper_tri, na.rm = TRUE)

CorrPlot <- ggplot(melted_CorrMat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white",show.legend=TRUE)+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
    name="Correlation") +
  theme_minimal()+ # minimal theme
 ggtitle("(c) Correlation Matrix") +
 theme(plot.title = element_text(size=12, face="bold",hjust = 0.5)) +
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 10, hjust = 1))+
# theme(axis.text.y = element_text(angle = 0, vjust = 1, 
#    size = 10, hjust = 1))+
 theme(axis.text.y = element_blank())+
 coord_fixed() + 
geom_text(aes(Var2, Var1, label = value), color = "black", size = 2.5) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.65, 0.7),
  legend.direction = "horizontal",
        plot.margin=unit(c(-0.1,-0.3,0.5,-0.9), "cm"))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))

namesR<-c("log(Bilirubin)","log(Albumin)","log(alk.phos)","log(Cholesterol)","log(ALT)","log(Clotting)","Platelets","Ascites","Hepatom","Spiders")
namesP<-c("Intercept","Time","Std.Dev")
BetaAcc<- matrix(accBetax,nrow=P[1],ncol=N,byrow=FALSE)
EpsAcc<- matrix(accSigsqEps,nrow=Eps,ncol=1,byrow=TRUE)
BetaAcc1<-rbind(BetaAcc,c(EpsAcc,NA,NA,NA))
colnames(BetaAcc1)<-namesR[1:N]
rownames(BetaAcc1)<-namesP[1:(P[1]+1)]
BetaAcc1<-BetaAcc1[nrow(BetaAcc1):1,]

BetaAccPlot <- ggplot(melt(BetaAcc1), aes(Var2, Var1, fill = value))+
 geom_tile(color = "white",show.legend=TRUE)+
 scale_fill_gradient2(low = "red", high = "green", mid = "white", 
   midpoint = 50, limit = c(0,100), space = "Lab", 
    name="Accuracy (%)") +
  theme_minimal()+ # minimal theme
 ggtitle("(a) Fixed Effects Accuracy (%)") +
 theme(plot.title = element_text(size=14, face="bold",hjust = 0.5)) +
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 theme(axis.text.y = element_text(angle = 0, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() + 
geom_text(aes(Var2, Var1, label = value), color = "black", size = 5) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  #plot.margin=unit(c(1,-0.5,1,1), "cm")
  plot.margin=unit(c(0.1,1,0.2,1), "cm"),
  axis.ticks = element_blank(),
  legend.position = "none",
  )

postscript(paste("PBC",N,"AccuracyPlotRE.eps",sep=""))
grid.arrange(AccPlot, CorrPlot, ncol=2,nrow=1,widths=c(1,1))
dev.off()

postscript(paste("PBC",N,"AccuracyPlotFE.eps",sep=""))
grid.arrange(BetaAccPlot)
dev.off()


g<-grid.arrange(arrangeGrob(BetaAccPlot,ncol=1),arrangeGrob(AccPlot,CorrPlot,ncol=2),nrow=2,ncol=1,widths=1,heights=c(0.9,2))
ggsave(file="PBC10MarkerAccuracy.eps", g) 

StreamlineAcc<-c(accBetax,accSigsqEps,accSigmaR,accSigmaRDiag)
Times<-c(MCMCTime,modVB$VBTime)
save(StreamlineAcc,Times,file=paste("PBC",N,"Comp.RData",sep=""))
}
   
      
########## End of PBCmod ##########

