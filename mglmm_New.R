##edit the input to the function to allow for different kinds of input.
mglmmVB <- function(y, x, z, id, responseType = "Gaussian",method="JJ", doStreamlined = TRUE, maxIter = 35, useMatForDg = TRUE) {
##First set up design matrices etc   
##Maybe need to sort and order by id (just in case)
m <- length(unique(id))  ##Number of patients
y<-as.data.frame(y)  ##TO allow for only 1 marker case.
R<-ncol(y)    ##Number of marker
RC<-sum(responseType=="Gaussian")
RP<-sum(responseType=="Poisson")
RB<-sum(responseType=="Bernoulli")
rc<-which(responseType=="Gaussian")
rp<-which(responseType=="Poisson")
rb<-which(responseType=="Bernoulli")

#lambda <- function(x) {
#       nzi <- (1:length(x))[x!=0]
#       ans <- rep(0.125,length(x))
#       ans[nzi] <- tanh(x[nzi]/2)/(4*x[nzi])
#       return(ans)
#     }
#     phi <- function(x) {
#       nzi <- (1:length(x))[x!=0]
#       ans <- rep(0.125,length(x))
#       ans[nzi] <- (x[nzi]/2) - log(1+exp(x[nzi])) + (0.25*x[nzi]*tanh(x[nzi]/2))
#       return(ans)
#     }

reBlockInds <-lapply(unique(id),function(u){which(id%in%u)})
nVec<-sapply(reBlockInds,function(u){length(u)})
numObs <- sum(nVec)  # Total number of observations.
CumNVec<-c(1,R*cumsum(nVec)+1)
 
##need to count the number of columns of x and z.
pr<-qr<-NULL
for(r in 1:R){
pr[r]<-ncol(as.data.frame(x[[r]]))
qr[r]<-ncol(as.data.frame(z[[r]]))
} 

Q<-sum(qr)  ##Total Number of random effects.
P<-sum(pr)  ##Total Number of fixed effects.
CumPr<-c(1,cumsum(pr)+1)  ##Identifiers for where each block of X and Z should start and finish (ie corresponds to patient and marker rows.
CumQr<-c(1,cumsum(qr)+1)
X<-matrix(0,ncol=P,nrow=numObs*R)
Z<-matrix(0,ncol=Q,nrow=numObs*R)
Y<-vector(length=R*numObs)

if(!doStreamlined){
Z1<-matrix(0,ncol=Q*m,nrow=numObs*R)
for(i in 1:m){
 for(r in 1:R){
  X[(CumNVec[i]+(r-1)*nVec[i]):(CumNVec[i]+(r*nVec[i])-1),CumPr[r]:(CumPr[(r+1)]-1)]<-as.matrix(x[[r]])[reBlockInds[[i]],]
  Z[(CumNVec[i]+(r-1)*nVec[i]):(CumNVec[i]+(r*nVec[i])-1),CumQr[r]:(CumQr[(r+1)]-1)]<-as.matrix(z[[r]])[reBlockInds[[i]],]
  Z1[(CumNVec[i]+(r-1)*nVec[i]):(CumNVec[i]+(r*nVec[i])-1),(Q*(i-1)+CumQr[r]):(Q*(i-1)+CumQr[(r+1)]-1)]<-as.matrix(z[[r]])[reBlockInds[[i]],]
  Y[(CumNVec[i]+(r-1)*nVec[i]):(CumNVec[i]+(r*nVec[i])-1)]<-as.numeric(y[reBlockInds[[i]],r])
  }
 }
}else{
for(i in 1:m){
 for(r in 1:R){
  X[(CumNVec[i]+(r-1)*nVec[i]):(CumNVec[i]+(r*nVec[i])-1),CumPr[r]:(CumPr[(r+1)]-1)]<-as.matrix(x[[r]])[reBlockInds[[i]],]
  Z[(CumNVec[i]+(r-1)*nVec[i]):(CumNVec[i]+(r*nVec[i])-1),CumQr[r]:(CumQr[(r+1)]-1)]<-as.matrix(z[[r]])[reBlockInds[[i]],]
  Y[(CumNVec[i]+(r-1)*nVec[i]):(CumNVec[i]+(r*nVec[i])-1)]<-as.numeric(y[reBlockInds[[i]],r])
 }
}
}

id1<-rep(unique(id),times=nVec*R)  ##Vector to keep track of id's
Riden<-unlist(sapply(nVec,function(u){rep(colnames(y),each=u)}))

DataCheck<-ifelse(!doStreamlined,dta<-data.frame(Y=Y,id1=id1,Riden=Riden,X=X,Z=Z,Z1=Z1),dta<-data.frame(Y=Y,id1=id1,Riden=Riden,X=X,Z=Z))
##Need to remove any rows with NA's
is.NA <- as.logical(apply(is.na(dta), 1, sum))
Y<-Y[!is.NA]
X<-X[!is.NA,]
Z<-Z[!is.NA,]
if(!doStreamlined){Z1<-Z1[!is.NA,]} 
id1<-id1[!is.NA]
Riden<-Riden[!is.NA]
m<-length(unique(id1))  ##Redefine m in case any patients have dropped out with NA's
##Need to retrack ID's and Marker rows...
  reBlockInds <-lapply(unique(id1),function(u){which(id1%in%u)})
  nVec<-sapply(reBlockInds,function(u){length(u)})
  numObs <- sum(nVec)  # Total number of observations.
  CumNVec<-c(1,R*cumsum(nVec)+1)

RBlockInds <- vector("list", length = R)
nRVec<-NULL
for(r in 1:R){
RBlockInds[r]<-list(which(Riden==colnames(y)[r]))
nRVec[r]<-sum(Riden==colnames(y)[r])
}

## Set up dimension variables for design matrices.
ncZ <- Q*m
##Set up hyperparameter choices. (In principle should have a line in the input to adjust these as necessary)
nuVal <- 2
sigsq.beta <- 1e5
A.epsr <- rep(1e5,RC); A.k <- rep(1e5,Q)

## Set adjusted response vector.
yAdj <- Y
Z_aslist<-lapply(1:m,function(u){Z[reBlockInds[[u]],,drop=F]})

   if (!doStreamlined) {
      C <- cbind(X,Z1);         CTC <- crossprod(C)
      ncC <- ncol(C)
   } else {
     ## Create lists of matrices required for MFVB.
     ZtX<-vector("list",length=m)
     ZtZ<-vector("list",length=m)
     for(i in 1:m){ZtX[[i]]<-crossprod(Z_aslist[[i]],X[reBlockInds[[i]],,drop=F])
                   ZtZ[[i]]<-crossprod(Z_aslist[[i]],Z_aslist[[i]])}
     G <- vector("list",length=m)
     H <- vector("list",length=m)
     GH<- vector("list",length=m)
     }

   ## Do MFVB initialisations. (Could have input in function to alter these if desirable)
   mu.q.recip.sigsq.epsr <- rep(1,RC)
   mu.q.recip.Sigma.eps<-vector(length=numObs)
   mu.q.recip.a.epsr <- rep(1,RC)
   M.q.inv.SigmaR <- diag(Q)
   mu.q.recip.ak <- rep(1,Q)
   if(!doStreamlined){mu.q.betau<-c(rep(1,P),rep(0,m*Q))}
   mu.q.betauG<-rep(1,P)
   mu.q.uR<-vector("list",length=m)
   for(i in 1:m){mu.q.uR[[i]]<-rep(0,Q)}
   if(!doStreamlined){Sigma.q.betau<-diag((P+m*Q))}
   Sigma.q.betauG<-diag(P)
   Sigma.q.uR<-vector("list",length=m)
   for(i in 1:m){Sigma.q.uR<-diag(Q)}
   A.q.sigsq.epsr<-vector(length=RC)
   for(r in rc){A.q.sigsq.epsr[r]<-0.5*(length(RBlockInds[[r]])+1)}
   B.q.sigsq.epsr<-vector(length=RC)
   wtVec1 <- rep(1,length(Y))
   wtVec2 <- rep(1,length(Y))
   #xiVec <- rep(1,length(Y))
   

##Set up Monahan and Stefanski p and s values..
pMS<-c(0.003246343272134,0.051517477033972,0.195077912673858,0.315569823632818,0.274149576158423,0.131076880695470,0.027912418727972,0.001449567805354) 
sMS<-c(1.365340806296348,1.059523971016916,0.830791313765644,0.650732166639391,0.508135425366489,0.396313345166341,0.308904252267995,0.238212616409306)

##Set up Ematrices to identify types of marker..
E1<-E2<-E3<-rep(0,numObs)
for(r in rc){E1[RBlockInds[[r]]]<-1}
for(r in rp){E2[RBlockInds[[r]]]<-1}
for(r in rb){E3[RBlockInds[[r]]]<-1}


   ## Set up function for the lower bound on the marginal log-likelihood.
   logML <- function(nuVal,R,P,ncZ,Q,m,numObs,mu1,E2,wtP,E3,B0muSig,B1muSig,xiVec,phiVec,
                     sigsq.beta,mu.q.betauG,Sigma.q.betauG,det.Sigma.q.betau,
                     B.q.SigmaR,A.k,M.q.inv.SigmaR,
                     mu.q.recip.ak,B.q.a.k,A.epsr,B.q.sigsq.epsr,B.q.a.epsr,
                     mu.q.recip.a.epsr,mu.q.recip.sigsq.epsr,mu.q.recip.Sigma.eps) {
      ASigma <- nuVal + Q - 1
      logCqRA <- (0.5*ASigma*Q*log(2) + 0.25*Q*(Q-1)*log(pi) +
                  sum(lgamma(0.5*(ASigma+1-(1:Q)))))
      logCqRAm <- (0.5*(ASigma+m)*Q*log(2) + 0.25*Q*(Q-1)*log(pi) +
                  sum(lgamma(0.5*(ASigma+m+1-(1:Q)))))
      ans1 <- (0.5*Q*ASigma*log(2*nuVal)
               - (0.5*Q+RC)*log(pi) - 0.5*P*log(sigsq.beta))
      if (!doStreamlined) {
         ans2 <- (-1/(2*sigsq.beta)*(sum(mu.q.betau[1:P]^2) + sum(diag(Sigma.q.betau[(1:P),(1:P)]))))
      } else {
         ans2 <- (-1/(2*sigsq.beta)*(sum(mu.q.betauG^2) + sum(diag(Sigma.q.betauG))))
      }
      ans3 <- (0.5*det.Sigma.q.betau + 0.5*(ncZ+P))  
      ans4 <- (-logCqRA + logCqRAm -0.5*(ASigma+m)*determinant(B.q.SigmaR)$modulus)
      ans5 <- (-sum(log(A.k)) + Q*lgamma(0.5*(nuVal+Q)))
      ans6 <- (sum(nuVal*diag(M.q.inv.SigmaR)*mu.q.recip.ak))
      ans7 <- (-0.5*(nuVal+Q)*sum(log(B.q.a.k)))
      ans8 <- (-0.5*sum((nRVec[1:RC]+1)*log(B.q.sigsq.epsr)) + sum(lgamma(0.5*(nRVec[1:RC]+1))) - 0.5*sum(nRVec[1:RC])*log(2*pi))
      ans9 <- (-sum(log(A.epsr)) - sum(log(B.q.a.epsr)) + sum(mu.q.recip.a.epsr*mu.q.recip.sigsq.epsr))
      ifelse(RP>0,ans10 <- sum(yAdj*E2*mu1) -sum(lgamma(E2*yAdj+1)) - sum(E2*wtP),ans10<-0) 
      ifelse(RB>0,ifelse(method=="KMW",ans11 <- sum(yAdj*E3*mu1) -sum(E3*B0muSig),ans11<-sum((yAdj-0.5)*E3*mu1) - E3*sum((B1muSig/2)*(xiVec^2)) + E3*sum(phiVec)),ans11<-0)
      ans <- ans1+ans2+ans3+ans4+ans5+ans6+ans7+ans8+ans9+ans10+ans11
         return(list(ans1,ans2,ans3,ans4,ans5,ans6,ans7,ans8,ans9,ans10,ans11,ans))
   }

   logMLcurr <- -1e20
   converged <- FALSE ; itnum <- 0
   tolerance <- 0.000001
   logMLgrid <- NULL
   ridgeVec <- rep((1/sigsq.beta),P)

StartTime<-Sys.time()
   while (!converged) {
     itnum <- itnum + 1
     ## Set current log(ML) value to previous one.
     if (itnum == 1) logMLprev <- logMLcurr
     if (itnum > 1)  logMLprev <- logMLcurr[[length(logMLcurr)]][1]
     ## Update mu.q.recip.Sigma.eps
       for(r in 1:R){ifelse(responseType[r]=="Gaussian",mu.q.recip.Sigma.eps[RBlockInds[[r]]]<-mu.q.recip.sigsq.epsr[r],mu.q.recip.Sigma.eps[RBlockInds[[r]]]<-1)}
       
    if (!doStreamlined) {
       ## Set up covariance matrix M.q.Sigma.
       if (Q == 1) {
         ridgeVec2 <- rep(M.q.inv.SigmaR,m)
         M.q.Sigma <- diag(c(ridgeVec,ridgeVec2))
       } else if (Q > 1) {
         ridgeMat <- kronecker(diag(m),M.q.inv.SigmaR)
         M.q.Sigma <- adiag(diag(ridgeVec),ridgeMat)
       }

       ## Update q*(beta,u) parameters.
nu.q.betau<-crossprod(C,(mu.q.recip.Sigma.eps*(yAdj-wtVec1))) - M.q.Sigma%*%mu.q.betau
Sigma.q.betau<-solve(crossprod(C,(mu.q.recip.Sigma.eps*wtVec2*C))+ M.q.Sigma)
mu.q.betau<- mu.q.betau + Sigma.q.betau%*%nu.q.betau

mu1<-C%*%mu.q.betau
##Need to calculate B0muSig and B1muSig if there are any Bernoulli variables using Monahan and Stefanski method.

if(method=="KMW"){
OmMS<-sqrt(tcrossprod(rep(1,numObs),rep(1,8)) + tcrossprod(diag(tcrossprod(C%*%Sigma.q.betau,C)),sMS^2))
B0muSig<-pnorm(tcrossprod(mu1,sMS)/OmMS)%*%pMS
B1muSig<-(dnorm(tcrossprod(mu1,sMS)/OmMS)/OmMS)%*%(pMS*sMS)
}

#if(method=="JJ"){
#B1muSig <- 2*lambda(xiVec)
#B0muSig <- B1muSig*mu1 -0.5
#phiVec <- phi(xiVec)
#}

##Need to define wtVec1 = diag(b'(sigmax+mu))
wtP<-exp(mu1 + .5*rowSums((C%*%Sigma.q.betau)*C))
wtVec1<-as.vector(E1*mu1 + E2*wtP + E3*B0muSig)
##Need to define wtVec2 = diag(b''(sigmax+mu))
wtVec2<-as.vector(E1 + E2*wtP + E3*B1muSig)
       ## Compute the determinant of Sigma.q.betau.
       det.Sigma.q.betau <- determinant(Sigma.q.betau)$modulus

       ##Update xiVec
       #if(method=="JJ"){
       #EsqMat <- Sigma.q.betau + tcrossprod(mu.q.betau)
       #xiVec <- sqrt(diag(C%*%EsqMat%*%t(C)))
       #                }
     }

    if (doStreamlined) {
       ## Update q*(beta,u) parameters.

       sVec <- rep(0,P); Smat <- matrix(0,P,P)
       for (i in 1:m) {
         G[[i]] <- crossprod(X[reBlockInds[[i]],,drop=F],(mu.q.recip.Sigma.eps[reBlockInds[[i]]]*wtVec2[reBlockInds[[i]]]*Z_aslist[[i]]))
         H[[i]] <- solve(crossprod(Z_aslist[[i]],(mu.q.recip.Sigma.eps[reBlockInds[[i]]]*wtVec2[reBlockInds[[i]]]*Z_aslist[[i]])) + M.q.inv.SigmaR)
         GH[[i]]<- G[[i]]%*%H[[i]]
         sVec <- sVec + as.vector(GH[[i]]%*%((M.q.inv.SigmaR%*%mu.q.uR[[i]])-crossprod(Z_aslist[[i]],(mu.q.recip.Sigma.eps[reBlockInds[[i]]]*(yAdj[reBlockInds[[i]]]-wtVec1[reBlockInds[[i]]])))))
         Smat <- Smat + GH[[i]]%*%t(G[[i]])
         }

       Sigma.q.betauG <- solve(crossprod(X,(mu.q.recip.Sigma.eps*wtVec2*X)) + diag(ridgeVec) - Smat)
       mu.q.betauG.OLD<-mu.q.betauG
       mu.q.uR.OLD<-mu.q.uR
       Sigma.q.betauG.OLD<-Sigma.q.betauG
       Sigma.q.uR.OLD<-Sigma.q.uR
       mu.q.betauG <- mu.q.betauG + as.vector(Sigma.q.betauG%*%(crossprod(X,(mu.q.recip.Sigma.eps*(yAdj-wtVec1))) - ridgeVec*mu.q.betauG + sVec))
        
       Sigma.q.uR <- vector("list",length=m)
       for (i in 1:m) {Sigma.q.uR[[i]] <- (H[[i]] + H[[i]]%*%crossprod(G[[i]],(Sigma.q.betauG%*%GH[[i]])))}

for (i in 1:m) {
         mu.q.uR[[i]] <- mu.q.uR[[i]] + as.vector(H[[i]]%*%(crossprod(Z_aslist[[i]],(mu.q.recip.Sigma.eps[reBlockInds[[i]]]*(yAdj[reBlockInds[[i]]]-wtVec1[reBlockInds[[i]]]))) 
                         - crossprod(G[[i]],(mu.q.betauG-mu.q.betauG.OLD)) - M.q.inv.SigmaR%*%mu.q.uR[[i]]))
       }


       f1<-function(x){!is.finite(x)}
       if(sum(f1(mu.q.betauG))!=0 | sum(sapply(mu.q.uR,f1))!=0){mu.q.betauG<-mu.q.betauG.OLD;mu.q.uR<-mu.q.uR.OLD;Sigma.q.betauG<-Sigma.q.betauG.OLD;Sigma.q.uR<-Sigma.q.uR.OLD;break}

## Compute the determinant of Sigma.q.betau.
      detAMat <- NULL
       for (i in 1:m)
         detAMat <- c(detAMat,determinant(solve(H[[i]]))$modulus)

       detBMat <- matrix(0,P,P)
       for (i in 1:m)
         detBMat <- detBMat + GH[[i]]%*%t(G[[i]])

       detBMat <- crossprod(X,(mu.q.recip.Sigma.eps*wtVec2*X)) + diag(ridgeVec) - detBMat
       det.Sigma.q.betau <- -(sum(detAMat) + determinant(detBMat)$modulus)

##Define mu1 and sigma1 for the change of variables
mu1a<-X%*%mu.q.betauG
mu1b<-NULL
mu1b<-do.call(rbind,lapply(1:m,function(u){Z_aslist[[u]]%*%mu.q.uR[[u]]}))
mu1<-mu1a+mu1b

sigma1<-rowSums((X%*%Sigma.q.betauG)*X)  ##need to check why this works. Much faster but what is the algebra?
for(i in 1:m){
             sigma1[reBlockInds[[i]]]<-sigma1[reBlockInds[[i]]] - 2*rowSums((Z_aslist[[i]]%*%t(GH[[i]])%*%Sigma.q.betauG)*X[reBlockInds[[i]],,drop=F])+ rowSums((Z_aslist[[i]]%*%Sigma.q.uR[[i]])*Z_aslist[[i]])
             } 

if(method=="KMW"){
OmMS<-sqrt(tcrossprod(rep(1,numObs),rep(1,8)) + tcrossprod(sigma1,sMS^2))
B0muSig<-pnorm(tcrossprod(mu1,sMS)/OmMS)%*%pMS
B1muSig<-(dnorm(tcrossprod(mu1,sMS)/OmMS)/OmMS)%*%(pMS*sMS)
}

#if(method=="JJ"){
#B1muSig <- 2*lambda(xiVec)
#B0muSig <- B1muSig*mu1 -0.5
#phiVec <- phi(xiVec)
#}
##Need to define wtVec1 = diag(b'(sigmax+mu))
mu2<-mu1+0.5*sigma1
Thresh<-2
wtP<-mu2
wtP[mu2>=Thresh]<-exp(Thresh)*(.5*mu2[mu2>=Thresh]^2 + (1-Thresh)*mu2[mu2>=Thresh] + (1-Thresh*(.5*Thresh -1)))
wtP[mu2<Thresh]<-exp(mu2[mu2<Thresh])

#wtP<-exp(mu1+0.5*sigma1)
wtVec1<-as.vector(E1*mu1 + E2*wtP + E3*B0muSig)
##Need to define wtVec2 = diag(b''(sigmax+mu))
wtVec2<-as.vector(E1 + E2*wtP + E3*B1muSig)

##Update xiVec
#if(method=="JJ"){
#xiSqd <- diag(X%*%(Sigma.q.betauG + tcrossprod(mu.q.betauG))%*%t(X))

#for (i in 1:m) {
#           EsqMatCurr <- (-Sigma.q.betauG)%*%G[[i]]%*%H[[i]]
#           + tcrossprod(mu.q.betauG,mu.q.uR[[i]])
#           xiSqd[reBlockInds[[i]]] <- (xiSqd[reBlockInds[[i]]]
#                                       + 2*diag(X[reBlockInds[[i]],]%*%EsqMatCurr%*%t(Z[reBlockInds[[i]],])))
#           EsqMatCurr <- Sigma.q.uR[[i]] + tcrossprod(mu.q.uR[[i]])
#           xiSqd[reBlockInds[[i]]] <- (xiSqd[reBlockInds[[i]]]
#                                       + diag(Z[reBlockInds[[i]],]%*%EsqMatCurr%*%t(Z[reBlockInds[[i]],])))
#         }
#
#xiVec <- sqrt(xiSqd)
#}
               ## Update q*(sigsq.eps) parameter.
residSS<-trTermTwo<-trTermThree<-vector(length=RC)
BiasY<-as.vector((yAdj-mu1)^2)
ZtZSigmaqu <- vector("list", length=m)
for(i in 1:m){ZtZSigmaqu[[i]]<-diag((ZtZ[[i]]%*%Sigma.q.uR[[i]]))}

for(r in 1:RC){
         residSS[r] <- sum(BiasY[RBlockInds[[r]]])  
         trTermTwotmp <- sum(diag(crossprod(X[RBlockInds[[r]],CumPr[r]:(CumPr[(r+1)]-1)])%*%Sigma.q.betauG[CumPr[r]:(CumPr[(r+1)]-1),CumPr[r]:(CumPr[(r+1)]-1)])) 
   for (i in 1:m){trTermTwotmp <- trTermTwotmp + sum(ZtZSigmaqu[[i]][CumQr[r]:(CumQr[(r+1)]-1)])}
         trTermTwo[r]<-trTermTwotmp  
         trTermThreetmp <- 0  
         for (i in 1:m){
           trTermThreetmp <- trTermThreetmp + sum(diag(ZtX[[i]]%*%Sigma.q.betauG%*%GH[[i]])[CumQr[r]:(CumQr[(r+1)]-1)])
           }
       trTermThree[r]<-trTermThreetmp
   }  ##Closes the R loop
 }   ##Closes the streamlined group 

       ## Update q*(1/sigsq.eps) parameters.
if(!doStreamlined){
BiasY<-as.vector((yAdj - C%*%mu.q.betau)^2)
CTCSigma<-CTC%*%Sigma.q.betau
}

       if (!doStreamlined){
         for(r in 1:RC){
         rblock<-c(CumPr[r]:(CumPr[(r+1)]-1),P+rep(seq(0,(ncZ-1),by=Q),each=qr[r])+CumQr[r]:(CumQr[(r+1)]-1))
         B.q.sigsq.epsr[r] <- mu.q.recip.a.epsr[r] + 0.5*(sum(BiasY[RBlockInds[[r]]]) + sum(diag(CTCSigma)[rblock]))
                       }
                          }  

      if (doStreamlined){B.q.sigsq.epsr <- mu.q.recip.a.epsr + 0.5*(residSS + trTermTwo - 2*trTermThree)}

      mu.q.recip.sigsq.epsr <- A.q.sigsq.epsr/B.q.sigsq.epsr
       
       ## Update q*(a.eps) parameter:

       B.q.a.epsr <- mu.q.recip.sigsq.epsr + (1/A.epsr^2)
       mu.q.recip.a.epsr <- 1/B.q.a.epsr

     ## Update q*(a.k) parameters.
     B.q.a.k <- nuVal*diag(M.q.inv.SigmaR) + (1/A.k^2)
     mu.q.recip.ak <- (0.5*(nuVal + Q))/B.q.a.k

     ## Update q*(SigmaR^{-1}) parameters.
     if (Q == 1) B.q.SigmaR <- 2*nuVal*mu.q.recip.ak
     if (Q > 1) B.q.SigmaR <- 2*nuVal*diag(mu.q.recip.ak)
                indsStt <- P + 1
     for (i in 1:m) {
       if (!doStreamlined) {
         indsEnd <- indsStt + Q - 1 ; inds <- indsStt:indsEnd
         B.q.SigmaR <- (B.q.SigmaR
                        + Sigma.q.betau[inds,inds]
                        + tcrossprod(mu.q.betau[inds]))
         indsStt <- indsStt + Q
       }

       if (doStreamlined)
         B.q.SigmaR <- B.q.SigmaR + Sigma.q.uR[[i]] + tcrossprod(mu.q.uR[[i]])
     }
     M.q.inv.SigmaR <- (nuVal + m + Q - 1)*solve(B.q.SigmaR)
    
     ## Obtain current log(ML).
     logMLcurr <- logML(nuVal,R,P,ncZ,Q,m,numObs,mu1,E2,wtP,E3,B0muSig,B1muSig,xiVec,phiVec,
                        sigsq.beta,mu.q.betauG,Sigma.q.betauG,det.Sigma.q.betau,
                        B.q.SigmaR,A.k,M.q.inv.SigmaR,
                        mu.q.recip.ak,B.q.a.k,A.epsr,B.q.sigsq.epsr,B.q.a.epsr,
                        mu.q.recip.a.epsr,mu.q.recip.sigsq.epsr,mu.q.recip.Sigma.eps)
     
     logMLgrid[itnum] <- logMLcurr[[length(logMLcurr)]][1]
          
     ## Compute relative error.
     relErr <- abs((logMLcurr[[length(logMLcurr)]][1]/logMLprev)-1)
     
     ## Check `converged' conditions.
     if (itnum >= maxIter) {
       converged <- TRUE
       print("WARNING: maximum number of iterations exceeded.")
    }
     if (relErr < tolerance) converged <- TRUE
     cat(itnum,sep="\n")
     print(cbind(itnum,(logMLcurr[[length(logMLcurr)]][1]-logMLprev)))
        }

EndTime<-Sys.time()
VBTime<-EndTime-StartTime

   if (!doStreamlined) {
       return(list(mu.q.betauG=mu.q.betau,
                   Sigma.q.betauG=Sigma.q.betau,
                   A.q.sigsq.epsr=A.q.sigsq.epsr,B.q.sigsq.epsr=B.q.sigsq.epsr,
                   B.q.SigmaR=B.q.SigmaR,logMLgrid=logMLgrid,id=unique(id1),VBTime=VBTime))
   }

   if (doStreamlined) {
       return(list(mu.q.betauG=mu.q.betauG,Sigma.q.betauG=Sigma.q.betauG,
                   A.q.sigsq.epsr=A.q.sigsq.epsr,B.q.sigsq.epsr=B.q.sigsq.epsr,
                   B.q.SigmaR=B.q.SigmaR,logMLgrid=logMLgrid,
                   mu.q.recip.sigsq.epsr=mu.q.recip.sigsq.epsr,
                   mu.q.recip.a.epsr=mu.q.recip.a.epsr,
                   M.q.inv.SigmaR=M.q.inv.SigmaR,
                   mu.q.recip.ak=mu.q.recip.ak,
                   mu.q.recip.Sigma.eps=mu.q.recip.Sigma.eps,
                   Sigma.q.uR=Sigma.q.uR,
                   mu.q.uR=mu.q.uR,
                   id=unique(id1),
                   VBTime=VBTime))
   }
}
