updateZ = function(Y,Z,Beta,iSigma,Eta,Lambda, X,Pi,dfPi,distr,rL, ind, EtaStar){
   ZPrev = Z
   ny = nrow(Y)
   ns = ncol(Y)
   nr = ncol(Pi)
   np = apply(Pi, 2, function(a) length(unique(a)))

   switch(class(X),
      matrix = {
         LFix = X%*%Beta
      },
      list = {
         LFix = matrix(NA,ny,ns)
         for(j in 1:ns)
            LFix[,j] = X[[j]]%*%Beta[,j]
      }
   )
   LRan = LRanStar = vector("list", nr)
   for(r in seq_len(nr)){
      if(rL[[r]]$xDim == 0){
         LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
         LRanStar[[r]] = EtaStar[[r]][Pi[,r],]%*%Lambda[[r]]
      } else{
         LRan[[r]] = LRanStar[[r]] = matrix(0,ny,ns)
         for(k in 1:rL[[r]]$xDim){
             LRan[[r]] = LRan[[r]] + (Eta[[r]][Pi[,r],]*rL[[r]]$x[as.character(dfPi[,r]),k]) %*% Lambda[[r]][,,k]
             LRanStar[[r]] = LRanStar[[r]] + (EtaStar[[r]][Pi[,r],]*rL[[r]]$x[as.character(dfPi[,r]),k]) %*% Lambda[[r]][,,k]

         }
      }
   }
   if(nr > 0){
      E = LFix + Reduce("+", LRan)
      EStar=LFix + Reduce("+", LRanStar)
   } else
      E = LFix

   Z = ZStar = matrix(NA,ny,ns)
   indNA = is.na(Y)
   std = matrix(iSigma^-0.5,ny,ns,byrow=TRUE)

   indColNormal = (distr[,1]==1)
   Z[,indColNormal] = Y[,indColNormal]
   ZStar[,indColNormal] = Y[,indColNormal]

   indColProbit = (distr[,1]==2)
   pN = sum(indColProbit)
   if(pN > 0){
      ZProbit = ZStarProbit = matrix(NA,ny,pN)
      YProbit = Y[,indColProbit]
      EProbit = E[,indColProbit]
      EStarProbit = EStar[,indColProbit]
      stdProbit = std[,indColProbit]
      indCellProbit = !indNA[,indColProbit]
      if(any(indCellProbit)){
         YProbit = as.logical(YProbit[indCellProbit])
         e = EProbit[indCellProbit]
         eSTar=EStarProbit[indCellProbit]
         s = stdProbit[indCellProbit]
         lB = rep(-Inf, length(YProbit))
         uB = rep(Inf, length(YProbit))
         lB[YProbit] = 0
         uB[!YProbit] = 0
         z = rtruncnorm(length(YProbit), a=lB, b=uB, mean=e, sd=s) # this is often the bottleneck for performance
         zStar = rtruncnorm(length(YProbit), a=lB, b=uB, mean=eStar, sd=s) # this is often the bottleneck for performance

         ZProbit[indCellProbit] = z
         ZProbitStar[indCellProbit] = zStar
         Z[,indColProbit] = ZProbit         
         ZStar[,indColProbit] = ZProbitStar

         
      }
   }

   indColPoisson = (distr[,1]==3)
   pN = sum(indColPoisson)
   if(pN > 0){
      r = 1e3# acquiring Poisson as limit of negative-binomial
      ZPoisson = matrix(NA,ny,pN)
      YPoisson = Y[,indColPoisson]
      EPoisson = E[,indColPoisson]
      stdPoisson = std[,indColPoisson]
      indCellPoisson = !indNA[,indColPoisson]
      if(any(indCellPoisson)){
         y = YPoisson[indCellPoisson]
         e = EPoisson[indCellPoisson]
         s = stdPoisson[indCellPoisson]
         zPrev = ZPrev[,indColPoisson][indCellPoisson]
         w = rpg(num=length(y), h=y+r, z=zPrev-1*log(r))
         prec = s^-2
         sigmaZ = (prec + w)^-1
         muZ = sigmaZ*((y-r)/(2) + prec*(e-log(r))) + 1*log(r)
         z = rnorm(length(y), muZ, sqrt(sigmaZ))
         if(any(is.na(z) | is.nan(z))){
            print("Fail in Poisson Z update")
         }
         ZPoisson[indCellPoisson] = z
         Z[,indColPoisson] = ZPoisson
      }
   }

   Z[indNA] = rnorm(sum(indNA), E[indNA], std[indNA])
   ZStar[indNA] = rnorm(sum(indNA), EStar[indNA], std[indNA])

   return(list(Z=Z, ZStar=ZStar)) ## need to make sure whatever calls this knows that I changed the output
}


