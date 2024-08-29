 MRanalysis<-function(whichIV){
 MR_single<-as.data.frame(matrix(,1,11))
 colnames(MR_single)<-c("IV","Exposure","Outcome",
 "ivw_Estimate","ivw_Pvalue","ivw_CILower","ivw_CIUpper",'Spearman_E','Spearman_P','Pearson_E','Pearson_P')
 oggetto = mr_input(bx = as.numeric(betaX[whichIV]),
 bxse = as.numeric(sebetaX[whichIV]),
 by = as.numeric(betaY[whichIV]),
 byse = as.numeric(sebetaY[whichIV]),
 exposure = "PhosphositesExp", outcome = "ProteinExp",
 snps = colnames(Z)[whichIV]) 
 MRresult_ivw<-mr_ivw(oggetto,model = "default",psi=cor.test(X,Y)$estimate)
 MR_single[1,1:7]<- c(paste(colnames(Z)[whichIV],collapse="; "),'selectedExposure','selectedOutcome',
          MRresult_ivw$Estimate,MRresult_ivw$Pvalue,MRresult_ivw$CILower,MRresult_ivw$CIUpper)
 MR_single[1,8:11]<-c(cor.test(X, Y, method = "spearman")$estimate,cor.test(X, Y, method = "spearman")$p.value,cor.test(X, Y, method = "pearson")$estimate,cor.test(X, Y, method = "pearson")$p.value)
 return(MR_single)
 }
 library(MendelianRandomization)
 library(EnvStats)
 Sample.N<-100000 
 theta = 0 # 0.6 
 thetaUx = 0.75
 thetaUy = 0.75
 invalidIVratio = 0 # 0.3 0.5
 IV.num<-30 
 SNP.num<-50
 Phoslink<-c()
 FDR_MR<-c()
 Min_MR<-c()
 GWAS_MR<-c()

 for(simulation in 1:2000){ 
 
 Gamma<-runif(IV.num, min=0.08, max=0.10)
 which.invalidIV<-sample(1:IV.num,IV.num*invalidIVratio)
 Alpha<-rep(0,IV.num)
 Alpha[which.invalidIV]<-rnorm(length(which.invalidIV),mean =0,sd=0.15)  
 Phi<-rep(0,IV.num)
 GG = sapply(1:IV.num, function(i){rbinom(Sample.N, 2, 0.3)})
 UU = GG %*% Phi + rnorm(Sample.N,mean = 0,sd=1)
 XX = GG %*% Gamma + thetaUx*UU+ rnorm(Sample.N,mean = 0,sd=1)
 
 GG_extra = sapply(1:(SNP.num-IV.num), function(i){rbinom(Sample.N, 2, 0.3)})
 ZZ<-cbind(GG,GG_extra)
  
 for(Sample.n in c(50,100,200,300)){
 selectSample=sample(1:Sample.N,Sample.n)
 U = UU[selectSample,,drop=F]
 G = GG[selectSample,] 
 X = XX[selectSample,,drop=F]
 Y = X %*% theta + G %*% Alpha + thetaUy*U+ rnorm(Sample.n,mean = 0,sd=1)
 Z<-cbind(G,GG_extra[selectSample,])
 Z<-Z[sample(1:Sample.n,replace = FALSE),]
  
 colnames(Z)<-paste0("SNP",1:ncol(Z))
 J = ncol(Z)
 N = nrow(Z)
 betaX <- array(NA, dim=J)
 betaY <- array(NA, dim=J)
 sebetaY <- array(NA, dim=J)
 sebetaX <- array(NA, dim=J)
 PrXX <- array(NA, dim=J)
 PrX <- array(NA, dim=J)
 PrY <- array(NA, dim=J)
 for(iIV in 1:J){
 if(sd(ZZ[,iIV])!=0&sd(Z[,iIV])!=0){
 regXX <- lm(XX ~ ZZ[,iIV])
 PrXX[iIV] <- summary(regXX)$coefficients[2,4]
 regX <- lm(X ~ Z[,iIV])
 betaX[iIV] <- summary(regX)$coefficients[2,1]
 sebetaX[iIV] <- summary(regX)$coefficients[2,2]
 PrX[iIV] <- summary(regX)$coefficients[2,4]
 regY <- lm(Y ~ Z[,iIV])
 betaY[iIV] <- summary(regY)$coefficients[2,1]
 sebetaY[iIV] <- summary(regY)$coefficients[2,2]
 PrY[iIV] <- summary(regY)$coefficients[2,4]
}
}
 #Phoslink
  whichIV<-intersect(which(p.adjust(PrXX,'bonferroni')<0.05),which(PrX<0.05));print(whichIV) 
 if(length(whichIV)!=0){
 MR_single<-MRanalysis(whichIV)
 Phoslink<-rbind(Phoslink,cbind(MR_single,Sample.n=Sample.n))
}
}
 simulation=simulation+1
 }
 write.table(Phoslink,paste0('Phoslink_theta',theta,'_invalidIVratio',invalidIVratio,'.txt'),row.names=F,col.names=T,sep="\t",quote=F)
