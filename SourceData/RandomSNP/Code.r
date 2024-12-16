 SNP<-read.table("SNP.txt", sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE,row.names=1)
 set.seed(123)
 rownames(SNP)<-rownames(SNP)[sample(1:nrow(SNP),nrow(SNP))]
 write.table(SNP,file="SNPrandom.txt",sep="\t",quote=F,row.name=T)
 
 PhosphositesExp<-read.table("PhosphositesExp.txt", sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE,row.names=1)
 ProteinExp<-read.table("ProteinExp.txt", sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE,row.names=1)
 SNP<-read.table("SNPrandom.txt", sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE,row.names=1)

 intersample<-intersect(intersect( names(PhosphositesExp), names(ProteinExp)),names(SNP))
 ProteinExp<-ProteinExp[,match(intersample,names(ProteinExp))]
 PhosphositesExp<-PhosphositesExp[,match(intersample,names(PhosphositesExp))]
 PhosphositesExp<-PhosphositesExp[apply(PhosphositesExp,1,function(x){length(which(is.na(x)))})<=71,]
 SNP<-SNP[,match(intersample,names(SNP))]

 CisPair<-read.table("CisPair.txt", sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE)
 PhosphositesID<-read.table("PhosphositesID.txt", sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE)

 phosSNPs <- read.table("phosSNPs.txt", sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE)[,c(1,2)]
 phosSNPs<-phosSNPs[phosSNPs[,2]%in%rownames(SNP),]
 phosSNPs<-phosSNPs[paste(phosSNPs[,1],phosSNPs[,2])%in%paste(CisPair[,1],CisPair[,2]),]

 IVs<-unique(na.omit(CisPair))
 uniqueExposure<-intersect(unique(rownames(PhosphositesID)[PhosphositesID[,4]%in%IVs[,1]]),rownames(PhosphositesExp))#7932
 length(uniqueExposure)
 ClinicalfromTableS1A <- t(read.table("ClinicalfromTableS1A.csv", sep = ",", row.names=2,header = TRUE,fill=T, stringsAsFactors = FALSE)[,c(11,12,29)])
 rownames(ClinicalfromTableS1A)<-c('Gender','Age','Smoking')
 Covariates<-ClinicalfromTableS1A[,match(intersample,colnames(ClinicalfromTableS1A))]

 AGE<-as.numeric(Covariates[2,,drop=F])
 Gender<-as.numeric(Covariates[1,,drop=F])
 SMOKING<-as.numeric(Covariates[3,,drop=F])
 
 library(MendelianRandomization)

 Phoslink_Result<-c()
 FDR_MR_Result<-c()
 GWAS_MR_Result<-c()
 Min_MR_Result<-c()
 GWAS_genome_MR_Result<-c()

StartTime<-Sys.time()
print(paste('StartTime:',StartTime))
print(whichExposure)

 for(whichExposure in 1:nrow(PhosphositesExp)) {
 selectedExposure<-uniqueExposure[whichExposure]
 IV<- IVs[IVs[,1]%in%PhosphositesID[rownames(PhosphositesID)==selectedExposure,4],2]
 XX = as.numeric(PhosphositesExp[rownames(PhosphositesExp)==selectedExposure,])

 for(whichOutcome in 1:nrow(ProteinExp)){ 
 if(whichOutcome%%1000==0){print(whichOutcome);print(Sys.time())}
 selectedOutcome<-rownames(ProteinExp)[whichOutcome]
 X=XX
 Y = as.numeric(ProteinExp[rownames(ProteinExp)==selectedOutcome,]) 
 Xisna<-which(is.na(XX))
 Yisna<-which(is.na(Y))
 isna<-unique(c(Xisna,Yisna))
 Z<-t(SNP[rownames(SNP)%in%IV,,drop=F])#行是样本，列是基因
 
XZ<-Z
if(length(Xisna)>0){
XZ<-Z[-Xisna,,drop=F]
}
YZ<-Z
if(length(Yisna)>0){
YZ<-Z[-Yisna,,drop=F]
}
Z<-Z[,intersect(which(apply(XZ,2,function(x){length(which(x==0))<0.9*nrow(XZ)})),which(apply(YZ,2,function(x){length(which(x==0))<0.9*nrow(YZ)}))),drop=F]

 age=AGE
 gender=Gender
 smoking=SMOKING

 if(length(isna)<=71&ncol(Z)!=0){
 J = ncol(Z)
 N = nrow(Z)
 betaX <- array(NA, dim=J)
 betaY <- array(NA, dim=J)
 sebetaY <- array(NA, dim=J)
 sebetaX <- array(NA, dim=J)
 PrX <- array(NA, dim=J)
 PrY <- array(NA, dim=J)
 R2X<-array(NA, dim=J)
 R2Y<-array(NA, dim=J)
 R2P<-array(NA, dim=J)
 for(iIV in 1:J){
 if(length(Xisna)>0){regX <- lm(X[-Xisna] ~ Z[-Xisna,iIV]+AGE[-Xisna]+Gender[-Xisna]+SMOKING[-Xisna])}
 else{regX <- lm(X ~ Z[,iIV]+AGE+Gender+SMOKING)}
 if(length(Yisna)>0){regY <- lm(Y[-Yisna] ~ Z[-Yisna,iIV]+AGE[-Yisna]+Gender[-Yisna]+SMOKING[-Yisna])}
 else{regY <- lm(Y ~ Z[,iIV]+AGE+Gender+SMOKING)}

 betaX[iIV] <- summary(regX)$coefficients[2,1]
 sebetaX[iIV] <- summary(regX)$coefficients[2,2]
 PrX[iIV] <- summary(regX)$coefficients[2,4]
 betaY[iIV] <- summary(regY)$coefficients[2,1]
 sebetaY[iIV] <- summary(regY)$coefficients[2,2]
 PrY[iIV] <- summary(regY)$coefficients[2,4]
 }
 
 Phoslink_IV<-which(PrX<0.05&colnames(Z)%in%phosSNPs[phosSNPs[,1]==PhosphositesID[rownames(PhosphositesID)==selectedExposure,4],2])
 FDR_IV<- which(p.adjust(PrX,'fdr')<0.05)
 GWAS_IV<- which(p.adjust(PrX,'bonferroni')<0.05)
 GWAS_genome_IV<- which(PrX<5e-8)
 Min_IV<-c()
 if(any(p.adjust(PrX,'fdr')<0.05)){ Min_IV<-which.min(PrX)}
 
 #Phoslink 
 if(length(Phoslink_IV)>=1){
 MR_single<-as.data.frame(matrix(,1,13))
 colnames(MR_single)<-c("IV","Exposure","Outcome",
 "ivw_Estimate","ivw_Pvalue","ivw_CILower","ivw_CIUpper",'Spearman_E','Spearman_P','Pearson_E','Pearson_P','TSLS_E','TSLS_P')
 oggetto = mr_input(bx = as.numeric(betaX[Phoslink_IV]),
 bxse = as.numeric(sebetaX[Phoslink_IV]),
 by = as.numeric(betaY[Phoslink_IV]),
 byse = as.numeric(sebetaY[Phoslink_IV]),
 exposure = "PhosphositesExp", outcome = "ProteinExp",
 snps = colnames(Z)[Phoslink_IV]) 
 MRresult_ivw<-mr_ivw(oggetto,model = "default",weights="delta",psi=cor.test(X,Y)$estimate)
 MR_single[1,1:7]<- c(paste(colnames(Z)[Phoslink_IV],collapse="; "),selectedExposure,selectedOutcome,
          MRresult_ivw$Estimate,MRresult_ivw$Pvalue,MRresult_ivw$CILower,MRresult_ivw$CIUpper) 
 MR_single[1,8:11]<-c(cor.test(X, Y, method = "spearman")$estimate,cor.test(X, Y, method = "spearman")$p.value,cor.test(X, Y, method = "pearson")$estimate,cor.test(X, Y, method = "pearson")$p.value)
 Phoslink_Result<-rbind(Phoslink_Result,MR_single)
				}						
		}
	}
}

EndTime<-Sys.time()
print(paste('EndTime:',EndTime))
write.table(Phoslink_Result,"Phosphosite2Protein_Phoslink.txt",row.names=F,col.names=F,sep="\t",quote=F)

 Type<-'Phoslink' #Phoslink_Result FDR_MR_Result GWAS_genome_MR_Result GWAS_MR_Result
 Time=1
 MRresult<-read.table("Phosphosite2Protein_Phoslink.txt", sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE)

 LDres<-read.table("TotalLD_DQ.tsv", sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE)

 PhosphositesExp<-read.table("PhosphositesExp.txt", sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE,row.names=1)
 ProteinExp<-read.table("ProteinExp.txt", sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE,row.names=1)
 SNP<-read.table("SNPrandom.txt", sep = "\t", header = T,fill=T,quote="", stringsAsFactors = FALSE,row.names=1)

 intersample<-intersect(intersect( names(PhosphositesExp), names(ProteinExp)),names(SNP))
 ProteinExp<-ProteinExp[,match(intersample,names(ProteinExp))]
 PhosphositesExp<-PhosphositesExp[,match(intersample,names(PhosphositesExp))]
 PhosphositesExp<-PhosphositesExp[apply(PhosphositesExp,1,function(x){length(which(is.na(x)))})<=71,]
 SNP<-SNP[,match(intersample,names(SNP))]

 ClinicalfromTableS1A <- t(read.table("ClinicalfromTableS1A.csv", sep = ",", row.names=2,header = TRUE,fill=T, stringsAsFactors = FALSE)[,c(11,12,29)])
 rownames(ClinicalfromTableS1A)<-c('Gender','Age','Smoking')
 Covariates<-ClinicalfromTableS1A[,match(intersample,colnames(ClinicalfromTableS1A))]

 AGE<-as.numeric(Covariates[2,,drop=F])
 Gender<-as.numeric(Covariates[1,,drop=F])
 SMOKING<-as.numeric(Covariates[3,,drop=F])
 
 library(MendelianRandomization)

 Phoslink_Result<-c()
 FDR_MR_Result<-c()
 GWAS_MR_Result<-c()
 Min_MR_Result<-c()
 GWAS_genome_MR_Result<-c()
 LDcutoff<-0.2
 
 for(i in 1:nrow(MRresult)){
 if(i%%1000==0){print(i)}
 
 selectedExposure<-MRresult[i,2]
 selectedOutcome<-MRresult[i,3]
 IV<- unlist(strsplit(MRresult[i,1],split='; '))
 XX = as.numeric(PhosphositesExp[rownames(PhosphositesExp)==selectedExposure,])

 X=XX
 Y = as.numeric(ProteinExp[rownames(ProteinExp)==selectedOutcome,]) 
 Xisna<-which(is.na(XX))
 Yisna<-which(is.na(Y))
 isna<-unique(c(Xisna,Yisna))
 Z<-t(SNP[rownames(SNP)%in%IV,,drop=F])#行是样本，列是基因

 age=AGE
 gender=Gender
 smoking=SMOKING

 if(length(isna)<=71&ncol(Z)!=0){
 J = ncol(Z)#候选IV数量
 N = nrow(Z)#样本数量
 betaX <- array(NA, dim=J)
 betaY <- array(NA, dim=J)
 sebetaY <- array(NA, dim=J)
 sebetaX <- array(NA, dim=J)
 PrX <- array(NA, dim=J)
 PrY <- array(NA, dim=J)
 R2X<-array(NA, dim=J)
 R2Y<-array(NA, dim=J)
 R2P<-array(NA, dim=J)
 for(iIV in 1:J){
 if(length(Xisna)>0){regX <- lm(X[-Xisna] ~ Z[-Xisna,iIV]+AGE[-Xisna]+Gender[-Xisna]+SMOKING[-Xisna])}
 else{regX <- lm(X ~ Z[,iIV]+AGE+Gender+SMOKING)}
 if(length(Yisna)>0){regY <- lm(Y[-Yisna] ~ Z[-Yisna,iIV]+AGE[-Yisna]+Gender[-Yisna]+SMOKING[-Yisna])}
 else{regY <- lm(Y ~ Z[,iIV]+AGE+Gender+SMOKING)}
 library(psych)

 R2X[iIV] <-summary(regX)$r.squared
 R2Y[iIV] <-summary(regY)$r.squared
 R2P[iIV] <-r.test(length(X[-Xisna]),summary(regX)$r.squared,summary(regY)$r.squared)$p
 
 betaX[iIV] <- summary(regX)$coefficients[2,1]
 sebetaX[iIV] <- summary(regX)$coefficients[2,2]
 PrX[iIV] <- summary(regX)$coefficients[2,4]
 betaY[iIV] <- summary(regY)$coefficients[2,1]
 sebetaY[iIV] <- summary(regY)$coefficients[2,2]
 PrY[iIV] <- summary(regY)$coefficients[2,4]
 }
 
 LD<-LDres[LDres[,1]%in%colnames(Z)&LDres[,2]%in%colnames(Z),]
 DeleteLD<-c()
 while(any(LD[,3]>LDcutoff)){
 whichld=which(LD[,3]>LDcutoff)[1]
 DeleteLD<-c(DeleteLD,LD[whichld,which.max(c(PrX[colnames(Z)==LD[whichld,1]],PrX[colnames(Z)==LD[whichld,2]]))])
 LD<-LD[LD[,1]!=LD[whichld,which.max(c(PrX[colnames(Z)==LD[whichld,1]],PrX[colnames(Z)==LD[whichld,2]]))]&LD[,2]!=LD[whichld,which.max(c(PrX[colnames(Z)==LD[whichld,1]],PrX[colnames(Z)==LD[whichld,2]]))],]
 }
 Phoslink_IV<-(1:length(PrX))
 if(length(DeleteLD)!=0){Phoslink_IV<-(1:length(PrX))[-which(colnames(Z)%in%DeleteLD)]}
 FDR_IV<-(1:length(PrX))
 if(length(DeleteLD)!=0){FDR_IV<-(1:length(PrX))[-which(colnames(Z)%in%DeleteLD)]}
 GWAS_IV<-(1:length(PrX))
 if(length(DeleteLD)!=0){GWAS_IV<-(1:length(PrX))[-which(colnames(Z)%in%DeleteLD)]}
 GWAS_genome_IV<-(1:length(PrX))
 if(length(DeleteLD)!=0){GWAS_genome_IV<-(1:length(PrX))[-which(colnames(Z)%in%DeleteLD)]}

 #Phoslink 
 if(Type=='Phoslink'){
 if(length(Phoslink_IV)>=1){
 MR_single<-as.data.frame(matrix(,1,13))
 colnames(MR_single)<-c("IV","Exposure","Outcome",
 "ivw_Estimate","ivw_Pvalue","ivw_CILower","ivw_CIUpper",'Spearman_E','Spearman_P','Pearson_E','Pearson_P','TSLS_E','TSLS_P')
 oggetto = mr_input(bx = as.numeric(betaX[Phoslink_IV]),
 bxse = as.numeric(sebetaX[Phoslink_IV]),
 by = as.numeric(betaY[Phoslink_IV]),
 byse = as.numeric(sebetaY[Phoslink_IV]),
 exposure = "PhosphositesExp", outcome = "ProteinExp",
 snps = colnames(Z)[Phoslink_IV]) 
 MRresult_ivw<-mr_ivw(oggetto,model = "default",weights="delta",psi=cor.test(X,Y)$estimate)
 MR_single[1,1:7]<- c(paste(colnames(Z)[Phoslink_IV],collapse="; "),selectedExposure,selectedOutcome,
          MRresult_ivw$Estimate,MRresult_ivw$Pvalue,MRresult_ivw$CILower,MRresult_ivw$CIUpper) 
 MR_single[1,8:11]<-c(cor.test(X, Y, method = "spearman")$estimate,cor.test(X, Y, method = "spearman")$p.value,cor.test(X, Y, method = "pearson")$estimate,cor.test(X, Y, method = "pearson")$p.value)
 Phoslink_Result<-rbind(Phoslink_Result,MR_single)
				}	
	}
								
		}
	}
if(!is.null(dim(Phoslink_Result)) ){write.table(Phoslink_Result,"Phoslink_multiIV_withoutLD_LDcutoff0.2.tsv"),row.names=F,col.names=T,sep="\t",quote=F)}



