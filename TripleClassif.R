

# DESCRIPTION:
# Classification using logistic regression based on first five principal components of beta values from 64 CpGs, 40 DCA consensus variables and 12 clinical variables 

# USAGE:
# tripleClassif(Clin, MethylationD, model, PCAMatr, DCAMatr, PICpG, MeanSD)



# REQUIRED ARGUMENTS
# Clin - a data.frame with 12 columns, where column names are: "Age","Gender","Tumor_diam","Grade","TNM","Hb","Albumin","Calcium","ALP","GGT","TPC","Creatinine"
# MethylationD - normalized beta values from Illumina EPIC methylation arrays 
# model - a logistic regression model 
# PCAMatr - Matrix containing the loadings from the principal component analysis 
# DCAMatr - Matrix that indicates which CpG sites that are included in the 40 DCA consensus variables 
# PICpG - Vector containing names of CpG sites that previously have been identified as relevent in ccRCC 
# MeanSD - Data.frame with two columns: Mean_v and SD_v, containing mean values and standard deviations that are used for standardizing variables before principal component analysis
# DETAILS:
# NAs are not allowed in Clin or in MethylationD for CpGs included in PICpG

# VALUE:
# Returns a numeric vector with classification label (1- high risk for progression, 0- low risk for progression) 


tripleClassif <- function(Clin, MethylationD, model, PCAMatr, DCAMatr, PICpG, MeanSD)
{
 	
	if(anyNA(Clin))
		stop("Clinical data contains NA")

	if(!all.equal(colnames(Clin),c("Age","Gender","Tumor_diam","Grade","TNM","Hb","Albumin","Calcium","ALP","GGT","TPC","Creatinine")))
		stop("Wrong column names in Clinical data")

	DCAvar <- ConstrDCA(MethylationD,DCAMatr)		

	#Sort out data priviously identified CpGs with connection to ccRCC 
	MethP <- MethylationD[PICpG,]	

	if(anyNA(MethP))
		stop("Sites within PICPGS contains NA")				

	CombData <- CombD(DCAvar,MethP,Clin)
	PCA_var <- NewPCA(CombData,PCAMatr, MeanSD)
	Class_data <- LGClassify(as.data.frame(PCA_var), model, CutOff=0.1909037)

	

	return(Class_data)

}



# DESCRIPTION:

# Calculates DCA consensus variables 

# USAGE:
# ConstrDCA(MethD,DCADM) 


# REQUIRED ARGUMENTS
# MethD - Normalized beta values from Illumina EPIC methylation arrays 
# DCADM - Matrix with CpGs as rownames and DCA consensus variables as colnames. Ones indicate which CpGs that are included in the DCA variables

# DETAILS:
# NAs are ingored when calculating mean vaules for DCA consensus variables


ConstrDCA <- function(MethD,DCADM) 
{
	DCA_new <- matrix(NA,nrow=40,ncol=dim(MethD)[2])

	for(i in 1:40)
	{
		Ind_temp <- DCADM[,i]==1
		RN <- rownames(DCADM)[Ind_temp]
		CM <- colMeans(MethD[RN,], na.rm=TRUE)
		DCA_new[i,] <- CM
		
	}
	colnames(DCA_new) <- colnames(MethD)
	rownames(DCA_new) <- paste("DCA", 1:40, sep="")

	return(DCA_new)

}


# DESCRIPTION:

# Combines DCA variables, Methylation panel and clinical variables to a data frame with samples as columns and variables as rows

# USAGE:
# CombD(DCA_var,Meth_p,Clin_d)

# REQUIRED ARGUMENTS:
# DCA_var - Matrix with patients as columns and DCA consensus variables as rows.
# Meth_p - Matrix with beta values, where each column corresponds to one patient and each row to a CpG site
# Clin_d - Clinical data containing the following columns: "Age","Gender","Tumor_diam","Grade","TNM","Hb","Albumin","Calcium","ALP","GGT","TPC","Creatinine"


CombD <- function(DCA_var,Meth_p,Clin_d)
{
	Ind_temp <- colnames(Meth_p) %in% rownames(Clin_d)
	Meth_p2 <- Meth_p[,Ind_temp]
	MatchD <- match(colnames(Meth_p2),colnames(DCA_var))
	MatchD2 <- match(colnames(Meth_p2),rownames(Clin_d))
	DCA_match <- DCA_var[,MatchD]
	Clin_match <- Clin_d[MatchD2,]

	CombiData <- rbind(t(Clin_match),DCA_match,Meth_p2)
	
	return(CombiData)

}




# DESCRIPTION: 
# Creates the first 5 principal components based on input rotation matrix

# USAGE:
#NewPCA(CData, RotMat, MeanSD_df)


# REQUIRED ARGUMENTS:
# CData - Dataframe with combined data containing clinical variables, DCA consensus variables and beta values from priviously identified CpGs connected to ccRCC
# RotMat - Matrix with loadings for the first five principal components for variables in CData
# MeanSD_df -  Data.frame with two columns: Mean_v and SD_v, containing mean values and standard deviations that are used for standardizing variables before principal component analysis


NewPCA <- function(CData,RotMat, MeanSD_df)
{
	Transp_data <- as.matrix(t(CData))
	Ind <- match(colnames(Transp_data),rownames(MeanSD_df))
	M_SD <- MeanSD_df[Ind,]

	ScD <- sweep(Transp_data, MARGIN=2, STATS=M_SD$Mean_v,FUN="-",check.margin=TRUE)
	ScD <- sweep(ScD, MARGIN=2, STATS=M_SD$SD_v,FUN="/",check.margin=TRUE)


	PC_new <- ScD%*%RotMat


	return(PC_new)

}


# DESCRIPTION: 
# Predicts data based on input model

# USAGE:
# LGClassify(FivePC, LG_Mod, CutOff)


# REQUIRED ARGUMENTS:
# FivePC - Data frame with 5 columns that corresponds to the first five principal components 
# LG_mod - glm object containing the model to use for prediction
# CutOff - The cutoff for posterior probability at which the patients should be classified as high risk

LGClassify <- function(FivePC, LG_Mod, CutOff)
{
	
 	PredV <- predict(LG_Mod, FivePC,type="response")
	ClassPat <- 1*(PredV>CutOff)

	return(ClassPat)


}


