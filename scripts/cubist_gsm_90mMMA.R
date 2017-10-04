#-------------------------------------------------------------------------------------------------------------------------------------------#
###SOIL ATTRIBUTE MODELLING USING CUBIST & K-FOLD CROSS VALIDATION INCLUDING UNCERTAINTY CALCULATION (ESTIMATION OF UPPER & LOWER LIMITS)####
#-------------------------------------------------------------------------------------------------------------------------------------------#

##SUMMARY: 
# This code will fit a CUBIST model, perform k-fold cross validation, calculate the uncertainties, and map the prediction and lower and upper limits at each soil depth for multiple soil attributes. 
# This is performed at each soil depth where the validation metrics and mapping procedures are performed at each k fold and written to the output directory.
# The CUBIST coefficients, variable usage, partitions and model/validation diagnostics are also writtten.
# After the 10th K-fold the predictions as well as the lower and upper limits are averaged to arrive at the final outputs. Validation metrics are also amalgamated.
# The code is capable of modelling multiple soil attributes (i.e. batch processing), therefore several soil training tables can be placed in the training directory (refer to instructions below). 
# Much of the operations are performed using parallel processing techniques. Ensure packages "doParrallel" and "foreach" are installed.

##INSTRUCTIONS:
# In your designated working directory you will need to create two folders (ensure names are spelt correctly): "TrainingData" & "Covariates" (e.g: D:/DSM/TrainingData & D:/DSM/Covariates)
# The "TrainingData" folder will house the soil attribute training data and must have the following structure:  X, Y, ID, soildepth1, soildepth2, soildepth3, etc..

# For example:  X       Y       ID    x0to5cm      x5to15cm     x15to30cm
#             540512  5377883   1     44.38902517  50.79748922  NA
#             540712  5376883   2     18.38972517  27.39758926  74.51045857
#             531512  5375983   3     8.389025179  29.19448962  NA
#             535112  5371183   4     51.98907517  65.294489234 84.91145859

# Multiple soil attribute training data can be inserted to the TrainingData folder as seperate ".csv" files to represent different soil properties. 
# These are processed sequentially and can be considered as processing soil attributes in batch mode. 
# Ensure appropriate names represent each training ".csv" file, e.g. pH, EC, OC, Drainage etc...
# Ensure Training and Covariate data directory is set correctly on line 71,72 of this script!!!!!
# Ensure ensure that the workspace is correctly saved line 714 of this script!!!!!
# The "Covariate" folder will house all available covariates (e.g. terrain derivatives, radiometrics, satellite imagery etc...).
# Ensure all covariates are of the same extent and coordinate system. You can check this using SAGA GIS, freely availble from: http://sourceforge.net/projects/saga-gis/files/
# Place covariates in folder with appropriate names. All covartiate grids should be in SAGA native files with the extension ".sdat". 
# If you wish to use another covariate raster type/extension you will need to change the raster extentions on lines 104, 105 and 289 to your preferred format.
# Ensure that the list of factors for categorical data is set correctly L217, 235,240!!!!!

# You will need the following modules: raster,sp,gstat,Cubist,snow,cvTools,rgdal,doParrallel,foreach. Ensure they are installed and working.

##INTERPRETATION OF FINAL DIAGNOSTIC TABLE OUTPUT
# Calib_RMSE = RMSE of the cubist training model predictions versus the actual training data 
# Calib_R2 = R2 of the cubist training model predictions versus the actual training data 
# Calib_Bias = Mean of the cubist training model predictions minus the mean of the actual training data 
# Calib_CC = Concordance of the cubist training model predictions versus the actual training data 
# Valid_RMSE = RMSE of the cubist prediction on the independent validation dataset versus the actual validation data 
# Valid_R2 = R2 of the cubist prediction on the independent validation dataset versus the actual validation data 
# Valid_Bias = Mean of the cubist prediction on the independent validation dataset minus the mean of the actual validation data 
# ValidCC = Concordance of the cubist prediction on the independent validation dataset versus the actual validation data 
# Perc.within.UPP.LOW.limits = Proportion of validation data within the confines of the estimated upper and lower limits. Should be ~90% to indicate good upper/lower predictions.

####################################################################################################################################################################################################
###CODE BEGINS HERE#################################################################################################################################################################################
####################################################################################################################################################################################################
rm(list = ls()) 
#Bring in Modules
library(raster)
library(sp)
#library(epiR)
library(gstat)
library(Cubist)
library(snow)
library(cvTools)
library(rgdal)
library(caret)
library(doParallel)
library(foreach)


#library(bigmemory)
#library(biganalytics)
####################################################################################################################################################################################################
memory.limit()
memory.size(max=T) #Adjust if necessary

DataDirectory="~/dev/GSM.net/data/TrainingData_log" #########################################Needs to be adjusted accordingly!!!!!###################
CovariateDirectory="~/dev/GSM.net/data/Covariates_500m"#######################################Needs to be adjusted accordingly!!!!!###################

#DataDirectory="/projet/orleans/mtitia/migale_transfer/GSM_500m/data/TrainingData_log" #########################################Needs to be adjusted accordingly!!!!!###################
#CovariateDirectory="/projet/orleans/mtitia/migale_transfer/GSM_500m/data/Covariates_500m"##

RootDirectory<-gsub(gsub("^.*/","",DataDirectory),"",DataDirectory)
dir.create(paste(RootDirectory,"temp",sep=""))
TempDirectory=paste(RootDirectory,"temp",sep="")
dir.create(paste(RootDirectory,"tiles",sep=""))
TileDirectory=paste(RootDirectory,"tiles",sep="")
rasterOptions(tmpdir=TempDirectory)
dir.create(paste(RootDirectory,"Outputs",sep=""))
OutputDirectory=paste(RootDirectory,"Outputs",sep="")

#concordance Function
ccc<- function(observed, predicted){
  mx=mean(observed)
  my=mean(predicted)
  s2x=var(observed)
  s2y=var(predicted)
  sxy=mean((observed-mx)*(predicted-my))
  ccc=2*sxy/(s2x+s2y+(mx-my)^2 )
  return(ccc)}

setwd(DataDirectory) 

#Obtain list of soil data files
#DataFiles<- list.files(getwd(),  pattern=".csv", full.names=FALSE)
DataFiles<- list.files(getwd(),  pattern=".RData", full.names=FALSE)

#Iterate over all soil data files in Data Directory
#for (g in c(8,5,3,6,7,2,4,1,9)){
  g=8
  #Read in data from soil data file 
  setwd(DataDirectory)
  #data<- read.table(DataFiles[g], sep=";",header=TRUE) # all data
  #data<- read.csv(DataFiles[g], sep=";",header=TRUE) # all data
  load(DataFiles[g]) 
  data <- get(substr(DataFiles[g], 1,nchar(DataFiles[g])-6))
  data <- data[complete.cases(data$X),]
  
  #Define soil data name
  SoilVar<- gsub(substr(DataFiles[g], nchar(DataFiles[g])-5, nchar(DataFiles[g])),"_",DataFiles[g])
  
  #Import covariates as a covariate stack
  setwd(CovariateDirectory) 
  list.files(getwd(),  pattern="tif$", full.names=FALSE)
  files<- list.files(getwd(), pattern='tif$')
  r1<- raster(files[1])
  for(i in 2:length(files)){
    r1<- stack(r1,files[i])
    #browser()
  }
   NAvalue(r1) <- 65535
 
  # create a list with tiles of rasters
  extent_rast <- extent(r1)
  ystep <- round((extent_rast@ymax - extent_rast@ymin)/30)
  y_ext <- c(seq(extent_rast@ymin, extent_rast@ymax, ystep),extent_rast@ymax)
  
  r1_list <- list()
  
  for (tile in 1:(length(y_ext)-1)) {
    print(tile)
    new_ext <- extent(extent_rast@xmin, extent_rast@xmax, y_ext[tile],y_ext[tile+1] )
    r1_tile <- crop(r1, new_ext, filename=paste(TileDirectory,"/raster_tile_",tile,".tif",sep=""), type="GTiff", overwrite=T )
    r1_list[[tile]] <- r1_tile
    rm(r1_tile)
    gc()
    }
  r1_list <- r1_list[-c(1,2)] # remove tiles without data, so in this case specificly 2
  
  save(r1_list, file=paste(OutputDirectory,"/raster_tiles.RData", sep=""))
  
  #Extract covariate values to the training data
  dat.df<-data.frame(data)
  dat.sdf<-SpatialPointsDataFrame(coords=dat.df[,1:2], data=dat.df)
  dat.sdf.ext<-extract(r1,dat.sdf)
  datacov<-cbind(data,dat.sdf.ext)
  
  
  if(file.exists(paste((OutputDirectory),"/",SoilVar,sep=""))){
    OutputDirectorySoilVar=paste((OutputDirectory),"/",SoilVar,sep="")
  }else{
    dir.create(paste((OutputDirectory),"/",SoilVar,sep=""))
    OutputDirectorySoilVar=paste((OutputDirectory),"/",SoilVar,sep="")
  }
  
  ##Create K-fold datasets to assess model uncertainty
  ## MMA : you don't need cvFolds for doing this? Or use folds$subsets[folds$which == i,1]
  ## to reach a specific i-fold
  k<-10
  folds<- cvFolds(n=nrow(data), K=k, R = 1, type = "random") 
  ksize<- trunc(nrow(data)/10)
  kmat<- matrix(NA,nrow=ksize,ncol=10)
  for (i in 0:(k-1)){
    st= (ksize*i)+1
    en= (ksize*i)+ksize
    kmat[,i+1]<- folds$subsets[st:en]
  }
  if(ksize != nrow(data)/10){
  #} else{
    kmat_ext<- matrix(folds$subsets[(en+1):(nrow(data))])
    kmat2<- matrix(NA,nrow=1,ncol=10)
    for (i in 1:nrow(kmat_ext)){
      kmat2[1,i]<-kmat_ext[i,1]
    }
    kmat<- rbind(kmat,kmat2)
  }
  
  #Create diagnostic matrix (to be written)
  diagnostic_mat<-matrix(NA,nrow=((ncol(data)-3)*k),ncol=11)
  colnames(diagnostic_mat)<-c("Kfold_level","VarName","Calib_RMSE","Calib_R2","Calib_Bias","Calib_CC","Valid_RMSE","Valid_R2","Valid_Bias","Valid_CC","Perc within UPP/LOW limits")
  cnt<-0

## iterates over the CV-Folds
#for(q in 1:10){
for(q in 1){
    
    valdat<- kmat[,q]
    valdat<-na.omit(valdat)
    
    tradat_mat<- kmat[,-q]
    tradat<- c(tradat_mat[,1],tradat_mat[,2],tradat_mat[,3],tradat_mat[,4],tradat_mat[,5],tradat_mat[,6],tradat_mat[,7],tradat_mat[,8],tradat_mat[,9])
    tradat<-na.omit(tradat)
    
    #Iterate over each soil depth in data file
    
#    for(k in 1:(ncol(data)-3)){
    for(k in 1:1){
      #k <- 1
         #for(k in 4:6){
      cnt<-cnt+1
      diagnostic_mat[cnt,1]<-q
      #Define soil depth name
      #Tvar<- names <- (data)[k+3]
      Tvar<- names((data)[k+3])
      diagnostic_mat[cnt,2]<-Tvar
      
      #Define soil data depth training table
      soil.dat <- datacov[,c(1:3,k+3,(ncol(data)+1):ncol(datacov))]
      soil.dat <- na.omit(soil.dat)
      soil.dat <- soil.dat[which(soil.dat[,4]>(-9999) | soil.dat[,4]<(9999)), ]   
      
      #Retrieves the training data/validation data
      mod.dat<- soil.dat[tradat,] #calibration set
      mod.dat <- na.omit(mod.dat)
      val.dat<- soil.dat[valdat,] #validation set
      val.dat <- na.omit(val.dat)
      
      for (kk in c(5,6,9,16,24,27)){ # use names instead
        mod.dat[,kk] <- as.factor(mod.dat[,kk])
        val.dat[,kk] <- as.factor(val.dat[,kk])
      }
      
      #Set target vaiable calls
      TvarC<-paste("mod.dat$",Tvar,sep="")
      target.C<- eval(parse(text=TvarC))
            
      TvarV<-paste("val.dat$",Tvar,sep="")
      target.V<-eval(parse(text=TvarV))
      
      
      covstackpnts_list <- list()
      uplo_mat_list <- list()
      
      # create a factor list for predictions on raster
      factor_list <- list()
      elements <- c(5,6,9,16,24,27)
      for (kk in 1:length(elements)){
      fact <- elements[kk]
      factor_list[[kk]] <- levels(mod.dat[,fact])
      }
      names(factor_list) <- names(mod.dat[,c(5,6,9,16,24,27)])
      
     
      # create a new raster to predict every training model on : raster_pred
      raster_na <- function(rast_stack) {
            # set the factor values to NA for those which are not in the cubist model
            for (vars in 1:length(factor_list)){
              message("replace", names(factor_list[vars]))
                  rast_temp <- rast_stack[[names(factor_list[vars])]]
                  rast_fact <- sort(unique(getValues(rast_temp)))
                  df_fact <- as.numeric(factor_list[[vars]])
                  
                  fact_replace <- matrix(NA, nrow=length(rast_fact),ncol=2)
                  fact_replace[,1] <- rast_fact
                  # for (variable in 1:length(df_fact)){
                    
                  for (variable in 1:length(fact_replace[,1])){  
                    if(rast_fact[variable] %in% df_fact){ fact_replace[variable,2] <-  rast_fact[variable]} else {fact_replace[variable,2] <- NA }
                  }
                              
                  rast_stack[[names(factor_list[vars])]] <- reclassify(rast_stack[[names(factor_list[vars])]], fact_replace, filename=paste(TempDirectory, "/",names(rast_temp), "_rast_pred.tif", sep=""), type="GTiff", overwrite=T)
                  rm(rast_temp)
                  rm(variable)
                  rm(df_fact)
                  rm(fact_replace)
                  gc()
                }
              
            # set the NA value of GRASS to NA
            
            for (vars in 1:dim(rast_stack)[3]){
              message("replace NA", names(rast_stack[[vars]]))
              NA_replace <- matrix(NA, nrow=1,ncol=2)
              NA_replace[1,] <- c(65535, NA)
              
              if(names(rast_stack[[vars]]) == "NPP_max") { 
              message("skip NPP_max")
              rast_stack[[vars]] <- rast_stack[[vars]]
              } else {
              rast_stack[[vars]] <- reclassify(rast_stack[[vars]], NA_replace, filename=paste(TempDirectory, "/",names(rast_stack[[vars]]), "_replace_NA.tif", sep=""), type="GTiff", overwrite=T)
              gc()             
              }
            }
          return(rast_stack)
      }
      
      rast_pred <- raster_na(r1)
      names(rast_pred) <- names(r1)
      
      #raster_na <- function(rast_stack) {
        # set the factor values to NA for those which are not in the cubist model
        #for (vars in 1:length(factor_list)){
        #  rast_stack[[names(factor_list[vars])]][!rast_stack[[names(factor_list[vars])]]%in% as.numeric(levels(mod.dat[,names(factor_list[vars])]))] <- NA
        #} 
        
        ## set the NA value of GRASS to NA
        #for (vars in 1:dim(rast_stack)[3]){
        #  rast_stack[[vars]][rast_stack[[vars]]== 65535] <- NA
        #}
        #return(rast_stack)
        #}
        #gc()
      
        #rast_pred <- raster_na(r1)
      #writeRaster(rast_pred, filename=paste(TempDirectory,"/raster_pred.grd",sep=""), type="raster", overwrite=T )
      
      #Remove all temporary files 
      #tempfiles<- list.files(TempDirectory, pattern='.tif')
      #for(u in 1: length(tempfiles)){
      #  delete<-paste(TempDirectory,"/",tempfiles[u],sep="")
      #  file.remove(delete)
      #}
      

      # create a list with tiles of rasters
      rast_pred_list <- list()
      for (tile in 1:(length(y_ext)-1)) {
        print(tile)
        new_ext <- extent(extent_rast@xmin, extent_rast@xmax, y_ext[tile],y_ext[tile+1] )
        rast_pred_tile <- crop(rast_pred, new_ext,filename=paste(TempDirectory,"/rasterpred_tile_",tile,".tif",sep=""), type="GTiff", overwrite=T)
        rast_pred_list[[tile]] <- rast_pred_tile
        rm(rast_pred_tile)
        gc()
      }
            
      rast_pred_list <- rast_pred_list[-c(1,2)]    

      #Produce cubist model and predict onto training and validation datasets
      cl<-makeCluster(detectCores()-1)
      #cubistPred_drain<-cubist(x= mod.dat[,-c(1:4,20)], y=target.C,cubistControl(unbiased = F,extrapolation = 0, sample = 0,seed = sample.int(4096, size = 1) - 1L,label = "outcome"),committees = 1)     # fit cubist model
      cubistPred_drain<-cubist(x= mod.dat[,5:ncol(mod.dat)], y=target.C,cubistControl(unbiased = F,extrapolation = 0, sample = 0,seed = sample.int(4096, size = 1) - 1L,label = "outcome"),committees = 1)     # fit cubist model
      #grd_caret <- expand.grid(committees=c(20,10,5,1),neighbors=c(9,5,1,0))
      # cubistPred_drain <- train(x= mod.dat[,5:ncol(mod.dat)], y=target.C, method="cubist", tuneGrid=grd_caret)
      stopCluster(cl)
      setwd(OutputDirectorySoilVar)
      mod.pred<- predict(cubistPred_drain, newdata = mod.dat)
      mod.pred.V<- predict(cubistPred_drain, newdata = val.dat)
      #mod.pred<- predict(cubistPred_drain, newdata = mod.dat,neighbors = 0, )
      #mod.pred.V<- predict(cubistPred_drain, newdata = val.dat,neighbors = 0, )
      
      #DIOGNOSTICS (calibration)
      diagnostic_mat[cnt,3]<- sqrt(mean((target.C -mod.pred)^2))
      R2.cal<- lm(mod.pred ~target.C)
      diagnostic_mat[cnt,4]<-as.matrix(summary(R2.cal)$adj.r.squared)
      diagnostic_mat[cnt,5]<-mean(mod.pred)- mean(target.C)
      diagnostic_mat[cnt,6]<- ccc(target.C,mod.pred)
      
      #DIOGNOSTICS (validation)
      diagnostic_mat[cnt,7]<- sqrt(mean((target.V -mod.pred.V)^2))
      R2.val<- lm(mod.pred.V ~target.V)
      diagnostic_mat[cnt,8]<-as.matrix(summary(R2.val)$adj.r.squared)
      diagnostic_mat[cnt,9]<-mean(mod.pred.V)- mean(target.V)
      diagnostic_mat[cnt,10]<- ccc(target.V,mod.pred.V)
      
      setwd(OutputDirectory)
      diagnostic_mat2<-matrix(NA,nrow=1,ncol=11)
      colnames(diagnostic_mat2)<-c("Kfold_level","VarName","Calib_RMSE","Calib_R2","Calib_Bias","Calib_CC","Valid_RMSE","Valid_R2","Valid_Bias","Valid_CC","Perc within UPP/LOW limits")
      diagnostic_mat2[1,]<-diagnostic_mat[cnt,]
      #write.table(diagnostic_mat2,file=paste(SoilVar,Tvar,"_Kfold_",q,"_Diagnostics.csv",sep=""),sep=",", col.names=T,row.names=F, append=T)
      write.table(diagnostic_mat2,file=paste(SoilVar,Tvar,"_Diagnostics_0.csv",sep=""),sep=",", col.names=F,row.names=F, append=T)
      save(cubistPred_drain, file=paste(SoilVar,Tvar,"_Kfold_",q,"_.RData", sep=""))
      pred_data <- cbind(mod.dat[,1:3],target.C, mod.pred)
      names(pred_data) <- c("x","y","id", "target","predict")
      val_data <- cbind(val.dat[,1:3],target.V,mod.pred.V)
      names(val_data) <- c("x","y","id","target", "predict")
      save(pred_data, file=paste(q,"_",Tvar,"_", SoilVar,"pred_data.RData", sep=""))
      save(val_data, file=paste(q,"_",Tvar,"_", SoilVar,"val_data.RData", sep=""))
      #}}}

      
      #Add residuals to predictions
      val.dat$FP<- mod.pred.V 
      
      #sep up cubist rule iterations in matrix
      #rules<-cubistPred_drain$splits[2:5]
      if(is.null(cubistPred_drain$splits[2:6])){
        print("1 rule")} else{
      rules<-cubistPred_drain$splits[2:6] # including the rules on categorical variables
      #rules<-cubistPred_drain$finalModel$splits[1:6] # including the rules on categorical variables
      #rules<-cubistPred_drain$finalModel$splits[2:5]
      
      # remove quotes around the categorical data
      rules_prepa <- as.data.frame(rules[,c("rule","variable")])
      rules_prepa1 <- as.data.frame(sapply(rules_prepa, function(rmquote) gsub("\"", "", rmquote)))
        
      rules$variable <- as.character(rules_prepa1$variable) 
      rules <- transform(rules, rule=as.numeric(rule),variable=as.character(variable), dir=as.character(dir), category=as.character(category))
      }
        
      #If 'rules' result in NULL value indicates only one cubist rule established and the cubist model can be run for the entire dataset during LOCV
      ##otherwise each rule will be sequentially run in the 'else' argument 
      #if(is.null(rules)){ 
      if(is.null(cubistPred_drain$splits)){ 
        #Set up validation matrix to where the upper and lower limit attributes will be written
        val.dat2<-val.dat
        
        #RULE-BASED LOCV - this will estimate the Upper/Lower limits
        ########## Leave-one-out cross validation #####################################################################################################
        # VLM: run cubist model based on the final model variables, committees and neighbors 
        # cubistPred_drain$finalModel$tuneValue$neighbors, cubistPred_drain$finalModel$tuneValue$committees
        tbl=mod.dat
        looresfunc<- function(j){
          loocubistPred<-cubist(x= tbl[-j,5:(ncol(tbl)-1)], y=tbl[-j,4],cubistControl(unbiased = F,extrapolation = 0, sample = 0,seed = sample.int(4096, size = 1) - 1L,label = "outcome"),committees = 1)     # fit cubist model
          looPred<-predict( loocubistPred, newdata = tbl[j,],neighbors = 0)
          looRes<- tbl[j,4] - (looPred)#+int.resids1.r1) 
          return(looRes)
        }
        print("Single Rule. Running parallel processing function on LOCV loop...")
        cl<-makeCluster(detectCores()-1)
        registerDoParallel(cl)
        #strt<-Sys.time()
        looResiduals<-foreach(j=1:nrow(tbl), .packages=library(Cubist), .combine='rbind') %dopar%{looresfunc(j)}
        #print(Sys.time()-strt)
        stopCluster(cl)
        looResiduals<-as.numeric(looResiduals)
        #################################################################################################################################################
        
        r.ulPI<- quantile(looResiduals, probs = c(0.05,0.95), na.rm = FALSE,names = F, type = 7) # rule lower and upper PI
        
        #Use the upper/lower limits on the validation matrix to assess its accuracy
        
        val.dat2$FP_lowerPI<- r.ulPI[1]
        val.dat2$FP_upperPI<- r.ulPI[2]
        
        val.dat2$FP_lower<-val.dat2$FP + val.dat2$FP_lowerPI
        val.dat2$FP_upper<-val.dat2$FP + val.dat2$FP_upperPI
        
        #Attach the upper/lower values to the XY covariate feature space.
        for (tile in 1:length(r1_list)) {
        #for (tile in 1:3) {
          print(tile)
          
          beginCluster()
          xyrast_tile <- subset(r1_list[[tile]], subset=1)
          covstackpnts<- data.frame(rasterToPoints(xyrast_tile))
          covstackpnts<-na.omit(covstackpnts)
          rownames(covstackpnts) <- 1:nrow(covstackpnts)
          endCluster()
          
          covstackpnts_list[[tile]] <- covstackpnts
          
          
          #Set up upper/lower designation matrices
          uplo_mat<- matrix(NA,ncol=4,nrow=nrow(covstackpnts))
          uplo_mat[is.na(uplo_mat)] <- 0
          colnames(uplo_mat)<- c("X","Y","lower","upper")
          
          uplo_mat[,1]<-covstackpnts[,1]
          uplo_mat[,2]<-covstackpnts[,2]    
          uplo_mat[,3]<-r.ulPI[1]
          uplo_mat[,4]<-r.ulPI[2]
                    
          uplo_mat_list[[tile]] <- uplo_mat
          
          rm(covstackpnts)
          rm(uplo_mat)
          gc()
        }  
      
        #Determine the proportion of predicted values that reside within the upper/lower limits within the validation dataset and write to diagnostic matrix.
        TvarV<-paste("val.dat2$",Tvar,sep="")
        target.V<-eval(parse(text=TvarV))  
        target.V.FP_lower<-eval(parse(text="val.dat2$FP_lower"))
        target.V.FP_upper<-eval(parse(text="val.dat2$FP_upper"))
        val.dat2$PICP<- as.numeric(target.V >= target.V.FP_lower & target.V <= target.V.FP_upper)
        #summary(as.factor(val.dat2$PICP))
        diagnostic_mat[cnt,11]<-summary(as.factor(val.dat2$PICP))[2]/nrow(val.dat2)
      }else{ #Each rule to be segmented and run sequentially
        
        #Form the rules matrix
        rulesdf<- data.frame(rules) # this is the original
      
        #Set up validation matrix to where the upper and lower limit attributes will be written
        val.dat2<-val.dat[1,]
        val.dat2$rule<- -9999
        val.dat2$FP_lowerPI=-9999
        val.dat2$FP_upperPI=-9999
        val.dat2$FP_lower<--9999
        val.dat2$FP_upper<--9999
        
        #set up Upper/Lower value matrix
        numrules<-max(rules[,1])
        r.ulPI_mat<-matrix(NA,nrow=numrules,ncol=2)
        
        #Derive XY covariate values used to form the Cubist rules and write as matrix
        setwd(CovariateDirectory) 
        impcov<-paste(unique(rulesdf[,2]))
        files<-""
        for(h in 1:length(impcov)){
          files<-c(files,paste(impcov[h],".tif",sep=""))
        }
        files<-files[-1] 
        
        
        ############### here we build in a new loop to run over the different spatial point dataframes - which results
        # in the uplo_mat which are used in the mapping for the upper and lower boundaries or maybe build a function?
            
        for (tile in 1:length(r1_list)) {
        #for (tile in 1:3) {
        covstack_vars <- c()
            for (tt in 1:length(files)) {
            subset_var <-  substr(files[tt], 1,nchar(files[tt])-4)  
            covstack_vars <- c(covstack_vars,subset_var)
            }
        names(r1_list[[tile]]) <- names(r1)
        covstack <- subset(r1_list[[tile]], subset=covstack_vars)
        
        beginCluster()
        print(tile)
        covstackpnts<- data.frame(rasterToPoints(covstack))
        covstackpnts<-na.omit(covstackpnts)
        rownames(covstackpnts) <- 1:nrow(covstackpnts)
        endCluster()
        
        covstackpnts_list[[tile]] <- covstackpnts
          
         
          #Set up upper/lower designation matrices
          uplo_mat<- matrix(NA,ncol=4,nrow=1)
          uplo_mat[is.na(uplo_mat)] <- 0
          colnames(uplo_mat)<- c("X","Y","lower","upper")
          uplo_mat_list[[tile]] <- uplo_mat
        
          rm(covstackpnts)
          rm(uplo_mat)
          gc()
         }  
                 
        
        #Iterate over each cubist rule to derive upper/lower limits and 
        #designate the upper/lower limits to the covariate feature space for mapping purposes.
        
        for(i in 1:numrules){
          #for(i in 1:10){
          #for(i in 1:c(1,2,4,6)){
          #Derive a single cubist rule from the cubist rule matrix and manipulate it into a function
          cubrule<-rulesdf[rulesdf$rule==i,]
          cubrule<-cubrule[order(cubrule$variable),]
          cubrule<-data.frame(lapply(cubrule[,-1], as.character), stringsAsFactors=FALSE) # the colomn with the rule no. is deleted
          
          #Ensure same name variables are together (eg. DEM<220 & DEM>230)
          cubrule_spl<-split(cubrule,cubrule$variable) # the dataframe is split into a list containing each row as a component
          cubrule_uni<-unique(cubrule$variable)# the unique variable names as character string
          cubrule_matr<-matrix(NA,nrow=length(cubrule_uni),ncol=1)
          
          #################### rewrite lines
        for(j in 1:length(cubrule_uni)){
            # so lets start with a continuous variable j=2
            uvar<- paste("cubrule_spl$",cubrule_uni[j],sep="") # it goes back to the list with the decision rule
            uvar1<-eval(parse(text=uvar)) # it picks up the decision rule within the list
            
          if(nrow(uvar1)>1&&is.na(uvar1$value)==F){
            uvar1_mat<- matrix(NA,nrow=nrow(uvar1),ncol=1)
            for (f in 1:nrow(uvar1_mat)){
              uvar1_mat[f,1]<- paste("tbl$",uvar1[f,1],uvar1[f,2],uvar1[f,3], sep="")
            }
            uvar2 <- ""
            for (f in 1:nrow(uvar1_mat)){
              uvar2 <- paste(uvar2,"&",uvar1_mat[f,1])
            }
            uvar2 <- substr(uvar2, 4, nchar(uvar2))
            cubrule_matr[j,1] <- paste("(",uvar2,")",sep="")
          }else if (is.na(uvar1$value)==F){
            cubrule_matr[j,1] <- paste("tbl$",uvar1[1,1],uvar1[1,2],uvar1[1,3], sep="")
          }else{
            cubrule_matr[j,1] <- paste("tbl$",uvar1[1,1],"%in%c(",uvar1[1,4],")==T", sep="")
          } 
      }
          
          # so what if is.na(uvar1$value)=F, j = categorical values are in -- uvar1[1,4], 
          #cubrule_matr[j,1]<-paste("tbl$",uvar1[1,1],uvar1[1,2],uvar1[1,4], sep="") original
          # or make a vector of the values and do something like %in%c(,,)? possible query: 1%in%c(1,2,3)..paste("%in%c(",uvar1[1,4],")=T",sep="")
          
          #Make cubist rule into 1 line of query to form the function
          cubrule1<- ""
          for (j in 1:nrow(cubrule_matr)){
            cubrule1<- paste(cubrule1,"&",cubrule_matr[j,1]) 
          }
          cubrule1<-substr(cubrule1, 4, nchar(cubrule1)) # changed from 4
          #Produce the function
          CubistRule<- function(tbl){
            cr <- eval(parse(text=cubrule1))
            return(cr)
          }
          
          #RULE-BASED LOCV - this will estimate the Upper/Lower limits
          ########## Leave-one-out cross validation #####################################################################################################
          tbl=mod.dat
          tbl$rule<-CubistRule(tbl)
          tbl<-tbl[tbl$rule==TRUE,]
          looresfunc<- function(j){
            loocubistPred<-cubist(x= tbl[-j,5:(ncol(tbl)-1)], y=tbl[-j,4],cubistControl(unbiased = F,extrapolation = 0, sample = 0,seed = sample.int(4096, size = 1) - 1L,label = "outcome"),committees = 1)     # fit cubist model
            looPred<-predict(loocubistPred, newdata = tbl[j,],neighbors = 0)
            looRes<- tbl[j,4] - (looPred)#+int.resids1.r1) 
            return(looRes)
            gc()
            
          }
                    
      if (length(tbl[,1]) > 3){  
          print("Multiple Rules. Running parallel processing function on LOCV loop...")
          cl<-makeCluster(detectCores()-1)
          registerDoParallel(cl)
          #strt<-Sys.time()
          looResiduals<-foreach(j=1:nrow(tbl), .packages=library(Cubist), .combine='rbind') %dopar%{looresfunc(j)}
          #print(Sys.time()-strt)
          stopCluster(cl)
          looResiduals<-as.numeric(looResiduals)
      }else{ 
        err_rule <- c("error_rule", i)
        write.table(err_rule, file=paste(OutputDirectory, "/", SoilVar, Tvar,"_kfold_",q,  "_error_rules.txt",sep=""), append=T, col.names=F, row.names=F)
      }
          #################################################################################################################################################
          r.ulPI<- quantile(looResiduals, probs = c(0.05,0.95), na.rm = FALSE,names = F, type = 7) # rule lower and upper PI
          r.ulPI_mat[i,1]<-r.ulPI[1]
          r.ulPI_mat[i,2]<-r.ulPI[2]
          
          #Use the upper/lower limits on the validation matrix to assess its accuracy
          tbl=val.dat
          tbl$rule<-CubistRule(tbl)
          tbl<-tbl[tbl$rule==TRUE,]
          tbl$rule[which(tbl$rule==TRUE)]<-i
          
          if(nrow(tbl)==0){
            
            #Determine which XYs of the covariates are applicable to this particular cubist rule. Attach the XY's to the upper/lower designation matrices.
            # for each tile
          for (tile in 1:length(covstackpnts_list)){  
            tbl=covstackpnts_list[[tile]]
            tbl$rule<-CubistRule(tbl)
            tbl1<-tbl[tbl$rule==TRUE,]
            
            uplo_mat2<- matrix(NA,ncol=4,nrow=nrow(tbl1))
            uplo_mat2[,1]<-tbl1[,1]
            uplo_mat2[,2]<-tbl1[,2]    
            uplo_mat2[,3]<-r.ulPI[1]
            uplo_mat2[,4]<-r.ulPI[2]
            uplo_mat<-rbind(uplo_mat_list[[tile]],uplo_mat2)
            uplo_mat_list[[tile]] <- uplo_mat  # store the boundaries per tile
            
            covstackpnts_list[[tile]]<-tbl[tbl$rule==FALSE,]
            covstackpnts_list[[tile]]<-covstackpnts_list[[tile]][,-ncol(covstackpnts_list[[tile]])]
            }
          
          } else {
            tbl$FP_lowerPI=0
            tbl$FP_upperPI=0 
            
            tbl$FP_lowerPI[which(tbl$rule==i)]<- r.ulPI[1]
            tbl$FP_upperPI[which(tbl$rule==i)]<- r.ulPI[2]
            
            tbl$FP_lower<-tbl$FP + tbl$FP_lowerPI
            tbl$FP_upper<-tbl$FP + tbl$FP_upperPI
            
            val.dat2<-rbind(val.dat2,tbl)
            
            #Determine which XYs of the covariates are applicable to this particular cubist rule. Attach the XY's to the upper/lower designation matrices.
            # for each tile
            for (tile in 1:length(covstackpnts_list)){ 
            tbl=covstackpnts_list[[tile]]
            tbl$rule<-CubistRule(tbl)
            tbl1<-tbl[tbl$rule==TRUE,]
            
            uplo_mat2<- matrix(NA,ncol=4,nrow=nrow(tbl1))
            uplo_mat2[,1]<-tbl1[,1]
            uplo_mat2[,2]<-tbl1[,2]    
            uplo_mat2[,3]<-r.ulPI[1]
            uplo_mat2[,4]<-r.ulPI[2]
            uplo_mat<-rbind(uplo_mat_list[[tile]],uplo_mat2) 
            
            uplo_mat_list[[tile]] <- uplo_mat  # store the boundaries per tile
            
            covstackpnts_list[[tile]]<-tbl[tbl$rule==FALSE,]
            covstackpnts_list[[tile]]<-covstackpnts_list[[tile]][,-ncol(covstackpnts_list[[tile]])]
            gc()
          }   # run over tiles
          gc()
         
        } # closes else statement
     } # closes numrules
        #Determine the proportion of predicted values that reside within the upper/lower limits within the validation dataset and write to diagnostic matrix.
        val.dat2<-val.dat2[-1,] 
        TvarV<-paste("val.dat2$",Tvar,sep="")
        target.V<-eval(parse(text=TvarV))  
        target.V.FP_lower<-eval(parse(text="val.dat2$FP_lower"))
        target.V.FP_upper<-eval(parse(text="val.dat2$FP_upper"))
        val.dat2$PICP<- as.numeric(target.V >= target.V.FP_lower & target.V <= target.V.FP_upper)
        #summary(as.factor(val.dat2$PICP))
        diagnostic_mat[cnt,11]<-summary(as.factor(val.dat2$PICP))[2]/nrow(val.dat2)
        setwd(OutputDirectorySoilVar)
        gc()
        
      }
      
      ##Create sub folders and write out model results and diagnostics
      if(file.exists(paste((OutputDirectorySoilVar),"/",Tvar,sep=""))){
        OutputDirectory2=paste((OutputDirectorySoilVar),"/",Tvar,sep="")
      }else{
        dir.create(paste((OutputDirectorySoilVar),"/",Tvar,sep=""))
        OutputDirectory2=paste((OutputDirectorySoilVar),"/",Tvar,sep="")
      }
      
      setwd(OutputDirectory2)
      if(is.null(cubistPred_drain$splits)){}else{write.table(rulesdf,file=paste(SoilVar,Tvar,"_Kfold_",q,"_CubistPartitions.csv",sep=""),sep=",", col.names=T,row.names=F)}
      write.table(cubistPred_drain$coefficients,file=paste(SoilVar,Tvar,"_Kfold_",q,"_CubistCoefficients.csv",sep=""),sep=",", col.names=T,row.names=F)
      write.table(cubistPred_drain$usage,file=paste(SoilVar,Tvar,"_Kfold_",q,"_CubistVariableUsage.csv",sep=""),sep=",", col.names=T,row.names=F)
      diagnostic_mat2<-matrix(NA,nrow=1,ncol=11)
      colnames(diagnostic_mat2)<-c("Kfold_level","VarName","Calib_RMSE","Calib_R2","Calib_Bias","Calib_CC","Valid_RMSE","Valid_R2","Valid_Bias","Valid_CC","Perc within UPP/LOW limits")
      diagnostic_mat2[1,]<-diagnostic_mat[cnt,]
      write.table(diagnostic_mat2,file=paste(SoilVar,Tvar,"_Kfold_",q,"_Diagnostics.csv",sep=""),sep=",", col.names=T,row.names=F)
      
     gc()
      
     save.image("D:/Vlmulder/project_GSM/GSM_soc/ws_soc_90m.RData")
      #Cubist Mapping#################################################################################################
      
      #setwd(PredictedDirectory)
      setwd(TempDirectory)
      # here loop over the list with rasters
      
     for (tile in 1:length(rast_pred_list)){
       # for (tile in 1:3){
         
      names(rast_pred_list[[tile]]) <- names(r1)
      beginCluster()
      #rast_pred_list[[tile]]
      drain_cubistMap <- clusterR(rast_pred_list[[tile]], predict, args=list(cubistPred_drain,neighbors = 0),factors=factor_list, na.rm=T,inf.rm=T,
                                  filename= paste(SoilVar,Tvar,"_Kfold_",q,"_predicted_tile", tile,".tif",sep=""),format="GTiff",progress="text",overwrite=T)
      endCluster()
      gc()
      
      
      #Add lower and upper prediction intervas (where there is more than 1 rule) and map.
      #uplo_mat<-uplo_mat[-1,]
      uplo_mat<-uplo_mat_list[[tile]][-1,]
      lower_xy<- data.frame(uplo_mat[,1:3])
      upper_xy<- data.frame(uplo_mat[,1:2],uplo_mat[,4])
      
      lower_rast<- rasterFromXYZ(lower_xy, res=res(r1), crs=NA, digits=2)
      upper_rast<- rasterFromXYZ(upper_xy, res=res(r1), crs=NA, digits=2)
      
      
      
      # create the lower prediction boundaries per tile
      setwd(TempDirectory)
      
      r_low <- extend(lower_rast,drain_cubistMap)
      crs(r_low) <- CRS("+init=epsg:2154")
      
      
      r2<- stack(r_low,drain_cubistMap)
      crs(r2) <- CRS("+init=epsg:2154")
      
      f1 <- function(x) calc(x, sum, na.rm=T)
      beginCluster()
      lower_pred <- clusterR(r2, fun=f1, filename =paste("temp_low_tile",tile, ".tif", sep=""),format="GTiff",progress="text",overwrite=T)
      endCluster()
      gc()
      
      # create the lower prediction boundaries per tile
      r_up <- extend(upper_rast,drain_cubistMap)
      crs(r_up) <- CRS("+init=epsg:2154")
      
      r2<- stack(r_up,drain_cubistMap)
      crs(r2) <- CRS("+init=epsg:2154")
      
      f1 <- function(x) calc(x, sum, na.rm=T)
      beginCluster()
      upper_pred <- clusterR(r2, fun=f1, filename =paste("temp_up_tile",tile, ".tif", sep=""),format="GTiff",progress="text",overwrite=T)
      endCluster()
      gc()
      
      }
     
     #Merge the predicted cubist values from each tile
     load_raster <- function (x) {
       maps <- list()
       for (rast in 1:length(x)) {  
         maps[[rast]] <- raster(x[rast])
       }
       return(maps)
     }
     
      ##Place in separate subfolder
      if(file.exists(paste((OutputDirectory2),"/predicted",sep=""))){
        PredictedDirectory=paste((OutputDirectory2),"/predicted",sep="")
      }else{
        dir.create(paste((OutputDirectory2),"/predicted",sep=""))
        PredictedDirectory=paste((OutputDirectory2),"/predicted",sep="")
      }
      
     setwd(TempDirectory)
     
     pred_files<- list.files(getwd(), pattern="predicted_tile", full.names=F)
     predictions <- load_raster(pred_files)
     gc()
     pred_merge <- predictions[[1]]
     writeRaster(pred_merge, filename=paste(TempDirectory,"/temp_pred_merge.tif",sep=""), format="GTiff", overwrite=T)
     
     for(i in 2:length(predictions)){
       message("merge predict",i)
       pred_merge <- merge(pred_merge,predictions[[i]], filename=paste(TempDirectory,"/temp_pred_merge_",i, ".tif",sep=""), format="GTiff", overwrite=T)
       gc()
        } 
     gc()
     
     setwd(PredictedDirectory)
     writeRaster(pred_merge, filename= paste(SoilVar,Tvar,"_Kfold_",q,"_predicted.tif",sep=""),format="GTiff",progress="text",overwrite=T)
     
     
     ##Place lower predictions in separate subfolder
      if(file.exists(paste((OutputDirectory2),"/lower",sep=""))){
        LowerPredDirectory=paste((OutputDirectory2),"/lower",sep="")
      }else{
        dir.create(paste((OutputDirectory2),"/lower",sep=""))
        LowerPredDirectory=paste((OutputDirectory2),"/lower",sep="")
      }
     
     setwd(TempDirectory)
     
     low_files<- list.files(getwd(), pattern="temp_low_tile", full.names=F)
     lowers <- load_raster(low_files)
     
     low_merge <- lowers[[1]]
     gc()
     for(i in 2:length(lowers)){
       message("lower merge", i)
       low_merge <- merge(low_merge,lowers[[i]], filename=paste(TempDirectory, "/temp_low_merge_",i, ".tif",sep=""), format="GTiff", overwrite=T)
       gc()
     } 
     gc()
     setwd(LowerPredDirectory)
     writeRaster(low_merge, filename= paste(SoilVar,Tvar,"_Kfold_",q,"_lowerpred.tif",sep=""),format="GTiff",progress="text",overwrite=T)
     
      ##Place upper predictions in separate subfolder
      if(file.exists(paste((OutputDirectory2),"/upper",sep=""))){
        UpperPredDirectory=paste((OutputDirectory2),"/upper",sep="")
      }else{
        dir.create(paste((OutputDirectory2),"/upper",sep=""))
        UpperPredDirectory=paste((OutputDirectory2),"/upper",sep="")
      }
     setwd(TempDirectory)
     
     up_files<- list.files(getwd(), pattern="temp_up_tile", full.names=F)
     uppers <- load_raster(up_files)
     gc()
     up_merge <- uppers[[1]]
     
     for(i in 2:length(uppers)){
       message("upper merge", i)
       up_merge <- merge(up_merge,uppers[[i]],filename=paste(TempDirectory,"/temp_up_merge_",i,".tif",sep=""), format="GTiff", overwrite=T)
       gc()
     } 
     gc()
     
     setwd(UpperPredDirectory)
     writeRaster(up_merge, filename= paste(SoilVar,Tvar,"_Kfold_",q,"_upperpred.tif",sep=""),format="GTiff",progress="text",overwrite=T)
     
     #Remove temp raster files 
     rasfiles<- list.files(TempDirectory, pattern='.tif')
     for(u in 1: length(rasfiles)){
     delete<-paste(TempDirectory,"/",rasfiles[u],sep="")
     file.remove(delete)
     }
      
      
    } # closes k
  gc()  
  } # closes q
  
  #Write out final Soil Diagnostics
  setwd(OutputDirectorySoilVar)
  diagnostic_mat1<-data.frame(diagnostic_mat[order(diagnostic_mat[,2], diagnostic_mat[,1]),]) 
  diagnostic_mat2<- split(diagnostic_mat1,diagnostic_mat1[,2])
  
  diagnam<-names(diagnostic_mat2)
  for(p in diagnam){
    print(p)
    diagvar<-paste("diagnostic_mat2$",p,sep="")
    diagmat<-eval(parse(text=diagvar))
    
    diagmat1<-cbind(diagmat$Kfold_level,diagmat[,3:ncol(diagmat)])
    diagmat2<-matrix(NA,nrow=nrow(diagmat1),ncol=ncol(diagmat1))
    for(s in 1:ncol(diagmat1)){
      for(t in 1:nrow(diagmat1)){
        diagmat2[t,s]<-as.numeric(as.character(diagmat1[t,s]))
      }
    }
    diagmat2<-rbind(diagmat2,NA)
    diagmat2[(nrow(diagmat2)),2:ncol(diagmat2)]<-colMeans(diagmat2[1:(nrow(diagmat2)-1),2:ncol(diagmat2)])
    diagmat2<-diagmat2[order(diagmat2[,1]),]
    
    diagmat2<-rbind(diagmat2,NA)
    for(o in 2:ncol(diagmat2)){
      diagmat2[nrow(diagmat2),o]<-sd(diagmat2[1:(nrow(diagmat2)-2),o])
    }
    diagmat2<-cbind(c("K1","K2","K3","K4","K5","K6","K7","K8","K9","K10","MEAN","STDEV"),diagmat2[,2:ncol(diagmat2)])
    colnames(diagmat2)<- c(names(diagmat[1]),names(diagmat[3:length(names(diagmat))]))
    write.table(diagmat2,file=paste(SoilVar,p,"_Diagnostics.csv",sep=""),sep=",", col.names=T,row.names=F)
  }
  gc()
  #Average out K-fold predictions
  
  #predicted
  for(p in diagnam){
    OutputDirectory3=paste((OutputDirectorySoilVar),"/",p,"/predicted",sep="")
    setwd(OutputDirectory3)
    files<- list.files(getwd(), pattern='.tif$')
    lowerstack<- raster(files[1])
    for(o in 2:length(files)){
      lowerstack<- stack(lowerstack,files[o])}
    f1 <- function(x) calc(x, mean, na.rm=T)
    setwd(OutputDirectorySoilVar)
    beginCluster()
    upper_pred <- clusterR(lowerstack, fun=f1, filename =paste(SoilVar,p,"_predicted_mean.tif",sep=""),format="GTiff",progress="text",overwrite=T)
    endCluster()
    #DeleteDirectory=paste((OutputDirectorySoilVar),"/",p,sep="")
    #setwd(DeleteDirectory)
    #unlink("predicted", recursive = TRUE, force = TRUE)
  } 
  gc()
  #lower
  for(p in diagnam){
    OutputDirectory3=paste((OutputDirectorySoilVar),"/",p,"/lower",sep="")
    setwd(OutputDirectory3)
    files<- list.files(getwd(), pattern='.tif$')
    lowerstack<- raster(files[1])
    for(o in 2:length(files)){
      lowerstack<- stack(lowerstack,files[o])}
    f1 <- function(x) calc(x, mean, na.rm=T)
    setwd(OutputDirectorySoilVar)
    beginCluster()
    upper_pred <- clusterR(lowerstack, fun=f1, filename =paste(SoilVar,p,"_lowerpred_mean.tif",sep=""),format="GTiff",progress="text",overwrite=T)
    endCluster()
    #DeleteDirectory=paste((OutputDirectorySoilVar),"/",p,sep="")
    #setwd(DeleteDirectory)
    #unlink("lower", recursive = TRUE, force = TRUE)
  }  
  gc()
  #upper
  for(p in diagnam){
    OutputDirectory3=paste((OutputDirectorySoilVar),"/",p,"/upper",sep="")
    setwd(OutputDirectory3)
    files<- list.files(getwd(), pattern='.tif$')
    lowerstack<- raster(files[1])
    for(o in 2:length(files)){
      lowerstack<- stack(lowerstack,files[o])}
    f1 <- function(x) calc(x, mean, na.rm=T)
    setwd(OutputDirectorySoilVar)
    beginCluster()
    upper_pred <- clusterR(lowerstack, fun=f1, filename =paste(SoilVar,p,"_upperpred_mean.tif",sep=""),format="GTiff",progress="text",overwrite=T)
    endCluster()
    #DeleteDirectory=paste((OutputDirectorySoilVar),"/",p,sep="")
    #setwd(DeleteDirectory)
    #unlink("upper", recursive = TRUE, force = TRUE)
    #Remove temp raster files 
    #rasfiles<- list.files(TempDirectory, pattern='$')
    #for(u in 1: length(rasfiles)){
    #  delete<-paste(TempDirectory,"/",rasfiles[u],sep="")
    #  file.remove(delete)
    #}
  }  # closes upper
  #Remove all temporary files 
#tempfiles<- list.files(TempDirectory, pattern='$')
#for(u in 1: length(tempfiles)){
#delete<-paste(TempDirectory,"/",tempfiles[u],sep="")
# file.remove(delete)
  }
  gc()
} # closes g
setwd(RootDirectory)
unlink("temp", recursive = TRUE, force = TRUE)
