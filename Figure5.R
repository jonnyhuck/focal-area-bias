
##
# Produce the 'squirrel model' chart for Figure 8
##

#Load libraries
library(terra)
library(sf)
library(raster)

#load in the squirrel data
sciurus<- read.csv("data/sciurus.csv")


#subset the data to only include points with complete coordinates.
sciurus<-subset(sciurus,!is.na(sciurus$Latitude))

#remove all points with uncertainty > 100m
sciurus<-sciurus[sciurus$Coordinate.uncertainty_m<100,]


#make spatial points layer

#create crs object
sciurus.latlong<-data.frame(x=sciurus$Longitude,y=sciurus$Latitude)



#Use coordinates object to create spatial points object

sciurus.sp<-vect(sciurus.latlong,geom=c("x","y"))

#check that the points now have our desired crs. 
crs(sciurus.sp)<-"epsg:4326"


studyExtent<-c(-4.2,-2.7,56.5,57.5) #list coordinates in the order: min Long, max Long, min Lat, max Lat

#now crop points to this area

C<-crop(sciurus.sp,studyExtent)

#load land-cover map
LCM<-rast("C:/SE/LCM25.tif")

#project sciurus points to LCM
sciurusFin<-project(C,crs(LCM))

#get coordinates
sciurusCoords<-crds(sciurusFin)

#buffer to include analysis extent (using 5km buffers later)
x.min <- min(sciurusCoords[,1]) - 5000
x.max <- max(sciurusCoords[,1]) + 5000
y.min <- min(sciurusCoords[,2]) - 5000
y.max <- max(sciurusCoords[,2]) + 5000

#make object from extent
extent.new <- ext(x.min, x.max, y.min, y.max)

#crop LCM to this extent
LCM <- crop(LCM$LCMUK_1, extent.new)

# make reproducible
set.seed(11)

#generate background (pseudo absence) points

back.xy <- spatSample(LCM, size=1000,as.points=TRUE,ext=ext(sciurusFin)) 




#make data frame from background coordinates
Abs<-data.frame(crds(back.xy),Pres=0)

#make data frame from presence coordinates
Pres<-data.frame(crds(sciurusFin),Pres=1)
colnames(Pres)<-c("x","y","Pres")

library(dplyr)

# bind the two data frames by row
sciurusData<-rbind(Pres,Abs) 


#access levels of the raster by treating them as categorical data

LCM<-as.factor(LCM)

#create an vector object with new values for a binary woodland layer
reclass <- c(0,1,rep(0,19))

# combine with the LCM categories into a matrix of old and new values.
RCmatrix<- cbind(levels(LCM)[[1]],reclass)

RCmatrix<-RCmatrix[,2:3]# only need columns 2 and 3
RCmatrix$LCMUK_1<-as.numeric(RCmatrix$LCMUK_1)#convert character to numeric


#Use the reclassify() function to asssign new values to LCM with reclassification matrix
broadleaf <- classify(LCM, RCmatrix)


# set up vector of budder distances
radii<-seq(100,5000,by=100)

#this function iteratres of all points, buffers over all radii and extracts woodland cover
landBuffer <- function(x){         
  
  #use terra to convert each row to a spatial point using the columns "x" and "y" as the geometry.
  SciurusPoints<-vect(sciurusData[x,],geom=c("x","y"),crs="epsg:27700")
  # initiate list of woodland cover % values
  cover<-list()
  for(i in radii){  #for loop to iterate over the radius distances
    
    sciurusBuffer <- buffer(SciurusPoints, width=i)           #buffer each point
    bufferlandcover <- crop(broadleaf, sciurusBuffer)         #crop the landcover layer to the buffer extent
    
    masklandcover <- mask(bufferlandcover, sciurusBuffer)            # mask the above to the buffer
    landcoverArea <- sum(values(masklandcover),na.rm = TRUE)*625      #and sum landcover area
    
   # convert to precentage cover
    percentcover <- landcoverArea/(pi*i^2)*100
    cover[[i]]<-percentcover # store in list
  }
  print(x) #(to gauge progress)
  
  return(unlist(cover))} #get the list of percentages ("unlist" makes sure we get a vector each time so that we don't end up with a list of lists). 

# get row numbers to iterate over each sampling point
nPoints<-1:nrow(sciurusData)

#run the function
woodlandBuffers<-lapply(nPoints,FUN=landBuffer)

# bind all results together
glmDataBuff<-do.call("rbind",woodlandBuffers)

#convert to data frame
glmDataBuff<-data.frame(glmDataBuff)

# set up column names
colnames(glmDataBuff)<-c(paste("w", radii, sep = ""))

#add in "presence" variable from original dataset
glmDataBuff$Pres<-sciurusData$Pres


#initiate lists for glm results and loglik values
glm_list <- list()
ll_values <- numeric()

#iterate over radii with glms
for (r in radii) {
  
  formula <- formula(paste("Pres ~ w", r, sep = ""))
  
  # Running the GLM to test the relationship between the broadleaf cover and the different buffer size (radius)
  glm_model <- glm(formula, family = "binomial", data = glmDataBuff)
  
  # Store the model in the list
  glm_list[[as.character(r)]] <- glm_model
  
  # Get the log-likelihood and store it
  ll_values <- c(ll_values, as.numeric(logLik(glm_model)))
}

# Create data frame of results
glmRes <- data.frame(dist = radii, ll = ll_values)


#####################################

#Now re-build the function above but with the FAB correction

landBufferFAB <- function(x){         
  
  
  SciurusPoints<-vect(sciurusData[x,],geom=c("x","y"),crs="epsg:27700")  
  coverFAB<-list()
  for(i in radii){  #for loop to iterate over the radius distances
    
    sciurusBuffer <- buffer(SciurusPoints, width=i)           #buffer each point
    bufferlandcover <- crop(broadleaf, sciurusBuffer)         #crop the landcover layer to the buffer extent
    
    masklandcover <- mask(bufferlandcover, sciurusBuffer)            # mask the above to the buffer
    sciurusD<-distance(rasterize(SciurusPoints,masklandcover,field=1))# get distance of each cell in the buffer from the focal point
    
    sciurusD[sciurusD==0]<-1 # this is just to avoid zeroes messing up calculations later on
    
    sciurusD <- crop(sciurusD, sciurusBuffer)#crop the broadleaf layer to the extent of the buffer 
    sciurusD <- mask(sciurusD, sciurusBuffer)# mask to buffer
    
    
    sciurusD<- i/sciurusD # to prepare the FAB correction (buffer radius divided by distance from focal point)
    
    bufferlandcover5kmFAB<-bufferlandcover*sciurusD # apply the correction to the data
    
    woodSum<-sum(values(bufferlandcover5kmFAB),na.rm=T)# sum the corrected values 
    weightSum<-sum(values(sciurusD),na.rm = T) #sum all the weights
    
    fabCorrectedPC<-woodSum/weightSum*100 # get the corrected result as a percentage
    
    coverFAB[[i]]<-fabCorrectedPC} # add to list of results
  
  print(x)  #(to gauge progress)
  return(unlist(coverFAB))} #get the list of percentages 


# run FAB correction
woodlandFAB<-lapply(nPoints,FUN=landBufferFAB)

# bind all results
glmDataFAB<-do.call("rbind",woodlandFAB)

#convert to data frame
glmDataFAB<-data.frame(glmDataFAB)

#set column names
colnames(glmDataFAB)<-c(paste("w", radii, sep = ""))

#add in presence data
glmDataFAB$Pres<-sciurusData$Pres

#initiate lists for model results
glm_list <- list()
ll_values <- numeric()

#iterate of radii with glms
for (r in radii) {
  
  formula <- formula(paste("Pres ~ w", r, sep = ""))
  
  # Running the GLM to test the relationship between the broadleaf cover and the different buffer size (radius)
  glm_model <- glm(formula, family = "binomial", data = glmDataFAB)
  
  # Store the model in the list
  glm_list[[as.character(r)]] <- glm_model
  
  # Get the log-likelihood and store it
  ll_values <- c(ll_values, as.numeric(logLik(glm_model)))
}

# Create data frame
glmResFAB <- data.frame(dist = radii, ll = ll_values)


################################################
#Gaussian kernel approach

#function for iterating over all points and radii and generating Gaussian kernel
gaussFun10<-function(x){
  
  #create spatial point from coordinates
  SciurusPoints<-vect(sciurusData[x,],geom=c("x","y"),crs="epsg:27700")  
  
  #initiate list to hold land-cover data
  coverGauss<-list()
  for(i in radii){  #for loop to iterate over the radius distances
    
    sciurusBuffer <- buffer(SciurusPoints, width=i)     #buffer each point
    bufferlandcover <- crop(broadleaf, sciurusBuffer)   #crop the landcover layer to the buffer extent
    
    
    fw.i <- focalMat(bufferlandcover, c(i/2,i),type= 'Gauss')  #buffer around each cell of the broadleaf layer, applying sigma and radius distance
    
    #this next line makes sure that Gauusian filter matrix and woodland raster always have same dimensions
    fw.i<-fw.i[-c(1:(nrow(fw.i)-nrow(bufferlandcover))),-c(1:(ncol(fw.i)-ncol(bufferlandcover)))]
    
   
    maskMat<-matrix(bufferlandcover,nrow=nrow(bufferlandcover),ncol=ncol(bufferlandcover))
   
    focalBroadleaf<-maskMat*fw.i #the focal() function now creates the Gaussian filter raster
    
    #sum all values from Gaussian raster
    sumG<-sum(focalBroadleaf)
    
    # get percentage by dividing by all weights and multiplying by 100
    pcG<-sumG/sum(fw.i)*100
    
    
    coverGauss[[i]]<-pcG} #store in a list
  
  
  print(x)  #(to gauge progress)
  
  return(unlist(coverGauss))} #get the list of percentages 


# run Gaussian function
woodlandGauss<-lapply(nPoints,FUN=gaussFun10)

#bind all results
glmDataGauss<-do.call("rbind",woodlandGauss)

#convert to data frame
glmDataGauss<-data.frame(glmDataGauss)

#set column name
colnames(glmDataGauss)<-c(paste("w", radii, sep = ""))

#add in presence data
glmDataGauss$Pres<-sciurusData$Pres

#initiate lists for modelling
glm_list <- list()
ll_values <- numeric()

#iterate over results with glms
for (r in radii) {
  
  formula <- formula(paste("Pres ~ w", r, sep = ""))
  
  # Run the GLM
  glm_model <- glm(formula, family = "binomial", data = glmDataGauss)
  
  # Store the model in the list
  glm_list[[as.character(r)]] <- glm_model
  
  # Get the log-likelihood and store it
  ll_values <- c(ll_values, as.numeric(logLik(glm_model)))
}

# Create dataframe
glmResGauss <- data.frame(dist = radii, ll = ll_values)

head(glmResGauss)
plot(glmResGauss)

############################################
#FAB-corrected Gaussian
############################################
gaussFunFAB<-function(x){
  #create spatiak point
  SciurusPoints<-vect(sciurusData[x,],geom=c("x","y"),crs="epsg:27700")  
  #initiate list to store land-cover results
  coverGauss<-list()
  for(i in radii){  #for loop to iterate over the radius distances
    
    
    sciurusBuffer <- raster::buffer(SciurusPoints, width=i)           #buffer each point
    bufferlandcover <- raster::crop(broadleaf, sciurusBuffer)         #crop the landcover layer to the buffer extent
    
    
    sciurusR<-terra::rasterize(SciurusPoints,bufferlandcover,field=1)
    sciurusD<-distance(sciurusR)# get distance of each cell in the buffer from the focal point
    
    sciurusD[sciurusD==0]<-1# this is just to avoid zeroes messing up calculations later on
    
    sciurusD<-(i/sciurusD) # to prepare the FAB correction (buffer radius divided by distance from focal point)
    
    #get woodland layer to right extent
    bufferlandcover<-crop(bufferlandcover,sciurusD)
    
    #convert landcover to matrix
    bufferlandcoverMat<-matrix(bufferlandcover,nrow=nrow(bufferlandcover),ncol=ncol(bufferlandcover))
    
    #convert distance raster to matrix
    fabMat<-matrix(sciurusD,nrow=nrow(sciurusD),ncol=ncol(sciurusD))
    
    #get FAB-corrected land-cover
    correctedFab<-fabMat*bufferlandcoverMat
    
    
    #buffer around each cell of the broadleaf layer, applying sigma and the buffer distance. 
    fw.i <- focalMat(sciurusD,c(i/2,i),type= 'Gauss')  
    
    #to ensure Gaussian filter and land-cover have same dimensions
    fw.i<-fw.i[-c(1:(nrow(fw.i)-nrow(sciurusD))),-c(1:(ncol(fw.i)-ncol(sciurusD)))]
    
    # 
    maskFab<-matrix(sciurusD,ncol=ncol(sciurusD),nrow=nrow(sciurusD))
    
    #FAB corrected Gaussian kernal on land-cover
    focalBroadleaf<-correctedFab*fw.i 
    
    #FAB corrected Gaussian weights
    focalFab<-maskFab*fw.i
    
    
    #sum of broadleaf values divided by all weights to get % cover
    pcG<-sum(focalBroadleaf)/sum(focalFab)*100
    
  
    coverGauss[[i]]<-pcG}# store in list
  
  
  print(x)  #(to gauge progress)
  
  return(unlist(coverGauss))} #return the list of percentages 



# run the function
woodlandGaussFab<-lapply(nPoints,FUN=gaussFunFAB)

#bind all results
glmDataGaussFAB<-do.call("rbind",woodlandGaussFab)

#convert to data frame
glmDataGaussFAB<-data.frame(glmDataGaussFAB)

#set up column names
colnames(glmDataGaussFAB)<-c(paste("w", radii, sep = ""))

#add in presence data
glmDataGaussFAB$Pres<-sciurusData$Pres


#initiate lists for model results
glm_list <- list()
ll_values <- numeric()

#iterate over radii and run glm models
for (r in radii) {
  
  formula <- formula(paste("Pres ~ w", r, sep = ""))
  
  # Run the GLM
  glm_model <- glm(formula, family = "binomial", data = glmDataGaussFAB)

  # Store model reulsts in a list
  glm_list[[as.character(r)]] <- glm_model
  # Getting the log-likelihood and storing it
  ll_values <- c(ll_values, as.numeric(logLik(glm_model)))
}

# Create dataframe
glmResGaussFAB <- data.frame(dist = radii, ll = ll_values)


################################Join all results to first model data frame
glmRes$FAB<-glmResFAB$ll # add FAB results
glmRes$Gauss<-glmResGauss$ll # add Gaussian results
glmRes$GaussFAB<-glmResGaussFAB$ll #add FAB-corrected Gaussian results

# for plotting
library(ggplot2)
library(ggpubr)

plotNet05<- ggplot(glmRes)+geom_smooth(aes(dist,ll,colour="buff"),
                                              method=lm,
                                              formula=y~poly(x,5), linewidth=1)+
  geom_point(aes(dist,ll,colour="buff"))+
  
  
  geom_smooth(aes(dist,FAB,colour="FAB"),
              method=lm,
              formula=y~poly(x,5), linewidth=1)+
  geom_point(aes(dist,FAB,colour="FAB"))+
  
  geom_smooth(aes(dist,Gauss,colour="Gauss"),
              method=lm,
              formula=y~poly(x,5), linewidth=1)+
  
  geom_point(aes(dist,Gauss,colour="Gauss"))+
  
  geom_smooth(aes(dist,GaussFAB,colour="GaussFab"),
              method=lm,
              formula=y~poly(x,5), linewidth=1)+
  
  geom_point(aes(dist,GaussFAB,colour="GaussFab"))+
  
  
  
  xlab("dist")+ylab("Log Lik")+scale_color_manual(name = "", 
                                                  
                                                  values = c("FAB" = "blue","buff"="red","Gauss"="black","GaussFab"="darkgreen"))+theme(text = element_text(size = 20)) +  
  theme_pubr(base_size = 15)+font("legend.text",size=18)

#produce plot
plotNet05




