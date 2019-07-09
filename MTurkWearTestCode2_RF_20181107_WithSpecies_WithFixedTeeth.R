#Here we are going to do the wear test, but we are going to combine dissertation and Mech Turk results to use in the training data set. 
#Hopefuly this gives us better results with more training data.  

#install_github("vbonhomme/Momocs")
library(geomorph)
library(EBImage)
library(shapes)
library(devtools)
library(Momocs)
library(plyr)
library(dplyr)
library(caret)
library(randomForest)
library(e1071)
library(caret)
set.seed(123456)
numHarm<-20

normSetting <- FALSE





#Using the older versions of these functions.  
efourier_i <- function(ef, nb.h, nb.pts = 120) {
  if (is.null(ef$ao))
    ef$ao <- 0
  if (is.null(ef$co))
    ef$co <- 0
  an <- ef$an
  bn <- ef$bn
  cn <- ef$cn
  dn <- ef$dn
  ao <- ef$ao
  co <- ef$co
  if (missing(nb.h)) {
    nb.h <- length(an)
  }
  theta <- seq(0, 2 * pi, length = nb.pts + 1)[-(nb.pts + 1)]
  #theta <- seq(0 - pi/2, 2 * pi - pi/2, length = nb.pts + 1)[-(nb.pts + 1)]
  hx <- matrix(NA, nb.h, nb.pts)
  hy <- matrix(NA, nb.h, nb.pts)
  for (i in 1:nb.h) {
    hx[i, ] <- an[i] * cos(i * theta) + bn[i] * sin(i * theta)
    hy[i, ] <- cn[i] * cos(i * theta) + dn[i] * sin(i * theta)
  }
  x <- (ao/2) + apply(hx, 2, sum)
  y <- (co/2) + apply(hy, 2, sum)
  coo <- cbind(x, y)
  colnames(coo) <- c("x", "y")
  return(coo)
}

#Note: 
##Before doing efourier and training the models and doing prediction.
##I need to rotate all shapes in training and test using procrustes so that they all have the same alignment.  



#This has been run on the Mech Turk data.  It only needs to be run once.  
# 
# #Import MTurk results
# folders<-list.files("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/results")
# 
# for (j in folders){print(j)
# #Pick out one folder
# folderTemp <- paste("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/results/",j,sep="")
# fff<-list.files(folderTemp,full=TRUE)
# 
# #Shell script
# #ImageMagick
# system(paste("cd ",folderTemp,'; for i in *.png ; do convert "$i" -background white -alpha remove "$i" ; done',sep=" "))
# #Convert to .jpg
# system(paste("cd ",folderTemp,'; for i in *.png ; do convert "$i" "${i%.*}.jpg" ; done',sep=" "))
# system(paste("cd ",folderTemp,'; for i in *.jpeg ; do convert "$i" "${i%.*}.jpg" ; done',sep=" "))
# #Convert to gray scale
# system(paste("cd ",folderTemp,'; for i in *.jpg; do convert "$i" -colorspace Gray "$i"; done',sep=" "))
# #trim
# system(paste("cd ",folderTemp,'; for i in *.jpg; do convert -trim "$i" "$i"; done',sep=" "))
# }

################################################
#Links DSCN numbers to Taxa
################################################
require(gdata)
keys<-list.files("./TeethWear/Extant teeth Images and Key/Key to bovid teeth/",full=TRUE)
keys<-keys[grep(".csv",keys)]
keysList<-list()
for (k in keys){
  temp<-read.csv(k)
  if(any(names(temp)=="sub")){
    temp<-subset(temp,select=-sub)
  }
  names(temp)[1:10]<-c("Tribe","Genus","Species","Institution","No","Origin","Element","M3","M2","M1")
  temp<-apply(temp,2,as.character)
  temp<-temp[temp[,"Tribe"]!="",]
  keysList[[k]]<-temp[,1:10]
}

keyFile<-as.data.frame(do.call(rbind,keysList))
keyFile$M3<-paste0("dscn",as.character(keyFile$M3))
keyFile$M2<-paste0("dscn",as.character(keyFile$M2))
keyFile$M1<-paste0("dscn",as.character(keyFile$M1))

keyFile$M1[keyFile$M1=="dscn"]<-NA
keyFile$M2[keyFile$M2=="dscn"]<-NA
keyFile$M3[keyFile$M3=="dscn"]<-NA

keyFile$M1<-toupper(keyFile$M1)
keyFile$M2<-toupper(keyFile$M2)
keyFile$M3<-toupper(keyFile$M3)

#Remove the L from the names
keyFile$M1 <- gsub("L","",keyFile$M1)
keyFile$M2 <- gsub("L","",keyFile$M2)
keyFile$M3 <- gsub("L","",keyFile$M3)

################################################
#test Data from mechanical turk.  
#This is the data to be classified
################################################

imgsList<-list()
dscn<-list.files("./TeethWear/resultsClean")

for (d in dscn){print(d)
  
  files<-list.files(paste0("./TeethWear/resultsClean/",d))
  files <- files[grep("jpg",files)]
  
  imgsList[[d]]<-import_jpg(paste0("./TeethWear/resultsClean/",d,"/",files,sep=""))
  
  #Manual removal of this MTurk image
  imgsList[[d]][which(names(imgsList[[d]])=="3A1COHJ8NJU7FB7Z1T98O6918IG8H5")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="3C8HJ7UOP7T8RL9X1GPYTVE1ODMZMA")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="3R0T90IZ1SBVRI21YZ7V5STJJ9PGCL")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="3NJM2BJS4W514VV01IXIZ17BK0BCPK")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="37C0GNLMHF23ZHJ9MITKD7YCASHD6T")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="3N1FSUEFL5ZPKUFV3U05G9EYE454DZ")]<-NULL
  
  #Remove images with too few observations.  These are clearly mistakes.  
  
  remove<-names(which(unlist(lapply(imgsList[[d]],nrow))<=50))
  if (length(remove)>0){
    for (q in remove){
      imgsList[[d]][[q]]<-NULL
    }
  }
  
  
}

dscn<-names(imgsList)
imgsListMean<-list()
for (d in dscn){print(d)
  
  outCoo <- Out(imgsList[[d]])
  outCoo <- coo_slidedirection(outCoo, "left")
  #"DSCN2650" "DSCN2672" "DSCN3727" "DSCN3924" "DSCN5454" "DSCN5692" "DSCN5696"
  if (d %in% c("dscn4808", "dscn6009", "dscn2931","dscn3019","dscn3289l","dscn3300l","dscn2616","dscn2628","dscn2631",
               "dscn2644", "dscn2624","dscn2640","dscn2879","dscn3080l","dscn3296l","dscn3950",
               "dscn4151","dscn4796","dscn4817","dscn7398","dscn7418",
               "dscn2650","dscn2672","dscn3727","dscn3924","dscn5454","dscn5692","dscn5696")){
    for (k in 1: length(outCoo$coo)){
      which.min(abs(outCoo$coo[[k]][,2]))
      
      start <- min(which.min(outCoo$coo[[k]][,1]))
      print(start)
      outCoo$coo[[k]] <- coo_slide(outCoo$coo[[k]],start)
    }
  }
  
  test <- efourier(outCoo,nb.h=numHarm, norm = normSetting)
  
  numPoints<-200
  ptsList<-list()
  arr<-array(NA, dim=c(numPoints, 2, dim(test$coe)[1]))
  for (j in 1:dim(test$coe)[1]){
    #arr[,,j]<-ptsList[[j]]<-efourier_shape(test$coe[j,grep("A",colnames(test$coe))],test$coe[j,grep("B",colnames(test$coe))],test$coe[j,grep("C",colnames(test$coe))],test$coe[j,grep("D",colnames(test$coe))],plot=TRUE,nb.pts=numPoints)
    efa<-list(an=test$coe[j,grep("A",colnames(test$coe))],
              bn=test$coe[j,grep("B",colnames(test$coe))],
              cn=test$coe[j,grep("C",colnames(test$coe))],
              dn=test$coe[j,grep("D",colnames(test$coe))],
              a0=0,c0=0)
    arr[,,j]<-ptsList[[j]]<-efourier_i(efa,nb.pts=numPoints)
  }
  
  # plot(arr[,,1],type="l", main = d,col = "white")
  # for (g in 1:2){
  # points(arr[1:10,,g],type="l")
  # }
  # points(imgsListMean[[d]],type="l",col="pink")
  
  if (length(ptsList)==1){
    imgsListMean[[d]] <- ptsList[[1]]
  }
  
  if (length(ptsList)>1){
    #Average across ptsList
    imgsListMean[[d]] <- ptsList[[1]]
    all <- do.call(cbind,ptsList)
    #Average across Juliet ands MTurk Workers
    xxx<-seq(1,dim(all)[2],2)
    yyy<-seq(1,dim(all)[2],2)+1
    
    keep <- abs(scale(all[1,xxx]))<1
    xxx <- xxx[keep]
    yyy <- yyy[keep]
    
    imgsListMean[[d]][,1]<-apply(all[,xxx],1,mean)
    imgsListMean[[d]][,2]<-apply(all[,yyy],1,mean)
  }
  
  
}

# #set of shapes to be predicted.
#imgsListMean[[d]]<-mshape<-as.matrix(procGPA(arr, scale = FALSE, reflect = FALSE)$mshape)
# 
# plot(imgsList[["dscn3468"]][[1]],main = "dscn3468",type="l")
# points(imgsList[["dscn3468"]][[2]],type="l",col = "red")
# points(imgsList[["dscn3468"]][[3]],type="l", col = "orange")
# points(imgsList[["dscn3468"]][[4]],type="l", col = "gold")
# points(imgsList[["dscn3468"]][[5]],type="l", col = "darkgreen")
# points(imgsList[["dscn3468"]][[6]],type="l", col = "blue")
# points(imgsList[["dscn3468"]][[7]],type="l", col = "pink")
# 
# plot(imgsList[["dscn3481"]][[1]],main = "dscn3481",type="l",lwd =3)
# points(imgsList[["dscn3481"]][[2]],type="l",col = "red")
# points(imgsList[["dscn3481"]][[3]],type="l", col = "orange")



names(imgsListMean) <- toupper(gsub("l","",names(imgsListMean)))


#Checking the teeth
# for (i in names(imgsListMean)){
#   png(paste0("/Users/gregorymatthews/Dropbox/brophyTeeth/WearTestCheck/",i,".png"))
#   plot(imgsListMean[[toupper(i)]], main = i, asp = 1)
#   dev.off()
# }




#save(imgsList,file="/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/imgsList_20181108.RData")
#load("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/imgsList_20181108.RData")
#save(imgsListMean,file="/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/imgsListMean_20181108.RData")
#load("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/imgsListMean_20181108.RData")
################################################
#Training Data
#Note: Delete tooth 5733 from the training data set entirely.  
################################################
ptsTrainList<-list()
refFile<-data.frame()

setwd("./TeethWear/Raw data/Point files/")

tribeVec <- list.files()
for (t in tribeVec){print(t)
  speciesVec<-list.files(paste0("./TeethWear/Raw data/Point files/",t))
  for (s in speciesVec){print(s)
    teethVec<-list.files(paste0("./TeethWear/Raw data/Point files/",t,"/",s))
    for (tooth in teethVec){
      filesVec<-list.files(paste0("./TeethWear/Raw data/Point files/",t,"/",s,"/",tooth))
      for (f in filesVec){
        file<-(paste0("./TeethWear/Raw data/Point files/",t,"/",s,"/",tooth,"/",f))
        
        if (file != "./TeethWear/Raw data/Point files/Antilopini/Antidorcas marsupialis/LM2 A marsupialis/DSCN3190"){
          temp<-matrix(unlist(read.table(file))[-c(1:20)],ncol=13,byrow=TRUE)[,-13]
          
          pts<-data.frame(x=c(t(temp[,seq(1,11,2)])),y=c(t(temp[,seq(1,11,2)+1])))
        }
        
        if (file == "./TeethWear/Raw data/Point files/Antilopini/Antidorcas marsupialis/LM2 A marsupialis/DSCN3190"){
          pts<-read.table(file)
          names(pts)<-c("x","y")
        }
        
        #Check if the tooth is in the test data set.  If it is don't add it to the training data set.
        if (!f%in%toupper(names(imgsListMean))){
          ptsTrainList[[substring(tooth,1,3)]][[f]]<-as.matrix(pts)
          refFile<-rbind(refFile,data.frame(ref=f,tooth=substring(tooth,1,3),tribe=t,species=s))
        }
        
      }
    }
  }
}


sp <- as.character(unique(keyFile$Species))[-c(17,22)]
for (s in sp){
  #Tooth check by species 
  ids <- keyFile[keyFile$Element == "maxilla" & keyFile$Species == s,]$M3
  temp <- ptsTrainList[["UM3"]][names(ptsTrainList[["UM3"]]) %in% ids]
  
  
  plot(scale(temp[[1]],scale = FALSE),type= "l", main = paste0("UM3",s), xlim = c(-500,500), ylim = c(-500,500))
  for (i in 2:length(temp)){
    points(scale(temp[[i]], scale = FALSE),type = "l")
  }
}

# efourier(ptsTrainList[[1]][[1]])
# efourier(imgsListMean[[1]])
# efourier(ptsTrainList[["LM1"]][["DSCN5598"]])
# 
# plot(scale(ptsTrainList[[3]][[1]], scale = FALSE),type="l")
# points(scale(imgsListMean[[2]], scale = FALSE),type="l")
# points(scale(ptsTrainList[["LM1"]][["DSCN5598"]], scale = FALSE), type = "l")

############################################################################
#These are teeth that were scaled wrong in the dissertation data.  They were re
#done by Juliet.    
#These are in the folder called "fixedTeeth"
############################################################################
fixedImgsList<-list()
for (tooth in c("LM1","LM2","LM3","UM1","UM2")){print(tooth)
  dscn<-list.files(paste0("./TeethWear/fixedTeeth/",tooth))
  for (d in dscn){
    id <- gsub(".jpg","",d)
    fixedImgsList[[tooth]][[id]] <- import_jpg(paste0("/Users/gregorymatthews/Dropbox/brophyTeeth/fixedTeeth/",tooth,"/",d,sep=""))[[1]]  
    ptsTrainList[[tooth]][[id]] <- fixedImgsList[[tooth]][[id]]
  }
}

#Manually remove a tooth per request of juliet
ptsTrainList[["UM2"]][["DSCN5733"]] <- NULL


#Also delete these teeth per request of Juliet: DSCN3284, DSCN4939, DSCN4967, DSCN5004, DSCN5008, DSCN5031
for (tooth in c("LM1","LM2","LM3","UM1","UM3")){
  for (i in c("DSCN3284","DSCN4939","DSCN4967","DSCN5004","DSCN5008","DSCN5031")){
    ptsTrainList[[tooth]][[i]] <- NULL
  }
}



# plot(ptsTrainList[[1]][[1]], xlim = c(0,2000), ylim = c(0,2000))
# for (i in 1:length(ptsTrainList[[1]])){
# points(ptsTrainList[[1]][[i]], type = "l")
# }
# 
# #plot(imgsListMean[[1]], xlim = c(-500,500),ylim = c(-500,500))
# for (i in 1:length(imgsListMean)){
#   points(imgsListMean[[i]]+c(500,500), type = "l",col = "red")
# }

#Find the missing observations
#setdiff(toupper(names(imgsListMean)),as.vector(unlist(lapply(outList,function(x){x$ref}))))


#Now train the model.  
set.seed(1234)
outList <- list()
ptsList <- mod <- modSpecies <- efTrainList <- list()
predTempProbsSpecies <- list()
rm(tooth)
#Remove observatiobs from the training set that are in the test set.  
#No UM2 in the test data set.  
for (tooth in c("LM1","LM2","LM3","UM1","UM3")){print(tooth)
  #Align the training and test sets using procrustes so that their Fourier coefficients all have the same interpretation
  modSpecies[[tooth]] <- list()
  
  #Subset the Test data set to include only the ones for the target tooth 
  if (substring(tooth,1,1)=="U"){
    ids<-which(toupper(names(imgsListMean))%in%toupper(subset(keyFile,Element=="maxilla")[[substring(tooth,2,3)]]))
  }
  
  if (substring(tooth,1,1)=="L"){
    ids<-which(toupper(names(imgsListMean))%in%toupper(subset(keyFile,Element=="mandible")[[substring(tooth,2,3)]]))
  }
  
  #Create the Coo object for Training and test sets. 
  outCooTrain <- Out(ptsTrainList[[tooth]])
  outCooTest <- Out(imgsListMean[ids])
  
  outCooTest <- coo_slidedirection(outCooTest, "right")
  outCooTrain <- coo_slidedirection(outCooTrain, "right")
  
  
  #Perform efourier for test and training data set.  
  efTrain <- as.data.frame(efourier(outCooTrain,nb.h=numHarm, norm = normSetting)$coe)
  efTest <- as.data.frame(efourier(outCooTest,nb.h=numHarm, norm = normSetting)$coe)
  
  #Perform PCA.  
  pc<-princomp(rbind(efTrain,efTest))
  efTrain<-cbind(efTrain,pc$scores[1:nrow(efTrain),1:20])
  efTest<-cbind(efTest,pc$scores[-c(1:nrow(efTrain)),1:20])
  
  efTrain$ref<-names(ptsTrainList[[tooth]])
  efTest$ref<-toupper(names(imgsListMean)[ids])
  
  efTrainList[[tooth]] <- efTrain<-merge(efTrain,refFile,by.x="ref",by.y="ref",all.x=TRUE)
  
  #####################
  #Train the model
  #####################
  form<-as.formula(paste0("tribe~",paste(c(paste0("A",1:numHarm),paste0("B",1:numHarm),paste0("C",1:numHarm),paste0("D",1:numHarm),paste0("Comp.",1:20)),collapse="+")))
  #Omit PCS
  #form<-as.formula(paste0("tribe~",paste(c(paste0("A",1:numHarm),paste0("B",1:numHarm),paste0("C",1:numHarm),paste0("D",1:numHarm)),collapse="+")))
  
  #tune.control <- tune.control(cross=6)
  #TribeInnerTune.svm<-tune.svm(y=efTrain$tribe,x=efTrain[,c(paste0("A",1:numHarm),paste0("B",2:numHarm),paste0("C",2:numHarm),paste0("D",1:numHarm),paste0("Comp.",1:20))],gamma=2^seq(-10,-3,length=15),cost=10^seq(-2,1.5,length=15),tune.control=tune.control)
  #print(c(TribeInnerTune.svm$best.parameters$gamma,TribeInnerTune.svm$best.parameters$cost))
  
  #Train a random forest
  mod[[tooth]]<-rfFitted<-randomForest(form,data=efTrain,ntree=5000)
  
  
  formula<-as.formula(paste0("species~",paste(c(paste0("A",1:numHarm),paste0("B",1:numHarm),paste0("C",1:numHarm),paste0("D",1:numHarm),paste0("Comp.",1:20)),collapse="+")))
  tribes <- as.character(sort(unique(efTrain$tribe)))
  for (j in c(1,4:7)){
    temp <- subset(efTrain, efTrain$tribe == tribes[j])
    temp$species<-as.factor(as.character(temp$species))
    
    #sampsize <- round(dim(temp)[1]*sampSizeVecSpecies[m])
    # #################################
    # #Tuning for Species
    # #################################
    # SpecInnerTune.rf <- tune.randomForest(formula,data = temp,ntree = 2000,mtry = seq(10,50,10))
    # print(SpecInnerTune.rf$best.parameters$mtry)
    # tuningParmList[[tooth[m]]][["Species"]][[j]][[i]]<-SpecInnerTune.rf$best.parameters
    
    modSpecies[[tooth]][[j]] <- cart <- randomForest(formula,data=temp,ntree=5000)
    
  }
  
  #################################
  #Fit with best tunes parameters
  #################################
  #svmFitted <- svm(y=efTrain$tribe,x=efTrain[,c(paste0("A",1:numHarm),paste0("B",2:numHarm),paste0("C",2:numHarm),paste0("D",1:numHarm),paste0("Comp.",1:20))],probability=TRUE,gamma = TribeInnerTune.svm$best.parameters$gamma,cost = TribeInnerTune.svm$best.parameters$cost)
  
  #Which observations in the test set are a specific tooth.  
  #ids<-which(toupper(names(imgsListMean))%in%toupper(keyFile$M1) |
  #             toupper(names(imgsListMean))%in%toupper(keyFile$M2) |
  #             toupper(names(imgsListMean))%in%toupper(keyFile$M3)) 
  
  
  
  #predict the tribes
  predTempProbs<-predict(rfFitted,efTest[,c(paste0("A",1:numHarm),paste0("B",1:numHarm),paste0("C",1:numHarm),paste0("D",1:numHarm),paste0("Comp.",1:20))],type="prob")
  predTempProbs<-t(apply(predTempProbs,1,function(x){
    x[x==1]<-0.999;x[x==0]<-0.001
    #Now normalize
    x <- x/sum(x)
    return(x)
  }))
  
  
  #predict the species
  predTempProbsSpecies[[tooth]]<-list()
  for (i in c(1:7)){
    
    if (i %in% c(2:3)){
      predTempProbsSpecies[[tooth]][[i]] <- as.matrix(predTempProbs[,tribes[i]])
    }
    
    if (i %in% c(1,4:7)){
      predTempProbsSpecies[[tooth]][[i]]<-predict(modSpecies[[tooth]][[i]] ,efTest[,c(paste0("A",1:numHarm),paste0("B",1:numHarm),paste0("C",1:numHarm),paste0("D",1:numHarm),paste0("Comp.",1:20))],type="prob")
      for (q in 1:nrow(predTempProbsSpecies[[tooth]][[i]])){  
        predTempProbsSpecies[[tooth]][[i]][q,] <- predTempProbsSpecies[[tooth]][[i]][q,] * predTempProbs[q,tribes[i]]
      }
    }
  }
  
  
  predSpecies <- do.call(cbind,predTempProbsSpecies[[tooth]])
  
  #use ids to pull out only the target tooth 
  outList[[tooth]]<-data.frame(ref=toupper(names(imgsListMean)[ids]),predTempProbs,predSpecies)
  
}


#getting species count 
#checkspecies <- do.call(rbind,efTrainList)

######################################################
#Pull in ID's on the chipps vs worn teeth.  
######################################################
fils<-toupper(list.files("./TeethWear/Wear Test/chipped/"))
fils<-fils[grep("JPG",fils)]
filsChipped<-substring(fils,1,nchar(fils)-4)

fils<-c(toupper(list.files("./TeethWear/Wear Test/notchipped/")),list.files("/Users/gregorymatthews/Dropbox/Wear Test/notchipped/Incorrect/"))
fils<-fils[grep("JPG",fils)]
filsNotChipped<-substring(fils,1,nchar(fils)-4)

######################################################################
#Now pull out only the teeth we need from each
#Prediction are made on all test teeth but we only need the ones for that specific model
outList[["LM1"]]<-subset(merge(outList[["LM1"]][,1:28],keyFile[,c("Tribe","Genus","Species","Element","M1")],by.x="ref",by.y="M1",all.x=TRUE),Element=="mandible")
outList[["LM2"]]<-subset(merge(outList[["LM2"]][,1:28],keyFile[,c("Tribe","Genus","Species","Element","M2")],by.x="ref",by.y="M2",all.x=TRUE),Element=="mandible")
outList[["LM3"]]<-subset(merge(outList[["LM3"]][,1:28],keyFile[,c("Tribe","Genus","Species","Element","M3")],by.x="ref",by.y="M3",all.x=TRUE),Element=="mandible")

outList[["UM1"]]<-subset(merge(outList[["UM1"]][,1:28],keyFile[,c("Tribe","Genus","Species","Element","M1")],by.x="ref",by.y="M1",all.x=TRUE),Element=="maxilla")
#outList[["UM2"]]<-subset(merge(outList[["UM2"]][,1:8],keyFile[,c("Tribe","Genus","Species","Element","M2")],by.x="ref",by.y="M2",all.x=TRUE),Element=="mandible")
outList[["UM3"]]<-subset(merge(outList[["UM3"]][,1:28],keyFile[,c("Tribe","Genus","Species","Element","M3")],by.x="ref",by.y="M3",all.x=TRUE),Element=="maxilla")


outList[["LM1"]]$tooth <- "LM1"
outList[["LM2"]]$tooth <- "LM2"
outList[["LM3"]]$tooth <- "LM3"

outList[["UM1"]]$tooth <- "UM1"
#outList[["UM2"]]$tooth <- "UM2"
outList[["UM3"]]$tooth <- "UM3"

preds<-do.call(rbind,outList)

names(preds)[13:14] <- c("Antidorcas.marsupialis","Syncerus.cafer")

names(preds)[9:28] <- do.call(rbind,strsplit(names(preds)[9:28],"[.]"))[,2]

preds$Species <- gsub(" ","",preds$Species)

################################################
#Merge on indicator for chipped vs not chipped.
################################################
preds$Chipped<-preds$ref%in%filsChipped+0
preds$Chipped[preds$ref%in%filsChipped]<-"Chipped"
preds$Chipped[!preds$ref%in%filsChipped]<-"Not Chipped"


preds$probTribe <- NA
for (tr in as.character(unique(preds$Tribe))){
  preds$probTribe[preds$Tribe==tr] <- preds[[tr]][preds$Tribe==tr]
}

preds$probSpecies <- NA
for (sp in as.character(unique(preds$Species))){
  preds$probSpecies[preds$Species==sp] <- preds[[sp]][preds$Species==sp]
}

preds$predictedTribe<-c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")[apply(preds[,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")],1,which.max)]
preds$correct<-(preds$predictedTribe==preds$Tribe)+0

#Check species correct
preds$predictedSpecies<-c("gnou","buselaphus","dorcas","oryx","niger","taurinus","campestris","gazella","capreolus","equinus" ,"arundinum","cafer")[apply(preds[,c("gnou","buselaphus","dorcas","oryx","niger","taurinus","campestris","gazella","capreolus","equinus" ,"arundinum","cafer")],1,which.max)]
preds$correctSpecies<-(preds$predictedSpecies==preds$Species)+0

# 
 write.csv(preds,file="/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/classProbs20181109.csv")
# write.csv(preds,file="/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/results20180604_RandomForest_withSpecies.csv")
 save(preds,file="/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/preds_results20181109_RandomForest_withSpecies.RData")
 save(mod,file="/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/mod_results20181109_RandomForest_withSpecies.RData")
save.image(file="/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/mod_results20181109_RandomForest_withSpecies_ALL.RData")




########################################################################
#SCRAP
########################################################################
# temp<-subset(preds,Tribe=="Alcelaphini")
# boxplot(temp$Alcelaphini~temp$Chipped)
# library(dplyr)
# group_by(preds,tooth,Tribe,Chipped)  %>% summarize(mn=mean(prob),n=n())
# group_by(preds,tooth,Chipped)  %>% summarize(mn=mean(prob),n=n())
# 
# boxplot(preds$prob~preds$tooth)
# boxplot(preds$prob~preds$Tribe)
# 
# preds$correct<-(preds$predictedTribe==preds$Tribe)+0
# 
# 
# library(ggplot2)
# ggplot(data=preds) + geom_boxplot(aes(x=Tribe,y=probTribe)) + geom_point(aes(y=probTribe,x=Tribe,colour=as.factor(correct))) 
# ggplot(data=preds) + geom_boxplot(aes(x=Tribe,y=prob)) + geom_point(aes(y=prob,x=Tribe,colour=as.factor(correct))) + facet_wrap(~Chipped)
# ggplot(data=preds) + geom_boxplot(aes(x=predictedTribe,y=prob)) + geom_point(aes(y=prob,x=predictedTribe,colour=as.factor(correct))) 
# ggplot(data=preds) + geom_boxplot(aes(x=predictedTribe,y=prob)) + geom_point(aes(y=prob,x=predictedTribe,colour=as.factor(correct))) + facet_wrap(~Chipped)
# 
# 
# 
# 
# ##############################
# #Overall
# ##############################
# #which.max(preds[1,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")])
# preds$predictedTribe<-c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")[apply(preds[,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")],1,which.max)]
# table(preds$predictedTribe,preds$Tribe)
# confusionMatrix(preds$predictedTribe,preds$Tribe)
# 
# 
# table(preds$predictedTribe[preds$tooth=="LM1"],preds$Tribe[preds$tooth=="LM1"])
# table(preds$predictedTribe[preds$tooth=="LM2"],preds$Tribe[preds$tooth=="LM2"])
# table(preds$predictedTribe[preds$tooth=="LM3"],preds$Tribe[preds$tooth=="LM3"])
# table(preds$predictedTribe[preds$tooth=="UM1"],preds$Tribe[preds$tooth=="UM1"])
# table(preds$predictedTribe[preds$tooth=="UM3"],preds$Tribe[preds$tooth=="UM3"])
# 
# ##############################
# #Chipped
# ##############################
# #which.max(preds[1,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")])
# preds$predictedTribe<-c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")[apply(preds[,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")],1,which.max)]
# 
# table(subset(preds,Chipped==1)$predictedTribe,subset(preds,Chipped==1)$Tribe)
# confusionMatrix(subset(preds,Chipped==1)$predictedTribe,subset(preds,Chipped==1)$Tribe)
# 
# 
# ##############################
# #Not Chipped
# ##############################
# #which.max(preds[1,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")])
# preds$predictedTribe<-c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")[apply(preds[,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")],1,which.max)]
# 
# table(subset(preds,Chipped==0)$predictedTribe,subset(preds,Chipped==0)$Tribe)
# confusionMatrix(subset(preds,Chipped==0)$predictedTribe,subset(preds,Chipped==0)$Tribe)
# 
# 
