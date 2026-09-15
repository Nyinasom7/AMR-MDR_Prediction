##########project files################
  
  dir.create("Data")
  
  ##### script files####
  dir.create("Script")
  
  ###### results#### 
  dir.create("Results")
  
  #################### SETUP- For reproducibility##########################
  library(tidyverse)
  library(readxl)
  library(caret)
  library(randomForest)
  library(e1071)
  library(pROC)
  #### for variable importance plot
  library(viridis) ## publication colors
  library(corrplot)
  library(ggpubr)
  
  set.seed(42)
  theme_set(theme_bw(base_size = 12)) ## publication theme
  sessionInfo()## save this for methods section
  
  da_raw<-read_excel("C:/Users/PROSPER/Desktop/AMR prediction/Data/da1.xlsx")
  
  df_raw<-read_excel("Data/da1.xlsx")
  
  data<-df_raw%>%rename(SPECIES_ISOLATED="SPECIES ISOLATED")
  ###replace NA with unknown Isolate####
  data<-data%>%mutate(SPECIES_ISOLATED=ifelse(is.na(SPECIES_ISOLATED),"Unknown isolate",SPECIES_ISOLATED))
  
  
  ##############PREPROCESSING############################
  ##clean data#####
  clean_data<-function(x){
    x<-as.character(x)
    x<-gsub("[>=<]","",x)
    x<-na_if(x, "-"); x <-na_if(x, "")
    as.numeric(x)
  }
  
  data<-data%>%mutate(across(ends_with("_MIC"),clean_data))
  data<-data%>%mutate(across(ends_with("_INT"),~toupper(trimws(.))))
  
  data$CEF_SCREEN<-factor(ifelse(df$CEF_SCREEN=="POS","R","S"),levels=c("S","R"))
  
  ##Impute by species median###
  data<-data%>%group_by(SPECIES_ISOLATED)%>%
    mutate(across(ends_with("_MIC"),~ifelse(is.na(.),median(.,na.rm=T),.)))%>%
    ungroup()
  
  data<-data%>%
    mutate(across(ends_with("_MIC"),~as.numeric(.)))%>%
    mutate(across(ends_with("_MIC"),~ifelse(is.na(.),median(.,na.rm=TRUE), .)))
  sum(is.na(data$CIP_MIC))
  #Log2 for ML
  data<-data%>%mutate(across(ends_with("_MIC") & !starts_with("LOG"), ~log2(.+0.125), .names="LOG_{.col}"))
  
  #######EDA##########################################################
  data%>%count(SPECIES_ISOLATED)%>%mutate(pct=n/sum(n)*100)
  data%>%count(SPECIMEN_TYPE)%>%mutate(pct=n/sum(n)*100)
  data%>%count(GENDER)%>%mutate(pct=n/sum(n)*100)
  
  
  ##species distribution##
  s1<-ggplot(data,aes(x=SPECIES_ISOLATED,fill=SPECIES_ISOLATED))+
    geom_bar()+geom_text(stat="count",aes(label=..count..),vjust=-0.5)+
    labs(title="Distribution of CoNS species(n=583) ",x="Species",y="Count")+scale_fill_viridis_d()+theme(legend.position ="none" )
  s1
  
  ggsave("Fig1A_Species_Dist.png",s1, width=6, height = 4,dpi=300)
  
  ## Resistance per antibiotic#####
  resist_r<-data%>%
    select(ends_with("_INT"))%>%
    pivot_longer(everything(),names_to="Antibiotic",values_to="Status")%>%
    group_by(Antibiotic)%>%
    summarise(R_rate =mean(Status=="R",na.rm=TRUE)*100, n=n())
  resist_r
  
  print(resist_r%>%arrange(desc(R_rate)))
  
  #### Resistance HeatMap ###
  s2<-resist_r%>%
    ggplot(aes(x=Antibiotic,y=R_rate,fill=R_rate))+geom_col()+coord_flip()
  
  s2

  
  #######MDR definition#####
  
  
  data<-data%>%
    mutate(c_beta=ifelse(BEN_INT=="R"|OXA_INT=="R"|CEF_SCREEN=="R",1,0),
           c_amin=ifelse(GEN_INT=="R",1,0),
           c_fluoro=ifelse(CIP_INT=="R",1,0),
           c_lin=ifelse(CLIN_INT=="R"|LIN_INT=="R",1,0),
           c_macro=ifelse(ERY_INT=="R",1,0),
           c_glyco=ifelse(VAN_INT=="R",1,0),
           c_tet=ifelse(TET_INT=="R",1,0),
           c_rif=ifelse(RIF_INT=="R",1,0),
           c_trim=ifelse(TRIM_INT=="R",1,0),
           c_nit=ifelse(NIT_INT=="R",1,0),
           n_classes=c_beta+c_amin+c_fluoro+c_macro+c_lin+c_glyco+c_rif+c_trim+c_nit,
           MDR=factor(ifelse(n_classes>=3,"MDR","non-MDR"),levels=c("non-MDR","MDR")))
       

  data2<-write.csv(data,"Clean_ready data")  
  print(table(data$MDR))  
  print(prop.table(table(data$MDR)))
#################FEATURE SELECTION#######################################
  
features<-c("GENDER","SPECIMEN_TYPE","SPECIES_ISOLATED",
            grep("^LOG_",names(data),value=TRUE))
  
  features
  
  df_ml<-data%>%select(all_of(features),MDR)%>%drop_na()
  
  df_ml
  write.csv(df_ml,"Clean_ready data.csv") 
  
  ##############Machine learning model #####################
  

  #### split####
  set.seed(123)
  Split_data<-createDataPartition(df_ml$MDR,p=0.8,list=FALSE)
train<-df_ml[Split_data,]
test<-df_ml[-Split_data,]
nrow(train)
nrow(test)

prop.table(table(train$MDR))
prop.table(table(test$MDR))

####Preprocess and Cross Validation###
pre<-preProcess(train%>%select(starts_with("LOG_")),method=c("center","scale"))
  
  #######Drop#############zero_variance LOG column first#####
nzv<-nearZeroVar(train%>%select(starts_with("LOG_")),saveMetrics = TRUE)

print(nzv[nzv$nzv,])#### see which ones
####keep only good Log cols#####
log_cols<-rownames(nzv)[!nzv$nzv]
print(log_cols)

pre<-preProcess(train%>%select(all_of(log_cols)),method=c("center","scale"))
  print(pre)
  
  train_scale<-predict(pre,train%>%select(all_of(log_cols)))
  test_scale<-predict(pre,test%>%select(all_of(log_cols)))

  train_scale$MDR<-as.factor(train$MDR)
                             
    test_scale$MDR<-as.factor(test$MDR)                         
    train_scale$MDR<-factor(train_scale$MDR)
    levels(train_scale$MDR)<-make.names(levels(train_scale$MDR))
    
    #check#
    print(levels(train_scale$MDR))
    
    
 ####balanced 96% problem####
train_bal<-upSample(x=train_scale%>%select(-MDR),
                    y=train_scale$MDR,yname="MDR")
table(train_bal$MDR)

###Ctrl###

ctrl<-trainControl(method="cv",number=5,
                   classProbs = TRUE,summaryFunction = twoClassSummary,
                   savePredictions = "final")


####model implementation######

#Random Forest#
model_rf<-train(MDR~.,data=train_bal,method="rf",trControl=ctrl,metric="ROC")

print(model_rf)
varImp(model_rf)##Import

####Logistic Regression###
model_glm<-train(MDR~.,data=train_bal,method="glm",family="binomial",trControl=ctrl,metric="ROC")
print(model_glm)


#####Support Vector machine######
model_svm<-train(MDR~.,train_bal,method="svmRadial",trControl=ctrl,metric="ROC")
print(model_svm)



####comparison####
results<-resamples(list(
  RF=model_rf,
  GLM=model_glm,
  SV=model_svm))

summary(results)
dotplot(results,metric="ROC")
dotplot(results,metric="Sens")
dotplot(results,metric="Spec")



##### Model performance validation#########################
test_scale$MDR<-factor(test_scale$MDR)
pred_rf<-factor(pred_rf,levels=levels(test_scale$MDR))
pred_SV<-factor(pred_SV,levels=levels(test_scale$MDR))


confusionMatrix(pred_rf,test_scale$MDR,positive="MDR")
