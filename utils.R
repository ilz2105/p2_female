CreateStudyFactorAndModel=function(StudyFactorModelCombo,gender){
  require(tidyverse)
  if(gender=="Male"){
    StudyFactorModelCombo=StudyFactorModelCombo[-6]#drop thyroid
    StudyFactorAndModel=NULL
    for (i_1 in 1:length(StudyFactorModelCombo[[1]])) {
      for (i_2 in 1:length(StudyFactorModelCombo[[2]])) {
        for (i_3 in 1:length(StudyFactorModelCombo[[3]])) {
          for (i_4 in 1:length(StudyFactorModelCombo[[4]])) {
            for (i_5 in 1:length(StudyFactorModelCombo[[5]])) {
              for (i_6 in 1:length(StudyFactorModelCombo[[6]])) {
                tmp=c(StudyFactorModelCombo[[1]][i_1],StudyFactorModelCombo[[2]][i_2],StudyFactorModelCombo[[3]][i_3],
                      StudyFactorModelCombo[[4]][i_4],StudyFactorModelCombo[[5]][i_5],StudyFactorModelCombo[[6]][i_6])
                StudyFactorAndModel=cbind(StudyFactorAndModel,tmp)
              }
            }
          }
        }
      }
    }
  }else{#female
    StudyFactorAndModel=NULL
    for (i_1 in 1:length(StudyFactorModelCombo[[1]])) {
      for (i_2 in 1:length(StudyFactorModelCombo[[2]])) {
        for (i_3 in 1:length(StudyFactorModelCombo[[3]])) {
          for (i_4 in 1:length(StudyFactorModelCombo[[4]])) {
            for (i_5 in 1:length(StudyFactorModelCombo[[5]])) {
              for (i_6 in 1:length(StudyFactorModelCombo[[6]])) {
                for (i_7 in 1:length(StudyFactorModelCombo[[7]])) {
                  tmp=c(StudyFactorModelCombo[[1]][i_1],StudyFactorModelCombo[[2]][i_2],StudyFactorModelCombo[[3]][i_3],
                        StudyFactorModelCombo[[4]][i_4],StudyFactorModelCombo[[5]][i_5],StudyFactorModelCombo[[6]][i_6],
                        StudyFactorModelCombo[[7]][i_7])
                  StudyFactorAndModel=cbind(StudyFactorAndModel,tmp)
                }
              }
            }
          }
        }
      }
    }
  }
  NumOfModel=ncol(StudyFactorAndModel)
  colnames(StudyFactorAndModel)=str_c("model",1:NumOfModel)
  
  StudyFactorMat=StudyFactorModelMat=matrix(NA_character_,nrow=nrow(StudyFactorAndModel),ncol=ncol(StudyFactorAndModel))
  for (i in 1:NumOfModel) {
    StudyFactorMat[,i]=as.data.frame(str_split(StudyFactorAndModel[,i],"-"))[1,] %>% as.vector %>% unlist %>% as.character
    StudyFactorModelMat[,i]=as.data.frame(str_split(StudyFactorAndModel[,i],"-"))[2,] %>% as.vector %>% unlist %>% as.character
  }
  colnames(StudyFactorMat)=colnames(StudyFactorModelMat)=str_c("model",1:NumOfModel)
  list(StudyFactorMat=StudyFactorMat,StudyFactorModelMat=StudyFactorModelMat)
}