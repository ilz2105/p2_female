library(tidyverse)
dir_resu=dir_data=dir_code="C:/Users/Lulu/OneDrive - cumc.columbia.edu/Documents/crc/git_cancer_internal/code_lz/p2 rf/p2_female_analyses/"
TargetGender="Female"
VersionTag="v1"
source(str_c(dir_code,"crc_mvm_v3.R"))
source(str_c(dir_code,"utils.R"))
seer_data=read_csv(str_c(dir_data,"crc_2074.csv")) %>% dplyr::select(period,age_grp,race,gender,count,pop,main.vs.sa)
load(str_c(dir_data,"p2_rf_list.Rdata"))
rf_data=p2_rf_list %>% as_tibble# %>% select(period,age_grp,race,gender,cs_m,adj_avg_alc_m,overweight_up_m,obese_m,thy_m,dia_m.brfss,hbp_m.brfss,calcium_m)
IncluPeriodEffect=F
gender=TargetGender#Male/Female
if(VersionTag=="v1"){
  StudyFactorModelCombo=list(c("cs_m-amc_cum10y5ylag2"),
                             c("adj_avg_alc_m-amc_cum10y5ylag2"),
                             c("overweight_up_m-no_lag","obese_m-no_lag",
                               "overweight_up_m-amc_no_lag","obese_m-amc_no_lag",
                               "overweight_up_m-lag","obese_m-lag",
                               "overweight_up_m-amc_lag2","obese_m-amc_lag2"),
                             c("hbp_m.brfss-no_lag","hbp_m.brfss-amc_no_lag","hbp_m.brfss-lag","hbp_m.brfss-amc_lag2"),
                             c("dia_m.brfss-no_lag","dia_m.brfss-amc_no_lag","dia_m.brfss-lag","dia_m.brfss-amc_lag2"),
                             c("thy_m-no_lag","thy_m-amc_no_lag","thy_m-lag","thy_m-amc_lag2"),
                             c("calcium_m-lag","calcium_m-amc_lag2","calcium_m-amc_cum10y5ylag2"))
}else if(VersionTag=="v2"){
  StudyFactorModelCombo=list(c("cs_m-lag","cs_m-amc_lag2","cs_m-amc_cum10y5ylag2"),
                             c("adj_avg_alc_m-lag","adj_avg_alc_m-amc_lag2","adj_avg_alc_m-amc_cum10y5ylag2"),
                             c("overweight_up_m-no_lag","obese_m-no_lag",
                               "overweight_up_m-amc_no_lag","obese_m-amc_no_lag",
                               "overweight_up_m-lag","obese_m-lag",
                               "overweight_up_m-amc_lag2","obese_m-amc_lag2"),
                             c("hbp_m.brfss-no_lag","hbp_m.brfss-amc_no_lag","hbp_m.brfss-lag","hbp_m.brfss-amc_lag2"),
                             c("dia_m.brfss-no_lag","dia_m.brfss-amc_no_lag","dia_m.brfss-lag","dia_m.brfss-amc_lag2"),
                             c("thy_m-no_lag","thy_m-amc_no_lag","thy_m-lag","thy_m-amc_lag2"),
                             c("calcium_m-lag","calcium_m-amc_lag2","calcium_m-amc_cum10y5ylag2"))
}
StudyFactorMat=CreateStudyFactorAndModel(StudyFactorModelCombo,gender)$StudyFactorMat
StudyFactorModelMat=CreateStudyFactorAndModel(StudyFactorModelCombo,gender)$StudyFactorModelMat

models=colnames(StudyFactorMat)
NumOfModel=length(models)
OutCombo=NULL
for (model_i in models) {
  out=crc_mvm_v3(rf_data=rf_data,seer_data=seer_data,gender=gender,
                 study_factors=StudyFactorMat[,model_i],
                 study_factor_models=StudyFactorModelMat[,model_i],
                 IncluPeriodEffect=IncluPeriodEffect)
  OutCombo=bind_rows(OutCombo,out$StudyFactorResult %>% mutate(model=model_i))
  if(as.numeric(str_sub(model_i,6,-1)) %% 50==0) print(str_c("Finished ",model_i," for ",gender,"; ",NumOfModel," in total."))
}
#OutCombo_v1_male=OutCombo
eval(parse(text=str_c("OutCombo_",VersionTag,"_",gender,"=OutCombo")))
#save(OutCombo_v1_male,file="OutComb_v1_male.RData")
eval(parse(text=str_c("save(OutCombo_",VersionTag,"_",gender,",file='",dir_resu,"OutComb_",VersionTag,"_",gender,".RData')")))
