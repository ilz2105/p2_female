#NEED TO install tidyverse, MASS, rsq, car
#v1: interaction among age-period-race is considered only when IncluPeriodEffect==T
#v2: further output, for each study factor, AIC, BIC, delta_AIC, delta_BIC
#v3: 1) fixed bug in Drop1StudyFactorThenCalcBICAIC---if IncluPeriodEffect==T; 2) can handle quintile models; 3) add Rsq; 4) add VIF

crc_mvm_v3=function(rf_data,seer_data,gender,
                    study_factors,study_factor_models,
                    age_grp_bound=c(5,9),
                    IncluPeriodEffect=T,
                    IncluAgePerInter=F,IncluRacPerInter=F,IncluAgeRacInter=F,P_THRESHOLD=0.05
){
  #misc prepare
  require(tidyverse)
  if(tolower(gender) %in% c("women","female")) TargetGender="Female" else if(tolower(gender) %in% c("men","male")) TargetGender="Male" else stop("Wrong gender value!")
  if(length(study_factors)!=length(study_factor_models)) stop("Study factors and corresponding models not the same length!")
  col_nam_PerAgeGenRac=c("period","age_grp","gender","race")
  col_nam_PerAgeRac=c("period","age_grp","race")
  std_fun=function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)
  
  ###########
  #get seer data for main and sensitivity analysis for one gender
  ###########
  seer_9=        seer_data %>% dplyr::filter(main.vs.sa=="main",gender==TargetGender) %>% dplyr::select(-gender) %>% as_tibble
  seer_9_13_18  =seer_data %>% dplyr::filter(main.vs.sa=="sa",gender==TargetGender) %>% dplyr::select(-gender) %>% as_tibble
  seer_ext=c("9","9_13_18")
  ###########
  #get targeted risk factor (RF) data for one gender
  ###########
  TargetRfData=rf_data %>% dplyr::select(all_of(col_nam_PerAgeGenRac),all_of(study_factors)) %>%
    dplyr::filter(gender==TargetGender) %>% dplyr::select(-gender) %>% as_tibble
  ########################
  #five regression models
  ########################
  #notes:
  #1) each fun below rtrn a dataframe containing cols of c("period","age_grp","race") and rf name
  #2) when age_grp_bound=c(5,9), the output$age_grp are in 5:9; when =c(10,14), the output$age_grp are in 10:14
  #3) output$period don't go beyong 3:10
  #4) within the above bounds, some rf contains NAs; or certain age-period-race combination (or row) may not exist
  fun_no_lag=function(data.agg,TargetRF,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    out=data.agg %>% dplyr::filter(age_grp %in% AgeForOutput,period %in% PeriodForOutput) %>%#age 20-49; exclude Nhanes1
      dplyr::select(all_of(col_nam_PerAgeRac),all_of(TargetRF))
    #build a function to get quintiles
    get_quintile=function(x){
      out=x
      non.na.id=!is.na(x)
      #consider quintiles or tertiles
      breaks_for_cut=quantile(x[non.na.id],probs=(0:5)/5) %>% {x=.;x[1]=x[1]-0.1;x}#deal with low boundary problem
      if(length(unique(breaks_for_cut))<6){#some quintile cut points are the same
        out[non.na.id]=cut(x[non.na.id],breaks=quantile(x[non.na.id],probs=(0:3)/3) %>% {x=.;x[1]=x[1]-0.1;x})
        return(factor(out,levels=paste0(1:3)))
      }else{#the six quintile cut points are all different
        out[non.na.id]=cut(x[non.na.id],breaks=breaks_for_cut)#assign fctr to numeric, out becomes numeric with 1:5, including NAs
        return(factor(out,levels=paste0(1:5)))
      }
    }
    out[,TargetRF]=get_quintile(unlist(out[,TargetRF]))
    out
  }
  fun_lag=function(data.agg,TargetRF,UseWithOthFun=F,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    out=data.agg %>% dplyr::select(all_of(col_nam_PerAgeRac),all_of(TargetRF))
    out[,TargetRF]=NA_real_
    
    #digress to deal with the below situation
    ###a rf don't have period 4 data (yet the rf has period 2 data); Thus we need a row for period 4 below to grab period 2 rf data;
    ###colitis have data for periods 2-6; Thus we need rows for period 7-8 below to grab period 5-6 data;
    period_miss=base::setdiff(2:10,unique(out$period))
    if(length(period_miss)>0){
      for (period_miss_i in period_miss) {
        out=bind_rows(out,
                      out %>% dplyr::filter(period==out$period[1]) %>% mutate(period=period_miss_i))
      }
      out=out %>% distinct#just be defensive
    }
    
    for (i in 1:nrow(out)) {
      #fix period, age, race group, 
      tmp_demo=out %>% dplyr::select(all_of(col_nam_PerAgeRac)) %>% slice(i)
      #find period-2 and age-2 but same race
      data_ear=data.agg %>% dplyr::filter(period==tmp_demo$period-2,age_grp==tmp_demo$age_grp-2,race==tmp_demo$race) %>% dplyr::select(all_of(TargetRF))
      if(nrow(data_ear)==1) out[i,TargetRF]=data_ear#sometimes, there are NAs
    }
    
    if(UseWithOthFun==F){
      out %>% dplyr::filter(age_grp %in% AgeForOutput,period %in% PeriodForOutput) %>% fun_no_lag(.,TargetRF,PeriodForOutput,AgeForOutput)
    }else{
      out# out will be given to amc_lag2; age_grp and period bounds will be considered in that fun;
    }
  }
  fun_amc_no_lag=function(data.agg,TargetRF,UseWithOthFun=F,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    #if UseWithOthFun==T, 1) don't standardize RF, 2) include age20-24 (agegrp4) and nhanes1 (period2)
    if(UseWithOthFun==F){
      out=data.agg %>% dplyr::filter(age_grp %in% AgeForOutput, period %in% PeriodForOutput) %>% dplyr::select(all_of(col_nam_PerAgeRac),all_of(TargetRF))
    }else{
      out=data.agg %>% dplyr::select(all_of(col_nam_PerAgeRac),all_of(TargetRF))
      # output will be given to amc_lag2 or amc_cum10y5ylag2; age_grp and period bounds will be considered in that fun
    }
    
    #group by age_grp and race
    #get mean
    age_race_mean=out %>% dplyr::select(-period) %>% group_by(age_grp,race) %>% summarise_all(mean,na.rm=T)
    #match the mean to out
    tmp=left_join(x=out,y=age_race_mean,by=c("age_grp"="age_grp","race"="race"))
    #minus mean for each col
    tmp.x=tmp %>% dplyr::select(contains(".x"))#don't contain demographic varaibles
    tmp.y=tmp %>% dplyr::select(contains(".y"))
    tmp=tmp.x-tmp.y#overwrite!
    if(UseWithOthFun==F) tmp=apply(tmp,2,std_fun)#standardize and overwrite!
    colnames(tmp)=TargetRF#b/c I changed colname earliers
    
    bind_cols(out[,col_nam_PerAgeRac],as_tibble(tmp))
  }
  fun_amc_lag2=function(data.agg,TargetRF,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){ # do amc before lag
    out=fun_lag(fun_amc_no_lag(data.agg,TargetRF,UseWithOthFun=T),TargetRF,UseWithOthFun=T)
    out %>% dplyr::filter(age_grp %in% AgeForOutput, period %in% PeriodForOutput) %>% mutate(across(.cols=all_of(TargetRF),.fns=std_fun))
  }
  fun_amc_cum10y5ylag2=function(data.agg,TargetRF,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    tmp=data.agg
    #fill in period 4 (linear interpolation)
    #identify if needed; if so, fill in; do this for all rf
    tmp_tmp=tmp %>% dplyr::select(all_of(col_nam_PerAgeRac),all_of(TargetRF))
    tmp_tmp=tmp_tmp[!is.na(unlist(tmp_tmp[,TargetRF])),]#overwrite and keep rows not NA
    period_vec=unique(tmp_tmp$period)#periods with data for rf i
    if(all( (c(3,4,5) %in% period_vec)==c(T,F,T) )){
      tmp_tmp_pad_4=tmp_tmp %>% dplyr::filter(period %in% c(3,5)) %>% dplyr::select(-period) %>% 
        group_by(age_grp,race) %>% summarize_all(mean,na.rm=F) %>% dplyr::mutate(period=4)#na.rm doesn't matter b/c no NA
      tmp_tmp=bind_rows(tmp_tmp,tmp_tmp_pad_4)#overwrite!
    }
    tmp=tmp_tmp#overwrite! tmp can't be overwriten while for-loop was run
    out=tmp=fun_amc_no_lag(tmp,TargetRF,UseWithOthFun=T)
    out[,TargetRF]=NA_real_
    
    #digress to deal with the below situation
    ###colitis have data for periods 2-6; Thus we need rows for period 7 below to grab period 5-6 data;
    period_miss=base::setdiff(2:10,unique(out$period))
    if(length(period_miss)>0){
      for (period_miss_i in period_miss) {
        out=bind_rows(out,
                      out %>% dplyr::filter(period==out$period[1]) %>% mutate(period=period_miss_i))
      }
      out=out %>% distinct
    }
    
    for (i in 1:nrow(out)) {
      #fix period, age, race group, 
      tmp_demo=out %>% dplyr::select(all_of(col_nam_PerAgeRac)) %>% slice(i)
      # 5 yr lag, 10 y cum
      data_cum=tmp %>% dplyr::filter(
        (period==tmp_demo$period-1 & age_grp==tmp_demo$age_grp-1) | 
          (period==tmp_demo$period-2 & age_grp==tmp_demo$age_grp-2), race==tmp_demo$race) %>% dplyr::select(all_of(TargetRF))
      if(nrow(data_cum)==2)  out[i,TargetRF]=colSums(data_cum,na.rm=F)#sometimes, there are NAs
    }
    out %>% dplyr::filter(age_grp %in% AgeForOutput, period %in% PeriodForOutput) %>% mutate(across(.cols=all_of(TargetRF),.fns=std_fun))
  }
  
  ##########################################################
  #create a function for dummy coding age, period, and race#
  ##########################################################
  dummy_code_AgePerRac=function(df){
    periods=sort(unique(df$period))[-1]#no dummy coding for reference level
    age_grps=sort(unique(df$age_grp))[-1]
    for (i_periods in periods) {
      #df$pi_periods=if_else(df$period==i_periods,'1','0') %>% factor
      eval(parse(text=str_c("df$p",i_periods,"=if_else(df$period==i_periods,'1','0') %>% factor")))
    }
    for (i_age_grps in age_grps) {
      #df$ai_age_grps=if_else(df$age_grp==i_age_grps,'1','0')
      eval(parse(text=str_c("df$a",i_age_grps,"=if_else(df$age_grp==i_age_grps,'1','0') %>% factor")))
    }
    df$rBlack=if_else(df$race=="Black",'1','0') %>% factor
    
    AllAgeTerm=str_c(str_c("a",age_grps),collapse="+")#works even if length(age_grps) is 1
    AllPerTerm=str_c(str_c("p",periods ),collapse="+")
    AllAgePerRacTerm=str_c(AllAgeTerm,"+",AllPerTerm,"+rBlack")
    attr(df,"AllAgePerRacTerm")=AllAgePerRacTerm
    attr(df,"AllAgeTerm")=AllAgeTerm
    attr(df,"AllPerTerm")=AllPerTerm
    attr(df,"AllRacTerm")="rBlack"
    df
  }
  
  ########################################################
  #create a function to get interaction terms for testing#
  ########################################################
  get_inter_term=function(df,IncluAgePerInter,IncluRacPerInter,IncluAgeRacInter){
    if(IncluAgePerInter==T & IncluRacPerInter==F & IncluAgeRacInter==F){
      InterTerm=str_c("(",attributes(df)$AllAgeTerm,")*(",attributes(df)$AllPerTerm,")")
    }else if(IncluAgePerInter==T & IncluRacPerInter==T & IncluAgeRacInter==F){
      InterTerm=str_c("(",attributes(df)$AllAgeTerm,")*(",attributes(df)$AllPerTerm,")+(",
                      attributes(df)$AllPerTerm,")*rBlack")
    }else if(IncluAgePerInter==T & IncluRacPerInter==T & IncluAgeRacInter==T){
      InterTerm=str_c("(",attributes(df)$AllAgeTerm,")*(",attributes(df)$AllPerTerm,")+(",
                      attributes(df)$AllPerTerm,")*rBlack+(",attributes(df)$AllAgeTerm,")*rBlack")
    }else if(IncluAgePerInter==F & IncluRacPerInter==F & IncluAgeRacInter==F){
      InterTerm=NULL
    }else{
      stop("Wrong value for IncluAgePerInter, IncluRacPerInter, or IncluAgeRacInter!")
    }
  }
  
  ###############################
  #create a function to calculate GOF, AIC&BIC&Rsq when one study factor is dropped; this allows me to get delta_aic/bic/rsq later
  ###############################
  Drop1StudyFactorThenCalcGOF=function(data,study_factor_model_nameS,IncluPeriodEffect,FullRegCoefNum){
    NumOfStudyFactor=length(study_factor_model_nameS)
    AIC_vec=BIC_vec=Rsq_vec=list()
    for (i in 1:NumOfStudyFactor) {
      if(IncluPeriodEffect==T){
        suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(DataForReg2)$AllAgePerRacTerm,SigInterTerm,
                                                "+",str_c(study_factor_model_nameS[-i],collapse="+"),"+offset(log(pop)),data=DataForReg2)"))))
      }else{
        suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(DataForReg2)$AllAgeTerm,"+",attributes(DataForReg2)$AllRacTerm,
                                                "+",str_c(study_factor_model_nameS[-i],collapse="+"),"+offset(log(pop)),data=DataForReg2)"))))
      }
      irr_out=exp(coef(res_tmp))
      sam_siz=nrow(DataForReg2)
      k_model=(length(irr_out)+1)#+1 free param for error
      RepNum=FullRegCoefNum-(k_model-1)#whether the study factor is amc or quintile or tertile
      BIC_vec[[i]]=rep(k_model*log(sam_siz)-summary(res_tmp)$twologlik,RepNum)
      AIC_vec[[i]]=rep(2*k_model-summary(res_tmp)$twologlik,RepNum)
      Rsq_vec[[i]]=rep(rsq::rsq(res_tmp,adj=T),RepNum)
    }
    list(AIC_vec=unlist(AIC_vec),BIC_vec=unlist(BIC_vec),Rsq_vec=unlist(Rsq_vec))
  }
  
  ######################
  #create a function to calculate vif for study factors, indicator for multi collinearity
  ######################
  CalcVIFs=function(data,study_factor_model_nameS,IncluPeriodEffect){
    if(IncluPeriodEffect==T){
      suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(DataForReg2)$AllAgePerRacTerm,SigInterTerm,
                                              "+",str_c(study_factor_model_nameS,collapse="+"),"+offset(log(pop)),data=DataForReg2)"))))
    }else{
      suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(DataForReg2)$AllAgeTerm,"+",attributes(DataForReg2)$AllRacTerm,
                                              "+",str_c(study_factor_model_nameS,collapse="+"),"+offset(log(pop)),data=DataForReg2)"))))
    }
    #vif_tmp=car::vif(res_tmp)#matrix or vector
    vif_tmp=try(car::vif(res_tmp))
    if(class(vif_tmp)=="try-error"){
      writeLines(geterrmessage())
      NA_real_
    }else{
      VIF_vec=list()
      if(is.matrix(vif_tmp)==F){#i.e., no quintile models included
        VIF_vec[[1]]=vif_tmp[grepl("lag",names(vif_tmp))]
      }else{
        for (i in 1:length(study_factor_model_nameS)) {
          VIF_vec[[i]]=rep(vif_tmp[study_factor_model_nameS[i],"GVIF^(1/(2*Df))"],vif_tmp[study_factor_model_nameS[i],"Df"])
        }
      }
      unlist(VIF_vec)
    }
  }
  
  ####################
  #apply diff models, get a combined dataset with all RFs, change RF names to include model name, drop rows with NA
  ####################
  #first get final output period/age_grp (available) for each study_factor-model combo
  #get the intersect of the above vectors (of period/age_grp), which is the final output period/age_grp (available) when joining all study_factor-model combo 
  NumOfStudyFactor=length(study_factors)
  for (i in 1:NumOfStudyFactor) {
    Fn_tmp=get(str_c("fun_",study_factor_models[i]))
    OneRfData_tmp=Fn_tmp(TargetRfData,study_factors[i]) %>% drop_na#too expensive here, improve later, perhaps
    PeriodForOutput_tmp=unique(OneRfData_tmp$period)
    AgeForOutput_tmp=unique(OneRfData_tmp$age_grp)
    if(i==1){
      PeriodForOutput=PeriodForOutput_tmp
      AgeForOutput=AgeForOutput_tmp
    }else{
      PeriodForOutput=base::intersect(PeriodForOutput,PeriodForOutput_tmp)
      AgeForOutput=base::intersect(AgeForOutput,AgeForOutput_tmp)
    }
  }
  
  study_factor_model_nameS=character(NumOfStudyFactor)
  for (i in 1:NumOfStudyFactor) {
    Fn_tmp=get(str_c("fun_",study_factor_models[i]))
    study_factor_model_nameS[i]=str_c(study_factors[i],"_",study_factor_models[i])
    #OneRfDatai=Fn_tmp(TargetRfData,study_factors[i],PeriodForOutput=PeriodForOutput,AgeForOutput=AgeForOutput)%>%rename(study_factor_model_nameS[i]=study_factors[i])
    eval(parse(text=str_c("OneRfData",i,"=Fn_tmp(TargetRfData,study_factors[i],PeriodForOutput=PeriodForOutput,AgeForOutput=AgeForOutput) %>% 
                          rename(",study_factor_model_nameS[i],"=",study_factors[i],")")))
    if(i==1){
      TargetRfDataUse=OneRfData1
    } else{
      #TargetRfDataUse=inner_join(TargetRfDataUse,OneRfDatai)
      suppressMessages(eval(parse(text=str_c("TargetRfDataUse=inner_join(TargetRfDataUse,OneRfData",i,")"))))
    }
  }
  TargetRfDataUse=TargetRfDataUse %>% drop_na
  
  #############
  #analysis for seer 9 and 9-18
  #############
  out_studyfactor=out=NULL
  for (seer_version in seer_ext) {
    #DataForReg1=suppressMessages(inner_join(seer_seer_version,TargetRfDataUse))
    eval(parse(text=str_c("DataForReg1=suppressMessages(inner_join(seer_",seer_version,",TargetRfDataUse))")))#join study factor data and seer data
    #dummy coding for age, period, and race
    DataForReg2=dummy_code_AgePerRac(DataForReg1)
    
    if(IncluPeriodEffect==T){
      #if do interaction: do regression to test significant interaction terms
      InterTerm=get_inter_term(DataForReg2,IncluAgePerInter,IncluRacPerInter,IncluAgeRacInter)#get all interaction terms for testing
      if(any(c(IncluAgePerInter,IncluAgeRacInter,IncluRacPerInter))==T){
        suppressWarnings(eval(parse(text=str_c("res_test_inter=MASS::glm.nb(count~",attributes(DataForReg2)$AllAgePerRacTerm,"+",InterTerm,
                                               "+offset(log(pop)),data=DataForReg2)"))))
        #select significant interaction terms
        #grab the matrix holding param names (as rownames) and p-values, etc
        SigInterTerm=summary(res_test_inter)$coefficients %>% {x=.;as_tibble(x) %>% mutate(param_name=rownames(x))} %>%
          filter(`Pr(>|z|)`<P_THRESHOLD) %>% filter(grepl(":",param_name)) %>% .$param_name %>%
          str_replace(.,"1:","*") %>% str_sub(.,1,-2) %>% str_c(.,collapse="+") %>% str_c("+",.)
      }else SigInterTerm=NULL
      
      #regression
      suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(DataForReg2)$AllAgePerRacTerm,
                                              SigInterTerm,
                                              "+",str_c(study_factor_model_nameS,collapse="+"),"+offset(log(pop)),data=DataForReg2)"))))
    }else{#IncluPeriodEffect==F
      #regression
      suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(DataForReg2)$AllAgeTerm,"+",attributes(DataForReg2)$AllRacTerm,
                                              "+",str_c(study_factor_model_nameS,collapse="+"),"+offset(log(pop)),data=DataForReg2)"))))
    }
    
    #output
    irr_out=exp(coef(res_tmp))
    irr_ci_out=suppressMessages(exp(confint(res_tmp)))
    irr_p_out=summary(res_tmp)$coefficients[,4]
    irr_p_out=if(length(irr_p_out)<length(irr_out)) NA_real_ else irr_p_out#when one pvalue is NA, the model is too saturated to be useful, so output NA for all pvalue
    
    GOF_vec=Drop1StudyFactorThenCalcGOF(data=DataForReg2,study_factor_model_nameS=study_factor_model_nameS,IncluPeriodEffect=IncluPeriodEffect,
                                        FullRegCoefNum=length(irr_out))
    VIF_vec=CalcVIFs(data=DataForReg2,study_factor_model_nameS=study_factor_model_nameS,IncluPeriodEffect=IncluPeriodEffect)
    
    sam_siz=nrow(DataForReg2)
    k_model=(length(irr_out)+1)#+1 free param for error
    BIC=k_model*log(sam_siz)-summary(res_tmp)$twologlik
    AIC=2*k_model-summary(res_tmp)$twologlik
    Rsq=rsq::rsq(res_tmp,adj=T)
    
    #periods_available=str_c(as.character(sort(unique(DataForReg2$period))),collapse="+")
    #num_periods_available=length(unique(DataForReg2$period))
    periods_available=str_c(PeriodForOutput,collapse="+")
    num_periods_available=length(PeriodForOutput)
    
    if(IncluPeriodEffect==T){
      if(any(c(IncluAgePerInter,IncluAgeRacInter,IncluRacPerInter))==F) SigInterTerm=" None"#overwrite NULL for output
    }else{
      SigInterTerm=" None"
    }
    out_tmp=tibble(gender=gender,seer=seer_version,covariate_name=names(irr_out),
                   IRR=irr_out,IRR_low=irr_ci_out[,1],IRR_up=irr_ci_out[,2],`p-value`=irr_p_out,BIC=BIC,AIC=AIC,Rsq=Rsq,
                   `Periods available`=periods_available,`Num of periods available`=num_periods_available,`Interaction term`=str_sub(SigInterTerm,2,-1))
    out_studyfactor_tmp=out_tmp %>% dplyr::filter(grepl("lag",covariate_name)) %>%#b/c all study factor name contain "lag"!!
      mutate(delta_AIC=GOF_vec$AIC_vec-AIC,delta_BIC=GOF_vec$BIC_vec-BIC,delta_Rsq=Rsq-GOF_vec$Rsq_vec,VIF=VIF_vec)
    
    out=bind_rows(out,out_tmp)
    out_studyfactor=bind_rows(out_studyfactor,out_studyfactor_tmp)
  }#end of for (seer_version in seer_ext)
  list(FullResult=out,
       StudyFactorResult=out_studyfactor %>% dplyr::select(gender,seer,covariate_name,IRR,IRR_low,IRR_up,`p-value`,contains("delta"),VIF,everything()))
}



