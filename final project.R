---
title: "Final Project"
author: "Yukun Shen"
output: word_document
---

setwd("/Users/shizhaoli/desktop/R/")
install.packages("ivprobit")
library(stargazer)
library(gdata)
library(ggplot2)
library(psych) 
library(ggeffects)
library(QuantPsyc)
library(usdm)
library(lmtest)
library(multiwayvcov)
library(sandwich)
library(foreign)
library(AER)
library(aod)
library(Rcpp)
library(mfx)
library(nnet)
library(reshape2)
library(VIF)
library(MASS)
library(readstata13)
library(ivprobit)
#==========================================================
## FIRST Q (OLS REGRESSION )
#==========================================================
```{r}
data1 = read.dta13("/Users/yukunshen/Desktop/BOPS/consumer level data.dta")
data2 = read.dta13("/Users/yukunshen/Desktop/BOPS/online daily prod_cat sales-returns data.dta")
data3 = read.dta13("/Users/yukunshen/Desktop/BOPS/online daily sales-returns data.dta")
data4 = read.dta13("/Users/yukunshen/Desktop/BOPS/transaction level data.dta")
stargazer(data3, type="text", median=TRUE, iqr=TRUE,digits=1, title="Descriptive Statistics")  
ggplot(data3, aes(x=salesvalue)) + geom_histogram(colour="green")
ggplot(data3, aes(x=log(salesvalue))) + geom_histogram(colour="green")
summary(data3)
mydata1 <- subset(data3,day<786)#generate a new data frame that exclude days beyond 783
mydata1$group<-ifelse(mydata1$store_number==5998,0,1)#seperate stores by creating a dummy variable
mydata1$time<-ifelse(mydata1$day<366,0,1)#creating a time dummy variable
mydata1$avg_femalemean<-ifelse(is.na(mydata1$avg_female),mean(mydata1$avg_female,na.rm=TRUE),mydata1$avg_female)
mydata1$avg_incomemean<-ifelse(is.na(mydata1$avg_income),mean(mydata1$avg_income,na.rm=TRUE),mydata1$avg_income)
mydata1$avg_homeownermean<-ifelse(is.na(mydata1$avg_homeowner),mean(mydata1$avg_homeowner,na.rm=TRUE),mydata1$avg_homeowner)
mydata1$avg_residencymean<-ifelse(is.na(mydata1$avg_residency),mean(mydata1$avg_residency,na.rm=TRUE),mydata1$avg_residency)
mydata1$avg_childownermean<-ifelse(is.na(mydata1$avg_childowner),mean(mydata1$avg_childowner,na.rm=TRUE),mydata1$avg_childowner)
mydata1$avg_agemean<-ifelse(is.na(mydata1$avg_age),mean(mydata1$avg_age,na.rm=TRUE),mydata1$avg_age)
df<-mydata1[c("time","avg_femalemean","avg_incomemean","avg_homeownermean","avg_residencymean","avg_childownermean","returnquantity","avg_agemean","group","returnvalue","returnquantity")]
cor(df) 
vifcor(df)
m1 <- lm(log(salesvalue+1)~time*group,data=mydata1)
stargazer(m1, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
#add control variables
is.factor(mydata1$month_dummy)
mydata1$month_dummy <- as.factor(mydata1$month_dummy)
m2 <- lm(log(salesvalue+1)~time*group+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=mydata1)
stargazer(m2, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
meffects1 <- ggpredict(m2, terms=c("time", "group")) # generates a tidy data frame  

ggplot(meffects1,aes(x, predicted, colour=group)) + geom_line(size=1.3) + 
  xlab("time") + ylab("predicted sales")

# Check for heteroscedasticity
gqtest(m2)
bptest(m2)
HWrobstder <- sqrt(diag(vcovHC(m2, type="HC1")))
stargazer(m2,m2,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))
#==========================================================
## Q1:POISSON REGRESSION & Negative Binomial
#==========================================================
poisson1 <- glm(salesquantity~time*group,family="poisson",data=mydata1)
stargazer(poisson1,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson1a <- glm(salesquantity~1, data=mydata1, family="poisson") # This is the command to run a logit on null model 
lrtest(poisson1, poisson1a)
#add control variables
poisson2 <- glm(salesquantity~time*group+avg_agemean+avg_femalemean+avg_incomemean+month_dummy+avg_childownermean+avg_homeownermean,family="poisson",data=mydata1)
stargazer(poisson1,  poisson2,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
lrtest(poisson2, poisson1a)

#negative binomial model
negbin1 <- glm.nb(salesquantity~time*group,data=mydata1)

stargazer(negbin1,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#add controp varaiables
negbin2 <- glm.nb(salesquantity~time*group+avg_femalemean+avg_incomemean+avg_agemean+month_dummy+avg_childownermean+avg_homeownermean,data=mydata1)

stargazer(negbin1,  negbin2,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#check for heterokadasticity
gqtest(negbin2) 
bptest(negbin2)
HWrobstder <- sqrt(diag(vcovHC(negbin2, type="HC1"))) # produces Huber-White robust standard errors 

stargazer(negbin2,negbin2,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))  # 
negbin1a <- glm.nb(log(salesquantity) ~ 1, data = mydata1) 

lrtest(negbin2, negbin1a) #significant P-value means model fit
#choose model
lrtest(poisson2,negbin2)
#==========================================================
##Second Question (OLS)
#==========================================================

ggplot(data3, aes(x=returnvalue)) + geom_histogram(colour="green")
ggplot(data3, aes(x=log(returnvalue))) + geom_histogram(colour="green")
m21 <- lm(log(returnvalue+1)~time*group,data=mydata1)
stargazer(m21, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
m22 <- lm(log(returnvalue+1)~time*group+avg_agemean+avg_femalemean+avg_incomemean+month_dummy+salesvalue+avg_childownermean+avg_homeownermean,data=mydata1)
stargazer(m21, m22,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))     
#check which model is better
anova(m21,m22)
# Check for heteroscedasticity
gqtest(m22)
bptest(m22)
HWrobstde <- sqrt(diag(vcovHC(m22, type="HC1")))
stargazer(m22,m22,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))
#==========================================================
##Second Question (POISSON REGRESSION & Negative Binomial)
#==========================================================

poisson21 <- glm(returnquantity~time*group,family="poisson",data=mydata1)
stargazer(poisson21,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson21a <- glm(returnquantity~1, data=mydata1, family="poisson") # This is the command to run a logit on null model 
lrtest(poisson1, poisson1a)
poisson22 <- glm(returnquantity~time*group+avg_agemean+avg_incomemean+avg_femalemean+salesquantity+month_dummy+avg_childownermean,family="poisson",data=mydata1)
stargazer(poisson21, poisson22,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
anova(poisson21,poisson22,test="Chisq")#indicates poisson 22 is better
poisson22a <- glm(returnquantity~1, data=mydata1, family="poisson") # This is the command to run a logit on null model 
lrtest(poisson22, poisson22a)

#negative binomial
negbin21 <- glm.nb(returnquantity~time*group,data=mydata1)

stargazer(negbin21,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))

negbin22 <- glm.nb(returnquantity~time*group+avg_femalemean+avg_incomemean+salesquantity+month_dummy+avg_agemean+avg_childownermean,data=mydata1)

stargazer(negbin21,  negbin22,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
anova(negbin21,negbin22,test="Chisq")
#generate an empty model
negbin21a <- glm.nb(returnquantity ~ 1, data = mydata1) 
lrtest(negbin22, negbin21a)  #P-value significant means model fit
#choose which model is better
lrtest(poisson22,negbin22)  #significant P-value indicates negtive binomial model is better


#==========================================================
##Third Question
#==========================================================

summary(data1)
ggplot(data1, aes(x=salesvalue)) + geom_histogram(colour="black")
ggplot(data1, aes(x=log(salesvalue))) + geom_histogram(colour="black")

#drop all NA values
mydata2 <- na.omit(data1)

#change child & homeowner_code into dummy variable
mydata2$childdummy <- ifelse(mydata2$child=="Y",1,0)
mydata2$homeowner_codedummy <- ifelse(mydata2$homeowner_code=="O",1,0)
#check multicollinearity
df3<-mydata2[c("bops_in_effect","bops_user","est_income_code","age_band","purchase_time_period","female","childdummy")]
cor(df3)


#generate initial model
m31 <- lm(log(salesvalue+1)~bops_user*bops_in_effect,data=mydata2)
stargazer(m31, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
m32 <- lm(log(salesvalue+1)~bops_user*bops_in_effect+est_income_code+age_band+childdummy+homeowner_codedummy+female,data=mydata2)
stargazer(m31, m32,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 
# Check for heteroscedasticity
gqtest(m32)
bptest(m32)
HWrobstder <- sqrt(diag(vcovHC(m32, type="HC1")))
stargazer(m32,m32,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))

#==========================================================
##Third Question Poisson & Negative Binomial
#==========================================================
poisson31 <- glm(salesquantity~bops_user*bops_in_effect,family="poisson",data=mydata2)
stargazer(poisson31,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson31a <- glm(salesquantity~1, data=mydata2, family="poisson") # This is the command to run a logit on null model 
lrtest(poisson1, poisson1a)
#add control variables
poisson32 <- glm(salesquantity~bops_user*bops_in_effect+est_income_code+age_band+female+childdummy+female,family="poisson",data=mydata2)
stargazer(poisson31,  poisson32,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
lrtest(poisson32, poisson31a)

#negative binomial model
negbin31 <- glm.nb(salesquantity~bops_user*bops_in_effect,data=mydata2)

stargazer(negbin31,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#add controp varaiables
negbin32 <- glm.nb(salesquantity~bops_user*bops_in_effect+est_income_code+age_band+female+childdummy+female,data=mydata2)

stargazer(negbin31,  negbin32,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))

# Check for heteroscedasticity
gqtest(negbin32)
bptest(negbin32)
HWrobstder <- sqrt(diag(vcovHC(negbin32, type="HC1")))
stargazer(negbin32,negbin32,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

negbin31a <- glm.nb(log(salesquantity) ~ 1, data = mydata2) #generate a perfect model for comparison

lrtest(negbin32, negbin31a) #significant P-value means model fit
#choose model
lrtest(poisson32,negbin32) #significant P-Value indicates negative binomial is better

#==========================================================
##Fourth Question
#==========================================================

mydata4 <- data4[!is.na(data4$bops), ] #exclude all NA variables in bops
mydata4$childdummy <- ifelse(mydata4$child=="Y",1,0)
mydata4$homeowner_codedummy <- ifelse(mydata4$homeowner_code=="o",1,0)
mydata4$est_income_codemedian<-ifelse(is.na(mydata4$est_income_code),median(mydata4$est_income_code,na.rm=TRUE),mydata4$est_income_code)
mydata4$age_bandmedian<-ifelse(is.na(mydata4$age_band),median(mydata4$age_band,na.rm=TRUE),mydata4$age_band)
mydata4$length_of_residencemedian<-ifelse(is.na(mydata4$length_of_residence),median(mydata4$length_of_residence,na.rm=TRUE),mydata4$length_of_residence)
mydata41 <- mydata4[!is.na(mydata4$female),]
mydata42<-mydata41[!mydata41$child=="", ]
m41 <- lm(return~bops+price+age_band+est_income_code+female,data=data4)
stargazer(m41, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
m41a <- lm(return~bops+price+age_bandmedian+est_income_codemedian+female,data=mydata42)
stargazer(m41a, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
logit41<- glm(return~bops+price+age_bandmedian+est_income_codemedian+female,data = mydata42, family="binomial")
stargazer(logit41, 
          title="Regression Results", type="text", 
          column.labels=c("Logit-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))

stargazer(logit41, 
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("OddsRatios"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))


#==========================================================
##Fivth Question
#==========================================================

#SALESVALUE
stargazer(data2, type="text", median=TRUE, iqr=TRUE,digits=1, title="Descriptive Statistics")  
ggplot(data2, aes(x=salesvalue)) + geom_histogram(colour="green")
ggplot(data2, aes(x=log(salesvalue))) + geom_histogram(colour="green")
summary(data2)
mydata5 <- subset(data2,day<786)#generate a new data frame that exclude days beyond 783
mydata5$group<-ifelse(mydata5$store_number==5998,0,1)#seperate stores by creating a dummy variable
mydata5$time<-ifelse(mydata5$day<366,0,1)#creating a time dummy variable
mydata5$avg_femalemean<-ifelse(is.na(mydata5$avg_female),mean(mydata5$avg_female,na.rm=TRUE),mydata5$avg_female)
mydata5$avg_incomemean<-ifelse(is.na(mydata5$avg_income),mean(mydata5$avg_income,na.rm=TRUE),mydata5$avg_income)
mydata5$avg_homeownermean<-ifelse(is.na(mydata5$avg_homeowner),mean(mydata5$avg_homeowner,na.rm=TRUE),mydata5$avg_homeowner)
mydata5$avg_residencymean<-ifelse(is.na(mydata5$avg_residency),mean(mydata5$avg_residency,na.rm=TRUE),mydata5$avg_residency)
mydata5$avg_childownermean<-ifelse(is.na(mydata5$avg_childowner),mean(mydata5$avg_childowner,na.rm=TRUE),mydata5$avg_childowner)
mydata5$avg_agemean<-ifelse(is.na(mydata5$avg_age),mean(mydata5$avg_age,na.rm=TRUE),mydata5$avg_age)
df5<-mydata5[c("time","avg_femalemean","avg_incomemean","avg_homeownermean","avg_residencymean","avg_childownermean","avg_agemean","group","returnvalue","returnquantity","product_category")]
cor(df5) 
vifcor(df5)
m51 <- lm(log(salesvalue+1)~time*group,data=mydata5)
stargazer(m51, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
#add control variables
is.factor(mydata5$month_dummy)
mydata5$month_dummy <- as.factor(mydata5$month_dummy)
m52 <- lm(log(salesvalue+1)~time*group+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=mydata5)
stargazer(m51,m52, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 

#Visiualization
meffects1 <- ggpredict(m52, terms=c("time", "group")) # generates a tidy data frame  

ggplot(meffects1,aes(x, predicted, colour=group)) + geom_line(size=1.3) + 
  xlab("time") + ylab("predicted sales")

# Check for heteroscedasticity
gqtest(m52)
bptest(m52)
HWrobstder <- sqrt(diag(vcovHC(m52, type="HC1")))
stargazer(m52,m52,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#=========================================================
##Fivth Question RETURNVALUE
#==========================================================
ggplot(mydata5, aes(x=returnvalue)) + geom_histogram(colour="green")
ggplot(mydata5, aes(x=log(returnvalue))) + geom_histogram(colour="green")

#Build initial model
m51r <- lm(log(returnvalue+1)~time*group,data=mydata5)
stargazer(m51r, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 

#add control variables
m52r <- lm(log(returnvalue+1)~time*group+salesvalue+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=mydata5)
stargazer(m51r,m52r, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
#Check for heterokedasticity
gqtest(m52r)
bptest(m52r)
HWrobstder <- sqrt(diag(vcovHC(m52r, type="HC1")))
stargazer(m52r,m52r,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))
#Visiualization
meffects1 <- ggpredict(m52r, terms=c("time", "group")) # generates a tidy data frame  

ggplot(meffects1,aes(x, predicted, colour=group)) + geom_line(size=1.3) + 
  xlab("time") + ylab("predicted sales")


#==========================================================
##Fivth Question SALESQUANTITY
#==========================================================
poisson51 <- glm(salesquantity~time*group,family="poisson",data=mydata5)
stargazer(poisson51,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson51a <- glm(salesquantity~1, data=mydata5, family="poisson") # This is the command to run a logit on null model 
lrtest(poisson51r, poisson51a)
#add control variables
poisson52 <- glm(salesquantity~time*group+avg_agemean+avg_femalemean+avg_incomemean+month_dummy+avg_childownermean+avg_homeownermean,family="poisson",data=mydata5)
stargazer(poisson51,  poisson52,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
lrtest(poisson52, poisson51a)

#negative binomial model
negbin51 <- glm.nb(salesquantity~time*group,data=mydata5)

stargazer(negbin51,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#add controp varaiables
negbin52 <- glm.nb(salesquantity~time*group+avg_femalemean+avg_incomemean+avg_agemean+month_dummy+avg_childownermean+avg_homeownermean,data=mydata5)

stargazer(negbin51,  negbin52,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#check for heterokadasticity
gqtest(negbin52) 
bptest(negbin52)
HWrobstder <- sqrt(diag(vcovHC(negbin52, type="HC1"))) # produces Huber-White robust standard errors 

stargazer(negbin52,negbin52,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))  # 
negbin51a <- glm.nb(salesquantity ~ 1, data = mydata5) #Genarate an empty model

lrtest(negbin52, negbin51a) #significant P-value means model fit
#choose model
lrtest(poisson52,negbin52)

#==========================================================
##Fivth Question RETURNQUANTITY
#==========================================================

poisson51r <- glm(returnquantity~time*group,family="poisson",data=mydata5)
stargazer(poisson51r,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson51ar <- glm(returnquantity~1, data=mydata5, family="poisson") # This is the command to run a logit on null model 
lrtest(poisson51r, poisson51ar)
#add control variables
poisson52r <- glm(returnquantity~time*group+avg_agemean+avg_femalemean+avg_incomemean+month_dummy+avg_childownermean+avg_homeownermean+salesquantity,family="poisson",data=mydata5)
stargazer(poisson51r,  poisson52r,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
lrtest(poisson52r, poisson51ar)

#negative binomial model
negbin51r <- glm.nb(returnquantity~time*group,data=mydata5)

stargazer(negbin51r,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#add controp varaiables
negbin52r <- glm.nb(returnquantity~time*group+avg_femalemean+avg_incomemean+avg_agemean+month_dummy+avg_childownermean+avg_homeownermean+salesquantity,data=mydata5)

stargazer(negbin51r,  negbin52r,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#check for heterokadasticity
gqtest(negbin52r) 
bptest(negbin52r)
HWrobstder <- sqrt(diag(vcovHC(negbin52r, type="HC1"))) # produces Huber-White robust standard errors 

stargazer(negbin52r,negbin52r,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))  # 
negbin51ar <- glm.nb(returnquantity ~ 1, data = mydata5) 

lrtest(negbin52r, negbin51ar) #significant P-value means model fit
#choose model
lrtest(poisson52r,negbin52r)


#==========================================================
##Sixth Question SALEVALUE
#==========================================================

m61 <- lm(log(salesvalue+1)~time*group*product_category,data=mydata5)
stargazer(m61, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
#add control variables
is.factor(mydata5$product_category)
mydata5$product_category <- as.factor(mydata5$product_category)
m62 <- lm(log(salesvalue+1)~time*group*product_category+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=mydata5)
stargazer(m61,m62, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 


# Check for heteroscedasticity
gqtest(m62)
bptest(m62)
HWrobstder <- sqrt(diag(vcovHC(m62, type="HC1")))
stargazer(m62,m62,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))


#==========================================================
##Sixth Question RETURNVALUE
#==========================================================

m61r <- lm(log(returnvalue+1)~time*group*product_category,data=mydata5)
stargazer(m61r, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
#add control variables
m62r <- lm(log(returnvalue+1)~time*group*product_category+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean+salesvalue,data=mydata5)
stargazer(m61r,m62r, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 


# Check for heteroscedasticity
gqtest(m62r)
bptest(m62r)
HWrobstder <- sqrt(diag(vcovHC(m62r, type="HC1")))
stargazer(m62r,m62r,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#==========================================================
##Sixth Question SALESQUANTITY
#==========================================================
poisson61 <- glm(salesquantity~time*group*product_category,family="poisson",data=mydata5)
stargazer(poisson61,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson61a <- glm(salesquantity~1, data=mydata5, family="poisson") # This is the command to run a logit on null model 
lrtest(poisson61, poisson61a)
#add control variables
poisson62 <- glm(salesquantity~time*group*product_category+avg_agemean+avg_femalemean+avg_incomemean+month_dummy+avg_childownermean+avg_homeownermean,family="poisson",data=mydata5)
stargazer(poisson61,  poisson62,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
lrtest(poisson62, poisson61a)

#negative binomial model
negbin61 <- glm.nb(salesquantity~time*group*product_category,data=mydata5)

stargazer(negbin61,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#add controp varaiables
negbin62 <- glm.nb(salesquantity~time*group*product_category+avg_femalemean+avg_incomemean+avg_agemean+month_dummy+avg_childownermean+avg_homeownermean,data=mydata5)

stargazer(negbin61,  negbin62,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#check for heterokadasticity
gqtest(negbin62) 
bptest(negbin62)
HWrobstder <- sqrt(diag(vcovHC(negbin62, type="HC1"))) # produces Huber-White robust standard errors 

stargazer(negbin62,negbin62,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))  # 
negbin61a <- glm.nb(salesquantity ~ 1, data = mydata5) #Genarate an empty model

lrtest(negbin62, negbin61a) #significant P-value means model fit
#choose model
lrtest(poisson62,negbin62)


#==========================================================
##Sixth Question RETURNQUANTITY
#==========================================================

poisson61r <- glm(returnquantity~time*group*product_category,family="poisson",data=mydata5)
stargazer(poisson61r,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson61ar <- glm(returnquantity~1, data=mydata5, family="poisson") # This is the command to run a logit on null model 
lrtest(poisson61r, poisson61ar)
#add control variables
poisson62r <- glm(returnquantity~time*group*product_category+avg_agemean+avg_femalemean+avg_incomemean+month_dummy+avg_childownermean+avg_homeownermean+salesquantity,family="poisson",data=mydata5)
stargazer(poisson61r,  poisson62r,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
lrtest(poisson62r, poisson61ar)

#negative binomial model
negbin61r <- glm.nb(returnquantity~time*group*product_category,data=mydata5)

stargazer(negbin61r,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#add controp varaiables
negbin62r <- glm.nb(returnquantity~time*group*product_category+avg_femalemean+avg_incomemean+avg_agemean+month_dummy+avg_childownermean+avg_homeownermean+salesquantity,data=mydata5)

stargazer(negbin61r,  negbin62r,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#check for heterokadasticity
gqtest(negbin62r) 
bptest(negbin62r)
HWrobstder <- sqrt(diag(vcovHC(negbin62r, type="HC1"))) # produces Huber-White robust standard errors 

stargazer(negbin62r,negbin62r,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))  # 
negbin61ar <- glm.nb(returnquantity ~ 1, data = mydata5) 

lrtest(negbin62r, negbin61ar) #significant P-value means model fit
#choose model
lrtest(poisson62r,negbin62r)

