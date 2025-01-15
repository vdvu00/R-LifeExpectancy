## Vanessa Vu
## BS806 Final Project
## -------------------------------------------------------------------------

## load in data
setwd("/Users/vanessavu/Downloads/BS806 multivariable analysis/final project")
life <- read.csv("lifeexpectancy.csv")

## 2010 data only
life2010 <- subset(life, Year == 2010)
summary(life2010)
## remove year variable
life2010 <- subset(life2010, select = -Year)
summary(life2010)
## make status a factor
life2010$Status = as.factor(life2010$Status)
summary(life2010)
## make country the row names
rownames(life2010) = life2010[,1]
life2010 <- life2010[,-1]
summary(life2010)
## remove NAs from all variables 
life2010 <- na.omit(life2010)
summary(life2010)
## went from 183 to 128 obs

## SUMMARY STATS - life expectancy
summary(life2010$Life.expectancy)
hist(life2010$Life.expectancy)

## CORRELATION
library(ggplot2)
library(GGally)
ggpairs(life2010)
## too big of a plot

## CONDITION INDEX - COLLINEARITY
scaled2010 <- as.data.frame(scale(life2010[,2:20]))
scaled2010 <- cbind(life2010[,1], scaled2010)
colnames(scaled2010)[1] <- "Status"
model0 <- lm(Life.expectancy ~., data = scaled2010)
library(olsrr)
ols_eigen_cindex(model0)

## VARIABLE SELECTION WITH AIC (forward)
modelall <- ~Status + Adult.Mortality + infant.deaths + Alcohol + 
  percentage.expenditure + Hepatitis.B + Measles + BMI + under.five.deaths +
  Polio + Total.expenditure + Diphtheria + HIV.AIDS + GDP +
  Population + thinness..1.19.years + thinness.5.9.years +
  Income.composition.of.resources + Schooling
modelint <- lm(Life.expectancy ~ 1, data = life2010)

m.forward <- step(modelint, scope = modelall, direction = "forward", k=2)
summary(m.forward)$coefficients

## REGRESSION MODEL
model0 <- lm(Life.expectancy ~ Income.composition.of.resources +
               HIV.AIDS + Adult.Mortality + Schooling + Total.expenditure +
               Alcohol + Status, data = life2010)
summary(model0)
model1 <- lm(Life.expectancy ~ Income.composition.of.resources +
               HIV.AIDS + Adult.Mortality + Schooling + Total.expenditure +
               Alcohol, data = life2010)
summary(model1)
plot(model1$residuals)

## MODEL DIAGNOSTICS
## checking normality
plot(fitted(model1), residuals(model1), xlab="Fitted", ylab="Residuals", 
     main = "Residuals vs Fitted")
abline(h=0)
qqnorm(residuals(model1), ylab="Residuals")
qqline(residuals(model1))
hist(model1$residuals, xlab = "Residuals", main = "Initial Model Histogram")
shapiro.test(model1$residuals)
## null of shapiro-wilk test = normally distributed
## with alpha = 0.05, we fail to reject null so it is normally distributed

## equal variance & linearity
plot(fitted(model1), residuals(model1), xlab="Fitted", ylab="Residuals")
abline(h=0)
## curvature test (linearity)
library(car)
library(alr4)
residualPlots(model1)

## UNUSUAL POINTS
## outliers:
stresid <- rstudent(model1) # studentized residual
n <- length(stresid)
pprime <- 7 # num of parameters
## threshold for lower tail
thresh <- qt(0.05/(n*2), df = n-pprime-1, lower.tail=TRUE)
thresh
## checking which points are greater than threshold 
outlier <- abs(stresid) > abs(thresh)
which((outlier)==1)
## no outliers according to this threshold

## influential points:
cooks <- cooks.distance(model1)
check <- cooks[cooks > 4/n] # rule of thumb
sort(check, decreasing = TRUE) 
## SHOULD USE THIS THRESHOLD BELOW
cooks[cooks>0.5] # check Di>0.5
influenceIndexPlot(model1)

## leverage points:
x <- as.matrix(data.frame(x0=rep(1,nrow(life2010)),
                          income=life2010$Income.composition.of.resources, 
                          hiv=life2010$HIV.AIDS, 
                          adultmort=life2010$Adult.Mortality, 
                          schooling=life2010$Schooling, 
                          totalexp=life2010$Total.expenditure,
                          alcohol=life2010$Alcohol))
hat <- x%*%solve(t(x)%*%x)%*%t(x)
lev <- diag(hat)
names(lev) = rownames(life2010)
## rule of thumb - leverage > 2p'/n
levrule <- lev[lev > 2*(pprime/n)]
levrule
sort(levrule, decreasing = TRUE) 

## REMOVE BHUTAN AND TURKMENISTAN
life2010new <- life2010[!rownames(life2010) %in% c("Bhutan", "Turkmenistan"), ]

## REMODEL WITH NEW DATASET FOR 2010
finallm <- lm(Life.expectancy ~ Income.composition.of.resources +
                HIV.AIDS + Adult.Mortality + Schooling + Total.expenditure +
                Alcohol, data = life2010new)
summary(finallm)

## MODEL DIAGNOSTICS ON UPDATED MODEL
## checking normality
plot(fitted(finallm), residuals(finallm), xlab="Fitted", ylab="Residuals", 
     main = "Residuals vs Fitted")
abline(h=0)
qqnorm(residuals(finallm), ylab="Residuals")
qqline(residuals(finallm))
hist(finallm$residuals, xlab = "Residual", main = "Final Model Histogram")
shapiro.test(finallm$residuals)
## null of shapiro-wilk test = normally distributed
## with alpha = 0.05, we fail to reject null so it is normally distributed

## equal variance & linearity
plot(fitted(finallm), residuals(finallm), xlab="Fitted", ylab="Residuals",
     main = "Residuals vs Fitted")
abline(h=0)
## curvature test (linearity)
residualPlots(finallm)

## CONFOUNDING
## income comp of resources is most significant
## test confounding by other predictors on income comp of resources
crudemod <- lm(Life.expectancy ~ Income.composition.of.resources, data = life2010new)
summary(crudemod)
## beta estimate of adj model = 53.466092
## beta estimate of crude model = 50.679 
## 10% rule:
rule <- abs((50.679-53.466092)/50.679)*100
rule
## no confounding by group of other covariates

## REGRESSION TREE
tree1 <- tree(Life.expectancy ~ Status + Adult.Mortality + infant.deaths + Alcohol + 
              percentage.expenditure + Hepatitis.B + Measles + BMI + under.five.deaths +
              Polio + Total.expenditure + Diphtheria + HIV.AIDS + GDP +
              Population + thinness..1.19.years + thinness.5.9.years +
              Income.composition.of.resources + Schooling, data = life2010new, 
              control=tree.control(nobs=nrow(life2010new)), mindev = 0.001)
summary(tree1)
set.seed(1)
tree.list <- cv.tree(tree1, FUN=prune.tree, K=10)
tree.list
## lowest deviance is K=9
## best tree
finaltree <- prune.tree(tree1, best=9)
summary(finaltree)
plot(finaltree)
text(finaltree, pretty=0, cex=0.8)
finaltree

## TREE VS LM 
## dummy variables
life2010new$branch1 <- 0;
life2010new$branch1[life2010new$HIV.AIDS < 0.7 & life2010new$Adult.Mortality < 363.5 &
                   life2010new$Income.composition.of.resources < 0.534] <- 1;
life2010new$branch2 <- 0;
life2010new$branch2[life2010new$HIV.AIDS >= 0.7 & life2010new$Adult.Mortality < 363.5 &
                   life2010new$Income.composition.of.resources < 0.534] <- 1;
life2010new$branch3 <- 0;
life2010new$branch3[life2010new$Adult.Mortality >= 363.5 &
                   life2010new$Income.composition.of.resources < 0.534] <- 1;
life2010new$branch4 <- 0;
life2010new$branch4[life2010new$Income.composition.of.resources < 0.6005 & 
                   life2010new$HIV.AIDS < 1.4 & 
                   life2010new$Income.composition.of.resources < 0.8025 &
                   life2010new$Income.composition.of.resources >= 0.534] <- 1;
life2010new$branch5 <- 0;
life2010new$branch5[life2010new$Adult.Mortality < 168.5 &
                  life2010new$Income.composition.of.resources >= 0.6005 & 
                   life2010new$HIV.AIDS < 1.4 & 
                   life2010new$Income.composition.of.resources < 0.8025 &
                   life2010new$Income.composition.of.resources >= 0.534] <- 1;
life2010new$branch6 <- 0;
life2010new$branch6[life2010new$Adult.Mortality >= 168.5 &
                   life2010new$Income.composition.of.resources >= 0.6005 & 
                   life2010new$HIV.AIDS < 1.4 & 
                   life2010new$Income.composition.of.resources < 0.8025 &
                   life2010new$Income.composition.of.resources >= 0.534] <- 1;
life2010new$branch7 <- 0;
life2010new$branch7[life2010new$HIV.AIDS >= 1.4 &
                   life2010new$Income.composition.of.resources < 0.8025 &
                   life2010new$Income.composition.of.resources >= 0.534] <- 1;
life2010new$branch8 <- 0;
life2010new$branch8[life2010new$Income.composition.of.resources < 0.856 &
                   life2010new$Income.composition.of.resources >= 0.8025 &
                   life2010new$Income.composition.of.resources >= 0.534] <- 1;
life2010new$branch9 <- 0;
life2010new$branch9[life2010new$Income.composition.of.resources >= 0.856 &
                   life2010new$Income.composition.of.resources >= 0.8025 &
                   life2010new$Income.composition.of.resources >= 0.534] <- 1;
## linear model of tree
treelm <- lm(Life.expectancy ~ branch2 + branch3 + branch4 + branch5 + branch6 + branch7 + 
               branch8 + branch9, data = life2010new)
summary(treelm)
## compare to lm 
summary(finallm)

## ANOVA & AIC
anova(finallm)
anova(treelm)
anova(finallm, treelm)

AIC(finallm)
AIC(treelm)

## INCORPORATE SHIV'S CODE
anova(treelm, finallm)

summary(treelm) #R2 = .924
summary(finallm) #R2 = .906
AIC(treelm) #AIC = 602.34
AIC(finallm) #AIC = 625.08

### to find a "best fit" - we will use AIC to compare and ANOVA to compare
### but we need to compare AIC with the same number of observations, so we will use
### the dataset with the outliers removed - wanna compare it to the model with the 
### outliers removed

#Hand calculate ANOVA from linear regression model - treelm
#residual standard error = 2.532
#s2 = (RSE^2) = (2.532^2)
#n = number of observations = 126
nrow(life2010new)

#p = number of predictors = 8
#RSS = s2*(n-p-1)
(2.532^2)*(126-8-1)
deviance(treelm)
df.residual(treelm)

summary(treelm)
#R2 = 0.924

RSS.reg.tree <- (2.532^2)*(126-8-1)
df.reg.tree <- (126-8-1)
MSE.reg.tree <- RSS.reg.tree/(126-8-1)
MSE.reg.tree
#SYY = RSS/(1-R2)
SYY.reg.tree = (RSS.reg.tree)/(1-0.924)
SYY.reg.tree
#SSreg - R2*SYY
SSreg.reg.tree <- (0.924)*(SYY.reg.tree)
SSreg.reg.tree
#MSReg = SSReg/p
MSreg.reg.tree <- SSreg.reg.tree/(8);
MSreg.reg.tree


#Hand calculate ANOVA from linear regression model - m.new
summary(finallm)
#R2 = .906
#RSE = 2.792
#(RSE)^2 = s2 = (2.792)^2
#n = 126
#p = 6
RSS.lin.reg <- ((2.792)^2)*(126-6-1)
deviance(finallm)
MSE.lin.reg <- (RSS.lin.reg)/(126-6-1)
MSE.lin.reg

SYY.lin.reg <- (RSS.lin.reg)/(1-.906)
SYY.lin.reg
SSreg.lin.reg <- SYY.lin.reg*.906
SSreg.lin.reg
MSReg.lin.reg <- (SSreg.lin.reg)/6
MSReg.lin.reg

#check with ANOVA
anova(treelm, finallm)



## -------------------------------------------------------------------------

## EXTRA CODE ON OTHER YEARS 
## 2000 data only
life2000 <- subset(life, Year == 2000)
summary(life2000)
## remove year variable
life2000 <- subset(life2000, select = -Year)
summary(life2000)
## make status a factor
life2000$Status = as.factor(life2000$Status)
summary(life2000)
## make country the row names
rownames(life2000) = life2000[,1]
life2000 <- life2000[,-1]
summary(life2000)
## make sure same countries are in 2000 data as in 2010
life2000 <- life2000[rownames(life2010),]
summary(life2000)
## TOTAL EXPENDITURE MISSING FOR IRAQ IN 2000 DATA
## SHOULD WE REMOVE FROM 2010 DATA TO COMPARE FAIRLY?

## REGRESSION FOR 2000
## using same predictors as 2010
model2 <- lm(Life.expectancy ~ Income.composition.of.resources +
               HIV.AIDS + Adult.Mortality + Schooling + Total.expenditure +
               Alcohol + Status, data = life2000)
summary(model2)
summary(model1)
## different significance

## remove NAs from 2000 data - only 60 countries left
life2000_nona <- na.omit(life2000)
summary(life2000_nona)


## 2012 data only
life2012 <- subset(life, Year == 2012)
summary(life2012)
## remove year variable
life2012 <- subset(life2012, select = -Year)
summary(life2012)
## make status a factor
life2012$Status = as.factor(life2012$Status)
summary(life2012)
## make country the row names
rownames(life2012) = life2012[,1]
life2012 <- life2012[,-1]
summary(life2012)
## make sure same countries are in 2012 data as in 2010
life2012 <- life2012[rownames(life2010),]
summary(life2012)

## REGRESSION FOR 2012
## using same predictors as 2010
model3 <- lm(Life.expectancy ~ Income.composition.of.resources +
               HIV.AIDS + Adult.Mortality + Schooling + Total.expenditure +
               Alcohol, data = life2012)
summary(model3)
summary(model1)

