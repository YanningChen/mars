source("mars.R")
source("plot.mars.R")
source("predict.mars.R")
source("print.mars.R")
source("summary.mars.R")
source("anova.mars.R")


library(ISLR)
data("Auto")
data("College")
data("Hitters")

# data1: Auto
mc=mars.control(Mmax=2)
Auto_mars<- mars(Auto$mpg~Auto$displacement+Auto$weight, data=Auto, control = mc) # test mars method
predict.mars(Auto_mars) # test predict method
print.mars(Auto_mars) # test print method
summary.mars(Auto_mars) # test summary method
anova.mars(Auto_mars) # test anova method
plot.mars(Auto_mars) # test plot method


# data2: College
College_mars=mars(College$Grad.Rate~College$F.Undergrad+College$P.Undergrad,data=College)


# data3: Hitters
Hitters_mars=mars(Hitters$Salary~Hitters$HmRun+Hitters$Hits
                  , data=Hitters)

