# Games-Howell test using Hatano_socialstatisticsbasic.R
#Modify "Dataset", "group" and "response" as yours.

#load dataset
summary(Dataset)
value <- Dataset$response
group <- factor(Dataset$group)

#one-way ANOVA and post-hoc multiple comparison
output.anova <- oneway.factorial.anova(value~group,boxplot=T)
output.compare <- multiple.comparison.test(value~group, method="G")

#Output results
output.anova
output.compare
