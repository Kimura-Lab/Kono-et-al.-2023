# Games-Howell test using Aoki_all.R
#Modify "Dataset", "group" and "response" as yours.

#####Games-Howell test#####
boxplot(response ~ factor(group), ylab="response", xlab="group", data=Dataset)
(res <- oneway.test(response ~ factor(group), data=Dataset, var.equal=FALSE))
(res <- tukey(Dataset$response, Dataset$group, method="G"))

#####Tukey-Kramer test#####
(res <- tukey(Dataset$response, Dataset$group)
