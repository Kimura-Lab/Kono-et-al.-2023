# Games-Howell and Tukey-Kramer test using Aoki_all.R
#Change "Dataset", "group" and "response" as yours.  
Modified EZR1.61 for Windows is also available.

#####Games-Howell test#####  
`boxplot(response ~ factor(group), ylab="response", xlab="group", data=Dataset)`  
`(res <- oneway.test(response ~ factor(group), data=Dataset, var.equal=FALSE))`  
`(res <- tukey(Dataset$response, Dataset$group, method='G'))` 

#####Tukey-Kramer test#####  
`(res <- oneway.test(response ~ factor(group), data=Dataset, var.equal=TRUE))`  
`(res <- tukey(Dataset$response, Dataset$group))`  
