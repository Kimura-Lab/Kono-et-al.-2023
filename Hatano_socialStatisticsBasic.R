# Functions for Social Statistics(Basic)
# 
# Functions List
# 
# descriptive(value,start,width,col)
# descriptive.grouping(formula,data,col)
# z.test(x,mu,sigma2,a)
# t.test.paired(x1,x2,a)
# t.test.independent(formula,group,no1,no2,a)
# oneway.factorial.anova(formula,data,boxplot,col)
# oneway.repeated.measures.anova(data)
# twoway.factorial.anova(formula,data,type,interaction)
# simple.main.effects(formula,data,anova)
# multiple.comparison.test(formula,data,repeated,anova,method,fact,summary)
# Levene.test(formula,data,center,boxplot,col)
# crosstab(row,column)
# correl(data,control,use,method)
# linest(formula,data,direction,step)
# samplingTest(x,N,num,col)
# convert.crosstable.to.weighted.data(matrix)
# convert.data.frame.crosstable.to.weighted.data(data.frame)
# convert.weighted.data.to.original(data.frame,freq.name)
# recode.dummy.variables(data, ref=NULL)
# write.output(result,filename)
# 
# 2016/10/18 Ver.1.0
# 2016/11/12 Ver.1.1 refine function descriptive.grouping(), linest(), oneway.factorial.anova(), add function Levene.test()
# 2016/11/21 Ver.1.2 add function twoway.factorial.anova()
# 2018/01/08 Ver.1.2.1 bug fixed in function linest()
# 2021/11/24 Ver.1.3 refine function t.test.independent()
# 

# Calculate descriptive statistics and output histogram.
# 
# @param value numeric vector.
# @param start numeric scalar.
# @param width numeric scalar.
# @param right boolean scalar.
# @param col char scalar.
# @return list matrix
descriptive <- function(value,start=min(value),width=(ceiling(max(value) - min(value))/ceiling(log(length(value)+1))),right=FALSE,col="royalblue"){
    origin.value <- value
    value <- na.omit(origin.value)
    end <- max(value) + width
        right <- right
    
    #histogram & mode
    cut <- seq(start,end,width)
    hist <- hist(value,cut,right=right,col=col)
    freq <- matrix(hist$counts)
    rownames(freq) <- format(cut[-length(cut)],nsmall=2)
    colnames(freq) <- c("freq")
                         
    mode <- mean(hist$mids[which(hist$counts == max(hist$counts))] )

    #fundamental statistics
    items <- c("N","mean","median","mode","max","min","range","V","sd","u2","u","Missing value")
    descriptiv <- matrix(0,12,1)
    descriptiv[,1] <- c(
        N <- length(value),
        mean.value <- mean(value),
        med <- median(value),
        mode,
        max.value <- max(value),
        min.value <- min(value),
        range <- max.value - min.value,
        V <- sum((value - mean.value)^2)/N,
        sd <- sqrt(V),
        u2 <- var(value),
        u <- sd(value),
        length(origin.value) - length(value)
    )
    colnames(descriptiv) <- c("statistics")
    rownames(descriptiv) <- items
    
    result <- list(statistics=data.matrix(descriptiv),freq=freq)
    class(result) <- c("listMatrix")

    return(result)
}

# Calculate descriptive statistics for grouped data and output boxplot.
# 
# @param formula formula object.
# @param data data.frame.
# @param col char scalar
# @return matrix
descriptive.grouping <- function(formula,data=NULL,col="darkorange"){
    if(is.vector(formula) && is.factor(data)){
        origin.dat <- data.frame(formula,data)
        dat <-  na.omit(data.frame(formula,data))
        formula <- formula~data
    }

    else if(is.null(data)){
        origin.dat <- model.frame(formula,na.action=na.pass)
        dat <- model.frame(formula,na.action=na.omit)
    }else{
        origin.dat <- model.frame(formula,data,na.action=na.pass)
        dat <- model.frame(formula,data,na.action=na.omit)
    }

    origin.value <- origin.dat[[1]]
    origin.group <- origin.dat[[2]]
    value <- dat[[1]]
    group <- factor(dat[[2]])

    #summary
    group.N <- tapply(value,group,length)
    origin.group.N <- tapply(origin.value,origin.group,length)
    group.mean <- tapply(value,group,mean)
    group.var <- tapply(value,group,var)
    group.min <- tapply(value,group,min)
    group.max <- tapply(value,group,max)
    group.median <- tapply(value,group,median)

    summary <- matrix(0,length(summary(group)) + 1,7)
    summary[,1] <- c(matrix(group.N),length(value))
    summary[,2] <- c(matrix(group.mean),mean(value))
    summary[,3] <- c(matrix(group.var),var(value))
    summary[,4] <- c(matrix(group.min),min(value))
    summary[,5] <- c(matrix(group.max),max(value))
    summary[,6] <- c(matrix(group.median),median(value))
    summary[,7] <- c(matrix(origin.group.N)-matrix(group.N),length(origin.value)-length(value))
    rownames(summary) <- c(names(summary(group)),"Sum")
    colnames(summary) <- c("N","mean","u2","min","max","median","Missing value")

    #boxplot
    boxplot(split(value,group),col=col)
    
    result <- list(summary=summary)
    class(result) <- c("listMatrix")

    return(result)
}

# Execute z-test or t-test(One sample) and output result of matrix type.
# 
# @param x numeric vector. Sampling Data.
# @param mu numeric scalar. Population Mean.
# @param sigma2 numeric scalar. Population Variance.
# @param a numeric scalar. Significance level.
# @return list matrix
z.test <- function(x,mu=0,sigma2=0,a=0.05){
    origin.x <- x
    x <- na.omit(origin.x)
    #statistics
    statistics <- matrix(0,6,1)
    statistics[,1] <- c(
        N <- length(x),
        df <- N - 1,
        mean.value <- mean(x),
        u2 <- var(x),
        u <- sd(x),
        length(origin.x) - length(x)
    )
    colnames(statistics) <- c("summary")
    rownames(statistics) <- c("N","df","mean","u2","u","Missing Value")
    
    if(sigma2 != 0){
        #Z-test
        z.result <- matrix(0,8,1)
        z.result[,1] <- c(
            SE <- sqrt(sigma2/N),
            Z <- (mean.value - mu)/SE,
            rejection <- qnorm(1 - a/2),
            p.2 <- (1 - pnorm(abs(Z)))*2,
            p.greater <- 1 - pnorm(Z),
            p.less <- pnorm(Z),
            upper <- mean.value + rejection*SE,
            lower <- mean.value - rejection*SE
        )
        colnames(z.result) <- c("z-test")
        rownames(z.result) <- c("SE","Z","rejection","p(!=)","p(>)","p(<)","Upper","Lower")

        #VAR-test
        var.result <- matrix(0,4,1)
        var.result[,1] <- c(
            chisq <- sum((x - mean.value)^2)/sigma2,
            p.2 <- if(u2 > sigma2){(1 - pchisq(chisq,N - 1))*2}else{pchisq(chisq,N - 1)*2},
            p.greater <- 1 - pchisq(chisq,N - 1),
            p.less <- pchisq(chisq,N - 1)
        )
        colnames(var.result) <- c("var-test")
        rownames(var.result) <- c("chisq","p(!=)","p(>)","p(<)")
        
        result <- list(statistics=statistics,z_test=z.result,var_test=var.result)
    }else{
        #t-test
        t.result <- matrix(0,10.1)
        t.result[,1] <- c(
            SE <- sqrt(var(x)/N),
            t <- (mean.value - mu)/SE,
            rejection <- qt(1 - a/2,df),
            p.2 <- (1 - pt(abs(t),df))*2,
            p.greater <- 1 - pt(t,df),
            p.less <- pt(t,df),
            upper <- mean.value + rejection*SE,
            lower <- mean.value - rejection*SE,
            d <- (mean.value - mu)/u,
            r <- sqrt(t^2/(t^2+N-1))
        )
        colnames(t.result) <- c("t-test")
        rownames(t.result) <- c("SE","t","rejection","p(!=)","p(>)","p(<)","Upper","Lower","d","r")
        
        result <- list(statistics=statistics,t_test=t.result)
    }
    class(result) <- c("listMatrix")

    return(result)
}

# Execute t-test(paired data) and output result of matrix type.
# 
# @param x1 numeric vector.
# @param x2 numeric vector.
# @param a numeric scalar.
# @return list matrix
t.test.paired <- function(x1,x2,a=0.05){
    omit.x <- na.omit(data.frame(x1,x2))
    omit.x1 <- omit.x$x1
    omit.x2 <- omit.x$x2
    diff <- omit.x1 - omit.x2
    N <- length(diff)

    #statistics
    statistics <- matrix(0,5,2)
    statistics[1,] <- c(length(omit.x1),length(omit.x2))
    statistics[2,] <- c(mean(omit.x1,na.rm=TRUE),mean(omit.x2,na.rm=TRUE))
    statistics[3,] <- c(sd(omit.x1,na.rm=TRUE),sd(omit.x2,na.rm=TRUE))
    statistics[4,] <- c(sqrt(var(omit.x1,na.rm=TRUE)/length(omit.x1)),sqrt(var(omit.x2,na.rm=TRUE)/length(omit.x2)))
    statistics[5,] <- c(length(x1)-length(omit.x1),length(x2)-length(omit.x2))
    colnames(statistics) <- c("x1","x2")
    rownames(statistics) <- c("N","mean","u","SE","Missing value")

    #t-test
    result <- matrix(0,12,1)
    result[,1] <- c(
        mean(diff),
        sd <- sd(diff),
        SE <- sqrt(var(diff)/length(diff)),
        t.value <- (mean(diff)) / SE,
        df <-  N - 1,
        p.2 <- (1 - pt(abs(t.value),df))*2,
        p.greater <- 1 - pt(t.value,df),
        p.less <- pt(t.value,df),
        upper <- mean(diff) + qt(1 - a/2,df)*SE,
        lower <- mean(diff) - qt(1 - a/2,df)*SE,
        d <- mean(diff)/sqrt(sum((diff - mean(diff))^2)/df),
        r <- sqrt(t.value^2/(t.value^2+df))
    )
    colnames(result) <- c("t-test")
    rownames(result) <- c("mean diff","sd","SE","t","df","p(!=)","p(>)","p(<)","Upper","Lower","d","r")
    
    result <- list(statistics=statistics,test.result=result)
    class(result) <- c("listMatrix")

    return(result)
}

# Execute t-test(independent data) and output result of matrix type.
# 
# @param formula formula object OR numeric vector.
# @param data data.frame OR numeric vector OR factor.
# @param no1 char scalar.
# @param no2 char scalar.
# @param a numeric scalar.
# @return list matrix.
t.test.independent <- function(formula,data=NULL,no1=NULL,no2=NULL,a=0.05,var.equal=FALSE){
    if(is.vector(formula) && is.vector(data) && is.numeric(data)){
        x.grouped <- list(formula,data)
        if(is.null(no1) || is.null(no2)){
            no1 <- "x"
            no2 <- "y"
        }
        is.both.vec <- TRUE
    }
    else{
        is.both.vec <- FALSE
        if(is.vector(formula) && (is.factor(data) || is.vector(data) && is.character(data))){
            origin.dat <- data.frame(formula,data)
            dat <-  na.omit(data.frame(formula,data))
            formula <- formula~data
        }
        else if(is.null(data)){
            origin.dat <- model.frame(formula,na.action=na.pass)
            dat <- model.frame(formula,na.action=na.omit)
        }
        else{
            origin.dat <- model.frame(formula,data,na.action=na.pass)
            dat <- model.frame(formula,data,na.action=na.omit)
        }
        origin.x <- origin.dat[[1]]
        origin.group <- factor(origin.dat[[2]])
        x <- dat[[1]]
        group <- factor(dat[[2]])

        if(is.null(no1) || is.null(no2)){
            origin.group.grouped <- split(origin.group,origin.group)
            group.grouped <- split(group,group)
            origin.x.grouped <- split(origin.x,origin.group)
            x.grouped <- split(x,group)
            no1 <- levels(group)[1]
            no2 <- levels(group)[2]
        }
        else{
            origin.group.grouped <- list(origin.group[origin.group==no1],origin.group[origin.group==no2])
            group.grouped <- list(group[group==no1],group[group==no2])
            origin.x.grouped <- list(origin.x[origin.group==no1],origin.x[origin.group==no2])
            x.grouped <- list(x[group==no1],x[group==no2])
        }
    }

    SE1 <- sqrt(var(x.grouped[[1]])/length(x.grouped[[1]]))
    SE2 <- sqrt(var(x.grouped[[2]])/length(x.grouped[[2]]))

    df1 <- length(x.grouped[[1]])-1
    df2 <- length(x.grouped[[2]])-1
    if(df1 < 2 || df2 < 2){
        stop("Probably, definition of grouping variable is incorrect.")
    }

    #statistics
    statistics <- matrix(0,5,2)
    statistics[1,] <- c(length(x.grouped[[1]]),length(x.grouped[[2]]))
    statistics[2,] <- c(mean(x.grouped[[1]]),mean(x.grouped[[2]]))
    statistics[3,] <- c(sd(x.grouped[[1]]),sd(x.grouped[[2]]))
    statistics[4,] <- c(sqrt(var(x.grouped[[1]])/length(x.grouped[[1]])),sqrt(var(x.grouped[[2]])/length(x.grouped[[2]])))
    if(is.both.vec){
        statistics[5,] <- c(0,0)
    }
    else{
        statistics[5,] <- c(length(na.omit(origin.group.grouped[[1]]))-length(group.grouped[[1]]),length(na.omit(origin.group.grouped[[2]]))-length(group.grouped[[2]]))
    }
    colnames(statistics) <- c(no1,no2)
    rownames(statistics) <- c("N","mean","u","SE","Missing value")

    #t-test
    mean.diff <- mean(x.grouped[[1]])-mean(x.grouped[[2]])
    if(var.equal){
        df <- length(x.grouped[[1]]) - 1 + length(x.grouped[[2]]) - 1
        v <- (sum((x.grouped[[1]] - mean(x.grouped[[1]]))^2) + sum((x.grouped[[2]] - mean(x.grouped[[2]]))^2))/df
        SE <- sqrt(v/length(x.grouped[[1]]) + v/length(x.grouped[[2]]))
    }
    else{
        SE <- sqrt(var(x.grouped[[1]])/length(x.grouped[[1]]) + var(x.grouped[[2]])/length(x.grouped[[2]]))
        df <- SE^4/(SE1^4/df1+SE2^4/df2)
    }
    
    result <- matrix(0,11,1)
    result[,1] <- c(
        t.value <- mean.diff/SE,
        df,
        p.2 <- (1 - pt(abs(t.value),df))*2,
        p.greater <- 1 - pt(t.value,df),
        p.less <- pt(t.value,df),
        mean.diff,
        SE,
        upper <- mean.diff + qt(1 - a/2,df)*SE,
        lower <- mean.diff - qt(1 - a/2,df)*SE,
        d <- (mean(x.grouped[[1]],na.rm=TRUE)-mean(x.grouped[[2]],na.rm=TRUE))/sqrt((df1*var(x.grouped[[1]],na.rm=TRUE)+df2*var(x.grouped[[2]],na.rm=TRUE))/(df1 + df2)),
        r <- sqrt(t.value^2/(t.value^2+df))
    )
    colnames(result) <- c("t-test")
    rownames(result) <- c("t","df","p(!=)","p(>)","p(<)","mean diff","SE","Upper","Lower","d","r")
    
    result <- list(statistics=statistics,test.result=result)
    class(result) <- c("listMatrix")

    return(result)
}

# Execute oneway-factorial-anova and output result of matrix type.
# 
# @param formula formula object.
# @param data data.frame.
# @param boxplot logical scalar.
# @param col char scalar.
# @return list matrix.
oneway.factorial.anova <- function(formula, data=NULL, boxplot=FALSE, col="darkorange"){
    if(is.vector(formula) && is.factor(data) || is.vector(data)){
        origin.dat <- data.frame(formula, data)
        dat <-  na.omit(data.frame(formula, data))
        formula <- formula~data
    }else if(is.null(data)){
        origin.dat <- model.frame(formula,na.action=na.pass)
        dat <- model.frame(formula,na.action=na.omit)
    }else{
        origin.dat <- model.frame(formula,data,na.action=na.pass)
        dat <- model.frame(formula,data,na.action=na.omit)
    }
    origin.value <- origin.dat[[1]]
    origin.group <- factor(origin.dat[[2]])
    value <- dat[[1]]
    group <- factor(dat[[2]])

    #summary
    group.N <- tapply(value,group,length)
    origin.group.N <- tapply(origin.value,group,length)
    group.mean <- tapply(value,group,mean)
    group.var <- tapply(value,group,var)

    summary <- matrix(0,length(summary(group)) + 1,5)
    summary[,1] <- c(matrix(group.N),length(value))
    summary[,2] <- c(matrix(group.N)-1,length(value)-length(summary(group)))
    summary[,3] <- c(matrix(group.mean),mean(value))
    summary[,4] <- c(matrix(group.var),var(value))
    summary[,5] <- c(matrix(origin.group.N - group.N),length(origin.value)-length(value))
    rownames(summary) <- c(names(summary(group)),"Sum")
    colnames(summary) <- c("N","df","mean","u2","Missing Value")
    if(boxplot){
        boxplot(split(value,group),col=col)
    }

    #ANOVA
    #degree of freedom
    factor.df <- length(group.mean) - 1
    residual.df <- sum(group.N -1)
    total.df <- length(value) - 1

    #Sum of Squares
    factor.ss <- sum((group.mean - mean(value))^2*group.N)
    residual.ss <- sum(tapply(value,group,function(x){sum((x -mean(x))^2)}))
    total.ss <- sum((value - mean(value))^2)

    #Mean of Squares
    factor.ms <- factor.ss / factor.df
    residual.ms <-  residual.ss / residual.df
    total.ms <- total.ss / total.df

    eta <- factor.ss/total.ss

    f.value <- factor.ms / residual.ms
    p.value <- 1 - pf(f.value,factor.df,residual.df)

    ANOVA <- matrix(0,2,6)
    ANOVA[,1] <- c(factor.ss,residual.ss)
    ANOVA[,2] <- c(factor.df,residual.df)
    ANOVA[,3] <- c(factor.ms,residual.ms)
    ANOVA[,4] <- c(f.value,NA)
    ANOVA[,5] <- c(p.value,NA)
    ANOVA[,6] <- c(eta,NA)

    colnames(ANOVA) <- c("Sum Sq","df","Mean Sq","F","P","eta2")
    rownames(ANOVA) <- c("Factor","Residual")

    #Welch
    group.rec.se2 <- group.N / group.var

    welch.f.value <- sum(group.rec.se2*(group.mean - sum(group.rec.se2 * group.mean)/sum(group.rec.se2))^2)/(factor.df*(1 + 2*(factor.df -1)/(length(group.mean)^2 - 1)*sum((1 - group.N/group.var/sum(group.rec.se2))^2/(group.N - 1))))
    welch.df <- (length(group.mean)^2 - 1)/(3*sum((1 - group.N/group.var/sum(group.rec.se2))^2/(group.N - 1)))
    welch.p.value <- 1 - pf(welch.f.value,factor.df,welch.df)

    test.result <- matrix(0,2,4)
    test.result[1,] <- c(f.value,factor.df,residual.df,p.value)
    test.result[2,] <- c(welch.f.value,factor.df,welch.df,welch.p.value)
    colnames(test.result) <- c("F","df1","df2","p")
    rownames(test.result) <- c("Fisher","Welch")

    result <- list(summary=summary,anova=ANOVA,test.result=test.result)

    class(result) <- c("listMatrix")
    
    return(result)
}

# Execute oneway-repeated-measures-anova and output result of matrix type.
# 
# @param data data.frame.
# @return list matrix.
oneway.repeated.measures.anova <- function(data){
    origin.data <- data
    data <- na.omit(data)

    samples.mean <- rowMeans(data)
    factor.mean <- colMeans(data)
    all.mean <- mean(samples.mean)

    samples.N <- length(samples.mean)
    factor.N <- length(factor.mean)
    N <- samples.N*factor.N
    missing.Value.N <- length(rowMeans(origin.data)) - samples.N
    samples.df <- samples.N - 1
    factor.df <- factor.N - 1
    res.df <- N - samples.df - factor.df -1

    summary <- matrix(0, factor.N + 1, 3)
    summary[,1] <- c(rep(samples.N,factor.N),N)
    summary[,2] <- c(factor.mean, all.mean)
    summary[,3] <- c(rep(missing.Value.N, factor.N), missing.Value.N*factor.N)
    rownames(summary) <- c(names(data),"Sum")
    colnames(summary) <- c("N","mean","Missing Value")

    samples.effects <- samples.mean - all.mean
    factor.effects <- factor.mean - all.mean

    factor.effects.matrix <- matrix(rep(factor.effects,samples.N),ncol=factor.N,byrow=T)

    samples.ss <- sum(samples.effects^2)*factor.N
    factor.ss <- sum(factor.effects^2)*samples.N
    res.ss <-sum((data - samples.effects - all.mean - factor.effects.matrix)^2)

    samples.ms <- samples.ss / samples.df
    factor.ms <- factor.ss / factor.df
    res.ms <- res.ss / res.df

    factor.F <- factor.ms / res.ms

    factor.p <- 1 - pf(factor.F,factor.df,res.df)
    
    ANOVA <- matrix(0,3,5)
    ANOVA[,1] <- c(factor.ss,samples.ss,res.ss)
    ANOVA[,2] <- c(factor.df,samples.df,res.df)
    ANOVA[,3] <- c(factor.ms,samples.ms,res.ms)
    ANOVA[,4] <- c(factor.F,NA,NA)
    ANOVA[,5] <- c(factor.p,NA,NA)

    colnames(ANOVA) <- c("Sum Sq","df","Mean Sq","F","P")
    rownames(ANOVA) <- c("factor","samples","Residual")
    
    result <- list(summary=summary,anova=ANOVA)
    class(result) <- c("listMatrix")
    return(result)
}

# Execute twoway-factorial-anova and output result of matrix type.
# 
# @param formula formula object.
# @param data data.frame.
# @param type 1 or 2.
# @param interaction  logical scalar.
# @return list matrix.
twoway.factorial.anova <- function(formula,data=NULL,type=2,interaction=FALSE){
    if(is.null(data)){
        origin.dat <- model.frame(formula,na.action=na.pass)
        dat <- model.frame(formula,na.action=na.omit)
    }else{
        origin.dat <- model.frame(formula,data,na.action=na.pass)
        dat <- model.frame(formula,data,na.action=na.omit)
    }
    origin.value <- origin.dat[[1]]
    origin.factor1 <- origin.dat[[2]]
    origin.factor2 <- origin.dat[[3]]
    value <- dat[[1]]
    factor1 <- factor(dat[[2]])
    factor2 <- factor(dat[[3]])

    factor1.mean <- tapply(value,factor1,mean)
    factor2.mean <- tapply(value,factor2,mean)

    factor1.N <- tapply(value,factor1,length)
    factor2.N <- tapply(value,factor2,length)

    factor1.num <- length(factor1.N)
    factor2.num <- length(factor2.N)
    factor1.df <- factor1.num - 1
    factor2.df <- factor2.num -1
    mix.df <- factor2.df*factor1.df
    res.df <- length(value) -1 - (factor2.df + factor1.df + mix.df)

    #summary
    origin.factor1.N <- tapply(origin.value,origin.factor1,length)
    origin.factor2.N <- tapply(origin.value,origin.factor2,length)
    factor1.var <- tapply(value,factor1,var)
    factor2.var <- tapply(value,factor2,var)

    summary.factor1 <- matrix(0,length(summary(factor1)) + 1,5)
    summary.factor1[,1] <- c(matrix(factor1.N),length(value))
    summary.factor1[,2] <- c(matrix(factor1.N)-1,length(value)-length(summary(factor1)))
    summary.factor1[,3] <- c(matrix(factor1.mean),mean(value))
    summary.factor1[,4] <- c(matrix(factor1.var),var(value))
    summary.factor1[,5] <- c(matrix(origin.factor1.N - factor1.N),length(origin.value)-length(value))
    rownames(summary.factor1) <- c(names(summary(factor1)),"Sum")
    colnames(summary.factor1) <- c("N","df","mean","u2","Missing Value")

    summary.factor2 <- matrix(0,length(summary(factor2)) + 1,5)
    summary.factor2[,1] <- c(matrix(factor2.N),length(value))
    summary.factor2[,2] <- c(matrix(factor2.N)-1,length(value)-length(summary(factor1)))
    summary.factor2[,3] <- c(matrix(factor2.mean),mean(value))
    summary.factor2[,4] <- c(matrix(factor2.var),var(value))
    summary.factor2[,5] <- c(matrix(origin.factor2.N - factor2.N),length(origin.value)-length(value))
    rownames(summary.factor2) <- c(names(summary(factor2)),"Sum")
    colnames(summary.factor2) <- c("N","df","mean","u2","Missing Value")

    mix.num <- table(factor1,factor2)
    mix.mean <- tapply(value,list(factor1,factor2),function(x)mean(x))

    #ANOVA
    if(type == 1){
        factor1.tmp.ss <- factor1.ss <- sum((sapply(factor1,function(x)mean(value[which(factor1==x)])) - mean(value))^2)
        factor2.tmp.ss <- factor2.ss <- drop1(lm(value~factor1+factor2),"factor2")$"Sum of Sq"[[2]]
    }else{
        factor1.tmp.ss <- drop1(lm(value~factor1),"factor1")$"Sum of Sq"[[2]]
        factor1.ss <- drop1(lm(value~factor1+factor2),"factor1")$"Sum of Sq"[[2]]
        factor2.tmp.ss <- factor2.ss <- drop1(lm(value~factor1+factor2),"factor2")$"Sum of Sq"[[2]]
    }
    
    factor1.ms <- factor1.ss / factor1.df
    factor2.ms <- factor2.ss / factor2.df

    mix.ss <- drop1(lm(value~factor1*factor2))$"Sum of Sq"[[2]]
    mix.ms <- mix.ss / mix.df

    total.ss <- (length(value) -1)*var(value)
    res.ss <- total.ss - factor1.tmp.ss - factor2.tmp.ss - mix.ss
    res.ms <- res.ss/res.df

    factor1.f <- factor1.ms/res.ms
    factor2.f <- factor2.ms/res.ms
    mix.f <- mix.ms/res.ms
    
    factor1.p <- 1 - pf(factor1.f,factor1.df,res.df)
    factor2.p <- 1 - pf(factor2.f,factor2.df,res.df)
    mix.p <- 1 - pf(mix.f,mix.df,res.df)

    ANOVA <- matrix(0,4,5)
    ANOVA[,1] <- c(factor1.ss,factor2.ss,mix.ss,res.ss)
    ANOVA[,2] <- c(factor1.df,factor2.df,mix.df,res.df)
    ANOVA[,3] <- c(factor1.ms,factor2.ms,mix.ms,res.ms)
    ANOVA[,4] <- c(factor1.f,factor2.f,mix.f,NA)
    ANOVA[,5] <- c(factor1.p,factor2.p,mix.p,NA)

    colnames(ANOVA) <- c("Sum Sq","df","Mean Sq","F","P")
    rownames(ANOVA) <- c(formula[[3]][[2]],formula[[3]][[3]],"mix","Residual")
    
    result <- list(factor1=summary.factor1,factor2=summary.factor2,mix.num=mix.num,mix.mean=mix.mean,anova=ANOVA)
    class(result) <- c("listMatrix")
    
    if(interaction){
        colors <- c(1:length(levels(factor2)))
        interaction.plot(factor1,factor2,value,col=colors,lwd=2)
    }

    return(result)
}

# Execute simple main effect test and output result of matrix type.
# 
# @param formula formula object.
# @param data data.frame.
# @param anova  listMatrix object(only twoway.factorial.anova() returns)
# @return list matrix.
simple.main.effects <- function(formula, data=NULL, anova){
    make.anova.matrix <- function(f1,f2){
        anova.matrix <- matrix(0,length(levels(f1))+1,5)
        for(x in 1:length(levels(f1))){
            value.tmp <- value[f1==levels(f1)[x]]
            factor.tmp <- f2[f1==levels(f1)[x]]

            factor.N <- tapply(value.tmp,factor.tmp,length)
            factor.mean <- tapply(value.tmp,factor.tmp,mean)
            factor.ss <- sum((factor.mean - mean(value.tmp))^2*factor.N)
            factor.df <- length(factor.mean) - 1
            factor.ms <- factor.ss / factor.df

            f.value <- factor.ms / residual.ms
            p.value <- 1 - pf(f.value,factor.df,residual.df)
            
            anova.matrix[x,] <- c(factor.ss,factor.df,factor.ms,f.value,p.value)
        }
        anova.matrix[length(levels(f1))+1,] <- c(residual.ss,residual.df,residual.ms,NA,NA)
        colnames(anova.matrix) <- c("Sum Sq","df","Mean Sq","F","P")
        rownames(anova.matrix) <- c(levels(f1),"Residual")
        return(anova.matrix)
    }
    if(is.vector(formula) && is.factor(data)){
        dat <-  na.omit(data.frame(formula,data))
        formula <- formula~data
    }else if(is.null(data)){
        dat <-model.frame(formula,na.action=na.omit)
    }else{
        dat <- model.frame(formula,data,na.action=na.omit)
    }
    value <- dat[[1]]
    factor1 <- factor(dat[[3]])
    factor2 <- factor(dat[[2]])

    result.anova <- anova$anova
    residual.ss <- result.anova[nrow(result.anova)]
    residual.df <- result.anova[nrow(result.anova)*2]
    residual.ms <- result.anova[nrow(result.anova)*3]

    anova.matrix1 <- make.anova.matrix(factor1,factor2)
    anova.matrix2 <- make.anova.matrix(factor2,factor1)
    result <- list(anova.matrix1,anova.matrix2)
    names(result) <- c(names(dat[2]),names(dat[3]))

    return(result)
}

# Execute multiple.comparison.test and output result of matrix type.
# 
# @param formula formula object.
# @param data data.frame.
# @param anova listMatrix object(only twoway.factorial.anova() returns)
# @param method char scalar.
# @param fact char scalar.
# @param summary boolean scalar.
# @return list matrix.
multiple.comparison.test <- function(formula, data=NULL, anova=NULL, repeated=F, method=c("Games-Howell","Tukey","bonferroni","holm","BH"), fact=NULL, summary=T){
    method <- match.arg(method)

    if(repeated){
        data <- formula
        origin.data <- data
        data <- na.omit(data)

        row.mean <- rowMeans(data)
        factor.mean <- colMeans(data)

        N <- length(row.mean)
        factor.N <- length(factor.mean)
        missing.Value.N <- length(rowMeans(origin.data)) - N
        df <- N - 1
    
        summary.matrix <- matrix(0, factor.N + 1, 3)
        summary.matrix[,1] <- c(rep(N, factor.N),length(data))
        summary.matrix[,2] <- c(factor.mean, mean(row.mean))
        summary.matrix[,3] <- c(rep(missing.Value.N, factor.N), missing.Value.N*factor.N)
        rownames(summary.matrix) <- c(names(data),"Sum")
        colnames(summary.matrix) <- c("N","mean","Missing Value")

        #t-test
        combn.diff <- combn(factor.N, 2, function(x) apply(data[x],1,diff))
        combn.var <- combn(choose(factor.N,2), 1, function(x) var(combn.diff[,x]))
        var <- mean(combn.var)

        mean.diff <- combn(choose(factor.N,2), 1, function(x) mean(combn.diff[,x]))
        se <- sqrt(var/N)
        t.values <- mean.diff / se

        if(method=="Games-Howell") method <- "Tukey"
        if(method == "Tukey"){
            p.values <- 1 - ptukey(t.values*sqrt(2), factor.N, df*2)
        }
        if(method == "bonferroni" || method == "holm" || method == "BH"){
            p.values <- (1 - pt(abs(t.values),df*2))*2
            p.values <- p.adjust(p.values, method = method)
        }
        factor1.name <- names(data)
    }else{
        if(is.vector(formula) && is.factor(data)){
            origin.dat <- data.frame(formula,data)
            dat <-  na.omit(data.frame(formula,data))
            formula <- formula~data
        }else if(is.null(data)){
            origin.dat <- model.frame(formula,na.action=na.pass)
            dat <-model.frame(formula,na.action=na.omit)
        }else{
            origin.dat <- model.frame(formula,data,na.action=na.pass)
            dat <- model.frame(formula,data,na.action=na.omit)
        }
        origin.value <- origin.dat[[1]]
        origin.factor1 <- factor(origin.dat[[2]])
        
        value <- dat[[1]]
        factor1 <- factor(dat[[2]])

        if(length(dat)>=3){
            origin.factor2 <- factor(origin.dat[[3]])
            factor2 <- factor(dat[[3]])
            if(!is.null(fact)){
                origin.value <- origin.dat[[1]][origin.factor2 == fact]
                origin.factor1 <- origin.dat[[2]][origin.factor2 == fact]
                value <- value[factor2 ==fact]
                factor1 <- factor1[factor2 == fact]
            }
        }

        #summary
        factor1.N <- tapply(value,factor1,length)
        rec.factor.N <- 1/factor1.N
        origin.factor1.N <- tapply(origin.value,origin.factor1,length)
        factor1.mean <- tapply(value,factor1,mean)
        factor1.var <- tapply(value,factor1,var)

        if(length(dat)>=3 && is.null(fact)){
            mix.mean <- tapply(value,list(factor1,factor2),function(x)mean(x))
            factor1.mean <- apply(mix.mean,1,mean)
            cont.N <- table(dat[, 2:3])
            factor2.N <- tapply(value,factor2,length)
            rec.factor.N <- apply(1/cont.N, 1, mean) / length(factor2.N)
        }

        summary.matrix <- matrix(0,length(summary(factor1)) + 1,5)
        summary.matrix[,1] <- c(matrix(factor1.N),length(value))
        summary.matrix[,2] <- c(matrix(factor1.N)-1,length(value)-length(summary(factor1)))
        summary.matrix[,3] <- c(matrix(factor1.mean),mean(value))
        summary.matrix[,4] <- c(matrix(factor1.var),var(value))
        summary.matrix[,5] <- c(matrix(origin.factor1.N - factor1.N),length(origin.value)-length(value))
        rownames(summary.matrix) <- c(names(summary(factor1)),"Sum")
        colnames(summary.matrix) <- c("N","df","mean","u2","Missing Value")

        #multiple comparison test
        if(method=="Tukey" || method=="bonferroni" || method=="holm" || method=="BH"){
            if(!is.null(anova)){
                result.anova <- anova$anova
                residual.ss <- result.anova[nrow(result.anova)]
                residual.df <- result.anova[nrow(result.anova)*2]
                residual.ms <- result.anova[nrow(result.anova)*3]
            }else{
                residual.ss <- sum(tapply(value,factor1,function(x){sum((x -mean(x))^2)}))
                residual.df <- sum(factor1.N -1)
                residual.ms <- residual.ss / residual.df
            }
            tmp <- combn(length(factor1.N),2,function(x){
                mean.diff <- -diff(factor1.mean[x])
                se <- sqrt(sum(residual.ms*rec.factor.N[x]))
                return(c(mean.diff,se))
            })
            mean.diff <- tmp[1,]
            se <- tmp[2,]
            df <- residual.df
            t.values <- abs(mean.diff)/se
        }
        if(method=="bonferroni"){
            p.values <- (1 - pt(t.values, df))*2*choose(length(factor1.N),2)
            p.values[which(p.values > 1)] <- 1
        }else if(method=="holm" || method=="BH"){
            p.values <- (1 - pt(t.values, df))*2
            p.values <- p.adjust(p.values, method = method)
        }else if(method=="Tukey"){
            p.values <- 1 - ptukey(t.values*sqrt(2), length(factor1.N), df)
        }else if(method=="Games-Howell"){
            tmp <- combn(length(factor1.N),2,function(x){
                mean.diff <- -diff(factor1.mean[x])
                df <- sum(factor1.var[x]/factor1.N[x])^2/sum((factor1.var[x]/factor1.N[x])^2/(factor1.N[x]-1))
                se <- sqrt(sum(factor1.var[x]/factor1.N[x]))
                return(c(mean.diff,df,se))
            })
            mean.diff <- tmp[1,]
            df <- tmp[2,]
            se <- tmp[3,]

            t.values <- abs(mean.diff)/se
            p.values <- 1 - ptukey(t.values*sqrt(2), length(factor1.N), df)
        }
        factor1.name <- levels(factor1)
    }
    result <- cbind(mean.diff,se,t.values, df, p.values) 
    colnames(result) <- c("mean diff","SE","t","df","p")
    rownames(result) <- combn(factor1.name,2,paste,collapse=" vs. ")
    
    if(summary){
        result <- list(summary.matrix,result)
        names(result) <- c("summary",method)
    }else{
        result <- list(result)
        names(result) <- c(method)
    }

    class(result) <- c("listMatrix")
    
    return(result)
}

# Execute Levene Test
# 
# @param formula formula object.
# @param data data.frame.
# @param boxplot logical scalar.
# @param col char scalar.
# @return list matrix.
Levene.test <- function(formula,data=NULL,center=mean,boxplot=FALSE,col="darkorange"){
    if(is.vector(formula) && is.factor(data)){
        origin.dat <- data.frame(formula,data)
        dat <-  na.omit(data.frame(formula,data))
        formula <- formula~data
    }

    else if(is.null(data)){
        origin.dat <- model.frame(formula,na.action=na.pass)
        dat <- model.frame(formula,na.action=na.omit)
    }else{
        origin.dat <- model.frame(formula,data,na.action=na.pass)
        dat <- model.frame(formula,data,na.action=na.omit)
    }
    origin.value <- origin.dat[[1]]
    origin.group <- origin.dat[[2]]
    value <- dat[[1]]
    group <- factor(dat[[2]])

    #summary
    group.N <- tapply(value,group,length)
    origin.group.N <- tapply(origin.value,origin.group,length)
    group.mean <- tapply(value,group,mean)
    group.var <- tapply(value,group,var)

    summary <- matrix(0,length(summary(group)) + 1,5)
    summary[,1] <- c(matrix(group.N),length(value))
    summary[,2] <- c(matrix(group.N)-1,length(value)-length(summary(group)))
    summary[,3] <- c(matrix(group.mean),mean(value))
    summary[,4] <- c(matrix(group.var),var(value))
    summary[,5] <- c(matrix(origin.group.N - group.N),length(origin.value)-length(value))
    rownames(summary) <- c(names(summary(group)),"Sum")
    colnames(summary) <- c("N","df","mean","u2","Missing Value")
    if(boxplot){
        boxplot(split(value,group),col=col)
    }

    #Levene
    factor.df <- length(group.mean) - 1
    residual.df <- sum(group.N - 1)

    group.mean.z <- tapply(value,group,function(x){mean(abs(x - center(x)))})
    all.mean.z <- sum(tapply(value,group,function(x){sum(abs(x - center(x)))}))/sum(group.N)

    factor.z.ss <- sum((group.mean.z - all.mean.z)^2*group.N)
    residual.z.ss <- sum(tapply(value,group,function(x){sum((abs(x - center(x))-mean(abs(x - center(x))))^2)}))

    factor.z.var <- factor.z.ss / factor.df
    residual.z.var <-  residual.z.ss / residual.df

    f.value <- factor.z.var / residual.z.var
    p.value <- 1 - pf(f.value,factor.df,residual.df)

    Levene <- matrix(0,2,5)
    Levene[,1] <- c(factor.z.ss,residual.z.ss)
    Levene[,2] <- c(factor.df,residual.df)
    Levene[,3] <- c(factor.z.var,residual.z.var)
    Levene[,4] <- c(f.value,NA)
    Levene[,5] <- c(p.value,NA)

    colnames(Levene) <- c("Sum Sq","df","Mean Sq","F","P")
    rownames(Levene) <- c("Factor","Residual")

    center.str <- paste("center = ", deparse(substitute(center)))

    result <- list(summary=summary,center=center.str,Levene=Levene)
    class(result) <- c("listMatrix")
    
    return(result)
}

# Make crosstable and Execute chi-square test.
# 
# @param row Factor vector.
# @param column Factor vector.
# @param data data.frame.
# @return table.
crosstab <- function(row, column, data=NULL){
    if(is.null(data)){
        cross.table <- table(row,column)
    }else{
        x <- substitute(row)
        y <- substitute(column)
        cross.table <- table(data[[x]],data[[y]])
    }
    names(dimnames(cross.table)) <- NULL
    class(cross.table) <- c("crosstab", "table")
    return(cross.table)
}

summary.crosstab <- function(object,upper=1000){
    upperN <- upper

    #CrossTable
    crossTable <- object
    class(crossTable) <- "table"
    crossTableMarin <- addmargins(crossTable)
    dim <- dim(crossTable)
    row <- dim[1]
    column <- dim[2]
    N <- sum(crossTable)
    df <- (nrow(crossTable) - 1)*(ncol(crossTable) - 1)

   if((df == 1) || min(crossTable)<=5){
        upperN <- 2000
   }

    #ratio
    row.ratio <- prop.table(crossTable,margin=1)
    column.ratio <- prop.table(crossTable,margin=2)

    #standardized residual
    expected.crossTable <- apply(crossTable,2,function(x){sum(x)*apply(crossTable,1,sum)/N})
    res.crossTable <- crossTable - expected.crossTable
    std.res.crossTable <- res.crossTable/sqrt(expected.crossTable)
    
    #Yates' correction
    yates.std.res.crossTable <- (abs(res.crossTable) - 0.5)/sqrt(expected.crossTable)
    yates.std.res.crossTable[which(yates.std.res.crossTable < 0)] <- 0

    #chisq
    chisq.value <- sum(std.res.crossTable^2)
    yates.chisq.value <- sum(yates.std.res.crossTable^2)

    #p
    p.value <- 1 - pchisq(chisq.value,df)
    yates.p.value <- 1 - pchisq(yates.chisq.value,df)
    if(N <= upperN){
        fisher <- fisher.test(crossTable, alternative = "two.sided", workspace=10000000,simulate.p.value=TRUE)
    }

    #chisq-test output
    if(df == 1 && N <= upperN){
        chi.result <- matrix(0,3,3)
    }else if(df == 1 && N > upperN || df > 1 && N <=upperN){
        chi.result <- matrix(0,2,3)
    }else{
        chi.result <- matrix(0,1,3)
    }
    chi.result[1,1] <- chisq.value
    chi.result[1,2] <- df
    chi.result[1,3] <- p.value
    if(df ==1){
        chi.result[2,1] <- yates.chisq.value
        chi.result[2,2] <- df
        chi.result[2,3] <- yates.p.value
        if(N <= upperN){
            chi.result[3,1] <- NA
            chi.result[3,2] <- NA
            chi.result[3,3] <- fisher[[1]]
        }
    }else if(N <= upperN){
        chi.result[2,1] <- NA
        chi.result[2,2] <- NA
        chi.result[2,3] <- fisher[[1]]
    }

    if(df == 1 && N <= upperN){
        rownames(chi.result) <- c("Peason","Yates","Fisher")
    }else if(df ==1 && N > upperN){
        rownames(chi.result) <- c("Peason","Yates")
    }else if(N <= upperN){
        rownames(chi.result) <- c("Peason","Fisher")
    }else{
        rownames(chi.result) <- c("Peason")
    }
    colnames(chi.result) <- c("chi sq","df","P")

    #residual analysis
    se.crossTable <- sqrt(apply(crossTable,2,function(x){(1 - sum(x)/N )*(1 - apply(crossTable,1,sum)/N)}))
    adj.res.crossTable <- std.res.crossTable/se.crossTable
    res.p.crossTable <- (1 - pnorm(abs(adj.res.crossTable)))*2
    res.p.crossTable <- round(res.p.crossTable, 4)

    #Cramer V
    Cramer <- matrix(0,1,1)
    Cramer[1,1] <- sqrt(chisq.value/(N*(min(nrow(crossTable),ncol(crossTable))-1)))
    rownames(Cramer) <- c("Cramer's V")
    
    result <- list(crossTable=crossTable,row.ratio=row.ratio,column.ratio=column.ratio,chisq.test=chi.result,residualAnalysis=adj.res.crossTable,res.p.value=res.p.crossTable,Cramer=Cramer)
    class(result) <- c("listMatrix")
    
    return(result)
}

# Correlation matrix or partial correlation matrix and t-test.
# 
# @param data numeric matrix
# @param control char vector
# @param use  pairwise or complete
# @param method  peason, kendall, spearman
# @return list matrix
correl <- function (data, control = NULL, use = c("pairwise","complete"),  method = c("pearson", "kendall", "spearman")) {
    use <- match.arg(use)
    method <- match.arg(method)
    n <- t(!is.na(data)) %*% (!is.na(data))

    if ((use == "complete") | (min(n) == max(n))){
        n <- min(n)
    }

    if(is.null(control)){
        r <- cor(data, use = use, method = method)
        t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
        t[t == Inf] <- NA
        p <- 2 * (1 - pt(abs(t), (n - 2)))
        p[p > 1] <- 1
        result <- list(n = n, r = r, t = t, p = p)
    }else{
        control.data <- data[,(colnames(data) %in% control)]
        control.labels <-control
        dat <- data[, !(colnames(data) %in% control)]
        var.labels <- colnames(dat)
        numVar <- length(dat)
        numAllVar <- numVar + length(control)
        data <- cbind(dat,control.data)
        numData <- dim(data)[1]

        r <- cor(data, use = use, method = method)
        t.matrix <- matrix(0, ncol = numAllVar, nrow = numAllVar)
        for (i in 1:numAllVar) {
            t.matrix[i, i] <- 1
        }
        r.matrix <- t.matrix %*% r %*% t(t.matrix)
        r.matrix[abs(r.matrix) > 1] <- NA
    
        var.r.matrix <- r.matrix[1:numVar, 1:numVar]
        control.var.r.matrix <- r.matrix[1:numVar, (numVar + 1):numAllVar]
        control.r.matrix <- r.matrix[(numVar + 1):numAllVar, (numVar + 1):numAllVar]
        control.r.inv <- solve(control.r.matrix)
    
        partial.r <- var.r.matrix - control.var.r.matrix %*% control.r.inv %*% t(control.var.r.matrix)
        partial.r <- cov2cor(partial.r)
        colnames(partial.r) <- rownames(partial.r) <- var.labels

        t <- partial.r/sqrt((1 - partial.r^2)/(numData - numAllVar + 2))
        t[t == Inf] <- NA
        p <- 2 * (1 - pt(abs(t), (numData - numAllVar + 2)))
        colnames(p) <- rownames(p) <- var.labels

        result <- list(n = n, control = control.labels, partial.r = partial.r, t = t, p = p)
    }
    class(result) <- c("listMatrix")
    return(result)
}

# Execute linear regression analysis and output result of matrix type.
# 
# @param formula formula object.
# @param data data.frame.
# @param direction  c("both", "backward", "forward"), .
# @param step logical scalar.
# @return list matrix.
linest <- function(formula, data=NULL, direction=c("enter", "both", "backward", "forward"), steps=FALSE){
    direction <- match.arg(direction)

    if(is.null(data)){
        data <- model.frame(formula)
    }else{
        data <- model.frame(formula,data)
    }
    attach(data)
    model <- lm(formula,data)
    out <- summary(model)

    if(steps || direction !="enter"){
        if(direction == "enter") direction <- "both"
        model <- step(model,data, direction=direction)
        out <- summary(model)
    }
    detach(data)
    
    coefficients <- data.frame(out$coefficients)
    label <- rownames(coefficients)

    #summary
    R <- sqrt(out$r.squared)
    R2 <- out$r.squared
    adjustedR2 <- out$adj.r.squared
    SE <- out$sigma
    summary <- matrix(0,4,1)
    summary[,1] <- c(R[[1]],R2[[1]],adjustedR2[[1]],SE[[1]])
    colnames(summary) <- c("Model Summary")
    rownames(summary) <- c("R","R2","Adjusted R2","SE")

    #ANOVA
    anova.table <- anova(model)
    df.R <- out$fstatistic[[2]]
    df.e <- out$fstatistic[[3]]
    df.T <- df.R + df.e

    SS.R <- sum(anova.table$`Sum Sq`[-length(anova.table$Df)])
    SS.e <- anova.table$`Sum Sq`[length(anova.table$Df)]
    SS.T <- sum(anova.table$`Sum Sq`)

    MS.R <- SS.R/df.R
    MS.e <- SS.e/df.e
    MS.T <- SS.T/df.T

    F <- out$fstatistic[[1]]
    anova.p <- 1 - pf(F,df.R,df.e)

    anova.matrix <- matrix(NA,3,5)
    anova.matrix[,1] <- c(SS.R,SS.e,SS.T)
    anova.matrix[,2] <- c(df.R,df.e,df.T)
    anova.matrix[,3] <- c(MS.R,MS.e,MS.T)
    anova.matrix[1,4] <- F
    anova.matrix[1,5] <- anova.p
    colnames(anova.matrix) <- c("SS","df","MS","F","P")
    rownames(anova.matrix) <- c("Regression","Residual","Total")

    #Coefficient
    num.org.terms <- length(labels(terms(model)))
    num.terms <- length(coefficients$Estimate) - 1
    B <- coefficients$Estimate
    SE2 <- coefficients$Std..Error
    beta <- matrix(NA, num.terms + 1, 1)
    uy <- sd(model$model[[1]])
    for(term.index in 2:num.terms){
        beta[term.index] <- sd(model.matrix(model)[,term.index])*B[term.index]/uy
    }
    t <- coefficients$t.value
    P <- coefficients$Pr...t..
    F <- t^2
    if(num.terms >= 2){
        partial.r <- cov2cor(vcov(model)[-1, -1])
        det.partial.r <- det(partial.r)
        assign <- 1:(num.terms)
        vif <- matrix(NA, num.terms, 1)
        for (term in 1:num.terms) {
            subs <- which(assign == term)
            vif[term, 1] <- det(as.matrix(partial.r[subs, subs])) * det(as.matrix(partial.r[-subs, -subs]))/det.partial.r
        }
    }else{
        vif <- matrix(NA, num.terms, 1)
    }
    VIF <- c(NA,vif)
    coefficient.matrix <- matrix(0,df.R+1,7)
    coefficient.matrix[,1] <- B
    coefficient.matrix[,2] <- SE2
    coefficient.matrix[,3] <- beta
    coefficient.matrix[,4] <- t
    coefficient.matrix[,5] <- P
    coefficient.matrix[,6] <- F
    coefficient.matrix[,7] <- VIF
    colnames(coefficient.matrix) <- c("b","SE","beta","t","P","F","VIF")
    rownames(coefficient.matrix) <- rownames(coefficients)
    
    result <- list(summary=summary,anova=anova.matrix,coefficient=coefficient.matrix)
    class(result) <- c("listMatrix")

    return(result)
}

# Execute random sampling and plot distribution of mean.
# 
# @param x numeric vector.
# @param N numeric scalar.
# @param num numeric scalar.
# @param col char scalar
samplingTest <- function(x,N=5,num=2500,col="royalblue"){
    x <- na.omit(x)
    if(N > length(x)) N <- length(x)
    mean.value <- matrix(0,num,1)
    for(i in 1:num){
      mean.value[i] <- mean(sample(x,N,replace=FALSE))
    }

    start <- min(mean.value)
    width <- ceiling(max(mean.value) - min(mean.value))/(ceiling(log(N)+1)*5+20)
    end <- max(mean.value)+width
    cut <- seq(start,end,width)
    hist(mean.value,cut,right=FALSE,freq=FALSE,col=col,main=paste("Distribution(N =",N,")"))
    curve(dnorm(x,mean=mean(mean.value),sd=sd(mean.value)),lwd=3,add=TRUE)
}

# Convert crosstable as matrix to weighted data.
# 
# @param matrix matrix
# @return data.frame
convert.crosstable.to.weighted.data <- function(matrix){
    vec <- vector()
    for(i in 1:ncol(matrix)){
        vec <- c(vec,rep(colnames(matrix)[i],nrow(matrix)))
    }
    result <- data.frame(var1=rep(rownames(matrix),ncol(matrix)),var2=vec,Freq=as.vector(matrix))

    return(result)
}

# Convert crosstable as data.frame to weighted data.
# 
# @param data.frame data.frame.
# @return data.frame
convert.data.frame.crosstable.to.weighted.data <- function(data.frame){
    row.names <- data.frame[[1]]
    data.rm.row.names <- data.frame[-1]
    matrix.data <- as.matrix(data.rm.row.names)
    row.names(matrix.data) <- row.names
    
    result <- data.frame(as.table(matrix.data))
    
    return(result)
}

# Convert weighted data to unweighted data.
# 
# @param data.frame data.frame.
# @param freq.name char scaler.
# @return data.frame
convert.weighted.data.to.original <- function(data.frame,freq.name="Freq"){
    data.frame <- na.omit(data.frame)
    freq <-  data.frame[colnames(data.frame) == freq.name][[1]]
    
    result <- data.frame(lapply(data.frame,function(listData){rep(listData,freq)}))[, colnames(data.frame) != freq.name]
    
    invisible(result)
}

# Record dummy variables from factor variables.
# 
# @param data.frame data.frame.
# @param ref char vector.
# @return data.frame
recode.dummy.variables <- function(data, ref=NULL){
    col.names <- colnames(data)
    result <- NULL
    for(col.num in 1:ncol(data)){
        recode.col.names <- NULL
        origin.col <- data[,col.names[col.num]]
        if(is.factor(origin.col)) {
            level.values <- levels(origin.col)
            level.n <- length(level.values)
            dummy.col <- matrix(0, length(origin.col), level.n)
            for(dummy.col.num in 1:level.n) dummy.col[which(level.values[dummy.col.num] == origin.col),dummy.col.num] <- 1
            recode.col.names <- paste(col.names[col.num], level.values, sep = ".")
            colnames(dummy.col) <- as.vector(recode.col.names)
            if(any(colnames(dummy.col) %in% ref)){
                dummy.col <- as.matrix(dummy.col[, !(colnames(dummy.col) %in% ref)])
                colnames(dummy.col) <- as.vector(recode.col.names[!recode.col.names %in% ref])
            }else if(is.null(ref)||ref !="none"){
                dummy.col <- as.matrix(dummy.col[,-1])
                colnames(dummy.col) <- as.vector(recode.col.names[-1])
            }
        }else{
            dummy.col <- as.matrix(origin.col)
            recode.col.names <- append(recode.col.names,col.names[col.num])
            colnames(dummy.col) <- as.vector(recode.col.names)
        }
        result <- cbind(result, dummy.col)
    }
    return(data.frame(result))
}

# write output data to csv.
# 
# @param listMatrix.
# @param filename char scaler.
write.output <- function(object, filename="output.csv"){
    result <- object
    class(result) <- "list"

    result[[1]]  <- as.matrix(result[[1]])
    suppressWarnings(write.table(result[[1]], filename, sep=",", row.names=T, col.names=NA, append=F))
    if(length(result) >= 2){
        for(i in 2:length(result)){
            result[[i]]  <- as.matrix(result[[i]])
            suppressWarnings(write.table(result[[i]], filename, sep=",", row.names=T, col.names=NA, append=T))
        }
    }
}