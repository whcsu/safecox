rm(list = ls())
library(glmnet)
library(survival)
library(MASS)
library(ggplot2)
library(reshape2)
library(survcomp)
#library(BiocManager)
library(microbenchmark)
####get data
get.data<- function(n,p,rho,q,c0){# n sample size; p dimension.
  beta <- rep(0,p);
  beta[1:q] <- c(rep(c(3,-3),q/2));
  sig <- matrix(0,p,p);
  sig <- rho^abs(row(sig)-col(sig));
  #sig <- matrix(rho,p, p);
  diag(sig)<-rep(1,p);
  status <- 0;
  while(status==0) 
  {
    Z_temp <- mvrnorm(5*n,rep(0,p),sig);
    haz <- drop(1+Z_temp%*%beta);
    ind <- (1:length(haz))[haz>0][1:n];#取haz中大于0的指标且扩充到n
    if(!is.na(sum(haz[ind])))
    {
      status <- 1;
    }
  }
  Z <- Z_temp[ind,];
  T <- rexp(n,rate=haz[ind]);
#c0 controls the censoring rate, the bigger c0,the smaller ratio
  #c0=3 produces approximately 25% censoring rate
  C <- runif(n,min=0,max=c0);
  de <- as.numeric(T<=C);
  X <- pmin(T,C);
  Y <- Surv(X,de);
  
  return(list(Z=Z,beta=beta,de=de,X=X,Y=Y));
}
lambda.sk<-function(x,y,data,k,g){
  p=ncol(x)
  y<-as.matrix(y)
  Q = t(x)%*%x
  B = t(x)%*%y
  f<-rep(0,p)
  tmpbeta<-rep(0,p)
  ss2 <- function(j,tmpbeta,Q,B)
  {
    a <- sum(tmpbeta*Q[,j])-tmpbeta[j]*Q[j,j]
    s <- 2*(a-B[j])
    return(abs(s))
  }#F0
  for(j in 1:p){
    f[j]<-ss2(j,tmpbeta,Q,B)
    #f[j]<-ss2(j,tmpbeta,Q,B)*abs(beta.ini[j])
  }
  max_lambda<-max(f)
  lambda_k<-rep(0,k)
  for(i in 1:k){
    min_lambda<-g*max_lambda
    lambda_k[i]<-max_lambda*(min_lambda/max_lambda)^(i/k)
  }
  return(list(lambda_k=lambda_k,max_lambda=max_lambda))
}
SER <- function(n,p,x,y,lambda,beta.ini){
  ####beta ini####
  #if(n<p){
  #fit<-cv.glmnet(x,Surv(y$time,y$status),family = "cox",alpha = 0)
  #beta.ini=coef(fit,s=fit$lambda.min)#岭回归 s=fit$lambda.min
  #}
  # else{beta.ini = coxph(Surv(y$time,y$status)~x)$coef}
  m=0#删除变量个数计数
  sf<-rep(0,p)
  y<-data.frame(y)
  x<-data.frame(x)
  ck = apply(x[which(y$status==1),],2,FUN=sum)#所有失效时间变量加和(p维向量)
  s1=rep(0,p)
  s2=rep(0,p)
  for(i in 1:n){
    minx=apply(x[i:n,],2,FUN =min)
    s1=s1+minx
    maxx=apply(x[i:n,],2,FUN=max)
    s2=s2+maxx}
  s=rbind(ck-s1,s2-ck)
  sf=apply(s,2,max)#最大值
  #safe rule
  lams=lambda/abs(beta.ini)
  index<-rep(0,p)#删除变量索引(1表示删除)
  for(i in 1:p){
    if(!is.na(sf[i]))
    {if(sf[i]<lams[i]){index[i]=1;m=m+1} 
    }  
  }
  x=x[,which(index==0)]
  return(list(index=index,z=x,m=m))
} 
eff_safe<-function(x,y,k,p,m,lambda_k,beta.ini){
  #fit<-cv.glmnet(x,Surv(y$time,y$status),family = "cox",alpha = 0,standardize = TRUE)
  #beta.ini=coef(fit,s=fit$lambda.min)
  coxadalafit=cv.glmnet(x,Surv(y$time,y$status),family="cox",alpha=1,penalty.factor = 1/abs(beta.ini),standardize = TRUE)
  rej<-rep(0,k)
  scre<-rep(0,k)
  tr<-rep(0,k)
  for(i in 1:k){
    print(i)
    tr[i]<-p-coef(coxadalafit,s=lambda_k[i])@p[2]
    scre[i]<-(p-m[i])/p
    if(tr[i]!=0){rej[i]<-m[i]/tr[i] }
    else{rej[i]=1}
  }
  return(list(rej=rej,scre=scre,tr=tr))
}
coxadala.perf=function(data,Rn,totalfold,lambda){
  n=dim(data)[1]
  testresult <-data.frame(row.names=1:Rn)    
  ci_coxlasso<-c(rep(0,Rn))
  ibs_coxlasso<-c(rep(0,Rn))
  for (ii in seq(1, Rn, by=totalfold)){
    print(ii)
    mydata<-as.data.frame(data[sample(nrow(data)),])
    y1<-mydata[,1:2]
    names(y1)<-c("time","status")
    mydata<-data.frame(y1,mydata[,-c(1:2)])
    #Create 5 equally size folds
    folds <- cut(seq(1,nrow(mydata)),breaks=totalfold,labels=FALSE)
    #Perform 2 fold cross validation
    for(k in 1:totalfold){
      #Segement your data by fold using the which() function
      testIndexes <- which(folds==k,arr.ind=TRUE)
      teset <- mydata[testIndexes, ]
      trset <- mydata[-testIndexes, ]
      #Use the test and train data partitions however you desire...
      tryCatch({
        coxridgefit=cv.glmnet(x=as.matrix(trset[,-c(1,2)]),y=Surv(trset$time, trset$status),family="cox",alpha=0)
        beta.ini<- coef(coxridgefit, s = coxridgefit$lambda.min)
        coxadalafit=cv.glmnet(x=as.matrix(trset[,-c(1,2)]),y=Surv(trset$time, trset$status),family="cox",alpha=1,penalty.factor = 1/abs(beta.ini))
        coxglmpre=predict(coxadalafit,as.matrix(teset[,-c(1,2)]), type="response",s=lambda)
        ci_coxlasso[ii+k-1]=unlist(concordance.index(coxglmpre,teset$time,teset$status)[1])
        
        if(!is.na(ci_coxlasso[ii+k-1])){
          #cal ibs
          uniquetimes=sort(unique(trset$time))
          coxglmpre2=as.numeric(predict(coxadalafit,as.matrix(teset[,-c(1,2)]),s=lambda, type="link"))
          basesurvresult <- basesurv(Surv(teset$time,teset$status), coxglmpre2,uniquetimes)
          cox_surv_prob  <- exp(exp(coxglmpre2) %*% -t(basesurvresult$cumBaseHaz))
          ibs_coxlasso[ii+k-1]=isbrier(Surv(teset$time,teset$status),t(cox_surv_prob),uniquetimes)
          
        }
        if (is.na(ci_coxlasso[ii+k-1]))
          ibs_coxlasso[ii+k-1]=NA
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
  }
  ci.m=mean(ci_coxlasso,na.rm=T)
  ibs.m=mean(ibs_coxlasso,na.rm=T)
  return(list(ci=ci_coxlasso,ibs=ibs_coxlasso,ci.m=ci.m,ibs.m=ibs.m))
}
isbrier <-function(obj, pred, btime = range(obj[,1])){
  if(!inherits(obj, "Surv"))
    stop("obj is not of class Surv")
  
  # check for right censoring
  
  # <FIXME>
  class(obj) <- NULL
  # </FIXME>
  if (attr(obj, "type") != "right")
    stop("only right-censoring allowed")
  N <- nrow(obj)	
  
  # get the times and censoring of the data, order them with resp. to time
  
  time <- obj[,1]
  ot <- order(time)
  cens <- obj[ot,2]
  time <- time[ot]
  
  # get the times to compute the (integrated) Brier score over
  
  if (is.null(btime)) stop("btime not given")
  if (length(btime) < 1) stop("btime not given")
  
  if (length(btime) == 2) {
    if (btime[1] < min(time)) warning("btime[1] is smaller than min(time)")
    if (btime[2] > max(time)) warning("btime[2] is larger than max(time)")
    btime <- time[time >= btime[1] & time <=
                    btime[2]]
  }
  
  ptype <- class(pred)
  # <begin> S3 workaround
  if (is.null(ptype)) {
    if (is.vector(pred)) ptype <- "vector"
    if (is.list(pred)) ptype <- "list"
  }
  # <end>
  if (ptype == "numeric" && is.vector(pred)) ptype <- "vector"
  
  survs <- NULL
  switch(ptype, survfit = {
    survs <- getsurvprob(pred, btime)
    survs <- matrix(rep(survs, N), nrow=length(btime))
  }, list = {
    if (!inherits(pred[[1]], "survfit")) stop("pred is not a list of survfit objects") 
    if (length(pred) != N) stop("pred must be of length(time)")
    pred <- pred[ot]
    survs <-  matrix(unlist(lapply(pred, getsurvprob, times = btime)),
                     nrow=length(btime), ncol=N)
  }, vector = {
    if (length(pred) != N) stop("pred must be of length(time)")
    if (length(btime) != 1) stop("cannot compute integrated Brier score with pred")
    survs <- pred[ot]
  }, matrix = {
    # <FIXME>
    if (all(dim(pred) == c(length(btime), N)))
      survs <- pred[,ot]
    else
      stop("wrong dimensions of pred")
    # </FIXME>
  })
  if (is.null(survs)) stop("unknown type of pred")
  
  # reverse Kaplan-Meier: estimate censoring distribution
  
  ### deal with ties
  hatcdist <- prodlim(Surv(time, cens) ~ 1,reverse = TRUE)
  csurv <- predict(hatcdist, times = time, type = "surv")
  
  if(length(csurv)<length(btime)){
    #adding more zero here
    addsurv=rep(0,length(btime)-length(csurv))
    csurv=c(csurv,addsurv)
  }
  
  csurv[csurv == 0] <- Inf
  # hatcdist <- survfit(Surv(time, 1 - cens) ~ 1)
  # csurv <- getsurv(hatcdist, time)
  # csurv[csurv == 0] <- Inf
  bsc <- rep(0, length(btime))
  
  # compute Lebesque-integrated Brier score
  
  if (length(btime) > 1) {
    for (j in 1:length(btime)) {
      help1 <- as.integer(time <= btime[j] & cens == 1)
      help2 <- as.integer(time > btime[j])
      bsc[j] <-  mean((0 - survs[j,])^2*help1*(1/csurv) +
                        (1-survs[j,])^2*help2*(1/csurv[j]))
      #  print(bsc[j])
      # print(j)
    }
    
    ### apply trapezoid rule
    idx <- 2:length(btime)
    RET <- diff(btime) %*% ((bsc[idx - 1] + bsc[idx]) / 2)
    RET <- RET / diff(range(btime))
    
    ### previously was
    #diffs <- c(btime[1], btime[2:length(btime)] -
    #                     btime[1:(length(btime)-1)])
    #RET <- sum(diffs*bsc)/max(btime)
    names(RET) <- "integrated Brier score"
    attr(RET, "time") <- range(btime)
    
    # compute Brier score at one single time `btime'
    
  } else {
    help1 <- as.integer(time <= btime & cens == 1)
    help2 <- as.integer(time > btime)
    cs <- predict(hatcdist, times=btime, type = "surv")
    ### cs <- getsurv(hatcdist, btime)
    if (cs == 0) cs <- Inf
    RET <-  mean((0 - survs)^2*help1*(1/csurv) +
                   (1-survs)^2*help2*(1/cs))
    names(RET) <- "Brier score"
    attr(RET, "time") <- btime
  }
  RET
}
basesurv <- function (response, lp, times.eval = NULL, centered = FALSE)
{
  if (is.null(times.eval)) times.eval <- sort(unique(response[,1]))
  
  t.unique <- sort(unique(response[,1][response[,2] == 1]))
  alpha    <- length(t.unique)
  
  for (i in 1:length(t.unique)) {
    alpha[i] <- sum(response[,1][response[,2] == 1] == t.unique[i])/sum(exp(lp[response[,1] >=  t.unique[i]]))
  }
  
  obj   <- approx(t.unique, cumsum(alpha), yleft=0, xout = times.eval, rule=2)
  
  if (centered) obj$y <- obj$y * exp(mean(lp))
  obj$z <- exp(-obj$y)
  
  names(obj) <- c("times","cumBaseHaz","BaseSurv")
  return(obj)
}

set.seed(123)
n=200;p=100;k=100;q=10;c0=1;rho=0.6
d<-get.data(n,p,rho,q,c0)  
ratio<-1-sum(d$de)/n #censoring ratio 30+%
print(ratio)
data<-data.frame(d$X,d$de,d$Z)
data<-data[order(data$d.X),]
y<-data[,1:2]
names(y)<-c("time","status")
x1<-as.matrix(data[,-c(1:2)] )
x<-scale(x1);data<-data.frame(y,x)
write.csv(data,file="S1.csv")
#ridge for initial beta
fit_ridge<-cv.glmnet(x,Surv(y$time,y$status),family = "cox",alpha = 0,standardize = TRUE)
beta.ini=coef(fit_ridge,s=fit_ridge$lambda.min)#岭回归

#screening for 100 lambda 
lam<-lambda.sk(x,y,data,k,0.01)
index<-matrix(0,ncol=k,nrow=p)
m<-rep(0,k)
for(i in 1:k){
  fit<-SER(n,p,x,y,lam$lambda_k[i],beta.ini)
  index[,i]<-fit$index
  m[i]<-fit$m
} 


####plot eff with screen ratio and reject ratio###

eff<-eff_safe(x,y,k,p,m,lam$lambda_k,beta.ini)

RR.data<-data.frame(lam$lambda_k/lam$max_lambda,eff$rej)
names(RR.data)<-c("lam","Rejection Ratio")
md<- melt(RR.data, id.vars="lam")
RR_plot<-ggplot(data=md,aes(x=lam,y=value))+
  scale_x_log10()+
  geom_line(aes(col=variable),size=1.2)+
  labs(x=expression(frac(lambda,lambda[max])),y="Rejection Ratio")+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=16),axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),legend.position = "non")#去除图例
SR.data<-data.frame(lam$lambda_k/lam$max_lambda,eff$scre)
names(SR.data)<-c("lam","Screen Ratio")
md<- melt(SR.data, id.vars="lam")
SR_plot<-ggplot(data=md,aes(x=lam,y=value))+
  scale_x_log10()+
  geom_line(aes(col=variable),size=1.2)+
  labs(x=expression(frac(lambda,lambda[max])),y="Screen Ratio")+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=16),axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),legend.position = "non")#去除图例
SR_plot
RR_plot


###choose the special lambda
Rn=10
totalfold=2
p1<-matrix(0,nrow=k,ncol=2)
ps1<-matrix(0,nrow=k,ncol=2)

for(i in 1:100){
  print(i)
  if(m[i]<=p-2)
  {
    cat("没有删完的lambda:",i)
    ss.data<-data.frame(y,x[,which(index[,i]!=1)])
    ci<-coxadala.perf(data,10,2,lam$lambda_k[i])$ci.m
    ib<-coxadala.perf(data,10,2,lam$lambda_k[i])$ibs.m
    p1[i,]<-c(ci,ib)
    ci_s<-coxadala.perf(ss.data,10,2,lam$lambda_k[i])$ci.m
    ib_s<-coxadala.perf(ss.data,10,2,lam$lambda_k[i])$ibs.m
    ps1[i,]<-c(ci_s,ib_s)
  }
}

lam_best_index<-51
lambda.best<-lam$lambda_k[lam_best_index]

##comparsion boxplot for three model with c-index and b-score
Rn=100#run 100 times
totalfold=2
ci_coxlasso<-c(rep(0,Rn))
ci_coxadaptive_lasso<-c(rep(0,Rn))
ci_coxsafe_adaptive<-c(rep(0,Rn))
ibs_coxlasso<-c(rep(0,Rn))
ibs_coxadaptive_lasso<-c(rep(0,Rn))
ibs_coxsafe_adaptive<-c(rep(0,Rn))

s.data<-data.frame(y,x[,which(index[,lam_best_index]!=1)])#after safe data

for (ii in seq(1, Rn, by=totalfold)){
  print(ii)
  data<-data[sample(nrow(data)),]
  s.data<-s.data[sample(nrow(s.data)),]
  folds <- cut(seq(1,nrow(data)),breaks=totalfold,labels=FALSE)
  for(k in 1:totalfold){
    testIndexes <- which(folds==k,arr.ind=TRUE)
    y1<-data[,1:2]
    y2<-s.data[,1:2]
    names(y1)<-c("time","status")
    names(y2)<-c("time","status")
    data<-data.frame(y1,data[,-c(1:2)])
    ss.data<-data.frame(y2,s.data[,-c(1:2)])
    teset <- data[testIndexes, ]
    trset <- data[-testIndexes, ]
    s.teset <- ss.data[testIndexes, ]
    s.trset <- ss.data[-testIndexes, ]
    print("coxlass")
    tryCatch({
      coxglmfit=cv.glmnet(x=as.matrix(trset[,-c(1:2)]),y=Surv(trset$time, trset$status),family="cox",alpha=1)
      coxglmpre=predict(coxglmfit,as.matrix(teset[,-c(1:2)]), type="response",s=coxglmfit$lambda.min)
      ci_coxlasso[ii+k-1]=unlist(concordance.index(coxglmpre,teset$time,teset$status)[1])#cox c 指数
      if(!is.na(ci_coxlasso[ii+k-1])){
        #cal ibs
        uniquetimes=sort(unique(trset$time))
        coxglmpre2=as.numeric(predict(coxglmfit,as.matrix(teset[,-c(1:2)]),s=lambda.best, type="link"))
        basesurvresult <- basesurv(Surv(teset$time,teset$status), coxglmpre2,uniquetimes)
        cox_surv_prob <- exp(exp(coxglmpre2) %*% -t(basesurvresult$cumBaseHaz))
        ibs_coxlasso[ii+k-1]=isbrier(Surv(teset$time,teset$status),t(cox_surv_prob),uniquetimes)
      }
      
      if (is.na(ci_coxlasso[ii+k-1]))
        ibs_coxlasso[ii+k-1]=NA
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    print("coxadaptive_lasso")
    tryCatch({
      coxridgefit=cv.glmnet(x=as.matrix(trset[,-c(1,2)]),y=Surv(trset$time, trset$status),family="cox",alpha=0)
      beta.ini<- coef(coxridgefit, s = coxridgefit$lambda.min)
      coxadalafit=cv.glmnet(x=as.matrix(trset[,-c(1,2)]),y=Surv(trset$time, trset$status),family="cox",alpha=1,penalty.factor = 1/abs(beta.ini))
      coxglmpre=predict(coxadalafit,as.matrix(teset[,-c(1,2)]), type="response",s=lambda.best)
      ci_coxadaptive_lasso[ii+k-1]=unlist(concordance.index(coxglmpre,teset$time,teset$status)[1])
      
      if(!is.na(ci_coxadaptive_lasso[ii+k-1])){
        #cal ibs
        uniquetimes=sort(unique(trset$time))
        coxglmpre2=as.numeric(predict(coxadalafit,as.matrix(teset[,-c(1,2)]),s=lambda.best, type="link"))
        basesurvresult <- basesurv(Surv(teset$time,teset$status), coxglmpre2,uniquetimes)
        cox_surv_prob  <- exp(exp(coxglmpre2) %*% -t(basesurvresult$cumBaseHaz))
        ibs_coxadaptive_lasso[ii+k-1]=isbrier(Surv(teset$time,teset$status),t(cox_surv_prob),uniquetimes)
        
      }
      if (is.na(ci_coxadaptive_lasso[ii+k-1]))
        ibs_coxadaptive_lasso[ii+k-1]=NA
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    print("coxsafe_adaptive")
    #logrank
    tryCatch({	
      coxridgefit=cv.glmnet(x=as.matrix(s.trset[,-c(1,2)]),y=Surv(s.trset$time, s.trset$status),family="cox",alpha=0)
      beta.ini<- coef(coxridgefit, s = coxridgefit$lambda.min)
      coxadalafit=cv.glmnet(x=as.matrix(s.trset[,-c(1,2)]),y=Surv(s.trset$time, s.trset$status),family="cox",alpha=1,penalty.factor = 1/abs(beta.ini))
      coxglmpre=predict(coxadalafit,as.matrix(s.teset[,-c(1,2)]), type="response",s=lambda.best)
      ci_coxsafe_adaptive[ii+k-1]=unlist(concordance.index(coxglmpre,s.teset$time,s.teset$status)[1])
      
      if(!is.na(ci_coxsafe_adaptive[ii+k-1])){
        #cal ibs
        uniquetimes=sort(unique(trset$time))
        coxglmpre2=as.numeric(predict(coxadalafit,as.matrix(s.teset[,-c(1,2)]),s=lambda.best, type="link"))
        basesurvresult <- basesurv(Surv(s.teset$time,s.teset$status), coxglmpre2,uniquetimes)
        cox_surv_prob  <- exp(exp(coxglmpre2) %*% -t(basesurvresult$cumBaseHaz))
        ibs_coxsafe_adaptive[ii+k-1]=isbrier(Surv(s.teset$time,s.teset$status),t(cox_surv_prob),uniquetimes)
        
      }
      if (is.na(ci_coxsafe_adaptive[ii+k-1]))
        ibs_coxsafe_adaptive[ii+k-1]=NA
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
}

##c-index 

CI<-c(ci_coxlasso,ci_coxsafe_adaptive,ci_coxadaptive_lasso)
CMethods<-c(rep("Lasso_Cox",100),rep("SAFE_aLasso_Cox",100),rep("aLasso_Cox",100))
CMethod<- factor(CMethods, levels=c("SAFE_aLasso_Cox","aLasso_Cox","Lasso_Cox"), ordered=TRUE)
xdataC<-data.frame(CI,CMethod)
##
mytheme<- theme(axis.title=element_text(face="plain",size=12, color="black"),
                  axis.text=element_text(face="plain", size=12,color="black"), 
                  panel.background=element_rect(color="white"),
                  panel.grid.minor.y=element_blank(),
                  panel.grid.minor.x=element_blank())

boxplot_CI<-ggplot(data=xdataC,aes(x=CMethod,y=CI,fill=CMethods))+
  geom_boxplot(width=0.5,outlier.colour="black", outlier.shape=1.5, outlier.size=3)+
  labs(title="",x="",y="C-index")+theme(legend.position="none")+mytheme+
  scale_fill_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  theme(axis.title.y=element_text(size=16),axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),legend.position="none")

###b-score
IBS<-c(ibs_coxlasso,ibs_coxsafe_adaptive,ibs_coxadaptive_lasso)
IBMethods<-c(rep("Lasso_Cox",100),rep("SAFE_aLasso_Cox",100),rep("aLasso_Cox",100))
IBMethod<- factor(CMethods, levels=c("SAFE_aLasso_Cox","aLasso_Cox","Lasso_Cox"), ordered=TRUE)
xdataIB<-data.frame(IBS,IBMethods)

boxplot_IBS<-ggplot(data=xdataIB,aes(x=IBMethod,y=IBS,fill=IBMethods))+
  geom_boxplot(width=0.5,outlier.colour="black", outlier.shape=1.5, outlier.size=3)+
  labs(title="",x="",y="B-score")+theme(legend.position="none")+mytheme+
  scale_fill_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  theme(axis.title.y=element_text(size=16),axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),legend.position="none")

boxplot_CI
boxplot_IBS

#####time
Rn=10
totalfold=2
for (ii in seq(1, Rn, by=totalfold)){
  print(ii)
  data<-data[sample(nrow(data)),]
  s.data<-s.data[sample(nrow(s.data)),]
  folds <- cut(seq(1,nrow(data)),breaks=totalfold,labels=FALSE)
  for(k in 1:totalfold){
    print(k)
    testIndexes <- which(folds==k,arr.ind=TRUE)
    y1<-data[,1:2]
    y2<-s.data[,1:2]
    names(y1)<-c("time","status")
    names(y2)<-c("time","status")
    data<-data.frame(y1,data[,-c(1:2)])
    ss.data<-data.frame(y2,s.data[,-c(1:2)])
    teset <- data[testIndexes, ]
    trset <- data[-testIndexes, ]
    s.teset <- ss.data[testIndexes, ]
    s.trset <- ss.data[-testIndexes, ]
    mbm<- microbenchmark("lasso_Cox"={
      coxglmfit=cv.glmnet(x=as.matrix(trset[,-c(1:2)]),y=Surv(trset$time, trset$status),family="cox",alpha=1)
      coxglmpre=predict(coxglmfit,as.matrix(teset[,-c(1:2)]), type="response",s=lambda.best)
    },
    "alasso_Cox"={
      ridge_cv <- cv.glmnet(x=as.matrix(trset[,-c(1,2)]),y=Surv(trset$time, trset$status),type.measure = "deviance",nfold = 10,family = "cox",alpha = 0)
      best_ridge_coef <- coef(ridge_cv, s = ridge_cv$lambda.min)
      coxglmfit=cv.glmnet(x=as.matrix(trset[,-c(1,2)]),Surv(trset$time, trset$status),type.measure = "deviance",nfold = 10,family = "cox",alpha = 1,penalty.factor = 1 / abs(best_ridge_coef))
      coxglmpre=predict(coxglmfit,as.matrix(teset[,-c(1,2)]), type="response",s=lambda.best)
    },
    "safe_alasso_Cox"={
      ridge_cv <- cv.glmnet(x=as.matrix(s.trset[,-c(1,2)]),y=Surv(s.trset$time, s.trset$status),type.measure = "deviance",nfold = 10,family = "cox",alpha = 0)
      best_ridge_coef <- coef(ridge_cv, s = ridge_cv$lambda.min)
      coxglmfit=cv.glmnet(x=as.matrix(s.trset[,-c(1,2)]),Surv(s.trset$time,s.trset$status),type.measure = "deviance",nfold = 10,family = "cox",alpha = 1,penalty.factor = 1 / abs(best_ridge_coef))
      coxglmpre=predict(coxglmfit,as.matrix(s.teset[,-c(1,2)]), type="response",s=lambda.best)
    }, times=100L)
  }
}
autoplot(mbm) 



