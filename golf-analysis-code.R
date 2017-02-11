# Read and clean data

# PGA2012R<-read.delim("/Users/William/BYU-Research/Dr.Fellingham-golf/PGATour/Exports/shot-2012riviera.TXT",sep=";",header=FALSE,as.is=TRUE)
# 
# names(PGA2012R)<-column.names
# 
# riv2012<-PGA2012R[,c(5,10,13:15,17,19:29,37)]
# names(riv2012)<-c("Player..","Round","Hole","Hole.Score","Par.Value","Shot","X..of.Strokes","From.Location.Scorer","From.Location.Laser.","To.Location.Scorer.","To.Location.Laser.","Distance","Distance.to.Pin","In.the.Hole.Flag","Around.the.Green.Flag","X1st.Putt.Flag","Distance.to.Hole.after.the.Shot","Distance.from.Center")
# riv2012.sub<-riv2012[riv2012$X..of.Strokes==1,]
# rep.times<-riv2012.sub$Hole.Score[riv2012.sub$Shot==1]
# riv2012.sub$obs.id<-rep(1:length(rep.times),times=rep.times)
# riv2012.sub$Shot.Real<-unlist(aggregate(X..of.Strokes~obs.id,data = riv2012.sub,FUN=cumsum)[,2])
# 
# riv2012.sub$ParMinusShot<-riv2012.sub$Par.Value-riv2012.sub$Shot.Real
# 
# holes<-riv2012.sub$Hole[riv2012.sub$Shot.Real==1]
# score.to.par<-riv2012.sub$Hole.Score[riv2012.sub$Shot.Real==1]-riv2012.sub$Par.Value[riv2012.sub$Shot.Real==1]
# 
# drive.dist<-riv2012.sub$Distance[riv2012.sub$Shot.Real==1]
# 
# drive.off<-riv2012.sub$Distance.from.Center[riv2012.sub$Shot.Real==2]
# hole.in.one.id<-as.numeric(riv2012.sub$obs.id[which(riv2012.sub$Hole.Score==1)])
# drive.off2<-drive.off[1:(hole.in.one.id-1)]
# drive.off3<-drive.off[hole.in.one.id:length(drive.off)]
# drive.off.full<-c(drive.off2,0,drive.off3)
# 
# sub.scramble<-riv2012.sub[riv2012.sub$ParMinusShot<=1  | (riv2012.sub$ParMinusShot==2 & riv2012.sub$In.the.Hole.Flag=="Y"),c("Distance.to.Pin","From.Location.Scorer.","obs.id")]
# sub.scramble$Distance.to.Pin[which(sub.scramble$From.Location.Scorer.=="Green")]<-0
# dist.scramble<-unlist(aggregate(Distance.to.Pin~obs.id,data = sub.scramble,FUN = sum)[,2])
# 
# first.putt<-riv2012.sub[riv2012.sub$X1st.Putt.Flag=="Y" | (riv2012.sub$X1st.Putt.Flag!="Y" & riv2012.sub$In.the.Hole.Flag=="Y" & riv2012.sub$From.Location.Scorer.!="Green"),c("Distance.to.Pin","obs.id")]
# first.putt2<-unlist(aggregate(Distance.to.Pin~obs.id,data = first.putt,FUN = sum)[,2])
# 
# putt.made<-riv2012.sub[riv2012.sub$In.the.Hole.Flag=="Y",c("Distance.to.Pin","From.Location.Scorer.")]
# putt.made$Distance.to.Pin[putt.made$From.Location.Scorer.!="Green"]<-0
# putt.made.dist<-putt.made$Distance.to.Pin
# 
# clean.dat2012<-data.frame(holes,score.to.par,(drive.dist/36),(drive.off.full/36),(dist.scramble/36),(first.putt2/12),(putt.made.dist/12))
# 
# names(clean.dat2012)<-c("hole","score.to.par","drive.dist","drive.off.full","dist.scramble","first.putt2","putt.made.dist")
# save(clean.dat2012,file = "~/STAT GRAD/Stat637Golf2012")
# load(file = "~/STAT GRAD/Stat637Golf2012")
# 
# PGA2011R<-read.delim("/Users/William/BYU-Research/Dr.Fellingham-golf/PGATour/Exports/shot-2011riviera.TXT",sep=";",header=FALSE,as.is=TRUE)
# 
# names(PGA2011R)<-column.names
# 
# PGA2010R<-read.delim("/Users/William/BYU-Research/Dr.Fellingham-golf/PGATour/Exports/shot-2010riviera.TXT",sep=";",header=FALSE,as.is=TRUE)
# 
# names(PGA2010R)<-column.names
# 
# PGA2009R<-read.delim("/Users/William/BYU-Research/Dr.Fellingham-golf/PGATour/Exports/shot-2009riviera.TXT",sep=";",header=FALSE,as.is=TRUE)
# 
# names(PGA2009R)<-column.names
# 
# PGA2008R<-read.delim("/Users/William/BYU-Research/Dr.Fellingham-golf/PGATour/Exports/shot-2008riviera.TXT",sep=";",header=FALSE,as.is=TRUE)
# 
# names(PGA2008R)<-column.names
# 
# PGA2007R<-read.delim("/Users/William/BYU-Research/Dr.Fellingham-golf/PGATour/Exports/shot-2007riviera.TXT",sep=";",header=FALSE,as.is=TRUE)
# 
# names(PGA2007R)<-column.names
# 
# PGA2006R<-read.delim("/Users/William/BYU-Research/Dr.Fellingham-golf/PGATour/Exports/shot-2006riviera.TXT",sep=";",header=FALSE,as.is=TRUE)
# 
# names(PGA2006R)<-column.names
# 
# PGA2005R<-read.delim("/Users/William/BYU-Research/Dr.Fellingham-golf/PGATour/Exports/shot-2005riviera.TXT",sep=";",header=FALSE,as.is=TRUE)
# 
# names(PGA2005R)<-column.names
# 
# PGA2004R<-read.delim("/Users/William/BYU-Research/Dr.Fellingham-golf/PGATour/Exports/shot-2004riviera.TXT",sep=";",header=FALSE,as.is=TRUE)
# 
# names(PGA2004R)<-column.names[-39]
# 
# PGA2003R<-read.delim("/Users/William/BYU-Research/Dr.Fellingham-golf/PGATour/Exports/shot-2003riviera.TXT",sep=";",header=FALSE,as.is=TRUE)
# 
# names(PGA2003R)<-column.names[-39]
# 
# PGAR<-rbind(PGA2003R,PGA2004R,PGA2005R[,-39],PGA2006R[,-39],PGA2007R[,-39],PGA2008R[,-39],PGA2009R[,-39],PGA2010R[,-39],PGA2011R[,-39],PGA2012R[,-39])
# 
# ####################################################################
# ############                 Model              ####################
# ####################################################################
# 
# library(msm)
# library(MCMCpack)
# 
# clean.dat<-clean.dat2012
# names(clean.dat)<-c("hole","score.to.par","drive.dist","drive.off.full","dist.scramble","first.putt2","putt.made.dist")
# clean.dat$drive.dist<-clean.dat$drive.dist/25
# clean.dat$long.drive<-clean.dat$drive.dist>300
# clean.dat$drive.off.full<-clean.dat$drive.off.full/5
# clean.dat$dist.scramble<-clean.dat$dist.scramble/25
# #clean.dat$dist.scramble[clean.dat$dist.scramble==0]<-NA
# clean.dat$first.putt2<-clean.dat$first.putt2/5
# clean.dat$putt.made.dist<-clean.dat$putt.made.dist/5

#prepare for mcmc

length<-9
burn<-1
#n<-nrow(clean.dat)
h<-18
obs<-nrow(clean.dat)/h
clean.dat$score.to.par[clean.dat$score.to.par>2]<-2
clean.dat$score.to.par[clean.dat$score.to.par<(-1)]<-(-1)
y<-matrix(nrow = obs,ncol = h)
for(k in 1:h){  
  y[,k]<-clean.dat$score.to.par[clean.dat$hole==k]
}

#exploratory scatter plots
par(mfrow=c(1,2))
plot(clean.dat$drive.dist,clean.dat$drive.off.full,col=(clean.dat$score.to.par+2),ylab="Off Center",xlab="Drive Dist")
legend("topleft",legend=c("Birdie","Par","Bogey","Double"),col=c(1:4),pch=1,cex=0.8)
plot(clean.dat$dist.scramble,clean.dat$first.putt2,col=(clean.dat$score.to.par+2),ylab="First Putt",xlab="Scramble Dist")
legend("topright",legend=c("Birdie","Par","Bogey","Double"),col=c(1:4),pch=1,cex=0.8)

X<-array(0,dim = c(obs,4,h))
for (k in 1:h){
  X[,,k]<-as.matrix(clean.dat[clean.dat$hole==k,c("long.drive","drive.off.full","dist.scramble","first.putt2")])
}

beta<-matrix(0,ncol=ncol(X),nrow = h)
beta.save<-array(dim = c(h,ncol(X),(burn+length)))
beta.save[,,1]<-beta

mu<-matrix(0,ncol=5,nrow=(burn+length))
mu[1,]<-rep(1,5)

m<-rep(0,5)
V<-10*diag(5)

w<-h+1
I<-diag(5)
Sigma<-array(dim = c(ncol(X),ncol(X),(burn+length)))
Sigma[,,1]<-I

z.o<-matrix(0,nrow=obs,ncol = h)
Z<-array(dim = c(obs,h,(burn+length)))
Z[,,1]<-z.o

gamma<-matrix(rep(c(0,1,2),each=h),ncol = 3)
gamma.save<-array(dim = c(h,ncol(gamma),(burn+length)))
gamma.save[,,1]<-gamma

birdie<-list()
even<-list()
bogey<-list()
dubbogey<-list()

#run mcmc
for(k in 1:h){
  birdie[[k]]<-which(y[,k]==-1)
  even[[k]]<-which(y[,k]==0)
  bogey[[k]]<-which(y[,k]==1)
  dubbogey[[k]]<-which(y[,k]==2)
}

for(d in 2:(length+burn)){
  # update beta (k = 1,...,n_h)
  for (k in 1:h){
    sigstar <- solve( solve(Sigma[,,(d-1)]) + t(X[,,k])%*%X[,,k] )
    mustar <- sigstar %*% (solve(Sigma[,,(d-1)])%*%mu[(d-1),] + t(X[,,k])%*%Z[,k,(d-1)])
    beta.save[k,,d] <- mvrnorm(1,mustar,sigstar)
  }
  
  # update mu
  Vstar <- solve( solve(V) +  solve(Sigma[,,(d-1)]) )
  mstar <- Vstar %*% (solve(V)%*%m + solve(Sigma[,,(d-1)])%*%colMeans(beta.save[,,d]))
  mu[d,] <- mvrnorm(1,mstar,Vstar)
  
  # update Sigma
  wstar<-w+h
  S.mu<-matrix(0,nrow = 5,ncol = 5)
  for(k in 1:h){
    S.mu<-S.mu+(beta.save[k,,d]-mu[d,])%*%t((beta.save[k,,d]-mu[d,]))
  }
  
  Istar<-solve(I+S.mu)
  Sigma[,,d] <- riwish(wstar,Istar)

  # update Z
  for(k in 1:h){ 
    for (i in 1:obs){
      Z[i,k,d]<-rtnorm(1,t(X[i,,k])%*%beta.save[k,,d],1,gamma.save[k,3,(d-1)],Inf)
      if(y[i,k]==1)
      {
        Z[i,k,d]<-rtnorm(1,t(X[i,,k])%*%beta.save[k,,d],1,gamma.save[k,2,(d-1)],gamma.save[k,3,(d-1)])
      }
      if(y[i,k]==0)
      {
        Z[i,k,d]<-rtnorm(1,t(X[i,,k])%*%beta.save[k,,d],1,gamma.save[k,1,(d-1)],gamma.save[k,2,(d-1)])
      }
      if(y[i,k]==-1)
      {
        Z[i,k,d]<-rtnorm(1,t(X[i,,k])%*%beta.save[k,,d],1,-Inf,gamma.save[k,1,(d-1)])
      }
    }
  }
  # update gamma
  for(k in 1:h){
    lower1<-max(Z[birdie[[k]],k,d])
    upper1<-min(min(Z[even[[k]],k,d]),gamma.save[k,2,(d-1)])
    gamma.save[k,1,d]<-runif(1,min(lower1,upper1),max(lower1,upper1))
  
    lower2<-max(Z[even[[k]],k,d],gamma.save[k,1,(d-1)])
    upper2<-min(Z[bogey[[k]],k,d],gamma.save[k,3,(d-1)])
    gamma.save[k,2,d]<-runif(1,min(lower2,upper2),max(lower2,upper2))
  
    lower3<-max(Z[bogey[[k]],k,d],gamma.save[k,2,(d-1)])
    upper3<-min(Z[dubbogey[[k]],k,d])
    gamma.save[k,3,d]<-runif(1,min(lower3,upper3),max(lower3,upper3))
  }
}

##### Frequentist models  #####
# library(VGAM)
# clean.dat2012.par4<-clean.dat2012[clean.dat2012$hole%in%c(2,3,5,7:10,12,13,15,18),]
# clean.dat2012.par4$drive.discrete<-NA
# clean.dat2012.par4$drive.discrete[clean.dat2012.par4$drive.dist<=250]<-"250-"
# clean.dat2012.par4$drive.discrete[clean.dat2012.par4$drive.dist>250 & clean.dat2012.par4$drive.dist<=275]<-"250-275"
# clean.dat2012.par4$drive.discrete[clean.dat2012.par4$drive.dist>275 & clean.dat2012.par4$drive.dist<=300]<-"275-300"
# clean.dat2012.par4$drive.discrete[clean.dat2012.par4$drive.dist>300 & clean.dat2012.par4$drive.dist<=325]<-"300-325"
# clean.dat2012.par4$drive.discrete[clean.dat2012.par4$drive.dist>325]<-"325+"
# 
# fit1<-vglm(score.to.par~drive.discrete,data = clean.dat2012.par4,family = cumulative(parallel = T))
# summary(fit1)
# 
# clean.dat2012.par4$drive.simple<-clean.dat2012.par4$drive.dist/25
# clean.dat2012.par4$first.putt5<-clean.dat2012.par4$first.putt2/5
# clean.dat2012.par4$putt.made.dist5<-clean.dat2012.par4$putt.made.dist/5
# clean.dat2012.par4$dist.scramble[clean.dat2012.par4$dist.scramble==0]<-NA
# clean.dat2012.par4$dist.scramble25<-clean.dat2012.par4$dist.scramble/25
# fit2<-vglm(score.to.par~drive.simple,data = clean.dat2012.par4,family = cumulative(parallel = T))
# summary(fit2)
# 
# fit2.1<-vglm(score.to.par~putt.made.dist10,data = clean.dat2012.par4,family = cumulative(parallel = T))
# summary(fit2.1)
# 
# fit2.2<-vglm(score.to.par~dist.scramble25,data = clean.dat2012.par4,family = cumulative(parallel = T))
# summary(fit2.2)
# 
# fit2.3<-vglm(score.to.par~putt.made.dist5+dist.scramble25+first.putt5,data = clean.dat2012.par4,family = cumulative(parallel = T))
# summary(fit2.3)
# 
# clean.dat2012.par4$long.drive<-clean.dat2012.par4$drive.dist>300
# fit2.4<-vglm(score.to.par~long.drive,data = clean.dat2012.par4,family = cumulative(parallel = T))
# summary(fit2.4)
# 
# clean.dat2012.par5<-clean.dat2012[clean.dat2012$hole%in%c(1,11,17),]
# clean.dat2012.par5$drive.simple<-clean.dat2012.par5$drive.dist/25
# fit3<-vglm(score.to.par~drive.simple,data = clean.dat2012.par5,family = cumulative(parallel = T))
# summary(fit3)
# 
# clean.dat2012.par5$long.drive<-clean.dat2012.par5$drive.dist>300
# fit4<-vglm(score.to.par~long.drive,data = clean.dat2012.par5,family = cumulative(parallel = T))
# summary(fit4)
# 

# Par 3s: 4, 6, 14, 16 
# Par 4s: 2, 3, 5, 7, 8, 9, 10, 12, 13, 15, 18
# Par 5s: 1, 11, 17

############################################################
#####     Results from long run         ##########
############################################################

# First run (variables = long drive, drive off center, scrambling, 1st putt length, putt length made)
load("~/Downloads/betas1.Rdata")
load("~/Downloads/Z1.Rdata")
length<-5000
burn<-1000
M<-length+burn
h<-18

par(mfrow=c(3,5))
for(t in 1:3){
  for(s in 1:5)
    plot(beta.save[t,s,-(1:burn)],type = "l")
}

for(t in 1:3){
  for(s in 1:5)
    plot(density(beta.save[t,s,-(1:burn)]),main = "")
}

lesszero<-matrix(nrow = h,ncol = 5)

for(t in 1:18){
  for(s in 1:5)
    lesszero[t,s]<-length(which(beta.save[t,s,-(1:burn)]<0))/length
}

beta.means<-matrix(nrow = h,ncol = 5)

for(t in 1:18){
  for(s in 1:5)
    beta.means[t,s]<-mean(beta.save[t,s,-(1:burn)])
}

# Second run (variables = drive length, drive off center, scrambling, 1st putt length, putt length made)
load("~/Downloads/betas2.Rdata")
length<-5000
burn<-1000
M<-length+burn
h<-18

par(mfrow=c(3,5))
for(t in 1:3){
  for(s in 1:5)
    plot(beta.save[t,s,-(1:burn)],type = "l")
}

for(t in 1:3){
  for(s in 1:5)
    plot(density(beta.save[t,s,-(1:burn)]),main="")
}

lesszero<-matrix(nrow = h,ncol = 5)

for(t in 1:18){
  for(s in 1:5)
    lesszero[t,s]<-length(which(beta.save[t,s,-(1:burn)]<0))/length
}

beta.means<-matrix(nrow = h,ncol = 5)

for(t in 1:18){
  for(s in 1:5)
    beta.means[t,s]<-mean(beta.save[t,s,-(1:burn)])
}

# Second run (variables = long drive, drive off center, scrambling, 1st putt length)
load("~/Downloads/betas3.Rdata")
length<-5000
burn<-1000
M<-length+burn
h<-18

par(mfrow=c(3,4))
for(t in 1:3){
  for(s in 1:4)
    plot(beta.save[t,s,-(1:burn)],type = "l")
}

for(t in 1:3){
  for(s in 1:4)
    plot(density(beta.save[t,s,-(1:burn)]),main="")
}

lesszero<-matrix(nrow = h,ncol = 4)

for(t in 1:18){
  for(s in 1:4)
    lesszero[t,s]<-length(which(beta.save[t,s,-(1:burn)]<0))/length
}


beta.means<-matrix(nrow = h,ncol = 4)

for(t in 1:18){
  for(s in 1:4)
    beta.means[t,s]<-mean(beta.save[t,s,-(1:burn)])
}

load("C:/Users/Stephen/Downloads/betas6.RData")
holes.of.interest <- c()
#5,4,4,3,4,3,4,4,4,4,5,4,4,3,4,3,5,4
colors = c(2,4,4,3,4,3,4,4,4,4,2,4,4,3,4,3,2,4)
#plot beta densities
par(mfrow=c(2,2))
plot(density(beta.save[1,1,-c(1:1000)]),main="Driving Distance Over 300 Yards",lwd=3,col=2,ylim=c(0,3.5))
for(i in 2:18) {
  lines(density(beta.save[i,1,-c(1:1000)]),col=colors[i],lwd=3)
}
legend("topleft",legend=c("Par 3","Par 4","Par 5"),lty=1,col=c(3,4,2),lwd=3,cex=0.5)

plot(density(beta.save[1,2,-c(1:1000)]),main="Distance Off From Center",lwd=3,col=2,xlim=c(-0.05,1),ylim=c(0,18))
for(i in 2:18) {
  lines(density(beta.save[i,2,-c(1:1000)]),col=colors[i],lwd=3)
}
legend("topright",legend=c("Par 3","Par 4","Par 5"),lty=1,col=c(3,4,2),lwd=3,cex=0.5)

plot(density(beta.save[2,3,-c(1:1000)]),main="Distance of Scramble",lwd=3,col=4,ylim=c(0,14),xlim=c(-0.05,1))
for(i in c(3:10,12:16,18)) {
  lines(density(beta.save[i,3,-c(1:1000)]),col=colors[i],lwd=3)
}
for(i in c(1,11,17)) {
  lines(density(beta.save[i,3,-c(1:1000)]),col=colors[i],lwd=5)
}
legend("topright",legend=c("Par 3","Par 4","Par 5"),lty=1,col=c(3,4,2),lwd=3,cex=0.5)

plot(density(beta.save[1,4,-c(1:1000)]),main="Distance of First Putt",lwd=3,col=2,xlim=c(-0.1,0.15),ylim=c(0,55))
for(i in 2:18) {
  lines(density(beta.save[i,4,-c(1:1000)]),col=colors[i],lwd=3)
}
legend("topright",legend=c("Par 3","Par 4","Par 5"),lty=1,col=c(3,4,2),lwd=3,cex=0.5)

load("C:/Users/Stephen/Downloads/Z6.RData")

Z.means<-matrix(0,ncol = 18,nrow = 437)
for(j in 1:437){
  Z.means[j,]<-rowMeans(Z[j,,])
}

par(mfrow=c(1,1))

plot(density(Z.means[,1]),lwd=4,col=colors[1],ylim=c(0,2),main="Scoring Distributions for Three Holes")
lines(density(Z.means[,4]),lwd=4,col=colors[4])
lines(density(Z.means[,2]),lwd=4,col=colors[2])
legend("topright",legend=c("Hole 1, Par 5","Hole 4, Par 3","Hole 2, Par 4"),col=c(2,3,4),lwd=4)

plot(density(Z.means[,2]),lwd=4)
#plot(density(Z.means[,3]),lwd=4)
plot(density(Z.means[,4]),lwd=4)
plot(density(Z.means[,5]),lwd=4)
#plot(density(Z.means[,6]),lwd=4)
#plot(density(Z.means[,7]),lwd=4)
plot(density(Z.means[,8]),lwd=4)
plot(density(Z.means[,9]),lwd=4)
#plot(density(Z.means[,10]),lwd=4)
plot(density(Z.means[,11]),lwd=4)
plot(density(Z.means[,12]),lwd=4)
#plot(density(Z.means[,13]),lwd=4)
plot(density(Z.means[,14]),lwd=4)
plot(density(Z.means[,15]),lwd=4)
#plot(density(Z.means[,16]),lwd=4)
#plot(density(Z.means[,17]),lwd=4)
plot(density(Z.means[,18]),lwd=4)

#posterior predictive
M <- 6000
n.h <- 2
pred.dist <- numeric(M)
for(i in 1:M) {
  index <- sample(1:M,1)
  beta.samp <- beta.save[n.h,,index]
  Z.samp <- numeric(437)
  for(j in 1:437) {
    X.samp <- numeric(4)
    X.samp[1] <- sample(X[,1,n.h],1)
    X.samp[2] <- sample(X[,2,n.h],1)
    X.samp[3] <- sample(X[,3,n.h],1)
    X.samp[4] <- sample(X[,4,n.h],1)
    Z.samp[j] <- rnorm(1,X.samp%*%beta.samp,1)
  }
  score <- numeric(4)
  score[1] <- sum(which(Z.samp < 0))
  score[2] <- sum(which(Z.samp > 0 && Z.samp < 1))
  score[3] <- sum(which(Z.samp > 1 && Z.samp < 2))
  score[4] <- sum(which(Z.samp > 2))
  pred.dist[i] <- which.max(score)-2
}


M<-length+burn

#prob.pred<-array(0,dim=c(n.test,4,h))
z.pred<-matrix(0,nrow = length,ncol = h)
#cutoffs<-c(0,1,2)

X.samp<-array(dim = c(18,4,5000))
for(j in 1:4){
  for(k in 1:h){
    X.samp[k,j,]<-sample(X[,j,k],5000,replace = TRUE)
  }
}

for(k in 1:h){
  index<-sample((burn+1):M,length,replace = T)
  for(d in 1:length){
    z.pred[d,k]<-rnorm(1,t(X.samp[k,,d])%*%beta.save[k,,index[d]],1)
  }
}

prob.pred<-matrix(0,nrow = 18,ncol = 4)

for(k in 1:18){
  prob.pred[k,]<-hist(z.pred[,k],plot = F,breaks=c(-Inf,0,1,2,Inf))$count/5000
}

pred.avg.score<-numeric(18)
values<-matrix(c(-1,0,1,2),ncol = 1)
pred.avg.score<-(prob.pred%*%values)


ranks<-data.frame(order(pred.avg.score,decreasing = TRUE),sort(pred.avg.score,decreasing = TRUE))
names(ranks)<-c("Hole","Predicted Mean Score to Par")


obs.avg.score<-numeric(18)
for(k in 1:18){
  obs.avg.score[k]<-mean(clean.dat$score.to.par[clean.dat$hole==k])
}

obs.ranks<-data.frame(order(obs.avg.score,decreasing = TRUE),sort(obs.avg.score,decreasing = TRUE))
names(obs.ranks)<-c("Hole","Observed Mean Score to Par")


comp.ranks<-cbind(ranks,obs.ranks)
names(comp.ranks)<-c("Predicted Hole","Predicted Mean Score to Par","Observed Hole","Observed Mean Score to Par")