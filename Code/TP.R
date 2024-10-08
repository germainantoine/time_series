library(readr)
library(zoo)
library(xts)
library(mgcv)
library(forecast)

cac <- read_delim("D:/Users/Antoine/Downloads/CAC40/350000.TXT", delim = "\t", escape_double = FALSE, col_types = cols(Date = col_date(format = "%d.%m.%Y")), trim_ws = TRUE)

head(cac)
str(cac)
summary(cac)

year<-format(cac$Date,"%Y")
mean.year<-tapply(cac$End, as.factor(year),mean)

plot(cac$Date,cac$Size, type='l',xlab="Temps",ylab="Encours du sous-jacent du CAC 40",col='royalblue3')

plot(cac$Date, cac$End, type = "l", xlab="Temps", ylab="Points du CAC 40", col='royalblue3')


plot(mean.year, type='b',axes=F,xlab="Temps",ylab="Points du CAC 40")
axis(1, c(1:length(mean.year)), labels=names(mean.year))
axis(2)
box()

cac$Date <- as.Date(cac$Date)
cac.xts<-xts(cac[,-1],order.by=cac$Date)
plot(cac.xts$End)

#####Estimation de la tendance
###Régression linéaire

linear <- lm(cac.xts$End ~ cac$Date)
plot(cac.xts$End, col='royalblue1')
lines(linear$fitted, col = "red")

plot(cac.xts$End-linear$fitted, col='royalblue1')

acf(cac.xts$End-linear$fitted,lag=7389)

hist(cac.xts$End-linear$fitted)

###moyenne mobile

fenetre <- 3
mobile <- stats::filter(cac.xts$End, filter = array(1/fenetre, dim = fenetre), method = c("convolution"), sides = 2, circular = FALSE)
mobile <-xts(mobile, order.by = cac$Date)
plot(cac.xts$End,col='royalblue1')
lines(mobile, col = "red")

plot(cac.xts$End-mobile, col='royalblue1')


acf(cac.xts$End,lag=7389)
#################noyau Gaussien

gaussien <- smooth.spline(cac.xts$End, spar = 0.5)
plot(cac.xts$End,col='royalblue1')
gaussien.xts<-xts(gaussien$y,order.by=cac$Date)
lines(gaussien.xts, col = "red")

plot(cac.xts$End-gaussien.xts, col='royalblue1')

#################polynomes locaux

poly <- loess(cac$End ~ index(cac), degree=1, span=0.9)
polyy<-xts(poly$fitted,order.by=cac$Date)
plot(cac.xts$End,type='l', col = 'royalblue1')
lines(polyy,col='red')

plot(cac.xts$End-polyy, col='royalblue1')

#################regression sur bases de splines

g<-gam(cac$End~s(index(cac), k=3))
ygam<-xts(g$fitted,order.by=cac$Date)
plot(cac.xts$End,type='l',col='royalblue1')
lines(ygam,col='red')

plot(cac.xts$End-ygam, col='royalblue1')

#####Meilleure estimation

plot(cac.xts$End,type='l', col = 'royalblue1')
lines(polyy,col='red')
lines(ygam,col='yellow')
lines(gaussien.xts, col = "green")
lines(mobile, col = "brown")
lines(linear$fitted, col = "violet")

moyerr <- (cac.xts$End-mobile)/7389
moyerr<-sum(moyerr[2,-2])

#####Lissage exponentiel
###Simple
expSmooth=function(x,alpha)
{
  xsmooth=x
  for(i in c(2:length(x)))
  {
    xsmooth[i]<-(1-alpha)*xsmooth[i-1]+alpha*x[i]
  }
  return(xsmooth)
}

alpha<-0.05
smooth<-expSmooth(cac$End,alpha)
plot(cac.xts$End, col='royalblue1')
smooth.xts<-xts(smooth,order.by=cac$Date)
lines(smooth.xts,col='red', lwd=2)

###Double
DoubleExpSmooth=function(x,alpha)
{
  xsmooth=x
  l<-array(x[1],dim=length(x))
  b<-array(x[2]-x[1],dim=length(x))
  
  for(i in c(2:length(x)))
  {
    l[i]<-xsmooth[i-1]+(1-(1-alpha)^2)*(x[i]-xsmooth[i-1])
    b[i]<-b[i-1]+alpha^2*(x[i]-xsmooth[i-1])
    xsmooth[i]<-l[i]+b[i]
  }
  
  res<-list()
  res$smooth<-xsmooth
  res$l<-l
  res$b<-b
  return(res)
}

alpha<-seq(10^-5,0.95,length=100)
forecast<-lapply(alpha,DoubleExpSmooth,x=cac$End)
erreur<-unlist(
  lapply(forecast,
         function(x){mean((tail(cac.xts$End,n-1)-head(x$smooth,n-1))^2)}))
plot(alpha,erreur,type='l')

smooth<-DoubleExpSmooth(cac$End,alpha[which.min(erreur)])
plot(cac.xts$End,type='l')
smooth_dates <- index(cac.xts)
smooth.xts <- xts(smooth$smooth, order.by = smooth_dates)
lines(smooth.xts,col='red')

#####Lissage exponentiel
###Simple
expSmooth=function(x,alpha)
{
  xsmooth=x
  for(i in c(2:length(x)))
  {
    xsmooth[i]<-(1-alpha)*xsmooth[i-1]+alpha*x[i]
  }
  return(xsmooth)
}

#alpha<-0.05
alpha<-0.8
smooth<-expSmooth(cac$End,alpha)
plot(cac.xts$End, col='royalblue1')
smooth.xts<-xts(smooth,order.by=cac$Date)
lines(smooth.xts,col='red', lwd=0.5)

###Double
cac$Date <- as.Date(cac$Date)
cac.xts<-xts(cac[,-1],order.by=cac$Date)

DoubleExpSmooth=function(x,alpha)
{
  xsmooth=x
  l<-array(x[1],dim=length(x))
  b<-array(x[2]-x[1],dim=length(x))
  
  for(i in c(2:length(x)))
  {
    l[i]<-xsmooth[i-1]+(1-(1-alpha)^2)*(x[i]-xsmooth[i-1])
    b[i]<-b[i-1]+alpha^2*(x[i]-xsmooth[i-1])
    xsmooth[i]<-l[i]+b[i]
  }
  
  res<-list()
  res$smooth<-xsmooth
  res$l<-l
  res$b<-b
  return(res)
}
n=length(cac$End)
alpha<-seq(10^-5,0.95,length=100)
forecast<-lapply(alpha,DoubleExpSmooth,x=cac$End)
erreur<-unlist(lapply(forecast,function(x){mean((tail(cac$End,n-1)-head(x$smooth,n-1))^2)}))
plot(alpha,erreur,type='l')

smooth<-DoubleExpSmooth(cac$End,alpha[which.min(erreur)])
plot(cac.xts$End,type='l')
smooth.xts<-xts(smooth$smooth,order.by=cac$Date)
lines(smooth.xts,col='red') 


#####Prévision 
predict.expSmooth<-function(Xsmooth,inst,horizon,smooth.type="double")
{
  
  if(smooth.type=="simple")
  {
    n<-length(Xsmooth)
    prev<-c(Xsmooth[1:inst],rep(Xsmooth[inst],horizon))
  }
  
  if(smooth.type=="double")
  {
    n<-length(Xsmooth$smooth)
    prev<-c(Xsmooth$smooth[1:inst],Xsmooth$l[inst]+Xsmooth$b[inst]*c(1:horizon))
  }
  return(prev)
}


cac_0<-cac$End[6589:7389]
cac_1<-cac$End[7343:7389]
cac_pr<-cac$End[6589:7389]

cac$Date <- as.Date(cac$Date)
cac.xts<-xts(cac[,-1],order.by=cac$Date)


alpha=0.2
X.d.exp.mooHW<-DoubleExpSmooth(cac_pr,alpha)
prevHW<-predict.expSmooth(X.d.exp.mooHW,
                          inst=754,horizon=46, smooth.type="double")
plot(cac_pr,pch=20,ylim=range(cac_pr,prevHW),type='l')
lines(prevHW,col='red',lwd=2,type='l')
abline(v=754,lty='dashed')

X.d.exp.mooHW<-DoubleExpSmooth(cac$End,alpha)
prevHW<-predict.expSmooth(X.d.exp.mooHW,
                          inst=754,horizon=46, smooth.type="double")
cac$Date <- as.Date(cac$Date)
cac.xts<-xts(cac[,-1],order.by=cac$Date)
plot(cac.xts$End,pch=20,xlab="Temps",ylab="Investissements du CAC40",type='l')
lines(prevHW,col='red',lwd=2,type='l')
abline(v=754,lty='dashed')


##Prévision par lissage simple 
cac_0.xts<-cac.xts$End[6589:7342]

cac_1.xts<-cac.xts$End[7343:7389]
ets1 <- ets(cac_0.xts, model="ANN")
ets2 <- ets(cac_0.xts, model="AAN")

ets1.forecast <- forecast(ets1, h=length(cac_1.xts))
ets2.forecast <- forecast(ets2, h=length(cac_1.xts))

plot(cac$Date[6589:7389],c(cac_0.xts,cac_1.xts),xlab="Temps", ylab="Investissements CAC40",type='l')
lines(cac$Date[6589:7342],cac_0.xts,type='l')
lines(cac$Date[7343:7389],cac_1.xts,col='red')
lines(cac$Date[7343:7389],  ets1.forecast$mean, col='blue')
lines(cac$Date[7343:7389],  ets2.forecast$mean, col='purple')

## Holt-Winters
cac$Date <- as.Date(cac$Date)
cac.xts<-xts(cac[,-1],order.by=cac$Date)

cac_0.xts<-cac.xts$End[6589:7342]
cac_1.xts<-cac.xts$End[7343:7389]

hw<-HoltWinters(cac.xts$End, gamma=F)
hw.forecast<-predict(hw, n.ahead=length(cac_1.xts)) 
plot(cac$Date[6589:7389],c(cac_0.xts,cac_1.xts),type='l', 
     xlab="Temps", ylab="Investissement CAC40",)
lines(cac$Date[7343:7389],cac_1.xts,col='red')
lines(cac$Date[7343:7389],  hw.forecast, col='blue')
plot(cac$Date, hw)
lines(cac$Date, cac.xts$End)



