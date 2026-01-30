lower3 <- function(x,m=1) {
  #runs the Wei and Stram disagg method
  #must be a time series (ts) or xts or zoo object
  #x is the series
  #m order of disaggregation/aggregation
  xdis <- c(3,4,12)
  if(length(index(x))==0 & is.ts(x)==FALSE)stop("Must be a ts, zoo, or xts series")
  if(is.na(match(m,xdis)))stop("invalid disagg level")
  if(is.ts(x)==FALSE) {
    temp <- ts_ts(x)
    if(frequency(temp)!=1 & frequency(temp)!=4)stop("must be an annual or quarterly series")
    x <- temp
  }
  nu <- length(x)
  nw <- nu*m
  y <- matrix(0,nrow=27,ncol=14)
  y[,1] <- 999999
  n <- 0

  for(i in 1:1) {
    for(j in 0:1) {
    for(ka in 0:0) {

        n <- n + 1
        r <- 0
        z <- try(arima(x,order=c(i,j,ka)),silent=TRUE)

        if(!inherits(z,"try-error")) {
          y[n,1] <- z$aic
        y[n,c(2:4)] <- c(i,j,ka)
        if(i!=0) {
          y[n,5:(5+i-1)] <- z$coef[1:i]
        }
        if(ka!=0) {
          y[n,8:(8+ka-1)] <- z$coef[(i+1):(ka+i)]
        }
         if(i==2) {
           if(abs(y[n,5])> 2)y[n,1] <- 999999
           if(abs(y[n,6]) >1)y[n,1] <- 999999
           my.poly <- polynomial(c(-y[n,6:5],1))
           my.poly.r <- polyroot(my.poly)

           if(abs(Re(my.poly.r[1]))>0.999 | abs(Re(my.poly.r[2]))>2.0)y[n,1] <- 999999
           if(abs(Im(my.poly.r[1]))> 0.05 & abs(Im(my.poly.r[2]))>0.05)y[n,1] <- 999999
           xphig <- my.poly.r[1]^(1/m)

           if(abs(Im(xphig))>0.05)y[n,1] <- 999999
         }
        if(i==1 & m %% 2==0) {
          if(y[n,5]<0)y[n,1] <- 999999
        }
        y[n,13] <- z$sigma2
        y[n,14] <- r
        if(i==1 & abs(y[n,5])>=0.99)y[n,1] <- 999999
        if(ka==1 & abs(y[n,8])>=0.99)y[n,1] <- 999999
        }
}

}
  }
  for(i in 1:2) {
      n <- n + 1
  z <- try(arima(x,order=c(0,1,i)),silent=TRUE)
  if(!inherits(z,"try-error")) {
        y[n,1] <- z$aic
    y[n,c(2:4)] <- c(0,1,i)
      y[n,8:(8+i-1)] <- z$coef[1:i]


    y[n,13] <- z$sigma2
    y[n,14] <- 1
    if(i==1 & abs(y[n,8])>=0.99)y[n,1] <- 999999
      }
  }

  for(i in 1:1) {
    n <- n + 1
    z <- try(arima(x,order=c(0,0,i)),silent=TRUE)
    if(!inherits(z,"try-error")) {
      y[n,1] <- z$aic
      y[n,c(2:4)] <- c(0,0,i)
      y[n,8:(8+i-1)] <- z$coef[1:i]

      y[n,13] <- z$sigma2
      y[n,14] <- i
      if(i==1 & abs(y[n,8])>=0.99)y[n,1] <- 999999
    }
  }




  yy <- which(y[,1]==min(y[,1],na.rm=TRUE))[1]


  od1 <- y[yy,c(2:4,14)]
  zz <- y[yy,]
   np <- od1[1]
   nd <- od1[2]
   nr <- ifelse(np==0,sum(od1[2:3]),sum(od1[1:2]))
   ns <- sum(np+nr)
   phit <- NULL
   thet <- NULL
   tack <- NULL
   sig <- zz[13]
   if(np!=0) {
     npl <- 1:length(np)
     phit <- zz[5:(5+np-1)]
   }
   if(od1[3]!=0) {
     thet <--zz[8:(8+od1[3]-1)]
   }
   if(np!=0) {
    if(od1[3]!=0) tack <- tacvfARMA(phi=phit,theta = thet,maxLag=length(x)-1,sigma2=sig)
    if(od1[3]==0) tack <- tacvfARMA(phi=phit,maxLag=length(x)-1,sigma2=sig)
   }
    if(np==0)tack <- tacvfARMA(theta=thet,maxLag=length(x)-1,sigma2=sig)


     xp1 <- polynomial(rep(1,length=m))

     d1 <- y[yy,3]
     xp2 <- xp1^(2*(d1+1))

    v1 <- -(d1+1)*(m-1)
    v2 <- -v1

    ord2 <- max(y[yy,2:4])

#    r1 <- y[yy,14]
    r1 <- ord2
    ell1 <- (m*r1) + length(xp2)
    ell2 <- ell1 - v2

    v1a <- tack[1:(r1+1)]
    v1b <- ell1 - m
    n2 <- length(xp2)
    ad1 <- matrix(0,nrow=(r1+1),ncol=ell1)
    for(i in 1:(r1+1)) {
      j1 <- (i-1)*m + 1
      j2 <- j1 + (n2-1)
      ad1[i,j1:j2] <- coef(xp2)
    }


    ad2 <- matrix(0,nrow=(r1+1),ncol=ell2)
    ad2[,1:ell2] <- ad1[,(v2+1):ell1]
    n3 <- (n2+1)/2
    for(i in 2:n3) {
      ad2[1,i] <- 2* ad2[1,i]
    }
    for(i in 2:(r1+1)) {
      k <- v2
      for(j in 2:(v2+1)) {
      ad2[i,j] <- ad2[i,j] + ad1[i,k]
      k <- k-1
      }
    }


    cm5 <- matrix(0,nrow=(nu-d1),ncol=(nw-d1))
  nc <- (d1+1)*(m-1) + 1
  cpoly <- polynomial(rep(1,m))
  cp2 <- cpoly^(d1+1)
  cp3 <- coef(cp2)
  ncp1 <- length(cp3)
  cm5[1,1:ncp1] <- cp3
  for(i in 2:(nu-d1)) {
    j1 <- (i-1)*m + 1
    j2 <- j1 + (ncp1-1)
    cm5[i,j1:j2] <- cp3
  }
ng <- rev(dim(ad2))
nin1 <- max(ng)-min(ng)
  gma <- matrix(0,nrow=ng[1],ncol=ng[2])
  for(i in 1:ng[2]) {
    gma[i,i] <- 1
  }

  if(od1[1]!=0) {

    xphia <- polynomial(c(-y[yy,6:5],1))
    xphir <- polyroot(xphia)

    xphir <- Re(xphir)
    xphib <- sign(xphir)*abs(xphir)^(1/m)

  }
    if(od1[1]==2) {
      gp1 <- polynomial(c(-Re(xphib[1]),1))
      gp2 <- polynomial(c(-Re(xphib[2]),1))
      gp3 <- gp1 * gp2

      gm2 <- -coef(gp3)[1]^(1:nin1)
      gm1 <- -coef(gp3)[2]^(1:nin1)
      gma[(ng[2]+1):ng[1],(ng[2]-1)] <- gm2
      gma[(ng[2]+1):ng[1],ng[2]] <- gm1
    }
    if(od1[1]==1) {
      gm1 <- xphib[2]^(1:nin1)
      gma[(ng[2]+1):ng[1],ng[2]] <- gm1
    }

    xu.vec <- numeric(length=(nu-d1))
    xu.vec <- tack[1:(nu-d1)]
    xu.mat <- matrix(0,nrow=(nu-d1),ncol=(nu-d1))
    xu.mat[1,] <- xu.vec
    j <- (nu-d1)
    for(i in 2:(nu-d1)) {
      j <- j-1
      xu.mat[i,i:(nu-d1)] <- xu.vec[1:j]
    }
    xu.mat <- xu.mat + t(xu.mat) - diag(diag(xu.mat))
    xui.mat <- solve(xu.mat)
    ag1a <- ad2 %*% gma
    xv <- solve(ag1a,tack[1:dim(ad2)[1]])
    nv1 <- length(xv)

    nv2 <- nv1-1
    xv.vec <- numeric(length=(nw-d1))
    xv.vec[1:nv1] <- xv

    if(od1[1]==2) {
      xv.vec <- c(xv,xv[nv1]*coef(gp3)[2]^(1:(nw-nv1-d1)) + xv[(nv1-1)]*coef(gp3)[1]^(1:(nw-nv1-d1)))

    }


    if(od1[1]==1) {
      xv.vec <- c(xv,xv[nv1]*xphib[2]^(1:(nw-nv1-d1)))
    }

    xv.mat <- matrix(0,nrow=(nw-d1),ncol=(nw-d1))
    xv.mat[1,] <- xv.vec
    j <- (nw-d1)
    for(i in 2:(nw-d1)) {
      j <- j-1
      xv.mat[i,i:(nw-d1)] <- xv.vec[1:j]
    }
    xv.mat <- xv.mat + t(xv.mat) - diag(diag(xv.mat))

    dely <- matrix(0,nrow=(nu-d1),ncol=nu)
    if(d1==0) {
      dely <- diag(nu)
    } else {
      for(i in 1:(nu-d1)) {
        dely[i,i:(i+1)] <- c(-1,1)
      }
    }
    delx <- matrix(0,nrow=(nw-d1),ncol=nw)
    if(d1==0) {
      delx <- diag(nw)
    } else {
      for(i in 1:(nw-d1)) {
        delx[i,i:(i+1)] <- c(-1,1)
      }
    }
    if(d1!=0) {
    top1 <- delx
    bot1 <- c(rep(0,(nw-m)),rep(1,m))

    big1 <- rbind(top1,bot1)
    big1s <- solve(big1)
    }
    top2 <- xv.mat %*% t(cm5) %*% xui.mat %*% dely
    if(d1!=0) {

    bot2 <- c(rep(0,(nu-d1)),rep(1,d1))
    big2 <- rbind(top2,bot2)
    fin1 <- big1s %*% big2 %*% x
    } else {
      top2 <- xv.mat %*% t(cm5) %*% xui.mat
      fin1 <- top2 %*% x
}
     zzz <- list(bigy=c(od1[1:2],y[yy,14]),fin1=fin1)

    return(zzz)
  }




