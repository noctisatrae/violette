library("lars")

convcheck<-function(beta1,beta2){
		a<-apply(abs(beta1+beta2),2,max)
		b<-apply(abs(beta1-beta2),2,max)		
	    d<-length(a)
	    x<-rep(1,d)
	    for (i in 1:d){
	       x[i]<-min(a[i],b[i])
	    }
	    max(x)
}
backsolvet <-
function(r, x, k=ncol(r))
{
  backsolve(r,x,k,transpose=TRUE)
}

delcol <-
function(r, z, k = p)
{
	p <- dim(r)[1]
	r <- r[,  - k, drop = FALSE]
	z <- as.matrix(z)
	dz <- dim(z)
	storage.mode(r) <- storage.mode(z) <- "double"
	.Fortran("delcol",
		r,
		as.integer(p),
		as.integer(k),
		z,
		as.integer(dz[1]),
		as.integer(dz[2]),
                PACKAGE="lars")[c(1, 4)]
}

downdateR <-
function(R, k = p)
{
	p <- dim(R)[1]
	if(p == 1)
		return(NULL)
	R <- delcol(R, rep(1, p), k)[[1]][ - p,  , drop = FALSE]
	attr(R, "rank") <- p - 1
	R	# Built-in Splus utility
}


updateRR<-function(xnew, R = NULL, xold, lambda,eps = .Machine$double.eps)
{
	xtx <- (sum(xnew^2)+lambda)/(1+lambda)
	norm.xnew <- sqrt(xtx)
	if(is.null(R)) {
		R <- matrix(norm.xnew, 1, 1)
		attr(R, "rank") <- 1
		return(R)
	}
	Xtx <- drop(t(xnew) %*% xold)/(1+lambda)
	r <- backsolvet(R, Xtx)
	rpp <- norm.xnew^2 - sum(r^2)
	rank <- attr(R, "rank")	
	if(rpp <= eps)
		rpp <- eps
	else {
		rpp <- sqrt(rpp)
		rank <- rank + 1
	}
	R <- cbind(rbind(R, 0), c(r, rpp))
	attr(R, "rank") <- rank
	R
}

solvebeta<-function(x, y, paras, max.steps, sparse=c("penalty","varnum"), eps = .Machine$double.eps)
{
       sparse <- match.arg(sparse)
       if (missing(sparse)) sparse <- "penalty"
	nm <- dim(x)
	n <- nm[1]
	m <- nm[2]
	im <- seq(m)
	one <- rep(1, n)
	vn <- dimnames(x)[[2]]
        lambda<-paras[1]
        if(lambda>0){
	   maxvars <- m
        }
        if (lambda==0) {
           maxvars <- min(m,n-1)
           if (m==n){
             maxvars<-m
           }
        }
	d1 <- sqrt(lambda)
	d2 <- 1/sqrt(1 + lambda)
	Cvec <- drop(t(y) %*% x) * d2
        ssy <- sum(y^2)
	residuals <- c(y, rep(0, m))
        if(missing(max.steps)) {max.steps <- 50 * maxvars}
        penalty<-max(abs(Cvec))
        if (sparse=="penalty" && penalty*2/d2<=paras[2])
          { 
           beta<-rep(0,m)
          }
        else {
               beta <- rep(0, m)
               first.in <- integer(m)
               active <- NULL
               ignores <- NULL
               actions <- as.list(seq(max.steps))
               drops <- FALSE
               Sign <- NULL
               R <- NULL
	       k <- 0
	       while((k < max.steps) & (length(active) < maxvars - length(ignores))) 
		{
		 action <- NULL
		 k <- k + 1
        	 inactive <- if(k == 1) im else im[ - c(active, ignores)]
		 C <- Cvec[inactive]
		 Cmax <- max(abs(C))
		 if(!any(drops)) {
			new <- abs(C) == Cmax 
			C <- C[!new]
			new <- inactive[new]       
			for(inew in new) {
                        R <- updateRR(x[, inew], R, x[, active], lambda) 
                               if(attr(R, "rank") == length(active)) {
					nR <- seq(length(active))
					R <- R[nR, nR, drop = FALSE]
					attr(R, "rank") <- length(active)
					ignores <- c(ignores, inew)
					action <- c(action,  - inew)
				}
				else {
					if(first.in[inew] == 0)
						first.in[inew] <- k
					active <- c(active, inew)
					Sign <- c(Sign, sign(Cvec[inew]))
					action <- c(action, inew)
				}
			}
		}
		else action <-  - dropid
                Gi1 <- backsolve(R, backsolvet(R, Sign))
            	A <- 1/sqrt(sum(Gi1 * Sign))
            	w <- A * Gi1
                u1<-drop(x[,active,drop=FALSE]%*%w*d2)  
                u2<-rep(0,m)
                u2[active]<-d1*d2*w
                u<-c(u1,u2)
                if(lambda>0){
	            maxvars <- m-length(ignores)
                           }
                if (lambda==0){
                    maxvars <- min(m-length(ignores),n-1)
                           }
       		if(length(active) == maxvars - length(ignores)) {
			gamhat <- Cmax/A
		}
		else {
			a <- (drop(u1 %*% x[,  - c(active, ignores)]) + d1 *
				u2[ - c(active, ignores)]) * d2
      			gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
				gamhat <- min(gam[gam > eps], Cmax/A)                        
			Cdrop <- c(C - gamhat * a,  - C + gamhat * a) - (Cmax -
				gamhat * A)
		}
		dropid <- NULL
		b1 <- beta[active]
		z1 <-  - b1/w
		zmin <- min(z1[z1 > eps], gamhat)
		if(zmin < gamhat) {
			gamhat <- zmin
			drops <- z1 == zmin
		}
               else drops <- FALSE
               beta2<-beta          
               beta[active] <- beta[active] + gamhat * w
               residuals <- residuals - (gamhat*u)
               Cvec <- (drop(t(residuals[1:n]) %*% x) + d1 * residuals[ - (1:n)]) * d2
               penalty <- c(penalty,penalty[k]-abs(gamhat*A))
               if(sparse=="penalty" && rev(penalty)[1]*2/d2<=paras[2]){
                   s1<-rev(penalty)[1]*2/d2
                   s2<-rev(penalty)[2]*2/d2
                   beta<-(s2-paras[2])/(s2-s1)*beta+(paras[2]-s1)/(s2-s1)*beta2
                   beta<-beta*d2
                   break
               }
		if(any(drops)) {
			dropid <- seq(drops)[drops]
			for(id in rev(dropid)) {
				R <- downdateR(R, id)
      			}
			dropid <- active[drops]
                        beta[dropid] <- 0
			active <- active[!drops]
			Sign <- Sign[!drops]
                      }
               if(sparse=="varnum" && length(active)>=paras[2]){
                 break
               }
	       if(!is.null(vn))
			names(action) <- vn[abs(action)]
                	actions[[k]] <- action 
	}
             }
     
        return (beta)
}