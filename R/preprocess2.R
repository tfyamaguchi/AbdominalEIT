############################################################ Kirchhoff's Voltage Law
Kirchhoff <- function(data, interval, fact) {
    average <- mean(rowSums(data))
    exdatanew <- data - average/fact  #35 #37
    # exdatanew <- data[1:32,]-rowMeans(data)
    return(list(expout = exdatanew, status = 0))
}
# 
KirchhoffRev <- function(data, interval) {
    interval <- as.numeric(interval)
    rawdata <- data[1:32, 1:32]
    shifteddata <- matrix(0, 32, 32)
    colmat <- matrix(0, 32, 32)
    for (i in 1:interval) {
        shifteddata <- shifteddata + rawdata[, (0:31 + (i - 1))%%32 + 1]
        colmat <- colmat + diag(32, 32, 32)[, (0:31 + (1 - i))%%32 + 1]
    }
    y <- matrix(c(rowSums(shifteddata), colSums(shifteddata)), 64, 1)
    designmat <- rbind(matrix(interval, 32, 32), colmat)
    x <- lm.fit(designmat, y)$coef
    exdatanew <- data[1:32, 1:32] <- rawdata - matrix(x, 32, 32, byrow = TRUE)
    return(list(expout = exdatanew, status = 0))
}
# 
KirchhoffRev2 <- function(data, interval, fact) {
    interval <- as.numeric(interval)
    rawdata <- data[1:32, 1:32]
    shifteddata <- rawdata
    colmat <- matrix(0, 32, 32)
    y <- matrix(rowSums(shifteddata), 32, 1)
    # fact <- 32
    designmat <- matrix(fact, 32, 1)
    x <- lm.fit(designmat, y)$coef
    exdatanew <- rawdata - matrix(x, 32, 32, byrow = TRUE)
    return(list(expout = exdatanew, status = 0))
}
# 
KirchhoffRec1 <- function(data, interval) {
    wdata <- data
    rmat <- matrix(0, 1024, 32)
    for (rc in 1:1024) {
        rmat[rc, ceiling(rc/32)] <- 1
    }
    wrmat <- rmat
    i <- 1
    while (i < interval) {
        wdata <- wdata + data[, c((i + 1):32, 1:i)]
        wrmat <- wrmat + rmat[, c((32 - i + 1):32, 1:(32 - i))]
        i <- i + 1
    }
    lmat <- matrix(0, 1024, 32)
    for (rc in 1:1024) {
        lmat[rc, (rc - 1)%%32 + 1] <- 1
    }
    wlmat <- lmat * interval
    designmat <<- cbind(as.vector(as.matrix(wdata)), wlmat, wrmat)
    designmat2 <- cbind(as.vector(as.matrix(data)), lmat, rmat)
    LL <- function(x) {
        datahat <- designmat %*% (c(1, x))
        datahatmat <- matrix(datahat, 32, 32)
        b <- c(rowMeans(datahatmat), colMeans(datahatmat))
        return(sum(b^2))
    }
    res <- optim(rep(0, 64), LL, , method = "SANN", control = list(trace = 5, maxit = 10000))
    datanew <- designmat2 %*% (c(1, res$par))
    rm(designmat, envir = .GlobalEnv)
    cat("Parameters", res$convergence, "\n")
    return(list(expout = matrix(datanew, 32, 32), status = 0))
}
# 
KirchhoffRec4 <- function(data, interval) {
    # data <- data[,c(32,1:31)]
    wdata <- data
    rmat <- matrix(0, 1024, 32)
    for (rc in 1:1024) {
        rmat[rc, ceiling(rc/32)] <- -1
    }
    wrmat <- rmat
    i <- 1
    while (i < interval) {
        wdata <- wdata + data[, c((i + 1):32, 1:i)]
        wrmat <- wrmat + rmat[, c((32 - i + 1):32, 1:(32 - i))]
        i <- i + 1
    }
    lmat <- matrix(0, 1024, 32)
    for (rc in 1:1024) {
        lmat[rc, (rc - 1)%%32 + 1] <- -1
    }
    wlmat <- lmat * interval
    designmat <<- cbind(as.vector(as.matrix(wdata)), wlmat, wrmat)
    designmat2 <<- cbind(as.vector(as.matrix(data)), lmat, 1 * rmat)
    LL <- function(x) {
        datahat <- designmat %*% (c(1, x[1:64]))
        datahatmat <- matrix(datahat, 32, 32)
        b1 <- rowSums(datahatmat)
        b2 <- colSums(datahatmat)
        ll <- c(dnorm(b1, mean = 0, sd = x[65], log = TRUE))
        return(sum(ll))
    }
    # res <- optim(c(rep(0,64),0,0.1),LL,control=list(fnscale=-1,maxit=10000))
    res <- optim(c(rep(0, 64), 10), LL, method = "Nelder-Mead", control = list(fnscale = -1, 
        parscale = c(rep(1, 64), 10), maxit = 50000, trace = 5))
    datanew <- designmat2 %*% (c(1, res$par[1:64]))
    datamat <- matrix(datanew, 32, 32)
    # datamat <- datamat[,c(2:32,1)] rm(designmat,envir=.GlobalEnv)
    cat("Parameters", res$par[65], res$convergence, "\n")
    cat("Parameters", res$par[1:64], "\n")
    return(list(expout = datamat, par = res$par, status = 0))
}
# 
KirchhoffRec3 <- function(data, interval) {
    wdata <- data
    rmat <- matrix(0, 1024, 32)
    for (rc in 1:1024) {
        rmat[rc, ceiling(rc/32)] <- -1
    }
    wrmat <- rmat
    i <- 1
    while (i < interval) {
        wdata <- wdata + data[, c((i + 1):32, 1:i)]
        wrmat <- wrmat + rmat[, c((32 - i + 1):32, 1:(32 - i))]
        i <- i + 1
    }
    lmat <- matrix(0, 1024, 32)
    for (rc in 1:1024) {
        lmat[rc, (rc - 1)%%32 + 1] <- -1
    }
    wlmat <- lmat * interval
    designmat <<- cbind(as.vector(as.matrix(wdata)), wlmat, wrmat)
    designmat2 <<- cbind(as.vector(as.matrix(data)), 1 * lmat, 1 * rmat)
    LL <- function(x) {
        datahat <- designmat %*% (c(1, x[1:64]))
        datahatmat <- matrix(datahat, 32, 32)
        b1 <- rowSums(datahatmat)
        b2 <- colSums(datahatmat)
        ll <- c(dnorm(b1, mean = x[65], sd = x[66], log = TRUE), 0.8 * dnorm(b2, 
            mean = x[65], sd = x[66], log = TRUE))
        return(sum(ll))
    }
    # res <- optim(c(rep(0,64),0,0.1),LL,control=list(fnscale=-1,maxit=10000))
    res <- optim(c(rep(0, 64), 0, 10), LL, method = "Nelder-Mead", control = list(fnscale = -1, 
        parscale = c(rep(1, 64), 1, 1), maxit = 90000, trace = 5))
    res <- optim(res$par, LL, method = "Nelder-Mead", control = list(fnscale = -1, 
        parscale = c(rep(0.5, 64), 1, 1), maxit = 90000, trace = 5))
    res <- optim(res$par, LL, method = "Nelder-Mead", control = list(fnscale = -1, 
        parscale = c(rep(1, 64), 1, 1), maxit = 90000, trace = 5))
    datanew <- designmat2 %*% (c(1, res$par[1:64]))
    datamat <- matrix(datanew, 32, 32)
    # rm(designmat,envir=.GlobalEnv)
    cat("Parameters", res$par[65:66], res$convergence, "\n")
    cat("Parameters", res$par[1:64], "\n")
    return(list(expout = datamat, par = res$par, status = 0))
}
# 
KirchhoffRec2 <- function(data, interval) {
    wdata <- data
    rmat <- diag(32, 32, 32)
    wrmat <- rmat
    i <- 1
    while (i < interval) {
        wrmat <- wrmat + rmat[, c(((-i)%%32 + 1):32, 1:(32 - i))]
        i <- i + 1
    }
    lmat <- matrix(1, 32, 32)
    wlmat <- lmat * interval
    dmat <- rbind(wlmat, wrmat)
    b <- c(rowSums(wdata), colSums(wdata))
    a <- ginv(t(dmat) %*% dmat) %*% t(dmat) %*% b
    datahatmat <- wdata + matrix(rep(a, 32), 32, 32, byrow = TRUE)
    return(list(expout = datahatmat, status = 0))
}

################################################# Exploration 32 electrode system to 64 system
lst <- function(data, bound32) {
    theta <- rep(0, 129)
    theta[2:66] <- bound32$boundry[2:66, 4]
    theta[65:129] <- theta[1:65] + 360
    expout <- NULL
    expout <- rbind(expout, theta[1:65])
    for (ll in 1:32) {
        vtemp <- data$vdatan32[ll, 1:32]
        z <- rep(0, 65)
        z[seq(1, 65, 2)] <- c(vtemp[32], vtemp[1:32])
        z[seq(2, 64, 2)] <- c(vtemp[1:32])
        expout <- rbind(expout, z)
    }
    return(list(expout = expout, status = 0))
}
############################################################ Reciplocal theorem
reciplocal <- function(realdata) {
    ex <- as.matrix(realdata$data[1:32, 1:32])
    interval <- realdata$interval
    # print(ex)
    diagmat <<- matrix(0, nrow(ex), ncol(ex))
    for (i in 1:interval) {
        diagmat <<- diagmat + (ex[, (i + (1:32) - 2)%%32 + 1])
    }
    exnew <- (diagmat + t(diagmat))/2
    if (interval == 1) 
        realdata$data[1:32, 1:32] <- exnew
    return(realdata)
}
############################################################ 
doit <- function(realdata, bound, abs, ang, fact) {
    localdata <- NULL
    for (i in 1:length(realdata)) {
        realpart <- realdata[[i]]$data[, 1:32]
        realpart <- Kirchhoff(data = realdata[[i]]$data[, 1:32], interval = realdata[[i]]$interval, 
            fact = fact)$expout
        if (abs) {
            imaginarypart <- Kirchhoff(data = realdata[[i]]$data[, 33:64], interval = realdata[[i]]$interval, 
                fact = fact)$expout
            zz <- complex(real = realpart, imaginary = imaginarypart)
            zz.shift <- complex(modulus = Mod(zz), argument = Arg(zz) + pi * ang[i]/180)
            realdata[[i]]$data <- matrix(Re(zz.shift), 32, 32)
        } else realdata[[i]]$data <- realpart
        colnames(realdata[[i]]$data) <- paste("V", 1:32, sep = "")
        localdata <- rbind(localdata, realdata[[i]]$data)
    }
    # browser()
    localdata <- realdata[[1]]$data  ## added 2015/4/8
    extemp <- localdata
    write.table(extemp, file = "VDATAN32.txt")
    ########################################################### Smoothing ret <- hosein(data=exdata,klim0=15,ksmooth=18,ksmooth2=27) ret$status
    exdatal <- NULL
    for (i in 1:length(realdata)) exdatal <- c(exdatal, list(list(data = extemp[((i - 
        1) * 32) + (1:32), ], interval = realdata[[i]]$interval)))
    # cl <- makeCluster(detectCores()) system.time(exppar <<-
    # parLapply(cl,exdatal,reciplocal)) system.time(exppar <<-
    # lapply(exdatal,reciplocal))
    retnew <- exdatal
    # stopCluster(cl) browser() retnew <- lapply(exppar,preprocess1,parallel=TRUE)
    plotelect(retnew[[1]]$data[1:32, 1:32], shift = TRUE, log = FALSE)
    expnew <- lapply(retnew, lstsqrsAR, bound32 = bound)
    # expnew <<- lapply(retnew,lst,bound32=boundG)
    return(list(cumdata = expnew, potdata = retnew))
    # return(list(cumdata=retnew,potdata=retnew))
}
################################################################# Preprocessing2
preprocess2 <- function(realdata, bound, abs, ang = c(9, 18), fact = 32) {
    res <- doit(realdata = realdata, bound = bound, abs = abs, ang = ang, fact = fact)
    # browser()
    return(list(cumdata = res$cumdata, potdata = res$potdata))
}
#################################### lstsqrsAR integration
lstsqrsAR <- function(data, bound32) {
    mbc <- 32
    vtemp <- rep(NA, length = 2 * mbc + 1)
    theta <- seq(0, 360, by = 360/64)
    # 
    int <- data$interval
    expout <- NULL
    expout <- rbind(expout, theta[1:65])
    zz <- rep(NA, 65)
    y <- rep(NA, 65)
    # browser()
    for (ll in 1:32) {
        vtemp[1:64] <- c(unlist(data$data[ll, 1:32]), unlist(data$data[ll, 1:32]))
        y[1] <- vtemp[1]
        for (i in 2:64) {
            y[i] <- y[i - 1] + vtemp[i]
        }
        z <- y[1:32]
        z <- z - mean(z)
        x <- seq(0, 2 * pi, length.out = 33)[1:32]
        pispl <- periodicSpline(x, z, period = 2 * pi, ord = 2)
        zy <- predict(pispl, seq(0, 2 * pi, length.out = 65))$y[c(seq(1, 63, 2))]
        # browser() zz[seq(2,64,2)] <- -zy[c(1:32)]
        zz[seq(2, 64, 2)] <- -z[c(1:32)]
        expout <- rbind(expout, zz)
    }
    # browser()
    return(list(expout = expout, status = 0))
}
####################################  
