############################################################################ 
#' Mesh generation
#' @param rad numeric vectors. Shape parameter vector (<1.0). 
#' @param zrange numeric. Length of the model.
#' @param mxbnd numeric. Maximum bandwidth.
#' @param mxnode numeric. Number of nodes.
#' @param bndd numeric.
#' @param mxn numeric.
#' @param mzn numeric. 
#' @param mxb numeric.
#' @param mze numeric.
#' @param mxe numeric.
#' @param zcod numeric.
#' @param nzin numeric.
#' @param nzot numeric.
#' @param zcalc numeric.
#' @param layer numeric vector.
meshgen <- function(rad, zrange, mxbnd, mxnode, bndd, mxn, mzn, mxb, mze, mxe, zcod, 
    nzin, nzot, zcalc = TRUE, layer = c(0.175, 0.295, 0.4, 0.5, 0.55)) {
    pb6 <- txtProgressBar(min = 1, max = 6, style = 3)
    #browser()
    setTxtProgressBar(pb6, 1)
    xxl <- vector("numeric", length = mxn)
    yyl <- vector("numeric", length = mxn)
    nex <- 0
    nez <- 0
    nodex <- array(0, c(4, mxe, mze))
    nodez <- array(0, c(4, mxe, mze))
    nnodex <- 0
    nnodez <- 0
    nb <- 0
    ibnd <- vector("numeric", length = mxb + 2)
    out <- .Fortran("b9c3d128", as.vector(as.double(rad)), xcod2 = as.vector(as.double(xxl)), 
        ycod2 = as.vector(as.double(yyl)), as.integer(nex), as.integer(nez), as.vector(as.integer(nodex)), 
        as.vector(as.integer(nodez)), as.integer(nnodex), as.integer(nnodez), as.integer(nb), 
        as.vector(as.integer(ibnd)), as.double(layer[1]), as.double(layer[2]), as.double(layer[3]), 
        as.double(layer[4]), PACKAGE = "AbdominalEIT")  #,as.double(layer[5]))
    # browser()
    setTxtProgressBar(pb6, 2)
    # h9m5 <- (as.matrix(read.table('H9M5DATA.QQQ',fill=TRUE,colClasses='numeric')))
    out2 <- inputdP2(zrange = zrange, mxbnd = mxbnd, mxnode = mxnode, bndd = bndd, 
        mxn = mxn, mzn = mzn, mxb = mxb, mze = mze, mxe = mxe, zcalc = zcalc, xcod = as.double(out[[2]]), 
        ycod = as.double(out[[3]]), nex = as.integer(out[[4]]), nez = as.integer(out[[5]]), 
        nodex = array(as.integer(out[[6]]), c(4, mxe, mze)), nodez = array(as.integer(out[[7]]), 
            c(4, mxe, mze)), nnodex = as.integer(out[[8]]), nnodez = as.integer(out[[9]]), 
        nb = as.integer(out[[10]]), ibnd = as.vector(as.integer(out[[11]])), zcod = zcod, 
        nzin = nzin, nzot = nzot, pb = pb6)
    setTxtProgressBar(pb6, 6)
    cat("\n")
    return(out2)
}
########################################################################### SUBROUTINE INPUTD
inputdP2 <- function(zrange, mxbnd, mxnode, bndd, mxn, mzn, mxb, mze, mxe, zcalc, 
    xcod, ycod, nex, nez, nodex, nodez, nnodex, nnodez, nb, ibnd, zcod, nzin, nzot, 
    pb) {
    # browser()
  library(parallel)
    nodex <- aperm(nodex, c(1, 3, 2))
    nodez <- aperm(nodez, c(1, 3, 2))
    theta <- vector("numeric", length = mxb + 1)
    # st <- array(0, c(4,4,mze,mxe)) stz <- array(0, c(4,4,mze,mxe)) vol <- array(0,
    # c(mze,mxe)) cex <- array(0, c(3,4,mze,mxe)) xav <- vector('numeric',
    # length=mxe) yav <- vector('numeric', length=mxe)
    arar <- vector("numeric", length = mxe)
    sumr <- vector("numeric", length = mxb + 1)
    xb <- vector("numeric", length = mxb + 1)
    yb <- vector("numeric", length = mxb + 1)
    # ibnd <- vector('numeric', length=mxb+2)
    bnd <- vector("numeric", length = mxb + 1)
    # nodex <- array(0, c(4,mze,mxe)) nodez <- array(0, c(4,mze,mxe))
    xcod2 <- xcod
    ycod2 <- ycod
    # nex <- round(mesh[1,1]) nez <- round(mesh[1,2])
    l <- 1
    # for(k in 1:nex){ for(j in 1:nez) { iex <- round(mesh[l+1,1]) iez <-
    # round(mesh[l+1,2]) nodex[1:4,j,k] <- round(mesh[l+1,1+2*(1:4)]) nodez[1:4,j,k]
    # <- round(mesh[l+1,2+2*(1:4)]) if((j != iez) | (k != iex)) { cat('ELEMENT NUMBER
    # INCONSISTEMT \n') cat('l,j,iez,k,iex',l,j,iez,k,iex,'\n') iret <- 4 stop()
    # return(iret) } l <- l+1 } } nnodex <- round(mesh[l+1,1]) #419 nnodez <-
    # round(mesh[l+1,2]) #13
    l <- l + 1
    nnodezx <- nnodex * nnodez
    if (nnodezx > mxnode) {
        cat(nnodezx, mxnode, "\n")
        cat("NNODEZX EXCEEDS MXNODE \n")
        iret <- 5
        stop()
        return(iret)
    }
    # for (j in 1:nnodex) { node <- round(mesh[l+1,1]) xcod[j] <-
    # as.double(mesh[l+1,2]) ycod[j] <- as.double(mesh[l+1,3]) l <- l+1 if(j != node)
    # { cat('X-NODE NUMBER INCONSISTENT \n') iret <- 6 stop() return(iret) } }
    # browser() for (j in 1:nnodez) { node <- round(mesh[l+1,1]) zcod[j] <-
    # mesh[10*l+2] l <- l+1 if (j != node) { cat('Z-NODE NUMBER INCONSISTENT \n')
    # iret <- 7 stop() return(iret) } } nb <- round(mesh[l+1,1]) l <- l+1 for (j in
    # 1:nb) { ibnd[j] <- round(mesh[l+1,1]) l <- l+1 }
    ibnd[nb + 1] <- ibnd[1]
    #### RECAPTURE TO THE REAL BOUNDARY #######################################
    nbound <- round(bndd[1])
    if (nb != 2 * nbound) {
        cat("INCONSISTENT BOUNDARY SIZE \n")
        cat(nb, "\n")
        cat(nbound, "\n")
        iret <- 8
        stop()
        return(iret)
    }
    sumxx <- 0
    sumyy <- 0
    xmax <- 0
    xmin <- 0
    ymax <- 0
    ymin <- 0
    for (j in 1:nbound) {
        numb <- bndd[j * 4 + 1]
        xb[2 * j - 1] <- bndd[j * 4 + 2]
        yb[2 * j - 1] <- bndd[j * 4 + 3]
        if (xmax < xb[2 * j - 1]) 
            xmax <- xb[2 * j - 1]
        if (xmin > xb[2 * j - 1]) 
            xmin <- xb[2 * j - 1]
        if (ymax < yb[2 * j - 1]) 
            ymax <- yb[2 * j - 1]
        if (ymin > yb[2 * j - 1]) 
            ymin <- yb[2 * j - 1]
        sumxx <- sumxx + xb[2 * j - 1]
        sumyy <- sumyy + yb[2 * j - 1]
    }
    size <- zrange * (xmax + ymax - xmin - ymin)/4
    sumxx <- sumxx/nbound
    sumyy <- sumyy/nbound
    xb[2 * (1:nbound) - 1] <- xb[2 * (1:nbound) - 1] - sumxx  #########
    yb[2 * (1:nbound) - 1] <- yb[2 * (1:nbound) - 1] - sumyy  #########
    xb[2 * nbound + 1] <- xb[1]
    yb[2 * nbound + 1] <- yb[1]
    # 
    # xb[2*(1:nbound)-1] <- (xb[(2*(1:nbound)-3)%%(2*nbound)+1]+xb[2*(1:nbound)+0])/2
    # yb[2*(1:nbound)-1] <- (yb[(2*(1:nbound)-3)%%(2*nbound)+1]+yb[2*(1:nbound)+0])/2
    xb[2 * (1:nbound)] <- (xb[2 * (1:nbound) - 1] + xb[2 * (1:nbound) + 1])/2
    yb[2 * (1:nbound)] <- (yb[2 * (1:nbound) - 1] + yb[2 * (1:nbound) + 1])/2
    # xb[2*nbound+1] <- xb[1] yb[2*nbound+1] <- yb[1]
    sumr[1] <- 0
    for (j in 1:nbound) {
        rr <- sqrt((xb[2 * j + 1] - xb[2 * j - 1])^2 + (yb[2 * j + 1] - yb[2 * j - 
            1])^2)
        sumr[2 * j + 1] <- sumr[2 * j - 1] + rr
    }
    dsize <- sumr[2 * nbound + 1]/nb
    f360 <- 2 * pi/sumr[2 * nbound + 1]
    c360 <- 360/sumr[2 * nbound + 1]
    # cat('f360 c360',f360,c360,'\n')
    sumr[2 * (1:nbound)] <- (sumr[2 * (1:nbound) - 1] + sumr[2 * (1:nbound) + 1])/2
    theta[1:(2 * nbound + 1)] <- sumr[1:(2 * nbound + 1)] * c360
    sumr[1:(2 * nbound + 1)] <- sumr[1:(2 * nbound + 1)] * f360
    setTxtProgressBar(pb, 3)
    ##### SUMR:
    ##### #########################################
    #browser()
    gc()
    cl5 <- parallel::makeCluster(parallel::detectCores(),type="PSOCK")  
    registerDoParallel(cl5)
    xycod <- foreach(j = 1:nnodex, .combine = rbind) %dopar% {
        # zz <- file('recapR.txt','w') for (j in 1:nnodex) {
        xxx <- xcod[j]
        yyy <- ycod[j]
        # if((abs(xxx) < 0.0001) & (abs(yyy) < 0.0001)) next
        if ((abs(xxx) < 1e-04) & (yyy < 0)) {
            thetat <- 0
        } else {
            thetat <- (atan2(-xxx, yyy) + pi)
        }
        rrr <- sqrt(xxx^2 + yyy^2)
        # browser() for (k in 1:nb) {
        k <- which((thetat >= sumr[1:nb]) & (thetat < sumr[1:nb + 1]))  #{
        if (length(k) != 1) 
            stop()
        # cat(thetat, sumr[k],sumr[k+1],'\n')
        xxt <- ((thetat - sumr[k]) * xb[k + 1] + (sumr[k + 1] - thetat) * xb[k])/(sumr[k + 
            1] - sumr[k])
        yyt <- ((thetat - sumr[k]) * yb[k + 1] + (sumr[k + 1] - thetat) * yb[k])/(sumr[k + 
            1] - sumr[k])
        xcod[j] <- rrr * xxt
        ycod[j] <- rrr * yyt
        # } } cat(j,xxt,yyt,xcod[j],ycod[j],rrr,'\n',file=zz)
        c(xcod[j], ycod[j])
    }
    # browser()
    stopCluster(cl5)
    xcod <- xycod[, 1]
    ycod <- xycod[, 2]
    zspace <- 15  #15.24
    nbw <- NULL
    if (zcalc) {
        ret <- zmesh(zrange = zrange, mzn = mzn, nnodez = nnodez, zspace = zspace, 
            dsize = dsize, nex = nex, nez = nez, nodex = nodex, nodez = nodez, mxbnd = mxbnd)
        zcod <- ret$zcod
        nzin <- ret$nzin
        nzot <- ret$nzot
        nbw <- ret$nbw
    }
    ########## BASIC AREA OF CURRENT POLES BND(1)
    ########## ################## tate <- zspace/8
    tate <- (zcod[nzin] - zcod[nzin - 1])/4
    xxs <- xcod[ibnd[nb]]
    yys <- ycod[ibnd[nb]]
    xxt <- xcod[ibnd[1]]
    yyt <- ycod[ibnd[1]]
    bnd[1] <- sqrt((xxs - xxt)^2 + (yys - yyt)^2) * tate
    for (j in 2:(nb + 1)) {
        xxs <- xxt
        yys <- yyt
        xxt <- xcod[ibnd[j]]
        yyt <- ycod[ibnd[j]]
        bnd[j] <- sqrt((xxs - xxt)^2 + (yys - yyt)^2) * tate
    }
    ########## REAL ELEMENT AREA #############################################
    setTxtProgressBar(pb, 4)
    #browser()
    gc()
    cl6 <- parallel::makeCluster(parallel::detectCores(),type="PSOCK")  
    registerDoParallel(cl6)
    arar <- foreach(i = 1:nex, .combine = "c") %dopar% {
        n1 <- nodex[1, 4, i]
        n2 <- nodex[2, 4, i]
        n3 <- nodex[3, 4, i]
        x1 <- xcod[n1]
        x2 <- xcod[n2]
        x3 <- xcod[n3]
        y1 <- ycod[n1]
        y2 <- ycod[n2]
        y3 <- ycod[n3]
        # xav[i] <- (x1+x2+x3)/3 yav[i] <- (y1+y2+y3)/3
        wmat <- matrix(c(x1, x2, x3, y1, y2, y3, 1, 1, 1), 3, 3)
        area <- det(wmat)
        # area <- x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3
        area/2
        # cat(i,area,'\n',file=zz)
    }
    #browser()
    stopCluster(cl6)  #close(zz)
    ########## STIFFNESS MATRIX ###############################################
    setTxtProgressBar(pb, 5)
    gc()
    cl7 <- parallel::makeCluster(parallel::detectCores(),type="PSOCK")  
    #ret <- apply(        cbind(rep(1:nez,nex),  rep(1:nex, each = nez)), 1, FUN = mcex,xcod = xcod, ycod = ycod, zcod = zcod, nodex = nodex, nodez = nodez)
    ret <- parApply(cl7, cbind(rep(1:nez, nex), rep(1:nex, each = nez)), 1, FUN = mcex,xcod = xcod, ycod = ycod, zcod = zcod, nodex = nodex, nodez = nodez)
    #plot(tm)
    stopCluster(cl7)
    # browser()
    ans <- unlist(ret, recursive = FALSE)
    cex <- array(unlist(ans[seq(1, (nez * nex * 5 - 4), 5)]), c(3, 4, mze, mxe))
    st <- array(unlist(ans[seq(2, (nez * nex * 5 - 3), 5)]), c(4, 4, mze, mxe))
    stxy <- array(unlist(ans[seq(3, (nez * nex * 5 - 2), 5)]), c(4, 4, mze, mxe))
    stz <- array(unlist(ans[seq(4, (nez * nex * 5 - 1), 5)]), c(4, 4, mze, mxe))
    vol <- array(unlist(ans[seq(5, (nez * nex * 5), 5)]), c(mze, mxe))
    # 
    return(list(theta = theta, cex = cex, st = st, stxy = stxy, stz = stz, vol = vol, 
        arar = arar, bnd = bnd, ibnd = ibnd, xcod = xcod, ycod = ycod, zcod = zcod, 
        xcod2 = xcod2, ycod2 = ycod2, nodex = nodex, nodez = nodez, nb = nb, nbw = nbw, 
        nnodex = nnodex, nnodez = nnodez, nnodezx = nnodezx, nex = nex, nez = nez, 
        nzin = nzin, nzot = nzot, xb = xb, yb = yb))
}
# 
zmesh <- function(zrange, mzn, nnodez, zspace, dsize, nex, nez, nodex, nodez, mxbnd) {
    #################### Z-COORDINATE INPUT ##################################
    zcod <- vector("numeric", length = mzn)
    nzin <- 8  #9
    nzot <- 6  #5
    # zmin <- -zrange zmax <- zrange REARRANGE OF Z CORDINATES
    # ###################################### if (nzin-nzot != 2) stop if (zmax-zmin <
    # 2*zspace) stop
    nzcent = (nzin + nzot)/2
    zcod[nzcent] <- 0
    zcod[nzin] <- zspace/2
    zcod[nzot] <- -zcod[nzin]
    if ((nnodez < nzin + 1) || (nzot < 2)) 
        stop
    zcod[nzin + 1] <- zcod[nzin] + zspace/2  #dsize
    zcod[nzot - 1] <- zcod[nzot] - zspace/2  #dsize
    if (nnodez > nzin + 1) {
        # if (zmax <= zspace) stop
        nk <- nnodez - nzin - 1
        summ <- sum(2^((1:nk) - 1))
        zbase <- (zrange - zspace)/summ
        cat("inn:zbase, zspace/2", zbase, zspace/2, "\n")
        for (j in 1:nk) {
            zcod[nzin + j + 1] <- zcod[nzin + j] + zbase * 2^(j - 1)
            zcod[nzot - j - 1] <- zcod[nzot - j] - zbase * 2^(j - 1)
        }
    }
    # if (nzot > 2) { if (zmin >= -zspace) stop nk <- nzot-2 summ <-
    # sum(3^((1:nk)-1)) zbase <- (zrange-15)/summ cat('out: zbase,
    # zspace/2',zbase,zspace/2,'\n') zcod[nzot-2] <- zcod[nzcent] - zbase*3^2 for (j
    # in ((2:nk)-1)) zcod[nzot-j-2] <- zcod[nzot-j-1] - zbase*3^(j-0) }
    for (j in 1:mzn) cat(j, zcod[j], "\n")
    ########## BANDWIDTH ##################################################### nbwmat <-
    ########## matrix(0,nez,nex) for (i in 1:nex){ for(ii in 1:nez){ nabs1 <-
    ########## nnodez*(nodex[1:4,ii,i]-1)+nodez[1:4,ii,i] nabs2 <-
    ########## nnodez*(nodex[1:4,ii,i]-1)+nodez[1:4,ii,i] nbwmat[ii,i] <-
    ########## max(abs(outer(nabs1,nabs2,FUN='-'))) } }
    bw <- function(x, nnodez, nodex, nodez) {
        i <- x[1]
        ii <- x[2]
        nabs1 <- nnodez * (nodex[1:4, ii, i] - 1) + nodez[1:4, ii, i]
        nabs2 <- nnodez * (nodex[1:4, ii, i] - 1) + nodez[1:4, ii, i]
        return(max(abs(outer(nabs1, nabs2, FUN = "-"))))
    }
    # browser()
    dmat <- cbind(rep(1:nex, nez), rep(1:nez, each = nex))
    #cl <- makeCluster(detectCores(),type="SOCK")
    nbwvec <- apply(X=dmat,MARGIN=1,FUN=bw,nnodez=nnodez,nodex=nodex,nodez=nodez)
    #nbwvec <- parApply(cl, X = dmat, MARGIN = 1, FUN = bw, nnodez = nnodez, nodex = nodex, nodez = nodez)
    #stopCluster(cl)
    nbw <- max(nbwvec) + 1
    cat("Half bandwidth= ", nbw, "\n")
    if (mxbnd < nbw) {
        cat("REWRITE MXBND=", mxbnd, "\n")
        iret <- 13
        stop()
    } else iret <- 0
    ######################################### 
    return(list(zcod = zcod, nzin = nzin, nzot = nzot, nbw = nbw))
}
# 
mcex <- function(X, xcod, ycod, zcod, nodex, nodez) {
    # x <- vector(length=4) y <- vector(length=4) z <- vector(length=4)
    cex <- matrix(0, 3, 4)
    # st <- stxy <- stz <- matrix(0,4,4)
    l <- X[1]
    n <- X[2]
    x <- xcod[nodex[1:4, l, n]]
    y <- ycod[nodex[1:4, l, n]]
    z <- zcod[nodez[1:4, l, n]]
    # browser()
    vol <- det(cbind(rep(1, 4), x, y, z))/6
    # d1 <-
    # x[2]*y[3]*z[4]+y[2]*z[3]*x[4]+z[2]*x[3]*y[4]-z[2]*y[3]*x[4]-x[2]*z[3]*y[4]-y[2]*x[3]*z[4]
    # d2 <-
    # x[3]*y[4]*z[1]+y[3]*z[4]*x[1]+z[3]*x[4]*y[1]-z[3]*y[4]*x[1]-x[3]*z[4]*y[1]-y[3]*x[4]*z[1]
    # d3 <-
    # x[4]*y[1]*z[2]+y[4]*z[1]*x[2]+z[4]*x[1]*y[2]-z[4]*y[1]*x[2]-x[4]*z[1]*y[2]-y[4]*x[1]*z[2]
    # d4 <-
    # x[1]*y[2]*z[3]+y[1]*z[2]*x[3]+z[1]*x[2]*y[3]-z[1]*y[2]*x[3]-x[1]*z[2]*y[3]-y[1]*x[2]*z[3]
    # vol <- d1-d2+d3-d4 cex[1,1] <-
    # -y[3]*z[4]-y[4]*z[2]-y[2]*z[3]+z[3]*y[4]+z[4]*y[2]+z[2]*y[3]
    cex[1, 1] <- -det(cbind(rep(1, 3), y[2:4], z[2:4]))
    # cex[1,2] <- +y[4]*z[1]+y[1]*z[3]+y[3]*z[4]-z[4]*y[1]-z[1]*y[3]-z[3]*y[4]
    cex[1, 2] <- det(cbind(rep(1, 3), y[c(3, 4, 1)], z[c(3, 4, 1)]))
    # cex[1,3] <- -y[1]*z[2]-y[2]*z[4]-y[4]*z[1]+z[1]*y[2]+z[2]*y[4]+z[4]*y[1]
    cex[1, 3] <- -det(cbind(rep(1, 3), y[c(4, 1, 2)], z[c(4, 1, 2)]))
    # cex[1,4] <- +y[2]*z[3]+y[3]*z[1]+y[1]*z[2]-z[2]*y[3]-z[3]*y[1]-z[1]*y[2]
    cex[1, 4] <- det(cbind(rep(1, 3), y[1:3], z[1:3]))
    # cex[2,1] <- -z[3]*x[4]-z[4]*x[2]-z[2]*x[3]+x[3]*z[4]+x[4]*z[2]+x[2]*z[3]
    cex[2, 1] <- -det(cbind(rep(1, 3), z[2:4], x[2:4]))
    # cex[2,2] <- +z[4]*x[1]+z[1]*x[3]+z[3]*x[4]-x[4]*z[1]-x[1]*z[3]-x[3]*z[4]
    cex[2, 2] <- det(cbind(rep(1, 3), z[c(3, 4, 1)], x[c(3, 4, 1)]))
    # cex[2,3] <- -z[1]*x[2]-z[2]*x[4]-z[4]*x[1]+x[1]*z[2]+x[2]*z[4]+x[4]*z[1]
    cex[2, 3] <- -det(cbind(rep(1, 3), z[c(4, 1, 2)], x[c(4, 1, 2)]))
    # cex[2,4] <- +z[2]*x[3]+z[3]*x[1]+z[1]*x[2]-x[2]*z[3]-x[3]*z[1]-x[1]*z[2]
    cex[2, 4] <- det(cbind(rep(1, 3), z[1:3], x[1:3]))
    # cex[3,1] <- -x[3]*y[4]-x[4]*y[2]-x[2]*y[3]+y[3]*x[4]+y[4]*x[2]+y[2]*x[3]
    cex[3, 1] <- -det(cbind(rep(1, 3), x[2:4], y[2:4]))
    # cex[3,2] <- +x[4]*y[1]+x[1]*y[3]+x[3]*y[4]-y[4]*x[1]-y[1]*x[3]-y[3]*x[4]
    cex[3, 2] <- det(cbind(rep(1, 3), x[c(3, 4, 1)], y[c(3, 4, 1)]))
    # cex[3,3] <- -x[1]*y[2]-x[2]*y[4]-x[4]*y[1]+y[1]*x[2]+y[2]*x[4]+y[4]*x[1]
    cex[3, 3] <- -det(cbind(rep(1, 3), x[c(4, 1, 2)], y[c(4, 1, 2)]))
    # cex[3,4] <- +x[2]*y[3]+x[3]*y[1]+x[1]*y[2]-y[2]*x[3]-y[3]*x[1]-y[1]*x[2]
    cex[3, 4] <- det(cbind(rep(1, 3), x[1:3], y[1:3]))
    # vol <- vol/6 #
    cex <- cex/vol/6  ### originally 6
    # vol <- vol/6 for(k in 1:4){ for(j in 1:4){ st[k,j] <-
    # vol*(cex[1,k]*cex[1,j]+cex[2,k]*cex[2,j]+cex[3,k]*cex[3,j]) stxy[k,j] <-
    # vol*(cex[1,k]*cex[1,j]+cex[2,k]*cex[2,j]) stz[k,j] <- vol*cex[3,k]*cex[3,j] } }
    st <- vol * (cex[1, 1:4] %o% cex[1, 1:4] + cex[2, 1:4] %o% cex[2, 1:4] + cex[3, 
        1:4] %o% cex[3, 1:4])
    stxy <- vol * (cex[1, 1:4] %o% cex[1, 1:4] + cex[2, 1:4] %o% cex[2, 1:4])
    stz <- vol * (cex[3, 1:4] %o% cex[3, 1:4])
    return(list(cex = cex, st = st, stxy = stxy, stz = stz, vol = vol))
} 
