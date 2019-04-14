###################### SETUP DEVICE
setupdev <- function(window = 2, clean = TRUE, anisotropy = TRUE) {
    while (length(dev.list()) > 0) dev.off()
    # while (rgl.cur() >0) rgl.close()
    if (length(dev.list()) == 0) {
        dev.new(title = "Convergent status", width = 4, height = 6.3, xpos = 390, 
            ypos = 0, bg = "white")
        par(mfrow = c(2, 1), oma = c(0, 0, 0, 0), mar = c(5, 4, 1, 1))
        if (window >= 1) {
            dev.new(title = "Conductivity image", width = 4, height = 3, xpos = 785, 
                ypos = 0, bg = "white")
            par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
        }
        if (window >= 2) {
            dev.new(title = "Sensitivity image", width = 4, height = 3, xpos = 785, 
                ypos = 330, bg = "white")
            par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
        }
        if (window >= 2) {
            dev.new(title = "Conductivity image (transverse)", width = 4, height = 3, 
                xpos = 785, ypos = 660, bg = "white")
            par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
        }
        if (window >= 2) {
            dev.new(title = "Conductivity image (longitudinal)", width = 4, height = 3, 
                xpos = 390, ypos = 660, bg = "white")
            par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
        }
        if (anisotropy) {
            dev.new(title = "Anisotropy image", width = 4, height = 3, xpos = 0, 
                ypos = 660, bg = "white")
            par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
        } else {
            dev.new(title = "Elements volume & area", width = 4, height = 3, xpos = 0, 
                ypos = 660, bg = "white")
            par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(4, 3, 1, 1))
        }
        if (window >= 2) {
            dev.new(title = "Smoothed image (color)", width = 4, height = 3, xpos = 1180, 
                ypos = 0, bg = "white")
            par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
        }
        if (window >= 2) {
            dev.new(title = "Sensitvity image (Transeverse)", width = 4, height = 3, 
                xpos = 1180, ypos = 330, bg = "white")
            par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
            dev.new(title = "Sensitvity image (Longitudinal)", width = 4, height = 3, 
                xpos = 1180, ypos = 660, bg = "white")
            par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
            dev.new(title = "Mesh image", width = 4, height = 3, xpos = 0, ypos = 330, 
                bg = "white")
            par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
        }
        if (window >= 2) {
            dev.new(title = "Experimental and computed data", width = 4, height = 11, 
                xpos = 1575, ypos = 0, bg = "white")
            par(mfrow = c(4, 1), oma = c(0, 0, 0, 0), mar = c(5, 4, 1, 1))
        }
    }
}
################################################################################## FUNCTION PLOTSTATUS
plotstatusA <- function(filename, devmat, alammat, dd, nee, iter2, sigma, sigmaxy, 
    sigmaz, xcod, ycod, zcod, xcod2, ycod2, nodex, nnodex, nodez, nnodez, nexr, nex, 
    nez, eleint, sense, sensexy, sensez, idval, arar, vol, bound32, aniso, nb, ibnd, 
    initial, estres, xb, yb, aic, annotate = annotate) {
    condmax <- log10(2)
    dev.set(2)
    dev <- devmat[length(devmat)]
    barplot(devmat, ylab = "Residual sum of squares",log="y")
    mtext(paste("Total square error (raw value):", format(devmat[length(devmat)], 
        digit = 1, nsmall = 1), sep = ""), side = 1)
    mtext(paste("Iteration:", iter2), side = 1, padj = 1.5)
    mtext(paste("Interval of current electrodes:", paste(eleint, collapse = ", ")), 
        side = 1, padj = 3)
    # browser()
    plot(dd[[1]], pch = 20, ylab = "Singularvalue", xlab = "", col = "blue", bty = "l")
    # ,ylim=c(0,20))
    if (length(dd) == 2) {
        points(dd[[2]], pch = 20, col = "red")
        legend("topright", c("Transverse", "Longitudinal"), pch = 20, col = c("red", 
            "blue"), bty = "n")
        abline(v = nee$tranc[1], col = "blue")
        abline(v = nee$tranc[2], col = "red")
    } else {
        legend("topright", c("Singular value"), pch = 20, col = c("blue"), bty = "n")
        abline(v = nee$tranc[1], col = "blue")
    }
    mtext(paste("Tikhonov paramter:", round(alammat[length(alammat)], 3)), side = 1, padj = 4, 
        cex = 0.75)
    if (length(nee$tranc) == 2) {
        mtext(paste("Number of sigularvalue:", nee$tranc[1], nee$rank[1], "/", nee$tranc[2], 
            nee$rank[2]), side = 1, padj = 5.5, cex = 0.75)
        mtext(paste("Cumulative proportion of singularvalue:", round(sum(dd$xy[1:nee$tranc[1]])/sum(dd$xy) * 
            100, 1), "% / ", round(sum(dd$z[1:nee$tranc[2]])/sum(dd$z) * 100, 1), 
            "%", sep = ""), side = 1, adj = 0, padj = 7, cex = 0.75)
    } else {
        mtext(paste("Number of sigularvalue:", nee$tranc[1], nee$rank[1]), side = 1, 
            padj = 5.5, cex = 0.75)
        mtext(paste("Cumulative proportion of singularvalue:", round(sum(dd[[1]][1:nee$tranc[1]])/sum(dd[[1]]) * 
            100, 1), "%", sep = ""), side = 1, adj = 0, padj = 7, cex = 0.75)
    }
    idval <- c(idval, list(NumberOfIterations = paste(iter2)))
    idval <- c(idval, list(ErrorComment = paste(round(dev, 1))))
    dev.set(3)
    ret <- condimgP(sigma = sqrt((2 * sigmaxy^2 + sigmaz^2)/3), xcod = xcod, ycod = ycod, 
        zcod = zcod, nodex = nodex, iter2 = iter2, nex = nex, nexr = nexr, arar = arar, 
        threeD = ifelse(iter2%%5 == 0, TRUE, FALSE), annotate = annotate, xb = xb, 
        yb = yb)
    # text(-120,143,paste('RSS,',format(dev,digits=2,nsmall=2),';
    # AIC,',format(aic,digits=1,nsmall=1)),cex=0.7)
    text(-120, 143, paste("Residual sum of squares,", format(dev, digits = 2, nsmall = 2)), 
        cex = 0.7)
    text(55, 143, paste("SFA conductivity,", format(sigma[nexr[9]] * 1000, nsmall = 4, 
        digits = 3), "S/m"), cex = 0.7)
    ani.record()
    # 
    timg <- captureImage()
    idval$InstanceNumber <- "1"
    idval$PhotometricInterpretation = "MONOCHROME1"
    ### writeDICOMFile(dd=timg,fname='cond.dcm',endian='little',value=idval)
    ch <- sample(1:32, 2)
    dev.set(4)
    ret <- senseimg(sense = sense, xcod = xcod, ycod = ycod, nodex = nodex, nex = nex, 
        nexr = nexr, ch = ch, eleint = eleint, arar = arar, barplot = TRUE)
    idval$InstanceNumber <- "2"
    idval$PhotometricInterpretation = "MONOCHROME1"
    # writeDICOMFile(dd=t((img$z-log10(0.01))/(condmax-log10(0.01))*255*3),fname='smoothCond.dcm',endian='little',value=idval)
    dev.set(5)
    if(aniso){
      ret <- condimgP(sigma = sigmaxy, xcod = xcod, ycod = ycod, zcod = zcod, nodex = nodex, 
        iter2 = iter2, nex = nex, nexr = nexr, annotate = FALSE, xb = xb, yb = yb)
    }
    dev.set(6)
    if(aniso){
      ret <- condimgP(sigma = sigmaz, xcod = xcod, ycod = ycod, zcod = zcod, nodex = nodex, 
        iter2 = iter2, nex = nex, nexr = nexr, annotate = FALSE, xb = xb, yb = yb)
    }
    dev.set(7)
    if (aniso) {
        ret <- condimgP(sigma = sigmaz/sigmaxy, xcod = xcod, ycod = ycod, zcod = zcod, 
            nodex = nodex, iter2 = iter2, nex = nex, nexr = nexr, ratio = TRUE, annotate = FALSE, 
            xb = xb, yb = yb)
    } else {
      #  plot(x = 1:nex, y = arar, col = "red", pch = 20, bty = "l", xlab = "Elements", 
      #      ylab = expression(paste("Area (", mm^2, ")")))
      #  points(x = 1:nex, y = colSums(vol)/500, col = "blue", pch = 20)
    }
    dev.set(8)
    filta06(file = filename, sigma = sqrt((2 * sigmaxy^2 + sigmaz^2)/3), bound32 = bound32, 
        nbun = sum(abs(range(bound32$boundry[, "x"], na.rm = TRUE))), rfilt = 0.05, 
        ratmap = 1, col = TRUE, clean = TRUE, nex = nex, nez = nez, nodex = nodex, 
        nodez = nodez, nnodex = nnodex, nnodez = nnodez, xcod = xcod2, ycod = ycod2, 
        nb = nb, ibnd = ibnd)
    # filta06(file=paste('smooth',iter,sep=''),bound32=bound32,nbun=sum(abs(range(xcod))),rfilt=0.1,ratmap=1,col=TRUE,clean=TRUE)
    # img <-
    # smoothimg(sigma=sqrt((2*sigmaxy^2+sigmaz^2)/3),xcod=xcod,ycod=ycod,nodex=nodex,nex=nex,ch=ch,eleint=eleint,color=TRUE)
    # ret <-
    # senseimg(sense=sense,xcod=xcod,ycod=ycod,nodex=nodex,nex=nex,ch=ch,eleint=eleint)
    dev.set(9)
    if (aniso){
      ret <- senseimg(sense = sensexy, xcod = xcod, ycod = ycod, nodex = nodex, nex = nex, 
        nexr = nexr, ch = ch, eleint = eleint, arar = arar)
    }
    dev.set(10)
    if (aniso){
      ret <- senseimg(sense = sensez, xcod = xcod, ycod = ycod, nodex = nodex, nex = nex, 
        nexr = nexr, ch = ch, eleint = eleint, arar = arar)
    }
    dev.set(11)
    # ret <-
    # condimgP(sigma=sqrt((2*sigmaxy^2+sigmaz^2)/3),xcod=xcod,ycod=ycod,zcod=zcod,nodex=nodex,iter=iter,nex=nex,nexr=nexr,arar=arar,mesh=TRUE)
    dev.set(7)
    # browser() par(mfrow=c(2,1),oma=c(0,0,0,0),mar=c(5,4,1,1)) browser() H <-
    # hist(log10(ddata),plot=F,xlim=c(0.01,1)) plot(H$mids,H$counts,type='n',
    # xaxt='n', xlab='Conductivity (S/m)',ylab='Counts', bg='lightgrey' )
    # abline(v=(H$breaks),col='lightgrey',lty=2) abline(v=(H$mids),col='lightgrey')
    # abline(h=pretty(H$counts),col='lightgrey') plot(H,add=T,freq=T,col='blue')
    # hist(ddata,xlim=c(0.01,1),xlog=TRUE,xlab='Conductivity (S/m)')
    # abline(v=log(1000*sigma[(nexr[6]+1):nexr[9]]),col='red')
    estret <- NULL
    if (dev < 1) {
        # browser()
        ddata <- NULL
        for (i in 1:nexr[6]) ddata <- c(ddata, rep(sigma[i] * 1000, each = round(1 * 
            arar[i])))
        dist <- "Beta"
        estret <- fatEstimation(ddata, fn = dist, initial = initial, estres = estres)
        initial <- estret$par
        mtext(paste("Deviation,", format(dev, digit = 1, nsmall = 1)), side = 3, 
            cex = 0.7, outer = FALSE)
        cat(estret$par, "\n")
        dev.set(3)
        cavityarea <- sum(arar[1:nexr[6]])
        text(-120, 130, paste("Cavity area,", format(cavityarea/100, nsmall = 1, 
            digits = 1)), cex = 0.7)
        text(-60, 133, expression(cm^2), cex = 0.7)
        text(55, 130, paste("Visceral fat area,", format(cavityarea/100 * estret$par[ifelse(dist=="Beta",3,5)], 
            nsmall = 1, digits = 1)), cex = 0.7)
        text(128, 133, expression(cm^2), cex = 0.7)
    }
    #################################### 
    return(list(initial = initial, estres = estret))
}
################################################################################# SUBROUTINE CONDIMG
condimgP <- function(sigma, xcod, ycod, zcod, nodex, iter2, nex, nexr, ratio = FALSE, 
    threeD = FALSE, arar = NULL, mesh = FALSE, annotate = annotate, xb, yb) {
    condmax <- log10(1)
    cond <- ifelse(rep(ratio, nex), log10(sigma[1:nex]), log10(1000 * sigma[1:nex]))
    elements <- cbind(COND = cond, P1X = xcod[nodex[1, 4, 1:nex]], P1Y = ycod[nodex[1, 
        4, 1:nex]], P2X = xcod[nodex[2, 4, 1:nex]], P2Y = ycod[nodex[2, 4, 1:nex]], 
        P3X = xcod[nodex[3, 4, 1:nex]], P3Y = ycod[nodex[3, 4, 1:nex]])
    if (ratio) {
        elements[, "COND"] <- ifelse(elements[, "COND"] > log10(100), rep(log10(100), 
            nex), elements[, "COND"])
        elements[, "COND"] <- ifelse(elements[, "COND"] < log10(0.01), rep(log10(0.01), 
            nex), elements[, "COND"])
        elements <- cbind(elements, COL = 256 - 256/(log10(100) - log10(0.01)) * 
            (elements[, "COND"] - log10(0.01)))
    } else {
        elements[, "COND"] <- ifelse(elements[, "COND"] > condmax, rep(condmax, nex), 
            elements[, "COND"])
        elements[, "COND"] <- ifelse(elements[, "COND"] < log10(0.01), rep(log10(0.01), 
            nex), elements[, "COND"])
        elements <- cbind(elements, COL = 256 - 256/(condmax - log10(0.01)) * (elements[, 
            "COND"] - log10(0.01)))
    }
    ccx <- NULL
    ccx2 <- NULL
    ccy <- NULL
    for (eid in 1:nrow(elements)) {
        cx <- c(elements[eid, "P1X"], elements[eid, "P2X"], elements[eid, "P3X"])
        cy <- c(elements[eid, "P1Y"], elements[eid, "P2Y"], elements[eid, "P3Y"])
        cz <- c(zcod[length(zcod)], zcod[length(zcod)], zcod[length(zcod)])
        ccx <- c(ccx, list(cbind(cx - 20, cy)))
        ccx2 <- c(ccx2, list(cbind(cx, cy, cz)))
    }
    # browser() ccolor <-
    # c('brown','red','orange','yellow','blue','green','violet','gray','white','black')
    plot(0, 0, xlim = c(-192, 192), ylim = c(-144, 144), type = "n", axes = FALSE, 
        asp = 1)
    mapply(x = ccx[1:nexr[1]], FUN = polygon, col = gray(elements[1:nexr[1], "COL"]/256), 
        border = ifelse(mesh, "brown", FALSE), lty = "solid")
    mapply(x = ccx[(nexr[1] + 1):nexr[2]], FUN = polygon, col = gray(elements[(nexr[1] + 
        1):nexr[2], "COL"]/256), border = ifelse(mesh, "red", FALSE), lty = "solid")
    mapply(x = ccx[(nexr[2] + 1):nexr[3]], FUN = polygon, col = gray(elements[(nexr[2] + 
        1):nexr[3], "COL"]/256), border = ifelse(mesh, "orange", FALSE), lty = "solid")
    mapply(x = ccx[(nexr[3] + 1):nexr[4]], FUN = polygon, col = gray(elements[(nexr[3] + 
        1):nexr[4], "COL"]/256), border = ifelse(mesh, "yellow", FALSE), lty = "solid")
    mapply(x = ccx[(nexr[4] + 1):nexr[5]], FUN = polygon, col = gray(elements[(nexr[4] + 
        1):nexr[5], "COL"]/256), border = ifelse(mesh, "blue", FALSE), lty = "solid")
    mapply(x = ccx[(nexr[5] + 1):nexr[6]], FUN = polygon, col = gray(elements[(nexr[5] + 
        1):nexr[6], "COL"]/256), border = ifelse(mesh, "green", FALSE), lty = "solid")
    mapply(x = ccx[(nexr[6] + 1):nexr[7]], FUN = polygon, col = gray(elements[(nexr[6] + 
        1):nexr[7], "COL"]/256), border = ifelse(mesh, "violet", FALSE), lty = "solid")
    mapply(x = ccx[(nexr[7] + 1):nexr[8]], FUN = polygon, col = gray(elements[(nexr[7] + 
        1):nexr[8], "COL"]/256), border = ifelse(mesh, "gray", FALSE), lty = "solid")
    mapply(x = ccx[(nexr[8] + 1):nexr[9]], FUN = polygon, col = gray(elements[(nexr[8] + 
        1):nexr[9], "COL"]/256), border = ifelse(mesh, "white", FALSE), lty = "solid")
    rect(xleft = 155, ybottom = -130 + (seq(256, 0, by = -32)[1:8]), xright = 165, 
        ytop = -130 + (seq(256, 0, by = -32)[2:9]), col = gray(seq(0, 1, by = 1/8)[1:8] + 
            1/16))
    if (threeD) {
        # mapply(x=ccx2,FUN=triangles3d,col=gray(elements[,'COL']/256))
    }
    makecolor <- function(x, minc = log10(0.01), maxc = condmax) {
        c <- ((256 - x) * (maxc - minc)/256 + minc)
        c <- 10^c
        cc <- round(c, 2)
        cc3 <- round(c, 0)
        dif <- (nchar(cc) - nchar(cc3))
        ccret <- ifelse(dif == 2, paste(cc, "0", sep = ""), ifelse(dif == 0, paste(cc, 
            ".00", sep = ""), cc))
        ccret <- ifelse(nchar(ccret) > 4, substr(ccret, 1, 4), ccret)
        ccret <- ifelse(substr(ccret, nchar(ccret), nchar(ccret)) == ".", substr(ccret, 
            1, (nchar(ccret) - 1)), ccret)
        return(ccret)
    }
    if (ratio) 
        text(x = 165, y = (seq(256, 0, by = -32) - 130), label = makecolor(seq(0, 
            256, by = 32), minc = log10(0.01), maxc = log10(100)), adj = -0.2, cex = 1) else text(x = 165, y = (seq(256, 0, by = -32) - 130), label = makecolor(seq(0, 
        256, by = 32)), adj = -0.2, cex = 1)
    # text(x=165,y=134,label=ifelse(ratio,expression(10^x),expression(10^x)),cex=0.7)
    text(x = 165, y = 140, label = ifelse(ratio, "", expression(S/m)), cex = 1)
    # lines(x=c(-190,-90),y=c(-130,-130))
    # text(100,-130,paste('iteration;',iter),cex=1.0)
    text(55, -147, paste("Subcutaneous fat area,", format(sum(arar[(nexr[6] + 1):nex])/100, 
        nsmall = 1, digits = 1)), cex = 0.7)
    text(144, -143, expression(cm^2), cex = 0.7)
    tempmat <- cbind(xb[2:129], yb[2:129]) - cbind(xb[1:128], yb[1:128])
    text(55, -135, paste("Circumference,", format(sum(apply(tempmat, 1, FUN = norm, 
        type = "2"))/10, nsmall = 1, digits = 1), "cm"), cex = 0.7)
    # 
    if (annotate) {
        text(-120, -147, paste("Cross-sectional area,", format(sum(arar[1:nex])/100, 
            nsmall = 1, digits = 1)), cex = 0.7)
        text(-36, -143, expression(cm^2), cex = 0.7)
        if (iter2 > 0) 
            text(-165, -115, paste("Iter;", iter2), cex = 1)
    }
    lines(x = c(-190, -190), y = c(-135, -125))
    lines(x = c(-90, -90), y = c(-135, -125))
    arrows(x0 = -170, x1 = -190, y0 = -130, y1 = -130, length = 0.1)
    arrows(x0 = -110, x1 = -90, y0 = -130, y1 = -130, length = 0.1)
    text(-140, -130, "100 mm", cex = 1)
    return(elements)
}
################################################################################# SUBROUTINE SENSEIMG
senseimg <- function(sense, xcod, ycod, nodex, nex, nexr, ch, eleint, arar, barplot = FALSE) {
    # elements
    # <-cbind(COND=sense[,ch[1],ch[2]],P1X=xcod[nodex[1,4,1:nex]],P1Y=ycod[nodex[1,4,1:nex]],P2X=xcod[nodex[2,4,1:nex]],P2Y=ycod[nodex[2,4,1:nex]],P3X=xcod[nodex[3,4,1:nex]],P3Y=ycod[nodex[3,4,1:nex]])
    # browser()
    elements <- cbind(COND = apply(sense, 1, sum)/arar, P1X = xcod[nodex[1, 4, 1:nex]], 
        P1Y = ycod[nodex[1, 4, 1:nex]], P2X = xcod[nodex[2, 4, 1:nex]], P2Y = ycod[nodex[2, 
            4, 1:nex]], P3X = xcod[nodex[3, 4, 1:nex]], P3Y = ycod[nodex[3, 4, 1:nex]])
    elements <- cbind(elements, COL = (1/(max(elements[, "COND"]) - min(elements[, 
        "COND"])) * (elements[, "COND"] - min(elements[, "COND"]))))
    if (barplot == FALSE) {
        ccx <- NULL
        ccy <- NULL
        for (eid in 1:nrow(elements)) {
            cx <- c(elements[eid, "P1X"], elements[eid, "P2X"], elements[eid, "P3X"])
            cy <- c(elements[eid, "P1Y"], elements[eid, "P2Y"], elements[eid, "P3Y"])
            ccx <- c(ccx, list(x = cx))
            ccy <- c(ccy, list(y = cy))
        }
        plot(0, 0, xlim = c(-192, 192), ylim = c(-144, 144), type = "n", axes = FALSE, 
            asp = 1)
        mapply(FUN = polygon, x = ccx, y = ccy, col = gray(elements[, "COL"]), border = FALSE, 
            lty = "blank")
        # text(-120,-130,paste('Currenr
        # electrods,',ch[2],'-',(ch[2]+eleint[[1]]-1)%%32+1),cex=0.8)
        # text(120,-130,paste('Voltage electrods,',ch[1],'-',(ch[1])%%32+1),cex=0.8)
    } else {
        # 
        #sensAvg <- apply(sense, 1, sum)
        sensAvg <- apply(sense, 1, norm,type="f")
        par(mar = c(5, 4, 1, 1))
        layerAvg <- c(mean(sensAvg[1:nexr[1]]) * 2, mean(sensAvg[(nexr[1] + 1):nexr[2]]) * 
            3, mean(sensAvg[(nexr[2] + 1):nexr[3]]) * 4, mean(sensAvg[(nexr[3] + 
            1):nexr[4]]) * 3, mean(sensAvg[(nexr[4] + 1):nexr[5]]) * 4, mean(sensAvg[(nexr[5] + 
            1):nexr[6]]) * 8)
        barplot(layerAvg, xlab = "Layers", ylab = "Sensitivity", names.arg = as.character(1:6), 
            sub = paste("Variance of sensitivities in layers,", round(var(layerAvg), 
                2)), cex.axis = 0.8, cex.sub = 0.8)
    }
    return(elements)
}
################################################################################# SUBROUTINE MESHIMG
meshimg <- function(xcod, ycod, nodex, nexr, nex, arar) {
    elements <- cbind(P1X = xcod[nodex[1, 4, 1:nex]], P1Y = ycod[nodex[1, 4, 1:nex]], 
        P2X = xcod[nodex[2, 4, 1:nex]], P2Y = ycod[nodex[2, 4, 1:nex]], P3X = xcod[nodex[3, 
            4, 1:nex]], P3Y = ycod[nodex[3, 4, 1:nex]])
    elements <- cbind(elements, COL = c(rep(191, nex0), rep(255, (nex - nex0))))
    ccx <- NULL
    ccy <- NULL
    for (eid in 1:nrow(elements)) {
        cx <- c(elements[eid, "P1X"], elements[eid, "P2X"], elements[eid, "P3X"])
        cy <- c(elements[eid, "P1Y"], elements[eid, "P2Y"], elements[eid, "P3Y"])
        ccx <- c(ccx, list(x = cx))
        ccy <- c(ccy, list(y = cy))
    }
    plot(0, 0, xlim = c(-192, 192), ylim = c(-144, 144), type = "n", axes = FALSE, 
        asp = 1)
    mapply(FUN = polygon, x = ccx, y = ccy, col = gray(elements[, "COL"]/255), border = TRUE, 
        lty = "solid")
    # eids <-
    # c(seq(nex0+1,nex0+121,8),seq(nex0+2,nex0+122,8),seq(nex0+8,nex0+128,8),seq(nex0+7,nex0+127,8))
    eids <- c(nex0 + 1, nex0 + 2, nex0 + 3, nex0 + 4)
    elements <- cbind(P1X = xcod[nodex[1, 1, eids]], P1Y = ycod[nodex[1, 1, eids]], 
        P2X = xcod[nodex[2, 1, eids]], P2Y = ycod[nodex[2, 1, eids]], P3X = xcod[nodex[3, 
            1, eids]], P3Y = ycod[nodex[3, 1, eids]])
    elements <- cbind(elements, COL = rep(128, nrow(elements)))
    ccx <- NULL
    ccy <- NULL
    for (eid in 1:nrow(elements)) {
        cx <- c(elements[eid, "P1X"], elements[eid, "P2X"], elements[eid, "P3X"])
        cy <- c(elements[eid, "P1Y"], elements[eid, "P2Y"], elements[eid, "P3Y"])
        ccx <- c(ccx, list(x = cx))
        ccy <- c(ccy, list(y = cy))
    }
    # mapply(FUN=polygon,x=ccx,y=ccy,col=gray(elements[,'COL']/255),border=TRUE,lty='solid')
    # eids <- c(seq(nex0-95,nex0-5,6),seq(nex0-90,nex0,6))
    eids <- c(nex0 - 90, nex0 - 89)
    elements <- cbind(P1X = xcod[nodex[1, 1, eids]], P1Y = ycod[nodex[1, 1, eids]], 
        P2X = xcod[nodex[2, 1, eids]], P2Y = ycod[nodex[2, 1, eids]], P3X = xcod[nodex[3, 
            1, eids]], P3Y = ycod[nodex[3, 1, eids]])
    elements <- cbind(elements, COL = rep(128, nrow(elements)))
    ccx <- NULL
    ccy <- NULL
    for (eid in 1:nrow(elements)) {
        cx <- c(elements[eid, "P1X"], elements[eid, "P2X"], elements[eid, "P3X"])
        cy <- c(elements[eid, "P1Y"], elements[eid, "P2Y"], elements[eid, "P3Y"])
        ccx <- c(ccx, list(x = cx))
        ccy <- c(ccy, list(y = cy))
    }
    # mapply(FUN=polygon,x=ccx,y=ccy,col=gray(elements[,'COL']/255),border=TRUE,lty='solid')
    # pp <- findnode(elements[1,],elements[16,])
    # points(pp[1],pp[2],col='black',pch=20)
    pp <- findnode(elements[1, ], elements[2, ])
    # points(pp[1],pp[2],col='black',pch=20) pp <-
    # findnode(elements[3,],elements[18,]) points(pp[1],pp[2],col='black',pch=20)
    text(-120, -130, paste("Cross-sectional area,", format(sum(arar[1:nex])/100, 
        nsmall = 1, digits = 1)), cex = 0.7)
    text(-43, -127, expression(cm^2), cex = 0.7)
    text(95, -130, paste("Subcutaneous fat area,", format(sum(arar[(nexr[6] + 1):nex])/100, 
        nsmall = 1, digits = 1)), cex = 0.7)
    text(177, -127, expression(cm^2), cex = 0.7)
    return(elements)
}
################################################################################# SUBROUTINE SMOOTHIMAGE
smoothimg <- function(sigma, xcod, ycod, nodex, nex, ch, eleint, color = FALSE) {
    condmax <- log10(2)
    nodesigma <- apply(cbind(1:419), 1, sigmanode, sigma, nodex, nex)
    smoothimg <- interp(x = xcod[1:419], y = ycod[1:419], z = nodesigma, xo = (-191:192), 
        yo = (-143:144), linear = TRUE)
    smoothimg$z[is.na(smoothimg$z)] <- log10(0.001)
    smoothimg$z[(smoothimg$z > condmax)] <- condmax
    smoothimg$x <- smoothimg$x - 20
    plot(0, 0, xlim = c(-192, 192), ylim = c(-144, 144), type = "n", axes = FALSE, 
        asp = 1)
    if (color) 
        localcolor <- jet.colors(128) else localcolor <- colorpanel(128, low = "white", high = "black")
    image(smoothimg, axes = FALSE, asp = ncol(smoothimg)/nrow(smoothimg), col = localcolor, 
        zlim = c(log10(0.01), condmax), bg = localcolor[1], add = TRUE)
    makecolor <- function(x, minc = log10(0.01), maxc = condmax) {
        c <- ((256 - x) * (maxc - minc)/256 + minc)
        c <- 10^c
        cc <- round(c, 2)
        cc3 <- round(c, 0)
        dif <- (nchar(cc) - nchar(cc3))
        ccret <- ifelse(dif == 2, paste(cc, "0", sep = ""), ifelse(dif == 0, paste(cc, 
            ".00", sep = ""), cc))
        ccret <- ifelse(nchar(ccret) > 4, substr(ccret, 1, 4), ccret)
        ccret <- ifelse(substr(ccret, nchar(ccret), nchar(ccret)) == ".", substr(ccret, 
            1, (nchar(ccret) - 1)), ccret)
        return(ccret)
    }
    # 
    rect(xleft = 155, ybottom = -130 + (seq(256, 0, by = -32)[1:8]), xright = 165, 
        ytop = -130 + (seq(256, 0, by = -32)[2:9]), col = localcolor[seq(0, 128, 
            by = 128/8)[8:1] + 128/16])
    text(x = 165, y = (seq(256, 0, by = -32) - 130), label = makecolor(seq(0, 256, 
        by = 32)), adj = -0.2, cex = 1, col = "black")
    text(x = 165, y = 140, label = expression(S/m), cex = 1, col = "black")
    lines(x = c(-190, -190), y = c(-135, -125), col = "black")
    lines(x = c(-90, -90), y = c(-135, -125), col = "black")
    arrows(x0 = -170, x1 = -190, y0 = -130, y1 = -130, length = 0.1, col = "black")
    arrows(x0 = -110, x1 = -90, y0 = -130, y1 = -130, length = 0.1, col = "black")
    text(-140, -130, "100 mm", cex = 1, col = "black")
    return(smoothimg)
}
##################################################################### 3D Plot
show3dmodel <- function(nodex, nodez, xcod, ycod, zcod, nex, nez, xb, yb, iinn, iott, 
    nzin, nzout) {
    # 
    nodex2 <- aperm(nodex, c(1, 3, 2))
    nodez2 <- aperm(nodez, c(1, 3, 2))
    # 
    p1 <- rbind(px = as.vector(xcod[nodex[1:4, 1:3, 1:nex]]), py = as.vector(ycod[nodex[1:4, 
        1:3, 1:nex]]), pz = as.vector(zcod[nodez[1:4, 1:3, 1:nex]]), 1)
    p2 <- rbind(px = as.vector(xcod[nodex[1:4, 4:6, 517:nex]]), py = as.vector(ycod[nodex[1:4, 
        4:6, 517:nex]]), pz = as.vector(zcod[nodez[1:4, 4:6, 517:nex]]), 1)
    p3 <- rbind(px = as.vector(xcod[nodex[1:4, 7:9, 517:nex]]), py = as.vector(ycod[nodex[1:4, 
        7:9, 517:nex]]), pz = as.vector(zcod[nodez[1:4, 7:9, 517:nex]]), 1)
    p4 <- rbind(px = as.vector(xcod[nodex[1:4, 10:12, 517:nex]]), py = as.vector(ycod[nodex[1:4, 
        10:12, 517:nex]]), pz = as.vector(zcod[nodez[1:4, 10:12, 517:nex]]), 1)
    p5 <- rbind(px = as.vector(xcod[nodex[1:4, 13:15, 517:nex]]), py = as.vector(ycod[nodex[1:4, 
        13:15, 517:nex]]), pz = as.vector(zcod[nodez[1:4, 13:15, 517:nex]]), 1)
    p6 <- rbind(px = as.vector(xcod[nodex[1:4, 16:18, 517:nex]]), py = as.vector(ycod[nodex[1:4, 
        16:18, 517:nex]]), pz = as.vector(zcod[nodez[1:4, 16:18, 517:nex]]), 1)
    p7 <- rbind(px = as.vector(xcod[nodex[1:4, 19:21, 517:nex]]), py = as.vector(ycod[nodex[1:4, 
        19:21, 517:nex]]), pz = as.vector(zcod[nodez[1:4, 19:21, 517:nex]]), 1)
    p8 <- rbind(px = as.vector(xcod[nodex[1:4, 22:24, 517:nex]]), py = as.vector(ycod[nodex[1:4, 
        22:24, 517:nex]]), pz = as.vector(zcod[nodez[1:4, 22:24, 517:nex]]), 1)
    p9 <- rbind(px = as.vector(xcod[nodex[1:4, 25:27, 517:nex]]), py = as.vector(ycod[nodex[1:4, 
        25:27, 517:nex]]), pz = as.vector(zcod[nodez[1:4, 25:27, 517:nex]]), 1)
    p10 <- rbind(px = as.vector(xcod[nodex[1:4, 28:30, 517:nex]]), py = as.vector(ycod[nodex[1:4, 
        28:30, 517:nex]]), pz = as.vector(zcod[nodez[1:4, 28:30, 517:nex]]), 1)
    p11 <- rbind(px = as.vector(xcod[nodex[1:4, 31:33, 517:nex]]), py = as.vector(ycod[nodex[1:4, 
        31:33, 517:nex]]), pz = as.vector(zcod[nodez[1:4, 31:33, 517:nex]]), 1)
    p12 <- rbind(px = as.vector(xcod[nodex[1:4, 34:36, 1:nex]]), py = as.vector(ycod[nodex[1:4, 
        34:36, 1:nex]]), pz = as.vector(zcod[nodez[1:4, 34:36, 1:nex]]), 1)
    # 
    vertices1 <- as.vector(p1)
    vertices2 <- as.vector(p2)
    vertices3 <- as.vector(p3)
    vertices4 <- as.vector(p4)
    vertices5 <- as.vector(p5)
    vertices6 <- as.vector(p6)
    vertices7 <- as.vector(p7)
    vertices8 <- as.vector(p8)
    vertices9 <- as.vector(p9)
    vertices10 <- as.vector(p10)
    vertices11 <- as.vector(p11)
    vertices12 <- as.vector(p12)
    # 
    indices1 <- indices2 <- NULL
    for (i in seq(1, (nex - 516) * 4 * 3, 4)) indices1 <- c(indices1, rep(i:(i + 
        3), 3))  # 
    for (i in seq(1, nex * 4 * 3, 4)) indices2 <- c(indices2, rep(i:(i + 3), 3))  # 
    # 
    open3d(windowRect = c(75, 1020, 535, 1520))
    bg3d("white")
    model1 <- tmesh3d(vertices1, indices2)
    model2 <- tmesh3d(vertices2, indices1)
    model3 <- tmesh3d(vertices3, indices1)
    model4 <- tmesh3d(vertices4, indices1)
    model5 <- tmesh3d(vertices5, indices1)
    model6 <- tmesh3d(vertices6, indices1)
    model7 <- tmesh3d(vertices7, indices1)
    model8 <- tmesh3d(vertices8, indices1)
    model9 <- tmesh3d(vertices9, indices1)
    model10 <- tmesh3d(vertices10, indices1)
    model11 <- tmesh3d(vertices11, indices1)
    model12 <- tmesh3d(vertices12, indices2)
    # 
    shade3d(model1, color = "gray")
    shade3d(model2, color = "gray")
    shade3d(model3, color = "gray")
    shade3d(model4, color = "gray")
    shade3d(model5, color = "gray")
    shade3d(model6, color = "gray")
    shade3d(model7, color = "gray")
    shade3d(model8, color = "gray")
    shade3d(model9, color = "gray")
    shade3d(model10, color = "gray")
    shade3d(model11, color = "gray")
    shade3d(model12, color = "gray")
    # 
    wire3d(model1, color = "black", size = 2)
    wire3d(model2, color = "black", size = 2)
    wire3d(model3, color = "black", size = 2)
    wire3d(model4, color = "black", size = 2)
    wire3d(model5, color = "black", size = 2)
    wire3d(model6, color = "black", size = 2)
    wire3d(model7, color = "black", size = 2)
    wire3d(model8, color = "black", size = 2)
    wire3d(model9, color = "black", size = 2)
    wire3d(model10, color = "black", size = 2)
    wire3d(model11, color = "black", size = 2)
    wire3d(model12, color = "black", size = 2)
    # 
    points3d(x = xb[iinn[, 1]], y = yb[iinn[, 1]], z = zcod[nzin], col = "red", size = 10)
    points3d(x = xb[(iinn[, 1] - 2)%%128 + 1], y = yb[(iinn[, 1] - 2)%%128 + 1], 
        z = zcod[nzin], col = "red", size = 5)
    points3d(x = xb[(iinn[, 1] + 0)%%128 + 1], y = yb[(iinn[, 1] - 0)%%128 + 1], 
        z = zcod[nzin], col = "red", size = 5)
    points3d(x = xb[iinn[, 1]], y = yb[iinn[, 1]], z = zcod[nzin + 1], col = "red", 
        size = 5)
    points3d(x = xb[iinn[, 1]], y = yb[iinn[, 1]], z = zcod[nzin - 1], col = "red", 
        size = 5)
    points3d(x = xb[iott[, 1]], y = yb[iott[, 1]], z = zcod[nzout], col = "blue", 
        size = 10)
    # rgl.postscript('3Dmodel1.pdf',fmt='pdf',drawText=FALSE)
    rgl.snapshot("3Dmodel1.png", fmt = "png")
    play3d(spin3d(axis = c(1, 0, 1)), duration = 2)
    rgl.snapshot("3Dmodel2.png", fmt = "png")
    # rgl.postscript('3Dmodel2.pdf',fmt='pdf',drawText=FALSE)
    p1 <- rbind(px = as.vector(xcod[nodex2[1:4, 101, 10]]), py = as.vector(ycod[nodex2[1:4, 
        101, 10]]), pz = as.vector(zcod[nodez2[1:4, 101, 10]]), 1)
    p2 <- rbind(px = as.vector(xcod[nodex2[1:4, 101, 11]]), py = as.vector(ycod[nodex2[1:4, 
        101, 11]]), pz = as.vector(zcod[nodez2[1:4, 101, 11]]), 1)
    p3 <- rbind(px = as.vector(xcod[nodex2[1:4, 101, 12]]), py = as.vector(ycod[nodex2[1:4, 
        101, 12]]), pz = as.vector(zcod[nodez2[1:4, 101, 12]]), 1)
    vertices1 <- as.vector(p1)
    vertices2 <- as.vector(p2)
    vertices3 <- as.vector(p3)
    indices <- c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4)
    # open3d(windowRect=c(1555,990,1955,1490))
    open3d(windowRect = c(555, 1020, 955, 1520))
    # msgWindow(type = 'minimize',which = dev.cur())
    bg3d("white")
    model1 <- tmesh3d(vertices1, indices)
    shade3d(model1, color = gray.colors(10)[2])
    model2 <- tmesh3d(vertices2, indices)
    shade3d(model2, color = gray.colors(10)[4])
    model3 <- tmesh3d(vertices3, indices)
    shade3d(model3, color = gray.colors(10)[6])
    text3d(x = xcod[nodex2[1, 101, 10]], y = ycod[nodex2[1, 101, 10]], z = zcod[nodez2[1, 
        101, 10]], "A1")
    text3d(x = xcod[nodex2[2, 101, 10]], y = ycod[nodex2[2, 101, 10]], z = zcod[nodez2[2, 
        101, 10]], "A2")
    text3d(x = xcod[nodex2[3, 101, 10]], y = ycod[nodex2[3, 101, 10]], z = zcod[nodez2[3, 
        101, 10]], "A3")
    text3d(x = xcod[nodex2[4, 101, 10]], y = ycod[nodex2[4, 101, 10]], z = zcod[nodez2[4, 
        101, 10]], "A4")
    text3d(x = xcod[nodex2[1, 101, 11]], y = ycod[nodex2[1, 101, 11]], z = zcod[nodez2[1, 
        101, 11]], "B1", adj = c(-1.5, -1.5))
    text3d(x = xcod[nodex2[2, 101, 11]], y = ycod[nodex2[2, 101, 11]], z = zcod[nodez2[2, 
        101, 11]], "B2", adj = c(-1.5, -1.5))
    text3d(x = xcod[nodex2[3, 101, 11]], y = ycod[nodex2[3, 101, 11]], z = zcod[nodez2[3, 
        101, 11]], "B3", adj = c(-1.5, -1.5))
    text3d(x = xcod[nodex2[4, 101, 11]], y = ycod[nodex2[4, 101, 11]], z = zcod[nodez2[4, 
        101, 11]], "B4", adj = c(-1.5, -1.5))
    text3d(x = xcod[nodex2[1, 101, 12]], y = ycod[nodex2[1, 101, 12]], z = zcod[nodez2[1, 
        101, 12]], "C1", adj = c(-0.5, -0.5))
    text3d(x = xcod[nodex2[2, 101, 12]], y = ycod[nodex2[2, 101, 12]], z = zcod[nodez2[2, 
        101, 12]], "C2", adj = c(-0.5, -0.5))
    text3d(x = xcod[nodex2[3, 101, 12]], y = ycod[nodex2[3, 101, 12]], z = zcod[nodez2[3, 
        101, 12]], "C3", adj = c(-0.5, -0.5))
    text3d(x = xcod[nodex2[4, 101, 12]], y = ycod[nodex2[4, 101, 12]], z = zcod[nodez2[4, 
        101, 12]], "C4", adj = c(-0.5, -0.5))
    writeWebGL(width = 500, height = 550)
}
################################################################### 
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", 
    "yellow", "#FF7F00", "red", "#7F0000"))
################################################################### 
sigmanode <- function(x, sigma, nodex, nex) {
    elements <- sigma[nodex[1, 1, 1:nex] == x | nodex[2, 1, 1:nex] == x | nodex[3, 
        1, 1:nex] == x]
    return(log10(1000 * mean(elements)))
}
################################################################### 
writeDICOMFile <- function(dd, fname, endian = "little", skip128 = TRUE, DICM = TRUE, 
    value) {
    x <- readDICOMFile(system.file("dcm/Abdo.dcm", package = "oro.dicom"))
    x$img <- (dd)
    obj <- nrow(x$img)
    x$hdr[x$hdr$name == "Rows", c("length", "value")] <- c("2", obj)
    obj <- ncol(x$img)
    x$hdr[x$hdr$name == "Columns", c("length", "value")] <- c("2", obj)
    obj <- "16"
    x$hdr[x$hdr$name == "BitsAllocated", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- "12"
    x$hdr[x$hdr$name == "BitsStored", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    obj <- "11"
    x$hdr[x$hdr$name == "HighBit", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    x$hdr[x$hdr$name == "PixelData", "length"] <- nrow(x$img) * ncol(x$img) * extractHeader(x$hdr, 
        "BitsAllocated")/8
    obj <- "1.0 1.0"
    x$hdr[x$hdr$name == "PixelSpacing", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    obj <- value$ImageComments
    x$hdr[x$hdr$name == "ImageComments", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- value$InstanceNumber
    x$hdr[x$hdr$name == "InstanceNumber", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- paste(x$hdr[x$hdr$name == "SOPInstanceUID", "value"], ".", value$InstanceNumber, 
        sep = "")
    x$hdr[x$hdr$name == "SOPInstanceUID", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- paste(x$hdr[x$hdr$name == "SeriesInstanceUID", "value"], ".", value$SeriesTime, 
        sep = "")
    x$hdr[x$hdr$name == "SeriesInstanceUID", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- paste(x$hdr[x$hdr$name == "MediaStorageSOPInstanceUID", "value"], ".", 
        value$InstanceNumber, sep = "")
    x$hdr[x$hdr$name == "MediaStorageSOPInstanceUID", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- value$PatientsName  #'Tohru Yamaguchi'
    if (nchar(obj)%%2 == 1) 
        obj <- paste(obj, " ", sep = "")
    x$hdr[x$hdr$name == "PatientsName", c("length", "value")] <- c(nchar(obj), obj)
    obj <- value$PatientID
    x$hdr[x$hdr$name == "PatientID", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    obj <- value$PatientsSex
    x$hdr[x$hdr$name == "PatientsSex", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    obj <- value$PatientsBirthDate
    x$hdr[x$hdr$name == "PatientsBirthDate", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- value$StudyDate
    x$hdr[x$hdr$name == "StudyDate", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    obj <- value$SeriesDate
    x$hdr[x$hdr$name == "SeriesDate", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    obj <- value$SeriesTime
    x$hdr[x$hdr$name == "SeriesTime", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    obj <- value$ProtocolName
    x$hdr[x$hdr$name == "ProtocolName", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    obj <- value$PatientsWeight
    x$hdr[x$hdr$name == "PatientsWeight", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- "TDU,Chiba,Japan"
    x$hdr[x$hdr$name == "InstitutionName", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- "Kao Corporation"
    x$hdr[x$hdr$name == "Manufacturer", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    obj <- value$StudyDescription
    x$hdr[x$hdr$name == "StudyDescription", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- value$SeriesDescription
    x$hdr[x$hdr$name == "SeriesDescription", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- value$NumberOfIterations
    x$hdr[x$hdr$name == "HeartRate", c("group", "element", "name", "code", "length", 
        "value")] <- c("0018", "9739", "NumberOfIterations", "US", "2", obj)
    obj <- value$ErrorComment
    x$hdr[x$hdr$name == "ScanningSequence", c("group", "element", "name", "code", 
        "length", "value")] <- c("300A", "0402", "SetupImageComment", "ST", nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- value$PatientSize
    x$hdr[x$hdr$name == "SequenceVariant", c("group", "element", "name", "code", 
        "length", "value")] <- c("0010", "1020", "PatientsSize", "DS", nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- "SC"
    x$hdr[x$hdr$name == "Modality", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    obj <- "EITsolve Version 1.50"
    x$hdr[x$hdr$name == "SoftwareVersions", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- value$PhotometricInterpretation
    x$hdr[x$hdr$name == "PhotometricInterpretation", c("length", "value")] <- c(nchar(obj) + 
        nchar(obj)%%2, obj)
    obj <- sum(as.integer(charToRaw(value$SeriesDescription)))
    x$hdr[x$hdr$name == "SeriesNumber", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    # x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name=='ScanningSequence',1],] x$hdr <-
    # x$hdr[-row(x$hdr)[x$hdr$name=='SequenceVariant',1],]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "ScanOptions", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "MRAcquisitionType", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "SliceThickness", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "RepetitionTime", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "NumberOfAverages", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "ImagingFrequency", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "ImagedNucleus", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "EchoNumbers", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "MagneticFieldStrength", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "SpacingBetweenSlices", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "EchoTime", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "NumberOfPhaseEncodingSteps", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "EchoTraInLength", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "PercentSampling", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "PercentPhaseFieldOfView", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "LowRRValue", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "HighRRValue", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "IntervalsAcquired", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "IntervalsRejected", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "ReceiveCoilName", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "TransmitCoilName", 1], ]
    x$hdr <- x$hdr[-row(x$hdr)[x$hdr$name == "InPlanePhaseEncodingDirection", 1], 
        ]
    obj <- (nrow(x$hdr) - 1) * 2
    x$hdr[x$hdr$name == "GroupLength", c("length", "value")] <- c(nchar(obj) + nchar(obj)%%2, 
        obj)
    ################################ 
    dicomGroup <- dicom.dic$group
    dicomElement <- dicom.dic$element
    dicomName <- dicom.dic$name
    dicomCode <- dicom.dic$code
    vrCode <- dicom.VR$code
    if (endian == "swap") {
        writeGroup <- writeElement <- function(x, con, endian) {
            xx <- 256 * (as.numeric(paste("0x", substr(x, 3, 4), sep = ""))) + as.numeric(paste("0x", 
                substr(x, 1, 2), sep = ""))
            writeBin(as.integer(xx), con, size = 2L, endian = endian)
        }
    } else {
        writeGroup <- writeElement <- function(x, con, endian) {
            xx <- 256 * (as.numeric(paste("0x", substr(x, 1, 2), sep = ""))) + as.numeric(paste("0x", 
                substr(x, 3, 4), sep = ""))
            writeBin(as.integer(xx), con, size = 2L, endian = endian)
        }
    }
    fid <- file(fname, "wb")
    if (skip128) {
        seek(fid, where = 128)
    }
    if (DICM) {
        writeChar("DICM", fid, nchars = 4, eos = NULL)
    }
    y <- x$hdr[order(x$hdr[, 2]), ]
    x$hdr <- y[order(y[, 1]), ]
    for (i in 1:nrow(x$hdr)) {
        # for (i in 1:90){
        seek.old <- seek(fid)
        writeGroup(x$hdr[i, ]$group, fid, endian = endian)
        writeElement(x$hdr[i, ]$element, fid, endian = endian)
        xx <- x$hdr[i, ]$code
        writeChar(xx, fid, nchar = 2, eos = NULL)
        if (xx == "OB" | xx == "OF" | xx == "SQ" | xx == "UT" | xx == "UN") {
            writeBin(as.integer("00"), fid, size = 2L, endian = endian)
            if (x$hdr[i, ]$value == "skipped") {
                writeBin(as.integer(as.integer(x$hdr[i, ]$length)%%256), fid, size = 4L, 
                  endian = endian)
                writeBin(as.integer(256), fid, size = 2L, endian = endian)
            } else {
                writeBin(as.integer(as.integer(x$hdr[i, ]$length)%%256), fid, size = 2L, 
                  endian = endian)
                writeBin(as.integer(floor(as.integer(x$hdr[i, ]$length)/256)), fid, 
                  size = 2L, endian = endian)
                writeBin(as.integer(x$hdr[i, ]$value), fid, size = 2L, endian = endian)
            }
        } else if (xx == "UL") {
            writeBin(as.integer(x$hdr[i, ]$length), fid, size = 2L, endian = endian)
            writeBin(as.integer(x$hdr[i, ]$value), fid, size = 4L, endian = endian)
        } else if (xx == "US") {
            writeBin(as.integer(x$hdr[i, ]$length), fid, size = 2L, endian = endian)
            writeBin(as.integer(x$hdr[i, ]$value), fid, size = 2L, endian = endian)
            seek(fid)
        } else if (xx == "OW") {
            writeBin(as.integer("00"), fid, size = 2L, endian = endian)
            writeBin(as.integer(x$hdr[i, ]$length), fid, size = 4L, endian = endian)
            img <- t(x$img)
            img <- img[, ncol(img):1]
            writeBin(as.integer(img), fid, endian = endian, size = 2)
        } else {
            writeBin(as.integer(x$hdr[i, ]$length), fid, size = 2L, endian = endian)
            if (as.integer(x$hdr[i, ]$length) != 0) 
                writeBin(charToRaw(x$hdr[i, ]$value)[1:(as.integer(x$hdr[i, ]$length))], 
                  fid, endian = endian)
            seek(fid)
        }
    }
    close(fid)
}

###################################################################### 
captureImage <- function() {
    dev.set(dev.prev())
    ras <- grid.cap()
    dev.set(dev.next())
    rast <- grid.cap()
    ras <- rast[(nrow(rast)):1, 1:ncol(rast)]
    img <- matrix(255 * 3 - colSums(col2rgb(as.vector(ras))), nrow = nrow(rast), 
        ncol = ncol(rast))
    return(img)
}
################################################################################ 
findnode <- function(x, y) {
    x1 <- as.double(c(x[1], x[2]))
    x2 <- as.double(c(x[3], x[4]))
    x3 <- as.double(c(x[5], x[6]))
    y1 <- as.double(c(y[1], y[2]))
    y2 <- as.double(c(y[3], y[4]))
    y3 <- as.double(c(y[5], y[6]))
    res <- NULL
    if (identical(x1, y1)) 
        res <- rbind(x1, res)
    if (identical(x1, y2)) 
        res <- rbind(x1, res)
    if (identical(x1, y3)) 
        res <- rbind(x1, res)
    if (identical(x2, y1)) 
        res <- rbind(x2, res)
    if (identical(x2, y2)) 
        res <- rbind(x2, res)
    if (identical(x2, y3)) 
        res <- rbind(x2, res)
    if (identical(x3, y1)) 
        res <- rbind(x3, res)
    if (identical(x3, y2)) 
        res <- rbind(x3, res)
    if (identical(x3, y3)) 
        res <- rbind(x3, res)
    res <- cbind(res, sqrt(res[, 1]^2 + res[, 2]^2))
    return(res[res[, 3] == max(res[, 3]), 1:2])
}
########################################## 
plotelect <- function(ex, shift = TRUE, log = TRUE) {
    exdata2 <- matrix(NA, 32, 32)
    for (j in 0:31) {
        for (i in 0:31) {
            if (shift) 
                exdata2[j + 1, i + 1] <- ex[j + 1, (i + j + 1)%%32 + 1] else exdata2[j + 1, i + 1] <- ex[j + 1, i + 1]
        }
    }
    plot(x = 1:32, exdata2[1, ], col = 1, log = ifelse(log, "y", ""), ylim = c(ifelse(log, 
        0.001, -200), ifelse(log, 50, 200)))
    for (i in 1:32) points(x = 1:32, exdata2[i, ], col = i, type = "b")
}
plotelect64 <- function(ex, shift = TRUE) {
    exdata2 <- matrix(NA, 32, 64)
    for (j in 1:1) {
        for (i in 0:63) {
            if (shift) 
                exdata2[j, i + 1] <- ex[j + 1, (i + j + 1)%%64 + 1] else exdata2[j, i + 1] <- ex[j + 1, i + 1]
        }
    }
    plot(x = 1:64, exdata2[1, ], col = 1, ylim = c(-250, 250), type = "n")
    for (i in 1:32) points(x = 1:64, exdata2[i, ], col = i, type = "b")
}
################################################################ 
geoMean <- function(x) {
    return(prod(x)^(1/length(x)))
} 
