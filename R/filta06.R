filta06 <- function(file, sigma, bound32, nbun = 100, kfilt = 1, rfilt = 0.1, strt = 2, 
    ratmap = 0.5, col = FALSE, clean = TRUE, element = 708, nex, nez, nodex, nodez, 
    nnodex, nnodez, xcod, ycod, nb, ibnd) {
    xy <- matrix(NA, 6, nex)
    zplxy <- matrix(0, 501, 501)
    for (j in 1:nex) {
        xy[, j] <- c(xcod[nodex[1, 4, j]], ycod[nodex[1, 4, j]], xcod[nodex[2, 4, 
            j]], ycod[nodex[2, 4, j]], xcod[nodex[3, 4, j]], ycod[nodex[3, 4, j]])
    }
    # browser()
    bound32$boundry[is.na(bound32$boundry)] <- 0
    # if (element==516) { fname <-
    # paste(Sys.getenv('R_LIBS_USER'),'/EITsolve/extdata/fem5data.RData',sep='') if
    # (file.exists(fname)) load(fname) else
    # load(paste(Sys.getenv('R_HOME'),'/library/EITsolve/extdata/fem5data.RData',sep=''))
    # fem5data[is.na(fem5data)] <- 99.99 mesh <- fem5data }else{ fname <-
    # paste(Sys.getenv('R_LIBS_USER'),'/EITsolve/extdata/f9m5data.RData',sep='') if
    # (file.exists(fname)) load(fname) else
    # load(paste(Sys.getenv('R_HOME'),'/library/EITsolve/extdata/f9m5data.RData',sep=''))
    # f9m5data[is.na(f9m5data)] <- 99.99 mesh <- f9m5data } h9m5data <-
    # read.table(paste(getwd(),'/h9m5data.qqq',sep=''),col.names=paste('V',1:10,sep=''),fill=TRUE)
    # mesh <- h9m5data browser() mesh <-
    # as.matrix(read.table(paste(getwd(),'/h9m5data.qqq',sep=''),fill=TRUE,colClasses='numeric'))[1:26052,1:10]
    # #### mesh[is.na(mesh)] <- 99.99 file.copy(file,'SGMOUT0.PRN',overwrite=TRUE)
    file <- basename(file)
    print(file)
    # browser()
    status <- .Fortran("filta06", as.character(substr(file, 1, 8)), as.double(as.vector(sigma)), 
        as.integer(nex), as.integer(nez), as.vector(as.integer((nodex[1:3, 4, ]))), 
        as.vector(as.integer(nodez[1:3, 1, 708])), as.integer(nnodex), as.integer(nnodez), 
        as.double(as.vector(xcod)), as.double(as.vector(ycod)), as.double(as.vector(xy)), 
        as.integer(nb), as.integer(as.vector(ibnd)), as.double(as.vector(t(as.matrix(bound32$boundry)))), 
        as.integer(length(as.matrix(bound32$boundry))), nbun = as.integer(nbun), 
        nbuny = as.integer(nbun), as.integer(kfilt), as.double(rfilt), as.double(strt), 
        as.double(ratmap), zplxy = as.double(as.vector(zplxy)), status = as.integer(1))
    # file2 <- paste(substr(file,1,8),'.DAT',sep='') nline <- readLines(file2,n=1)
    # nline <- unlist(strsplit(nline,' ')) nline <- as.numeric(nline[length(nline)])
    # dd<- read.csv(file2,skip=8,header=FALSE,colClasses='numeric',nrow=nline) dd2 <-
    # as.matrix(dd) browser()
    dd2 <- matrix(as.double(as.vector(status$zplxy)), nrow = 501, ncol = 501)
    dd2 <- dd2[1:as.integer(status$nbun), 1:as.integer(status$nbuny)]
    dd3 <- t(dd2[, 1:(ncol(dd2) - 1)])
    # p2 <- pixmapGrey(-dd3)
    # write.pnm(p2,file=paste(substr(file,1,8),'PGM',sep='.'),forceplain=TRUE)
    dd4 <- t(dd3[nrow(dd3):1, ])
    dd5 <- dd4
    # browser()
    dd4 <- 10^dd4
    dd4 <- 1000 * dd4
    dd4 <- log10(dd4)
    dd5 <- dd4
    dd5[dd4 < log10(0.01)] <- log10(0.001)
    dd5[dd4 > log10(1)] <- log10(1)
    dd6 <- -t(dd5)
    p2 <- pixmapGrey(dd6[nrow(dd6):1, ])
    write.pnm(p2, file = paste(substr(file, 1, 8), "PGM", sep = "."), forceplain = TRUE)
    if (col) {
        # jet.colors <- colorRampPalette(c('#00007F', 'blue', '#007FFF', 'cyan',
        # '#7FFF7F', 'yellow', '#FF7F00', 'red', '#7F0000'))
        color <- jet.colors(128)
    } else color <- colorpanel(128, low = "white", high = "black")
    par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), bg = "white")
    plot(0, 0, xlim = c(-192, 192), ylim = c(-144, 144), type = "n", axes = FALSE, 
        asp = 1)
    image(x = (1:nrow(dd5)) - (nrow(dd5)/2) - 20, y = (1:ncol(dd5)) - (ncol(dd5)/2), 
        z = dd5, axes = FALSE, asp = ncol(dd5)/nrow(dd5), col = color, zlim = c(log10(0.01), 
            log10(1)), bg = "white", add = TRUE)
    makecolor <- function(x, minc = log10(0.01), maxc = log10(1)) {
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
        ytop = -130 + (seq(256, 0, by = -32)[2:9]), col = color[seq(0, 128, by = 128/8)[8:1] + 
            128/16])
    text(x = 165, y = (seq(256, 0, by = -32) - 130), label = makecolor(seq(0, 256, 
        by = 32)), adj = -0.2, cex = 1, col = "black")
    text(x = 165, y = 140, label = expression(S/m), cex = 1, col = "black")
    lines(x = c(-190, -190), y = c(-135, -125), col = "black")
    lines(x = c(-90, -90), y = c(-135, -125), col = "black")
    arrows(x0 = -170, x1 = -190, y0 = -130, y1 = -130, length = 0.1, col = "black")
    arrows(x0 = -110, x1 = -90, y0 = -130, y1 = -130, length = 0.1, col = "black")
    text(-140, -130, "100 mm", cex = 1, col = "black")
    # if (clean==TRUE) { file.remove(file2) file.remove('sgmout0.prn')
    # ile.remove('sgmout1.prn') file.remove('thetest.txt') }
    return(list(z = dd5, status = status))
} 
