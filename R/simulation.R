######################### Start point ###################################################
#' Execute numerical simulation
#' @param model type of model
#' @param interval The interval of current electrodes
#' @param resEnh Do resolution enhancement
#' @param anisotropy Consider the conductive anisotropy
#' @param est64 shape parameter is 64 or 32
#' @return normaly zero
simulationStudy <- function(model, interval, resEnh, anisotropy, est64, pattern) {
    #browser()
  #Rhpc_initialize()
  #
    filename <- paste("Model", model + 1, sep = "")  #file.choose()
    # setwd(paste(Sys.getenv('R_HOME'),'/library/EITsolve/extdata/simulation/',sep=''))
    old_dir <- setwd(paste(.libPaths()[1], "/AbdominalEIT/data/", sep = ""))
    direc <- getwd()  # dirname(filename)
    # if (substr(direc,nchar(direc)-9,nchar(direc)) != 'Simulation')
    # setwd(paste(direc,'Simulation',sep='/')) interval <- c(1:15) interval <- c(1) #
    # for test purpose anisotropy <- FALSE
    if (resEnh <= 2) 
        heterogeneity <- FALSE else heterogeneity <- TRUE
    if ((resEnh%%2) == 1) 
        resEnh <- TRUE else resEnh <- FALSE
    #parallel <- TRUE
    #bypass <- FALSE
    absolute <- FALSE
    # 
    idval <- list(PatientID = "301602")
    idval <- c(idval, list(ProtocolName = paste("100kHz,", "1.0mA")))
    s <- strsplit(getwd(), "/")[[1]]
    ss <- substr(s[length(s)], 1, 8)
    idval <- c(idval, list(PatientsName = "Tohru^F^Yamaguchi", PatientsSex = "M ", 
        PatientsBirthDate = "19630505", PatientsWeight = "62.5", PatientSize = "174.0 ", 
        StudyDate = ss, SeriesDate = format(Sys.time(), "%Y%m%d"), SeriesTime = format(Sys.time(), 
            "%H%M%S"), SeriesDescription = s[length(s)], InstanceNumber = "1", PhotometricInterpretation = "MONOCHROME1"))
    ### 
    #mode <- "parallel"
    ################################################################### make a model
    graphics.off()
    while (rgl.cur()) rgl.close()
    # rgl.quit() library(rgl)
    mxe <- 708
    mze <- 36
    mxn <- 419
    mzn <- 13
    mxb <- 128
    maxcur <- 1024
    jcur <- 32
    kvol <- 32
    resEnhRatio <- list(c(2,2,4,4,4,6,1,1,1),c(2,2,4,4,4,6,1,1,1))  
    nexr <- c(6, 24, 52, 100, 164, 260, 388, 516, 708)
    nexc <- list(c(nexr[1],(nexr[2]-nexr[1]),(nexr[3]-nexr[2]),(nexr[4]-nexr[3]),(nexr[5]-nexr[4]),(nexr[6]-nexr[5]),(nexr[7]-nexr[6]),(nexr[8]-nexr[7]),(nexr[9]-nexr[8]))/resEnhRatio[[1]],
                 c(nexr[1],(nexr[2]-nexr[1]),(nexr[3]-nexr[2]),(nexr[4]-nexr[3]),(nexr[5]-nexr[4]),(nexr[6]-nexr[5]),(nexr[7]-nexr[6]),(nexr[8]-nexr[7]),(nexr[9]-nexr[8]))/resEnhRatio[[2]])
    nex0 <- 260
    nexx <- 261
    vb <- array(0, c(mxb + 1, jcur))
    # vb0 <- array(0, c(mxb+1,jcur))
    pot <- vector("numeric", length = mxn)
    acqmode <- rep(list(ifelse(pattern==0,"parallel","bypass")), 1)
    #browser()
    exyout <- rep(list(array(0, c(2, mxe))), 1)
    exyz0 <- rep(list(array(0, c(3, mze, mxe, jcur,kvol))), 1)
    exyz1 <- rep(list(array(0, c(3, mze, mxe, jcur,kvol))), 1)
    #exyz2 <- rep(list(array(0, c(3, mze, mxe, jcur,kvol))), 1)
    #exyz3 <- rep(list(array(0, c(3, mze, mxe, jcur,kvol))), 1)
    sensematxy <- rep(list(array(0, c(mxe, kvol, jcur))), 1)
    sensematz <- rep(list(array(0, c(mxe, kvol, jcur))), 1)
    sensemat <- rep(list(array(0, c(mxe, kvol, jcur))), 1)
    pod <- array(0, c(kvol, jcur))
    podex <- array(0, c(kvol, jcur))
    sigma <- vector("numeric", length = mxe)
    sigmaxy <- vector("numeric", length = mxe)
    sigmaz <- vector("numeric", length = mxe)
    iinn <- array(0, c(32, 4))
    iott <- array(0, c(32, 4))
    vbex <- array(0, c(mxb + 1, jcur))
    difp <- rep(list(array(0, c(kvol, jcur))), 1)
    difp0 <- rep(list(array(0, c(kvol, jcur))), 1)
    bx0xy <- vector("numeric", mxe)  # for anisotropy
    bx0z <- vector("numeric", mxe)
    numpol <- array(0, c(kvol, jcur))
    numpol2 <- array(0, c(kvol, jcur))
    mxbnd <- 520
    mxnode <- 5447
    podmax <- vector("numeric", 16)
    # ----------------START--------------------------------------- aa <-
    # rep(list(matrix(0,mxnode,mxbnd)),1) fname <-
    # paste(Sys.getenv('R_LIBS_USER'),'/EITsolve/extdata/h9m5data.RData',sep='') if
    # (file.exists(fname)) load(fname) else
    # load(paste(Sys.getenv('R_HOME'),'/library/EITsolve/extdata/h9m5data.RData',sep=''))
    # h9m5data[is.na(h9m5data)] <- 99.99 h9m5 <- h9m5data
    filenm <- paste("img_____")
    filenmd <- paste(filenm, ".txt", sep = "")
    substr(filenmd, 7, 8) <- "dv"
    cat("a9l5qxfa \n")
    cat("Generating mesh for h9l5qxfa... \n")
    aniso <- anisotropy
    anisovis <- 1
    if (model == 1) {
        notch <- 0.85
        flat <- 0.75
        vis <- 5e-04
        mus <- 5e-04
        sub <- 5e-05
        ra <- 100
        rb <- 100
    } else if (model == 2) {
        notch <- 0.75
        flat <- 0.85
        vis <- 5e-04
        mus <- 5e-04
        sub <- 5e-05
        ra <- 100
        rb <- 100
    } else if (model == 3) {
        notch <- 0.9
        flat <- 0.7
        vis <- 5e-04
        mus <- 5e-04
        sub <- 5e-05
        ra <- 100
        rb <- 100
    } else if (model == 4) {
        notch <- 0.7
        flat <- 0.9
        vis <- 5e-04
        mus <- 5e-04
        sub <- 5e-05
        ra <- 100
        rb <- 100
    } else if (model == 5) {
        notch <- 0.85
        flat <- 0.75
        vis <- 1e-04
        mus <- 5e-04
        sub <- 5e-05
        ra <- 140
        rb <- 100
    } else if (model == 6) {
        notch <- 0.75
        flat <- 0.85
        vis <- 1e-04
        mus <- 5e-04
        sub <- 5e-05
        ra <- 140
        rb <- 100
    } else if (model == 7) {
        notch <- 0.9
        flat <- 0.7
        vis <- 5e-05
        mus <- 5e-04
        sub <- 5e-05
        ra <- 140
        rb <- 100
    } else if (model == 8) {
        notch <- 0.7
        flat <- 0.9
        vis <- 5e-05
        mus <- 5e-04
        sub <- 5e-05
        ra <- 140
        rb <- 100
    } else if (model == 9) {
        notch <- 0.9
        flat <- 0.7
        aniso <- 10
        anisovis <- 10
        vis <- 2e-04
        mus <- 2e-04
        sub <- 5e-05
        ra <- 100
        rb <- 100
    } else if (model == 10) {
        notch <- 0.7
        flat <- 0.9
        vis <- 1e-04
        mus <- 2e-04
        sub <- 5e-05
        ra <- 140
        rb <- 100
        aniso <- 10
    } else if (model == 11 | model == 14 | model == 17) {
        # notch <- 0.80 flat <- 0.9
        vis <- 5e-05
        mus <- 5e-04
        sub <- 5e-05
        ra <- 130
        rb <- 80
        aniso <- 1
    } else if (model == 12 | model == 15 | model == 18) {
        # notch <- 0.75 flat <- 0.9
        vis <- 5e-05
        mus <- 5e-04
        sub <- 5e-05
        ra <- 140
        rb <- 100
        aniso <- 1
    } else if (model == 13 | model == 16 | model == 19) {
        # notch <- 0.7 flat <- 0.9
        vis <- 5e-05
        mus <- 5e-04
        sub <- 5e-05
        ra <- 150
        rb <- 120
        aniso <- 1
    } else {
        stop()
    }
    # rad64 <- rep(flat,64) rad64[33] <- notch rad64[c(32,34)] <- (notch+flat)/2
    rad64 <- rep(NA, 64)
    if (model == 17 | model == 18 | model == 19) {
        rad64[1] <- 0.82
        rad64[2] <- rad64[64] <- 0.86
        rad64[3] <- rad64[63] <- 0.85
        rad64[4] <- rad64[62] <- 0.86
        rad64[5] <- rad64[61] <- 0.81
        rad64[6] <- rad64[60] <- 0.77
        rad64[7] <- rad64[59] <- 0.72
        rad64[8] <- rad64[58] <- 0.695
        rad64[9] <- rad64[57] <- 0.67
        rad64[10] <- rad64[56] <- 0.695
        rad64[11] <- rad64[55] <- 0.72
        rad64[12] <- rad64[54] <- 0.725
        rad64[13] <- rad64[53] <- 0.745
        rad64[14] <- rad64[52] <- 0.795
        rad64[15] <- rad64[51] <- 0.82
        rad64[16] <- rad64[50] <- 0.84
        rad64[33] <- 0.68
        rad64[32] <- rad64[34] <- 0.72
        rad64[31] <- rad64[35] <- 0.72
        rad64[30] <- rad64[36] <- 0.73
        rad64[29] <- rad64[37] <- 0.745
        rad64[28] <- rad64[38] <- 0.745
        rad64[27] <- rad64[39] <- 0.745
        rad64[26] <- rad64[40] <- 0.745
        rad64[25] <- rad64[41] <- 0.745
        rad64[24] <- rad64[42] <- 0.72
        rad64[23] <- rad64[43] <- 0.72
        rad64[22] <- rad64[44] <- 0.745
        rad64[21] <- rad64[45] <- 0.795
        rad64[20] <- rad64[46] <- 0.85
        rad64[19] <- rad64[47] <- 0.87
        rad64[18] <- rad64[48] <- 0.87
        rad64[17] <- rad64[49] <- 0.86
    } else if (model == 14 | model == 15 | model == 16) {
        rad64[1] <- 0.85
        rad64[2] <- rad64[64] <- 0.9
        rad64[3] <- rad64[63] <- 0.875
        rad64[4] <- rad64[62] <- 0.85
        rad64[5] <- rad64[61] <- 0.825
        rad64[6] <- rad64[60] <- 0.8
        rad64[7] <- rad64[59] <- 0.75
        rad64[8] <- rad64[58] <- 0.725
        rad64[9] <- rad64[57] <- 0.7
        rad64[10] <- rad64[56] <- 0.725
        rad64[11] <- rad64[55] <- 0.75
        rad64[12] <- rad64[54] <- 0.75
        rad64[13] <- rad64[53] <- 0.775
        rad64[14] <- rad64[52] <- 0.825
        rad64[15] <- rad64[51] <- 0.85
        rad64[16] <- rad64[50] <- 0.875
        rad64[33] <- 0.7
        rad64[32] <- rad64[34] <- 0.75
        rad64[31] <- rad64[35] <- 0.75
        rad64[30] <- rad64[36] <- 0.75
        rad64[29] <- rad64[37] <- 0.775
        rad64[28] <- rad64[38] <- 0.775
        rad64[27] <- rad64[39] <- 0.775
        rad64[26] <- rad64[40] <- 0.775
        rad64[25] <- rad64[41] <- 0.775
        rad64[24] <- rad64[42] <- 0.75
        rad64[23] <- rad64[43] <- 0.75
        rad64[22] <- rad64[44] <- 0.775
        rad64[21] <- rad64[45] <- 0.825
        rad64[20] <- rad64[46] <- 0.9
        rad64[19] <- rad64[47] <- 0.9
        rad64[18] <- rad64[48] <- 0.9
        rad64[17] <- rad64[49] <- 0.9
    } else if (model == 11 | model == 12 | model == 13) {
        rad64[1] <- 0.86
        rad64[2] <- rad64[64] <- 0.91
        rad64[3] <- rad64[63] <- 0.895
        rad64[4] <- rad64[62] <- 0.88
        rad64[5] <- rad64[61] <- 0.855
        rad64[6] <- rad64[60] <- 0.83
        rad64[7] <- rad64[59] <- 0.78
        rad64[8] <- rad64[58] <- 0.755
        rad64[9] <- rad64[57] <- 0.73
        rad64[10] <- rad64[56] <- 0.755
        rad64[11] <- rad64[55] <- 0.78
        rad64[12] <- rad64[54] <- 0.805
        rad64[13] <- rad64[53] <- 0.83
        rad64[14] <- rad64[52] <- 0.855
        rad64[15] <- rad64[51] <- 0.88
        rad64[16] <- rad64[50] <- 0.905
        rad64[33] <- 0.8
        rad64[32] <- rad64[34] <- 0.83
        rad64[31] <- rad64[35] <- 0.835
        rad64[30] <- rad64[36] <- 0.84
        rad64[29] <- rad64[37] <- 0.845
        rad64[28] <- rad64[38] <- 0.84
        rad64[27] <- rad64[39] <- 0.835
        rad64[26] <- rad64[40] <- 0.83
        rad64[25] <- rad64[41] <- 0.825
        rad64[24] <- rad64[42] <- 0.83
        rad64[23] <- rad64[43] <- 0.85
        rad64[22] <- rad64[44] <- 0.9
        rad64[21] <- rad64[45] <- 0.915
        rad64[20] <- rad64[46] <- 0.93
        rad64[19] <- rad64[47] <- 0.93
        rad64[18] <- rad64[48] <- 0.93
        rad64[17] <- rad64[49] <- 0.93
    }
    ry <- -rb * cos(seq(0, 2 * pi, length.out = 65)[1:64])
    rx <- ra * sin(seq(0, 2 * pi, length.out = 65)[1:64])
    geodata <- sqrt(rx^2 + ry^2)[seq(2, 64, 2)]
    # write.table(file='abc_3.csv',sep='\n',data.frame(matrix(geodata,1,32)),row.names=FALSE,col.names=FALSE)
    bound <- bound32R(geodata)
    bndd <- t(bound$boundry)
    #browser()
    zcod <- NULL
    nzin <- NULL
    nzot <- NULL
    zrange <- 250
    layer <- c(0.152, 0.265, 0.371, 0.505, 0.73)  # c(0.13,0.24,0.34,0.49,0.55)
    ret <- meshgen(rad = rad64, zrange = zrange, mxbnd = mxbnd, mxnode = mxnode, 
        bndd = bndd, mxn = mxn, mzn = mzn, mxb = mxb, mze = mze, mxe = mxe, zcod, 
        nzin, nzot, layer = layer)  #c(0.14,0.24,0.34,0.48))  #c(0.152,0.265,0.371,0.505))
    theta <- ret$theta
    st <- ret$st
    stxy <- ret$stxy
    stz <- ret$stz
    cex <- ret$cex
    vol <- ret$vol
    xav <- ret$xav
    yav <- ret$yav
    arar <- ret$arar
    bnd <- ret$bnd
    ibnd <- ret$ibnd
    zcod <- ret$zcod
    nb <- ret$nb
    nbw <- ret$nbw
    nodex <- ret$nodex
    nodez <- ret$nodez
    nnodex <- ret$nnodex
    nnodez <- ret$nnodez
    nnodezx <- ret$nnodezx
    nex <- ret$nex
    nez <- ret$nez
    nzin <- ret$nzin
    nzot <- ret$nzot
    xcod <- ret$xcod
    ycod <- ret$ycod
    zcod <- ret$zcod
    xb <- ret$xb
    yb <- ret$yb
    cat("... finished\n")
    eleint <- interval
    # browser()
    
    sigmaxy[1:(nex0 - 96)] <- vis
    sigmaz[1:(nex0 - 96)] <- vis * anisovis  #0.0005
    sigma[1:(nex0 - 96)] <- sqrt((2 * (sigmaxy[1:(nex0 - 96)])^2 + (sigmaz[1:(nex0 - 
        96)])^2)/3)
    # 
    sigmaxy[(nex0 - 96 + 1):nex0] <- mus
    sigmaz[(nex0 - 96 + 1):nex0] <- mus * aniso  #0.0005
    sigma[(nex0 - 96 + 1):nex0] <- sqrt((2 * (sigmaxy[(nex0 - 96 + 1):nex0])^2 + 
        (sigmaz[(nex0 - 96 + 1):nex0])^2)/3)
    # 
    muscle <- 1:(nex0 - 96)
    if (model == 11 | model == 12 | model == 13) {
        # sigma[(nex0-96+15):(nex0-96+16)] <- sigmaz[(nex0-96+15):(nex0-96+16)] <-
        # sigmaxy[(nex0-96+15):(nex0-96+16)] <- vis sigma[(nex0-96+33):(nex0-96+34)] <-
        # sigmaz[(nex0-96+33):(nex0-96+34)] <- sigmaxy[(nex0-96+33):(nex0-96+34)] <- vis
        # sigma[(nex0-96+63):(nex0-96+64)] <- sigmaz[(nex0-96+63):(nex0-96+64)] <-
        # sigmaxy[(nex0-96+63):(nex0-96+64)] <- vis sigma[(nex0-96+81):(nex0-96+82)] <-
        # sigmaz[(nex0-96+81):(nex0-96+82)] <- sigmaxy[(nex0-96+81):(nex0-96+82)] <- vis
        muscle <- NULL
        muscle <- c(muscle, c(1:3, 4:6))  #
        muscle <- c(muscle, c((6 + 1):(6 + 6), (6 + 13):(6 + 18)))  #
        muscle <- c(muscle, c((24 + 1):(24 + 8), (24 + 21):(24 + 28)))  #
        muscle <- c(muscle, c((52 + 1):(52 + 14), (52 + 35):(52 + 48)))
        #muscle <- c(muscle, c((52 + 1):(52 + 18), (52 + 31):(52 + 48)))
        muscle <- c(muscle, c((100 + 1):(100 + 10), (100 + 55):(100 + 64),(100 + 27):(100 + 38)))
        sigmaxy[muscle] <- sigmaz[muscle] <- mus
        sigma[muscle] <- sqrt((2 * (sigmaxy[muscle])^2 + (sigmaz[muscle])^2)/3)
    } else if (model == 14 | model == 15 | model == 16) {
        # sigma[(nex0-96+15):(nex0-96+16)] <- sigmaz[(nex0-96+15):(nex0-96+16)] <-
        # sigmaxy[(nex0-96+15):(nex0-96+16)] <- vis sigma[(nex0-96+33):(nex0-96+34)] <-
        # sigmaz[(nex0-96+33):(nex0-96+34)] <- sigmaxy[(nex0-96+33):(nex0-96+34)] <- vis
        # sigma[(nex0-96+63):(nex0-96+64)] <- sigmaz[(nex0-96+63):(nex0-96+64)] <-
        # sigmaxy[(nex0-96+63):(nex0-96+64)] <- vis sigma[(nex0-96+81):(nex0-96+82)] <-
        # sigmaz[(nex0-96+81):(nex0-96+82)] <- sigmaxy[(nex0-96+81):(nex0-96+82)] <- vis
        muscle <- NULL
        muscle <- c(muscle, c(1:3, 4:6))  #
        muscle <- c(muscle, c((6 + 1):(6 + 6), (6 + 13):(6 + 18)))
        muscle <- c(muscle, c((6 + 7):(6 + 12)))  #
        muscle <- c(muscle, c((24 + 1):(24 + 8), (24 + 21):(24 + 28)))  #
        muscle <- c(muscle, c((52 + 1):(52 + 14), (52 + 35):(52 + 48)))
        #muscle <- c(muscle, c((52 + 1):(52 + 18), (52 + 31):(52 + 48)))
        muscle <- c(muscle, c((100 + 1):(100 + 10), (100 + 55):(100 + 64),(100 + 27):(100 + 38)))
        sigmaxy[muscle] <- sigmaz[muscle] <- mus
        sigma[muscle] <- sqrt((2 * (sigmaxy[muscle])^2 + (sigmaz[muscle])^2)/3)
    } else if (model == 17 | model == 18 | model == 19) {
        # sigma[(nex0-96+15):(nex0-96+16)] <- sigmaz[(nex0-96+15):(nex0-96+16)] <-
        # sigmaxy[(nex0-96+15):(nex0-96+16)] <- vis sigma[(nex0-96+33):(nex0-96+34)] <-
        # sigmaz[(nex0-96+33):(nex0-96+34)] <- sigmaxy[(nex0-96+33):(nex0-96+34)] <- vis
        # sigma[(nex0-96+63):(nex0-96+64)] <- sigmaz[(nex0-96+63):(nex0-96+64)] <-
        # sigmaxy[(nex0-96+63):(nex0-96+64)] <- vis sigma[(nex0-96+81):(nex0-96+82)] <-
        # sigmaz[(nex0-96+81):(nex0-96+82)] <- sigmaxy[(nex0-96+81):(nex0-96+82)] <- vis
        muscle <- NULL
        muscle <- c(muscle, c(1:3, 4:6))  #
        muscle <- c(muscle, c((6 + 1):(6 + 6), (6 + 13):(6 + 18)))
        muscle <- c(muscle, c((6 + 7):(6 + 12)))  ##
        muscle <- c(muscle, c((24 + 1):(24 + 8), (24 + 21):(24 + 28)))
        muscle <- c(muscle, c((24 + 9):(24 + 20)))  #
        #muscle <- c(muscle, c((24 + 9):(24 + 20)))  #
        muscle <- c(muscle, c((52 + 1):(52 + 14), (52 + 35):(52 + 48)))
        #muscle <- c(muscle, c((52 + 1):(52 + 18), (52 + 31):(52 + 48)))
        muscle <- c(muscle, c((100 + 1):(100 + 10), (100 + 55):(100 + 64),(100 + 27):(100 + 38)))
        sigmaxy[muscle] <- sigmaz[muscle] <- mus
        sigma[muscle] <- sqrt((2 * (sigmaxy[muscle])^2 + (sigmaz[muscle])^2)/3)
    }
    cavityarea <- sum(arar[1:nex0])
    visfatarea <- cavityarea - sum(arar[(nex0 - 96 + 1):nex0]) - sum(arar[muscle])
    # 
    sigma[nexx:nex] <- sigmaxy[nexx:nex] <- sigmaz[nexx:nex] <- sub  #0.00005
    # for (multi in 1:1) { # for multiple cat(eleint,'multi:',multi,'\n')
    # iinn[1:32,1] <- (4*(1:32)-2)%%128+1 #3 iott[1:32,1] <-
    # (4*((1:32)+eleint)-2)%%128+1 #3 iinn[1:32,2] <- (4*(1:32)-2)%%128+1 #3
    # iott[1:32,2] <- (4*(1:32)+2)%%128+1 #1 iinn[1:32,1] <- (4*(1:32)-2)%%128+1 #3
    # iott[1:32,1] <- (4*((1:32)+eleint)-2)%%128+1 #3 iinn[1:32,2] <-
    # (4*(1:32)-2)%%128+1 #3 iott[1:32,2] <- (4*(1:32)+2)%%128+1 #1
    #iinn[1:32, 1] <- (4 * (1:32) - 6)%%128 + 1
    #iott[1:32, 1] <- (4 * ((1:32) + eleint) - 6)%%128 + 1
    iinn[1:32, 1] <- (4 * (1:32) - 2)%%128 + 1
    iott[1:32, 1] <- (4 * (1:32) + 2)%%128 + 1
    iinn[1:32, 2] <- (4 * (1:32) - 2)%%128 + 1
    iott[1:32, 2] <- (4 * (1:32) + 2)%%128 + 1
    iinn[1:32, 3] <- (4 * (1:32) - 2)%%128 + 1
    iott[1:32, 3] <- (4 * (1:32) + 2)%%128 + 1
    iinn[1:32, 4] <- (4 * (1:32) - 2)%%128 + 1
    iott[1:32, 4] <- (4 * (1:32) + 2)%%128 + 1
    #
    numb <- 0
    for (j in 1:32) {
        for (k in 1:32) {
            numb <- numb + 1
            numpol[k, j] <- numb
        }
        for (k in 1:32) {
            numpol2[k, j] <- numpol[k, j]
            kp1 <- k%%32 + 1
            km1 <- (k + 30)%%32 + 1
            if (numpol[kp1, j] == 0 | numpol[km1, j] == 0) 
                numpol2[k, j] <- 0
        }
    }
    # 
    if (numb != maxcur) {
        cat("maxcur inconsitent \n")
        iret <- 1
        stop()
        return(iret)
    }
    nbound <- nb/2
    # } Forward Problem #########################################
    bb <- vector(mode = "numeric", length = mxnode)
    oneb2 <- 0.5
    oneb4 <- 0.25
    oneb8 <- 0.125
    oneb16 <- 1/16
    oneb32 <- 1/32
    oneb64 <- 1/64
    #aa0 <- matrix(0, mxnode, mxbnd)
    nvecout <- 1
    ########################################################################### 
    cat("-JUN1-\n")
    #aa0[1:mxnode, 1:mxbnd] <- 0
    aa0 <- setupaaPAfort(sigmaxy = sigmaxy, sigmaz = sigmaz, nodex = nodex, 
        nodez = nodez, st = st, stxy = stxy, stz = stz, nnodez = nnodez, nbw = nbw, 
        nez = nez, nex = nex, mxbnd=mxbnd,mxnode=mxnode,mze=mze,mxe=mxe)
    cat("-JUN2-\n")
    aa0 <- systemafort(numnp = nnodezx, mband = nbw, a = aa0, mxbnd = mxbnd, mxnode = mxnode)
    cat("-JUN3-\n")
    ########### Current poles ################################################
    mxcurp <- jcur
    mxvolp <- kvol
    nzcur <- nzin
    ############ loop ########################################################
    nmulti <- 1
    # 
    rm(ret)
    #browser()
    for (multi in 1:nmulti) {
        cat(eleint[[multi]], "multi:", multi, "acqmode; ", acqmode[[multi]], "\n")
        #ret <- junPA(npat = 1, iter = 0, aa0 = aa0, mxbnd = mxbnd, mxnode = mxnode, 
        #    maxcur = maxcur, numpol = numpol, sigma = sigma, sigmaxy = sigmaxy, sigmaz = sigmaz, 
        #    iinn = iinn, iott = iott, pod = pod, podex = podex, vb = vb, difp = difp[[multi]], 
        #    difp0 = difp0[[multi]], nnodex = nnodex, nodex = nodex, nodez = nodez, 
        #    nnodez = nnodez, nnodezx = nnodezx, st = st, stxy = stxy, stz = stz, 
        #    nbw = nbw, nb = nb, ibnd = ibnd, bnd = bnd, vol = vol, cex = cex, 
        #    nex = nex, nez = nez, nzin = nzin, nzot = nzot, jcur = jcur, kvol = kvol, 
        #    mxe = mxe, mze = mze, mxb = mxb, acqmode=acqmode[[multi]],useCluster=FALSE,nCPU=nCPU)
        ret <- junPA(npat = 1, iter = 0, mxbnd = mxbnd, mxnode = mxnode, 
                     maxcur = maxcur, numpol = numpol, sigma = sigma, sigmaxy = sigmaxy, sigmaz = sigmaz, 
                     iinn = iinn, iott = iott, pod = pod, podex = podex, vb = vb, difp = difp[[multi]], 
                     difp0 = difp0[[multi]], nnodex = nnodex, nodex = nodex, nodez = nodez, 
                     nnodez = nnodez, nnodezx = nnodezx, st = st, stxy = stxy, stz = stz, 
                     nbw = nbw, nb = nb, ibnd = ibnd, bnd = bnd, vol = vol, cex = cex, 
                     nex = nex, nez = nez, nzin = nzin, nzot = nzot, jcur = jcur, kvol = kvol, 
                     mxe = mxe, mze = mze, mxb = mxb, acqmode=acqmode[[multi]],useCluster=FALSE,nCPU=nCPU)
        # aa[[multi]] <- ret$aa
        difp[[multi]] <- ret$difp
        difp0[[multi]] <- ret$difp0
        # difpraw[[multi]] <- ret$difpraw
        exyz0[[multi]] <- ret$exyz0
        exyz1[[multi]] <- ret$exyz1
        #exyz2[[multi]] <- ret$exyz2
        #exyz3[[multi]] <- ret$exyz3
        # exyout[[multi]] <- ret$exyout
        pod <- ret$pod
        pot <- ret$pot
        vb <- ret$vb
        vb0 <- ret$vb0
    }
    dev <- sum(sapply(difp, sumup))/nmulti
    # devraw <- sum(sapply(difpraw,sumup))/nmulti
    cat("dev,sigmaxy[1]", dev, sigmaxy[1], "\n")
    # 
    pb <- txtProgressBar(min = 1, max = 4, style = 3)
    # 
    #browser()
    #
    for (i in 1:4) {  # 1:4
        
        cat("--- Computing the sensitivity of conductivity ---\n")
        for (multi in 1:nmulti) {
            cat(eleint[[multi]], "multi:", multi, "acqmode; ", acqmode[[multi]], "\n")
            ret <- senseCalc(mxnode = mxnode, numpol = numpol, sigma = sigma, sigmaxy = sigmaxy, 
                sigmaz = sigmaz, pod = pod, jcur = jcur, kvol = kvol, mxe = mxe, 
                mze = mze, nex = nex, nez = nez, vol = vol, exyz0 = exyz0[[multi]], 
                exyz1 = exyz1[[multi]])
            sensemat[[multi]] <- ret$sensemat
            sensematxy[[multi]] <- ret$sensematxy
            sensematz[[multi]] <- ret$sensematz
        }
        cat("--- Finished ---\n")
        cat("Radius of layers,", layer, "\n")
        layer <- adjLayer(layer, sensemat[[1]], nexr,nexc=nexc[[1]], rad64, arar, resEnhRatio=resEnhRatio[[1]])
        cat("Radius of layers,", layer, "\n")
        # 
        cat("---Generating a mesh--- \n")
        ret <- meshgen(rad = rad64, zrange = zrange, mxbnd = mxbnd, mxnode = mxnode, 
            bndd = bndd, mxn = mxn, mzn = mzn, mxb = mxb, mze = mze, mxe = mxe, zcod, 
            nzin, nzot, zcalc = FALSE, layer = layer)
        theta <- ret$theta
        # rad32 <- ret$rad
        xcod <- ret$xcod
        ycod <- ret$ycod
        xcod2 <- ret$xcod2
        ycod2 <- ret$ycod2
        # zcod <- ret$zcod
        arar <- ret$arar
        xb <- ret$xb
        yb <- ret$yb
        bnd <- ret$bnd
        xav <- ret$xav
        yav <- ret$yav
        st <- ret$st
        stxy <- ret$stxy
        stz <- ret$stz
        cex <- ret$cex
        vol <- ret$vol
        nez <- ret$nez
        cat("--- FINISHED, meshgen.---\n")
        for (multi in 1:nmulti) {
            cat(eleint[[multi]], "multi:", multi, "acqmode; ", acqmode[[multi]], "\n")
            #ret <- junPA(npat = 1, iter = 0, aa0 = aa0, mxbnd = mxbnd, mxnode = mxnode, 
            #    maxcur = maxcur, numpol = numpol, sigma = sigma, sigmaxy = sigmaxy, 
            #    sigmaz = sigmaz, iinn = iinn, iott = iott, pod = pod, podex = podex, 
            #    vb = vb, difp = difp[[multi]], difp0 = difp0[[multi]], nnodex = nnodex, 
            #    nodex = nodex, nodez = nodez, nnodez = nnodez, nnodezx = nnodezx, 
            #    st = st, stxy = stxy, stz = stz, nbw = nbw, nb = nb, ibnd = ibnd, 
            #    bnd = bnd, vol = vol, cex = cex, nex = nex, nez = nez, 
            #    nzin = nzin, nzot = nzot, jcur = jcur, kvol = kvol, mxe = mxe, mze = mze, 
            #    mxb = mxb, acqmode=acqmode[[multi]],useCluster=FALSE,nCPU=nCPU)
            ret <- junPA(npat = 1, iter = 0, mxbnd = mxbnd, mxnode = mxnode, 
                         maxcur = maxcur, numpol = numpol, sigma = sigma, sigmaxy = sigmaxy, 
                         sigmaz = sigmaz, iinn = iinn, iott = iott, pod = pod, podex = podex, 
                         vb = vb, difp = difp[[multi]], difp0 = difp0[[multi]], nnodex = nnodex, 
                         nodex = nodex, nodez = nodez, nnodez = nnodez, nnodezx = nnodezx, 
                         st = st, stxy = stxy, stz = stz, nbw = nbw, nb = nb, ibnd = ibnd, 
                         bnd = bnd, vol = vol, cex = cex, nex = nex, nez = nez, 
                         nzin = nzin, nzot = nzot, jcur = jcur, kvol = kvol, mxe = mxe, mze = mze, 
                         mxb = mxb, acqmode=acqmode[[multi]],useCluster=FALSE,nCPU=nCPU)
            # aa[[multi]] <- ret$aa
            difp[[multi]] <- ret$difp
            difp0[[multi]] <- ret$difp0
            exyz0[[multi]] <- ret$exyz0
            exyz1[[multi]] <- ret$exyz1
            #exyz2[[multi]] <- ret$exyz2
            #exyz3[[multi]] <- ret$exyz3
            pod <- ret$pod
            pot <- ret$pot
            vb <- ret$vb
            vb0 <- ret$vb0
        }
        dev <- sum(sapply(difp, sumup))/nmulti
        # devraw <- sum(sapply(difpraw,sumup))/nmulti
        cat("dev,sigmaxy[1]", dev, sigmaxy[1], "\n")
        # 
        setTxtProgressBar(pb, i)
        cat("\n")
    }
    cat("--- Computing the sensitivity of conductivity again ---\n")
    #browser()
    for (multi in 1:nmulti) {
        cat(eleint[[multi]], "multi:", multi, "acqmode; ", acqmode[[multi]], "\n")
        ret <- senseCalc(mxnode = mxnode, numpol = numpol, sigma = sigma, sigmaxy = sigmaxy, 
            sigmaz = sigmaz, pod = pod, jcur = jcur, kvol = kvol, mxe = mxe, mze = mze, 
            nex = nex, nez = nez, vol = vol, exyz0 = exyz0[[multi]], exyz1 = exyz1[[multi]])
        sensemat[[multi]] <- ret$sensemat
        sensematxy[[multi]] <- ret$sensematxy
        sensematz[[multi]] <- ret$sensematz
    }
    cat("--- Finished ---\n")
    layer <- adjLayer(layer, sensemat[[1]], nexr,nexc=nexc[[1]], rad64, arar, resEnhRatio=resEnhRatio[[1]])
    cat("Radius of layers,", layer, "\n")
    # 
    cat("---Generating a mesh--- \n")
    ret <- meshgen(rad = rad64, zrange = zrange, mxbnd = mxbnd, mxnode = mxnode, 
        bndd = bndd, mxn = mxn, mzn = mzn, mxb = mxb, mze = mze, mxe = mxe, zcod, 
        nzin, nzot, zcalc = FALSE, layer = layer)
    theta <- ret$theta
    xcod <- ret$xcod
    ycod <- ret$ycod
    xcod2 <- ret$xcod2
    ycod2 <- ret$ycod2
    arar <- ret$arar
    xb <- ret$xb
    yb <- ret$yb
    bnd <- ret$bnd
    xav <- ret$xav
    yav <- ret$yav
    st <- ret$st
    stxy <- ret$stxy
    stz <- ret$stz
    cex <- ret$cex
    vol <- ret$vol
    nez <- ret$nez
    cat("--- FINISHED, meshgen.---\n")
    
    ########### Loop #############################################################
    for (ncur in 1:jcur) {
        iin <- iinn[ncur, 1]
        iot <- iott[ncur, 1]
        cat("iin,iot", iin, iot, " ")
        ############ BB ##########################################################
        bb <- vector(mode = "numeric", length = nnodezx)  ###
        bb[1:nnodezx] <- 0
        # iinm1 <- (ibnd[(iin-2)%%nb+1])*nnodez + nzcur iinc <- (ibnd[iin ])*nnodez +
        # nzcur iinp1 <- (ibnd[(iin)%%nb+1 ])*nnodez + nzcur iotp1 <- (ibnd[(iot)%%nb+1
        # ])*nnodez + nzcur iotc <- (ibnd[iot ])*nnodez + nzcur iotm1 <-
        # (ibnd[(iot-2)%%nb+1])*nnodez + nzcur
        
        iinm1 <- (ibnd[(iin - 2)%%nb + 1] - 1) * nnodez + nzcur
        iinc <- (ibnd[iin] - 1) * nnodez + nzcur
        iinp1 <- (ibnd[(iin)%%nb + 1] - 1) * nnodez + nzcur
        iotp1 <- (ibnd[(iot)%%nb + 1] - 1) * nnodez + nzcur
        iotc <- (ibnd[iot] - 1) * nnodez + nzcur
        iotm1 <- (ibnd[(iot - 2)%%nb + 1] - 1) * nnodez + nzcur
        bb[iinc] <- oneb2
        bb[iotc] <- -oneb2
        if (iin%%2 == 1) {
            bb[iinc + 1] <- oneb8
            bb[iinc - 1] <- oneb8
            bb[iinm1] <- oneb4 * bnd[iin]/(bnd[iin] + bnd[iin + 1])
            bb[iinp1] <- oneb4 * bnd[iin + 1]/(bnd[iin] + bnd[iin + 1])
        } else {
            stop
        }
        if (iot%%2 == 1) {
            bb[iotc + 1] <- -oneb8
            bb[iotc - 1] <- -oneb8
            bb[iotm1] <- -oneb4 * bnd[iot]/(bnd[iot] + bnd[iot + 1])
            bb[iotp1] <- -oneb4 * bnd[iot + 1]/(bnd[iot] + bnd[iot + 1])
        } else {
            stop
        }
        # zz <- file('tstboundR.txt','w') for (j in 1:nb)
        # cat(j,bb[(j-1)*nnodez+1],bb[(j-1)*nnodez+2],bb[j],'\n',file=zz) close(zz)
        cat("-JUN4-\n")
        bb <- systembfort(numnp = nnodezx, mband = nbw, b = bb, a = aa0, mxbnd = mxbnd, 
            mxnode = mxnode)
        cat("-JUN5-\n")
        summ <- mean(bb[1:nnodezx])
        bb[1:nnodezx] <- bb[1:nnodezx] - summ
        ############## BOUNDARY POTENTIAL############
        summ <- 0
        for (j in 1:nb) {
            # vb0[1:nb,ncur] <- vb[1:nb,ncur] depending on acquisitin mode; START
            if (acqmode[[1]] == "parallel") {
                vb[j, ncur] <- bb[nnodez * (ibnd[j] - 1) + nzot]
                # vb[j,ncur] <- bb[nnodez*(ibnd[j])+nzot]
            } else if ((ibnd[j]) == (ibnd[iin]) | (ibnd[j]) == (ibnd[iot])) {
                vb[j, ncur] <- bb[nnodez * (ibnd[j] - 1) + nzot]
            } else if ((ibnd[j]) == (ibnd[iin] - 1) | (ibnd[j]) == (ibnd[iin] + 1) | 
                (ibnd[j]) == (ibnd[iot] - 1) | (ibnd[j]) == (ibnd[iot] + 1)) {
                vb[j, ncur] <- bb[nnodez * (ibnd[j] - 1) + nzot]
            } else if ((ibnd[j]) == (ibnd[iin] - 2) | (ibnd[j]) == (ibnd[iin] + 2) | 
                (ibnd[j]) == (ibnd[iot] - 2) | (ibnd[j]) == (ibnd[iot] + 2)) {
                vb[j, ncur] <- bb[nnodez * (ibnd[j] - 1) + (nzot + nzin)/2]
            } else vb[j, ncur] <- bb[nnodez * (ibnd[j] - 1) + nzin]
        }
        # browser() difp0s <- difp[,ncur]
        summ <- mean(vb[1:nb, ncur])
        # vb0[nb+1,ncur] <- vb0[1,ncur]
        vb[nb + 1, ncur] <- vb[1, ncur]
        vb[1:(nb + 1), ncur] <- vb[1:(nb + 1), ncur] - summ
        for (j in 1:kvol) {
            pod[j, ncur] <- (vb[iinn[j, 2], ncur] - vb[iott[j, 2], ncur])  ######### modified!!!
        }
        ######################################################################## exyouts <- exyz[1:2,3*(ifelse(mode=='parallel',nzot,nzin))+1,1:nex]
        for (n in 1:nnodex) pot[n] <- bb[nnodez * (n - 1) + ifelse((acqmode[[1]] != "parallel" & 
            iin != (n) & iot != (n)), nzin, nzot)]
        # return(list(difps=difps,difp0s=difp0s,pods=pods,exyz0s=exyz0s,exyouts=exyouts,pots=pots,vbs=vbs,vb0s=vb0s))
    }
    
    ############################ 
    bound <- list(boundry = bound)
    intd <- list(eleint)
    # geodata <- rep(100,32) ##### for test
    write.table(file = "abc_3.csv", sep = "\n", geodata, row.names = FALSE, col.names = FALSE)  # for test
    bound <- bound32R(geodata)
    bound$crossarea
    bound$circumference
    ################################ 
    file.remove("abc_3.csv")
    ldata <- NULL
    # pod <- pod[c(32,1:31),]
    exdata <- list(datan32 = t(pod))  ######!!!
    # 
    ldata <- NULL
    ldata <- c(ldata, list(list(data = exdata$datan32, interval = intd)))
    # browser()
    idval <- c(idval, list(StudyDescription = paste("Electrical Impedance Tomography", 
        ifelse(aniso > 1, "considering anisotropy", "not considering anisotropy")), 
        ImageComments = paste("Interval of current electrodes:", paste(paste(intd, 
            acqmode[[1]]), collapse = ", "))))
    # dev.new(width=4,height=6) browser()
    graphics.off()
    setupdev(anisotropy = ifelse(aniso == 1, FALSE, TRUE))
    dev.set(5)
    ret <- condimgP(sigma = sigmaxy, xcod = xcod, ycod = ycod, zcod = zcod, nodex = nodex, 
        iter2 = 0, nex = nex, nexr = nexr, arar = arar, annotate = TRUE, xb = xb, 
        yb = yb)
    dev.set(6)
    ret <- condimgP(sigma = sigmaz, xcod = xcod, ycod = ycod, zcod = zcod, nodex = nodex, 
        iter2 = 0, nex = nex, nexr = nexr, arar = arar, annotate = TRUE, xb = xb, 
        yb = yb)
    dev.set(11)
    ret <- condimgP(sigma = sigma, xcod = xcod, ycod = ycod, zcod = zcod, nodex = nodex, 
        iter2 = 0, nex = nex, nexr = nexr, arar = arar, mesh = FALSE, annotate = FALSE, 
        xb = xb, yb = yb)
    # 
    refimage <- captureImage()
    if (TRUE) {
        text(-120, 143, paste("Cavity area,", round(cavityarea/100, 1)), cex = 0.7)
        text(-60, 146, expression(cm^2), cex = 0.7)
        text(55, 143, paste("Visceral fat area,", round(visfatarea/100, 1)), cex = 0.7)
        text(128, 146, expression(cm^2), cex = 0.7)
    }
    boundry <- bound$boundry[seq(2, 64, 2), ]
    # boundr <- sqrt(boundry[,'x']^2+boundry[,'y']^2) for (i in 1:32){ theta <-
    # 2*pi/32*i-2*pi/64
    # text(x=(boundr[i]+12)*sin(theta)-20,y=-(boundr[i]+10)*cos(theta),labels=i) }
    # browser()
    dev.copy2eps(file = paste("meshModel", model, ".eps", sep = ""))
    # dev.off()
    dev.set(2)
    par(mfrow = c(3, 1), mar = c(5, 3, 3, 1))
    plot(x = 1:32, y = ldata[[1]]$data[1, 1:32], type = "b", pch = 20, xlab = "Electrode", 
        ylab = "Voltage (mV)", sub = paste("(a) interval =", ldata[[1]]$interval))
    # plot(x=1:32,y=ldata[[2]]$data[5,1:32],type='b',pch=20,xlab='Electrode',ylab='Voltage
    # (mV)',sub=paste('(b) interval =',ldata[[2]]$interval))
    # dev.copy2eps(file=paste(substr(filename,1,nchar(filename)-4),'rawdata','.eps',sep=''))
    # embedFonts(file=paste(substr(filename,1,nchar(filename)-4),'rawdata','.eps',sep=''))
    # extending to 64 electrode system applying simulation noise
    ldata[[1]]$data <- ldata[[1]]$data + matrix(rnorm(1024, mean = -0, sd = 0), 32, 
        32)
    ############################################################################ 
    expnew <- preprocess2(realdata = ldata, bound = bound, abs = absolute, ang = c(9, 
        18), fact = 32)  #9,18
    ############################################################################ 
    expf <- NULL
    for (i in 1:length(intd)) {
        temp <- list(expout = expnew$cumdata[[i]]$expout, ldata = expnew$potdata[[i]]$data)
        expf <- c(expf, list(temp))
        plotelect64(expnew$cumdata[[i]]$expout, shift = FALSE)
    }
    ########################################################################################## Solveing the inverse problemexpf of the electrical impedance tomography
    #browser()
    rm(list = c("ret","sensemat", "sensematxy", "sensematz", "exyz0", "exyz1", "exyout", 
        "aa0", "cex", "st", "stxy", "stz","difp","difp0","bx0xy","bx0z","sigma","sigmaxy","sigmaz"))
    gc();    gc()
    dev.set(2)
    par(mfrow = c(2, 1), oma = c(0, 0, 0, 0), mar = c(5, 4, 1, 1))
    options(device = "x11")
    system.time(ret <- a9l5qxfa(bound32 = bound, expout = expf, maxiter = 150, 
        acqmode = acqmode, eleint = intd, anisotropy = anisotropy, redalam = 0.75, alam0 = 900, 
        sigini = rnorm(708, 5e-05, 1e-14), alamlim = 0.1, zrange = zrange, idval = idval, 
        heterogeneity = heterogeneity, est64 = TRUE, resEnh = resEnh,
        Marquardt = FALSE))
    setwd(old_dir)
    return(0)
} 
