##################################### preProcess function
preProcess <- function(object, anisotropy) {
    filename <- file.choose()
    dir <- dirname(filename)
    dir
    filename <- basename(filename)
    filename
    setwd(dir)
    # 
    interval <- c(1:15)
    # interval <- c(2) # for test purpose anisotropy <<- TRUE
    parallel <- TRUE
    bypass <- FALSE
    absolute <- FALSE
    # data(exdata)
    geodata <- as.vector(read.csv(filename, header = FALSE, skip = 9, nrows = 1))[1:32]
    # write.table(file='abc_3.csv',sep='\n',geodata,row.names=FALSE,col.names=FALSE)
    # #???????
    iddata <- read.csv(filename, header = FALSE, skip = 0, nrows = 8, blank.lines.skip = FALSE)[1:8, 
        1:2]
    idval <<- list(PatientID = as.character(iddata[1, 2]))
    idval <<- c(idval, list(ProtocolName = paste(as.character(iddata[2, 2]), "kHz,", 
        as.character(iddata[3, 2]), "mA")))
    s <- strsplit(getwd(), "/")[[1]]
    ss <- substr(s[length(s)], 1, 8)
    idval <<- c(idval, list(PatientsName = "Tohru^F^Yamaguchi", PatientsSex = "M ", 
        PatientsBirthDate = "19630505", PatientsWeight = "62.5", PatientSize = "174.0 ", 
        StudyDate = ss, SeriesDate = format(Sys.time(), "%Y%m%d"), SeriesTime = format(Sys.time(), 
            "%H%M%S"), SeriesDescription = s[length(s)], InstanceNumber = "1", PhotometricInterpretation = "MONOCHROME1"))
    dd <- read.csv(filename, header = FALSE, skip = 10)
    write.table(file = "abc_2.csv", sep = ",", dd, row.names = FALSE, col.names = FALSE)
    ### 
    bound <<- bound32R(kei = geodata)
    bound$crossarea
    bound$circumference
    ################################ 
    for (fact in 32:40) {
        ldata <- NULL
        intd <<- NULL
        amode <<- NULL
        for (int in interval) {
            if (parallel) {
                exdata <- dtconv("abc_2.csv", eleint = int, mode = "parallel")
                exdata$datan32[, 1:32] <- -exdata$datan32[c(32, 1:31), c(1:32)]  #### depending on system
                len <- (length(exdata$datan32[exdata$datan32 != 0]))
                if (len >= 1024) {
                  amode <<- c(amode, "parallel")
                  intd <<- c(intd, int)
                  ldata <- c(ldata, list(list(data = exdata$datan32, interval = int)))
                }
            }
            if (ncol(dd) == 128 & bypass) {
                exdata <- dtconv("abc_2.csv", eleint = int, mode = "bypass")
                exdata$datan32[, 1:32] <- -exdata$datan32[, c(2:32, 1)]  #### depending on system
                len <- (length(exdata$datan32[exdata$datan32 != 0]))
                if (len >= 1024) {
                  amode <<- c(amode, "bypass")
                  intd <<- c(intd, int)
                  ldata <- c(ldata, list(list(data = exdata$datan32, interval = int)))
                }
            }
        }
        # file.remove('abc_2.csv') browser()
        idval <<- c(idval, list(StudyDescription = paste("Electrical Impedance Tomography", 
            ifelse(anisotropy, "considering anisotropy", "not considering anisotropy")), 
            ImageComments = paste("Interval of current electrodes:", paste(paste(intd, 
                amode), collapse = ", "))))
        # dev.new(width=4,height=6)
        graphics.off()
        # rgl.quit()
        dev.new(widreth = 4, height = 6)
        par(mfrow = c(3, 1), mar = c(4, 3, 3, 1))
        plot(x = 1:32, y = ldata[[1]]$data[1, 1:32], type = "b", pch = 20, xlab = "Electrode", 
            ylab = "Voltage (mV)", sub = paste("(a) interval =", ldata[[1]]$interval))
        # plot(x=1:32,y=ldata[[2]]$data[5,1:32],type='b',pch=20,xlab='Electrode',ylab='Voltage
        # (mV)',sub=paste('(b) interval =',ldata[[2]]$interval)) data extraction
        dtconv <- function(filename, eleint = 2, mode) {
            table <- read.table(file = filename, header = FALSE, sep = ",")
            if (mode == "parallel") 
                table <- table[, 1:64] else table <- table[, 65:128]
            esvec <- NULL
            for (i in 0:31) {
                esvec <- c(esvec, i * 32 + (i + eleint)%%32 + 1)
            }
            estable <- table[esvec, ]
            return(list(status = 0, datan32 = estable))
        }
        ############################################################################ extending to 64 electrode system
        trimcut <- function(x) {
            x <- x[x != max(x)]
            # x <- x[x != max(x)] x <- x[x != max(x)]
            return(x)
        }
        expnew <- preprocess2(realdata = ldata, bound = bound, abs = absolute, ang = c(9, 
            18), fact = fact)  #9,18
        testdata <- expnew$potdata[[1]]$data
        tres <- (max(apply(testdata, 2, FUN = trimcut)) < 0)
        cat("division factor", fact, ifelse(tres, "OK", "NG"), "\n")
        if (tres) 
            break
        # if (TRUE) break
    }
    ############################################################################ 
    expf <<- NULL
    for (i in 1:length(intd)) {
        temp <- list(expout = expnew$cumdata[[i]]$expout * 2, ldata = expnew$potdata[[i]]$data * 
            2)  #???????
        # temp
        # <-list(expout=expnew$cumdata[[i]]$expout*1,ldata=expnew$potdata[[i]]$data*1)
        # #???????
        expf <<- c(expf, list(temp))
        plotelect64(expnew$cumdata[[i]]$expout, shift = FALSE)
    }
}
########################################################################### Solveing the inverse problem of the electrical impedance tomography
imageReconstruction <- function(obj, resEnh, anisotropy) {
    # rgl.quit()
    graphics.off()
    setupdev(anisotropy = ifelse(anisotropy == 1, FALSE, TRUE))
    filename <<- basename(getwd())
    if (resEnh <= 2) 
        heterogeneity <- FALSE else heterogeneity <- TRUE
    if ((resEnh%%2) == 1) 
        resEnh <- TRUE else resEnh <- FALSE
    system.time(ret <- a9l5qxfa(bound32 = bound, expout = expf, maxiter = 250, 
        mode = amode, eleint = intd, anisotropy = anisotropy, redalam = 0.8, alam0 = 900, 
        sigini = rnorm(708, 8e-05, 1e-14), alamlim = 0.2, zrange = 250, idval = idval, 
        heterogeneity = heterogeneity, est64 = TRUE, resEnh = resEnh, 
        Marquardt = FALSE))
}
########################################################################### Post process
postProcess <- function(obj) {
    filename <- basename(getwd())
    anisotropy <- TRUE
    # Saving the images
    scr <- 2
    dev.set(scr)
    dev.copy2eps(file = paste(substr(filename, 1, nchar(filename) - 4), "status", 
        ".eps", sep = ""))
    # embedFonts(paste(substr(filename,1,nchar(filename)-4),'status','.eps',sep=''))
    scr <- scr + 1
    dev.set(scr)
    dev.copy2eps(file = paste(substr(filename, 1, nchar(filename) - 4), "average", 
        ".eps", sep = ""))
    # embedFonts(file=paste(substr(filename,1,nchar(filename)-4),'average','.eps',sep=''))
    scr <- scr + 1
    dev.set(scr)
    dev.copy2eps(file = paste(substr(filename, 1, nchar(filename) - 4), "smooth", 
        ".eps", sep = ""))
    # embedFonts(file=paste(substr(filename,1,nchar(filename)-4),'smooth','.eps',sep=''))
    if (anisotropy == TRUE) {
        scr <- scr + 1
        dev.set(scr)
        dev.copy2eps(file = paste(substr(filename, 1, nchar(filename) - 4), "trans", 
            ".eps", sep = ""))
        # embedFonts
        # (file=paste(substr(filename,1,nchar(filename)-4),'trans','.eps',sep=''))
        scr <- scr + 1
        dev.set(scr)
        dev.copy2eps(file = paste(substr(filename, 1, nchar(filename) - 4), "long", 
            ".eps", sep = ""))
        # embedFonts(file=paste(substr(filename,1,nchar(filename)-4),'long','.eps',sep=''))
        scr <- scr + 1
        dev.set(scr)
        dev.copy2eps(file = paste(substr(filename, 1, nchar(filename) - 4), "aniso", 
            ".eps", sep = ""))
        # embedFonts(file=paste(substr(filename,1,nchar(filename)-4),'aniso','.eps',sep=''))
        dev.set(11)
        dev.copy2eps(file = paste(substr(filename, 1, nchar(filename) - 4), "mesh", 
            ".eps", sep = ""))
        # embedFonts(file=paste(substr(filename,1,nchar(filename)-4),'mesh','.eps',sep=''))
    }
    scr <- scr + 1
    dev.set(scr)
    dev.copy2eps(file = paste(substr(filename, 1, nchar(filename) - 4), "smoothcol", 
        ".eps", sep = ""))
    # embedFonts(file=paste(substr(filename,1,nchar(filename)-4),'smoothcol','.eps',sep=''))
    if (FALSE) {
        scr <- scr + 1
        dev.set(scr)
        dev.copy2eps(file = paste(substr(filename, 1, nchar(filename) - 4), "sense", 
            ".eps", sep = ""))
        # embedFonts(file=paste(substr(filename,1,nchar(filename)-4),'sense','.eps',sep=''))
    }
    rgl.snapshot(paste(substr(filename, 1, nchar(filename) - 4), "3d", ".png", sep = ""), 
        fmt = "png", top = TRUE)
}

############################################  
