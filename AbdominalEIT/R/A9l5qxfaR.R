#############################################################################
# Three-dimensional inverse problem
# This program was originally written in FORTRAN77 by Kazuo Maki.
# A translated version was written by Tohru F. Yamaguchi.
#############################################################################
#############################################################################
#############################################################################y

# Main routine a9l5qxfa
# ' @useDynLib AbdominalEITextra
#' @useDynLib AbdominalEIT
#############################################################################
a9l5qxfa <- function(bound32,expout,zrange,sigini,maxiter,redalam,alam0,acqmode,eleint,alamlim,idval,anisotropy,heterogeneity,est64,resEnh,alpha=1,Marquardt=TRUE)
{
  while (rgl.cur()) rgl.close()
  if(FALSE) rgl.init()
  #
  #browser()
  snow::setDefaultClusterOptions(type="MPI",homogeneous=FALSE)
  #
  useGPU <- FALSE
  useCluster <- TRUE
  useRhpc <- FALSE
  nCPU <- 17 #32 #Rmpi::mpi.universe.size()-1 # 35
  #
  if (useGPU){
    library(gpuR)
  }
  #
  if (useCluster & useRhpc){
    requireNamespace("Rhpc") 
    Rhpc::Rhpc_initialize()
    cl2 <<- Rhpc::Rhpc_getHandle(nCPU)
    #cl3 <<- Rhpc::Rhpc_getHandle(2)
    # Rhpc set to options
    opstr=list("Rhpc.mpi.rank","Rhpc.mpi.procs","Rhpc.mpi.c.comm","Rhpc.mpi.f.comm")
    do.call("options",opstr)
    Rhpc::Rhpc_worker_call(cl2, "do.call","options", opstr)
    # warning! : pointer not export, worker Rhpc.mpi.c.comm is (nil) on master.
  }#
  mxe <- 708
  mze <- 36
  mxn <- 419
  mzn <- 13
  mxb <- 128
  maxcur <- 1024
  jcur <- 32
  kvol <- 32
  #  resEnhRatio <- list(c(2,3,4,3,4,6,1,1,1)
  resEnhRatio <- list(c(2,2,4,4,4,6,1,1,1),c(2,3,4,3,4,6,1,1,1))
  resEnhRatio2 <- list(c(1,1,2,2,2,3,1,1,1),c(1,1,1,1,2,3,1,1,1))
  doResEnh <- c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE)
  #resEnhRatio[[2]]<- ifelse(doResEnh,resEnhRatio[[2]],resEnhRatio[[1]])
  #
  if(resEnh){
    perm <- list(prod(ifelse(doResEnh,resEnhRatio[[1]]/resEnhRatio2[[1]],1)),prod(ifelse(doResEnh,resEnhRatio[[2]]/resEnhRatio2[[2]],1)))
  } else{
    perm <- list(1,1)
  }
  #
  nexr <- c(6,24,52,100,164,260,388,516,708)
  nexl <- c(nexr[1],diff(nexr))
  nexc <- list(c(nexr[1],(nexr[2]-nexr[1]),(nexr[3]-nexr[2]),(nexr[4]-nexr[3]),(nexr[5]-nexr[4]),(nexr[6]-nexr[5]),(nexr[7]-nexr[6]),(nexr[8]-nexr[7]),(nexr[9]-nexr[8]))/resEnhRatio[[1]],
               c(nexr[1],(nexr[2]-nexr[1]),(nexr[3]-nexr[2]),(nexr[4]-nexr[3]),(nexr[5]-nexr[4]),(nexr[6]-nexr[5]),(nexr[7]-nexr[6]),(nexr[8]-nexr[7]),(nexr[9]-nexr[8]))/resEnhRatio[[2]])
  nexp <- list(nexc[[1]][1:6]*c(1,1,1,1/2,1/2,1/2),nexc[[2]][1:6]*c(1,1,1,1/2,1/2,1/2))
  nexx <- list()
  for (j in 1:2){
    nexxtemp <- NULL
    for(i in 1:9) {
      nexxtemp <- c(nexxtemp,sum(nexc[[j]][1:i]))
    }
    nexx <- c(nexx,list(nexxtemp))
  }
  cat("Unknown parameters",nexp[[1]],"\n")
  cat("Unknown parameters",nexp[[2]],"\n")
  nexe <- list(ifelse(heterogeneity,nexx[[1]][9],nexx[[1]][6]+1),ifelse(heterogeneity,nexx[[2]][9],nexx[[2]][6]+1) )
  #
  #
  #
  #
  rndele <- list()
  for (j in 1:2){
  temp <- list(
    first=ifelse(rep(doResEnh[1],perm[[j]]),(1:perm[[j]]-1)%%resEnhRatio[[j]][1],rep(0,perm[[j]])),
    second=ifelse(rep(doResEnh[2],perm[[j]]),rep((1:resEnhRatio[[j]][2]-1),each=prod(ifelse(doResEnh[1],  resEnhRatio[[j]][1]/resEnhRatio2[[1]],  1))),rep(0,perm[[j]])),
    third =ifelse(rep(doResEnh[3],perm[[j]]),rep((seq(1,resEnhRatio[[j]][3],resEnhRatio2[[j]][3])-1),each=prod(ifelse(doResEnh[1:2],resEnhRatio[[j]][1:2]/resEnhRatio2[[2]][1:2],1))),rep(0,perm[[j]])),
    fourth=ifelse(rep(doResEnh[4],perm[[j]]),rep((seq(1,resEnhRatio[[j]][4],resEnhRatio2[[j]][4])-1),each=prod(ifelse(doResEnh[1:3],resEnhRatio[[j]][1:3]/resEnhRatio2[[2]][1:3],1))),rep(0,perm[[j]])),
    fifth =ifelse(rep(doResEnh[5],perm[[j]]),rep((seq(1,resEnhRatio[[j]][5],resEnhRatio2[[j]][5])-1),each=prod(ifelse(doResEnh[1:4],resEnhRatio[[j]][1:4]/resEnhRatio2[[2]][1:4],1))),rep(0,perm[[j]])),
    sixth =ifelse(rep(doResEnh[6],perm[[j]]),rep((seq(1,resEnhRatio[[j]][6],resEnhRatio2[[j]][6])-1),each=prod(ifelse(doResEnh[1:5],resEnhRatio[[j]][1:5]/resEnhRatio2[[2]][1:5],1))),rep(0,perm[[j]])))
  rndele <- c(rndele,list(temp))
  #
  }
  #
  #browser()
  #
  vb <- rep(list(array(0, c(mxb+1,jcur))),length(expout))
  vb0 <- rep(list(array(0, c(mxb+1,jcur))),length(expout))
  pot <- rep(list(vector("numeric", length=mxn)),length(expout))
  exyz0 <- rep(list(array(0, c(3,mze,mxe,jcur,kvol))),length(expout))
  exyz1 <- rep(list(array(0, c(3,mze,mxe,jcur,kvol))),length(expout))
  sensemat <- rep(list(array(0, c(mxe,kvol,jcur))),length(expout))
  sensematxy <- rep(list(array(0, c(mxe,kvol,jcur))),length(expout))
  sensematz <- rep(list(array(0, c(mxe,kvol,jcur))),length(expout))
  pod <- rep(list(array(0, c(kvol,jcur))),length(expout))
  podex <- rep(list(array(0, c(kvol,jcur))),length(expout))
  sigma <- sigmaxy <- sigmaz <- sigini  #vector("numeric", length=mxe)
  iinn <- rep(list(array(0, c(32,4))),length(expout))
  iott <- rep(list(array(0, c(32,4))),length(expout))
  vbex <- array(0, c(mxb+1,jcur))
  difp <- rep(list(array(0, c(kvol,jcur))),length(expout))
  difp0 <- rep(list(array(0, c(kvol,jcur))),length(expout))
  bx0xy <- vector("numeric",mxe) # for anisotropy
  bx0z <-  vector("numeric",mxe)
  numpol <- array(0, c(kvol,jcur))
  numpol2 <- array(0, c(kvol,jcur))
  mxbnd <- 520
  mxnode <- 5447
  podmax <- vector("numeric",16)
  layer <- c(0.155,0.295,0.4,0.55,0.6)
  # ----------------START---------------------------------------
  #aa <- rep(list(matrix(0,mxnode,mxbnd)),length(expout))
  bndd <- as.vector(t(as.matrix(bound32$boundry)))
  filenm <- paste("img_____")
  filenmd <- paste(filenm,".txt",sep="")
  substr(filenmd,7,8) <- "dv"
  cat("a9l5qxfa \n")
  ani.record(reset = TRUE)
  cat("zrange,maxiter \n")
  cat(zrange,maxiter,"\n")
  cat("redalam,alam0 \n")
  cat(redalam,alam0,"\n")
  cat("Generating mesh for h9l5qxfa...")
  rad64 <- rep(0.85,64)
  #browser()
  ret <- meshgen(rad=rad64,zrange=zrange,mxbnd=mxbnd,mxnode=mxnode,bndd=bndd,mxn=mxn,mzn=mzn,mxb=mxb,mze=mze,mxe=mxe,zcod,nzin,nzot,layer=layer)
  #theta <- ret$theta
  st <- ret$st
  stxy <- ret$stxy
  stz <- ret$stz
  cex <- ret$cex
  vol <- ret$vol
  #xav <- ret$xav
  #yav <- ret$yav
  arar <- ret$arar
  bnd <- ret$bnd
  ibnd <-  ret$ibnd
  nb <- ret$nb
  nodex <- ret$nodex
  nodez <- ret$nodez
  nnodex <- ret$nnodex
  nnodez <- ret$nnodez
  nnodezx <- ret$nnodezx
  nex <- ret$nex
  nez <- ret$nez
  xb <- ret$xb
  yb <- ret$yb
  xcod <- ret$xcod
  ycod <- ret$ycod
  xcod2 <- ret$xcod2
  ycod2 <- ret$ycod2
  zcod <- ret$zcod
  nzin <- ret$nzin
  nzot <- ret$nzot
  nbw <- ret$nbw
  #
  cat("... finished\n")
  #rm(list=c("ret"))
  gc();    gc()
  #browser()
  for (multi in 1:length(expout)) { # for multiple
    cat(eleint[[multi]],"multi:",multi,"\n")
    #iinn[[multi]][1:32,1] <- (4*(1:32)-6)%%128+1 #3
    #iott[[multi]][1:32,1] <- (4*((1:32)+eleint[[multi]])-6)%%128+1 #3
    iinn[[multi]][1:32,1] <- (4*(1:32)-2)%%128+1 #3
    iott[[multi]][1:32,1] <- (4*(1:32)+2)%%128+1 #1
    iinn[[multi]][1:32,2] <- (4*(1:32)-2)%%128+1 #3
    iott[[multi]][1:32,2] <- (4*(1:32)+2)%%128+1 #1
    iinn[[multi]][1:32,3] <- (4*(1:32)-2)%%128+1 #3
    iott[[multi]][1:32,3] <- (4*(1:32)+2)%%128+1 #1
    iinn[[multi]][1:32,4] <- (4*(1:32)-2)%%128+1 #3
    iott[[multi]][1:32,4] <- (4*(1:32)+2)%%128+1 #1
    #
    numb <- 0
    for (j in 1:32)  {
      for (k in 1:32)	{
        numb <- numb+1
        numpol[k,j] <- numb
      }
      for (k in 1:32)	{
        numpol2[k,j] <- numpol[k,j]
        kp1 <- k%%32+1
        km1 <- (k+30)%%32+1
        if (numpol[kp1,j] ==0 | numpol[km1,j]==0)	numpol2[k,j] <- 0
      }
    }
    #browser()
    #if (iclean != TRUE)	{
    #  zz <- file("numchkR.txt","w")
    #  cat (jcur,kvol,"\n",file=zz)
    #  for (j in 1:jcur)	cat(numpol[1:kvol,j],"\n",file=zz)
    #  for (j in 1:jcur)	cat(numpol2[1:kvol,j],"\n",file=zz)
    #  close(zz)
    #}
    if(numb != maxcur)	{
      cat("maxcur inconsitent \n")
      iret <- 1
      stop()
      return (iret)
    }
    nbound <- nb/2
    #
    #expout1 <- expout[[multi]]$expout[2:33,1:65] ##
    expout1 <- expout[[multi]]$expout[1:32,1:65] ##
    ldata1  <- expout[[multi]]$ldata[1:32,1:32]
    #
    for (j in 1:jcur){
      vbex[2*(1:(nbound))-1,j] <- expout1[j,c((nbound-1):nbound,1:(nbound-2))]
      #vbex[2*(1:nbound),j] <- (vbex[2*(1:nbound)-1,j] + vbex[2*(1:nbound)+1,j])/2
      #vbex[,j] <- vbex[c((mxb-3):mxb,1:(mxb-3)),j]
    }                                # cat("EXPOUT VBEX 2 \n")
    #if (iclean != TRUE) zz <- file("podexR.txt","w")
    #if (iclean != TRUE) cat(jcur,kvol,"\n",file=zz)
    #podmax[multi] <- 0
    for (ncur in 1:jcur){
      for (j in 1:kvol) {
        #podex[[multi]][j,ncur] <- vbex[iinn[[multi]][j,2],ncur] - vbex[iott[[multi]][j,2],ncur]
        podex[[multi]][j,ncur] <- ldata1[ncur,j] ####
      }
      #if (iclean != TRUE) cat(podex[[multi]][1:kvol,ncur],"\n",file=zz) 
    }
    #if (iclean != TRUE) close(zz)
  }
  rm(vbex)
  #sigma[1:nex] <- sigini
  #sigmaxy[1:nex] <- sigini
  #sigmaz[1:nex] <- sigini
  cat("Uniform conductvity is assumed.\n")
  if(FALSE) show3dmodel(nodex=nodex,nodez=nodez,xcod=xcod,ycod=ycod,zcod=zcod,nex=ret$nex,nez=nez,xb=xb,yb=yb,iinn=iinn[[1]],iott=iott[[1]],nzin=nzin,nzout=nzot)
    bxxy <- vector("numeric", length=mxe) # for anisotropy
    bxz  <- vector("numeric", length=mxe)
    bxtxy <- vector("numeric",mxe) #
    bxtz  <- vector("numeric",mxe)
    p <- vector("numeric",length=2*mxe) # for anisotorpy
    sigma0 <- sigma
    sigma0xy <- sigmaxy
    sigma0z <- sigmaz
    dd <- vector("numeric",length=nexx[[1]][6]+2)
    rsense <- rep(list(matrix(0,maxcur,32)), length(difp))
    ch <- c(as.character(1:9),LETTERS[1:21])
    alam <- alam0
    dev <- 0
    dev0 <- 0 ####### modified ###########
    nrest <- 0
    nquit <- 0
    nquit0 <- 0
    nquit00 <- 0
    iter <- 0
    nmulti <- length(difp)
    eps <- 0 #.0001 #0.03/nmulti
    #initial <- c(0.05,  0.10,    0, 0.5, 0.10,    0, 0.5) #skewNormal
    #initial <- c(0.05, 0.1,0.5,0.1,0.5) #Normal
    initial <- c(1, 3, 0.5) # Beta
    estres=NULL
##########################################################
    lmm <- 5
    dev <- sum(sapply(difp,sumup))/nmulti
                                        #devraw <- sum(sapply(difpraw,sumup))/nmulti
    facl <- 1
    cat("facl,dev,sigmaxy[1]",facl,dev,sigmaxy[1],"\n")
    dev1 <- dev
    cat("dev1",dev1,"\n")
    #if (nreasig == 0) {
#######################################################
        #browser()
        rm(ret)
        facl <- 10
        sigmaxy[1:nex] <- sigma0xy[1:nex]*facl
        sigmaz[1:nex] <- sigma0z[1:nex]*facl
        sigma[1:nex] <- sqrt((2*sigmaxy^2+sigmaz^2)/3)
        for (multi in 1:nmulti) {
            cat(eleint[[multi]],"multi:",multi,"acqmode; ",acqmode[[multi]],"\n")
            #retJun <- junPA(npat=1,iter=0,aa0=aa[[multi]],mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster=useCluster,nCPU=nCPU,useRhpc=useRhpc)
            retJun <- junPA(npat=1,iter=0,mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster=useCluster,nCPU=nCPU,useRhpc=useRhpc)
            #aa[[multi]] <- retJun$aa
            difp[[multi]] <- retJun$difp
            difp0[[multi]] <- retJun$difp0
            #difpraw[[multi]] <- ret$difpraw
            exyz0[[multi]] <- retJun$exyz0
            exyz1[[multi]] <- retJun$exyz1
            #exyout[[multi]] <- ret$exyout
            pod[[multi]] <- retJun$pod
            pot[[multi]] <- retJun$pot
            vb[[multi]] <- retJun$vb
            vb0[[multi]] <- retJun$vb0
      }
      rm(retJun) 
      dev <- sum(sapply(difp,sumup))/nmulti
                                        #devraw <- sum(sapply(difpraw,sumup))/nmulti
        cat("facl,dev,sigmaxy[1]",facl,dev,sigmaxy[1],"\n")
        dev2 <- dev
####################################################################################
        cat("dev1,dev2",dev1,dev2,"\n")
        facl <- 0.1
        #browser()
        sigmaxy[1:nex] <- sigma0xy[1:nex]*facl
        sigmaz[1:nex] <- sigma0z[1:nex]*facl
        sigma[1:nex] <- sqrt((2*sigmaxy^2+sigmaz^2)/3)
        for (multi in 1:nmulti) {
            cat(eleint[[multi]],"multi:",multi,"acqmode; ",acqmode[[multi]],"\n")
            #retJun <- junPA(npat=1,iter=0,aa0=aa[[multi]],mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster = useCluster,nCPU=nCPU,useRhpc=useRhpc)
            retJun <- junPA(npat=1,iter=0,mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster = useCluster,nCPU=nCPU,useRhpc=useRhpc)
            #aa[[multi]] <- retJun$aa
            difp[[multi]] <- retJun$difp
            difp0[[multi]] <- retJun$difp0
            #difpraw[[multi]] <- ret$difpraw
            exyz0[[multi]] <- retJun$exyz0
            exyz1[[multi]] <- retJun$exyz1
            #exyout[[multi]] <- ret$exyout
            pod[[multi]] <- retJun$pod
            pot[[multi]] <- retJun$pot
            vb[[multi]] <- retJun$vb
            vb0[[multi]] <- retJun$vb0
        }
        rm(retJun)
        dev <- sum(sapply(difp,sumup))/nmulti
        dev0 <- dev
###################################################################################
        cat("dev1,dev2,dev0",dev1,dev2,dev0,"\n")
        if((dev0+dev2) <= (2*dev1)) {
            cat("RESET APPROPRIATE INITIAL SIGMA \n")
            iret <- 3
            stop("iret",iret,"\n")
            return(iret)
        } else {
            facl <- (dev0-dev2)/(dev2+dev0-2*dev1)/2
            facl <- exp(facl*log(10))/2
        }
        sigmaxy[1:nex] <- sigma0xy[1:nex]*facl
        sigmaz[1:nex] <- sigma0z[1:nex]*facl
        sigma[1:nex] <- sqrt((2*sigmaxy^2+sigmaz^2)/3)
        cat("--- Solving forward problem ---\n")
        #browser()
        for (multi in 1:nmulti) {
            cat(eleint[[multi]],"multi:",multi,"acqmode; ",acqmode[[multi]],"\n")
            #retJun <- junPA(npat=1,iter=0,aa0=aa[[multi]],mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster = useCluster,nCPU=nCPU,useRhpc=useRhpc)
            retJun <- junPA(npat=1,iter=0,mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster = useCluster,nCPU=nCPU,useRhpc=useRhpc)
            #aa[[multi]] <- retJun$aa
            difp[[multi]] <- retJun$difp
            difp0[[multi]] <- retJun$difp0
            #difpraw[[multi]] <- ret$difpraw
            exyz0[[multi]] <- retJun$exyz0
            exyz1[[multi]] <- retJun$exyz1
            #exyout[[multi]] <- ret$exyout
            pod[[multi]] <- retJun$pod
            pot[[multi]] <- retJun$pot
            vb[[multi]] <- retJun$vb
            vb0[[multi]] <- retJun$vb0
            #
            if (acqmode=="bypass"){
              exyz1[[multi]] <- exyz0[[multi]]
              for (j in 1:ncur){
                exyz1[[multi]][ , , , ,(j-1-1)%%32+1] <- retJun$exyz2[ , , , ,(j-1-1)%%32+1]
                exyz1[[multi]][ , , , ,(j-0-1)%%32+1] <- retJun$exyz1[ , , , ,(j-0-1)%%32+1]
                exyz1[[multi]][ , , , ,(j+1-1)%%32+1] <- retJun$exyz3[ , , , ,(j+1-1)%%32+1]
              }
            } 
        }
        rm(retJun)
        
        dev <- sum(sapply(difp,sumup))/nmulti
        cat("--- Finished ---\n")
        cat("facl,dev,sigmaxy[1]",facl,dev,sigmaxy[1],"\n")
    #}
    iter <- 0
    nfix <- 0
    cdev <-NULL
    nee <- ifelse(heterogeneity,sum(nexc[[1]]),sum(nexc[[1]][1:6])+64)
    bb <- NULL
    dd <- NULL
    rank <- NULL
    aic0 <- 100000000
    #sigmaxy [1:nex0] <- sigmaz[1:nex0] <- 0.0005       ###
    #sigmaxy [nexx:nex] <- sigmaz[nexx:nex] <- 0.00005 # 0.00002  ###
    dev0 <- dev  ###
    sigma <- sqrt((2*sigmaxy^2+sigmaz^2)/3)
    etime <- (Sys.time())
    radChg <- 1
    boostFlag <- FALSE
    resoAlpha <- 0.01
    etime <- Sys.time()
    Conds0 <- NULL
    for (i in 1:sum(unlist(perm))) Conds0 <- c(Conds0,list(sigma0)) 
    #Rhpc::Rhpc_finalize()
##################### 3000 ###############################################
    repeat {
      #Rprof()
      #browser()
      if (useCluster & useRhpc){
        #Rhpc::Rhpc_initialize()
        #cl2 <<- Rhpc::Rhpc_getHandle(nCPU)
        gc(reset=TRUE)
        gc(reset=TRUE)
      }
      if (resEnh & dev<1500) {
        j <- 2
        resEnhLocal <- TRUE 
      } else {
        j <- 2
        resEnhLocal <- FALSE
      }
      if (resEnhLocal) {
        rndelelocal <- rndele
      }else{
        rndelelocal <- list(list(first=rep(0,1),second=rep(0,1),third=rep(0,1),fourth=rep(0,1),fifth=rep(0,1),sixth=rep(3,1)),
                           list(first=rep(0,1),second=rep(0,1),third=rep(0,1),fourth=rep(0,1),fifth=rep(0,1),sixth=rep(3,1)))
      }
      if (resEnhLocal == TRUE) cat("with Resolutoin enhancement \n")
      #browser()
      cdev <- c(cdev,dev)
      cat("--- Computing the sensitivity of conductivity --- ")
      gc();    gc()
      #browser()
      for (multi in 1:nmulti) {
          cat(eleint[[multi]],"multi:",multi,"acqmode; ",acqmode[[multi]],"\n")
          retSens <- senseCalc(mxnode=mxnode,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,pod=pod[[multi]],jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,nex=nex,nez=nez,vol=vol,exyz0=exyz0[[multi]],exyz1=exyz1[[multi]])
          sensemat[[multi]] <- retSens$sensemat
          sensematxy[[multi]] <- retSens$sensematxy
          sensematz[[multi]] <- retSens$sensematz
      }
      rm(list=c("retSens"))
      gc();    gc()
      cat("--- Computing the sensitivity of subcutanious fat area ---\n")
      ######################## modification of subcutaneous fat area ###########
      if(heterogeneity==FALSE) {
        retRad <- updateradius(sigmaxy=sigmaxy,sigmaz=sigmaz,exyz0=exyz0,exyz1=exyz1,bound32=bound32,difp=difp,nmulti=nmulti,nez=nez,nexx=nexx[[1]],nexr=nexr,xcod=xcod,ycod=ycod,zcod=zcod,nodex=nodex,nodez=nodez,zrange=zrange,mxbnd=mxbnd,mxnode=mxnode,mxn=mxn,mzn=mzn,mxb=mxb,mze=mze,mxe=mxe,nzin=nzin,nzot=nzot,iter=iter,est64=est64,useRhpc=useRhpc,useCluster = useCluster)
        rsens <- retRad$rsens #* log(10) ######################################
        boundr <- retRad$boundr
        rm(list=c("retRad"))
      }
      rm(exyz0,exyz1)
      gc();    gc()
      #browser()
      if (anisotropy==TRUE){
	    senselist <- list(sensematxy,sensematz,sensemat)
	    } else if (anisotropy==FALSE) {
	    senselist <- list(sensemat)
      }
      if (heterogeneity==FALSE) senselist <- c(senselist,list(rsens))
####################################################################
		  cat("--- Updateing conductivity --- \n")
      if(resEnhLocal)    {
        if (useGPU){
          cat("Resolution enhancement using GPU \n")
        } else {
        cat("Resolution enhancement using CPU \n")
        }
      }
      #
      #designMat <-designMatsqrt <- NULL
      retConds <-  NULL
      normbxs  <- NULL
      Conds <- NULL
########################### loop start ##################
      pb8 <- txtProgressBar(min = 1, max = 8, style = 3)
      setTxtProgressBar(pb8, 1)
      #    
      ################## L-curve method ############################
      #browser()
      #
      if (resEnhLocal){
        cp <- cd <- cq <- NULL
        lambdas <- alam*c(0.5,0.75,1,2,5,10)
        for (alam2 in lambdas){
          offset <- list()
          for (i in 1:perm[[j]]){
            offset <- c(offset,list(list(offset=c(rndelelocal[[j]]$first[i],rndelelocal[[j]]$second[i],rndelelocal[[j]]$third[i],rndelelocal[[j]]$fourth[i],rndelelocal[[j]]$fifth[i],rndelelocal[[j]]$sixth[i]),
                                     cond=TRUE,alam=alam2)))
          }
         
          offset <- c(offset, list(list(offset=rep(0,6),cond=FALSE,alam=alam2)))
          #    
          if (useRhpc & useCluster){
            #Rhpc::Rhpc_initialize()
            #cl2 <-Rhpc::Rhpc_getHandle(nCPU)
            tm <- proc.time()
            retCond <- parLapply2(cl2,x=offset,fun=updatecond,param=senselist,nee=nee,bb=bb,dd=dd,rank=rank,nmulti=nmulti,nexe=nexe[[j]],nexx=nexx[[j]],nex=nex,nexr=nexr,nexc=nexc[[j]],
                                 resEnhRatio=resEnhRatio[[j]],numpol=numpol,kvol=kvol,jcur=jcur,maxcur=maxcur,difp=difp,vol=vol,aniso=anisotropy,
                                 hetero=heterogeneity,est64=est64,iter=iter,alpha=alpha,radChg=radChg,dev=dev0,Marquardt=Marquardt,sigma0=sigma0)
            print(proc.time()-tm)
            #Rhpc::Rhpc_finalize()
          } else {
            #Rmpi::mpi.universe.size()-1
            cl <- parallel::makeCluster(parallel::detectCores(),type="PSOCK")
            tm <-snow.time(
            retCond <- parallel::parLapply(cl,X=offset,fun=updatecond,param=senselist,nee=nee,bb=bb,dd=dd,rank=rank,nmulti=nmulti,nexe=nexe[[j]],nexx=nexx[[j]],nex=nex,nexr=nexr,nexc=nexc[[j]],
                      resEnhRatio=resEnhRatio[[j]],numpol=numpol,kvol=kvol,jcur=jcur,maxcur=maxcur,difp=difp,vol=vol,aniso=anisotropy,
                      hetero=heterogeneity,est64=est64,iter=iter,alpha=alpha,radChg=radChg,dev=dev0,Marquardt=Marquardt,sigma0=sigma0)
            )
            dev.set(9)
            par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(4,4,2,1))
            plot(tm)
            stopCluster(cl)
          }
          #retCond <- lapply(X=offset,FUN=updatecond,param=senselist,nee=nee,bb=bb,dd=dd,rank=rank,alam=alam,nmulti=nmulti,nexe=nexe[[j]],nexx=nexx[[j]],nex=nex,nexr=nexr,nexc=nexc[[j]],
        #                          resEnhRatio=resEnhRatio[[j]],numpol=numpol,kvol=kvol,jcur=jcur,maxcur=maxcur,difp=difp,vol=vol,aniso=anisotropy,
        #                          hetero=heterogeneity,est64=est64,iter=iter,alpha=alpha,radChg=radChg,dev=dev0,Marquardt=Marquardt,sigma0=sigma0)
		    #cat("--- FINISHED --- \n")
        setTxtProgressBar(pb8, 2)
########################### Resolution enhancement ##############################
        #
        #
        setTxtProgressBar(pb8, 4)
        #
        #browser()
        cumcond <- NULL
        for(i in 1:perm[[j]]){
            cumcond <- cbind(cumcond,retCond[[i]]$bx)
        }
        q2 <- apply(abs(cumcond),1,which.max)
        q3 <- NULL
        for (i in 1:708) {
            q3 <- c(q3,cumcond[i,q2[i]])
        }
        sigma <- sigmaxy <- sigmaz <- sigma0*exp(q3)
        #
        cat("--- Solving forward problem ---\n")
        difptemp <- difp
        for (multi in 1:nmulti) {
            cat(eleint[[multi]],"multi:",multi,"acqmode; ",acqmode[[multi]],"\n")
            #retjun <- junPA(npat=1,iter=0,aa0=aa[[multi]],mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster = useCluster,nCPU=nCPU,useRhpc=useRhpc)
            retjun <- junPA(npat=1,iter=0,mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster = useCluster,nCPU=nCPU,useRhpc=useRhpc)
            difptemp[[multi]] <- retjun$difp
        }
        rm(list=c("retjun"))
        gc();    gc()
        cat("--- Finished forward problem---\n")
        dev <- sum(sapply(difptemp,sumup2))/nmulti
        cp <- c(cp,alam2)
        cd <- c(cd,dev)
        cq <- c(cq,1*norm(as.matrix(q3),type="f"))
        }
        #
        #browser()
        dev.set(7)
        plot(x=(cd),y=log10(cq),type="p",bty="l",ylab="model",xlab="devience",pch=20)
        text((cd),log10(cq),labels=round(cp,3),adj=c(1,2),cex=0.75)
        f <- interpSpline((cd),log10(cq))
        g <- smooth.spline((cd),log10(cp))
        xx <- seq(min((cd)), max((cd)),length.out = 101)
        lines(predict(f,xx),col = "yellow", lwd = 1.5)
        fx <- predict(f,xx,deriv=1)$y
        fxx <- predict(f,xx,deriv=2)$y
        #browser()
        curvetures <- curveture(fx=fx,fxx=fxx)
        curvetures[1:50] <- NA
        curvetures[length(curvetures)] <- NA
        newAlpha <- 10^(predict(g,xx[which.min(curvetures)])$y)
        cat("alpha=",round(alam,3),round(newAlpha,3),"\n")
        #
        offset <- list()
        for (i in 1:perm[[j]]){
            offset <- c(offset,list(list(offset=c(rndelelocal[[j]]$first[i],rndelelocal[[j]]$second[i],rndelelocal[[j]]$third[i],rndelelocal[[j]]$fourth[i],rndelelocal[[j]]$fifth[i],rndelelocal[[j]]$sixth[i]),
                                         cond=TRUE,alam=newAlpha)))
        }
        offset <- c(offset, list(list(offset=rep(0,6),cond=FALSE,alam=newAlpha)))
        #
        if (useRhpc & useCluster){
          #Rhpc::Rhpc_initialize()
          #cl2 <-Rhpc::Rhpc_getHandle(nCPU)
          retCond <- parLapply2(cl2,x=offset,fun=updatecond,param=senselist,nee=nee,bb=bb,dd=dd,rank=rank,nmulti=nmulti,nexe=nexe[[j]],nexx=nexx[[j]],nex=nex,nexr=nexr,nexc=nexc[[j]],
                               resEnhRatio=resEnhRatio[[j]],numpol=numpol,kvol=kvol,jcur=jcur,maxcur=maxcur,difp=difp,vol=vol,aniso=anisotropy,
                               hetero=heterogeneity,est64=est64,iter=iter,alpha=alpha,radChg=radChg,dev=dev0,Marquardt=Marquardt,sigma0=sigma0)
          #Rhpc::Rhpc_finalize()
        } else {
          cl <- parallel::makeCluster(parallel::detectCores(),type="PSOCK")
          retCond <- parallel::parLapply(cl,X=offset,fun=updatecond,param=senselist,nee=nee,bb=bb,dd=dd,rank=rank,nmulti=nmulti,nexe=nexe[[j]],nexx=nexx[[j]],nex=nex,nexr=nexr,nexc=nexc[[j]],
                             resEnhRatio=resEnhRatio[[j]],numpol=numpol,kvol=kvol,jcur=jcur,maxcur=maxcur,difp=difp,vol=vol,aniso=anisotropy,
                             hetero=heterogeneity,est64=est64,iter=iter,alpha=alpha,radChg=radChg,dev=dev0,Marquardt=Marquardt,sigma0=sigma0)
          stopCluster(cl)
        }
        #retCond <- lapply(X=offset,FUN=updatecond,param=senselist,nee=nee,bb=bb,dd=dd,rank=rank,alam=alam,nmulti=nmulti,nexe=nexe[[j]],nexx=nexx[[j]],nex=nex,nexr=nexr,nexc=nexc[[j]],
        #                          resEnhRatio=resEnhRatio[[j]],numpol=numpol,kvol=kvol,jcur=jcur,maxcur=maxcur,difp=difp,vol=vol,aniso=anisotropy,
        #                          hetero=heterogeneity,est64=est64,iter=iter,alpha=alpha,radChg=radChg,dev=dev0,Marquardt=Marquardt,sigma0=sigma0)
        #cat("--- FINISHED --- \n")
        setTxtProgressBar(pb8, 2)
        ########################### Resolution enhancement ##############################
        #
          #
        setTxtProgressBar(pb8, 4)
          #
          #browser()
        cumcond <- NULL
        for(i in 1:perm[[j]]){
            cumcond <- cbind(cumcond,retCond[[i]]$bx)
        }
        #
        q2 <- apply(abs(cumcond),1,which.max)
        q3 <- NULL
        for (i in 1:708) {
            q3 <- c(q3,cumcond[i,q2[i]])
        }
        sigma <- sigmaxy <- sigmaz <- sigma0*exp(q3)
        #
        cat("--- Solving forward problem ---\n")
        difptemp <- difp
        for (multi in 1:nmulti) {
          cat(eleint[[multi]],"multi:",multi,"acqmode; ",acqmode[[multi]],"\n")
          #retjun <- junPA(npat=1,iter=0,aa0=aa[[multi]],mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster = useCluster,nCPU=nCPU,useRhpc=useRhpc)
          retjun <- junPA(npat=1,iter=0,mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster = useCluster,nCPU=nCPU,useRhpc=useRhpc)
          difptemp[[multi]] <- retjun$difp
        }
        rm(list=c("retjun"))
        gc();    gc()
        cat("--- Finished forward problem---\n")
        dev <- sum(sapply(difptemp,sumup2))/nmulti
        cp <- c(newAlpha,cp)
        cd <- c(dev,cd)
        cq <- c(norm(as.matrix(q3),type="f"),cq)
        #browser()
        #
        dev.set(7)
        #browser()
        #plot(x=log10(cd),y=log10(cq),type="p",bty="l",ylab="model",xlab="devience",ylim=c(-9,5),pch=20,col="green")
        points(x=(cd[1]),y=log10(cq[1]),pch=20,col="blue")
        text((cd[1]),log10(cq[1]),labels=round(cp[1],3),adj=c(1,2),cex=0.75,col="blue")
        ord <- order(cp)
        cp <- cp[ord]
        cd <- cd[ord]
        cq <- cq[ord]
        f <- interpSpline((cd),log10(cq))
        g <- smooth.spline((cd),log10(cp))
        xx <- seq(min((cd)), max((cd)),length.out=101)
        lines(predict(f,xx), col = "green", lwd = 1.5)
        #xx <- log10(cd)
        fx <- predict(f,xx,deriv=1)$y
        fxx <-predict(f,xx,deriv=2)$y
        #browser()
        curvetures <- curveture(fx=fx,fxx=fxx)
        curvetures[1:50] <- NA
        curvetures[length(curvetures)] <- NA
        print(curvetures[51:length(curvetures)])
        newAlpha <- 10^(predict(g,xx[which.min(curvetures)])$y)
        #
        text(x=xx[which.min(curvetures)], y=predict(f,xx[which.min(curvetures)])$y,labels=round(newAlpha,3),adj=c(1,-2),cex=0.75,col="red")
        points(x=xx[which.min(curvetures)], y=predict(f,xx[which.min(curvetures)])$y, pch=20, col="red")
        #
        cat("alpha=",round(alam,3),round(newAlpha,3),"\n")
        #
        alam <- ifelse(newAlpha<5*alam & newAlpha>0.75*alam,newAlpha,alam)
        #
      }
      ############# L-curve criterion END ######################################################################################## ***
      #
      #
      setTxtProgressBar(pb8, 4)
          #
      #browser()
          #
        offset <- list()
        if (resEnhLocal){
          for (i in 1:perm[[j]]){
            offset <- c(offset,list(list(offset=c(rndelelocal[[j]]$first[i],rndelelocal[[j]]$second[i],rndelelocal[[j]]$third[i],rndelelocal[[j]]$fourth[i],rndelelocal[[j]]$fifth[i],rndelelocal[[j]]$sixth[i]),
                                         cond=TRUE,alam=alam)))
          }
        } else {
          offset <- c(offset,list(list(offset=c(rndelelocal[[j]]$first[1],rndelelocal[[j]]$second[1],rndelelocal[[j]]$third[1],rndelelocal[[j]]$fourth[1],rndelelocal[[j]]$fifth[1],rndelelocal[[j]]$sixth[1]),
                                       cond=TRUE,alam=alam)))
        }
        offset <- c(offset, list(list(offset=rep(0,6),cond=FALSE,alam=alam)))
        #  
        if (useRhpc & useCluster){
          #Rhpc::Rhpc_initialize()
          #cl2 <-Rhpc::Rhpc_getHandle(nCPU)
          retCond <- parLapply2(cl2,x=offset,fun=updatecond,param=senselist,nee=nee,bb=bb,dd=dd,rank=rank,nmulti=nmulti,nexe=nexe[[j]],nexx=nexx[[j]],nex=nex,nexr=nexr,nexc=nexc[[j]],
                               resEnhRatio=resEnhRatio[[j]],numpol=numpol,kvol=kvol,jcur=jcur,maxcur=maxcur,difp=difp,vol=vol,aniso=anisotropy,
                               hetero=heterogeneity,est64=est64,iter=iter,alpha=alpha,radChg=radChg,dev=dev0,Marquardt=Marquardt,sigma0=sigma0)
          #Rhpc::Rhpc_finalize()
        } else {
          cl <- makeCluster(parallel::detectCores(),type="PSOCK")
          retCond <- parallel::parLapply(cl,X=offset,fun=updatecond,param=senselist,nee=nee,bb=bb,dd=dd,rank=rank,nmulti=nmulti,nexe=nexe[[j]],nexx=nexx[[j]],nex=nex,nexr=nexr,nexc=nexc[[j]],
                               resEnhRatio=resEnhRatio[[j]],numpol=numpol,kvol=kvol,jcur=jcur,maxcur=maxcur,difp=difp,vol=vol,aniso=anisotropy,
                               hetero=heterogeneity,est64=est64,iter=iter,alpha=alpha,radChg=radChg,dev=dev0,Marquardt=Marquardt,sigma0=sigma0)
          stopCluster(cl)
        }
        #
        if (resEnhLocal){   
          cumcond <- NULL
          for(i in 1:perm[[j]]){
            cumcond <- cbind(cumcond,retCond[[i]]$bx)
          }
          q2 <- apply(abs(cumcond),1,which.max)
          q3 <- NULL
          for (i in 1:708) {
            q3 <- c(q3,cumcond[i,q2[i]])
          }
          sigma <- sigmaxy <- sigmaz <- sigma0*exp(q3)
          #
        } else {
          sigma <- sigmaxy <- sigmaz <- sigma0*exp(retCond[[1]]$bx)
        }
        cat("updated the conductivities \n")
        #cat("selected",sel,"\n")
        #if (sel > perm[[1]]) {
        #  k <- 2
        #  sel <- sel - perm[[1]]
        #} else {
        #  k <- 1
        #}
        #cat("Pattern",j,k,"\n")
        #cat(rndele[[k]]$first[sel],rndele[[k]]$second[sel],rndele[[k]]$third[sel],rndele[[k]]$fourth[sel],rndele[[k]]$fifth[sel],rndele[[k]]$sixth[sel],"\n")
#########################################################################################
      dd <- retCond[[1]]$dd
      rank<- retCond[[1]]$rank
      #nee <- retCond$nee
      bb <- retCond[[1]]$bb
      aic <- (32*log(dev/32) + nee*2)
      cat("\n AIC;",aic,"\n")
      #if (aic<aic0) nee <- min(nexe+64,nee+1)
      #else nee <- nee-1
      aic0 <- aic
      if (heterogeneity==FALSE){
          #browser()
          #dd <- r$d[1:16]*sum(r$d[1:16])/sum(r$d[1:16])
          #dR <- -r$v%*%diag(1/dd)%*%t(r$u)%*%difpmul*1/1
          #cat("rbx",round(ret$rbx,4),"\n")
          #if(iter>0) dRm <- (ret$rbx/20)
          #if(iter>4) dRm <- exp(retCond$rbx*boundr/10000/log(100))
          if(iter>2) dRm <- exp(retCond[[length(offset)]]$rbx/boundr)
          #else if(iter > 1)      dRm <- exp(retCond$rbx*boundr/10000/10/log(10)) #
          else dRm <- rep(1.0,length(retCond[[length(offset)]]$rbx))
          #rm(list=c("retCond"))
          #dRm <- dR
          rad0 <- rad64
          if(est64==TRUE){
          rad64 <- as.vector(rad0*dRm)
          rad64 <- ifelse(rad64>0.95,0.95,rad64)
          } else {
          rad32 <- rad0[seq(1,63,2)]*dRm
          #rad32 <- rad0[seq(1,63,2)]+dRm
          xx <- seq(0, 2*pi, length.out = 33 )[1:32]
          yy <- rad32[1:32] #c(rad32[32],rad32[1:31])
          pispl <- periodicSpline( xx, yy, period = 2 * pi,ord=2)
          rad64 <- predict( pispl ,seq(0, 2*pi, length.out = 65 ))$y
          }
          #B <- t(S)%*%(as.vector(difpmul))
          #dR <- -solve(a=t(S)%*%S,b=B)*1
          #dR <- dR[c(2:32,1)]
          #cat("dR",dR,"\n")
          #browser()
          radChg <- norm(matrix(rad64-rad0),type="f")
          cat("Change of radius;",radChg, "\n")
          cat("conductivity of subcutaneous fat;", sigmaxy[mxe], "\n")
          #cat(sigmaxy[mxe],"\n")
          #cat("modification of radius","\n")
          #cat(ifelse(dRm > 1,"+","-"),"\n")
          #cat("renewal; rad64\n")
          #cat(" 1-16:",format(rad64[1:16],nsmall=2,digits=2),"\n")
          #cat("17-32:",format(rad64[17:32],nsmall=2,digits=2),"\n")
          #cat("33-48:",format(rad64[33:48],nsmall=2,digits=2),"\n")
          #cat("49-64:",format(rad64[49:64],nsmall=2,digits=2),"\n")
          #
          if (dev < 100 & boostFlag==FALSE) boostFlag <- TRUE
          layer <- adjLayer(layer=layer,sensemat=sensemat[[1]],nexr=nexr,nexc=nexc[[1]],rad64=rad64,arar=arar,resEnhRatio=resEnhRatio[[1]])
          cat("Radius of layers,",layer,"\n")
          #
          cat("---Generating a mesh--- \n")
          ret <- meshgen(rad=rad64,zrange=zrange,mxbnd=mxbnd,mxnode=mxnode,bndd=bndd,mxn=mxn,mzn=mzn,mxb=mxb,mze=mze,mxe=mxe,zcod,nzin,nzot,zcalc=FALSE,layer=layer)
          #theta <- ret$theta
          rad32 <- ret$rad
          xcod <- ret$xcod
          ycod <- ret$ycod
          xcod2 <- ret$xcod2
          ycod2 <- ret$ycod2
          zcod <- ret$zcod
          arar <- ret$arar
          xb <- ret$xb
          yb <- ret$yb
          bnd <- ret$bnd
          #xav <- ret$xav
          #yav <- ret$yav
          st <- ret$st
          stxy <- ret$stxy
          stz <- ret$stz
          cex <- ret$cex
          vol <- ret$vol
          nez <- ret$nez
          rm(list=c("ret"))
          gc();    gc()
          #cat("--- FINISHED, updateradius.---\n")
      }
      setTxtProgressBar(pb8, 6)
		###################################################################
      rm(retCond)
      gc(); gc()
    #browser()
    substring(filenm,6,8) <- paste(substring("00",nchar(iter),2),iter,sep="")
    iter <- iter +1
    sigma0 <- sigma
    currentTime <- Sys.time()
    cat("Cycle time;",difftime(currentTime,etime),"\n")
    etime <- currentTime
    #if(iter%%10==1 | iter==maxiter) {
    if(TRUE) {
      resPlot <- plotstatusA(filename=filenm,devmat=cdev,alammat=c(0,alam0,alam),dd=list(dd=c(dd,NULL)),nee=list(tranc=c(nee),rank=c(rank)),iter=iter,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,xcod=xcod,ycod=ycod,zcod=zcod,xcod2=xcod2,ycod2=ycod2,nodex=nodex,nodez=nodez,nnodex=nnodex,nnodez=nnodez,nex=nex,nez=nez,nexr=nexr,eleint=eleint,sense=sensemat[[1]],sensexy=sensematxy[[1]],sensez=sensematz[[1]],idval=idval,arar=arar,vol=vol,bound32=bound32,aniso=anisotropy,nb=nb,ibnd=ibnd,initial=initial,estres=estres,xb=xb,yb=yb,aic=aic,annotate=FALSE)
      if(!is.null(resPlot)) {
        estres <- resPlot$estres
        initial <- resPlot$initial
      }
      rm(list="resPlot")
    }
    for(nn in 1:2){
            if(FALSE) rgl.set(nn)
            #play3d(spin3d(axis=c(0,1,1)),duration=1)
    }
            #nik <- nik +1
            #nskip <- 0
###### 2000 #################################################################
            #repeat{
                qflag <- 0
                #dev0 <- dev
                cat("--- Solving forward problem --- ")
                exyz0 <- rep(list(array(0, c(3,mze,mxe,jcur,kvol))),length(expout))
                exyz1 <- rep(list(array(0, c(3,mze,mxe,jcur,kvol))),length(expout))
                for (multi in 1:nmulti) {
                    cat(eleint[[multi]],"multi:",multi,"acqmode; ",acqmode[[multi]],"\n")
                    #retJun <- junPA(npat=1,iter=0,aa0=aa[[multi]],mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster = useCluster,nCPU=nCPU,useRhpc=useRhpc)
                    retJun <- junPA(npat=1,iter=0,mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster = useCluster,nCPU=nCPU,useRhpc=useRhpc)
                    #aa[[multi]] <- retJun$aa
                    difp[[multi]] <- retJun$difp
                    difp0[[multi]] <- retJun$difp0
                                        #difpraw[[multi]] <- ret$difpraw
                    exyz0[[multi]] <- retJun$exyz0
                    exyz1[[multi]] <- retJun$exyz1
                                        #exyout[[multi]] <- ret$exyout
                    pod[[multi]] <- retJun$pod
                    pot[[multi]] <- retJun$pot
                    vb[[multi]] <- retJun$vb
                    vb0[[multi]] <- retJun$vb0
                    #
                    if (acqmode=="bypass"){
                      exyz1[[multi]] <- exyz0[[multi]]
                      for (j in 1:ncur){
                        exyz1[[multi]][ , , , ,(j-1-1)%%32+1] <- retJun$exyz2[ , , , ,(j-1-1)%%32+1]
                        exyz1[[multi]][ , , , ,(j-0-1)%%32+1] <- retJun$exyz1[ , , , ,(j-0-1)%%32+1]
                        exyz1[[multi]][ , , , ,(j+1-1)%%32+1] <- retJun$exyz3[ , , , ,(j+1-1)%%32+1]
                      }
                    }
                }
                rm(retJun)
                cat("--- Finished forward problem--- iteration,",iter,"\n")
                dev <- sum(sapply(difp,sumup))/nmulti
                cat("(dev0-dev)/dev0,dev0,dev,alam",(dev0-dev)/dev0,dev0,dev,alam,"\n")
                #
                if ((dev0-dev)/dev0 <= eps & iter >5) {
                    nfix <- nfix+1  #3
                    cat("(dev0-dev)/dev <= eps  !!!\n")
                    #if (nfix==1)  alam <- min(alpha*alam0,alpha*dev^(1/2))
                    #else 
                    alam <- min(2*alam0, alam/redalam^3)  #^3
                    cat("Beta will be modified to;",alam,"\n")
                    #nee <- nee-10
                    cat("retreive conductivity \n")
                    sigmaxy[1:nex] <- sigma0xy[1:nex]
                    sigmaz[1:nex] <- sigma0z[1:nex]
                    sigma <- sqrt((2*sigmaxy^2+sigmaz^2)/3)
                    for (numlti in 1:nmulti) vb[[multi]][1:(nb+1),1:jcur] <- vb0[[multi]][1:(nb+1),1:jcur]
                    if (heterogeneity==FALSE) {
                    rad64 <- rad0
                    rad0 <- rad64
                    cat("retreive radius \n")
                    ret <- meshgen(rad=rad64,zrange=zrange,mxbnd=mxbnd,mxnode=mxnode,bndd=bndd,mxn=mxn,mzn=mzn,mxb=mxb,mze=mze,mxe=mxe,zcod,nzin,nzot,zcalc=FALSE,layer=layer)
                    #theta <- ret$theta
                    xcod <- ret$xcod
                    ycod <- ret$ycod
                    xcod2 <- ret$xcod2
                    ycod2 <- ret$ycod2
                    arar <- ret$arar
                    xb <- ret$xb
                    yb <- ret$yb
                    bnd <- ret$bnd
                    #xav <- ret$xav
                    #yav <- ret$yav
                    st <- ret$st
                    stxy <- ret$stxy
                    stz <- ret$stz
                    cex <- ret$cex
                    vol <- ret$vol
                    rm(list=c("ret"))
                    gc();    gc()
                    #
                    cat("--- FINISHED, retrieve radius.---\n")
                    }
                    cat("--- Solving forward problem ---\n")
                    for (multi in 1:nmulti) {
                        cat(eleint[[multi]],"multi:",multi,"acqmode; ",acqmode[[multi]],"\n")
                        #retjun <- junPA(npat=1,iter=0,aa0=aa[[multi]],mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster = useCluster,nCPU=nCPU,useRhpc=useRhpc)
                        retjun <- junPA(npat=1,iter=0,mxbnd=mxbnd,mxnode=mxnode,maxcur=maxcur,numpol=numpol,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,iinn=iinn[[multi]],iott=iott[[multi]],pod=pod[[multi]],podex=podex[[multi]],vb=vb[[multi]],difp=difp[[multi]],difp0=difp0[[multi]],nnodex=nnodex,nodex=nodex,nodez=nodez,nnodez=nnodez,nnodezx=nnodezx,st=st,stxy=stxy,stz=stz,nbw=nbw,nb=nb,ibnd=ibnd,bnd=bnd,vol=vol,cex=cex,nex=nex,nez=nez,nzin=nzin,nzot=nzot,jcur=jcur,kvol=kvol,mxe=mxe,mze=mze,mxb=mxb,acqmode=acqmode[[multi]],useCluster = useCluster,nCPU=nCPU,useRhpc=useRhpc)
                        #aa[[multi]] <- retjun$aa
                        difp[[multi]] <- retjun$difp
                        difp0[[multi]] <- retjun$difp0
                                        #difpraw[[multi]] <- ret$difpraw
                        exyz0[[multi]] <- retjun$exyz0
                        exyz1[[multi]] <- retjun$exyz1
                                        #exyout[[multi]] <- ret$exyout
                        pod[[multi]] <- retjun$pod
                        pot[[multi]] <- retjun$pot
                        vb[[multi]] <- retjun$vb
                        vb0[[multi]] <- retjun$vb0
                    }
                    rm(list=c("retjun"))
                    gc();    gc()
                    cat("--- Finished forward problem---\n")
                    dev <- sum(sapply(difp,sumup))/nmulti
                    cat("iter,dev,alam",iter,dev,alam,"\n")
                    dev0 <- dev
                    qflag <- 5000
                    #alam <- min(alam0, alam/redalam/redalam)
                    cat("nfix=",nfix, "\n")
                } else {
                    #browser()
                    #nfix <- 0
                    dev0 <- dev
                    if(Marquardt) alam <- max(alam * redalam, alpha*alamlim)
                    else {
                      if (nfix==0) alam <- max(alpha*norm(as.matrix(na.omit(do.call(c,as.vector(difp)))),type="f")^(2/3),alpha*alamlim)
                      #if (nfix==0) alam <- max(alpha*norm(as.matrix(na.omit(do.call(c,as.vector(difp)))),type="f")^(2/3)*8,alpha*alamlim)
                      else  {
                        #alam <- alam*redalam
                        nfix <- nfix-1
                        cat("nfix=",nfix, "\n")
                      }
                      #alam <- ifelse(dev>1,max(alam*redalam, alpha*dev^(1/2)),max(alpha*dev,alpha*alamlim))
                      #
                      #A <- retCond$A
                      #x <- retCond$x
                      #a <- as.vector(t(A)%*%A%*%x) 
                      #b <- as.vector(x) #norm(sensmulti[,-(nexx[6]+1)],type="2")
                      #alam <- ((-sum(a*b)+sqrt(sum(a*b)^2-sum(b^2)*(sum(a^2)-sum((t(A)%*%(alpha*do.call(c,as.vector(difp))))^2))))/sum(b^2))
                      #alam <- ((-sum(a*b)+sqrt(sum(a*b)^2-sum(b^2)*(sum(a^2)-(alpha*norm(t(A)%*%do.call(c,as.vector(difp0)),type="f"))^2)))/sum(b^2))
                    }
                    cat("alpha=",alam,"\n")
                    if (alam<0 | is.na(alam)) alam <- 10
                    #
                    cat("Beta will be modified to;",alam,"\n")
                    for (numlti in 1:nmulti) vb0[[multi]][1:(nb+1),1:jcur] <- vb[[multi]][1:(nb+1),1:jcur]
                    sigma0xy <- sigmaxy
                    sigma0z <- sigmaz
                    sigma0 <- sqrt((2*sigma0xy^2+sigma0z^2)/3) # average conductivity
                    rad0 <- rad64
                    qflag <- 5000
                    #cat("will modify the alam and radius, goto 3000 \n")
                }
############## 5000 #############################################################
     repeat{
                #Rprof(NULL)
                #browser()
                #if (useCluster & useRhpc) Rhpc::Rhpc_finalize()
                #summaryRprof()
                if ((iter >= maxiter) | ((dev0-dev)/dev <= eps & nfix > 30)){
                    devf <- dev
                    qflag <- 10000
                    break
                }
                if (qflag==6000) break
                #if (nee > nexx[6]) nee <- nexx[6]
                break
        } # 5000
############## 6000 ###############################
            if (qflag == 10000) break
            qflag <- ifelse(dev > dev0,3000,3000)
        if (qflag == 10000) break
    } # 3000
    devf <- dev
#
    ani.options(nmax = 240, ani.width=320,ani.height=300)
    saveGIF(ani.replay(), movie.name= "movie.gif", img.name= "record_plot", interval = 0.5, movietype = "gif", outdir = getwd())
    #
    saveVideo(ani.replay(), video.name= "movie.mp4", img.name= "record_plot", other.opts = "-pix_fmt yuv420p -b 300k",interval = 0.45)
    #ani.options(convert = "C:/Program Files/ImageMagick-6.9.2-Q16/convert.exe")
    files <- list.files(pattern="*.PGM")
    demon <- function(files) { 
    for (sfile in files) {
      pgmimg <- read.pnm(sfile)
      plot(pgmimg)
    }
  }
  saveGIF(demon(files), movie.name= "movie_smooth.gif", interval = 0.5, movietype = "gif", outdir = getwd())
#
cat("---------------NEWTON END---------------------\n")
cat("---------------ITERATION END------------------\n")
iret <- 0
if (useCluster & useRhpc) Rhpc::Rhpc_finalize()
return(list(status=iret))
}
#########################################################################
#
#   END OF FUNCTION
#
#########################################################################
#########################################################################

#######################################################################################
findBorder <- function(param,height,xcod,ycod,nodes,element,sigmaxy,sigmaz,face,iter,est64,useRhpc,useCluster){
  exyz0 <- param[[1]]
  exyz1 <- param[[2]]
  #browser()
  rm(param)
  #cl5 <- makeCluster(detectCores())
  #registerDoParallel(detectCores())
  #if (useRhpc & useCluster){
  if (FALSE){
    #Rhpc::Rhpc_initialize()
    #cl3 <-Rhpc::Rhpc_getHandle(12)
    ignore <- Rhpc::Rhpc_EvalQ(cl=cl3, expr={library(AbdominalEITextra,quietly=TRUE); NULL},envir=as.environment(globalenv()))
    #tm <- proc.time()  
    interresl <- Rhpc::Rhpc_sapplyLB(cl=cl3,X=1:12,FUN=resist,simplify="array",exyz0=exyz0,exyz1=exyz1,height=height,xcod=xcod,ycod=ycod,nodes=nodes,element=element,sigmaxy=sigmaxy,sigmaz=sigmaz,face=face,iter=iter,est64=est64)
    #
    #Rhpc::Rhpc_finalize()
  #} else if (!useRhpc & useCluster){
  } else if (FALSE){
    cl1 <- snow::makeMPIcluster(12)
    ignore <- snow::clusterEvalQ(cl=cl1, expr={library(AbdominalEITextra,quietly=TRUE); NULL})
    tm <- snow.time( 
      interresl <- snow::parSapply(cl=cl1,X=1:12,FUN=resist,simplify="array",exyz0=exyz0,exyz1=exyz1,height=height,xcod=xcod,ycod=ycod,nodes=nodes,element=element,sigmaxy=sigmaxy,sigmaz=sigmaz,face=face,iter=iter,est64=est64)
    )
    #interresl <- sapply(X=1:12,FUN=resist,simplify="array",exyz0=exyz0,exyz1=exyz1,height=height,xcod=xcod,ycod=ycod,nodes=nodes,element=element,sigmaxy=sigmaxy,sigmaz=sigmaz,face=face,iter=iter,est64=est64)
    stopCluster(cl1)
    cat("\n")
    print(tm$elapsed)
  } else {
    #
    #browser()
    cl2 <- parallel::makeCluster(parallel::detectCores(),type="PSOCK")
    ignore <- parallel::clusterEvalQ(cl=cl2, expr={library(AbdominalEITextra,quietly=TRUE); NULL})
    tm <- snow.time( 
      interresl <- parallel::parSapply(cl=cl2,X=1:12,FUN=resist,simplify="array",exyz0=exyz0,exyz1=exyz1,height=height,xcod=xcod,ycod=ycod,nodes=nodes,element=element,sigmaxy=sigmaxy,sigmaz=sigmaz,face=face,iter=iter,est64=est64)
    )
    #interresl <- sapply(X=1:12,FUN=resist,simplify="array",exyz0=exyz0,exyz1=exyz1,height=height,xcod=xcod,ycod=ycod,nodes=nodes,element=element,sigmaxy=sigmaxy,sigmaz=sigmaz,face=face,iter=iter,est64=est64)
    stopCluster(cl2)
    cat("\n")
    print(tm$elapsed)
  }
  #
  #interresl <- sapply(X=1:12,FUN=resist,simplify="array",exyz0=exyz0,exyz1=exyz1,height=height,xcod=xcod,ycod=ycod,nodes=nodes,element=element,sigmaxy=sigmaxy,sigmaz=sigmaz,face=face,iter=iter,est64=est64)
  #interresl <- parSapply(cl=cl5,X=1:12,FUN=resist,simplify="array",exyz0=exyz0,exyz1=exyz1,height=height,xcod=xcod,ycod=ycod,nodes=nodes,element=element,sigmaxy=sigmaxy,sigmaz=sigmaz,face=face,iter=iter,est64=est64)
  #interresl <- foreach(j=1:12,.combine="+") %dopar% resist(j,exyz0=exyz0,exyz1=exyz1,height=height,xcod=xcod,ycod=ycod,nodes=nodes,element=element,sigmaxy=sigmaxy,sigmaz=sigmaz,face=face,iter=iter,est64=est64)
  #stopCluster(cl5)
  #stopImplicitCluster()
  return(list(interres=interresl))
}
####################################################################################
sharednodes <- function(e1,e2,z,nodex,xcod,ycod){
    snodes <- NULL
    e1nodes <- nodex[1:3,z,e1]
    e2nodes <- nodex[1:3,z,e2]
    if (identical(e1nodes[1],e2nodes[1])) snodes <- c(snodes,e1nodes[1])
    if (identical(e1nodes[1],e2nodes[2])) snodes <- c(snodes,e1nodes[1])
    if (identical(e1nodes[1],e2nodes[3])) snodes <- c(snodes,e1nodes[1])
    if (identical(e1nodes[2],e2nodes[1])) snodes <- c(snodes,e1nodes[2])
    if (identical(e1nodes[2],e2nodes[2])) snodes <- c(snodes,e1nodes[2])
    if (identical(e1nodes[2],e2nodes[3])) snodes <- c(snodes,e1nodes[2])
    if (identical(e1nodes[3],e2nodes[1])) snodes <- c(snodes,e1nodes[3])
    if (identical(e1nodes[3],e2nodes[2])) snodes <- c(snodes,e1nodes[3])
    if (identical(e1nodes[3],e2nodes[3])) snodes <- c(snodes,e1nodes[3])
    snode1 <- c(xcod[snodes[1]],ycod[snodes[1]])
    snode2 <- c(xcod[snodes[2]],ycod[snodes[2]])
    inner <- (snode2-snode1)%*%c(-snode1[2],snode1[1])
    #cat(inner)
    if(is.na(inner)) {
        browser()
        browser()
        cat(snode1,snode2,e1,e2)
    }
    if(inner<0) snodes <- rev(snodes)
    return(snodes)
}
#####################################################################################
FindElements <- function(nodes,elements,z,nodex){
    e1 <- rep(NA,2)
    e2 <- rep(NA,2)
    ele  <- elements[1]
    for (n in 1:2){
        node <- nodes[n]
        for (zz in (z-1)*3+1:3){
            if(sum(match(nodex[1:4,zz,ele],node),na.rm=TRUE)==2) e1[n] <- zz
        }
    }
    ele <- elements[2]
    for (n in 1:2){
        node <- nodes[n]
        for (zz in (z-1)*3+1:3){
            if(sum(match(nodex[1:4,zz,ele],node),na.rm=TRUE)==2) e2[n] <- zz
        }
    }
    return(c(e1,e2))
}

############################
# updateradius
#############################
updateradius <- function(sigmaxy,sigmaz,exyz0,exyz1,bound32,difp,nmulti,nez,nexx,nexr,xcod,ycod,zcod,nodex,nodez,zrange,mxbnd,mxnode,mxn,mzn,mxb,mze,mxe,nzin,nzot,iter,est64,useRhpc,useCluster){
    pb5 <- txtProgressBar(min = 1, max = 5, style = 3)
    setTxtProgressBar(pb5, 1) 
    #cat("optimization has been started.\n")
    interres <- rep(list(NA),length(difp))
    #dRm <- rep(NA,32)
    height <- zcod[2:13]-zcod[1:12]
    face <- array(NA,c(12,32,4,12))
    element <- matrix(NA,7,32)
    nodes <- array(NA,c(12,12,32,2))
    element[1,] <- (nexr[4]+1)+(0:31)*2 #A
    element[2,] <- (nexr[4]+2)+(0:31)*2 #B
    element[3,] <- (nexr[5]+2)+(0:31)*3 #C
    element[4,] <- (nexr[5]+1)+(0:31)*3 #D
    element[5,] <- (nexr[5]+3)+(0:31)*3 #E
    element[6,] <- (nexr[6]+2)+(0:31)*4 #F
    element[7,] <- (nexr[6]+3)+(0:31)*4 #G
    #
    ccolor <- c("brown","red","orange","yellow","blue","green","violet","gray","white","black")
    #
    bound <- bound32$boundry[seq(2,64,2),]
    boundr <- sqrt(bound[,"x"]^2+bound[,"y"]^2)
    if(FALSE)  {
      dev.set(11)
      abline(h=0,lty=3,col="red")
      abline(v=-20,lty=3,col="blue")
      for (i in 1:32){
        theta <- 2*pi/32*i-2*pi/64
        text(x=(boundr[i]+12)*sin(theta)-20,y=-(boundr[i]+10)*cos(theta),labels=i)
      }
      for(i in 1:7){
        eids <- element[i,]
        elements <-cbind(P1X=xcod[nodex[1,1,eids]],P1Y=ycod[nodex[1,1,eids]],P2X=xcod[nodex[2,1,eids]],P2Y=ycod[nodex[2,1,eids]],P3X=xcod[nodex[3,1,eids]],P3Y=ycod[nodex[3,1,eids]])
                                        #elements <- cbind(elements, COL=rep(128,nrow(elements)))
        ccx <- NULL
        ccy <- NULL
        for (eid in 1:nrow(elements)) {
            cx <- c(elements[eid,"P1X"],elements[eid,"P2X"],elements[eid,"P3X"])
            cy <- c(elements[eid,"P1Y"],elements[eid,"P2Y"],elements[eid,"P3Y"])
            ccx <- c(ccx,list(x=cx))
            ccy <- c(ccy,list(y=cy))
        }
        #mapply(FUN=polygon,x=ccx,y=ccy,col=ccolor[i],border=TRUE,lty="solid")
      }
    }
    for(j in 1:12){                          #
      for(i in seq(1,32,2)) nodes[1,j,i,] <- sharednodes(element[2,i],element[1,i],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #BA
      for(i in seq(2,32,2)) nodes[1,j,i,] <- sharednodes(element[1,i],element[2,i],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #AB
      #
      for(i in seq(1,32,2)) nodes[2,j,i,] <- sharednodes(element[1,i],element[3,i],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #AC
      for(i in seq(2,32,2)) nodes[2,j,i,] <- sharednodes(element[2,i],element[3,i],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #BC
      #
      for(i in 1:32) nodes[3,j,i,]        <- sharednodes(element[3,i],element[4,i],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #CD
      #
      for(i in 1:32) nodes[4,j,i,]        <- sharednodes(element[3,i],element[5,i],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #CE
      #
      for(i in 1:32) nodes[5,j,i,]        <- sharednodes(element[4,i],element[6,i],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #DF
      for(i in 1:32) nodes[6,j,i,]        <- sharednodes(element[5,i],element[7,i],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #EG
      for(i in seq(1,32,2)) nodes[7,j,i,] <- sharednodes(element[1,(30+i)%%32+1],element[2,(30+i)%%32+1],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #AB-
      for(i in seq(2,32,2)) nodes[7,j,i,] <- sharednodes(element[2,(30+i)%%32+1],element[1,(30+i)%%32+1],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #BA-
      for(i in seq(2,32,2)) nodes[8,j,i,] <- sharednodes(element[1,(30+i)%%32+1],element[3,(30+i)%%32+1],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #AC-
      for(i in seq(1,32,2)) nodes[8,j,i,] <- sharednodes(element[2,(30+i)%%32+1],element[3,(30+i)%%32+1],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #BC-
      for(i in 1:32) nodes[9,j,i,]        <- sharednodes(element[3,(30+i)%%32+1],element[4,(30+i)%%32+1],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #CD-
      for(i in 1:32) nodes[10,j,i,]       <- sharednodes(element[3,(30+i)%%32+1],element[5,(30+i)%%32+1],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #CE-
      for(i in 1:32) nodes[11,j,i,]       <- sharednodes(element[4,(30+i)%%32+1],element[6,(30+i)%%32+1],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #DF-
      for(i in 1:32) nodes[12,j,i,]       <- sharednodes(element[5,(30+i)%%32+1],element[7,(30+i)%%32+1],3*j-2,nodex=nodex,xcod=xcod,ycod=ycod) #EG-
    }#
    setTxtProgressBar(pb5, 2) 
    if (iter >1){
      dev.set(11)
      for (j in 5:5){
          for (i in 1:32){
           # lines(x=c(xcod[nodes[j,6,i,1]]-20,xcod[nodes[j,6,i,2]]-20),y=c(ycod[nodes[j,6,i,1]],ycod[nodes[j,6,i,2]]),lwd=3,col=ccolor[(i-1)%%10+1])
          }
      }                                   #
      for (i in 1:32){
      #  points(x=xcod[nodes[5,6,i,1]]-20,y=ycod[nodes[5,6,i,1]],pch=19,cex=1.0,col=ccolor[(i-1)%%10+1])
      }
    }                       #
    #browser()
    setTxtProgressBar(pb5, 3) 
    for(j in 1:12){
        for(i in seq(1,32,2)) face[1,i,1:4,j] <- FindElements(nodes=nodes[1,j,i,],elements=c(element[2,i],element[1,i]),z=j,nodex=nodex) #BA
        for(i in seq(2,32,2)) face[1,i,1:4,j] <- FindElements(nodes=nodes[1,j,i,],elements=c(element[1,i],element[2,i]),z=j,nodex=nodex) #AB
        for(i in seq(1,32,2)) face[2,i,1:4,j] <- FindElements(nodes=nodes[2,j,i,],elements=c(element[1,i],element[3,i]),z=j,nodex=nodex) #AC
        for(i in seq(2,32,2)) face[2,i,1:4,j] <- FindElements(nodes=nodes[2,j,i,],elements=c(element[2,i],element[3,i]),z=j,nodex=nodex) #BC
        for(i in 1:32)        face[3,i,1:4,j] <- FindElements(nodes=nodes[3,j,i,],elements=c(element[3,i],element[4,i]),z=j,nodex=nodex) #CD
        for(i in 1:32)        face[4,i,1:4,j] <- FindElements(nodes=nodes[4,j,i,],elements=c(element[3,i],element[5,i]),z=j,nodex=nodex) #CE
        for(i in 1:32)        face[5,i,1:4,j] <- FindElements(nodes=nodes[5,j,i,],elements=c(element[4,i],element[6,i]),z=j,nodex=nodex) #DF
        for(i in 1:32)        face[6,i,1:4,j] <- FindElements(nodes=nodes[6,j,i,],elements=c(element[5,i],element[7,i]),z=j,nodex=nodex) #EG
        for(i in seq(1,32,2)) face[7,i,1:4,j] <- FindElements(nodes=nodes[7,j,i,],elements=c(element[1,(30+i)%%32+1],element[2,(30+i)%%32+1]),z=j,nodex=nodex) #AB-
        for(i in seq(2,32,2)) face[7,i,1:4,j] <- FindElements(nodes=nodes[7,j,i,],elements=c(element[2,(30+i)%%32+1],element[1,(30+i)%%32+1]),z=j,nodex=nodex) #BA-
        for(i in seq(2,32,2)) face[8,i,1:4,j] <- FindElements(nodes=nodes[8,j,i,],elements=c(element[1,(30+i)%%32+1],element[3,(30+i)%%32+1]),z=j,nodex=nodex) #AC-
        for(i in seq(1,32,2)) face[8,i,1:4,j] <- FindElements(nodes=nodes[8,j,i,],elements=c(element[2,(30+i)%%32+1],element[3,(30+i)%%32+1]),z=j,nodex=nodex) #BC-
        for(i in 1:32)        face[9,i,1:4,j] <- FindElements(nodes=nodes[9,j,i,],elements=c(element[3,(30+i)%%32+1],element[4,(30+i)%%32+1]),z=j,nodex=nodex) #CD-
        for(i in 1:32)       face[10,i,1:4,j] <- FindElements(nodes=nodes[10,j,i,],elements=c(element[3,(30+i)%%32+1],element[5,(30+i)%%32+1]),z=j,nodex=nodex) #CE-
        for(i in 1:32)       face[11,i,1:4,j] <- FindElements(nodes=nodes[11,j,i,],elements=c(element[4,(30+i)%%32+1],element[6,(30+i)%%32+1]),z=j,nodex=nodex) #DF-
        for(i in 1:32)       face[12,i,1:4,j] <- FindElements(nodes=nodes[12,j,i,],elements=c(element[5,(30+i)%%32+1],element[7,(30+i)%%32+1]),z=j,nodex=nodex) #EG-
    }
    #
    setTxtProgressBar(pb5, 4) 
    if (nmulti==1){
        param <- list(list(exyz0[[1]],exyz1[[1]]))
    } else{
        param <- list(list(exyz0[[1]],exyz1[[1]]),list(exyz0[[2]],exyz1[[2]]))
    }
    tlist <-  NULL
    for(multi in 1:nmulti){
        tlist <-  append(tlist,list(exyz0[[multi]],exyz1[[multi]]))
    }
    #param <- list(tlist)
    #
    now <- Sys.time()
    #if (useRhpc & useCluster){
    #cl2 <- makeCluster(2)#param,height,xcod,ycod,nodes,element,sigmaxy,sigmaz,face
    #res <- parLapply(cl2,param,fun=findBorder,height=height,xcod=xcod,ycod=ycod,nodes=nodes,element=element,sigmaxy=sigmaxy,sigmaz=sigmaz,face=face,iter=iter,est64=est64)
    res <- lapply(param,FUN=findBorder,height=height,xcod=xcod,ycod=ycod,nodes=nodes,element=element,sigmaxy=sigmaxy,sigmaz=sigmaz,face=face,iter=iter,est64=est64,useRhpc=useRhpc,useCluster=useCluster)
    #stopCluster(cl2)
    cat("\n","findborder",Sys.time()-now,"\n")
    for (multi in 1:nmulti) interres[[multi]] <- res[[multi]]$interres
    #
    interres2 <- lapply(interres,FUN=apply,MARGIN=2:4,sum)
    interres3 <- lapply(interres2,FUN=aperm,perm=c(3,2,1))
    interres4 <- lapply(interres3,FUN=apply,MARGIN=3,as.vector) # margin=3
    interresmulti <- do.call(rbind,(interres4))
#########
    setTxtProgressBar(pb5, 5) 
    cat("\n")
    return(list(rsens=interresmulti,boundr=boundr))
}
###################################################################################
#  Conductivity update subroutine
##################### 4000 begin subroutine #######################################
updatecond <- function(X,param,nee,bb,dd,rank,nmulti,nexe,nexx,nexr,nexc,resEnhRatio,nex,numpol,kvol,jcur,maxcur,difp,vol,aniso,hetero,est64,iter,alpha,radChg,dev,Marquardt,sigma0)	 {
    #browser()
    pb4 <- txtProgressBar(min = 1, max = 6, style = 3)
    condEst <- X$cond
    alam <- X$alam
    #sigma0 <- X$sigma0
    sensemat <- param
    sens  <- rep(list(matrix(0,maxcur,nexe)),nmulti)
    sensz <- rep(list(matrix(0,maxcur,nexe-1)),nmulti)
    bx <- vector("numeric", length=ifelse(aniso,2*nex,nex)) # for anisotropy
    #rbx <- vector("numeric", length=32) # for anisotropy
    setTxtProgressBar(pb4, 1) 
    offset1 <- X$offset[1]
    offset2 <- X$offset[2]
    offset3 <- X$offset[3]
    offset4 <- X$offset[4]
    offset5 <- X$offset[5]
    offset6 <- X$offset[6]
    offset7 <- offset8 <- offset9 <- 0
    for (multi in 1:nmulti) {
      if (aniso==TRUE) {
        if (hetero==TRUE) {
          sens[[multi]][numpol[1:kvol,1:jcur],1:nexe] <-                       2*aperm(sensemat[[1]][[multi]][1:nexr[9],1:kvol,1:jcur],perm=c(2,3,1))
          sens[[multi]][numpol[1:kvol,1:jcur],(nexe+1):(2*nexe)] <-            1*aperm(sensemat[[2]][[multi]][1:nexr[9],1:kvol,1:jcur],perm=c(2,3,1))
        } else {
          #sens[[multi]][numpol[1:kvol,1:jcur],1:nexx[2]] <-                  2*aperm(sensemat[[1]][[multi]][1:nexr[2],1:kvol,1:jcur],perm=c(2,3,1))
          sens[[multi]][numpol[1:kvol,1:jcur],          1:nexx[1]]<-2/2*aperm(apply(array(sensemat[[1]][[multi]][      0+(1:(nexr[1]-      0)+offset1-1)%%(nexr[1]  -0      )+1,1:kvol,1:jcur],c(resEnhRatio[1],(nexr[1]-0    )/resEnhRatio[1],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          sens[[multi]][numpol[1:kvol,1:jcur],(nexx[1]+1):nexx[2]]<-2/2*aperm(apply(array(sensemat[[1]][[multi]][nexr[1]+(1:(nexr[2]-nexr[1])+offset2-1)%%(nexr[2]-nexr[1])+1,1:kvol,1:jcur],c(resEnhRatio[2],(nexr[2]-nexr[1])/resEnhRatio[2],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          sens[[multi]][numpol[1:kvol,1:jcur],(nexx[2]+1):nexx[3]]<-2/2*aperm(apply(array(sensemat[[1]][[multi]][nexr[2]+(1:(nexr[3]-nexr[2])+offset3-1)%%(nexr[3]-nexr[2])+1,1:kvol,1:jcur],c(resEnhRatio[3],(nexr[3]-nexr[2])/resEnhRatio[3],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          sens[[multi]][numpol[1:kvol,1:jcur],(nexx[3]+1):nexx[4]]<-2/2*aperm(apply(array(sensemat[[1]][[multi]][nexr[3]+(1:(nexr[4]-nexr[3])+offset4-1)%%(nexr[4]-nexr[3])+1,1:kvol,1:jcur],c(resEnhRatio[4],(nexr[4]-nexr[3])/resEnhRatio[4],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          sens[[multi]][numpol[1:kvol,1:jcur],(nexx[4]+1):nexx[5]]<-2/2*aperm(apply(array(sensemat[[1]][[multi]][nexr[4]+(1:(nexr[5]-nexr[4])+offset5-1)%%(nexr[5]-nexr[4])+1,1:kvol,1:jcur],c(resEnhRatio[5],(nexr[5]-nexr[4])/resEnhRatio[5],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          sens[[multi]][numpol[1:kvol,1:jcur],(nexx[5]+1):nexx[6]]<-2/2*aperm(apply(array(sensemat[[1]][[multi]][nexr[5]+(1:(nexr[6]-nexr[5])+offset6-1)%%(nexr[6]-nexr[5])+1,1:kvol,1:jcur],c(resEnhRatio[6],(nexr[6]-nexr[5])/resEnhRatio[6],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          #
          sens[[multi]][numpol[1:kvol,1:jcur],nexx[6]+1] <- apply(sensemat[[3]][[multi]][(nexr[6]+1):nexr[9],1:kvol,1:jcur],c(2,3),sum) # subcutaneous fat
          # longitudinal direction
          sensz[[multi]][numpol[1:kvol,1:jcur],          1:nexx[1]]<- aperm(apply(array(sensemat[[2]][[multi]][      0+(1:(nexr[1]-      0)+offset1-1)%%(nexr[1]-0      )+1,1:kvol,1:jcur],c(resEnhRatio[1],(nexr[1]-0      )/resEnhRatio[1],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          sensz[[multi]][numpol[1:kvol,1:jcur],(nexx[1]+1):nexx[2]]<- aperm(apply(array(sensemat[[2]][[multi]][nexr[1]+(1:(nexr[2]-nexr[1])+offset2-1)%%(nexr[2]-nexr[1])+1,1:kvol,1:jcur],c(resEnhRatio[2],(nexr[2]-nexr[1])/resEnhRatio[2],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          sensz[[multi]][numpol[1:kvol,1:jcur],(nexx[2]+1):nexx[3]]<- aperm(apply(array(sensemat[[2]][[multi]][nexr[2]+(1:(nexr[3]-nexr[2])+offset3-1)%%(nexr[3]-nexr[2])+1,1:kvol,1:jcur],c(resEnhRatio[3],(nexr[3]-nexr[2])/resEnhRatio[3],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          sensz[[multi]][numpol[1:kvol,1:jcur],(nexx[3]+1):nexx[4]]<- aperm(apply(array(sensemat[[2]][[multi]][nexr[3]+(1:(nexr[4]-nexr[3])+offset4-1)%%(nexr[4]-nexr[3])+1,1:kvol,1:jcur],c(resEnhRatio[4],(nexr[4]-nexr[3])/resEnhRatio[4],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          sensz[[multi]][numpol[1:kvol,1:jcur],(nexx[4]+1):nexx[5]]<- aperm(apply(array(sensemat[[2]][[multi]][nexr[4]+(1:(nexr[5]-nexr[4])+offset5-1)%%(nexr[5]-nexr[4])+1,1:kvol,1:jcur],c(resEnhRatio[5],(nexr[5]-nexr[4])/resEnhRatio[5],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          sensz[[multi]][numpol[1:kvol,1:jcur],(nexx[5]+1):nexx[6]]<- aperm(apply(array(sensemat[[2]][[multi]][nexr[5]+(1:(nexr[6]-nexr[5])+offset6-1)%%(nexr[6]-nexr[5])+1,1:kvol,1:jcur],c(resEnhRatio[6],(nexr[6]-nexr[5])/resEnhRatio[6],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
        }
      } else {
        #if (hetero==TRUE) {
        #  sens[[multi]][numpol[1:kvol,1:jcur],1:nexe] <-      aperm(sensemat[[1]][[multi]][1:nexr[9],1:kvol,1:jcur],perm=c(2,3,1))
        #} else {
        #sens[[multi]][numpol[1:kvol,1:jcur],1:nexx[2]] <-     aperm(sensemat[[1]][[multi]][1:nexr[2],1:kvol,1:jcur],perm=c(2,3,1))
        sens[[multi]][numpol[1:kvol,1:jcur],          1:nexx[1]]<-aperm(apply(array(sensemat[[1]][[multi]][      0+(1:(nexr[1]-      0)+offset1-1)%%(nexr[1]-0      )+1,1:kvol,1:jcur],c(resEnhRatio[1],(nexr[1]-0      )/resEnhRatio[1],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
        sens[[multi]][numpol[1:kvol,1:jcur],(nexx[1]+1):nexx[2]]<-aperm(apply(array(sensemat[[1]][[multi]][nexr[1]+(1:(nexr[2]-nexr[1])+offset2-1)%%(nexr[2]-nexr[1])+1,1:kvol,1:jcur],c(resEnhRatio[2],(nexr[2]-nexr[1])/resEnhRatio[2],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
        sens[[multi]][numpol[1:kvol,1:jcur],(nexx[2]+1):nexx[3]]<-aperm(apply(array(sensemat[[1]][[multi]][nexr[2]+(1:(nexr[3]-nexr[2])+offset3-1)%%(nexr[3]-nexr[2])+1,1:kvol,1:jcur],c(resEnhRatio[3],(nexr[3]-nexr[2])/resEnhRatio[3],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
        sens[[multi]][numpol[1:kvol,1:jcur],(nexx[3]+1):nexx[4]]<-aperm(apply(array(sensemat[[1]][[multi]][nexr[3]+(1:(nexr[4]-nexr[3])+offset4-1)%%(nexr[4]-nexr[3])+1,1:kvol,1:jcur],c(resEnhRatio[4],(nexr[4]-nexr[3])/resEnhRatio[4],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
        sens[[multi]][numpol[1:kvol,1:jcur],(nexx[4]+1):nexx[5]]<-aperm(apply(array(sensemat[[1]][[multi]][nexr[4]+(1:(nexr[5]-nexr[4])+offset5-1)%%(nexr[5]-nexr[4])+1,1:kvol,1:jcur],c(resEnhRatio[5],(nexr[5]-nexr[4])/resEnhRatio[5],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
        sens[[multi]][numpol[1:kvol,1:jcur],(nexx[5]+1):nexx[6]]<-aperm(apply(array(sensemat[[1]][[multi]][nexr[5]+(1:(nexr[6]-nexr[5])+offset6-1)%%(nexr[6]-nexr[5])+1,1:kvol,1:jcur],c(resEnhRatio[6],(nexr[6]-nexr[5])/resEnhRatio[6],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
        #}
        #
        if (hetero) {
          sens[[multi]][numpol[1:kvol,1:jcur],(nexx[6]+1):nexx[7]]<-aperm(apply(array(sensemat[[1]][[multi]][nexr[6]+(1:(nexr[7]-nexr[6])+offset7-1)%%(nexr[7]-nexr[6])+1,1:kvol,1:jcur],c(resEnhRatio[7],(nexr[7]-nexr[6])/resEnhRatio[7],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          sens[[multi]][numpol[1:kvol,1:jcur],(nexx[7]+1):nexx[8]]<-aperm(apply(array(sensemat[[1]][[multi]][nexr[7]+(1:(nexr[8]-nexr[7])+offset8-1)%%(nexr[8]-nexr[7])+1,1:kvol,1:jcur],c(resEnhRatio[8],(nexr[8]-nexr[7])/resEnhRatio[8],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
          sens[[multi]][numpol[1:kvol,1:jcur],(nexx[8]+1):nexx[9]]<-aperm(apply(array(sensemat[[1]][[multi]][nexr[8]+(1:(nexr[9]-nexr[8])+offset9-1)%%(nexr[9]-nexr[8])+1,1:kvol,1:jcur],c(resEnhRatio[9],(nexr[9]-nexr[8])/resEnhRatio[9],kvol,jcur)),c(2,3,4),sum),c(2,3,1))
        } else sens[[multi]][numpol[1:kvol,1:jcur],nexx[6]+1] <-   apply(sensemat[[1]][[multi]][(nexr[6]+1):nexr[9],1:kvol,1:jcur],c(2,3),sum)
      }
    }
    setTxtProgressBar(pb4, 2) 
###########################################################################
    #cat("--- singular value decomposition ---\n")
    #cat("--- START --- \n")
    sensmulti <- do.call(rbind,sens)
    sensmultiz <- do.call(rbind,sensz)
    #browser()
    if (aniso) {
        if(hetero){
          normvec1 <- normvec2 <- normvec3 <- normvec4 <- 1
          sensmulti <- cbind(sensmulti,sensmultiz)
        } else{
        rsens <- sensemat[[4]] 
        sensall <- cbind(sensmulti,sensmultiz,rsens)
        normvec1 <- (norm(sensall[,1:nexx[6]],type="f")                /norm(sensall,type="f")) *185/60
        normvec2 <- (norm(sensall[,(2*nexx[6]+2):185],type="f")        /norm(sensall,type="f")) *185/64
        normvec3 <- (norm(as.matrix(sensall[,nexx[6]+1]),type="f")     /norm(sensall,type="f")) *185/1
        normvec4 <- (norm(sensall[,(nexx[6]+2):(2*nexx[6]+1)],type="f")/norm(sensall,type="f")) *185/60
        #normvec1 <- norm(sensmulti[,-(nexx[6]+1)],type="2")/norm(sensmultiz,type="2")
        #normvec2 <- norm(sensmulti[,-(nexx[6]+1)],type="2")/norm(sensemat[[4]],type="2")
        #normvec2 <- svd(sensmulti[,-(nexx[6]+1)])$d[1]/svd(sensemat[[4]])$d[1]
        #sensmulti <- cbind(sensmulti,diag(rep(normvec1,1024))%*%sensmultiz,diag(rep(normvec2,1024))%*%sensemat[[4]])
        sensmulti <- cbind(sensmulti[,-(nexx[6]+1)]/normvec1,sensmulti[,nexx[6]+1]/normvec3,sensmultiz/normvec4,rsens/normvec2)
        }
    } else {
        if(hetero){
        normvec1 <- normvec2 <- normvec3 <- 1
        } else{
          rsens <- sensemat[[2]]
          sensall <- cbind(sensmulti,rsens)
          normvec1 <- (norm(sensall[,1:nexx[6]],type="f")            /norm(sensall,type="f"))*(sum(nexc[1:6])+64+1)/sum(nexc[1:6])
          normvec2 <- (norm(sensall[,(nexx[6]+2):(sum(nexc[1:6])+64+1)],type="f")      /norm(sensall,type="f"))*(sum(nexc[1:6])+64+1)/64
          normvec3 <- (norm(as.matrix(sensall[,nexx[6]+1]),type="f")/norm(sensall,type="f"))*(sum(nexc[1:6])+64+1)/1
          #normvec2 <- norm(sensmulti,type="2")/norm(sensemat[[2]],type="2")
          #normvec2 <- svd(sensmulti[,-(nexx[6]+1)])$d[1]/svd(sensemat[[2]])$d[1]
          #sensmulti <- cbind(sensmulti[,-(nexx[6]+1)],sensmulti[,nexx[6]+1],diag(rep(normvec2,1024))%*%rsens)
          if(condEst){
            sensmulti <- cbind(sensmulti[,-(nexx[6]+1)]/normvec1,sensmulti[,nexx[6]+1]/normvec3) 
          }else{
            sensmulti <- rsens/normvec2
          }
        }
      
    }
    #cat("normvec2, normvec3;",normvec2,normvec3,"\n")
    #cat("size of sensmulti", nrow(sensmulti),ncol(sensmulti), "\n")
    #browser()
    #if(alpha == 1) normvec <- rep(5,ncol(sensmulti))
    #else           normvec <- rep(1,ncol(sensmulti))
    #normvec <- c(rep(1,ncol(sensmulti)-64),rep(normvec2,64))
    #else  normvec <- apply(sensmulti,2,FUN=function(x) norm(matrix(x),type="F"))
    #normvec <- apply(sensmulti,2,FUN=function(x) ((norm(matrix(x),type="F"))))
    #
    #sensmulti <- sensmulti%*%diag(normvec) #normalization
    ret <- svd(sensmulti) ###
    rank <- length(ret$d[ret$d > 1e-16])  #rankMatrix(sensmulti)
    senss <- ret$u
    dd <- ret$d
    bb <- ret$v
    #cat("size of bb", dim(bb), "\n")
    #cat("minimum of dd", dd[length(dd)], "\n")
    setTxtProgressBar(pb4, 3) 
    #cat("--- FINISHED ---- \n")
###########################################################################
    k1 <- ((alam/20) <= (dd^2)) ##
    tmat <- (rep(1,length(dd)))%o%dd
    tmat[upper.tri(tmat)] <- 0
    k2 <- (apply(tmat,1,sum)/sum(dd) < 0.995)
    nee1 <- max(length(k1[k1==TRUE]),2)
    #nee2 <- max(length(k2[k2==TRUE]),2)
    #dd <- dd[dd > 0.20]  # 0.24]
    #cat("maximum of dd", dd[1],"\n")
    #cat("minimum of dd",dd[nee],"\n")
    #cat("check point 3\n")
    #cat("nee:", nee,"\n")
    #cat("check point 4\n")
######## sensivility matrix #####################
    difpmulti <- do.call(c,as.vector(difp))
    #cat("size of difpmulti", length(difpmulti), "\n")
    #cat("size of senss", dim(senss),"\n")
############## Tikonove method #############################################
    bx0 <- -t(senss)%*%(difpmulti)
    #browser()
    bxt <- bx0*dd/(dd^2 + alam^2)
    #bxt <- bx0*dd/(dd^2)
    #normvec3 <- c(rep(1,ncol(sensmulti)-64),rep(100*normvec2,64))
    #bb <- bb%*%diag(1/normvec3) #normalization
############## Marquardt method ############################################
    if (Marquardt){
      r0 <- senss%*%bb
      cat("check point 5\n")
      ss <- t(r0)%*%(r0)
      bx0 <- -t(senss)%*%(difpmulti)
      bxt <- backsolve(ss+alam/20*diag(norm(sensmulti[,-(nexx[6]+1)],type="2"),nrow(ss)),x=bx0)
    }
############################################################################
    #nee <- min(nee1, length(dd[dd > 0.2]))
    bx1 <- bxt
#   bxt <- bxt[1:nee]
#    bb <- bb[,1:nee]
############################################################################
    setTxtProgressBar(pb4, 4) 
    if(condEst){
      if(aniso){
        if (hetero==TRUE){
          bx[1:nexr[9]] <- 1/2*bb[1:nexr[9],]%*%bxt
          bx[(nexr[9]+1):(nexr[9]+nexr[9])] <- bb[(nexr[9]+1):(nexr[9]+nexr[9]),]%*%bxt
        } else {
          #bx[1:(nexr[2])] <- 1/2*bb[1:nexx[2],]%*%bxt
          bx[          0+(1:(nexr[1]-      0)+offset1-1)%%(nexr[1]-      0)+1] <- 2/2*rep(bb[(      0+1):nexx[1],]%*%bxt,each=resEnhRatio[1])/normvec1
          bx[    nexr[1]+(1:(nexr[2]-nexr[1])+offset2-1)%%(nexr[2]-nexr[1])+1] <- 2/2*rep(bb[(nexx[1]+1):nexx[2],]%*%bxt,each=resEnhRatio[2])/normvec1
          bx[    nexr[2]+(1:(nexr[3]-nexr[2])+offset3-1)%%(nexr[3]-nexr[2])+1] <- 2/2*rep(bb[(nexx[2]+1):nexx[3],]%*%bxt,each=resEnhRatio[3])/normvec1
          bx[    nexr[3]+(1:(nexr[4]-nexr[3])+offset4-1)%%(nexr[4]-nexr[3])+1] <- 2/2*rep(bb[(nexx[3]+1):nexx[4],]%*%bxt,each=resEnhRatio[4])/normvec1
          bx[    nexr[4]+(1:(nexr[5]-nexr[4])+offset5-1)%%(nexr[5]-nexr[4])+1] <- 2/2*rep(bb[(nexx[4]+1):nexx[5],]%*%bxt,each=resEnhRatio[5])/normvec1
          bx[    nexr[5]+(1:(nexr[6]-nexr[5])+offset6-1)%%(nexr[6]-nexr[5])+1] <- 2/2*rep(bb[(nexx[5]+1):nexx[6],]%*%bxt,each=resEnhRatio[6])/normvec1
          bx[   (nexr[6]+1):nexr[9]] <- 1/normvec3*2/2*bb[nexx[6]+1,]%*%bxt
          #
          #bx[nex+1:(nexr[2])] <- bb[nexe + 1:nexx[2],]%*%bxt
          bx[nexr[9]+      0+(1:(nexr[1]-      0)+offset1-1)%%(nexr[1]-      0)+1] <- rep(bb[nexe + (      0+1):nexx[1],]%*%bxt,each=resEnhRatio[1])/normvec4
          bx[nexr[9]+nexr[1]+(1:(nexr[2]-nexr[1])+offset2-1)%%(nexr[2]-nexr[1])+1] <- rep(bb[nexe + (nexx[1]+1):nexx[2],]%*%bxt,each=resEnhRatio[2])/normvec4
          bx[nexr[9]+nexr[2]+(1:(nexr[3]-nexr[2])+offset3-1)%%(nexr[3]-nexr[2])+1] <- rep(bb[nexe + (nexx[2]+1):nexx[3],]%*%bxt,each=resEnhRatio[3])/normvec4
          bx[nexr[9]+nexr[3]+(1:(nexr[4]-nexr[3])+offset4-1)%%(nexr[4]-nexr[3])+1] <- rep(bb[nexe + (nexx[3]+1):nexx[4],]%*%bxt,each=resEnhRatio[4])/normvec4
          bx[nexr[9]+nexr[4]+(1:(nexr[5]-nexr[4])+offset5-1)%%(nexr[5]-nexr[4])+1] <- rep(bb[nexe + (nexx[4]+1):nexx[5],]%*%bxt,each=resEnhRatio[5])/normvec4
          bx[nexr[9]+nexr[5]+(1:(nexr[6]-nexr[5])+offset6-1)%%(nexr[6]-nexr[5])+1] <- rep(bb[nexe + (nexx[5]+1):nexx[6],]%*%bxt,each=resEnhRatio[6])/normvec4
          bx[nexr[9]+(nexr[6]+1):nexr[9]] <- 1/normvec3*bb[nexx[6]+1,]%*%bxt
          #
          #rbx <- bb[(2*nexe-1)+(1:ifelse(est64,64,32)),]%*%diag(rep(1/normvec2*1024/64,length(bxt)))%*%bxt
          #rbx <- bb[(2*nexe-1)+(1:ifelse(est64,64,32)),]%*%diag(rep(1/normvec2,length(bxt)))%*%bxt
          #rbx <- rbx/normvec2 ############
        }
      } else {
          #if (hetero==TRUE){
          #  bx[1:nex] <- bb[1:nex,]%*%bxt
          #} else {
          #bx[1:(nexr[2])] <- bb[1:nexx[2],]%*%bxt
          #browser()
          bx[      0+(1:(nexr[1]-      0)+offset1-1)%%(nexr[1]-      0)+1] <- rep(bb[(      0+1):nexx[1],]%*%bxt,each=resEnhRatio[1])/normvec1
          bx[nexr[1]+(1:(nexr[2]-nexr[1])+offset2-1)%%(nexr[2]-nexr[1])+1] <- rep(bb[(nexx[1]+1):nexx[2],]%*%bxt,each=resEnhRatio[2])/normvec1
          bx[nexr[2]+(1:(nexr[3]-nexr[2])+offset3-1)%%(nexr[3]-nexr[2])+1] <- rep(bb[(nexx[2]+1):nexx[3],]%*%bxt,each=resEnhRatio[3])/normvec1
          bx[nexr[3]+(1:(nexr[4]-nexr[3])+offset4-1)%%(nexr[4]-nexr[3])+1] <- rep(bb[(nexx[3]+1):nexx[4],]%*%bxt,each=resEnhRatio[4])/normvec1
          bx[nexr[4]+(1:(nexr[5]-nexr[4])+offset5-1)%%(nexr[5]-nexr[4])+1] <- rep(bb[(nexx[4]+1):nexx[5],]%*%bxt,each=resEnhRatio[5])/normvec1
          bx[nexr[5]+(1:(nexr[6]-nexr[5])+offset6-1)%%(nexr[6]-nexr[5])+1] <- rep(bb[(nexx[5]+1):nexx[6],]%*%bxt,each=resEnhRatio[6])/normvec1
          #}
          if (hetero) {
            bx[nexr[6]+(1:(nexr[7]-nexr[6])+offset7-1)%%(nexr[7]-nexr[6])+1] <- rep(bb[(nexx[6]+1):nexx[7],]%*%bxt,each=resEnhRatio[7])
            bx[nexr[7]+(1:(nexr[8]-nexr[7])+offset9-1)%%(nexr[8]-nexr[7])+1] <- rep(bb[(nexx[7]+1):nexx[8],]%*%bxt,each=resEnhRatio[8])
            bx[nexr[8]+(1:(nexr[9]-nexr[8])+offset9-1)%%(nexr[9]-nexr[8])+1] <- rep(bb[(nexx[8]+1):nexx[9],]%*%bxt,each=resEnhRatio[9])
            #rbx <- NULL
          } else { 
            bx[(nexr[6]+1):nexr[9]] <- 1/normvec3*bb[nexx[6]+1,]%*%bxt
            #bx[(nexr[6]+1):nexr[9]] <- bb[nexx[6]+1,]%*%bxt
            #rbx <- bb[nexe+(1:ifelse(est64,64,32)),]%*%bxt
          }
        #
      }
      rbx <- NULL
    } else {
            rbx <- bb[1:ifelse(est64,64,32),]%*%diag(rep(1/normvec2,length(bxt)))%*%bxt
    }
    setTxtProgressBar(pb4, 5) 
    #browser()
      #cMatrix <- matrix(0,nrow=(61+64),ncol=(708+64))
      #cMatrix[      1:(nexx[1]-1),      1:(nexr[1])] <- 1/resEnhRatio[1]
      #cMatrix[(nexx[1]+1):nexx[2],(nexr[1]+1):nexr[2]] <- 1/resEnhRatio[2]
      #cMatrix[(nexx[2]+1):nexx[3],(nexr[2]+1):nexr[3]] <- 1/resEnhRatio[3]
      #cMatrix[(nexx[3]+1):nexx[4],(nexr[3]+1):nexr[4]] <- 1/resEnhRatio[4]
      #cMatrix[(nexx[4]+1):nexx[5],(nexr[4]+1):nexr[5]] <- 1/resEnhRatio[5]
      #cMatrix[(nexx[5]+1):nexx[6],(nexr[5]+1):nexr[6]] <- 1/resEnhRatio[6]
      #cMatrix[(nexx[6]+1):(nexx[6]+1),(nexr[5]+1):nexr[9]] <- 1/resEnhRatio[7]/resEnhRatio[8]/resEnhRatio[9]
      #cMatrix[(nexx[6]+1+1:64),(nexr[9]+1:64)] <- 1
      #A <- sensmulti
      #x <- bb%*%(bx0*dd/(dd^2 + alam))
      #
    
    ############################ Conductivity update ######################################################################################
    if(condEst) {
      if (aniso==TRUE) {
        sigmaxy <- sigma0xy*exp(bx[1:nex])
        sigmaz  <- sigma0z*exp(bx[nex+1:nex]) # for anisotropy
        ######################################################################
        #sigmaxy[which(sigmaxy < sigmaxy[nex])[(which(sigmaxy < sigmaxy[nex])>nexr[5] & which(sigmaxy < sigmaxy[nex]) <=nexr[6])]] <- sigmaxy[nex]*1.5
        #sigmaz[which(sigmaz < sigmaz[nex])[(which(sigmaz < sigmaz[nex])>nexr[5] & which(sigmaz < sigmaz[nex]) <=nexr[6])]] <- sigmaz[nex]*1.5
        ######################################################################
        sigma <- sqrt((2*sigmaxy^2+sigmaz^2)/3)
      } else if (aniso==FALSE) {
        sigmaxy <- sigmaz <- sigma <- sigma0[1:nex]* exp(bx[1:nex])
        ##############################################################################
        #sigma[which(sigma < 1.1*sigma[nex])[(which(sigma < 1.1*sigma[nex])>nexr[5] & which(sigma < 1.1*sigma[nex]) <=nexr[6])]] <- sigma[nex]*1.1
        #sigmaxy[1:nex] <- sigmaz[1:nex] <- sigma[1:nex]
        ##############################################################################
      }
      normbx <- norm(as.matrix(na.omit(bx[1:nex])),type="f")
    } else {
      sigma <- sigmaxy <- sigmaz <- NULL
      normbx <- NULL
    }
#######################################################################
#######################################################################
  #cat("size of bx", length(bx), "\n")
  setTxtProgressBar(pb4, 6) 
  cat("\n")
  return(list(bx=bx,dd=dd,rank=rank,nee=nee,bb=bb,rbx=rbx,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,normbx=normbx))
}
#####################################################
#################################################################################
# SUBROUTINE systema
#################################################################################
#
systemaR <- function(numnp=nnodezx,mband=nbw,a=aa0,mxbnd=mxbnd,mxnode=mxnode){
  #browser()
  for (n in 1:(numnp-1)){
    for (l in 2:mband){
      c <- a[n,l] / a[(n-1)%%mxnode+1,n%/%mxnode+1]
      i <- n + l - 1
      if (i > numnp) next()
      #j <- 0
      #for (k in l:mband){
      k <- l:mband
       # j <- j + 1
        j <- k - l + 1
        a[i,j] <- a[i,j]-c*a[n,k]
      #}
      a[n,l] <- c
    }
    #cat("n=",n,"/",numnp,"\n")
  }
  #browser()
  return(a)
}
#
#
systemaRGPU <- function(numnp=nnodezx,mband=nbw,a=aa0,mxbnd=mxbnd,mxnode=mxnode){
  browser()
  gpua <- gpuMatrix(a,type="float")
  for (n in 1:(numnp-1)){
    for (l in 2:mband){
      c <- gpua[n,l] / gpua[(n-1)%%mxnode+1,n%/%mxnode+1]
      i <- n + l - 1
      if (i > numnp) next()
      #j <- 0
      for (k in l:mband){
        #k <- l:mband
        #gpuk <- gpuVector(k, type = NULL)
        # j <- j + 1
        j <- k - l + 1
        #gpuj <- gpuVector(j, type = "integer")
        gpua[i,j] <- gpua[i,j]-c*gpua[n,k]
      }
      gpua[n,l] <- c
      #
      #gpua <- vclMatrix(a)
    }
    cat("n=",n,"/",numnp,"\n")
  }
  #browser()
  a <- gpua[]
  return(a)
}
#################################################################################
# SUBROUTINE setupaaPAfort
#################################################################################
setupaaPAfort <- function(sigmaxy,sigmaz,nodex,nodez,st,stxy,stz,nnodez,nbw,nez,nex,mxbnd,mxnode,mze,mxe){
	#nr <- nrow(aal)
	#nc <- ncol(aal)
	aal <- matrix(0,nrow=mxnode,ncol=mxbnd)
	out <- .Fortran("setaaPA",aal=as.double(as.vector(aal)),laal=as.integer(length(as.vector(aal))),
			sigmaxy=as.double(as.vector(sigmaxy)),lsigmaxy=as.integer(length(as.vector(sigmaxy))),
			sigmaz=as.double(as.vector(sigmaz)),lsigmaz=as.integer(length(as.vector(sigmaz))),
			nodex=as.integer(as.vector(nodex)),lnodex=as.integer(length(as.vector(nodex))),
			nodez=as.integer(as.vector(nodez)),lnodez=as.integer(length(as.vector(nodez))),
			st=as.double(as.vector(stxy)),lst=as.integer(length(as.vector(stxy))),
			stz=as.double(as.vector(stz)),lstz=as.integer(length(as.vector(stz))),
			nnodez=as.integer(nnodez),nbw=as.integer(nbw),nez=as.integer(nez),nex=as.integer(nex),
			mxnode=as.integer(mxnode)
			)
	ret <- matrix(as.double(as.vector(out$aal)),mxnode,mxbnd)
	return(ret)
}
#
################
#  function sumup
################
sumup <- function(x) {
	y <- sum(x^2,na.rm=TRUE)
	return(y)
}
sumup2 <- function(x) {
  y <- sum(sqrt(x^2),na.rm=TRUE)
  return(y)
}
##############################################
# fieldR
##############################################
fieldR <- function(nex, nez, cex, bb, nnodez, nodex, nodez) {
  browser()
  exyzl1 <- matrix(0, nez, nex)
  exyzl2 <- matrix(0, nez, nex)
  exyzl3 <- matrix(0, nez, nex)
  for (n in 1:nex) {
    for (l in 1:nez) {
      #for (j in 1:3) {
        #exyzl[j, l, n] <- -sum(cex[j, 1:4, l, n] * (bb[nnodez * (nodex[1:4, l, n] - 1) + nodez[1:4, l, n]]))
      exyzl1[l, n] <- -sum(cex[1, 1:4, l, n] * (bb[nnodez * (nodex[1:4, l, n] - 1) + nodez[1:4, l, n]]))
      exyzl2[l, n] <- -sum(cex[2, 1:4, l, n] * (bb[nnodez * (nodex[1:4, l, n] - 1) + nodez[1:4, l, n]]))
      exyzl3[l, n] <- -sum(cex[3, 1:4, l, n] * (bb[nnodez * (nodex[1:4, l, n] - 1) + nodez[1:4, l, n]]))
      #}
    }
  }
  exyzl <- array(c(as.vector(exyzl1),as.vactor(exyzl2),as.vector(exyzl3)),dim=c(3,nez,nex))
  return(exyzl)
}
#
##########################################################################
# fwdjunPA
##########################################################################
fwdjunPA <- function(X,iinn,iott,nnodezx,ibnd,nnodez,nb,vb,oneb2,oneb8,oneb4,bnd,nbw,bb,aa0,mxbnd,mxnode,nex,nez,nodex,nodez,mze,mxe,vol,pod,jcur,kvol,cex,nzot,nzin,difp,podex,nnodex,acqmode) {
    param <- unlist(X)
    ncur  <- param[1]
    npat  <- param[2]
    pods <- vector("numeric",length=kvol)
    pots <- vector("numeric",length=kvol)
    vbs <- vector("numeric",length=nb)
    vb0s <- vector("numeric",length=nb)
    iin <- iinn[ncur,npat]
    iot <- iott[ncur,npat]
    #cat("iin,iot",iin,iot," ")
    #browser()
    #
############# BB ##########################################################
    bb <- vector(mode="numeric",length=nnodezx) ###
    bb[1:nnodezx] <- 0
    #
    if (npat==1) {
      iinm1 <- (ibnd[(iin-2)%%nb+1]-1)*nnodez + nzin
      iinc <-  (ibnd[iin          ]-1)*nnodez + nzin
      iinp1 <- (ibnd[(iin)%%nb+1  ]-1)*nnodez + nzin
      iotp1 <- (ibnd[(iot)%%nb+1  ]-1)*nnodez + nzin
      iotc <-  (ibnd[iot          ]-1)*nnodez + nzin
      iotm1 <- (ibnd[(iot-2)%%nb+1]-1)*nnodez + nzin
    } else if(npat==2){
      iinm1 <- (ibnd[(iin-2)%%nb+1]-1)*nnodez + nzot
      iinc <-  (ibnd[iin          ]-1)*nnodez + nzot
      iinp1 <- (ibnd[(iin)%%nb+1  ]-1)*nnodez + nzot
      iotp1 <- (ibnd[(iot)%%nb+1  ]-1)*nnodez + nzot
      iotc <-  (ibnd[iot          ]-1)*nnodez + nzot
      iotm1 <- (ibnd[(iot-2)%%nb+1]-1)*nnodez + nzot
    } else if(npat==3){
      iinm1 <- (ibnd[(iin-2)%%nb+1]-1)*nnodez + nzin
      iinc <-  (ibnd[iin          ]-1)*nnodez + nzin
      iinp1 <- (ibnd[(iin)%%nb+1  ]-1)*nnodez + nzin
      iotp1 <- (ibnd[(iot)%%nb+1  ]-1)*nnodez + nzot
      iotc <-  (ibnd[iot          ]-1)*nnodez + nzot
      iotm1 <- (ibnd[(iot-2)%%nb+1]-1)*nnodez + nzot
    }else if(npat==4){
      iinm1 <- (ibnd[(iin-2)%%nb+1]-1)*nnodez + nzot
      iinc <-  (ibnd[iin          ]-1)*nnodez + nzot
      iinp1 <- (ibnd[(iin)%%nb+1  ]-1)*nnodez + nzot
      iotp1 <- (ibnd[(iot)%%nb+1  ]-1)*nnodez + nzin
      iotc <-  (ibnd[iot          ]-1)*nnodez + nzin
      iotm1 <- (ibnd[(iot-2)%%nb+1]-1)*nnodez + nzin
    }
    #
    bb[iinc] <- oneb2
    bb[iotc] <- -oneb2
    if(iin%%2 == 1) {
        bb[iinc+1] <- oneb8
        bb[iinc-1] <- oneb8
        bb[iinm1] <- oneb4*bnd[iin]/(bnd[iin]+bnd[iin+1])
        bb[iinp1] <- oneb4*bnd[iin+1]/(bnd[iin]+bnd[iin+1])
    } else {
        stop
    }
    if(iot%%2 == 1) {
        bb[iotc+1] <- -oneb8
        bb[iotc-1] <- -oneb8
        bb[iotm1] <- -oneb4*bnd[iot]/(bnd[iot]+bnd[iot+1])
        bb[iotp1] <- -oneb4*bnd[iot+1]/(bnd[iot]+bnd[iot+1])
    } else {
        stop
    }
      #cat("-JUN4-\n")
    bb <- suppressWarnings(systembfort(numnp=nnodezx,mband=nbw,b=bb,a=aa0,mxbnd=mxbnd,mxnode=mxnode))
      #cat("-JUN5-\n")
    summ <- mean(bb[1:nnodezx])
    bb[1:nnodezx] <- bb[1:nnodezx]-summ
########### FIELD ########################################################
    exyz0s <- suppressWarnings(fieldfort(nex=nex,nez=nez,cex=cex,bb=bb,nnodez=nnodez,nodex=nodex,nodez=nodez,mxnode=mxnode,mze=mze,mxe=mxe)[1:3,1:nez,1:nex])
    #exyz0s <- fieldR(nex=nex,nez=nez,cex=cex,bb=bb,nnodez=nnodez,nodex=nodex,nodez=nodez)[1:3,1:nez,1:nex]
    #
    #
############## BOUNDARY POTENTIAL############
    difps <- NULL
    difp0s <- difp[,ncur]
    if(npat==1){
        summ <- 0
        for(j in 1:nb) {
            vb0s[1:nb] <- vb[1:nb,ncur]
            # depending on acquisitin mode; START
            if (acqmode=="parallel") {
                vbs[j] <- bb[nnodez*(ibnd[j]-1)+nzot]
            }	else if ( (ibnd[j])==(ibnd[iin]) | (ibnd[j])==(ibnd[iot]) )	{
                vbs[j] <- bb[nnodez*(ibnd[j]-1)+nzot]
            }	else if ( (ibnd[j])==(ibnd[iin]-1) | (ibnd[j])==(ibnd[iin]+1) | (ibnd[j])==(ibnd[iot]-1) | (ibnd[j])==(ibnd[iot]+1) ) {
                vbs[j] <- bb[nnodez*(ibnd[j]-1)+nzot]
            }	else if ( (ibnd[j])==(ibnd[iin]-2) | (ibnd[j])==(ibnd[iin]+2) | (ibnd[j])==(ibnd[iot]-2) | (ibnd[j])==(ibnd[iot]+2) ) {
                vbs[j] <- bb[nnodez*(ibnd[j]-1)+(nzot+nzin)/2]
            }	else vbs[j] <- bb[nnodez*(ibnd[j]-1)+nzin]
        }
                                        #browser()
        summ <- mean(vbs[1:nb])
        vb0s[nb+1] <- vb0s[1]
        vbs[nb+1] <- vbs[1]
        vbs[1:(nb+1)] <- vbs[1:(nb+1)]-summ
                                        #cat(iin,vb[iin:nb,ncur],vb[1:iin,ncur],"\n",file=xx) #############
                                        #cat(iin,vb[iin:nb,ncur],vb[1:iin,ncur],"\n") #############
                                        #exyz0s <- exyz[1:3,1:nez,1:nex]
                                        #difpsraw <- NULL
        for (j in 1:kvol) {
            pods[j] <- (vbs[iinn[j,2]]-vbs[iott[j,2]]) ######### modified!!!
            if(is.nan(log(pods[j]/podex[j,ncur]))){
                cat("Bad data ! j,ncur,pod[j,ncur],podex[j,ncur]",j,ncur, pod[j,ncur],podex[j,ncur],"\n")
                podex[j,ncur] <- (podex[((j-1)+kvol-1)%%kvol+1,ncur]+podex[((j+1)+kvol-1)%%kvol+1,ncur])/2
                if(is.nan(log(pod[j,ncur]/podex[j,ncur]))){
                    podex[j,ncur] <- (podex[((j-1)+kvol-1)%%kvol+1,((ncur-1)+32-1)%%32+1]+podex[((j+1)+kvol-1)%%kvol+1,((ncur+1)+32-1)%%32+1])/2
                }
                cat("inter polated ! j,ncur,pod[j,ncur],podex[j,ncur]",j,ncur, pod[j,ncur],podex[j,ncur],"\n")
            }
            difps <- c(difps,log(pods[j]/podex[j,ncur]))
        }
                                        #browser()
#########################################################################
                                        #exyouts <- exyz0s[1:2,3*(ifelse(mode=="parallel",nzot,nzin))+1,1:nex]
        for (n in 1:nnodex) pots[n] <- bb[nnodez*(n-1) + ifelse((acqmode!="parallel" & iin!=(n) & iot!=(n)),nzin,nzot)]
    }
    return(list(difps=difps,difp0s=difp0s,pods=pods,exyz0s=exyz0s,pots=pots,vbs=vbs,vb0s=vb0s))
}
##########################################################################
#########################################################################

###########################################################################
#
#    SUBROUTINE JUN
#
###########################################################################
#junPA <- function(npat,iter,aa0,mxbnd,mxnode,maxcur,numpol,sigma,sigmaxy,sigmaz,iinn,iott,pod,podex,vb,difp,difp0,nodex,nodez,nnodex,nnodez,nnodezx,st,stxy,stz,ibnd,bnd,vol,cex,nbw,nb,nex,nez,nzin,nzot,jcur,kvol,mxe,mze,mxb,acqmode,useCluster,nCPU,useRhpc=FALSE) {
junPA <- function(npat,iter,mxbnd,mxnode,maxcur,numpol,sigma,sigmaxy,sigmaz,iinn,iott,pod,podex,vb,difp,difp0,nodex,nodez,nnodex,nnodez,nnodezx,st,stxy,stz,ibnd,bnd,vol,cex,nbw,nb,nex,nez,nzin,nzot,jcur,kvol,mxe,mze,mxb,acqmode,useCluster,nCPU,useRhpc=FALSE) {
    bb <- vector(mode="numeric",length=mxnode)
    exyz0l <- exyz1l <- array(0,c(3,mze,mxe,jcur,kvol))
    vbl <-   array(0, c(mxb+1,jcur))
    vb0l <-  array(0, c(mxb+1,jcur))
    oneb2 <- 1/2
    oneb4 <- 1/4
    oneb8 <- 1/8
    oneb16 <- 1/16
    oneb32 <- 1/32
    oneb64 <- 1/64
###########################################################################
    #cat("npat,,iter\n",npat,iter,"\n")
    pb2 <- txtProgressBar(min = 1, max = 5, style = 3)
    setTxtProgressBar(pb2, 1) 
    #cat("-JUN1-")
    #browser()
    #aa0[1:mxnode,1:mxbnd] <- 0
    #now <-Sys.time()
    #aa0 <- setupaaPA(    aal=aa0,sigmaxy=sigmaxy,sigmaz=sigmaz,nodex=nodex,nodez=nodez,st=st,stxy=stxy,stz=stz,nnodez=nnodez,nbw=nbw,nez=nez,nex=nex)
    #cat("setupaaPA",difftime(Sys.time(),now,units="sec"),"\n")
    now <-Sys.time()
    aa0 <- setupaaPAfort(sigmaxy=sigmaxy,sigmaz=sigmaz,nodex=nodex,nodez=nodez,st=st,stxy=stxy,stz=stz,nnodez=nnodez,nbw=nbw,nez=nez,nex=nex,mxbnd=mxbnd,mxnode=mxnode,mze=mze,mxe=mxe)
    #setupaaPAfort <- function(aal,sigmaxy,sigmaz,nodex,nodez,st,stxy,stz,nnodez,nbw,nez,nex,mxbnd,mxnode,mze,mxe){
    cat("setupaaPAfort",difftime(Sys.time(),now,units="sec"),"\n")
    #
    setTxtProgressBar(pb2, 2) 
    #cat("-JUN2-")
    #browser()
    now <-Sys.time()
    aa0 <- suppressWarnings(systemafort(numnp=nnodezx,mband=nbw,a=aa0,mxbnd=mxbnd,mxnode=mxnode))
    cat("systemafort",difftime(Sys.time(),now,units="sec"),"\n")
    #now <-Sys.time()
    #aa0 <- systemaRGPU(numnp=nnodezx,mband=nbw,a=aa0,mxbnd=mxbnd,mxnode=mxnode)
    cat("systemaGPU",difftime(Sys.time(),now,units="sec"),"\n")
    #
    #browser()
    setTxtProgressBar(pb2, 3) 
    #cat("-JUN3-")
########### Current poles ################################################
############ loop ########################################################
    z <- apply(matrix(c(rep(1:jcur, 4), rep(1:4, each = jcur)), 4 * jcur, 
                      2),1,FUN=function(x) return(list(x)))
    #browser()
    #Rhpc_initialize()
    #cl<-Rhpc_getHandle(4)
    #ignore <- Rhpc_EvalQ(cl=cl, expr={library(AbdominalEITextra,quietly=TRUE); NULL},envir=as.environment(globalenv()))
    #
    gc();    gc()
    now <- Sys.time()
    if (useRhpc & useCluster){
      #Rhpc::Rhpc_initialize()
      #cl2 <-Rhpc::Rhpc_getHandle(nCPU)
      ignore <- Rhpc::Rhpc_EvalQ(cl=cl2, expr={library(AbdominalEITextra,quietly=TRUE); NULL},envir=as.environment(globalenv()))
      tm <- proc.time()  
        ret <- parLapply2(cl2, x = z[1:jcur], fun = fwdjunPA, iinn = iinn, iott = iott, nnodezx = nnodezx, ibnd = ibnd, nnodez = nnodez, nb = nb, 
                         vb = vb, oneb2 = oneb2, oneb8 = oneb8, oneb4 = oneb4, bnd = bnd, nbw = nbw, bb = bb, aa0 = aa0, mxbnd = mxbnd, mxnode = mxnode, 
                         nex = nex, nez = nez, nodex = nodex, nodez = nodez, mze = mze, mxe = mxe, vol = vol, pod = pod, jcur = jcur, kvol = kvol, 
                         cex = cex, nzot = nzot, nzin = nzin, difp = difp, podex = podex, nnodex = nnodex, acqmode = acqmode) 
      print(proc.time() - tm)
      #Rhpc::Rhpc_finalize()
    } else {
      #
      #library(AbdominalEITextra)
      cl <- parallel::makeCluster(parallel::detectCores(),type="PSOCK")
      #cl <- snow::makeCluster(ifelse(useCluster,nCPU,detectCores()),type=ifelse(useCluster,"MPI","SOCK"))
      ignore <- parallel::clusterEvalQ(cl=cl, expr={library(AbdominalEITextra,quietly=TRUE); NULL})
      tm <- snow.time(  
        ret <- parallel::parLapply(cl, X = z[1:jcur], fun = fwdjunPA, iinn = iinn, iott = iott, nnodezx = nnodezx, ibnd = ibnd, nnodez = nnodez, nb = nb, 
                         vb = vb, oneb2 = oneb2, oneb8 = oneb8, oneb4 = oneb4, bnd = bnd, nbw = nbw, bb = bb, aa0 = aa0, mxbnd = mxbnd, mxnode = mxnode, 
                         nex = nex, nez = nez, nodex = nodex, nodez = nodez, mze = mze, mxe = mxe, vol = vol, pod = pod, jcur = jcur, kvol = kvol, 
                         cex = cex, nzot = nzot, nzin = nzin, difp = difp, podex = podex, nnodex = nnodex, acqmode = acqmode) 
      )
    #
      #stopCluster(cl)
      dev.set(5)
      par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(4,4,2,1))
      plot(tm)
    }
    # 
      for (i in 1:jcur) {
        difp[,i] <- ret[[i]]$difps
        difp0[,i] <- ret[[i]]$difp0s
                                        #difpraw[,i] <- ret[[i]]$difpsraw
        pod[,i] <- ret[[i]]$pods
        exyz0l[,,,i,] <- ret[[i]]$exyz0s
        vbl[,i] <- ret[[i]]$vbs
        vb0l[,i] <- ret[[i]]$vb0s
    }
    #
    #
    if(useRhpc & useCluster){
    tm <- proc.time()  
        ret <- parLapply2(cl2, x = z[(jcur+1):(2*jcur)], fun = fwdjunPA, iinn = iinn, iott = iott, nnodezx = nnodezx, ibnd = ibnd, nnodez = nnodez, nb = nb, 
                          vb = vb, oneb2 = oneb2, oneb8 = oneb8, oneb4 = oneb4, bnd = bnd, nbw = nbw, bb = bb, aa0 = aa0, mxbnd = mxbnd, mxnode = mxnode, 
                          nex = nex, nez = nez, nodex = nodex, nodez = nodez, mze = mze, mxe = mxe, vol = vol, pod = pod, jcur = jcur, kvol = kvol, 
                          cex = cex, nzot = nzot, nzin = nzin, difp = difp, podex = podex, nnodex = nnodex, acqmode = acqmode) 
    print(proc.time()-tm)
    } else {
    #browser()
    tm <- snow.time(  
        ret <- parallel::parLapply(cl,X = z[(jcur+1):(2*jcur)], fun = fwdjunPA, iinn = iinn, iott = iott, nnodezx = nnodezx, ibnd = ibnd, nnodez = nnodez, nb = nb, 
                           vb = vb, oneb2 = oneb2, oneb8 = oneb8, oneb4 = oneb4, bnd = bnd, nbw = nbw, bb = bb, aa0 = aa0, mxbnd = mxbnd, mxnode = mxnode, 
                           nex = nex, nez = nez, nodex = nodex, nodez = nodez, mze = mze, mxe = mxe, vol = vol, pod = pod, jcur = jcur, kvol = kvol, 
                           cex = cex, nzot = nzot, nzin = nzin, difp = difp, podex = podex, nnodex = nnodex, acqmode = acqmode) 
    )
      #
    dev.set(6)
    par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(4,4,2,1))
    plot(tm)
    }
      #
    for (i in 1:jcur) {
        exyz1l[,,,i,] <- ret[[i]]$exyz0s
    }
    if (acqmode=="bypass") {
      exyz2l <- exyz3l <- array(0,c(3,mze,mxe,jcur,kvol))
      tm <- system.time(  
      ret <- parLapply(cl, X = z[(2*jcur+1):(3*jcur)], fun = fwdjunPA, iinn = iinn, iott = iott, nnodezx = nnodezx, ibnd = ibnd, nnodez = nnodez, nb = nb, 
                         vb = vb, oneb2 = oneb2, oneb8 = oneb8, oneb4 = oneb4, bnd = bnd, nbw = nbw, bb = bb, aa0 = aa0, mxbnd = mxbnd, mxnode = mxnode, 
                         nex = nex, nez = nez, nodex = nodex, nodez = nodez, mze = mze, mxe = mxe, vol = vol, pod = pod, jcur = jcur, kvol = kvol, 
                         cex = cex, nzot = nzot, nzin = nzin, difp = difp, podex = podex, nnodex = nnodex, acqmode = acqmode) 
      )
      #
      print(tm)
      #
      for (i in 1:jcur) {
        exyz2l[,,,i,] <- ret[[i]]$exyz0s
      }
      tm <- system.time(  
      ret <- parLapply(cl, X = z[(3*jcur+1):(4*jcur)], fun = fwdjunPA, iinn = iinn, iott = iott, nnodezx = nnodezx, ibnd = ibnd, nnodez = nnodez, nb = nb, 
                         vb = vb, oneb2 = oneb2, oneb8 = oneb8, oneb4 = oneb4, bnd = bnd, nbw = nbw, bb = bb, aa0 = aa0, mxbnd = mxbnd, mxnode = mxnode, 
                         nex = nex, nez = nez, nodex = nodex, nodez = nodez, mze = mze, mxe = mxe, vol = vol, pod = pod, jcur = jcur, kvol = kvol, 
                         cex = cex, nzot = nzot, nzin = nzin, difp = difp, podex = podex, nnodex = nnodex, acqmode = acqmode) 
      )
      #
      print(tm)
      #
      for (i in 1:jcur) {
        exyz3l[,,,i,] <- ret[[i]]$exyz0s
      }
    } else {
      exyz2l <- exyz3l <- NULL
    }
                                        #browser()
                                        #exyout <- ret[[1]]$exyouts #??? no need
                                        #pot  <- ret[[1]]$pots #??? no need
                                        #browser()
                                        #pod <- pod[c(32,1:31),]  ###
    setTxtProgressBar(pb2, 4)
    if(useRhpc & useCluster){
      #Rhpc::Rhpc_finalize()
    } else {
      stopCluster(cl)
    }
    setTxtProgressBar(pb2, 5) 
    cat("\n")
    cat("\n foward prblem=",Sys.time()-now,"\n")
    dev.set(12)
    par(mfrow=c(4,1),oma=c(0,0,0,0),mar=c(5,4,1,1))
    plot(x=1:32,y=podex[1:32,1],type="b",pch=20,xlab="Electrode",ylab="Voltage (mV)",sub=paste("current electrode 1,"),col="blue",ylim=c(-200,200))
    legend("topright",c("Experimental","Caluculated"),pch=c(20,10),bty="n",col=c("blue","red"))
    points(x=1:32,y=pod[1:32,1],type="b",pch=10,col="red")
    plot(x=1:32,y=podex[1:32,9],type="b",pch=20,xlab="Electrode",ylab="Voltage (mV)",sub=paste("current electrode 9,"),col="blue",ylim=c(-200,200))
    points(x=1:32,y=pod[1:32,9],type="b",pch=10,col="red")
    plot(x=1:32,y=podex[1:32,17],type="b",pch=20,xlab="Electrode",ylab="Voltage (mV)",sub=paste("current electrode 17,"),col="blue",ylim=c(-200,200))
    points(x=1:32,y=pod[1:32,17],type="b",pch=10,col="red")
    plot(x=1:32,y=podex[1:32,25],type="b",pch=20,xlab="Electrode",ylab="Voltage (mV)",sub=paste("current electrode 25,"),col="blue",ylim=c(-200,200))
    points(x=1:32,y=pod[1:32,25],type="b",pch=10,col="red")
    #return(list(status=0,aa=aa0,difp=difp,difp0=difp0,pod=pod,exyz0=exyz0l,exyz1=exyz1l,exyz2=exyz2l, exyz3=exyz3l, vb=vbl,vb0=vb0l))
    #return(list(status=0,aa=aa0,difp=difp,difp0=difp0,pod=pod,exyz0=exyz0l,exyz1=exyz1l,exyz2=exyz2l, exyz3=exyz3l, vb=vbl,vb0=vb0l))
    return(list(status=0,difp=difp,difp0=difp0,pod=pod,exyz0=exyz0l,exyz1=exyz1l,exyz2=exyz2l, exyz3=exyz3l, vb=vbl,vb0=vb0l))
}
#
#########################################################################################################################################################
# Sensitivity calculation
#########################################################################################################################################################
senseCalc <- function(mxnode,numpol,sigma,sigmaxy,sigmaz,pod,jcur,kvol,mxe,mze,nex,nez,vol,exyz0,exyz1) {
    sensematl   <- array(0,c(mxe,jcur,kvol))
    sensematxyl <- array(0,c(mxe,jcur,kvol))
    sensematzl  <- array(0,c(mxe,jcur,kvol))
#
    pb3 <- txtProgressBar(min = 1, max = 4, style = 3)
    setTxtProgressBar(pb3, 1) 
    #cl <- makeCluster(detectCores())
    #pb <- txtProgressBar(min = 1, max = jcur, style = 3)
    #foreach (j = 1:jcur) %do% {
    #    for (k in 1:kvol) {
    #        sensematl[1:nex,j,k] <- -sigma[1:nex]/pod[k,j]*apply(2*vol[1:nez,1:nex]*(exyz1[1,1:nez,1:nex,k]*exyz0[1,1:nez,1:nex,j]+exyz1[2,1:nez,1:nex,k]*exyz0[2,1:nez,1:nex,j]+exyz1[3,1:nez,1:nex,k]*exyz0[3,1:nez,1:nex,j]),MARGIN=2,FUN=sum)
    #        sensematxyl[1:nex,j,k] <- -sigma[1:nex]/pod[k,j]*apply(2*vol[1:nez,1:nex]*(exyz1[1,1:nez,1:nex,k]*exyz0[1,1:nez,1:nex,j]+exyz1[2,1:nez,1:nex,k]*exyz0[2,1:nez,1:nex,j]),MARGIN=2,FUN=sum)
    #        sensematzl[1:nex,j,k] <- -sigma[1:nex]/pod[k,j]*apply(2*vol[1:nez,1:nex]*(exyz1[3,1:nez,1:nex,k]*exyz0[3,1:nez,1:nex,j]),MARGIN=2,FUN=sum)
    #    }
    #    setTxtProgressBar(pb, j)
    #}
    #cat("\n")
    #cl4 <- makeCluster(detectCores())
    #browser()
    registerDoParallel(detectCores())
    sensematl[1:mxe,1:jcur,1:kvol]   <- foreach(j=1:jcur, .combine=rbind) %do% invcurPA(ncur=j,nex=nex,nez=nez,vol=vol,exyz0=exyz0,exyz1=exyz1,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,pod=pod,mxe=mxe,mze=mze,mxnode=mxnode,jcur=jcur,kvol=kvol)
    setTxtProgressBar(pb3, 2)
    sensematxyl[1:mxe,1:jcur,1:kvol] <- foreach(j=1:jcur, .combine=rbind) %do% invcurPAxy(ncur=j,nex=nex,nez=nez,vol=vol,exyz1=exyz1,exyz0=exyz0,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,pod=pod,mxe=mxe,mze=mze,mxnode=mxnode,jcur=jcur,kvol=kvol)
    setTxtProgressBar(pb3, 3)
    sensematzl[1:mxe,1:jcur,1:kvol]  <- foreach(j=1:jcur, .combine=rbind) %do% invcurPAz(ncur=j,nex=nex,nez=nez,vol=vol,exyz1=exyz1,exyz0=exyz0,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,pod=pod,mxe=mxe,mze=mze,mxnode=mxnode,jcur=jcur,kvol=kvol)
    #stopCluster(cl4)     
    stopImplicitCluster()
    #senseret <- apply(array(1:jcur,jcur),1,invcurPA,nex=nex,nez=nez,vol=vol,exyz1=exyz1,exyz0=exyz0,sigma=sigma,sigmaxy=sigmaxy,sigmaz=sigmaz,pod=pod,mxe=mxe,mze=mze,mxnode=mxnode,jcur=jcur,kvol=kvol)
                                        #browser()
    #sensematl[1:mxe,1:jcur,1:kvol] <- senseret
    #sensematxyl[1:mxe,1:jcur,1:kvol] <- senseretxy
    #sensematzl[1:mxe,1:jcur,1:kvol]  <- senseretz
    #
    setTxtProgressBar(pb3, 4)
    cat("\n")
    return(list(status=0,sensemat=aperm(sensematl,c(1,2,3)),sensematxy=aperm(sensematxyl,c(1,3,2)),sensematz=aperm(sensematzl,c(1,3,2))))
}
################################################################################################
#
# FUNCTION invcurPA
#
invcurPA <- function(ncur,nex,nez,vol,exyz0,exyz1,sigma,sigmaxy,sigmaz,pod,mxe,mze,mxnode,jcur,kvol) {
    sensematl <- array(0,c(mxe,jcur))
    #browser()
    for (j in 1:kvol) {
        sensematl[1:nex,j] <- -sigma[1:nex]/pod[ncur,j]*apply(2*vol[1:nez,1:nex]*(exyz1[1,1:nez,1:nex,j,ncur]*exyz0[1,1:nez,1:nex,ncur,j]+exyz1[2,1:nez,1:nex,j,ncur]*exyz0[2,1:nez,1:nex,ncur,j]+exyz1[3,1:nez,1:nex,j,ncur]*exyz0[3,1:nez,1:nex,ncur,j]),MARGIN=2,FUN=sum)
     	  #sensematxyl[1:nex,j] <- -sigmaxy[1:nex]/pod[ncur,j]*apply(2*vol[1:nez,1:nex]*(exyz1[1,1:nez,1:nex,j]*exyz0[1,1:nez,1:nex,ncur]+exyz1[2,1:nez,1:nex,j]*exyz0[2,1:nez,1:nex,ncur]),MARGIN=2,FUN=sum)
       	#sensematzl[1:nex,j] <- -sigmaz[1:nex]/pod[ncur,j]*apply(2*vol[1:nez,1:nex]*(exyz1[3,1:nez,1:nex,j]*exyz0[3,1:nez,1:nex,ncur]),MARGIN=2,FUN=sum)
      	#fieldmatxyl[1:nex,j] <- colSums(vol[1:nez,1:nex]*(exyz[1,1:nez,1:nex]*exyz0[1,1:nez,1:nex,j]+exyz[2,1:nez,1:nex]*exyz0[2,1:nez,1:nex,j]))
        #fieldmatzl[1:nex,j]  <- colSums(vol[1:nez,1:nex]*(exyz[3,1:nez,1:nex]*exyz0[3,1:nez,1:nex,j]))
    }
    return(senses=sensematl) #,fieldsxy=fieldmatxyl,fieldsz=fieldmatzl))
}
invcurPAxy <- function(ncur,nex,nez,vol,exyz1,exyz0,sigma,sigmaxy,sigmaz,pod,mxe,mze,mxnode,jcur,kvol) {
    #sensematl <- array(0,c(mxe,jcur))
    sensematxyl <- array(0,c(mxe,jcur))
    #sensematzl <- array(0,c(mxe,jcur))
    #fieldmatxyl <- array(0,c(mxe,jcur))
    #fieldmatzl <- array(0,c(mxe,jcur))
    for (j in 1:kvol) {
        #sensematl[1:nex,j] <- -sigma[1:nex]/pod[ncur,j]*apply(2*vol[1:nez,1:nex]*(exyz1[1,1:nez,1:nex,j]*exyz0[1,1:nez,1:nex,ncur]+exyz1[2,1:nez,1:nex,j]*exyz0[2,1:nez,1:nex,ncur]+exyz1[3,1:nez,1:nex,j]*exyz0[3,1:nez,1:nex,ncur]),MARGIN=2,FUN=sum)
     	sensematxyl[1:nex,j] <- -sigmaxy[1:nex]/pod[ncur,j]*apply(2*vol[1:nez,1:nex]*(exyz1[1,1:nez,1:nex,j,ncur]*exyz0[1,1:nez,1:nex,ncur,j]+exyz1[2,1:nez,1:nex,j,ncur]*exyz0[2,1:nez,1:nex,ncur,j]),MARGIN=2,FUN=sum)
       	#sensematzl[1:nex,j] <- -sigmaz[1:nex]/pod[ncur,j]*apply(2*vol[1:nez,1:nex]*(exyz1[3,1:nez,1:nex,j]*exyz0[3,1:nez,1:nex,ncur]),MARGIN=2,FUN=sum)
      	#fieldmatxyl[1:nex,j] <- colSums(vol[1:nez,1:nex]*(exyz[1,1:nez,1:nex]*exyz0[1,1:nez,1:nex,j]+exyz[2,1:nez,1:nex]*exyz0[2,1:nez,1:nex,j]))
        #fieldmatzl[1:nex,j]  <- colSums(vol[1:nez,1:nex]*(exyz[3,1:nez,1:nex]*exyz0[3,1:nez,1:nex,j]))
    }
    return(sensesxy=sensematxyl) #,fieldsxy=fieldmatxyl,fieldsz=fieldmatzl))
}
invcurPAz <- function(ncur,nex,nez,vol,exyz1,exyz0,sigma,sigmaxy,sigmaz,pod,mxe,mze,mxnode,jcur,kvol) {
    #sensematl <- array(0,c(mxe,jcur))
    #sensematxyl <- array(0,c(mxe,jcur))
    sensematzl <- array(0,c(mxe,jcur))
    #fieldmatxyl <- array(0,c(mxe,jcur))
    #fieldmatzl <- array(0,c(mxe,jcur))
    for (j in 1:kvol) {
        #sensematl[1:nex,j] <- -sigma[1:nex]/pod[ncur,j]*apply(2*vol[1:nez,1:nex]*(exyz1[1,1:nez,1:nex,j]*exyz0[1,1:nez,1:nex,ncur]+exyz1[2,1:nez,1:nex,j]*exyz0[2,1:nez,1:nex,ncur]+exyz1[3,1:nez,1:nex,j]*exyz0[3,1:nez,1:nex,ncur]),MARGIN=2,FUN=sum)
     	#sensematxyl[1:nex,j] <- -sigmaxy[1:nex]/pod[ncur,j]*apply(2*vol[1:nez,1:nex]*(exyz1[1,1:nez,1:nex,j]*exyz0[1,1:nez,1:nex,ncur]+exyz1[2,1:nez,1:nex,j]*exyz0[2,1:nez,1:nex,ncur]),MARGIN=2,FUN=sum)
       	sensematzl[1:nex,j] <- -sigmaz[1:nex]/pod[ncur,j]*apply(2*vol[1:nez,1:nex]*(exyz1[3,1:nez,1:nex,j,ncur]*exyz0[3,1:nez,1:nex,ncur,j]),MARGIN=2,FUN=sum)
      	#fieldmatxyl[1:nex,j] <- colSums(vol[1:nez,1:nex]*(exyz[1,1:nez,1:nex]*exyz0[1,1:nez,1:nex,j]+exyz[2,1:nez,1:nex]*exyz0[2,1:nez,1:nex,j]))
        #fieldmatzl[1:nex,j]  <- colSums(vol[1:nez,1:nex]*(exyz[3,1:nez,1:nex]*exyz0[3,1:nez,1:nex,j]))
    }
    return(sensesz=sensematzl) #,fieldsxy=fieldmatxyl,fieldsz=fieldmatzl))
}
#
adjLayer <- function(layer,sensemat,nexr,nexc,rad64,arar,resEnhRatio){
  #
  #browser()
  radAvg <- mean(rad64[1:64])
  lambda <- 0.25+radAvg/2
  #sensAvg <- apply(sensemat,1,norm,type="F")
  sensAvg <- apply(sensemat,1,sum)
  sensavg<- c(-sum(sensAvg[1:nexr[1]]),-sum(sensAvg[(nexr[1]+1):nexr[2]]),-sum(sensAvg[(nexr[2]+1):nexr[3]]),-sum(sensAvg[(nexr[3]+1):nexr[4]]),-sum(sensAvg[(nexr[4]+1):nexr[5]]),-sum(sensAvg[(nexr[5]+1):nexr[6]])) #/resEnhRatio[1:6]
  #sensavg<- c(-mean(sensAvg[1:nexr[1]]),-mean(sensAvg[(nexr[1]+1):nexr[2]]),-mean(sensAvg[(nexr[2]+1):nexr[3]]),-mean(sensAvg[(nexr[3]+1):nexr[4]]),-mean(sensAvg[(nexr[4]+1):nexr[5]]),-mean(sensAvg[(nexr[5]+1):nexr[6]]))*nexc[1:6]
  sSensavg <- sensavg[1:5]/c(layer[1],(layer[2]-layer[1]),(layer[3]-layer[2]),(layer[4]-layer[3]),(lambda-layer[4]))
  design <- sSensavg%o%sSensavg
  design <- -5*design/(c(3,6,7,16,16)%o%c(3,6,7,16,16))
  diag(design) <- -4*diag(design)
  design <- rbind(design,rep(1,5))
  design <- cbind(design,c(rep(-25/2,5),0))
  ans <- solve(design,c(rep(0,5),lambda))
  layerNew <- c(ans[1],sum(ans[1:2]),sum(ans[1:3]),sum(ans[1:4]))
  return(layerNew)
}
#
adjLayer2 <- function(layer,sensemat,nexr,rad64,arar,boost=TRUE){
  #browser()
  radAvg <- mean(rad64[1:64])
  rho <- radAvg
  sensAvg <- apply(sensemat,1,sum)
  #sensavg<- c(-mean(sensAvg[1:nexr[1]]),-mean(sensAvg[(nexr[1]+1):nexr[2]]),-mean(sensAvg[(nexr[2]+1):nexr[3]]),-mean(sensAvg[(nexr[3]+1):nexr[4]]),-mean(sensAvg[(nexr[4]+1):nexr[5]]),-mean(sensAvg[(nexr[5]+1):nexr[6]]))
  sensavg<- c(-sum(sensAvg[1:nexr[1]]),-sum(sensAvg[(nexr[1]+1):nexr[2]]),-sum(sensAvg[(nexr[2]+1):nexr[3]]),-sum(sensAvg[(nexr[3]+1):nexr[4]]),-sum(sensAvg[(nexr[4]+1):nexr[5]]),-sum(sensAvg[(nexr[5]+1):nexr[6]]))
  #sensavg<- sensavg/c(sum(arar[1:nexr[1]]),sum(arar[(nexr[1]+1):nexr[2]]),sum(arar[(nexr[2]+1):nexr[3]]),sum(arar[(nexr[3]+1):nexr[4]]),sum(arar[(nexr[4]+1):nexr[5]]),sum(arar[(nexr[5]+1):nexr[6]]))  ###
  sSensavg <- sensavg[1:6]/c(layer[1],(layer[2]-layer[1]),(layer[3]-layer[2]),(layer[4]-layer[3]),(layer[5]-layer[4]),(rho-layer[5]))
  design <- sSensavg%o%sSensavg
  if (boost) weight <- c(3,6,7,16,16,24)
  else  weight <- c(3,6,7,16,16,12)
  design <- -5*design/(weight%o%weight)
  diag(design) <- -5*diag(design)
  design <- rbind(design,rep(1,6))
  design <- cbind(design,c(rep(-18,6),0))
  ans <- solve(design,c(rep(0,6),rho))
  layerNew <- c(ans[1],sum(ans[1:2]),sum(ans[1:3]),sum(ans[1:4]),sum(ans[1:5]))
  return(layerNew)
}
#
# 
parLapply <- function(cl, x, fun, ...){
  clusterCall(cl, LB.init, fun, ...)
  r <- snow::clusterApplyLB(cl, x, fun=LB.worker)
  clusterEvalQ(cl, rm('.LB.fun', '.LB.args', pos=globalenv()))
  r
}
parLapply2 <- function(cl, x, fun, ...){
  Rhpc::Rhpc_worker_call(cl, LB.init, fun, ...)
  r <- Rhpc::Rhpc_lapplyLB(cl, x, LB.worker)
  Rhpc::Rhpc_EvalQ(cl, rm('.LB.fun', '.LB.args', pos=globalenv()),envir=as.environment(.GlobalEnv))
  r
}
#
parSapply2 <- function(cl, x, fun, ...){
  Rhpc::Rhpc_worker_call(cl, LB.init, fun, ...)
  r <- Rhpc::Rhpc_sapplyLB(cl, x, LB.worker)
  Rhpc::Rhpc_EvalQ(cl, rm('.LB.fun', '.LB.args', pos=globalenv()),envir=as.environment(.GlobalEnv))
  r
}
#
LB.init <- function(fun, ...){
  assign('.LB.fun', fun, pos=globalenv())
  assign('.LB.args', list(...), pos=globalenv())
  NULL
}
#
LB.worker <- function(x){
  do.call('.LB.fun', c(list(x), .LB.args))
}
#
ginv2 <- function (X, tol = sqrt(.Machine$double.eps), damper) {
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Xsvd <- svd(X)
  #if (is.complex(X)) 
  #  Xsvd$u <- Conj(Xsvd$u)
  #Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  #browser()
  #if (all(Positive)){ 
    #invXsvd <- 1/Xsvd$d
    Xsvd$d[Xsvd$d < 1e-12] <- 0
    invXsvd <- Xsvd$d/(Xsvd$d^2 + damper^2)
    Xsvd$v %*% (invXsvd * t(Xsvd$u))
  #}
  #else if (!any(Positive)) 
  #  array(0, dim(X)[2L:1L])
  #else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
  #                                             t(Xsvd$u[, Positive, drop = FALSE]))
}
#
curveture <- function(fx,fxx) {
  r <- ((1+fx^2)^(3/2))/(fxx)
  return(1/r)
}
#
#################################################################################################
#      END OF PROGRAM
#################################################################################################