############################################## 
invsigmafun <- function(x,sigma1,sigma2){
  s1 <- sigma1
  s2 <- sigma2
  z <- (x-s2)*(s1+2*x)/((s1-x)*(s2+2*x)-(s2-x)*(s1+2*x))
  return(z)
}
#
#
#
fatEstimation <- function(dd, fn = "Normal", initial, estres = NULL) {
    #browser()
    if (fn == "Beta") {
      npara <- 2
      # initial <- c(0.05, 0.1,0.5,0.1,0.5)
      llimits <- c(1,4, 0.08)
      ulimits <- c(6, 8, 0.6)
      dfn1 <- function(z, x) dbeta(invsigmafun(z,max(dd)+0.0001,min(dd)-0.0001), shape1 = x[1], shape2 = x[2], log = FALSE)
      dfn2 <- function(z, x) dbeta(invsigmafun(z,max(dd)+0.0001,min(dd)-0.0001), shape1 = x[2], shape2 = x[1], log = FALSE)
    } else if (fn == "Normal") {
        npara <- 2
        # initial <- c(0.05, 0.1,0.5,0.1,0.5)
        llimits <- c(0.02, 0.01, 0.4, 0.03, 0.08)
        ulimits <- c(0.35, 0.05, 0.8, 0.2, 0.8)
        dfn1 <- function(z, x) dnorm(z, mean = x[1], sd = x[2], log = FALSE)/(pnorm(max(dd), 
            mean = x[1], sd = x[2], log = FALSE) - pnorm(min(dd), mean = x[1], sd = x[2], 
            log = FALSE))
        dfn2 <- function(z, x) dnorm(z, mean = x[1], sd = x[2], log = FALSE)/(pnorm(max(dd), 
            mean = x[1], sd = x[2], log = FALSE) - pnorm(min(dd), mean = x[1], sd = x[2], 
            log = FALSE))
    } else if (fn == "skewNormal") {
        npara <- 3
        # initial <- c(0.05, 0.10, 0, 0.5, 0.10, 0, 0.5)
        llimits <- c(0, 0.05, -5, 0.4, 0.05, -5, 0.1)
        ulimits <- c(0.4, 0.15, 5, 1, 0.15, 5, 0.9)
        dfn1 <- function(z, x) dsn(z, dp = x[1:npara], log = FALSE)/(psn(max(dd), 
            dp = x[1:npara], log = FALSE) - psn(min(dd), dp = x[1:npara], log = FALSE))
        dfn2 <- function(z, x) dsn(z, dp = x[1:npara], log = FALSE)/(psn(max(dd), 
            dp = x[1:npara], log = FALSE) - psn(min(dd), dp = x[1:npara], log = FALSE))
    } else if (fn == "logNormal") {
        npara <- 2
        # initial <- c(-3, 0.2, -0.1, 0.2, 0.5)
        llimits <- c(-5, 0.05, -1, 0.05, 0.1)
        ulimits <- c(-0.8, 0.8, 1, 0.8, 0.9)
        # dfn1 <- function(z,x)
        # (dlnorm((z-min(dd)),meanlog=x[1],sdlog=x[2],log=FALSE))/(plnorm((max(dd)-min(dd)),meanlog=x[1],sdlog=x[2],log=FALSE))
        # dfn2 <- function(z,x)
        # (dlnorm((max(dd)-z),meanlog=x[1],sdlog=x[2],log=FALSE))/(plnorm((max(dd)-min(dd)),meanlog=x[1],sdlog=x[2],log=FALSE))
        dfn1 <- function(z, x) dlnorm(z, meanlog = x[1], sdlog = x[2], log = FALSE)/(plnorm(max(dd), 
            meanlog = x[1], sdlog = x[2], log = FALSE) - plnorm(min(dd), meanlog = x[1], 
            sdlog = x[2], log = FALSE))
        dfn2 <- function(z, x) dlnorm(z, meanlog = x[1], sdlog = x[2], log = FALSE)/(plnorm(max(dd), 
            meanlog = x[1], sdlog = x[2], log = FALSE) - plnorm(min(dd), meanlog = x[1], 
            sdlog = x[2], log = FALSE))
    } else if (fn == "Rayleigh") {
        npara <- 1
        initial <- c(1, 1, 0.5)
        llimits <- c(0.1, 0.1, 0.1)
        ulimits <- c(0.5, 0.5, 0.9)
        dfn1 <- function(z, x) (drayleigh((z - min(dd)), scale = x[1], log = FALSE))/(prayleigh(max(dd) - 
            min(dd), scale = x[1]))
        dfn2 <- function(z, x) (drayleigh((max(dd) - z), scale = x[1], log = FALSE))/(prayleigh(max(dd) - 
            min(dd), scale = x[1]))
    } else if (fn == "Nakagami") {
        npara <- 2
        llimits <- c(1, 1, 1, 1, 0.1)
        initial <- c(2, 5, 2, 5, 0.5)
        ulimits <- c(5, 9, 5, 9, 0.9)
        dfn1 <- function(z, x) (dnaka((z - min(dd)), shape = x[1], scale = x[2], 
            log = FALSE))/(pnaka((max(dd) - min(dd)), shape = x[1], scale = x[2]))
        dfn2 <- function(z, x) (dnaka((max(dd) - z), shape = x[1], scale = x[2], 
            log = FALSE))/(pnaka((max(dd) - min(dd)), shape = x[1], scale = x[2]))
    } else if (fn == "Rice") {
        npara <- 2
        initial <- c(20, 10, 20, 10, 0.5)
        llimits <- c(0, 0, 0, 0, 0.2)
        ulimits <- c(10, 10, 10.1, 10, 0.9)
        pars <- c(1, 1, 1, 1, 1)
        dfn1 <- function(z, x) (drice((z - min(dd)), vee = x[1], sigma = x[2], log = FALSE))/(price(q = (max(dd) - 
            min(dd)), vee = x[1], sigma = x[2]))
        dfn2 <- function(z, x) (drice((max(dd) - z), vee = x[1], sigma = x[2], log = FALSE))/(price(q = (max(dd) - 
            min(dd)), vee = x[1], sigma = x[2]))
        dfn3 <- function(z, x) (drice((z - min(dd)), vee = x[1], sigma = x[2], log = FALSE))/(price(q = (max(dd) - 
            min(dd)), vee = x[1], sigma = x[2]))
        dfn4 <- function(z, x) (drice((max(dd) - z), vee = x[1], sigma = x[2], log = FALSE))/(price(q = (max(dd) - 
            min(dd)), vee = x[1], sigma = x[2]))
    }
    price <- function(q, vee, sigma = 1) {
      if (!is.Numeric(q)) 
      stop("bad input for argument 'q'")
      if (!is.Numeric(vee)) 
      stop("bad input for argument 'vee'")
      if (!is.Numeric(sigma)) 
      stop("bad input for argument 'sigma'")
    integrate(f = drice, lower = 0, upper = q, vee = vee, sigma = sigma, stop.on.error = FALSE, 
              subdivisions = 150)$value
    }
    #browser()
    mle <- function(x, z) {
        ll <- log(x[2 * npara + 1] * dfn1(z, x[1:npara]) + (1 - x[2 * npara + 1]) * 
            dfn2(z, x[(npara + 1):(2 * npara)]))
        return(ll)
    }
    mle2 <- function(x, z) {
      ll <- log(x[npara + 1] * dfn1(z, x[1:npara]) + (1 - x[npara + 1]) * dfn2(z, x[1:npara]))
      return(ll)
    }
    LL <- function(x) {
        ml <- mle(x[1:(2 * npara + 1)], dd)
        cumll <- sum(ml, na.rm = FALSE)
        return(cumll)
    }
    LL2 <- function(x) {
      ml <- mle2(x[1:(npara + 1)], dd)
      cumll <- sum(ml, na.rm = TRUE)
      return(cumll)
    }
    
    estres0 <- estres
    if (fn=="Beta") {
      pars <- c(rep(0.1, npara), 0.1)
      e <- try(estres <- optim(par = initial, fn = LL2, method = "L-BFGS-B", lower = llimits, 
        upper = ulimits, control = list(fnscale = -1, trace = 2, maxit = 150, parscale = pars)))
    } else{
      pars <- c(rep(0.1, npara * 2), 1)
      e <- try(estres <- optim(par = initial, fn = LL, method = "L-BFGS-B", lower = llimits, 
        upper = ulimits, control = list(fnscale = -1, trace = 1, maxit = 150, parscale = pars)), silent = TRUE)
    }
    if (class(e) == "try-error")  estres <- estres0
    estres$par
    estres$value
    estres$convergence
    # dev.off() dev.new(title='Conductivity
    # distribution',width=4,height=3,xpos=90,ypos=320)
    # par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(5,1,0,1))
    # dd2 <- log10(dd)
    truehist(dd, prob = TRUE, bty = "n", lwd = 2, cex = 1.2, xlab = "", 
        ylab = "Density", xlim = c((0.05), (1.5)), axes = FALSE, main = "", cex.sub = 1,log="x")
    mtext(side = 1, adj = NA, padj = 3.5, cex = 1, "Conductivity (S/m)")
    # mtext(side=1,adj=NA,padj=5,cex=1.0,expression(lambda))
    # mtext(side=1,adj=NA,padj=5.4,cex=1.0,paste(' , ',round(res$par[2*npara+1],2),';
    # Fat area,',round(bound$crossarea*(res$par[2*npara+1])/100,0),sep=''))
    # mtext(side=1,adj=1,padj=3.9,cex=1.0,expression(paste(cm^2)))
    axis(side = 1, labels = (c(0.1,0.2,0.3,0.5,1)), at = (c(0.1,0.2,0.3,0.5,1)))
    if(!is.null(estres)){
    if (fn=="Beta") {
    curve(     estres$par[npara + 1]  * dfn1(x, x = estres$par[1:npara]), from = min(dd), to = max(dd), col = "red", add = TRUE, lty = 2, lwd = 2)
    curve((1 - estres$par[npara + 1]) * dfn2(x, x = estres$par[1:npara]), from = min(dd), to = max(dd), col = "green", add = TRUE, lty = 6, lwd = 2)
    curve(estres$par[npara + 1] * dfn1(x, x = estres$par[1:npara]) + (1 - estres$par[npara + 1]) * dfn2(x, x = estres$par[1:npara]), from = min(dd), 
          to = max(dd), col = "black", add = TRUE, lty = 1, lwd = 2)
    } else{
      curve(estres$par[2 * npara + 1] * dfn1(x, x = estres$par[1:npara]), from = min(dd), 
            to = max(dd), col = "red", add = TRUE, lty = 2, lwd = 2)
      curve((1 - estres$par[2 * npara + 1]) * dfn2(x, x = estres$par[(npara + 1):(2 * npara)]), from = min(dd), to = max(dd), col = "green", add = TRUE, lty = 6, 
            lwd = 2)
      curve(estres$par[2 * npara + 1] * dfn1(x, x = estres$par[1:npara]) + (1 - estres$par[2 * npara + 1]) * dfn2(x, x = estres$par[(npara + 1):(2 * npara)]), from = min(dd), 
            to = max(dd), col = "black", add = TRUE, lty = 1, lwd = 2)
    }
    }
    # legend('topleft',legend=c('Fat component','Non-fat
    # component'),lty=c(2,6),bty='n',col=c('red','green'),lwd=2)
    # dev.copy2eps(file=paste(unlist(strsplit(filename,'\\.'))[1],fn,'.eps',sep=''))
    return(estres)
} 
