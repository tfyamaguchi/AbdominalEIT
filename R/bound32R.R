bound32ROLD <- function(filename) {
    kei <- as.vector(t(read.table(filename, header = FALSE)))[1:32]
    # xx <- seq( 0+360/64, 360+360/64, length.out = 33 )[1:32]
    xx <- seq(0, 360, length.out = 33)[1:32]
    pispl <- periodicSpline(xx, kei, period = 360, ord = 2)
    boundary <- predict(pispl, seq(0, 360, length.out = 65))
    bound32 <- cbind(1:65, boundary$y * sin(boundary$x/180 * pi), -1 * boundary$y * 
        cos(boundary$x/180 * pi), boundary$x)
    colnames(bound32) <- c("sn", "x", "y", "deg")
    bound1 <- bound32
    bound2 <- as.data.frame(bound1[-1, ])
    bound2 <- cbind(bound2, radius = sqrt(bound2$x^2 + bound2$y^2))
    rang <- 360/64 * pi/180
    area <- sum(bound2$radius^2/2 * rang)
    circ <- sum(bound2$radius * rang)
    plot(bound2$x, bound2$y, xlim = c(-150, 150), ylim = c(-150, 150))
    vecb <- matrix(NA, 32, 2)
    vecb[1:32, ] <- c(bound2$x[(1:32 - 1) * 2 + 2] - bound2$x[((1:32 - 1 - 1) * 2)%%64 + 
        2], bound2$y[(1:32 - 1) * 2 + 2] - bound2$y[((1:32 - 1 - 1) * 2)%%64 + 2])
    for (i in 1:32) {
        vecb[i, ] <- c(vecb[i, 1]/sqrt(vecb[i, 1]^2 + vecb[i, 2]^2), vecb[i, 2]/sqrt(vecb[i, 
            1]^2 + vecb[i, 2]^2))
    }
    veca <- cbind(bound2$x[2 * (1:32) - 1], bound2$y[2 * (1:32) - 1])
    angle <- function(x) x[1:2] %*% x[3:4]/sqrt(x[1:2] %*% x[1:2])/sqrt(x[3:4] %*% 
        x[3:4])
    ct <- apply(cbind(veca, vecb), 1, angle)
    st <- sqrt(1 - ct^2)
    theta <- acos(ct) * 180/pi
    sph <- 4
    text(bound2$x[seq(1, 64, 2)] - 10, bound2$y[seq(1, 64, 2)] + 10, round(theta, 
        0))
    text(bound2$x[seq(1, 64, 2)] - 20, bound2$y[seq(1, 64, 2)] + 20, round(sph * 
        ct, 1))
    bound3 <- bound2[seq(1, 64, 2), ]
    nx <- NULL
    ny <- NULL
    for (i in 1:32) {
        nx <- c(nx, bound3[i, "x"] + ct[i] * sph * vecb[i, 1])
        ny <- c(ny, bound3[i, "y"] + ct[i] * sph * vecb[i, 2])
    }
    bound3 <- cbind(bound3, nx, ny)
    points(bound3$nx, bound3$ny, col = "red")
    # pispl <-
    # periodicSpline(ifelse(nx>=0,atan(ny/nx)*180/pi,180+atan(ny/nx)*180/pi)+90,sqrt(nx^2+ny^2),
    # period = 360 )
    pispl <- periodicSpline(ifelse(nx >= 0, atan(ny/nx) * 180/pi, 180 + atan(ny/nx) * 
        180/pi) + 90, sqrt(nx^2 + ny^2), period = 360)
    nnx <- nx[c(2:32, 1)]
    nny <- ny[c(2:32, 1)]
    nnxm <- (nx + nnx)/2
    nnym <- (ny + nny)/2
    nxn <- vector("numeric", 64)
    nyn <- vector("numeric", 64)
    nxn[seq(1, 63, 2)] <- nx
    nxn[seq(2, 64, 2)] <- nnxm
    nyn[seq(1, 63, 2)] <- ny
    nyn[seq(2, 64, 2)] <- nnym
    nxn <- c(nxn, nxn[1])
    nyn <- c(nyn, nyn[1])
    boundary2 <- predict(pispl, ifelse(nxn >= 0, atan(nyn/nxn) * 180/pi, 180 + atan(nyn/nxn) * 
        180/pi) + 90, sqrt(nxn^2 + nyn^2))
    bound32n <- cbind(1:65, boundary2$y * sin(boundary2$x/180 * pi), -1 * boundary2$y * 
        cos(boundary2$x/180 * pi), boundary2$x)
    bound32n <- rbind(c(64, NA, NA, NA), bound32n)
    colnames(bound32n) <- c("sn", "x", "y", "deg")
    status <- 0
    return(list(boundry = bound32n, status = status, crossarea = area, circumference = circ))
}
### bound32R <- function(filename){
bound32R <- function(kei) {
    # browser() kei <- as.vector(t(read.table(filename,header=FALSE)))[1:32]
    xx <- seq(2 * pi/64, 2 * pi + 2 * pi/64, length.out = 33)[1:32]
    pispl <- periodicSpline(xx, kei, period = 2 * pi, ord = 2)
    # boundary <- predict( pispl, seq(2*pi/64, 2*pi+2*pi/64, length.out = 65) )
    boundary <- predict(pispl, seq(0, 2 * pi, length.out = 65))
    # boundary$x <- boundary$x[c(64,1:63)] boundary$y <- boundary$y[c(64,1:63)]
    bound32 <- cbind(1:65, boundary$y * sin(boundary$x), -1 * boundary$y * cos(boundary$x), 
        (boundary$x * 180/pi))
    bound32 <- rbind(c(64, NA, NA, NA), bound32)
    colnames(bound32) <- c("sn", "x", "y", "deg")
    bound1 <- bound32
    bound2 <- as.data.frame(bound1[2:65, ])
    bound2 <- cbind(bound2, radius = sqrt(bound2$x^2 + bound2$y^2))
    rang <- 360/64 * pi/180
    area <- sum(bound2$radius^2/2 * rang)
    circ <- sum(bound2$radius * rang)
    status <- 0
    return(list(boundry = bound32, status = status, crossarea = area, circumference = circ))
} 
