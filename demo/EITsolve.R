data(bound)
data(exdata)
exdata$datan32[,1:32] <- Kirchhoff(data=exdata$datan32[,1:32],dev=34)$expout
ret <- hosein(data=exdata)
exp <- lstsqrs(data=ret,bound32=bound)
ret <- f9l5qxfa(bound32=bound,expout=exp)
ret <- convsmt()
ret <- filta06("IMG__&_G.PRN",bound32=bound,nbun=255,rfilt=0.1,ratmap=1,col=rich.colors(256)[128:256])
