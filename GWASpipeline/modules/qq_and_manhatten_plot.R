# Extract parameters from snakemake object.
method <- snakemake@wildcards$method
gwas <- snakemake@wildcards$gwas
goodvarsfile <- snakemake@params$goodvars
model <- snakemake@wildcards$model
assocfile <- snakemake@input$assoc
outfile <- snakemake@output$qqplot

library(data.table)
library(dplyr)
qqp<-function(x,ci=TRUE,add=FALSE,ylab="Observed log10(p-value)",xlab="Expected log10(p-value)",maxLogP,col=1,...){
    if(length(col)>1){
        col<-col[!is.na(x)]
    }
    x<-x[!is.na(x)]
    if(!missing(maxLogP)){
        x[x<10^-maxLogP]<-10^-maxLogP
    }
    N<-length(x)
    ord<-order(x)
    x<-x[ord]
    if(length(col)>1){
        col<-col[ord]
    }
    e<- -log((1:N-0.5)/N,10)
    if(add){
        points(e,-log(x,10),col=col,...)
    } else{
        plot(e,-log(x,10),ylab=ylab,xlab=xlab,col=col,...)
        abline(0,1,col=2,lwd=2)
    }
    if(ci){
        c97.5<-qbeta(0.975,1:N,N-(1:N)+1)
        c02.5<-qbeta(0.025,1:N,N-(1:N)+1) 
        lines(e,-log(c97.5,10))
        lines(e,-log(c02.5,10))
    }
}

if(method %in% c('gemma', 'gemmaloco')){
    pname <- 'p_score'
    chrname <- 'chr'
    idname <- 'rs'
} else {
    stop(paste0("Method '", method, "' not implemented."))
}

assoc <- fread(assocfile, data.table = F, select = c(chrname, pname, idname))
colnames(assoc) <- c('chr', 'p', 'id')

# Only plot good vars.
goodvars <- fread(goodvarsfile, da=F, header=F)[,1]
assoc <- assoc[assoc$id %in% goodvars,]

palette(c("#67a9cf","#2166ac"))
png(outfile, width = 12, height = 4, res=300, units = 'in')
par(family = 'arial', mar=c(3,3,2,1), mgp=c(1.6,0.6,0))
layout(matrix(1:2, nrow = 1), widths = c(0.3,0.7))

# QQ plot
qqp(assoc$p, main=paste(gwas, method, model, 'goodvars', sep=' '),
    las=1, bty='L', pch=19, cex=0.6, col='black')

# Manhatten plot
chrpos <- assoc %>% group_by(chr) %>%
    summarise( n=n() ) %>% as.data.frame()
chrpos$nsum <- cumsum(chrpos$n)
chrpos$half <- round(chrpos$nsum-(chrpos$n/2))

plot(0,0,col="transparent", xaxt="n", las=1, xaxt='n', bty='L',
    xlim=c(0, nrow(assoc)), ylim=c(0, -log10(min(assoc$p))),
    xlab='Chromosome',ylab='Observed log10(p-value)',pch=19, cex=0.4,
    main=paste(gwas, method, model, 'goodvars', sep=' '))
axis(1, at=chrpos$half, labels = chrpos$chr, cex.axis=0.8)
abline(h=-log10(5e-8), lty=2)
points(-log10(assoc$p), col=assoc$chr, pch=19, cex=0.5)
dev.off()

