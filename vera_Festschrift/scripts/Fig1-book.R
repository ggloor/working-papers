################################################################################
##### Fig1.R
##### Author: Jia Rong Wu, modified by Greg Gloor
##### jwu424 (at) gmail.com
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures
##### for the paper IQLR. Code for Figure 1.
#####
##### USAGE: Rscript --vanilla Fig1.R
#####
##### LICENSE
##### Copyright (c) 2016,2019 Jia Rong Wu, GBG
#####
##### Permission is hereby granted, free of charge, to any person obtaining a
##### copy of this software and associated documentation files (the "Software"),
##### to deal in the Software without restriction, including without limitation
##### the rights to use, copy, modify, merge, publish, distribute, sublicense,
##### and/or sell copies of the Software, and to permit persons to whom the
##### Software is furnished to do so, subject to the following conditions:
#####
##### The above copyright notice and this permission notice shall be included in
##### all copies or substantial portions of the Software.
#####
##### THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
##### IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
##### FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
##### THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
##### LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
##### FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
##### DEALINGS IN THE SOFTWARE.
################################################################################

##### Load Required Libraries
source("scripts/Variables.R")

################################### FIGURE 1a ##################################
##### This figure is generated off dataset reads
##### It is the original unmodified base dataset with 40 positive controls
##### This is original ALDEx2 denom="all"
################################### FIGURE 1 ##################################
conds <- c(rep("A", 10), rep("B", 10))

##### Read table and generate conditions
x.clr <- aldex.clr(reads, conds, mc.samples, verbose=FALSE, denom="all")
x.e <- aldex.effect(x.clr)
x.t  <- aldex.ttest(x.clr)
x.all <- data.frame(x.e, x.t)

##### Use assymmetric datasets

x.clr.0 <- aldex.clr(reads.0, conds, mc.samples, verbose=FALSE, denom="all")
x.e.0 <- aldex.effect(x.clr.0)
x.t.0  <- aldex.ttest(x.clr.0)
x.all.0 <- data.frame(x.e.0, x.t.0)


all.col=rgb("black")
all.pch=19
all.cex=0.4
called.col="red"
called.pch=20
called.cex=0.6
true.col="blue"
true.pch=19
true.cex=0.6
asym.col="cyan"

thres.line.col="darkgrey"
thres.lwd=1.5
test="wilcox"
cutoff=0.1
rare.col="black"
rare=0
rare.pch=20
rare.cex=0.2



f1 <- paste(figs.dir, "Fig_1-book.eps",sep="")
setEPS()
postscript(f1, height=7, width=3.5)
par(fig=c(0,1,0,1), new=TRUE, cex=0.8)
par(fig=c(0,1,0.5,1), new=TRUE)
	called <- x.all$wi.eBH <= cutoff

	plot(x.all$diff.win, x.all$diff.btw, xlab=xlab, ylab=ylab,
	  col=all.col, pch=all.pch, cex=all.cex, main="Symmetric",
	  ylim=c(ymin,ymax))

	points(x.all$diff.win[x.all$rab.all < rare], x.all$diff.btw[x.all$rab.all < rare],
	  col=rare.col,
	  pch=rare.pch, cex=rare.cex)

	points(x.all$diff.win[called], x.all$diff.btw[called],
	  col=called.col, pch=called.pch, cex=called.cex)

	points(x.all[true.set,"diff.win"], x.all[true.set,"diff.btw"],
	  col=true.col, pch=true.pch, cex=true.cex)

	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)
	legend("topright", inset= 0.1, legend=c("TP", "FP", "asymmetric", "NS"),
	  pch=19, col=c("blue", "red","cyan", "black"))


par(fig=c(0.1,0.6, 0.675,0.975), new=TRUE)
	hist(x.all$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "" ),
	  xlab="", ylab="", cex.main=0.6)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)


par(fig=c(0,1,0,.50), new=TRUE)
	called <- x.all.0$wi.eBH <= cutoff
	lv.set <- rownames(x.all.0[called,] & x.all.0$diff.win < 1.6)
	FP.set <- setdiff(lv.set, true.set)
    asym <- x.all.0$wi.eBH <= cutoff & x.all.0$diff.win > 1.6

	plot(x.all.0$diff.win, x.all.0$diff.btw, xlab=xlab, ylab=ylab,
	  col=all.col, pch=all.pch, cex=all.cex, main="2% asymmetric", ylim=c(-4,ymax))

	points(x.all.0$diff.win[x.all.0$rab.all < rare], x.all.0$diff.btw[x.all.0$rab.all < rare],
	  col=rare.col, pch=rare.pch, cex=rare.cex)

	points(x.all.0[FP.set,"diff.win"], x.all.0[FP.set,"diff.btw"], col=called.col,
	  pch=called.pch, cex=called.cex)

	points(x.all.0[asym,"diff.win"], x.all.0[asym,"diff.btw"], col=asym.col,
	  pch=called.pch, cex=called.cex)

	points(x.all.0[true.set,"diff.win"], x.all.0[true.set,"diff.btw"], col=true.col,
	  pch=true.pch, cex=true.cex)

	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.1,0.6, 0.175,0.475), new=TRUE)
	hist(x.all.0$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "" ), xlab="", ylab="", cex.main=0.6)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)

dev.off()

