# plots of relative abundance vs. dispersion in different types of datasets

library(ALDEx2)
# metatranscriptome
t4 <-  read.table("data/twntyfr.txt", header=T, sep="\t", quote="", comment.char="", row.names=1)

conds.t <-c("H","H","H","H","B","H","B","B","H","B","H","B","B","B","B","B","B","B","H","B","H","H")

t <- aldex.clr(t4, conds.t)
t.e <- aldex.effect(t)
t.lvha <- aldex.clr(t4, conds.t, denom="lvha")
t.iqlr <- aldex.clr(t4, conds.t, denom="iqlr")

# tiyani
tiyani <- read.table("data/pup_mage.txt", header=T, row.names=1)
conds.ty <- c(rep("P", 161),rep("M",86))
x.ty <- aldex.clr(tiyani, conds.ty)
x.ty.lvha <- aldex.clr(tiyani, conds.ty, denom="lvha")
x.ty.iqlr <- aldex.clr(tiyani, conds.ty, denom="iqlr")
x.ty.e <- aldex.effect(x.ty)

# barton yeast
yeast <- read.table("data/barton_yeast.txt", header=T, row.names=1)
conds.yeast <- c(rep("S",41),rep("W",43))
x.yeast <- aldex.clr(yeast, conds.yeast)
x.yeast.iqlr <- aldex.clr(yeast, conds.yeast, denom="iqlr")
x.yeast.lvha <- aldex.clr(yeast, conds.yeast, denom="lvha")
x.e.yeast <- aldex.effect(x.yeast)

data(selex)
conds.s <- c(rep("A", 7),rep("B",7))
s <- aldex.clr(selex, conds.s)
s.iqlr <- aldex.clr(selex, conds.s, denom="iqlr")
s.lvha <- aldex.clr(selex, conds.s, denom="lvha")
s.e <- aldex.effect(s)

setEPS()
postscript("vera-book/Disp-v-Abund.eps", height=7, width=6, horizontal=F)
par(mfrow=c(2,2))

plot(x.e.yeast$diff.win, x.e.yeast$rab.all, col="grey80", pch=19,
   xlim=c(0,8), ylim=c(-10,15), xlab="Dispersion", ylab="log2 Relative Abundance",
   main="transcriptome")
points(x.e.yeast$diff.win[x.yeast.iqlr@denom], x.e.yeast$rab.all[x.yeast.iqlr@denom],
    col="blue", pch=19)
points(x.e.yeast$diff.win[x.yeast.lvha@denom], x.e.yeast$rab.all[x.yeast.lvha@denom],
    col="red", pch=19)
legend("topright", inset=0.2, legend=c("all", "LVHA", "IQLR"), col=c("grey80", "red", "blue"),
    pch=19)

plot(t.e$diff.win, t.e$rab.all, col="grey80", pch=19,
   xlim=c(0,8), ylim=c(-10,15), xlab="Dispersion", ylab="log2 Relative Abundance",
   main="meta-transcriptome")
points(t.e$diff.win[t.iqlr@denom], t.e$rab.all[t.iqlr@denom],
    col="blue", pch=19)
points(t.e$diff.win[t.lvha@denom], t.e$rab.all[t.lvha@denom],
    col="red", pch=19)

plot(x.ty.e$diff.win, x.ty.e$rab.all, col="grey80", pch=19,
   xlim=c(0,8), ylim=c(-10,15), xlab="Dispersion", ylab="log2 Relative Abundance",
   main="microbiome")
points(x.ty.e$diff.win[x.ty.iqlr@denom], x.ty.e$rab.all[x.ty.iqlr@denom],
    col="blue", pch=19)
points(x.ty.e$diff.win[x.ty.lvha@denom], x.ty.e$rab.all[x.ty.lvha@denom],
    col="red", pch=19)

plot(s.e$diff.win, s.e$rab.all, col="grey80", pch=19,
   xlim=c(0,8), ylim=c(-10,15), xlab="Dispersion", ylab="log2 Relative Abundance",
   main="selex")
points(s.e$diff.win[s.iqlr@denom], s.e$rab.all[s.iqlr@denom],
    col="blue", pch=19)
points(s.e$diff.win[s.lvha@denom], s.e$rab.all[s.lvha@denom],
    col="red", pch=19)
dev.off()

