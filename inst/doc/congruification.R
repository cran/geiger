### R code from vignette source 'congruification.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: congruification.Rnw:26-29
###################################################
par(mfrow=c(1,1))
par(mar=c(3,3,1,0)+.1)
options(width=80)


###################################################
### code chunk number 2: congruification.Rnw:44-47
###################################################
require(geiger)
amph=get(data(amphibia)) ## load data and reassign to an object 'amph'
rl=amph$roelants ## Roelants et al. (2007) chronogram


###################################################
### code chunk number 3: congruification.Rnw:53-54
###################################################
gbresolve(c("Liua", "Hynobius", "Ranodon", "Andrias", "Siren"), rank="order", within="Caudata")


###################################################
### code chunk number 4: congruification.Rnw:62-67
###################################################
trl=gbresolve(rl, within="Amphibia", rank=c("genus", "family"))$tax
ex=subset(rl, tax=trl, rank="genus")

## show original and new tip labels of the exemplar phylogeny
print(cbind(original=rl$tip.label[1:6], exemplar=ex$tip.label[1:6]))


###################################################
### code chunk number 5: congruification.Rnw:73-75
###################################################
plot.phylo(ladderize(ex, right=FALSE), cex=0.25, label.offset=3, x.lim=c(-10, 450))
axisPhylo(cex.axis=0.75)


###################################################
### code chunk number 6: congruification.Rnw:82-84
###################################################
plot.phylo(ladderize(ex, right=FALSE), type="phylogram", cex=0.25, label.offset=3, no.margin=FALSE, x.lim=c(-10, 450))
axisPhylo(cex.axis=0.75)


###################################################
### code chunk number 7: congruification.Rnw:92-97
###################################################
pw=ladderize(amph$pyronwiens, right=FALSE) ## Pyron and Wiens (2011) phylogram

## resolve taxonomy of tips (restricted to within amphibians)
tmp=gbresolve(pw, within="Amphibia", rank=c("genus", "class"))
tax=tmp$tax


###################################################
### code chunk number 8: congruification.Rnw:101-102
###################################################
head(tax[,c("family", "order")])


###################################################
### code chunk number 9: congruification.Rnw:108-109
###################################################
plot.phylo(pw, type="fan", show.tip=FALSE, edge.width=0.2, no.margin=TRUE)


###################################################
### code chunk number 10: congruification.Rnw:114-115
###################################################
plot.phylo(pw, type="fan", show.tip=FALSE, edge.width=0.01, no.margin=TRUE)


###################################################
### code chunk number 11: congruification.Rnw:123-124
###################################################
res=congruify.phylo(reference=ex, target=pw, taxonomy=tax, scale=NA)


###################################################
### code chunk number 12: congruification.Rnw:130-135
###################################################
cal=res$calibrations

## options ('opts') can be user-specific 
write.treePL(pw, cal, base="amphibia", nsites=12712, opts=list(smooth=0.1, nthreads=2, 
	opt=1, optad=1, thorough=TRUE)) 


###################################################
### code chunk number 13: congruification.Rnw:141-142
###################################################
head(cal[,-which(names(cal)=="MRCA")]) ## excluding the hash keys in the first column


###################################################
### code chunk number 14: congruification.Rnw:146-153
###################################################
phy=amph$congruified

## exclude genus, tribe, and subfamily from the labels (the first three columns of 'tax')
phy=nodelabel.phylo(phy, tax[,-match(c("genus", "tribe", "subfamily"), colnames(tax))], strict=TRUE) 
plot.phylo(phy, show.tip=FALSE, edge.width=0.1, no.margin=TRUE, x.lim=c(-10,450))
axisPhylo(cex.axis=0.75)
nodelabels(phy$node.label, frame="none", col="lightskyblue", cex=0.85)


###################################################
### code chunk number 15: congruification.Rnw:160-166
###################################################
oo=phy
oo$tip.label=rep(paste(rep("",5),collapse=" "),Ntip(phy))
plot.phylo(oo, type="phylogram", show.tip=TRUE, label.offset=0, edge.width=0.1, x.lim=c(-10, 450))
axisPhylo()
xx=substring(phy$node.label, nchar(phy$node.label)-3, nchar(phy$node.label))=="idae"
nodelabels(phy$node.label, frame="none", col="lightskyblue", cex=ifelse(xx==TRUE, 0.25, 0.85))


###################################################
### code chunk number 16: congruification.Rnw:175-195
###################################################
mrcaID=function(phy, cal){
	cal=as.matrix(cal)
	res=sapply(1:nrow(cal), function(idx){ ## loop over rows
			tips=cal[idx, c("taxonA", "taxonB")] ## fetch spanning taxa
			return(geiger:::.mrca(tips, phy)) ## MRCA of spanning taxa (node ID)
	})
	N=Ntip(phy)
	n=Nnode(phy)
	nn=integer(N+n) ## create empty vector of same length as branches in tree
	nn[res]=1 ## identify nodes that appear within calibrations
	nn=nn[-c(1:N)] ## exclude tip branches
	return(nn) ## return vector (ordered from first to last internal node in the tree)
}

vec=mrcaID(phy, cal) ## get vector for calibrated nodes
sum(vec)==nrow(cal) ## check on whether the function is working appropriately
plot.phylo(phy, type="fan", show.tip=FALSE, edge.width=0.1)

## plot box at node only if calibrated
nodelabels(text=NULL, cex=ifelse(vec==0, NA, 2), frame="n", bg="white", pch=21)


###################################################
### code chunk number 17: congruification.Rnw:201-204
###################################################
vec=mrcaID(phy, cal)
plot.phylo(phy, type="fan", show.tip=FALSE, edge.width=0.05, no.margin=TRUE) 
nodelabels(text=NULL, cex=ifelse(vec==0, NA, 2), frame="n", bg="white", pch=21) 


