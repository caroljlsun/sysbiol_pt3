##### methylKit for converting the methylation data, and quality control #####

##### Reading methylation files #####

# Here, we are going to read in the methylation data into an object called
# 'methylRawList', which stores methyl info per sample for each base
# requires a minimum coverage of 10 reads

# methylkit is also used for comparing methylation of control vs treatment
# which we aren't really doing here, so some of the methods will not be
# relevant to us here

# all we're going to do here now is go through their provided examples!


# load in the library, ofc
library(methylKit)

# find the example files using 'system.files'. "extdata" refers to the 
# additional external data in packages, usually they are for examples or
# vignettes

# Here we will read in high throughput bisulfite sequencing call files
file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

what_format= system.file("extdata", "test1.myCpG.txt", package = "methylKit")

# so we have made a list of sample data, all of which are txt files
# with test and control conditions, and 2 samples per condition

# now we will read the files into 'myobjDB'
# and save into a database, in a folder called 'methylDB'

myobjDB=methRead(file.list,
                 sample.id=list("test1","test2","ctrl1","ctrl2"),
                 assembly="hg18",
                 treatment=c(1,1,0,0),
                 context="CpG",
                 dbtype = "tabix",
                 dbdir = "methylDB")

# arguments of methRead()
# (location, sample.id, assembly, dbtype = NA, pipeline = "amp",
# header = TRUE, skip = 0, sep = "\t", context = "CpG",
# resolution = "base", treatment, dbdir = getwd(), mincov = 10)
# 
# assembly is a string which tells you which reference genome it came from.
# context, is it CpG, CpH, CHH etc?
# resolution can be base or region
# dbdir is the directory where the flat file db will be saved. Default is wd
# mincov is the minimum read coverage. Default is 10. Anything lower will be
# ignored.

# Here, we'll get where the files are saved
print(myobjDB[[1]]@dbpath)
# the @ can get you other info, e.g. use @sample.id or @num.records

# we can also read methylation calls from sorted Bismark alignments, e.g.
# BAM or SAM files
# I'll ignore this for now, but the info is there on the website
# https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html

##### Basic stats #####

# So we have read in the methylation data. Now we will check some basic stats
# such as coverage and percent methylation

getMethylationStats(myobjDB[[2]],plot=FALSE,both.strands=FALSE)

# what does 'plot=' do? Well:

getMethylationStats(myobjDB[[2]],plot=TRUE,both.strands=FALSE)
# it plots a cute little histogram for the %methylation per base

getMethylationStats(myobjDB[[2]],plot=TRUE,both.strands=TRUE)
# this gives you the histogram,, for both strands

# to get read coverage per base, use:

getCoverageStats(myobjDB[[2]],plot=TRUE,both.strands=FALSE)

##### Filtering samples #####
# based on read coverage

filtered.myobjDB = filterByCoverage(myobjDB, lo.count = 10, 
                                    lo.perc = NULL,
                                    hi.count = NULL,
                                    hi.perc = 99.9)
# parsing arguments
# lo.count = bases with lower coverage than this will be discarded. 
# good for better statistical tests
# hi.count = bases w/ higher coverage than this will be discarded. Good for
# eliminating PCR effects
# lo.perc = [0-100] percentile of read counts. Bases w/ lower coverage than this is 
# discarded
# and vice versa for hi.perc

##### Comparative analysis #####

# May be less relevant to us but oh well let's learn

# We will need to 'get' the bases covered in all samples.
# the following function will merge all samples to one object for base pair locations
# covered in all samples...
# 'destrand = TRUE' will merge reads on both strands, which provides better coverage
# but only suitable for CpG methylation, and only works if you have base pair res.
# 'unite' will return a methylBase object for comparative analysis, and contains
# methyl info for bases covered in all samples
# so I guess it ignores bases that are covered in only some of the samples

meth=unite(myobjDB, destrand = F)

# let's have a look at the object

head(meth)
# for some reason I only get one row, and they get 6 in the vignette

# to relax the conditions (aka base exists in 3 of 4 samples or whatever)

meth.min=unite(myobjDB, min.per.group = 1L)
head(meth.min)
# again I only get one row for some strange reason

##### Sample correlation #####

# I think this would be good to use to test the samples for correlation, for 
#  m, f, b, t, and mincov_10

getCorrelation(meth, plot = T)
# we get a correlation matrix, as well as pearson's correlation coeff.
# and maybe a scatter plot (?)

##### Clustering samples #####

# cluster according to similarity of methylation profiles

clusterSamples(meth, dist = "correlation", method = "ward", plot = T)

# This function gives you a nice dendrogram of your clusters
# parsing parameters: 'dist' = distance measure. incl. euclidean, maximum, manhattan
# 'method'= agglomeration method , incl. single, complete, average, mcquitty
# agglomeration is a bottom up approach to hierarchical clustering, where each obs
# starts in its own cluster, and pairs are merged as they go up the hierarchy

clusterSamples(meth, dist = "correlation", method = "ward", plot = F)
# this can be used to return a dendrogram object which can be further manipulated by users

hc = clusterSamples(meth, dist = "correlation", method = "ward", plot = T)
# ?

# PCA analysis, scree plot then PC1 vs PC2

PCASamples(meth, screeplot = T)
PCASamples(meth)

##### Batch effects! #####

# make some batch data frame
# this is a bogus data frame
# we don't have batch information
# for the example data
sampleAnnotation=data.frame(batch_id=c("a","a","b","b"),
                            age=c(19,34,23,40))

as=assocComp(mBase=meth,sampleAnnotation)
as

# mbase is the methylBase object, must not have NA values
# sampleAnnotation is a data frame where columns are different annotations 
# and rows are the samples, in the same order as in the methylBase object.

# construct a new object by removing the first principal component
# from percent methylation value matrix
newObj=removeComp(meth,comp=1)

# Error in fread(methMat, header = FALSE, nrows = chunk) :
# input= must be a single character string containing a file name,
# a system command containing at least one space, a URL starting 
# 'http[s]://', 'ftp[s]://' or 'file://', or, the input data itself 
# containing at least one \n or \r

#afadjskfsad

##### Tiling window analysis #####

# for doing region/window resolution instead of base-pair resolution analysis
# e.g. the following code will tile the genome into windows of 1000bps
# and the methylation info will be summarised per 1000bps

myobj_lowCov = methRead(file.list,
                        sample.id=list("test1","test2","ctrl1","ctrl2"),
                        assembly="hg18",
                        treatment=c(1,1,0,0),
                        context="CpG",
                        mincov = 3
)


# notice mincov = 3, we can filter later based on the number of C's per 1000bp

tiles = tileMethylCounts(myobj_lowCov,win.size=1000,step.size=1000,cov.bases = 10)

head(tiles[[1]],3)

##### Finding differentially methylated bases or regions #####

# use the function 'calculateDiffMeth()'
# can use fisher's exact or logistic regression to calculate p-values, dep. on sample
# size
# if you have replicates, will use fishers
# the logistic regression essentially tests whether b^0+b^1*treament does a better job
# of predicting the log-odds ratio of methylation than just b^0

myDiff=calculateDiffMeth(meth)

# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
#
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
#
#
# get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

# the q-value gives the expected pFDR (positive false discovery rate) obtained by 
# rejecting the null hypothesis for any result with an equal or smaller q-value.

# to visualise the distribution of the differentially methylated regions

diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

##### Correcting for overdispersion #####

# overdispersion is when there is more variability in the data than assumed by the 
# distribution
# the logistic regression models the response variable meth as a binomial dist.
# therefore the variance = n*pi(1-pi) and the mean = n*pi
# where n = coverage for the base/region, pi = underlying methylation proportion

# to deal with overdispersion, add the scaling parameter 'phi' to the var formula:
# phi*n*pi(1-pi)
# to calculate overdispersion, use 'overdispersion = "MN"'

sim.methylBase1<-dataSim(replicates=6,sites=1000,
                         treatment=c(rep(1,3),rep(0,3)),
                         sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
)
# just some simulated data
# 'dataSim' simulates all kinds of DNA methylation data

my.diffMeth<-calculateDiffMeth(sim.methylBase1[1:3], 
                               overdispersion="MN",test="Chisq",mc.cores=1)

##### Accounting for covariation #####
# I assume this means, effects of sex, of time of collection etc.
# the package can test if a logistic regression model with covariates and treatment
# ( a full model), does better than one without treatment
 
# you'll have to supply the covariation in the form of a dataframe
# e.g. here the dataframe is of age
covariates=data.frame(age=c(30,80,34,30,80,40))

sim.methylBase<-dataSim(replicates=6,sites=1000,
                        treatment=c(rep(1,3),rep(0,3)),
                        covariates=covariates,
                        sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
)

my.diffMeth3<-calculateDiffMeth(sim.methylBase,
                                covariates=covariates,
                                overdispersion="MN",test="Chisq",mc.cores=1)

##### Annotating differentially methylated bases or regions #####
# by using the 'genomation' package and reading in gene annotations from a BED file
# you will have to force the methylKit objects into a GRanges object to do this

library(genomation)
# read the gene BED file
gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", 
                                            package = "methylKit"))
# annotate differentially methylated CpGs with 
# promoter/exon/intron using annotation data

annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
# The function annotates GRangesList or GRanges object as overlapping with
# promoter,exon,intron or intergenic regions. In this case, methylated bases


# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
                                     package = "methylKit"),
                         feature.flank.name=c("CpGi","shores"))

# convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

##### Regional analysis #####

# for when you have regions of interest to summarise methylation data for
# such as promoters, CpG islands
# 
promoters=regionCounts(myobjDB,gene.obj$promoters)

head(promoters)
# lol I get nothing

##### Convenience function #####
# get the distance from a TSS and or neared gene name with

diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

head(getAssociationWithTSS(diffAnn))
# target.row is the row number in myDiff25p

#It is also desirable to get percentage/number of differentially methylated regions 
#that overlap with intron/exon/promoters

getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

#We can also plot the percentage of differentially methylated bases overlapping
# with exon/intron/promoters

plotTargetAnnotation(diffAnn,precedence=TRUE,
                     main="differential methylation annotation")

# percentage of differentially methylated bases are on CpG islands,
# CpG island shores and other regions.

plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
                     main="differential methylation annotation")

# to get percentage of intron/exon/promoters that overlap with differentially
# methylated bases.

getFeatsWithTargetsStats(diffAnn,percentage=TRUE)

##### methylKit convenience functions #####

class(meth)
as(meth, "GRanges")

class(myDiff)
as(myDiff, "GRanges")

class(myobjDB[[1]])
as(myobjDB[[1]], "methylRaw")
