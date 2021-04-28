library(goProfiles)
library(org.Mm.eg.db)
library(ChIPpeakAnno)

load("/home/cjls4/feature_vectors/top_100_names.RData")
load("/home/cjls4/feature_vectors/top_100_names_rf.RData")
load("/home/cjls4/feature_vectors/top_100_names_svm.RData")

setwd("/home/cjls4/feature_vectors/")

top_100_names_cnn_entrez <- convert2EntrezID(IDs = top_100_names, orgAnn = "org.Mm.eg.db", ID_type = "gene_symbol")
top_100_names_rf_entrez <- convert2EntrezID(IDs = top_100_names_rf, orgAnn = "org.Mm.eg.db", ID_type = "gene_symbol")
top_100_names_svm_entrez <- convert2EntrezID(IDs = top_100_names_svm, orgAnn = "org.Mm.eg.db", ID_type = "gene_symbol")


# goProfiles::basicProfile(top_100_names_cnn_entrez, orgPackage = "org.Mm.eg.db")
# 
# a1 <- compareGeneLists(top_100_names_cnn_entrez, top_100_names_rf_entrez, orgPackage = "org.Mm.eg.db")
# 
# print(compSummary(a1))

cnn.MF <- basicProfile(top_100_names_cnn_entrez, onto = "MF", level = 2, orgPackage = "org.Mm.eg.db")
rf.MF <- basicProfile(top_100_names_rf_entrez, onto = "MF", level = 2, orgPackage = "org.Mm.eg.db")
svm.MF <- basicProfile(top_100_names_svm_entrez, onto = "MF", level = 2, orgPackage = "org.Mm.eg.db")

#setwd

setwd("/home/cjls4/Plots/")

#cnns and rfs
cnn.rf.MF <- mergeProfilesLists(cnn.MF, rf.MF, profNames = c("CNN", "RF"))

plotProfiles(cnn.rf.MF, aTitle = "CNN vs RF GO profile", percentage = T, colores = c("white", "black"),
             labelWidth = 70, legend = T)
#cnns and svms

cnn.svm.MF <- mergeProfilesLists(cnn.MF, svm.MF, profNames = c("CNN", "SVM"))

plotProfiles(cnn.svm.MF, aTitle = "CNN vs SVM GO profile", percentage = T, colores = c("white", "black"),
             labelWidth = 70, legend = T)
#rfs and svms

rf.svm.MF <- mergeProfilesLists(rf.MF, svm.MF, profNames = c("RF", "SVM"))

plotProfiles(rf.svm.MF, aTitle = "RF vs SVM GO profile", percentage = T, colores = c("white", "black"),
             labelWidth = 70, legend = T)



plotProfiles(cnn.MF)

cnn.BP <- basicProfile(top_100_names_cnn_entrez, onto = "BP", level = 2, orgPackage = "org.Mm.eg.db")
plotProfiles(cnn.BP)

cnn.CC <- basicProfile(top_100_names_cnn_entrez, onto = "CC", level = 2, orgPackage = "org.Mm.eg.db")
plotProfiles(cnn.CC)

#
data(kidneyGeneLists)
kidneyGeneLists
genListsClusters <- iterEquivClust(kidneyGeneLists, ontoLevels = 2:3,
                                   jobName = "Kidney Gene Lists_Equivalence Clustering (complete)",
                                   ylab = "Equivalence threshold distance",
                                   orgPackage="org.Hs.eg.db", method = "complete")
genListsClusters[["BP"]][["Level 3"]]
class(genListsClusters[["BP"]][["Level 3"]])

plot(genListsClusters$MF$Level2)
#equivalentGOProfiles(goObject, ...)
