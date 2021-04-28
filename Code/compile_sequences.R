
setwd("/home/cjls4/feature_vectors/")

load("MINUS_STRAND_G4_NEGATIVE_SEQ.RData")
load("MINUS_STRAND_G4_POSITIVE_SEQ.RData")

load("PLUS_STRAND_G4_NEGATIVE_SEQ.RData")
load("PLUS_STRAND_G4_POSITIVE_SEQ.RData")


ALL_STRAND_ALL_G4_SEQ <- rbind(MINUS_STRAND_G4_NEGATIVE_SEQ, 
                               MINUS_STRAND_G4_POSITIVE_SEQ, 
                               PLUS_STRAND_G4_NEGATIVE_SEQ, 
                               PLUS_STRAND_G4_POSITIVE_SEQ)


save(ALL_STRAND_ALL_G4_SEQ, file = "ALL_STRAND_ALL_G4_SEQ.RData")
