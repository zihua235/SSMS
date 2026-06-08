source("/nfs/home/users/z.huang/revision_net_construction/fastCCLasso.R")
library(NetCoMi)

# split datasets
mgs_ec <- read.delim("/nfs/home/users/z.huang/revision_net_construction/mgs_level4ec_relab_unstratified.tsv", row.names = 1)
Austria_mgs_ec <- mgs_ec[, grepl("ERR", names(mgs_ec))]
Italy1_id <- read.table("/nfs/home/users/z.huang/revision_net_construction/id_ThomasAM_Italian1.txt", header = FALSE)
Italy2_id <- read.table("/nfs/home/users/z.huang/revision_net_construction/id_ThomasAM_Italian2.txt", header = FALSE)
pattern_1 <- paste(Italy1_id$V1, collapse = "|")
pattern_2 <- paste(Italy2_id$V1, collapse = "|")
Italy1_mgs_ec <- mgs_ec[, grepl(pattern_1, names(mgs_ec))]
Italy2_mgs_ec <- mgs_ec[, grepl(pattern_2, names(mgs_ec))]

# pre-processing
Austria_mgs_ec_mt <- as.matrix(Austria_mgs_ec)
Austria_pseudo <- min(Austria_mgs_ec[Austria_mgs_ec>0]) * 0.01
Austria_mgs_ec_mt<- Austria_mgs_ec_mt + Austria_pseudo
rownames(Austria_mgs_ec_mt) <- rownames(Austria_mgs_ec)

Italy1_mgs_ec_mt <- as.matrix(Italy1_mgs_ec)
Italy1_pseudo <- min(Italy1_mgs_ec[Italy1_mgs_ec>0]) * 0.01
Italy1_mgs_ec_mt <- Italy1_mgs_ec_mt +  Italy1_pseudo
rownames(Italy1_mgs_ec_mt) <- rownames(Italy1_mgs_ec)

Italy2_mgs_ec_mt <- as.matrix(Italy2_mgs_ec)
Italy2_pseudo <- min(Italy2_mgs_ec[Italy2_mgs_ec>0]) * 0.01
Italy2_mgs_ec_mt<- Italy2_mgs_ec_mt + Italy2_pseudo
rownames(Italy2_mgs_ec_mt) <- rownames(Italy2_mgs_ec)

Austria_mgs_ec_mt <- t(Austria_mgs_ec_mt)
Italy1_mgs_ec_mt <- t(Italy1_mgs_ec_mt)
Italy2_mgs_ec_mt <- t(Italy2_mgs_ec_mt)


Austria_fastcclasso_net <- fastCCLasso(Austria_mgs_ec_mt, isCnt = FALSE, k_cv = 3, 
                               lam_min_ratio = 1E-4, k_max = 20, n_boot=100) 
print("austia done")
Italy1_fastcclasso_net <- fastCCLasso(Italy1_mgs_ec_mt, isCnt = FALSE, k_cv = 3, 
                                       lam_min_ratio = 1E-4, k_max = 20, n_boot=100) 
print("italt1 done")
Italy2_fastcclasso_net <- fastCCLasso(Italy2_mgs_ec_mt, isCnt = FALSE, k_cv = 3, 
                                       lam_min_ratio = 1E-4, k_max = 20, n_boot=100) 
print("italy2 fone")
save.image(file = '/nfs/home/users/z.huang/revision_net_construction/cclasso.RData')

