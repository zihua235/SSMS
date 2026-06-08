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
Austria_mgs_ec <- as.matrix(Austria_mgs_ec)
t_Austria_mgs_ec <- t(Austria_mgs_ec)
Austria_pseudo <- min(t_Austria_mgs_ec[t_Austria_mgs_ec>0]) * 0.1

Italy1_mgs_ec <- as.matrix(Italy1_mgs_ec)
t_Italy1_mgs_ec <- t(Italy1_mgs_ec)
Italy1_pseudo <- min(t_Italy1_mgs_ec[t_Italy1_mgs_ec>0]) * 0.1

Italy2_mgs_ec <- as.matrix(Italy2_mgs_ec)
t_Italy2_mgs_ec <- t(Italy2_mgs_ec)
Italy2_pseudo <- min(t_Italy2_mgs_ec[t_Italy2_mgs_ec>0]) * 0.1

# construct network
Austria_spearman_net <- netConstruct(t_Austria_mgs_ec,
                                 measure = "spearman", # 
                                 normMethod = "clr", 
                                 zeroMethod = "pseudo",
                                 zeroPar = list(pseudocount = Austria_pseudo),
                                 sparsMethod = "t-test",         
                                 adjust = "adaptBH",
                                 verbose = 3)


Italy1_spearman_net <- netConstruct(t_Italy1_mgs_ec,
                                    measure = "spearman", # 
                                    normMethod = "clr", 
                                    zeroMethod = "pseudo",
                                    zeroPar = list(pseudocount = Italy1_pseudo),
                                    sparsMethod = "t-test",         
                                    adjust = "adaptBH",
                                    verbose = 3)

Italy2_spearman_net <- netConstruct(t_Italy2_mgs_ec,
                                    measure = "spearman", # 
                                    normMethod = "clr", 
                                    zeroMethod = "pseudo",
                                    zeroPar = list(pseudocount = Italy2_pseudo),
                                    sparsMethod = "t-test",         
                                    adjust = "adaptBH",
                                    verbose = 3)

write.table(Italy1_spearman_net$edgelist1, file = "/nfs/home/users/z.huang/revision_net_construction/Italy1_spearman_net.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(Italy2_spearman_net$edgelist1, file = "/nfs/home/users/z.huang/revision_net_construction/Italy2_spearman_net.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(Austria_spearman_net$edgelist1, file = "/nfs/home/users/z.huang/revision_net_construction/Austria_spearman_net.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

