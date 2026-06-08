library(GENIE3)
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
Austria_mgs_ec_mt <- clr(Austria_mgs_ec_mt)
rownames(Austria_mgs_ec_mt) <- rownames(Austria_mgs_ec)

Italy1_mgs_ec_mt <- as.matrix(Italy1_mgs_ec)
Italy1_pseudo <- min(Italy1_mgs_ec[Italy1_mgs_ec>0]) * 0.01
Italy1_mgs_ec_mt <- Italy1_mgs_ec_mt +  Italy1_pseudo
Italy1_mgs_ec_mt <- clr(Italy1_mgs_ec_mt)
rownames(Italy1_mgs_ec_mt) <- rownames(Italy1_mgs_ec)

Italy2_mgs_ec_mt <- as.matrix(Italy2_mgs_ec)
Italy2_pseudo <- min(Italy2_mgs_ec[Italy2_mgs_ec>0]) * 0.01
Italy2_mgs_ec_mt<- Italy2_mgs_ec_mt + Italy2_pseudo
Italy2_mgs_ec_mt <- clr(Italy2_mgs_ec_mt)
rownames(Italy2_mgs_ec_mt) <- rownames(Italy2_mgs_ec)

# construct network
Austria_genie3_net <- GENIE3(Austria_mgs_ec_mt,
                    nCores = 8)

write.table(Austria_genie3_net,
            file = '/nfs/home/users/z.huang/revision_net_construction/Austria_genie3_net.txt',
            row.names = T,
            col.names = T,
            quote = F)

Italy1_genie3_net <- GENIE3(Italy1_mgs_ec_mt,
                             nCores = 8)

write.table(Italy1_genie3_net, 
            file = '/nfs/home/users/z.huang/revision_net_construction/Italy1_genie3_net.txt',    
            row.names = T, 
            col.names = T,
            quote = F)

Italy2_genie3_net <- GENIE3(Italy2_mgs_ec_mt,
                             nCores = 8)


write.table(Austria_genie3_net, 
            file = '/nfs/home/users/z.huang/revision_net_construction/Austria_genie3_net.txt',    
            row.names = T, 
            col.names = T,
            quote = F)
write.table(Italy1_genie3_net, 
            file = '/nfs/home/users/z.huang/revision_net_construction/Italy1_genie3_net.txt',    
            row.names = T, 
            col.names = T,
            quote = F)
write.table(Italy2_genie3_net, 
            file = '/nfs/home/users/z.huang/revision_net_construction/Italy2_genie3_net.txt',    
            row.names = T, 
            col.names = T,
            quote = F)

