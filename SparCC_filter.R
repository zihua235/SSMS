# 合并相关矩阵的P值矩阵
library(igraph)


# 观测值d额相关矩阵
cor_sparcc <- read.csv("/Users/zihuahuang/tools/SparCC3/example/basis_corr/cor_sparcc.out", sep = "\t", row.names = 1, check.names = FALSE)
table(abs(cor_sparcc) > 0.05)
#str(cor_sparcc)
#P值矩阵
pvals <- read.csv('/Users/zihuahuang/tools/SparCC3/example/pvals/cov_mat_SparCC.out', sep = "\t",row.names = 1, check.names = FALSE)
#str(pvals)
#保留 p <= 0.05 + cor >= 0.05的值
cor_sparcc[abs(cor_sparcc) <= 0.05 | pvals >= 0.05] <- 0
diag(cor_sparcc) <- 0
table(cor_sparcc != 0)
cor_sparcc[abs(cor_sparcc) > 0] <- 1

write.table(cor_sparcc, '/Users/zihuahuang/Documents/SSMS_network/SparCC/network.adj.txt', col.names = NA, sep = '\t', quote = FALSE)

# 统计图的信息
### 对于无向图，节点的度是邻接矩阵中每行或每列的非零元素的和
node_degrees <- rowSums(cor_sparcc)
table(node_degrees == 0)
### 对于无向图，可以直接从邻接矩阵中获取每个节点的邻居
node_neighbors <- lapply(seq_len(nrow(cor_sparcc)), function(i) which(cor_sparcc[i, ] != 0))
str(node_neighbors)





