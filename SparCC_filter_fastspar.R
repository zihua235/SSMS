rm(list=ls())

library(ggplot2)


# 观测值d额相关矩阵
cor_sparcc <- read.csv("/Users/zihuahuang/Documents/SSMS/SparCC/fastspar/DeFilippisF_2019_median_correlation.tsv", sep = "\t", row.names = 1, check.names = FALSE)
dim(cor_sparcc)
table(abs(cor_sparcc) > 0.03)
#str(cor_sparcc)
#P值矩阵
pvals <- read.csv('/Users/zihuahuang/Documents/SSMS/SparCC/fastspar/DeFilippisF_2019_pvalues.tsv', sep = "\t",row.names = 1, check.names = FALSE)

#保留 p <= 0.05 + cor >= 0.05的值
cor_sparcc[abs(cor_sparcc) <= 0.06 | pvals >= 0.05] <- 0    #replace
diag(cor_sparcc) <- 0
table(cor_sparcc != 0)
cor_sparcc[abs(cor_sparcc) > 0] <- 1
table(abs(cor_sparcc) > 0)
write.table(cor_sparcc, file = '/Users/zihuahuang/Documents/SSMS/SparCC/DeFilippisF_2019_0.06_0.05_network.adj.txt', col.names = NA, sep = '\t', quote = FALSE)

# 统计图的信息
### 对于无向图，节点的度是邻接矩阵中每行或每列的非零元素的和
node_degrees <- rowSums(cor_sparcc)
table(node_degrees == 0)
### 对于无向图，可以直接从邻接矩阵中获取每个节点的邻居
node_neighbors <- lapply(seq_len(nrow(cor_sparcc)), function(i) which(cor_sparcc[i, ] != 0))
degree_distribution <- sapply(node_neighbors, length)
str(as.data.frame(degree_distribution))
max(as.data.frame(degree_distribution))
min(as.data.frame(degree_distribution))
# plot of degree
ggplot(as.data.frame(degree_distribution), aes(x = degree_distribution)) +
geom_bar(fill = "lightblue", color = "black") +
  xlab("number of Neighbor") +
  ylab("Frequency")

