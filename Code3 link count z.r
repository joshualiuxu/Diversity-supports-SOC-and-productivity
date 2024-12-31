
# Load necessary libraries
library(igraph)
library(tidyverse)

# Import data
otu <- read.csv("rarefied_otu_table_euk.csv", row.names=1)
node_types <- read.csv("a.csv", row.names=1)
edgelist <- read.csv("biotic_edge.csv", row.names=1)

# Prepare biotic graph with absolute weights
ig_biotic_abs <- graph_from_edgelist(as.matrix(edgelist[, 1:2]), directed = FALSE)
ig_biotic_abs <- set_edge_attr(ig_biotic_abs, 'weight', index = E(ig_biotic_abs), abs(as.numeric(edgelist[, 3])))

# Assign node types to the graph
ig_biotic_abs <- set_vertex_attr(ig_biotic_abs, "type", index = V(ig_biotic_abs), node_types$type)

# Filter OTU table based on nodes in the graph
otutab <- otu
node <- vertex_attr(ig_biotic_abs)[[1]]
name <- intersect(row.names(otutab), node)
otu <- otutab[name,]

# Display number of nodes and edges
print(length(V(ig_biotic_abs)))  # Number of nodes
print(length(E(ig_biotic_abs)))  # Number of edges

# Construct subgraphs based on OTU presence
sub_graph_ig_biotic_linktotal <- lapply(names(otu), function(i) {
    sample_i <- otu[i]
    select_node <- rownames(sample_i)[which(sample_i > 0)]
    induced_subgraph(ig_biotic_abs, select_node)
})

# Calculate interaction counts between node types
sub_link_count <- lapply(sub_graph_ig_biotic_linktotal, function(ig) {
    count_matrix <- matrix(0, nrow = 4, ncol = 4, dimnames = list(c("Fungi", "Protists", "Metazoa", "Plants"), c("Fungi", "Protists", "Metazoa", "Plants")))
    
    for (e in E(ig)) {
        endpoints <- ends(ig, e)
        from_type <- V(ig)[endpoints[1]]$type
        to_type <- V(ig)[endpoints[2]]$type
        count_matrix[from_type, to_type] <- count_matrix[from_type, to_type] + 1
    }
    
    count_matrix
})

rerow = c("funfun","funpro","funmet","funpla",
           "profun","propro","promet","propla",
           "metfun","metpro","metmet","metpla",
           "plafun","plapro","plamet","plapla")
# Format and export link count data
df_total <- matrix(NA, nrow = 16, ncol = 260, dimnames = list(rerow, colnames(otu)))

for (i in seq_len(260L)) {
    # Convert count matrix to long format and store in 'df_total'
    df <- as.data.frame(sub_link_count[[i]])
    df$Row <- rownames(df)
    long_df <- tidyr::pivot_longer(df, cols = -Row, names_to = "Column", values_to = "Value")
    df_total[, i] <- long_df$Value
}

write.csv(df_total, "df_total.csv")

# 假设df_total是你原始的dataframe

# 选择前四行，并保留原始行名
funfun <- df_total[1, ]  # 第一行
propro <- df_total[6, ]  # 第六行
metmet <- df_total[11, ] # 第十一行
plapla <- df_total[16, ] # 第十六行

# 计算组合行，并为新行指定行名
funpro <- df_total[2, ] + df_total[5, ]   # 第二行加第五行
funmet <- df_total[3, ] + df_total[9, ]   # 第三行加第九行
funpla <- df_total[4, ] + df_total[13, ]  # 第四行加第十三行
promet <- df_total[7, ] + df_total[10, ]  # 第七行加第十行
propla <- df_total[8, ] + df_total[14, ]  # 第八行加第十四行
metpla <- df_total[12, ] + df_total[15, ] # 第十二行加第十五行

# 合并成新的dataframe并设置行名
data_re <- rbind(funfun, propro, metmet, plapla, funpro, funmet, funpla, promet, propla, metpla)

# 为new_df设置行名
rownames(data_re) <- c("funfun", "propro", "metmet", "plapla", "funpro", "funmet", "funpla", "promet", "propla", "metpla")

# 查看结果
print(data_re)

# 如果需要导出为CSV文件，同时保留行名
write.csv(data_re, "data_re.csv", row.names = TRUE)




###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################



###########################################################################################
###########################################################################################
# 第一部分：真实网络分析
# 计算 link count（连接计数）

# Load necessary libraries
library(igraph)      # 图网络分析包
library(tidyverse)   # 包含数据处理、可视化等功能的综合包

# Import data
otu <- read.csv("rarefied_otu_table_euk.csv", row.names = 1)  # 导入 OTU 表格数据
node_types <- read.csv("a.csv", row.names = 1)                # 导入节点类型数据
edgelist <- read.csv("biotic_edge.csv", row.names = 1)        # 导入边列表数据

# Prepare biotic graph with absolute weights
ig_biotic_abs <- graph_from_edgelist(as.matrix(edgelist[, 1:2]), directed = FALSE)  # 使用边列表创建无向图
ig_biotic_abs <- set_edge_attr(ig_biotic_abs, 'weight', index = E(ig_biotic_abs), abs(as.numeric(edgelist[, 3])))  # 设置边的权重为绝对值

# Assign node types to the graph
ig_biotic_abs <- set_vertex_attr(ig_biotic_abs, "type", index = V(ig_biotic_abs), node_types$type)  # 设置节点类型属性

# Filter OTU table based on nodes in the graph
otutab <- otu  # 备份 OTU 表
node <- vertex_attr(ig_biotic_abs)[[1]]  # 获取图中的节点名称
name <- intersect(row.names(otutab), node)  # 找出 OTU 表中存在的节点
otu <- otutab[name, ]  # 过滤 OTU 表格，保留在图中的节点

# Display number of nodes and edges
print(length(V(ig_biotic_abs)))  # 打印节点数量
print(length(E(ig_biotic_abs)))  # 打印边数量

# Construct subgraphs based on OTU presence
sub_graph_ig_biotic_linktotal <- lapply(names(otu), function(i) {
    sample_i <- otu[i]  # 获取当前样本
    select_node <- rownames(sample_i)[which(sample_i > 0)]  # 获取存在 OTU 的节点
    induced_subgraph(ig_biotic_abs, select_node)  # 根据节点构建子图
})

# Calculate interaction counts between node types
sub_link_count <- lapply(sub_graph_ig_biotic_linktotal, function(ig) {
    count_matrix <- matrix(0, nrow = 4, ncol = 4, dimnames = list(c("Fungi", "Protists", "Metazoa", "Plants"), c("Fungi", "Protists", "Metazoa", "Plants")))
    
    # 遍历图中的每一条边，计算节点类型之间的连接计数
    for (e in E(ig)) {
        endpoints <- ends(ig, e)
        from_type <- V(ig)[endpoints[1]]$type
        to_type <- V(ig)[endpoints[2]]$type
        count_matrix[from_type, to_type] <- count_matrix[from_type, to_type] + 1
    }
    
    count_matrix
})

# 行名重新编码，用于后续统计
rerow = c("funfun","funpro","funmet","funpla",
          "profun","propro","promet","propla",
          "metfun","metpro","metmet","metpla",
          "plafun","plapro","plamet","plapla")

# Format and export link count data
df_total <- matrix(NA, nrow = 16, ncol = 260, dimnames = list(rerow, colnames(otu)))

# 遍历每个样本，将 link count 结果格式化并存储在 df_total 中
for (i in seq_len(260L)) {
    df <- as.data.frame(sub_link_count[[i]])
    df$Row <- rownames(df)
    long_df <- tidyr::pivot_longer(df, cols = -Row, names_to = "Column", values_to = "Value")
    df_total[, i] <- long_df$Value
}

write.csv(df_total, "df_total.csv")  # 导出为 CSV 文件

# 选择特定行并进行组合
funfun <- df_total[1, ]
propro <- df_total[6, ]
metmet <- df_total[11, ]
plapla <- df_total[16, ]

# 组合特定行并创建新数据框
funpro <- df_total[2, ] + df_total[5, ]
funmet <- df_total[3, ] + df_total[9, ]
funpla <- df_total[4, ] + df_total[13, ]
promet <- df_total[7, ] + df_total[10, ]
propla <- df_total[8, ] + df_total[14, ]
metpla <- df_total[12, ] + df_total[15, ]

# 合并行结果并设置行名
data_re <- rbind(funfun, propro, metmet, plapla, funpro, funmet, funpla, promet, propla, metpla)
rownames(data_re) <- c("funfun", "propro", "metmet", "plapla", "funpro", "funmet", "funpla", "promet", "propla", "metpla")

# 导出组合结果为 CSV
write.csv(data_re, "data_re.csv", row.names = TRUE)

###########################################################################################
###########################################################################################
# 第二部分：随机网络分析
# 计算随机网络的 link count

# 创建一个空列表以存储随机网络的结果
df_list = list()

# 遍历每个样本的子图
for (i in seq_along(sub_graph_ig_biotic_linktotal)) {
    ig = sub_graph_ig_biotic_linktotal[[i]]
    endpoints <- ends(ig, E(ig))
    
    count_matrix_list = list()

    # 进行 99 次随机化循环
    for (j in 1:99) {
        count_matrix <- matrix(0, nrow = 4, ncol = 4, dimnames = list(c("Fungi", "Protists", "Metazoa", "Plants"), c("Fungi", "Protists", "Metazoa", "Plants")))

        # 随机打乱两列
        shuffled_endpoints <- endpoints
        shuffled_endpoints[, 1] <- sample(endpoints[, 1])  # 随机打乱第一列
        shuffled_endpoints[, 2] <- sample(endpoints[, 2])  # 随机打乱第二列

        # 遍历每个随机边，更新 count_matrix
        for (y in 1:nrow(shuffled_endpoints)) {
            from_node <- shuffled_endpoints[y, 1]
            to_node <- shuffled_endpoints[y, 2]
            from_type <- V(ig)[from_node]$type
            to_type <- V(ig)[to_node]$type
            count_matrix[from_type, to_type] <- count_matrix[from_type, to_type] + 1
        }

        count_matrix_list[[j]] = count_matrix
    }

    # 格式化并存储随机化的 count_matrix 数据
    df_total <- matrix(NA, nrow = 16, ncol = 99, dimnames = list(rerow, c(1:99)))

    for (h in seq_len(99L)) {
        df <- as.data.frame(count_matrix_list[[h]])
        df$Row <- rownames(df)
        long_df <- tidyr::pivot_longer(df, cols = -Row, names_to = "Column", values_to = "Value")
        df_total[, h] <- long_df$Value
    }

    # 选择特定行并进行组合
    funfun <- df_total[1, ]
    propro <- df_total[6, ]
    metmet <- df_total[11, ]
    plapla <- df_total[16, ]

    funpro <- df_total[2, ] + df_total[5, ]
    funmet <- df_total[3, ] + df_total[9, ]
    funpla <- df_total[4, ] + df_total[13, ]
    promet <- df_total[7, ] + df_total[10, ]
    propla <- df_total[8, ] + df_total[14, ]
    metpla <- df_total[12, ] + df_total[15, ]

    # 合并行结果并设置行名
    new_df <- rbind(funfun, propro, metmet, plapla, funpro, funmet, funpla, promet, propla, metpla)
    rownames(new_df) <- c("funfun", "propro", "metmet", "plapla", "funpro", "funmet", "funpla", "promet", "propla", "metpla")

    df_list[[i]] = new_df

    print(i)  # 打印当前样本索引以监控进度
}

###########################################################################################
###########################################################################################
# 第三部分：计算 Z 值和 P 值

# 假设 data_re 是 10 行 * 260 列的真实值矩阵
# df_list 是包含 260 个元素的列表，每个元素是 10 行 * 99 列的随机网络数据框

# 初始化结果存储列表
z_scores_list <- list()
p_values_list <- list()

# 遍历每个样本
for (i in 1:260) {
    # 获取真实值的第 i 列
    real_values <- data_re[, i]

    # 获取对应的随机网络数据，形状为 10 行 * 99 列
    random_values <- df_list[[i]]

    # 计算每行的随机网络均值和标准差
    random_means <- rowMeans(random_values)
    random_sds <- apply(random_values, 1, sd)

    # 计算 Z 值（标准化）
    z_scores <- (real_values - random_means) / random_sds

    # 计算 p 值
    p_values <- mapply(function(observed, random_dist) {
        t.test(random_dist, mu = observed, alternative = "two.sided")$p.value
    }, real_values, as.data.frame(t(random_values)))  # 需要转置

    # 存储结果
    z_scores_list[[i]] <- z_scores
    p_values_list[[i]] <- p_values

    print(i)  # 打印当前样本索引以监控进度
}

# 将列表结果合并为数据框
df_z_scores <- as.data.frame(do.call(cbind, z_scores_list))
df_p_values <- as.data.frame(do.call(cbind, p_values_list))

# 保存为 CSV 文件
write.csv(df_z_scores, "z_scores.csv")
write.csv(df_p_values, "p_values.csv")


###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
