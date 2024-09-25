library(umap)
library(uwot)
library(ggplot2)
library(DDoutlier)
library(diffudist)
library(fitdistrplus)
library(igraph)
library(MixGHD)
library(BiocParallel)
library(loe)
library(mclust)
library(cluster)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggrepel)
library(aricode)
library(Rtsne)
library(BiocParallel)



set.seed(2024) 

get_mapping <- function(dataset) {
  mapping_Oihane <- c(
    "Astrocytes" = "Astro",
    "Endothelial" = "Endo",
    "Ependymal" = "Epend",
    "Hybrid" = "Hyb",
    "Microglia" = "Micro",
    "Neurons" = "Neur",
    "Oligodendrocytes" = "Oligo"
  )

  mapping_Quake <- c(
    "B cell" = "BC",
    "ciliated columnar cell of tracheobronchial tree" = "CCCT",
    "classical monocyte" = "ClMono",
    "epithelial cell of lung" = "EpLung",
    "leukocyte" = "Leuk",
    "lung endothelial cell" = "LungEnd",
    "monocyte" = "Mono",
    "myeloid cell" = "Myel",
    "natural killer cell" = "NK",
    "stromal cell" = "Strom",
    "T cell" = "TC"
  )

  mapping_mir <- c(
    "AUTONOMIC_GANGLIA" = "AG",
    "BONE" = "Bo",
    "BREAST" = "Br",
    "CENTRAL_NERVOUS_SYSTEM" = "CNS",
    "ENDOMETRIUM" = "En",
    "FIBROBLAST" = "Fb",
    "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE" = "HL",
    "KIDNEY" = "Ki",
    "LARGE_INTESTINE" = "LIN",
    "LIVER" = "LIV",
    "LUNG" = "LUN",
    "OESOPHAGUS" = "Oe",
    "OVARY" = "Ov",
    "PANCREAS" = "Pa",
    "SKIN" = "Sk",
    "SOFT_TISSUE" = "STI",
    "STOMACH" = "STO",
    "THYROID" = "Th",
    "UPPER_AERODIGESTIVE_TRACT" = "UAT",
    "URINARY_TRACT" = "UT"
  )
    mapping_gene <- c(
    "AUTONOMIC_GANGLIA" = "AG",
    "BONE" = "Bo",
    "BREAST" = "Br",
    "CENTRAL_NERVOUS_SYSTEM" = "CNS",
    "ENDOMETRIUM" = "En",
    "FIBROBLAST" = "Fb",
    "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE" = "HL",
    "KIDNEY" = "Ki",
    "LARGE_INTESTINE" = "LIN",
    "LIVER" = "LIV",
    "LUNG" = "LUN",
    "OESOPHAGUS" = "Oe",
    "OVARY" = "Ov",
    "PANCREAS" = "Pa",
    "SKIN" = "Sk",
    "SOFT_TISSUE" = "STI",
    "STOMACH" = "STO",
    "THYROID" = "Th",
    "UPPER_AERODIGESTIVE_TRACT" = "UAT",
    "URINARY_TRACT" = "UT"
  )
    

#   mapping_Gutierrez <- NULL  #(CD4  CD8  iNKT NK  MAIT Vd1  Vd2 )

  switch(dataset,
         "Oihane" = mapping_Oihane,
         "mir" = mapping_mir,
         "Quake" = mapping_Quake,
         NULL)
}

GLOBAL_PARAMS <- list(
  n_components = 3,
  k_values = c(1, 2, 5, 10, 20),
  seed = 2024,
  figure_path = "./figures/",
  font_sizes = list(
    title = 36,
    axis = 22,
    axis_cci = 26,
    legend_title = 24,
    legend_text = 24,
      
    caption = 12
  ) 
)



color_list =c("#FB6A4A","#54278F","#006635","#3182BD","#DE2D26","#72A34F","#5D7AD3", "#756BB1","#FCAE91","#fe87ac","#AFABAB","#67A9CF","#CBC9E2","#4d982e","#E6873E","#545454","#aa3474","#ee8c7d","#2e5fa1","#FDD0A3","#C22F2F","#036f73")
names_list = c("PCA","MDS","iMDS","Sammon", "HLLE","Isomap","kPCA1","kPCA2","LEIM" , "UMAP1" , "UMAP2" , "tSNE1", "tSNE2", "PHATE1", "PHATE2","KEF")



# 设置通用主题
theme_custom <- function() {
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = GLOBAL_PARAMS$font_sizes$title, face = "bold"), ## ARI NMI CCI bold   visulization without
    axis.title = element_text(size = GLOBAL_PARAMS$font_sizes$axis),
    axis.text = element_text(size = GLOBAL_PARAMS$font_sizes$axis),
      axis.text.x = element_text(face = "bold"),   ## ARI NMI CCI bold   visulization without
       legend.key.size = unit(1.8, "cm"), # 调整图例键的大小

    legend.title = element_text(size = GLOBAL_PARAMS$font_sizes$legend_title),    # 调整图例标题的大小
    legend.text = element_text(size = GLOBAL_PARAMS$font_sizes$legend_text),  # 调整图例文本的大小
    panel.grid = element_blank(),
    legend.position = "",
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1.1)
  )
}


cal_cci <- function(ensemble.out, mev.out){
    n_components = 3
    cci_list = list()
    cci_viz<-dist.grp(ensemble.out$ensemble.dist.mat, info, k=c(1,2,5,10,20))
    cci_viz<-cbind(cci_viz, umap5(x=as.matrix(ensemble.out$ensemble.dist.mat), info=info,  n_components = n_components, do.plot=FALSE))
    cci_viz<-cbind(cci_viz, dist.grp(as.matrix(dist(eigen(max(ensemble.out$ensemble.dist.mat)-ensemble.out$ensemble.dist.mat)$vec[,1:n_components])), info, k=c(1,2,5,10,20))[,2])
    colnames(cci_viz)<-c("n_neighbors","viz_raw","viz_umap","viz_pca")
    cci_list[[1]] <- cci_viz
    cci_adm<-dist.grp(mev.out$diffu.dist, info, k=c(1,2,5,10,20))
    cci_adm<-cbind(cci_adm, umap5(x=mev.out$diffu.dist, info=info,  n_components = n_components, do.plot=FALSE))
    cci_adm<-cbind(cci_adm, dist.grp(as.matrix(dist(eigen(max(mev.out$diffu.dist)-mev.out$diffu.dist)$vec[,1:n_components])), info, k=c(1,2,5,10,20))[,2])
    colnames(cci_adm)<-c("n_neighbors","adm_raw","adm_umap","adm_pca")
    cci_list[[2]] <- cci_adm
    save(cci_list, file=paste("./Data/",dataset,"CCI.Rdata"))
    return (cci_list)
}



umap5<-function(x, info, do.plot=TRUE, n_components, k=c(1,2,5,10,20))
{
    library(uwot)

    for(n in 1:5)
    {
        ensemble.data<-umap(as.dist(x), n_components = n_components)
        ds<-dist.grp(as.matrix(dist(ensemble.data)), info, k=k)[,2]
        if(n==1) rec<-ds
        else rec<-cbind(rec,ds)
    }
    if(do.plot) pairs(ensemble.data, col=rainbow(length(unique(info)))[as.numeric(info)], pch=(1:length(unique(info)))[as.numeric(info)])
    to.return<-apply(rec,1,mean)
    return(to.return)
}

dist.grp<-function(distmat, grp, k=1:20)
{
    library(BiocParallel)
    grpmat<-matrix(0,nrow=length(grp),ncol=length(grp))
    for(i in 1:length(grp)) for(j in 1:length(grp)) if(grp[i]==grp[j]) grpmat[i,j]<-1

    rec<-cbind(k,k)
    diag(distmat)<-Inf
    
    dist.fun<-function(distmat, this.k, grpmat)
    {
        r<-NULL


        for(i in 1:nrow(distmat))
        {
            sel<-which(distmat[i,]<=quantile(distmat[i,], this.k/ncol(distmat)))
            r<-c(r,grpmat[i,sel])
        }
        r<-sum(r==1)/length(r)
        r
    }
    
    rec[,2]<-unlist(bplapply(k, dist.fun, distmat=distmat, grpmat=grpmat))
    rec
}



### dotplot Visualization 
visulization_fuc <- function(data, method_name,color_list,info){
    data <- as.data.frame(data)
#     print("Length of color_list:")
#     print(length(color_list))
    data$info = info
   plot <- ggplot(data, aes(x = data[[1]], y = data[[2]], color = info)) +
    geom_point(alpha = 0.75) +
    scale_color_manual(values = color_list) +
    labs(title = method_name,
         x = "",
         y = "") +
    theme_custom() +
     theme( legend.key.height = unit(1, "cm"),
            axis.text.x = element_blank(),  # 去掉 x 轴文字
            axis.text.y = element_blank(),  # 去掉 y 轴文字
            axis.ticks = element_blank(),    # 去掉刻度线
            axis.title.x = element_blank(), # 去掉 x 轴标题（刻度值）
            axis.title.y = element_blank()  # 去掉 y 轴标题（刻度值）
        )+
    guides(color = guide_legend(override.aes = list(size = 6, shape = 15, ncol = 2))) +
    theme(legend.key.height = unit(1, "cm"))
    

    target_dir <- "./figures/vizMethods"
    # 检查目录是否存在，如果不存在则创建
    if (!dir.exists(target_dir)) {
      dir.create(target_dir, recursive = TRUE)
    }
    ggsave(filename = file.path(target_dir, paste0(method_name, ".png")), plot = plot, width = 10, height = 8) 
  print(plot)
}



plot.bar.func <- function(data_long, dataset){
  p <- ggplot(data_long, aes(x = Metric, y = Value, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    labs(title = dataset,
         x = NULL, y = "") +
    scale_fill_manual(values = c("#4B8537", "#8669A9")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    theme_custom()
     ggsave(filename = paste0("./figures/arinmi.png"), plot = p, width = 6, height = 5)
    plot(p)
}





visualization_with_label <- function(umap_viz, method, info, dataset, color_list) {
  set.seed(2024)
  umap_df <- as.data.frame(umap_viz)
  umap_df$cluster <- info
  
  # 计算聚类中心和标签位置
  cluster_centers <- umap_df %>%
    group_by(cluster) %>%
    summarise(center_x = mean(V1),
              center_y = mean(V2)) %>%
    mutate(
      label_x = center_x + (center_x - mean(center_x)) * 0.2,
      label_y = center_y + (center_y - mean(center_y)) * 0.2
    )
  
  p <- ggplot(umap_df, aes(x = V1, y = V2, color = cluster)) +
    geom_point(alpha = 0.75) +
    scale_color_manual(values = color_list) +
    geom_segment(data = cluster_centers,
                 aes(x = center_x, y = center_y, xend = label_x, yend = label_y),
                 color = "grey50", linetype = "dashed") +
    geom_text_repel(data = cluster_centers, 
                    aes(x = label_x, y = label_y, label = cluster),
                    color = "black",
                    fontface = "bold",
                    box.padding = 0.5,
                    point.padding = 0.5,
                    segment.color = NA,size = 8) +
    labs(title = paste(method),
         x = "",
         y = "") +
    theme_custom() +
     theme( legend.key.height = unit(1, "cm"),
            axis.text.x = element_blank(),  # 去掉 x 轴文字
            axis.text.y = element_blank(),  # 去掉 y 轴文字
            axis.ticks = element_blank(),    # 去掉刻度线
            axis.title.x = element_blank(), # 去掉 x 轴标题（刻度值）
            axis.title.y = element_blank()  # 去掉 y 轴标题（刻度值）
        )+
    guides(color = guide_legend(override.aes = list(size = 6, shape = 15, ncol = 2))) +
    theme(legend.key.height = unit(0.8, "cm"))
 
  ggsave(filename = paste0("./figures/vizMethods/", method, "_labeled.png"), 
         plot = p, width = 10, height = 8)
  return(p)
}



plot_legend <- function(umap_viz, method_name, info, dataset, color_list, label_mapping = NULL) {
  set.seed(2024)
  # 准备数据
  umap_df <- as.data.frame(umap_viz)
  simplified_info <- simplify_labels(info, label_mapping)
  umap_df$cluster <- simplified_info
 p <- ggplot(umap_df, aes(x = umap_df[[1]], y = umap_df[[2]], color = cluster)) +
    geom_point(alpha = 0.75) +
    scale_color_manual(values = color_list) +
    labs(title = method_name,
         x = "",
         y = "") +
    theme_custom() +
    guides(color = guide_legend(override.aes = list(size = 8, shape = 15, ncol = 2))) +
     theme(legend.position = "right",  # 覆盖 theme_custom 中的 legend.position
#           legend.key.height = unit(1.3, "cm"),legend.spacing.x = unit(0.15, "cm"),
           legend.key.height = unit(1.3, "cm"),legend.spacing.x = unit(1, "cm"),
          legend.text = element_text(margin = margin(l = -0.5, unit = "cm")))
  ggsave(filename = paste0("./figures/vizMethods/legend.png"), plot = p, width = 10, height = 8)
  return(p)
}


plot_cci_results_cal <- function(cci_list, dataset) {  ## 如果没有保存CCI
  # 创建一个函数来处理每个矩阵并转换为长格式
  process_matrix <- function(matrix, method) {
    df <- as.data.frame(matrix)
    df$n_neighbors <- rownames(df)
    types <- colnames(df)[-1]
    values <- unlist(df[, -1])
    n_neighbors <- rep(df$n_neighbors, each = length(types))

    df_long <- data.frame(
      n_neighbors = n_neighbors,
      type = rep(types, times = nrow(df)),
      value = values,
      method = method
    )
    df_long$type <- sub(paste0(method, "_"), "", df_long$type)
    return(df_long)
  }
  viz_data <- process_matrix(cci_list[[1]], "viz")
  adm_data <- process_matrix(cci_list[[2]], "adm")
  all_data <- rbind(viz_data, adm_data)
  p <- ggplot(all_data, aes(x = type, y = value, fill = method)) +
    geom_boxplot() +
       scale_fill_manual(values = c("viz" = "#66C2A5", "adm" = "#FC8D62"))+
    labs(x = "Type", y = "Value", fill = "Method") +
    ggtitle(dataset) +
     them(
     axis.title = element_text(size = GLOBAL_PARAMS$font_sizes$axis_cci),
    axis.text = element_text(size = GLOBAL_PARAMS$font_sizes$axis_cci)
     )
    theme_custom()
   ggsave(filename = paste0("./figures/CCI.png"), plot = p, width = 12, height = 8)
  print(p)
}



plot_cci_results <- function(cci_list, dataset) {
  # 创建一个函数来处理每个矩阵并转换为长格式
  process_matrix <- function(matrix, method) {
    df <- as.data.frame(matrix)
    df$n_neighbors <- rownames(df)
    colnames(df) <- gsub(paste0(method, "_"), "", colnames(df))
    df_long <- pivot_longer(df, cols = c("raw", "umap", "pca"), 
                            names_to = "type", values_to = "value")
    df_long$method <- method
    return(df_long)
  }
  
  viz_data <- process_matrix(cci_list[[1]], "viz")
  diffu_data <- process_matrix(cci_list[[2]], "diffu")
  all_data <- rbind(viz_data, diffu_data)

  # 确保type的顺序
  all_data$type <- factor(all_data$type, levels = c("raw", "umap", "pca"))

  p <- ggplot(all_data, aes(x = type, y = value, fill = method)) +
    geom_boxplot() +
    scale_fill_manual(values = c("viz" = "#8669A9", "diffu" = "#4B8537"),    #c("viz" = "#7F7F7F", "diffu" = "#71C1A5")   绿灰配色
                      labels = c("viz" = "metaspec", "diffu" = "ADM")) +
    labs(x = " ", y = " ", fill = "Method") +
    ggtitle(dataset) +
    theme_custom() 
#     theme(legend.position = "right",  # 覆盖 theme_custom 中的 legend.position
#           legend.key.height = unit(1.2, "cm"))
  
  ggsave(filename = paste0("./figures/CCI.png"), plot = p, width = 12, height = 8)  #original 10
  print(p)
}



plot.SilhouetteWidth.func <- function(data_df, data_name, label_mapping = NULL){
    simplified_info <- simplify_labels(info, label_mapping)
    data_df$SimplifiedCluster <- simplified_info[match(data_df$Cluster, names(simplified_info))]
    write.csv(data_df, paste("./Data/",data_name,"ASW.csv"))
  p <- ggplot(data_df, aes(x = SimplifiedCluster, y = SilhouetteWidth, fill = Method)) + 
    scale_fill_manual(values = c("#4B8537", "#8669A9"))+
    geom_boxplot(position = position_dodge(0.8)) +
    labs(title = data_name,
         x = "",
         y = "Silhouette Width") +
       theme_custom()
#     theme(legend.position = "right",  # 覆盖 theme_custom 中的 legend.position
#           legend.key.height = unit(1.2, "cm"))
    ggsave(filename = paste0("./figures/Silhouette.png"), plot = p, width = 12, height = 8)
  return(p)
}



## 简化label
simplify_labels <- function(info, mapping = NULL) {
  if (is.null(mapping)) {
    simplified <- setNames(info, info)  
  } else {
    simplified <- mapping[as.character(info)]  
    simplified[is.na(simplified)] <- info[is.na(simplified)]
    return(simplified)
  }
}                
       


                         

  
#   Process each method in e
visualize_individual_methods <- function(e, names_list, seed = 2024) {
    for (i in 1:length(e)) {
        data <- data.frame(e[[i]])
        method <- names_list[[i]]
        set.seed(seed)
        visulization_fuc(data, method,color_list,info)
    }
}

                         
### calculating ARI and NMI
cal_ari_nmi <- function(lowDim_data, k, method_name, seed){
    set.seed(seed)
    cluster_viz <- kmeans(lowDim_data, centers = k)
    ARI <- adjustedRandIndex(info, cluster_viz$cluster)
    NMI <- NMI(info, cluster_viz$cluster)
    rec_ari_nmi <- data.frame(ARI = ARI, NMI = NMI)
    cat("******", method_name, "******\n")
    print(rec_ari_nmi) 
    return(rec_ari_nmi)
}

                         
# Process meta-spec method
process_and_visualize_meta_methods <- function(ensemble.out, mev.out, info, k,color_list, seed = 2024) {
  set.seed(seed)
  ARI_list <- list()
  ASW_list <- list()
  print(paste("Running R version:", R.version$major, ".", R.version$minor, sep = ""))
  method <- "meta-spec"
    
   print(find("umap"))
    if ("package:umap" %in% search()) {
  print("Using umap package")
  print(packageVersion("umap"))
} else if ("package:uwot" %in% search()) {
  print("Using uwot package")
  print(packageVersion("uwot"))
} else {
  print("Neither umap nor uwot package is loaded")
}

  umap_viz0 <- umap(ensemble.out$ensemble.dist.mat)
  

kmeans_source <- find("kmeans")
print(paste("kmeans function source:", kmeans_source))


if ("package:stats" %in% search()) {
  print("kmeans is available from stats package (R base)")
  print(paste("R version:", getRversion()))
} else if ("package:cluster" %in% search()) {
  print("cluster package is loaded, but note that it doesn't export a kmeans function")
  print(paste("cluster package version:", packageVersion("cluster")))
} else {
  print("Neither stats nor cluster package is explicitly loaded")
}


  set.seed(seed)
  cluster_viz <- kmeans(umap_viz0, centers = k)
  ARI <- adjustedRandIndex(info, cluster_viz$cluster)
  NMI <- NMI(info, cluster_viz$cluster)
  ARI_list[[1]] <- data.frame(ARI = ARI, NMI = NMI)
  print(ARI_list[[1]])

  ASW_list[[1]] <- silhouette(as.numeric(factor(info)), dist = dist(umap_viz0))[, 3]
  umap_viz <- as.data.frame(umap_viz0)
  visulization_fuc(umap_viz, method,color_list,info)
  
  # Process ADM method
  method <- "ADM"
   set.seed(seed)
  umap_adm0 <- umap(mev.out$diffu.dist) 
#     set.seed(2024)
  ARI_list[[2]] <- cal_ari_nmi(umap_adm0, k, method, seed)
  ASW_list[[2]] <- silhouette(as.numeric(factor(info)), dist = dist(umap_adm0))[, 3]
  umap_adm <- as.data.frame(umap_adm0)
  visulization_fuc(umap_adm, method,color_list,info)
  
  return(list(
    ARI_list = ARI_list,
    ASW_list = ASW_list,
    umap_viz = umap_viz0,
    umap_adm = umap_adm0
  ))
}
  

  
  # Plot Silhouette Width
 visualize_silhouette_width <- function(ASW_list, info, dataset,label_mapping = NULL) {
  # Extract silhouette width data
  asw_metaspec <- ASW_list[[1]]
  asw_adm <- ASW_list[[2]]
  
  # Create a data frame for plotting
  df_Oihane <- data.frame(
    Cluster = factor(info),
    SilhouetteWidth = c(asw_metaspec, asw_adm),
    Method = rep(c("meta-spec", "ADM"), each = length(asw_metaspec)),
    dataset = rep(dataset, each = length(info) * 2)
  )
  # Create and display the plot
  p1 <- plot.SilhouetteWidth.func(df_Oihane,dataset,label_mapping)
  print(p1)
  
}

                  
  # Plot ARI bar chart
visualize_ari_nmi <- function(ARI_list, dataset) {
  # Extract ARI and NMI data
  ARI_viz <- ARI_list[[1]]
  ARI_adm <- ARI_list[[2]]
  
  # Combine data
  data_combined <- rbind(
    cbind(ARI_viz, group = "meta-spec"),
    cbind(ARI_adm, group = "ADM")
  )
  
  # Convert to long format
  data_long <- data_combined %>%
    pivot_longer(cols = c(ARI, NMI), names_to = "Metric", values_to = "Value")
 
  # Create and display the plot
  plot.bar.func(data_long, dataset)
}




dataloader <- function(dataset, base_dir = ".") {
    path = paste0("/home/fanruibo/xiaobao/MEV/real data/",dataset,"/")
  setwd(path)
  
  load_common <- function() {
    load(paste0(path,"all embeded dim 3.bin"))
    load(file.path("Data", paste(dataset, "metaspec_out.Rdata")))  # ensemble.out
    load(file.path("Data", paste(dataset, "adm_out.Rdata")))  # mev.out
  }
  
  result <- switch(dataset,
    "Oihane" = {
      load("Oihane.Rdata")
      list(
        dat = Oihane,
        cell.type = as.factor(Oihane.info.cellType),
        info = as.factor(Oihane.info.cellType)
      )
    },
    "Gutierrez" = {
      load("Gutierrez.Rdata")
      list(
        dat = Gutierrez$data,
        cell.type = as.factor(Gutierrez$data.cellType),
        info = as.factor(Gutierrez$data.cellType)
      )
    },
     "Spleen" = {
       dat <- read.csv("./Spleen_pro_data.csv", header = TRUE, row.names = 1)
       lab <- read.csv("./Spleen_pro_label.csv", header = TRUE, row.names = 1)
#          print(lab)
       lab <- lab$SpatialGlue
#          print(lab)
      list(
         dat = dat,
         cell.type = lab,
         info = lab
      )
    },            
    "Quake" = {
      load("Quake_Smartseq2_Lung.Rdata")
      list(
        dat = Quake_Smartseq2_Lung$data,
        cell.type = as.factor(Quake_Smartseq2_Lung$data.cellType),
        info = as.factor(Quake_Smartseq2_Lung$data.cellType)
      )
    },
    "Brain5k" = {
      dat <- read.csv("Brain5k_data.csv", header = TRUE, row.names = 1)
      lab <- read.csv("Brain5k_metadata.csv", header = TRUE, row.names = 1)[,1]
      list(
        dat = dat,
        cell.type = as.factor(lab),
        info = as.factor(lab)
      )
    },
     "metabolism" = {
      dat <- read.csv("./metabolism_data.csv", header = TRUE, row.names = 1)
      lab <- read.csv("./metabolism_label.csv", header = TRUE, row.names = 1)
      lab<-lab$subtype  
      list(
        dat = dat,
        cell.type = as.factor(lab),
        info = as.factor(lab)
      )
    },
                        
    "kidney" = {
      dat <- read.csv("kidney_data.csv", header = TRUE, row.names = 1)
      lab <- read.csv("kidney_metadata.csv", header = TRUE, row.names = 1)$cell_type
      list(
        dat = dat,
        cell.type = as.factor(lab),
        info = as.factor(lab)
      )
    },
     "mir" = {load("CCLE_miRNA_20181103.bin")
    load("CCLE_RNAseq_genes_rpkm_20180929.bin")

    mir<-log10(1+mir)
    gene<-log10(1+gene)
    gene<-as.matrix(gene)
    mir<-as.matrix(mir)

    mir<-mir[,which(colnames(mir) %in% colnames(gene))]
    gene<-gene[,which(colnames(gene) %in% colnames(mir))]

    mir<-mir[,order(colnames(mir))]
    gene<-gene[,order(colnames(gene))]

    sum(colnames(mir)==colnames(gene))
    dim(mir)

    mir<-t(mir)
    gene<-t(gene)

    cv.mir<-apply(mir,2,sd)/apply(mir,2,mean)
    mir<-mir[,cv.mir>=0.1]

    ze<-apply(gene==0,2,sum)
    gene<-gene[,which(ze<=0.25*nrow(gene))]
    cv.gene<-apply(gene,2,sd)/apply(gene,2,mean)
    gene<-gene[,cv.gene>=0.5]

    ############################

    cells<-rownames(gene)
    cell.type<-cells
    for(i in 1:length(cells))
    {
        this<-strsplit(cells[i], "_")[[1]][-1]
        this<-paste(this, collapse="_")
        cell.type[i]<-this
    }

    ttt<-table(cell.type)
    sel<-which(cell.type %in% names(ttt)[ttt<10])
    cell.type[sel]<-NA

    anno<-read.table("Cell_lines_annotations_20181226.txt",header=T,sep="\t")
    anno<-anno[which(anno[,1] %in% cells),]
    anno<-anno[order(anno[,1]),]

    ##### limit to cells with annotations
    s.anno<-which(cells %in% anno[,1] & !is.na(cell.type))
    cells<-cells[s.anno]
    cell.type<-cell.type[s.anno]
    mir<-mir[s.anno,]
    gene<-gene[s.anno,]

    mir2<-mir
    gene2<-gene

    for(i in 1:ncol(mir2)) mir2[,i]<-(mir2[,i]-mean(mir2[,i]))/sd(mir2[,i])
    for(i in 1:ncol(gene2)) gene2[,i]<-(gene2[,i]-mean(gene2[,i]))/sd(gene2[,i])

    all.colors<-c("grey50","green","blue","cyan", "yellow","orange","red","black", "wheat","purple","darkblue","dodgerblue4","darkred","darkorange","darkcyan","magenta","firebrick","khaki4","cornsilk3","darkgoldenrod4")

    ################### 
    dat<-mir
    info<-as.factor(cell.type)
 
     list(
        dat = dat,
        cell.type = cell.type,
        info = info
      )
    },
     "gene" = {load("CCLE_miRNA_20181103.bin")
    load("CCLE_RNAseq_genes_rpkm_20180929.bin")

    mir<-log10(1+mir)
    gene<-log10(1+gene)
    gene<-as.matrix(gene)
    mir<-as.matrix(mir)

    mir<-mir[,which(colnames(mir) %in% colnames(gene))]
    gene<-gene[,which(colnames(gene) %in% colnames(mir))]

    mir<-mir[,order(colnames(mir))]
    gene<-gene[,order(colnames(gene))]

    sum(colnames(mir)==colnames(gene))
    dim(mir)

    mir<-t(mir)
    gene<-t(gene)

    cv.mir<-apply(mir,2,sd)/apply(mir,2,mean)
    mir<-mir[,cv.mir>=0.1]

    ze<-apply(gene==0,2,sum)
    gene<-gene[,which(ze<=0.25*nrow(gene))]
    cv.gene<-apply(gene,2,sd)/apply(gene,2,mean)
    gene<-gene[,cv.gene>=0.5]

    ############################

    cells<-rownames(gene)
    cell.type<-cells
    for(i in 1:length(cells))
    {
        this<-strsplit(cells[i], "_")[[1]][-1]
        this<-paste(this, collapse="_")
        cell.type[i]<-this
    }

    ttt<-table(cell.type)
    sel<-which(cell.type %in% names(ttt)[ttt<10])
    cell.type[sel]<-NA

    anno<-read.table("Cell_lines_annotations_20181226.txt",header=T,sep="\t")
    anno<-anno[which(anno[,1] %in% cells),]
    anno<-anno[order(anno[,1]),]

    ##### limit to cells with annotations
    s.anno<-which(cells %in% anno[,1] & !is.na(cell.type))
    cells<-cells[s.anno]
    cell.type<-cell.type[s.anno]
    mir<-mir[s.anno,]
    gene<-gene[s.anno,]

    mir2<-mir
    gene2<-gene

    for(i in 1:ncol(mir2)) mir2[,i]<-(mir2[,i]-mean(mir2[,i]))/sd(mir2[,i])
    for(i in 1:ncol(gene2)) gene2[,i]<-(gene2[,i]-mean(gene2[,i]))/sd(gene2[,i])

    all.colors<-c("grey50","green","blue","cyan", "yellow","orange","red","black", "wheat","purple","darkblue","dodgerblue4","darkred","darkorange","darkcyan","magenta","firebrick","khaki4","cornsilk3","darkgoldenrod4")

    ################### 
    dat<-gene
    info<-as.factor(cell.type)
 
     list(
        dat = dat,
        cell.type = cell.type,
        info = info
      )
    },
    stop("Unsupported dataset")
  )
  

 
  return(result)
}


mev<-function(e, k.dim=NULL, dist.power=0.5, conn.prop=0.003, raw.d.pwr=.5, diffu.steps=NA, diffu.factor=3.5, distr.template="combine", gamma.shape=3, gamma.rate=3, scale.dist=TRUE, symmetrize="mean", dist.quantile=0.25)
{
	########################
	if(is.null(k.dim)) k.dim=ncol(e[[1]])

	#this.e is one dimension reduction result
	print("working on mev")
	fake.fun<-function(this.e, raw.d.pwr, diffu.steps, conn.prop, scale.dist, symmetrize, diffu.factor, dist.quantile)
	{
	
		library(DDoutlier)
		library(diffudist)
		library(MixGHD)
		library(BiocParallel)
		library(fitdistrplus)
		
		zp.quantile <- function(x,y)
		{
		  o.x<-order(x)

		  o.y<-order(y)
		  r.y<-rank(y, ties.method = "average")

		  x<-x[o.x]
		  y<-y[o.y]


		  z.x<-seq(0,1,length.out=length(x))
		  z.y<-seq(0,1,length.out=length(y))

		  new.y<-stats::approx(x=z.x, y=x, xout=z.y)$y
		  #y[y>0]<-new.y2
		  y<-new.y[r.y]
		}

		estimate_mode <- function(x) {
		  d <- density(x)
		  d$x[which.max(d$y)]
		}
	
	
		move.outlier<-function(x, d=NULL, fraction=0.01)
		{
			x.is.vector<-F
			if(is.null(nrow(x))) 
			{
				x<-matrix(x, ncol=1)
				x.is.vector<-TRUE
			}
			
			for(n in 1:ncol(x))
			{
				this.x<-x[,n]
				if(is.null(d))
				{
					this.d<-diff(quantile(this.x,c(0.25,0.75)))/3
				}else{
					this.d<-d
				}
				
				this.out<-DB(matrix(this.x,ncol=1), this.d, fraction)
				this.sel<-which(this.out$classification == "Outlier")
				
				if(length(this.sel)>0)
				{
					new.x<-zp.quantile(this.x[-this.sel], this.x)
					x[,n]<-new.x
				}
			}
			x
		}
		
		this.e<-move.outlier(this.e)
		this.d<-as.matrix(dist(this.e))

		n<-nrow(this.d)
		
		diag(this.d)<-Inf
		
		######## global
		
		conn.cutoff<-quantile(as.dist(this.d), conn.prop)
		this.conn<-1*(this.d <= conn.cutoff)

		########## local
		
		for(i in 1:nrow(this.d))
		{
			sel<-which(this.d[i,]<=quantile(this.d[i,], conn.prop))
			this.conn[i,sel]<-1
		}

		this.graph<-graph_from_adjacency_matrix(this.conn, mode="undirected")

		
		pmat<-1/this.d^raw.d.pwr
		diag(pmat)=0
		pmat[pmat == Inf]<-max(pmat[pmat != Inf], na.rm=TRUE)
		pmat[is.na(pmat)]<-max(pmat, na.rm=TRUE)
		
		pmat<-pmat*this.conn
		for(i in 1:nrow(pmat)) pmat[i,]<-pmat[i,]/sum(pmat[i,])
#      write.csv(pmat,"./Data/pmat.csv")
		if(is.na(diffu.steps[1])) 
		{
			sp.mat<-shortest.paths(this.graph)
			sp.mat[sp.mat==Inf]<-NA
			mean.steps<-apply(sp.mat,1,quantile, na.rm=TRUE, probs=dist.quantile)
			diffu.steps<-quantile(mean.steps, c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm=TRUE)*diffu.factor
		}
		
		# sp.mat is shortest distance
		# diffu.steps is all the steps we try
		
		all.d<-new("list")
		for(i in 1:length(diffu.steps))
		{
			this.d<-as.matrix(suppressMessages(get_distance_matrix_from_T(pmat, diffu.steps[i])))
			if(scale.dist) 
			{
				this.d<- 1-cor(this.d, method="spearman")
				this.d[this.d>1]<-1
			}
       
			all.d[[i]]<-this.d
		}
		
		# merge.d is the matrix in which each point takes different diffusion time
		# the steps for each point is based on its mean distance to other points
		
		merge.d<-all.d[[1]]
#      write.csv(pmat,"./Data/pmat.csv")
		for(i in 1:nrow(merge.d))
		{
			this.dist<-abs(mean.steps[i]*diffu.factor-diffu.steps)
			this.closest<-which(this.dist==min(this.dist))[1]
			merge.d[i,]<-all.d[[this.closest]][i,]
		}
		this.d<-(merge.d+t(merge.d))/2
		
		if(symmetrize=="mean") {
			this.d<-this.d+t(this.d)
		}else if(symmetrize=="min") {
			this.d.2<-t(this.d)
			this.d[this.d>this.d.2]<-this.d.2[this.d>this.d.2]
		}else if(symmetrize=="max") {
			this.d.2<-t(this.d)
			this.d[this.d<this.d.2]<-this.d.2[this.d<this.d.2]			
		}

		to.return<-list(this.d=this.d)
		return(to.return)
	}

	# e is a list object. Each item is a dimension reduction result, an nXk matrix. K is the dimension.
	
	d<-bplapply(e, fake.fun, raw.d.pwr=raw.d.pwr, diffu.steps=diffu.steps, conn.prop=conn.prop, scale.dist=scale.dist, symmetrize=symmetrize, diffu.factor=diffu.factor, dist.quantile=dist.quantile)
	
	################### for Junning to change #######################
	
	# the code below combines the diffusion distance matrices. A quantile normalization is taken
	# combine: merge all into a long vector as template
	# gamma: find gamma parameter from each, and take average of the parameters
	# parametric: specify gamma parameters by user
	
	if(distr.template != "none")
	{
		
		if(distr.template == "combine")
		{
			d.template<-NULL
			for(i in 1:length(d)) d.template<-c(d.template, as.dist(d[[i]][[1]]))
		}else if(distr.template == "gamma"){
			for(m in 1:length(d))
			{
				this.fit<-fitdist(as.vector(as.dist(d[[m]][[1]])),"gamma",method="mme")$estimate
				if(m == 1) all.fit=this.fit
				else all.fit=rbind(all.fit, this.fit)
			}
			ave.fit= apply(all.fit,2,mean)
			
			d.template<-rgamma(nrow(d[[1]][[1]])*ncol(d[[1]][[1]]), shape=ave.fit[1], rate=ave.fit[2])
		}else if(distr.template=="parametric"){
			d.template<-rgamma(nrow(d[[1]][[1]])*ncol(d[[1]][[1]]), shape=gamma.shape, rate=gamma.rate)
		}

		for(m in 1:length(d))
		{
			this.d<-as.dist(d[[m]][[1]])
			this.d2<-zp.quantile(d.template, this.d)
			attributes(this.d2)<-attributes(this.d)			
			this.d2<-as.matrix(this.d2)
			d[[m]][[1]]<-this.d2
		}
	}		
	for(m in 1:length(d)){
#     write.csv(d[[m]],paste("./Data/diffu_dist",m, ".csv"))
    }
	dd<-d[[1]][[1]]*0
	for(m in 1:length(d)) dd<-dd+d[[m]][[1]]^dist.power

	to.return=list(diffu.dist=dd)
#    write.csv(dd,"./Data/diffu_dist_merged.csv")
	return(to.return)
}




dist.grp<-function(distmat, grp, k=1:20)
{
	library(BiocParallel)
	grpmat<-matrix(0,nrow=length(grp),ncol=length(grp))
	for(i in 1:length(grp)) for(j in 1:length(grp)) if(grp[i]==grp[j]) grpmat[i,j]<-1

	rec<-cbind(k,k)
	diag(distmat)<-Inf
	
	dist.fun<-function(distmat, this.k, grpmat)
	{
		r<-NULL
		
		
		for(i in 1:nrow(distmat))
		{
			sel<-which(distmat[i,]<=quantile(distmat[i,], this.k/ncol(distmat)))
			r<-c(r,grpmat[i,sel])
		}
		
		r<-sum(r==1)/length(r)
		r
	}
	
	rec[,2]<-unlist(bplapply(k, dist.fun, distmat=distmat, grpmat=grpmat))
	rec
}



candidate.visual <- function(data, dim=2, methods= c("PCA", "MDS", "iMDS", "Sammon", "LLE", "HLLE","Isomap", 
                                                     "kPCA", "LEIM", "UMAP", "tSNE", "PHATE"), 
                             kpca.sigma = c(0.001, 0.002), 
                             umap.k= c(30, 50), 
                             tsne.perplexity = c(30, 50),
                             phate.k = c(30,50),
                             cal_dist = TRUE){

  n=dim(data)[1]
  dim.red.data = list()
  if(cal_dist){
    dist.data = dist(data)
  }
  print(dim(data))
  name.method =c()
  ####################
  ###### PCA
  ####################
  
  #Original in meta-spec
  i=0
  if(sum(methods == "PCA")>0){
    i=i+1
    print("PCA calculating...")
    pc.data = rARPACK::svds(as.matrix(data), k =dim)
    dim.red.data[[i]] = pc.data$u[,1:2]
    name.method =  c(name.method, "PCA")
  }

    
  #####################
  ####### classical MDS
  #####################
  
  if(sum(methods == "MDS")>0){
    i=i+1
    print("MDS calculating...")
    dim.red.data[[i]] = cmdscale(dist.data, k=dim)
    name.method =  c(name.method, "MDS")
  }
  #####################
  ####### isoMDS
  #####################
  
  if(sum(methods == "iMDS")>0){
    i=i+1
     print("iMDS calculating...")
    imds.data = isoMDS(dist.data, k=dim)
    dim.red.data[[i]] = imds.data$points
    name.method =  c(name.method, "iMDS")
  }
  
  #####################
  ####### Sammon's nonlinear mapping
  #####################
  
  if(sum(methods == "Sammon")>0){
    i=i+1
    print("Sammon calculating...")
    sam.data = sammon(dist.data, k=dim)
    dim.red.data[[i]] = sam.data$points
    name.method =  c(name.method, "Sammon")
  }
  
  #####################
  ####### LLE
  #####################
  
  if(sum(methods == "LLE")>0){
    i=i+1
     print("LLE calculating...")
    #k.sel = calc_k(data, m = 2)
    lle.data = lle(data, m=dim, k=20, reg=2)
    dim.red.data[[i]] =  lle.data$Y
    name.method =  c(name.method, "LLE")
  }
  
  #####################
  ####### HLLE
  #####################
  
  if(sum(methods == "HLLE")>0){
    i=i+1
     print("HLLE calculating...")
    hlle.data <- embed(data, "HLLE", knn =20, ndim=dim)
    dim.red.data[[i]] = hlle.data@data@data
    name.method =  c(name.method, "HLLE")
  }
  
  #####################
  ####### isomap
  #####################
  
  if(sum(methods == "Isomap")>0){
    i=i+1
     print("Isomap calculating...")
    imp.data <- embed(data, "Isomap", knn = 20, ndim=dim)
    dim.red.data[[i]] = imp.data@data@data
    name.method =  c(name.method, "Isomap")
  }
  
  #####################
  ####### kPCA
  #####################
  
  if(sum(methods == "kPCA")>0){
    for(j in 1:length(kpca.sigma)){
      i=i+1
      print("kPCA calculating...")
      kpca.data <- embed(data,  "kPCA", kpar=list(sigma=kpca.sigma[j]), features=dim, ndim=dim)
      dim.red.data[[i]] = kpca.data@data@data
      name.method =  c(name.method, paste0("kPCA", j))
    }
  }
  
  #####################
  ####### Laplacian Eigenmap
  #####################
  
  if(sum(methods == "LEIM")>0){
    i=i+1
     print("LEIM calculating...")
    lem.data <- embed(data,  "LaplacianEigenmaps", ndim=dim)
    dim.red.data[[i]] = lem.data@data@data
    name.method =  c(name.method, "LEIM")
  }
  
  #####################
  ####### UMAP
  #####################
  
  if(sum(methods == "UMAP")>0){
    for(j in 1:length(umap.k)){
      i=i+1
      print("UMAP calculating...")
      umap.data <- umap(data,  n_neighbors = umap.k[j], n_components = dim)
      dim.red.data[[i]] = umap.data
      name.method =  c(name.method, paste0("UMAP", j))
    }
  }
  #####################
  ####### tSNE
  #####################
  
  if(sum(methods == "tSNE")>0){
    for(j in 1:length(tsne.perplexity)){
      i=i+1
      print("tSNE calculating...")
      tsne.data <- embed(data,  "tSNE", perplexity = tsne.perplexity[j], ndim=dim)
      dim.red.data[[i]] = tsne.data@data@data
      name.method =  c(name.method, paste0("tSNE", j))
    }
  }
  
  #####################
  ####### PHATE
  #####################
  
  if(sum(methods == "PHATE")>0){
    for(j in 1:length(phate.k)){
      i=i+1
      print("PHATE calculating...")
      dim.red.data[[i]] = phate(data, knn=phate.k[j], ndim=dim)$embedding
      name.method =  c(name.method, paste0("PHATE", j))
    }
  }
  
  #####################
  ####### KEF
  #####################
  if(sum(methods == "KEF") > 0){
    i = i + 1
    print("KEF calculating...")
    dist.mat = dist(data)
    K.mat = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.5)^2)
    eigen.K.mat = eigen(K.mat)
    kef_u1 = eigen.K.mat$vectors[,2] * eigen.K.mat$values[2]
    kef_u2 = eigen.K.mat$vectors[,3] * eigen.K.mat$values[3]
    print(dim(kef_u1))
    dim.red.data[[i]] = cbind(kef_u1, kef_u2)
    name.method =  c(name.method, "KEF")
  }
  
   
  
  
  return(list(embed.list=dim.red.data, method_name=name.method))
}


#################################################
######## meta-visualization function
#################################################

# ensemble.viz: using spectral method for quantifying visualization quality and generating meta-visualization.
# data.list: a list of 2-dimensionoal embeddings, which is created by candidate.visual.
# name.method: the names of the candidate visualizations.
# original.data: an option to use the original data for the quality assessment, instead of using eigenscores.
# Output: a list containing (1) a meta-distance for meta-visualization, (2) eigenscores for candidate visualizations,
#         (3) names of the method ordered by averaged eigenscores.



#this was used for running time evaluation for n>10000 data
ensemble.viz <- function(data.list, name.method=NA, original.data=NA){
  
  
  n=dim(data.list[[1]])[1]
  K=length(data.list)
  
  if(is.na(original.data)){
    ########## obtain weights
    
    ensemble.mat = matrix(ncol=n,nrow=n)
    weight = matrix(ncol=K,nrow=n)
    for(j in 1:n){
      
      local.dist=matrix(ncol=n,nrow=K)
      for(i in 1:K){
        local.dist[i,]= sqrt(rowSums(t(data.list[[i]][j,]-t(data.list[[i]]))^2))
      }
      
      comp.mat = matrix(ncol=K, nrow=K)
      embed.mat.norm = list()
      for(i in 1:K){
        for(k in 1:K){
          comp.mat[i,k] = sum(local.dist[i,]*local.dist[k,])/sqrt(sum(local.dist[k,]^2))/sqrt(sum(local.dist[i,]^2))
        }
        embed.mat.norm[[i]] = local.dist[i,]/sqrt(sum(local.dist[i,]^2))
      }
      weight[j,] = abs(eigen(comp.mat)$vectors[,1])
      
      
      ensemble.mat[,j] = apply(do.call(cbind, embed.mat.norm), 1, weighted.mean, w = weight[j,])*sum(weight[j,])
      if(j/1000==floor(j/1000)){
        print(paste0(j," samples done!"))
      }
    }
    
    
    return(list(ensemble.dist.mat=(ensemble.mat+t(ensemble.mat))/2, eigenscore=weight, 
                name.method=name.method[order(colMeans(weight), decreasing = T)]))
    
  }else{
    
    ensemble.mat = matrix(ncol=n,nrow=n)
    weight = matrix(ncol=K,nrow=n)
    data.mat= as.matrix(dist(original.data))
    for(j in 1:n){
      
      local.dist=matrix(ncol=n,nrow=K)
      for(i in 1:K){
        local.dist[i,]= sqrt(rowSums(t(data.list[[i]][j,]-t(data.list[[i]]))^2))
      }
      
      eigen.score = c()
      embed.mat.norm = list()
      for(i in 1:K){
        eigen.score[i]=sum(data.mat[,j]*local.dist[i,])/sqrt(sum(data.mat[,j]^2))/sqrt(sum(local.dist[i,]^2))
        embed.mat.norm[[i]] = local.dist[i,]/sqrt(sum(local.dist[i,]^2)) 
      }
      
      eigen.score[which(eigen.score<0)]=0
      weight[j,] = eigen.score^2/sum(eigen.score^2)
      
      ensemble.mat[,j] = apply(do.call(cbind, embed.mat.norm), 1, weighted.mean, w = weight[j,])*sum(weight[j,])
      
    }
    
    return(list(ensemble.dist.mat=(ensemble.mat+t(ensemble.mat))/2, eigenscore=weight, 
                name.method=name.method[order(colMeans(weight), decreasing = T)]))
    
  }  
}


ensemble.v.local <- function(data.list, name.method=NA, original.data=NA){

  
  n=dim(data.list[[1]])[1]
  K=length(data.list)
  
  
  dist.list=array(dim=c(K,n,n))
  for(i in 1:K){
    dist.list[i,,]=as.matrix(dist(data.list[[i]]))
  }
  
  if(is.na(original.data)){
    ########## obtain weights
    
    ensemble.mat = matrix(ncol=n,nrow=n)
    weight = matrix(ncol=K,nrow=n)
    embed.mat.norm = list()
    for(j in 1:n){
      comp.mat = matrix(ncol=K, nrow=K)
      for(i in 1:K){
        embed.mat.norm[[i]] = dist.list[i,,j]/sqrt(sum(dist.list[i,,j]^2))
        for(k in 1:K){
          comp.mat[i,k] = sum(dist.list[k,,j]*dist.list[i,,j])/sqrt(sum(dist.list[k,,j]^2))/sqrt(sum(dist.list[i,,j]^2))
        }
      }
      weight[j,] = abs(eigen(comp.mat)$vectors[,1])
      
      ensemble.mat[,j] = apply(do.call(cbind, embed.mat.norm), 1, weighted.mean, w = weight[j,])*sum(weight[j,])
      if(j/1000==floor(j/1000)){
        print(j)
      }
    }
    
    
    return(list(ensemble.dist.mat=(ensemble.mat+t(ensemble.mat))/2, eigenscore=weight, 
                name.method=name.method[order(colMeans(weight), decreasing = T)], dist.list=dist.list))
    
  }else{
    
    ensemble.mat = matrix(ncol=n,nrow=n)
    weight = matrix(ncol=K,nrow=n)
    data.mat= as.matrix(dist(original.data))
    for(j in 1:n){
      
      eigen.score = c()
      embed.mat.norm = list()
      for(i in 1:K){
        eigen.score[i]=sum(data.mat[,j]*dist.list[i,,j])/sqrt(sum(data.mat[,j]^2))/sqrt(sum(dist.list[i,,j]^2))
        embed.mat.norm[[i]] = dist.list[i,,j]/sqrt(sum(dist.list[i,,j]^2)) 
      }
      
      eigen.score[which(eigen.score<0)]=0
      weight[j,] = eigen.score^2/sum(eigen.score^2)
      
      ensemble.mat[,j] = apply(do.call(cbind, embed.mat.norm), 1, weighted.mean, w = weight[j,])*sum(weight[j,])
      
    }
    
    return(list(ensemble.dist.mat=(ensemble.mat+t(ensemble.mat))/2, eigenscore=weight, 
                name.method=name.method[order(colMeans(weight), decreasing = T)], dist.list=dist.list))
    
  }

}
