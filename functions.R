library(ggplot2)
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





dataloader <- function(dataset, base_dir = ".") {
    path = paste0("../dataset/")
    setwd（path）
  result <- switch(dataset,
    "Oihane" = {
      dat <- read.csv("./Oihane_data.csv", header = TRUE, row.names = 1)
      lab <- read.csv("./Oihane_metadata.csv", header = TRUE, row.names = 1)
      lab <- lab$cell_type1
      list(
         dat = dat,
         cell.type = lab,
         info = lab
      )
    },
    "Gutierrez" = {
      dat <- read.csv("./Gutierrez_data.csv", header = TRUE, row.names = 1)
      lab <- read.csv("./Gutierrez_metadata.csv", header = TRUE, row.names = 1)
      lab <- lab$cell_type1
      list(
         dat = dat,
         cell.type = lab,
         info = lab
      )
    },
     "Spleen" = {
       dat <- read.csv("./Spleen_pro_data.csv", header = TRUE, row.names = 1)
       lab <- read.csv("./Spleen_pro_label.csv", header = TRUE, row.names = 1)
       lab <- lab$SpatialGlue
      list(
         dat = dat,
         cell.type = lab,
         info = lab
      )
    },            
    "Quake" = {
      dat <- read.csv("./Quake_data.csv", header = TRUE, row.names = 1)
      lab <- read.csv("./Quake_metadata.csv", header = TRUE, row.names = 1)
      lab <- lab$cell_type1
      list(
         dat = dat,
         cell.type = lab,
         info = lab
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



cal_ari_nmi <- function(lowDim_data, k, method_name, seed, info = NULL){
  if (is.null(info)) {
    if (!exists("info")) {
      stop("'info' not provided and not found in global environment")
    }
    info <- get("info", envir = .GlobalEnv)
  }
  set.seed(seed)
  cluster_viz <- stats::kmeans(lowDim_data, centers = k)
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("aricode", quietly = TRUE)) {
    stop("Package 'aricode' is needed for this function to work. Please install it.", call. = FALSE)
  }
  ARI <- mclust::adjustedRandIndex(info, cluster_viz$cluster)
  NMI <- aricode::NMI(info, cluster_viz$cluster)
  rec_ari_nmi <- data.frame(ARI = ARI, NMI = NMI)
  cat("******", method_name, "******\n")
  print(rec_ari_nmi)
  return(rec_ari_nmi)
}
visualization_func <- function(data, method_name, color_list = NULL, info = NULL) {
  # Convert data to data frame and set column names
  data <- as.data.frame(data)
  colnames(data)[1:2] <- c("x", "y")
  data$info <- info
  
  # Create the base plot
  plot <- ggplot(data, aes(x = .data$x, y = .data$y)) +
    labs(title = method_name, x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"  # Remove the legend
    )

  # Add points to the plot
  if (!is.null(info) && !is.null(color_list)) {
    data$info <- info
    plot <- plot +
      geom_point(aes(color = .data$info), alpha = 0.75) +
      scale_color_manual(values = color_list)
  } else {
    plot <- plot + geom_point(alpha = 0.75)
  }

  return(plot)
}  

process_and_visualize_meta_methods <- function(mev.out, ensemble.out = NULL, info, k, color_list, seed = 2024) {
  # Input validation
  if (!is.null(ensemble.out) && !all(c("ensemble.dist.mat") %in% names(ensemble.out))) {
    stop("When provided, ensemble.out must contain 'ensemble.dist.mat'")
  }
  if (!all(c("diffu.dist") %in% names(mev.out))) {
    stop("mev.out must contain 'diffu.dist'")
  }
  if (!is.null(ensemble.out) && length(info) != nrow(ensemble.out$ensemble.dist.mat)) {
    stop("Length of info must match the number of rows in ensemble.dist.mat")
  }


  ARI_list <- list()
  ASW_list <- list()
  plots <- list()
  print(paste("Running R version:", R.version$major, ".", R.version$minor, sep = ""))
  # set.seed(seed)
  # Meta-spec method (only if ensemble.out is not NULL)
  if (!is.null(ensemble.out)) {
    method <- "meta-spec"
    umap_viz0 <- uwot::umap(ensemble.out$ensemble.dist.mat)
    ARI_list[[1]] <- cal_ari_nmi(umap_viz0, k, method, seed, info)
    ASW_list[[1]] <- cluster::silhouette(as.numeric(factor(info)), dist = stats::dist(umap_viz0))[, 3]
    umap_viz <- as.data.frame(umap_viz0)
    plots[[1]] = visualization_func(umap_viz, method, color_list, info)
  }

  # ADM method
  method <- "ADM"
  set.seed(seed)
  umap_adm0 <- uwot::umap(mev.out$diffu.dist)
  ARI_list[[length(ARI_list) + 1]] <- cal_ari_nmi(umap_adm0, k, method, seed, info)
  ASW_list[[length(ASW_list) + 1]] <- cluster::silhouette(as.numeric(factor(info)), dist = stats::dist(umap_adm0))[, 3]
  umap_adm <- as.data.frame(umap_adm0)
  plots[[length(ASW_list) + 1]] = visualization_func(umap_adm, method, color_list, info)

  result <- list(
    ARI_list = ARI_list,
    ASW_list = ASW_list,
    umap_adm = umap_adm0,
    plot = plots
  )

  if (!is.null(ensemble.out)) {
    result$umap_viz <- umap_viz0
  }

  return(result)
}


