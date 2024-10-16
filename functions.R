
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
    path = paste0("./dataset/")
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
       lab <- lab$SpatialGlue
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

