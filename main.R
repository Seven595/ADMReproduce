

source("./plot_functions.R")

dataset = "Oihane"  ## Gutierrez  Oihane Quake  pbmc  mir  kidney  Spleen metabolism  gene

dataload = dataloader(dataset)
dat = dataload$dat
info = dataload$info
path = paste0("../dataset/",dataset,"/")

##########################################
##
## Calculating  indivadual outputs
##
##########################################
# candidate.out = candidate.visual(dat, dim = 3, method=c("PCA", "MDS", "iMDS", "Sammon", "HLLE", "Isomap", 
#                                                   "kPCA", "LEIM", "UMAP", "tSNE","PHATE","KEF"),tsne.perplexity = c(10, 30))
# e<-candidate.out[[1]]
# name = candidate.out[[2]]
# save(e, file=paste0(path,"all embeded dim 3.bin"))  


##########################################
##
##  Calculating  meta outputs
##
##########################################
# load(paste0(path,"all embeded dim 3.bin"))
# ensemble.out = ensemble.viz(e, names(e))
# mev.out = mev(e, dist.power=dist.power, conn.prop=conn.prop,diffu.factor = diffu.factor)
# CCI <- cal_cci(ensemble.out, mev.out)
# save(ensemble.out, file=paste("./Data/",dataset,"metaspec_out.Rdata"))
# save(mev.out, file=paste("./Data/",dataset, "adm_out.Rdata"))
# save(CCI, file=paste("./Data/",dataset,"CCI.Rdata"))


##########################################
##
##  visulization
##
##########################################


load(paste0(path,"all embeded dim 3.bin"))
load(paste0(path, "Data/", " ", dataset, " metaspec_out.Rdata"))
load(paste0(path,"Data/ ", dataset, " adm_out.Rdata"))
load(paste0(path,"Data/ ", dataset, " CCI.Rdata"))
k = length(unique(info))
label_mapping <- get_mapping(dataset)
visualize_individual_methods(e, names_list, seed = 2024)
results = process_and_visualize_meta_methods(ensemble.out, mev.out, info, k, color_list,seed = 2024)
# plot_legend(results[[4]], "ADM", info, dataset, color_list, label_mapping = label_mapping) ##从该屠图上拿legend


visualization_with_label(results[[4]], "ADM", info, dataset, color_list)
visualization_with_label(results[[3]], "meta-spec", info, dataset, color_list)  
visualize_silhouette_width(results[[2]], info,dataset,label_mapping)
visualize_ari_nmi(results[[1]], dataset)
plot_cci_results(rec, dataset)

