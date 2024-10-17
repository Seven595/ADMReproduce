# Load necessary libraries
library(ADM)
source("./functions.R")

# Define color list for visualization
color_list <- c("#FB6A4A", "#54278F", "#006635", "#3182BD", "#DE2D26", "#72A34F", "#5D7AD3", 
                "#756BB1", "#FCAE91", "#fe87ac", "#AFABAB", "#67A9CF", "#CBC9E2", "#4d982e", 
                "#E6873E", "#545454", "#aa3474", "#ee8c7d", "#2e5fa1", "#FDD0A3", "#C22F2F", "#036f73")

# Set dataset name
dataset <- "Quake"  # Options: Gutierrez, Oihane, Quake, Brain5k, mir, Spleen, metabolism, gene

# Load data
dataload <- dataloader(dataset)
dat <- dataload$dat  # Data matrix
info <- dataload$info  # Label information
k <- length(unique(info))  # Number of classes
label_mapping <- get_mapping(dataset)  # Get label mapping
print("Data loaded successfully!")

# Set data path
path <- file.path("../dataset", dataset)

# Set random seed for reproducibility
set.seed(2024)

# Execute candidate visualization methods
candidate.out <- candidate.visual(
  dat, 
  dim = 3, 
  method = c("PCA", "MDS", "iMDS", "Sammon", "HLLE", "Isomap", "kPCA", "LEIM", "UMAP", "tSNE", "PHATE", "KEF"),
  tsne.perplexity = c(10, 30)
)
print("Individual methods completed!")

# Extract results
e <- candidate.out[[1]]  # Visualization results
names_list <- candidate.out[[2]]  # List of method names

# Execute ensemble visualization
ensemble.out <- ensemble.viz(e, names(e))
print("Meta-spec completed!")

# Execute ADM method
adm.out <- adm(e, distr.template = "combine")
print("ADM completed!")

# Process and visualize meta-method results
result <- process_and_visualize_meta_methods(adm.out, ensemble.out, info, k, color_list)

# Calculate the Category consistency index (CCI)
cci = cal_cci(ensemble.out, adm.out, info)
print(cci)
# Visualize individual method results
ind_result <- visualize_individual_methods(e, names_list, info, color_list, k)

# View results
# Use the following commands to view numerical results and plots for each method:
# ind_result[[1]]$plot  # Plot
# ind_result[[1]]$ari   # ARI value
# ind_result[[1]]$nmi   # NMI value
# ind_result[[1]]$silhouette  # Silhouette coefficient

# Output results summary
print("Results summary:")
for (i in seq_along(ind_result)) {
  method_name <- names(ind_result)[i]
  cat(sprintf("%s: ARI = %.4f, NMI = %.4f, Silhouette = %.4f\n",
              method_name,
              ind_result[[i]]$ari,
              ind_result[[i]]$nmi,
              ind_result[[i]]$silhouette))
}
