
################################################# DATA ORGANIZATION #######################################

library(officer) 
library(magrittr) 

# Define the base directory containing the image 
base_results_dir <- "./results_after_sva_pca_all_tissues" 
analysis_dir <- "./analysis_sva_pca_tissues" 

# Create analysis directory if it doesn't exist
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)  

# Define the path for the final PowerPoint presentation 
output_pptx_all <- file.path(analysis_dir, "All_Tissues_Results_Presentation.pptx")  

ppt_all <- read_pptx()  

# Get the PNG files in the directory 
image_files <- list.files(path = base_results_dir, pattern = "\\.png$", full.names = TRUE) 

# Iterate through each image file and add to the PowerPoint
for (image in image_files) {
  print(paste("Processing image:", image))   
  
  ppt_all <- ppt_all %>%   
    add_slide(layout = "Title and Content", master = "Office Theme") %>%   
    ph_with(value = "PCA Plot", location = ph_location_type(type = "title")) %>%   
    ph_with(external_img(image, width = 6, height = 4.5), location = ph_location(left = 1, top = 1.5, width = 6.2, height = 6))   
}  

# Save the PowerPoint presentation
print(ppt_all, target = output_pptx_all)  
cat("PowerPoint presentation created successfully at:", output_pptx_all, "\n")
