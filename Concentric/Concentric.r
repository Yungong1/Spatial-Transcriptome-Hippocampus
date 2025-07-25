library(Seurat)
library(dplyr)
library(ggplot2)

# Function to find spots within a certain circular radius
find_circular_neighbors <- function(data, vas_spot, radius_min, radius_max) {
  neighbors <- data %>%
    filter(
      sqrt((col - vas_spot$col)^2 + (row - vas_spot$row)^2) > radius_min &  # Outside inner radius
      sqrt((col - vas_spot$col)^2 + (row - vas_spot$row)^2) <= radius_max   # Inside outer radius
    )
  return(neighbors)
}

concentric = function(meta_data, distance1, distance2, distance3)
{
    # Identify VAS spots
    vas_spots <- meta_data %>%
      filter(major_annotation == "VAS") %>%
      select(col, row)

    # Step 1: Get first-layer neighbors (r = 0 to 1)
    neighbors_r1 <- lapply(1:nrow(vas_spots), function(i) {
      find_circular_neighbors(meta_data, vas_spots[i, ], radius_min = 0, radius_max = distance1)
    }) %>% bind_rows() %>% distinct(col, row, .keep_all = TRUE)

    # Step 2: Get second-layer neighbors (r = 1 to 2)
    neighbors_r2 <- lapply(1:nrow(vas_spots), function(i) {
      find_circular_neighbors(meta_data, vas_spots[i, ], radius_min = distance1 + 1, radius_max = distance2)
    }) %>% bind_rows() %>% distinct(col, row, .keep_all = TRUE)

    # Step 3: Get third-layer neighbors (r = 2 to 3)
    neighbors_r3 <- lapply(1:nrow(vas_spots), function(i) {
      find_circular_neighbors(meta_data, vas_spots[i, ], radius_min = distance2 + 1, radius_max = distance3)
    }) %>% bind_rows() %>% distinct(col, row, .keep_all = TRUE)

    # Find unique in layer 2
    unique_r2_barcode = setdiff(rownames(neighbors_r2), rownames(neighbors_r1)) 
    neighbors_r2 = neighbors_r2[unique_r2_barcode, ]

    # Find unique in layer 3
    unique_r3_barcode = setdiff(rownames(neighbors_r3), rownames(neighbors_r2))
    unique_r3_barcode = setdiff(unique_r3_barcode, rownames(neighbors_r1))
    neighbors_r3 = neighbors_r3[unique_r3_barcode, ]
    
    # Add the concentric information into the original meta data
    vas_spots <- meta_data %>%
      filter(major_annotation == "VAS")

    vas_spots$distance = "layer0"
    neighbors_r1$distance = "layer1"
    neighbors_r2$distance = "layer2"
    neighbors_r3$distance = "layer3"

    meta_data_dis = rbind(vas_spots, neighbors_r1, neighbors_r2, neighbors_r3)
    
    # Ensure both datasets have rownames as a column
    meta_data_dis$spot_id = rownames(meta_data_dis)

    meta_data$spot_id = rownames(meta_data)

    # Merge datasets using left_join
    merged_data <- meta_data %>%
      left_join(meta_data_dis %>% select(spot_id, distance), by = "spot_id") %>%
      mutate(distance = ifelse(is.na(distance), "Blank", distance))  # Assign "Blank" if no match
    
    return(merged_data)
}

meta_all = spot_data@meta.data
meta_all_dis = data.frame()
for (sample in unique(meta_all$ID))
{
    sub_data = subset(spot_data, subset = ID == sample)
    meta_data = sub_data@meta.data
    meta_data = concentric(meta_data, distance1 = 200, distance2 = 400, distance3 = 600)
    meta_all_dis = rbind(meta_all_dis, meta_data)
}

# Visulization
sub_data = subset(spot_data, subset = ID == "0618")
meta_data = sub_data@meta.data
meta_data = concentric(meta_data, distance1 =200, distance2 = 400, distance3 = 600)

# Define colors for each layer
color_mapping <- c("Blank" = "gray", 
                   "layer0" = "red", 
                   "layer1" = "pink", 
                   "layer2" = "blue", 
                   "layer3" = "purple")

# Plot the spatial distribution of spots with color mapping
pdf("Results/Region_plot/concentric_PART.pdf", width = 6, height = 6)
ggplot(meta_data, aes(x = col, y = row, color = distance)) +
  geom_point(size = 1.5, alpha = 0.8) +  # Plot spots with assigned colors
  scale_color_manual(values = color_mapping) +  # Use predefined colors
  labs(title = "Spatial Distribution of Spots",
       x = "Column",
       y = "Row",
       color = "Distance Layer") +
  coord_fixed() +  # Maintain spatial proportions
  theme_minimal()
dev.off()