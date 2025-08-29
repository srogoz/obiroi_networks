###########
###########
#import interpolated (antrax_xy data)
library(arrow)
library(R.matlab)
library(BiocManager)
library(rhdf5)
library(dplyr)
#setwd("~/ulr-lab/Users/Sarah_saro7615/EXP1_GROUPSIZEandCOMPOSITION")
df <-read_parquet("processed/b16-1_interpo_20250731.parquet", as_data_frame = TRUE)
df$frame <-as.numeric(rownames(df))

interaction_threshold <- 0.002 #interaction registered if distance closer than bodylength
df<- df %>% select(-colony) #drop tag, because all ants in one colony
#######################
#######################
#plot heatmap of (w,y)coordinates of ants


get_network <- function(df, interaction_threshold = 0.002) {
  # Ensure frame is integer-like
  df$frame <- as.integer(df$frame)
  
  total_frames <- max(df$frame, na.rm = TRUE)
  present_ants <- unique(df$ant)
  number_ants <- length(present_ants)
  
  # Create a mapping from ant name to index
  ant_to_idx <- setNames(seq_along(present_ants), present_ants)
  
  # Initialize a 3D array of zeros
  interaction_matrix <- array(0, dim = c(number_ants, number_ants, total_frames))
  
  for (t in seq_len(total_frames)) {
    frame_data <- subset(df, frame == t)
    
    for (ant1 in present_ants) {
      for (ant2 in present_ants) {
        if (ant1 == ant2) next
        
        row1 <- frame_data[frame_data$ant == ant1, ]
        row2 <- frame_data[frame_data$ant == ant2, ]
        
        # Skip if any data is missing
        if (nrow(row1) == 0 || nrow(row2) == 0) next
        if (any(is.na(c(row1$x, row1$y, row2$x, row2$y)))) next
        
        # Calculate Euclidean distance
        dist <- sqrt((row1$x - row2$x)^2 + (row1$y - row2$y)^2)
        
        if (dist <= interaction_threshold) {
          idx1 <- ant_to_idx[[ant1]]
          idx2 <- ant_to_idx[[ant2]]
          interaction_matrix[idx1, idx2, t] <- 1
        }
      }
    }
  }
  
  return(interaction_matrix)
}

###########
#apply function

interaction_mat <- get_network(df, interaction_threshold)