# Function to write tree scructure for folders

# Install fs package if not already installed
if (!requireNamespace("fs", quietly = TRUE)) {
      install.packages("fs")
}

# Load the fs package


# Define a function to get and save the folder tree with symbols
save_folder_tree_with_symbols <- function(path = ".", output_file = "folder_tree.txt") {
      library(fs)
      # Use fs::dir_tree to create the folder tree and get individual lines
      folder_tree <- dir_tree(path, recurse = TRUE)
      
      # Define a function to add symbols based on depth
      add_symbols <- function(lines) {
            result <- c()
            depth <- 0
            
            for (i in seq_along(lines)) {
                  # Calculate the depth based on the number of slashes in the path
                  depth <- sum(strsplit(lines[i], "/")[[1]] != "")
                  
                  # Determine if it's a file or folder and if it's the last in its group
                  is_last <- (i == length(lines)) || 
                        (sum(strsplit(lines[i + 1], "/")[[1]] != "") <= depth)
                  
                  # Add appropriate symbols for the depth level
                  prefix <- ifelse(depth > 1, paste0(strrep("│  ", depth - 2), ifelse(is_last, "└─ ", "├─ ")), "")
                  
                  result <- c(result, paste0(prefix, basename(lines[i])))
            }
            return(result)
      }
      
      # Create a visual representation of the folder tree with symbols
      visual_tree <- add_symbols(folder_tree)
      
      # Print the visual folder tree to the console
      cat(visual_tree, sep = "\n")
      
      # Save the visual folder tree to a text file
      writeLines(visual_tree, output_file)
}

# Call the function on the

save_folder_tree_with_symbols(path = here::here("SCRIPTS"), output_file = "folder_tree_SCRIPTS.txt")
save_folder_tree_with_symbols(path = here::here("ANALYSES"), output_file = "folder_tree_ANALYSES.txt")
save_folder_tree_with_symbols(path = here::here("DATASETS"), output_file = "folder_tree_DATASETS.txt")
