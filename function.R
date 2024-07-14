### define function
# 
extract_id <- function(input_string) {
  # 使用sub函数将第一个下划线后面的内容替换为空字符串
  cleaned_string <- sub("_.*", "", input_string)
  return(cleaned_string)
}

#
get_taxonomy <- function(text) {
  parts <- unlist(strsplit(text, "\\|"))
  species_col <- parts[7]
  modified_parts <- gsub('s__', '', species_col)
  return(modified_parts)
}

#
cosine_similarity <- function(x, y) {
  return(sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2))))
}                                                        #consine similarity  (compare each paired sample)

#
spearman_correlation <- function(x, y) {
  return(cor(x, y, method = "spearman"))
}                                                        #spearman correlation  (compare each paired sample)


#
calculate_jaccard_similarity <- function(mat1, mat2) {
  if (is.null(names(mat1)) || is.null(names(mat2))) {
    stop("Both matrices must have row names.")
  }
  
  # extract rownames of non-zero row
  non_zero_rows_mat1 <- names(mat1[mat1 != 0])
  non_zero_rows_mat2 <- names(mat2[mat2 != 0])
  
  # calculate intersect and union set
  intersection <- length(intersect(non_zero_rows_mat1, non_zero_rows_mat2))
  union <- length(union(non_zero_rows_mat1, non_zero_rows_mat2))
  
  # culculate jaccard similarity
  jaccard_similarity <- intersection / union
  
  return(jaccard_similarity)
}                                                      # jaccard similarity  (compare each paired sample)


#

#euclidean_distance <- function(x, y) {
#  return(Rfast::Dist(x, y, method = "euclidean"))               
#} 

euclidean_distance <- function(x, y) {
  sqrt(sum((x - y)^2))                                  # Euclidean_distance (compare each paired sample)
}

# RMES (compare each paired paired sample/EC number)
rmse <- function(x, y) {
  return(RMSE(x, y))
}   #spearman correlation  (compare each paired sample)

get_group_name3 <- function(string) {
  
  result <- sub("^[^.]+\\.[^.]+\\.", "", string)
  return(result)
  
} # for genefamily output






