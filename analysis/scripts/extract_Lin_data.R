library(readxl)
library(dplyr)

# 1. 设置路径
data_path <- "H:/SPHARM_analysis/analysis/data/raw_data/Lin_2024_scar_data"

# 3. 获取 .txt 文件列表
file_list <- list.files(path = data_path, pattern = "\\.txt$", full.names = TRUE)

process_rhino_file <- function(file_path) {
  lines <- readLines(file_path)
  specimen_name <- gsub("\\.txt$", "", basename(file_path))

  object_breaks <- grep("Rhino object:", lines)
  if (length(object_breaks) == 0) return(NULL)

  specimen_data <- data.frame()
  scar_counter <- 1
  
  for (i in 1:length(object_breaks)) {
    start_idx <- object_breaks[i]
    end_idx <- if (i < length(object_breaks)) object_breaks[i+1] - 1 else length(lines)
    block <- lines[start_idx:end_idx]
    
    coord_regex <- "\\((-?[0-9.e+-]+),\\s*(-?[0-9.e+-]+),\\s*(-?[0-9.e+-]+)\\)"
    
    start_line <- block[grep("CV\\[ 0\\]|start\\s*=", block)]
    end_line <- block[grep("CV\\[ 1\\]|end\\s*=", block)]
    
    if (length(start_line) > 0 && length(end_line) > 0) {
      s_match <- regmatches(start_line, regexec(coord_regex, start_line))[[1]]
      e_match <- regmatches(end_line, regexec(coord_regex, end_line))[[1]]
      
      if (length(s_match) >= 4 && length(e_match) >= 4) {
        temp_df <- data.frame(
          Specimen_ID = specimen_name,
          Scar_ID = scar_counter,
          X1 = as.numeric(s_match[2]), Y1 = as.numeric(s_match[3]), Z1 = as.numeric(s_match[4]),
          X2 = as.numeric(e_match[2]), Y2 = as.numeric(e_match[3]), Z2 = as.numeric(e_match[4])
        )
        specimen_data <- rbind(specimen_data, temp_df)
        scar_counter <- scar_counter + 1
      }
    }
  }
  
  if (nrow(specimen_data) > 0) {
    specimen_data$Length <- sqrt((specimen_data$X2 - specimen_data$X1)^2 + 
                                   (specimen_data$Y2 - specimen_data$Y1)^2 + 
                                   (specimen_data$Z2 - specimen_data$Z1)^2)
  }
  
  return(specimen_data)
}

# 5. 执行处理
all_results_list <- lapply(file_list, process_rhino_file)
final_combined_df <- do.call(rbind, all_results_list)

# 7. 格式转换为标准管道格式
format_for_pipeline <- function(df) {
  
  # 提取 Typology（文件名中下划线后的部分）
  # 例如：EXP01_Levallois preferential → Levallois preferential
  df$ID       <- df$Specimen_ID
  df$Typology <- sub("^[^_]+_", "", df$Specimen_ID)
  
  # 重命名坐标列
  df$Start_X <- df$X1
  df$Start_Y <- df$Y1
  df$Start_Z <- df$Z1
  df$End_X   <- df$X2
  df$End_Y   <- df$Y2
  df$End_Z   <- df$Z2
  
  # 计算 Direct 向量（End - Start 的单位向量）
  dv <- data.frame(
    dx = df$End_X - df$Start_X,
    dy = df$End_Y - df$Start_Y,
    dz = df$End_Z - df$Start_Z
  )
  norm <- sqrt(dv$dx^2 + dv$dy^2 + dv$dz^2)
  norm <- ifelse(norm < 1e-10, 1, norm)
  
  df$Direct_X <- dv$dx / norm
  df$Direct_Y <- dv$dy / norm
  df$Direct_Z <- dv$dz / norm
  
  # 用每个标本所有起点的质心作为参考点 Pos_X/Y/Z
  centroid <- df %>%
    group_by(ID) %>%
    summarise(
      Pos_X = mean(Start_X),
      Pos_Y = mean(Start_Y),
      Pos_Z = mean(Start_Z),
      .groups = "drop"
    )
  df <- df %>% left_join(centroid, by = "ID")
  
  # 只保留管道需要的列，顺序和现有数据一致
  df <- df %>%
    select(ID, Pos_X, Pos_Y, Pos_Z,
           Scar_ID,
           Start_X, Start_Y, Start_Z,
           End_X,   End_Y,   End_Z,
           Length,
           Direct_X, Direct_Y, Direct_Z,
           Typology)
  
  return(df)
}

# 执行转换
pipeline_df <- format_for_pipeline(final_combined_df)

# 检查
cat("标本数量：", n_distinct(pipeline_df$ID), "\n")
cat("片疤总数：", nrow(pipeline_df), "\n")
cat("类型分布：\n")
print(pipeline_df %>% distinct(ID, Typology) %>% count(Typology))

# 8. 导出原始版本
write.csv(final_combined_df, output_file, row.names = FALSE)

# 9. 导出管道版本
pipeline_file <- file.path(data_path, "Scar_Vectors_Lin2024_pipeline.csv")
write.csv(pipeline_df, pipeline_file, row.names = FALSE)
cat("\n管道格式已保存：", pipeline_file, "\n")

# 6. 导出
output_file <- file.path(data_path, "Scar_Vectors_Lin2024.csv")
write.csv(final_combined_df, output_file, row.names = FALSE)
