#
# Note: the code of these function was originally written by S.Mella
#

factor_resolution <- function(so){
  tmp_df <- so@meta.data
  ind <- grep("res.", colnames(tmp_df))
  for(i in ind){
    tmp_df[,i] <- factor(tmp_df[,i])
    tmp_df[,i] <- factor(tmp_df[,i],
                         levels = 0 : (length(levels(tmp_df[,i])) - 1) )
  }
  so@meta.data <- tmp_df
  return(so)
}
