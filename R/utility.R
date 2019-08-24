# utility functions

#' Make sparse annotation matrix from annotation dataframe
#'
#' @param anno_df dataframe with columns 
#' @param row_total total number of SNPs/loci
#'
#' @return dgCMatrix with features as columns and loci as rows
#' @export 
#'
spread_sparse <- function(anno_df,SNP="SNP",feature="feature",row_total=length(unique(anno_df$SNP))){
  SNP_var <- as_string(ensym(SNP))
  feature_var <- as_string(ensym(feature))
  sv <- anno_df[[SNP_var]]
  fv <- anno_df[[feature_var]]
  stopifnot(!is.null(sv),
            !is.null(fv),
            is.integer(sv))
  
  fv <- as.factor(fv)
  Matrix::sparseMatrix(i = sv,
                       j = as.integer(fv),
                       x = rep(1.0,length(sv)),
                       dims = c(row_total,length(levels(fv))),
                       dimnames = list(NULL,levels(fv)))
}


