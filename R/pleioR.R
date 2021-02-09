#' @title Fit a multitrait model to test for genetic pleiotropy
#' @description The main function of the pleiotest package. Fits a seemingly unrelated regression with possible unbalanced data, and covariates. It returns a pleio_class object to use as argument of the pleio_test function.
#' @param pheno A dataframe with phenotypic data with columns 'id', 'trait', and 'y'. The column 'y' contains the observations for the corresponding trait and id.
#' @param geno a matrix like object with the variants in columns and IDs in row names.
#' @param i vector with integers indicating rows from geno to use in the model.
#' @param j vector with integers indicating columns from geno to use in the model.
#' @param covariates (optional) a data frame or matrix containing covariates in columns and with IDs matching those in geno.
#' @param drop_subsets minimum sub-set sample size to consider for analysis.
#' @author Fernando M. Aguate (fmaguate@gmail.com)
#'
pleioR <- function(pheno, geno, i = NULL, j = NULL, covariates = NULL, drop_subsets = 10){

  if (is.null(rownames(geno)))
    stop('geno must have row names with ids')
  if (!class(pheno) %in% c('list', 'data.frame'))
    stop('pheno must be a list or a data frame')
  if (class(pheno) == 'data.frame'){
    pheno <- as.list(pheno)
  } else {
    if (var(sapply(pheno, length)) != 0)
      stop('If pheno is a list: id, trait and y should have the same length')
  }
  names(pheno)[grep('id|Var1', names(pheno), ignore.case = T)] <- 'id'
  names(pheno)[grep('trait|variable|Var2', names(pheno), ignore.case = T)] <- 'trait'
  names(pheno)[grep('y|value', names(pheno), ignore.case = T)] <- 'y'
  if (is.null(names(pheno)) | !all(names(pheno) %in% c('id', 'trait', 'y')))
    stop('pheno must be named with ids, variables and values. Try reshape2::melt().')
  if (!is.null(covariates)) {
    if (class(covariates) != 'matrix')
      covariates <- as.matrix(covariates)
    if (is.null(rownames(covariates)))
      stop('covariates must have row names with ids matching ids in geno')
  }

  if (is.null(i))
    i <- seq_len(nrow(geno))
  if (is.null(j))
    j <- seq_len(ncol(geno))

  y <- as.vector(pheno$y)
  y_na <- is.na(y)
  y <- y[!y_na]
  trait <- as.vector(pheno$trait)[!y_na]
  id <- as.vector(pheno$id)[!y_na]

  if(!any(rownames(geno)[i] %in% id))
    stop('At least some ids in pheno should match row names of geno')
  if(any(is.na(id)))
    stop('It cannot be NAs in ids. Remove them and try again.')
  subsets_list <- identify_subsets(trait, id)
  id_matrix <- subsets_list[[1]]
  subsets <- subsets_list[[2]]

  if(any(id_matrix[,ncol(id_matrix)] > drop_subsets)){
    subsets <- subsets[id_matrix[, ncol(id_matrix)] > drop_subsets]
    id_matrix <- id_matrix[id_matrix[, ncol(id_matrix)] > drop_subsets, , drop = F]
  }

  if(any(colSums(id_matrix[, -ncol(id_matrix), drop = F]) == 0))
    stop('Not enough data after dropping. Set a lower value for drop_subset')

  sets_xrows <- lapply(subsets, function(id_i)
    which(rownames(geno)[i] %in% id_i))
  if (!is.null(covariates))
    sets_cov <- lapply(subsets, function(id_i)
      which(rownames(covariates) %in% id_i))

  sets_rs <- list()
  sets_y <- list()
  for (set_i in 1:length(subsets)){
    if (length(subsets[[set_i]]) > 1){
      y_i <- y[id %in% subsets[[set_i]]]
      trait_i <- trait[id %in% subsets[[set_i]]]
      sets_y[[set_i]] <- do.call(cbind, split(y_i, trait_i))
      if (!is.null(covariates)){
        cov_i <- covariates[sets_cov[[set_i]], , drop = F]
        cov_i <- apply(cov_i, 2, function(x){ x[is.na(x)] <- mean(x, na.rm = T) ; x})
        for (set_j in 1:ncol(sets_y[[set_i]])){
          y_tmp <- sets_y[[set_i]][, set_j]
          sets_y[[set_i]][, set_j] <- y_tmp - predict(lm(y_tmp ~ cov_i), newdata = as.data.frame(cov_i))
        }
      }
      sets_rs[[set_i]] <- solve(cov(sets_y[[set_i]], use = 'complete.obs'))
    } else {
      sets_rs[[set_i]] <- 0
      sets_y[[set_i]] <- 0
    }
  }

  sets_which_trait <- lapply(sets_y, function(col_i) match(colnames(col_i), colnames(id_matrix)))

  pleio_res <- list()

  for (j_ in j){
    xj_ <- geno[i, j_]
    allele_freq <- mean(xj_, na.rm = T) / 2

    if(allele_freq != 0){
      xx_xy_list <- getXX_XYc(xj_, sets_xrows, sets_y, id_matrix, sets_which_trait)
      xrx_list <- xrsx_xrsy(id_matrix, sets_rs, xx_xy_list$xx, xx_xy_list$xy)
      lhs_solved <- solve(xrx_list$xrsx)
      betas <- lhs_solved %*% xrx_list$xrsy
      pleio_res[[length(pleio_res) + 1]] <- list('lhss' = lhs_solved,
                                                 'rhs' = xrx_list$xrsy,
                                                 'betas' = betas,
                                                 'n' = sum(!is.na(xj_)),
                                                 'allele_freq' = allele_freq)
      names(pleio_res)[length(pleio_res)] <- colnames(geno)[j_]
    }
  }

  class(pleio_res) <- append(class(pleio_res), 'pleio_class')
  return(pleio_res)
}
