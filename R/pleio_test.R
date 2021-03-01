#' @title Sequential Wald test for pleiotropy
#' @description Performs the sequential test to test pleiotropic effects from results of pleioR.
#' @param pleio_results object of class pleio_class (see function pleioR).
#' @param loop_breaker A numerical value for a maximum p-value used to stop the sequence if they pass this threshold.
#' @param save_at string with directory and/or file name (.rdata) to save the results.
#' @author Original code by Fernando Aguate.
#'
pleio_test <- function(pleio_results, loop_breaker = 0.99, save_at = NULL){
  if (!'pleio_class' %in% class(pleio_results))
    stop('pleio_results should be a pleio_class object')

  n_traits <- length(pleio_results[[1]]$rhs)
  traits_names <- row.names(pleio_results[[1]]$betas)
  indices <- list()
  for (i in 1:(n_traits - 1))
    indices[[length(indices) + 1]] <- combn(n_traits, i)

  contrast_matrices <- c(list(list(diag(n_traits))), lapply(indices, function(x)
    lapply(as.data.frame(x), function(w)
      diag(n_traits)[-w, , drop = F])))

  contrast_matrices_indices <- lapply(contrast_matrices[-1], function(x)
    lapply(x, function(w)
      paste(which(colSums(w) == 0), collapse = '_')))

  tests_res <- PleioSeqTestc(pleio_results, contrast_matrices, contrast_matrices_indices, loop_breaker)

  p_values_res <- do.call(rbind, lapply(tests_res, function(x) t(x$pValues)))
  indices_res <- as.data.frame(do.call(rbind, lapply(tests_res, function(x) t(x$index))))

  colnames(p_values_res) <- paste0('p', 1:ncol(p_values_res))
  colnames(indices_res) <- paste0('ind', 1:ncol(indices_res))
  if (length(names(pleio_results)) > 0)
    rownames(p_values_res) <- rownames(indices_res) <- names(pleio_results)

  names(traits_names) <- 1:length(traits_names)

  tests_res <- list('pValues' = p_values_res, 'Index' = indices_res[, -ncol(indices_res)], 'traits' = traits_names)

  if (!is.null(save_at)){
    if(!is.character(save_at))
      stop('Save at must be a character string with a directory or file name')
    if(length(grep('.rdata', save_at, ignore.case = T)) > 0){
      file_name <- strsplit(casefold(basename(save_at)), '.rdata')[[1]]
      dir_name <- dirname(save_at)
    }else{
      if(length(grep('\\.', save_at)) > 0)
        stop('The file must be .rdata')
      dir_name <- save_at
      file_name <- 'pleio_test_result'
    }

    if(!substring(dir_name, nchar(dir_name)) == '/')
      dir_name <- paste0(dir_name, '/')

    if(!dir.exists(dir_name))
      dir.create(dir_name)

    i <- 1
    while (file.exists(paste0(dir_name, file_name, '_', i, '.rdata')))
      i <- i + 1
    save(tests_res, file = paste0(dir_name, file_name, '_', i, '.rdata'))
  }

  return(tests_res)
}

