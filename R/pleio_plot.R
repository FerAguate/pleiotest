#' @title Multi-trait manhattan plot
#' @description Plots the p-values that test the hypothesis of pleiotropic effects on n_traits.
#' @param pleio_res Object of class list that results of function pleio_test.
#' @param alpha either 'bonferroni05' or a numeric threshold for significance level.
#' @param n_traits number indicating the level of pleiotropy to plot.
#' @param bp_positions dataframe with chromosomes and basepair positions for SNPs matching pleio_res. rownames should contain SNP names.
#' @param set_colors string with 3 colors to use in the plot (by default: c('goldenrod4', 'brown4', 'royalblue2')).
#' @param set_text dataframe or matrix with text to add for some p-values. Row names should contain SNP names. Only the first column (with strings) will be considered.
#' @param set_plot boolean (TRUE by default) indicating whether to return the manhattan plot.
#' @param chr_spacing integer indicating the spacing (in basepair positions) between chromosomes. By default = 1e5.
#' @param ... Additional graphic parameters for the plot.
#' @author Original code by Fernando Aguate.
#'
pleio_plot <- function(pleio_res, alpha = 'bonferroni05', n_traits = 2, bp_positions = NULL, set_colors = NULL, set_text = NULL, set_plot = TRUE, chr_spacing = 1e5, ...){
  p_values <- apply(pleio_res[[1]][, 1:n_traits, drop = F], 1, max)
  p_values[p_values == 0] <- 1e-300
  p_notna <- !is.na(p_values)
  p_values <- p_values[p_notna]

  if (!n_traits > ncol(pleio_res[[2]])){
    indices <- as.character(pleio_res[[2]][, n_traits])
  } else {
    indices <- rep('', length(p_values))
  }
  indices <- indices[p_notna]

  if (alpha == 'bonferroni05')
    alpha <- 0.05 / length(p_values)

  if (is.null(set_colors)){
    set_colors <- c('#1F968BFF', '#39568CFF', '#440154FF')
  } else {
    if (length(set_colors) < 3)
      set_colors <- set_colors[1:3]
  }

  if (!is.null(bp_positions)){
    if (class(bp_positions) != 'data.frame'){
      warning ('bp_positions must be a data frame')
      bp_positions <- as.data.frame(bp_positions)
    }
    if (is.null(rownames(bp_positions)))
      stop ('bp_positions must have names matching names with pleio_res')
    if (any(is.na(bp_positions)))
      stop('bp_positions cannot have NAs')

    chr_col <- grep('chr', colnames(bp_positions), ignore.case = T, value = T)[1]
    bpp_col <- grep('pos', colnames(bp_positions), ignore.case = T, value = T)[1]
    if(length(chr_col) == 0 | length(bpp_col) == 0)
      stop('bp_positions must have col names for chromosomes (chr) and positions (pos)')

    snp_names <- rownames(bp_positions)
    chr_integers <- as.integer(bp_positions[, chr_col])

    if (length(bpp_col) > 0){
      pos_integers <- as.integer(bp_positions[, bpp_col])
      chr_pos <- data.frame(chr_integers, pos_integers - min(pos_integers), row.names = snp_names)
    } else {
      chr_pos <- data.frame(chr_integers, pos_integers = 1:length(chr_integers), row.names = snp_names)
    }

    chr_pos <- chr_pos[do.call(order, chr_pos), , drop = F]
    order_values <- match(rownames(chr_pos), names(p_values))
    chr_pos <- chr_pos[!is.na(order_values),]

    if (all(is.na(order_values)))
      stop ('bp_positions must have names matching names with pleio_res')

    my_colors <- rep(set_colors[3], length(p_values))
    my_colors[chr_pos[,1] %% 2 == 1] <- set_colors[2]

    if(length(unique(chr_pos[,1])) > 1){
      x_lab <- 'Chromosome'
    } else {
      x_lab <- 'Position'
    }

  } else {
    order_values <- 1:length(p_values)
    my_colors <- rep(set_colors[3], length(p_values))
    x_lab <- 'Position'
  }

  order_values <- order_values[!is.na(order_values)]
  p_values <- p_values[order_values]
  indices <- indices[order_values]

  p_significant <- p_values < alpha
  my_colors[p_significant] <- set_colors[1]

  if (!is.null(bp_positions)){
    if (is.unsorted(chr_pos[, 2])){
      for (i in unique(chr_pos[,1])){
        chr_pos[chr_pos[,1] == i, 2] <- chr_pos[chr_pos[,1] == i, 2] - min(chr_pos[chr_pos[,1] == i, 2]) + 1
        if (i != unique(chr_pos[,1])[1])
          chr_pos[chr_pos[,1] == i, 2] <- chr_pos[chr_pos[,1] == i, 2] + tmp_max + chr_spacing
        tmp_max <- max(chr_pos[chr_pos[,1] == i, 2])
      }
      pos_max <- chr_pos[nrow(chr_pos), 2]
    } else {
      pos_max <- length(p_values)
    }
  } else {
    pos_max <- length(p_values)
  }

  if (set_plot){
    plot(NULL, cex = .5,
         xlim = c(0, pos_max), ylim = c(0, ceiling(-log10(min(p_values, na.rm = T)))),
         ylab = paste0('-log10(p value)'), xlab = x_lab,
         main = paste0('Testing for association with ', n_traits, ' traits'), xaxt = 'n', ...)

    if (!is.null(bp_positions)){
      if(length(unique(chr_pos[,1])) > 1){
        x_axis <- aggregate(chr_pos[,2] ~ chr_pos[,1], FUN = function(x) round(mean(x)))
        axis(1, x_axis[,2], x_axis[,1])
        end <- cumsum(rle(chr_pos[,1])$lengths)
        start <- c(1, end[-length(end)] + 1)
        for (i in which(1:length(end) %% 2 == 1))
          rect(chr_pos[start[i], 2], -15, chr_pos[end[i], 2], 325, col = 'grey88', border = NA)
        points(chr_pos[, 2], -log10(p_values), col = my_colors, cex = .5)
      } else {
        axis(1, at = which(p_significant), labels = which(p_significant))
        points(1:length(p_values), -log10(p_values), col = my_colors, cex = .5)
      }
    } else {
      axis(1, at = which(p_significant), labels = which(p_significant))
      points(1:length(p_values), -log10(p_values), col = my_colors, cex = .5)
    }

    abline(h = -log10(alpha), lty = 2)

    if (!is.null(set_text)){
      pos_Text <- match(rownames(set_text), names(p_values))
      text(pos_Text, -log10(p_values[pos_Text]), labels = set_text[,1])
    }
  }
  w_sig <- unname(which(p_significant))
  result <- data.frame('p_value' = p_values[w_sig],
                       'index' = indices[w_sig],
                       row.names = names(p_values)[w_sig])
  if (!is.null(bp_positions))
    result$chr <- chr_pos[w_sig, 1]
    result$bp <- chr_pos[w_sig, 2]
  return(result)
}
