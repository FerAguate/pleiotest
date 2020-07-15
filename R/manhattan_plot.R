#' @title Single Trait Manhattan plot
#' @description plots p-values with results from mt_gwas.
#' @param mt_gwas_results Results from function mt_gwas
#' @param trait integer number indicating the position of the trait to be plotted.
#' @param bp_positions dataframe with SNPs basepair positions. colnames should be 'chr' and 'position', rownames should be SNP identifiers matching names in mt_gwas.
#' @param ... Further graphical parameters to customize the plot. Options include: title, bty, pch, cex.lab and/or cex.main.
#' @author Original code by Fernando Aguate.
#'
manhattan_plot <- function(mt_gwas_results, trait, bp_positions, ...){
  colnames(bp_positions) <- c('chr', 'position')
  result_trait <- mt_gwas_results[[trait]]

  if(length(unlist(strsplit(rownames(result_trait)[1], '_'))) == 2){
    rownames(result_trait) <- sapply(strsplit(rownames(result_trait), '_'), function(x) x[1])
  }
  matched_bp <- match(rownames(result_trait), rownames(bp_positions))
  result_trait <- cbind(result_trait, bp_positions[matched_bp,])
  result_trait <- result_trait[order(result_trait$chr, result_trait$position),]
  my_my_colors <- as.factor(result_trait$chr)
  levels(my_colors) <- rep(c('blue1', 'darkorange3'), length.out = 26)
  plot(1:nrow(result_trait), -log10(result_trait$`p value`), col = as.character(my_colors),
       cex = .8, xaxt = 'n', xlab = 'chromosome', ylab = '-log10(p value)', ...)
  x_axis <- aggregate(1:nrow(result_trait) ~ result_trait$chr, FUN = function(x) round(mean(x)))
  axis(1, x_axis[,2], x_axis[,1])
}
