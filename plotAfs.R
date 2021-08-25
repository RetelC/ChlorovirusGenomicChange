plotAfs <- function(
  af_mat, tps, add_v=FALSE, col_v=NA, alpha_v=1, tps_hline=numeric(), 
  xlim_v=numeric(), ylim_v=c(0, 1), lines_lwd=1, dots_cex=1.2, bg_col=1, 
  title_v="Derived allele frequencies", xlab_v="Time", 
  ylab_v="Derived allele frequency", 
  vline_sampletimes=TRUE, include_dots=TRUE, interpolate_lines=TRUE, 
  pdf_name=character(), seed_v=42
){
  ## Convenient wrapper functions to plot allele frequency values using little
  ## code. There are a lot of arguments to control what the figure looks like, 
  ## but except for af_mat and tps, they're all optional and have sensible 
  ## default values. 
  
  ## af_mat = matrix with allele frequency values. Rows correspond to 
  ##   loci, columns to time points (samples)
  ## tps = numeric vector of time points (same length as ncol(af_mat))
  ## add_v = logical of length 1. If TRUE, lines are added to existing 
  ##   plot using lines(), instead of calling a new plot. 
  ## col_v = colour vector of length 1 or of same length as nrow(af_mat). 
  ##   If NA, every locus gets a randomly assigned colour. 
  ## lines_lwd = line width
  ## dots_cex = dots size character expansion
  ## bg_col = colour of axis labels
  ## tps_hline = vector of numeric values where vertical lines 
  ##   should be drawn. If unspecified and vline_sampletimes=TRUE, 
  ##   a vertical line is drawn at every time point specified in {tps}
  ## xlim_v, ylim_v = numeric of length 2 specifying x resp. y limits
  ## title_v = character of length 1 specifying title
  ## xlab_v = character of length 1 specifying x axis label
  ## vline_sampletimes = logical. If TRUE, vertical lines are drawn
  ##   at time points specified in {tps}
  ## include_dots = logical of length 1. If TRUE, draw dots using 
  ##   points() function at observed allele frequency values
  ## interpolate_lines = logical of length 1. If TRUE, draw lines
  ##   between all pairs of adjacent non-NA values, instead of only 
  ##   between neighbouring ones. 
  ## pdf_name = optional name of file where the plot should be stored
  ## seed_v = seed value; only needed when col_v is set to NA
  ## * added ylab_v argument *
  ## * now able to handle one locus instead of a matrix of multiple rows *
  ## * now able to set alpha value (useful for overlapping lines) *
  ## * now able to call lines() to add to an existing plot
  ## instead of creating a new one. *
  
  ## check input
  if(!is.numeric(tps)){
    stop("plotAfs(): Arguments tps is not numeric")
  }
  if((length(add_v) != 1) | !is.logical(add_v)){
    stop("plotAfs(): Something is wrong with the add_v argument")
  }
  if((length(alpha_v) != 1) | !is.numeric(alpha_v) | 
     !(alpha_v >= 0 & alpha_v <= 1)){
    stop("plotAfs(): Something is wrong with the alpha_v argument")
  }
  if(!(is.numeric(lines_lwd) & is.numeric(xlim_v) & is.numeric(ylim_v) & 
       is.numeric(tps_hline) & (length(ylim_v) == 2))){
    stop("plotAfs(): One of the arguments lwd, xlim, ylim, tps_hline is wrong")
  }
  if( !((length(lines_lwd) == 1) & is.numeric(lines_lwd) & 
        (length(dots_cex)==1) & is.numeric(dots_cex)) ){
    stop("plotAfs(): One of the arguments lines_lwd or dots_cex is wrong")
  }
  if(!(is.logical(include_dots) & is.logical(interpolate_lines) & 
       is.logical(vline_sampletimes))){
    stop("plotAfs(): One of the arguments include_dots, interpolate_lines or vline_sampletimes is wrong")
  }
  if(!(is.character(title_v) & is.character(xlab_v) & is.character(ylab_v))){
    stop("plotAfs(): One of the arguments title, xlab or ylab is wrong")
  }
  if(!is.numeric(seed_v)){
    stop("plotAfs(): Invalid seed value")
  }
  
  if(!(is.matrix(af_mat) | is.numeric(af_mat))){
    stop("plotAfs(): Something wrong with af_mat argument")
  }
  if(is.matrix(af_mat)){
    if(length(tps) != ncol(af_mat)){
      stop("plotAfs(): Length of tps doesn't match dimension of af_mat")
    }
  }else if(length(tps) != length(af_mat)){
    stop("plotAfs(): Length of tps doesn't match length of af_mat")
  }
  
  ## find number of alleles and samples
  na <- ifelse(is.matrix(af_mat), nrow(af_mat), 1)
  ns <- ifelse(is.matrix(af_mat), ncol(af_mat), length(af_mat))
  
  ## if time points are not in order, re-order the dataset
  ## and the time points vector. 
  if(!all(order(tps) == 1:ns)){
    print("{tps} is not ordered by increasing value; reordered to produce sensible lines")
    if(is.matrix(af_mat)){
      af_mat <- af_mat[, order(tps)]
    } else{
      af_mat <- af_mat[order(tps)]
    }
    tps <- tps[order(tps)]
  }
    
  ## if not specified, set xlim to range of tps
  if(length(xlim_v) == 0){
    xlim_v <- c(min(tps, na.rm=TRUE), max(tps, na.rm=TRUE))
  }
  
  ## if no colours are set, randomly assign them
  if(is.na(col_v)[1]){
    print("No input colour vector; set with rainbow() and random assignment")
    set.seed(seed_v); col_v <- rainbow(max(na, 6))[sample(1:max(na, 6))[1:na]]
  }else if(length(col_v) == 1){
    col_v <- rep(col_v, na)
  }else if(length(col_v) < na){
    stop("Colour vector is shorter than the number of frequencies")
  }
  
  ## if alpha not equal to 1, change col_v accordingly
  if(alpha_v != 1){
    col_v <- apply(
      col2rgb(col_v) / 255, 2, 
      (function(x) rgb(x[1], x[2], x[3], alpha=alpha_v))
    )
  }
  
  ## open graphics device
  if(length(pdf_name) != 0){
    pdf(pdf_name, width=9, heigh=4)
  }
  
  ## create empty plot
  if(!add_v){
    par(mar=c(2, 4, 4, 2))
    plot(xlim_v, ylim_v, type='n', ylab="", xlab="")
    ## add labels
    mtext(side=1, text=xlab_v, 2.2, col=bg_col)
    mtext(side=2, text=ylab_v, 2.4, col=bg_col)
    mtext(side=3, text=title_v, col=bg_col)
  }
  ## add horizontal lines at sampling times
  if(vline_sampletimes & length(tps_hline)!=0){ 
    abline(v=tps_hline, col=2, lty=2) 
  }else if(vline_sampletimes & length(tps_hline)==0){
    abline(v=tps, col=2, lty=2)
  }
  
  ## plot allele frequencies
  if(is.matrix(af_mat)){
    for(i in 1:na){
      if(include_dots){ 
        points(tps, af_mat[i, ], col=col_v[i], pch=16, cex=dots_cex)
      }
      if(!interpolate_lines){
        lines(tps, af_mat[i, ], col=col_v[i], lwd=lines_lwd)
      }else{
        lines(tps[!is.na(af_mat[i, ])], af_mat[i, !is.na(af_mat[i, ])], 
              col=col_v[i], lwd=lines_lwd)
      }
    }
  }else{
    if(include_dots){
      points(tps, af_mat, col=col_v, pch=16, cex=dots_cex)
    }
    if(!interpolate_lines){
      lines(tps, af_mat, col=col_v, lwd=lines_lwd)
    }else{
      lines(tps[!is.na(af_mat)], af_mat[!is.na(af_mat)], 
            col=col_v, lwd=lines_lwd)
    }
  }
  
  ## close graphics device
  if(length(pdf_name) != 0){
    dev.off()
  }
}



plotAfs_pres <- function(
  af_mat, tps, add_v=FALSE, col_v=NA, alpha_v=1, tps_hline=numeric(), 
  xlim_v=numeric(), ylim_v=c(0, 1), lines_lwd=1, dots_cex=1.2, bg_col=1, 
  title_v="Derived allele frequencies", xlab_v="Time", 
  ylab_v="Derived allele frequency", 
  vline_sampletimes=TRUE, include_dots=TRUE, interpolate_lines=TRUE, 
  pdf_name=character(), seed_v=42
){
  ## same as above, but with little whitespace around edges. Used for 
  ## presentations. 
  
  ## check input
  if(!is.numeric(tps)){
    stop("plotAfs(): Arguments tps is not numeric")
  }
  if((length(add_v) != 1) | !is.logical(add_v)){
    stop("plotAfs(): Something is wrong with the add_v argument")
  }
  if((length(alpha_v) != 1) | !is.numeric(alpha_v) | 
     !(alpha_v >= 0 & alpha_v <= 1)){
    stop("plotAfs(): Something is wrong with the alpha_v argument")
  }
  if(!(is.numeric(lines_lwd) & is.numeric(xlim_v) & is.numeric(ylim_v) & 
       is.numeric(tps_hline) & (length(ylim_v) == 2))){
    stop("plotAfs(): One of the arguments lwd, xlim, ylim, tps_hline is wrong")
  }
  if( !((length(lines_lwd) == 1) & is.numeric(lines_lwd) & 
        (length(dots_cex)==1) & is.numeric(dots_cex)) ){
    stop("plotAfs(): One of the arguments lines_lwd or dots_cex is wrong")
  }
  if(!(is.logical(include_dots) & is.logical(interpolate_lines) & 
       is.logical(vline_sampletimes))){
    stop("plotAfs(): One of the arguments include_dots, interpolate_lines or vline_sampletimes is wrong")
  }
  if(!(is.character(title_v) & is.character(xlab_v) & is.character(ylab_v))){
    stop("plotAfs(): One of the arguments title, xlab or ylab is wrong")
  }
  if(!is.numeric(seed_v)){
    stop("plotAfs(): Invalid seed value")
  }
  
  if(!(is.matrix(af_mat) | is.numeric(af_mat))){
    stop("plotAfs(): Something wrong with af_mat argument")
  }
  if(is.matrix(af_mat)){
    if(length(tps) != ncol(af_mat)){
      stop("plotAfs(): Length of tps doesn't match dimension of af_mat")
    }
  }else if(length(tps) != length(af_mat)){
    stop("plotAfs(): Length of tps doesn't match length of af_mat")
  }
  
  ## find number of alleles and samples
  na <- ifelse(is.matrix(af_mat), nrow(af_mat), 1)
  ns <- ifelse(is.matrix(af_mat), ncol(af_mat), length(af_mat))
  
  ## if time points are not in order, re-order the dataset
  ## and the time points vector. 
  if(!all(order(tps) == 1:ns)){
    print("{tps} is not ordered by increasing value; reordered to produce sensible lines")
    if(is.matrix(af_mat)){
      af_mat <- af_mat[, order(tps)]
    } else{
      af_mat <- af_mat[order(tps)]
    }
    tps <- tps[order(tps)]
  }
  
  ## if not specified, set xlim to range of tps
  if(length(xlim_v) == 0){
    xlim_v <- c(min(tps, na.rm=TRUE), max(tps, na.rm=TRUE))
  }
  
  ## if no colours are set, randomly assign them
  if(is.na(col_v)[1]){
    print("No input colour vector; set with rainbow() and random assignment")
    set.seed(seed_v); col_v <- rainbow(max(na, 6))[sample(1:max(na, 6))[1:na]]
  }else if(length(col_v) == 1){
    col_v <- rep(col_v, na)
  }else if(length(col_v) < na){
    stop("Colour vector is shorter than the number of frequencies")
  }
  
  ## if alpha not equal to 1, change col_v accordingly
  if(alpha_v != 1){
    col_v <- apply(
      col2rgb(col_v) / 255, 2, 
      (function(x) rgb(x[1], x[2], x[3], alpha=alpha_v))
    )
  }
  
  ## open graphics device
  if(length(pdf_name) != 0){
    pdf(pdf_name, width=9, heigh=4)
  }
  
  ## create empty plot
  if(!add_v){
    par(mar=c(3, 4, 1, 1))
    plot(xlim_v, ylim_v, xaxt='n', yaxt='n', type='n',
         ylab="", xlab="")
    axis(1, at=20*(0:5), cex.axis=2, las=1, padj=.5)
    axis(2, at=.2*(0:5), cex.axis=2, las=1)
    ## add labels
    mtext(side=1, text=xlab_v, 2.2, col=bg_col)
    mtext(side=2, text=ylab_v, 2.4, col=bg_col)
    mtext(side=3, text=title_v, col=bg_col)
  }
  ## add horizontal lines at sampling times
  if(vline_sampletimes & length(tps_hline)!=0){ 
    abline(v=tps_hline, col=2, lty=2) 
  }else if(vline_sampletimes & length(tps_hline)==0){
    abline(v=tps, col=2, lty=2)
  }
  
  ## plot allele frequencies
  if(is.matrix(af_mat)){
    for(i in 1:na){
      if(include_dots){ 
        points(tps, af_mat[i, ], col=col_v[i], pch=16, cex=dots_cex)
      }
      if(!interpolate_lines){
        lines(tps, af_mat[i, ], col=col_v[i], lwd=lines_lwd)
      }else{
        lines(tps[!is.na(af_mat[i, ])], af_mat[i, !is.na(af_mat[i, ])], 
              col=col_v[i], lwd=lines_lwd)
      }
    }
  }else{
    if(include_dots){
      points(tps, af_mat, col=col_v, pch=16, cex=dots_cex)
    }
    if(!interpolate_lines){
      lines(tps, af_mat, col=col_v, lwd=lines_lwd)
    }else{
      lines(tps[!is.na(af_mat)], af_mat[!is.na(af_mat)], 
            col=col_v, lwd=lines_lwd)
    }
  }
  
  ## close graphics device
  if(length(pdf_name) != 0){
    dev.off()
  }
}


plotAfs_pub <- function(
  af_mat, tps, add_v=FALSE, col_v=NA, alpha_v=1, 
  xlim_v=numeric(), ylim_v=c(0, 1), lines_lwd=1, dots_cex=1.2, bg_col=1, 
  xaxis_ticks = 10*(0:10), yaxis_ticks = .2*(0:5), 
  include_dots=TRUE, interpolate_lines=TRUE, 
  pdf_name=character(), seed_v=42
){
  ## same as above, but with very slim margins, without titles and 
  ## axis labels, with user-specified axis ticks. 
  ## Used for publication (I add the text labels in Keynote)
  
  ## check input
  if(!is.numeric(tps)){
    stop("plotAfs(): Arguments tps is not numeric")
  }
  if((length(add_v) != 1) | !is.logical(add_v)){
    stop("plotAfs(): Something is wrong with the add_v argument")
  }
  if((length(alpha_v) != 1) | !is.numeric(alpha_v) | 
     !(alpha_v >= 0 & alpha_v <= 1)){
    stop("plotAfs(): Something is wrong with the alpha_v argument")
  }
  if(!(is.numeric(lines_lwd) & is.numeric(xlim_v) & is.numeric(ylim_v) & 
       is.numeric(xaxis_ticks) & is.numeric(yaxis_ticks) &
       (length(xlim_v) == 2) & (length(ylim_v) == 2))){
    stop("plotAfs(): One of the arguments lwd, xlim, ylim, xaxis_ticks or yaxis_ticks is wrong")
  }
  if( !((length(lines_lwd) == 1) & is.numeric(lines_lwd) & 
        (length(dots_cex)==1) & is.numeric(dots_cex)) ){
    stop("plotAfs(): One of the arguments lines_lwd or dots_cex is wrong")
  }
  if(!(is.logical(include_dots) & is.logical(interpolate_lines))){
    stop("plotAfs(): One of the arguments include_dots, or interpolate_lines is wrong")
  }
  if(!is.numeric(seed_v)){
    stop("plotAfs(): Invalid seed value")
  }
  
  if(!(is.matrix(af_mat) | is.numeric(af_mat))){
    stop("plotAfs(): Something wrong with af_mat argument")
  }
  if(is.matrix(af_mat)){
    if(length(tps) != ncol(af_mat)){
      stop("plotAfs(): Length of tps doesn't match dimension of af_mat")
    }
  }else if(length(tps) != length(af_mat)){
    stop("plotAfs(): Length of tps doesn't match length of af_mat")
  }
  
  ## find number of alleles and samples
  na <- ifelse(is.matrix(af_mat), nrow(af_mat), 1)
  ns <- ifelse(is.matrix(af_mat), ncol(af_mat), length(af_mat))
  
  ## if time points are not in order, re-order the dataset
  ## and the time points vector. 
  if(!all(order(tps) == 1:ns)){
    print("{tps} is not ordered by increasing value; reordered to produce sensible lines")
    if(is.matrix(af_mat)){
      af_mat <- af_mat[, order(tps)]
    } else{
      af_mat <- af_mat[order(tps)]
    }
    tps <- tps[order(tps)]
  }
  
  ## if not specified, set xlim to range of tps
  if(length(xlim_v) == 0){
    xlim_v <- c(min(tps, na.rm=TRUE), max(tps, na.rm=TRUE))
  }
  
  ## if no colours are set, randomly assign them
  if(is.na(col_v)[1]){
    print("No input colour vector; set with rainbow() and random assignment")
    set.seed(seed_v); col_v <- rainbow(max(na, 6))[sample(1:max(na, 6))[1:na]]
  }else if(length(col_v) == 1){
    col_v <- rep(col_v, na)
  }else if(length(col_v) < na){
    stop("Colour vector is shorter than the number of frequencies")
  }
  
  ## if alpha not equal to 1, change col_v accordingly
  if(alpha_v != 1){
    col_v <- apply(
      col2rgb(col_v) / 255, 2, 
      (function(x) rgb(x[1], x[2], x[3], alpha=alpha_v))
    )
  }
  
  ## open graphics device
  if(length(pdf_name) != 0){
    pdf(pdf_name, width=9, heigh=4)
  }
  
  ## create empty plot
  if(!add_v){
    par(mai=c(0.12, 0.12, 0.12, 0.12))
    plot(xlim_v, ylim_v, xaxt='n', yaxt='n', type='n', ylab="", xlab="")
    ## add ticks
    axis(side=1, at=xaxis_ticks, cex.axis=1.66, labels=F); 
    axis(side=2, at=yaxis_ticks, las=1, cex.axis=1.33, labels=F)
    
  }
  
  ## plot allele frequencies
  if(is.matrix(af_mat)){
    for(i in 1:na){
      if(include_dots){ 
        points(tps, af_mat[i, ], col=col_v[i], pch=16, cex=dots_cex)
      }
      if(!interpolate_lines){
        lines(tps, af_mat[i, ], col=col_v[i], lwd=lines_lwd)
      }else{
        lines(tps[!is.na(af_mat[i, ])], af_mat[i, !is.na(af_mat[i, ])], 
              col=col_v[i], lwd=lines_lwd)
      }
    }
  }else{
    if(include_dots){
      points(tps, af_mat, col=col_v, pch=16, cex=dots_cex)
    }
    if(!interpolate_lines){
      lines(tps, af_mat, col=col_v, lwd=lines_lwd)
    }else{
      lines(tps[!is.na(af_mat)], af_mat[!is.na(af_mat)], 
            col=col_v, lwd=lines_lwd)
    }
  }
  
  ## close graphics device
  if(length(pdf_name) != 0){
    dev.off()
  }
}


