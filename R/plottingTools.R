#' @include N-class.R
NULL

#' Extract ggplot legend.
#' 
#' Returns the legend from a ggplot plotting object.
#' 
#' @importFrom ggplot2 ggplot_gtable
#' @importFrom ggplot2 ggplot_build
#' 
#' @param gplot ggplot object from which to pull legend
#' @return legend ggplot object.
#' @export
g_legend = function(gplot){
  tmp <- ggplot_gtable(ggplot_build(gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


#' Get y-axis plotting margins
#'
#' Gets flexible y-axis margins for basic N, Shape, View, Strand, and Round plots.
#'
#' @param x Value
#' @param dig Number of digits to use.
yMarg = function(x, dig =2) {
  if (x >= .01) {
    yM = round(x, dig)
  } else {
    yM = signif(x, dig)
  }
  return (yM)
}


#' Get data.frame demarcating regions for a ggplot object.
#' 
#' Builds data.frame of regions used to shade different regions of the basic
#' plots for object of class 'N' and 'Shape'.
#' 
#' @importFrom RColorBrewer brewer.pal
#' 
#' @param labels labels for each section to be assigned a rectangle
#' @param lengths lengths for each of the new sections
getRegions = function(labels, lengths) {
  numReg = length(labels[[1]])
  xstart = rep(0, numReg)
  xend = rep(0, numReg)
  col = labels[[1]]
  
  xstart[1] = .5
  xend[numReg] = .5+sum(lengths[[1]])
  for (i in 1:(numReg-1)) {
    xstart[(i+1)] = .5+sum(lengths[[1]][1:i])
    xend[i] = .5+sum(lengths[[1]][1:i])
  }
  rects = data.frame(xstart = xstart,
                     xend = xend,
                     col = col)
  rects$col = factor(rects$col, levels = unique(rects$col))
  regionColors = c(RColorBrewer::brewer.pal(7, "Greys"))
  rects$hexCol = ""
  if (length(levels(rects$col))  >=  1) {
    rects$hexCol[rects$col == levels(rects$col)[1]] = regionColors[1]
  }
  if (length(levels(rects$col)) >= 2) {
    rects$hexCol[rects$col == levels(rects$col)[2]] = regionColors[4]
  }
  if (length(levels(rects$col)) >= 3) {
    rects$hexCol[rects$col == levels(rects$col)[3]] = regionColors[7]
  }
  if (length(levels(rects$col)) >= 4) {
    rects$hexCol[rects$col == levels(rects$col)[4]] = regionColors[3]
  }
  if (length(levels(rects$col)) >= 5) {
    rects$hexCol[rects$col == levels(rects$col)[5]] = regionColors[5]
  }
  if (length(levels(rects$col)) >= 6) {
    rects$hexCol[rects$col == levels(rects$col)[6]] = regionColors[2]
  }
  if (length(levels(rects$col)) >= 7) {
    rects$hexCol[rects$col == levels(rects$col)[7]] = regionColors[6]
  }
  return (rects)
}

#' Plot PSAM
#' 
#' Returns standardized PSAM plot
#' 
#' @param model model from which to take nucleotide feature values
#' @param title title for plot
#' @param maxErr maximum height of error bars to use in plotting
#' @param iter Iteration of model to use for mono nucleotide values
#' @param regs Optional parameter giving labels to use for region shading. If used, lengths parameter must also be specified.
#' @param lengths lengths for each of the new sections
#' @export
plotPSAM = function(model, title="", maxErr=1, iter=NULL, regs=NULL, lengths = NULL) {
  values.Plot = getPlotValues(model@features@N, maxErr, iter)
  lab.Plot = getTopSeqLabels(model@features@N)
  nc = ncol(getPSAM(model))
  if ((is.null(regs)) & (is.null(lengths))) {
    psamPlot = ggplot2::ggplot(values.Plot)+
      geom_point(data = values.Plot, aes(x=PosID, y=Affinity, colour=N), size = 1, position = position_dodge(width = .2)) +
      geom_errorbar(data = values.Plot, aes(x=PosID, ymax = PlotErrorMax, ymin=PlotErrorMin, colour = N), width = 1, position = position_dodge(width = .2))+
      scale_colour_manual(name = "", values=c("green","blue","#FFC000","red"), breaks=c("A", "C", "G", "T"),
                          labels=c("A           ", "C   ", "G   ", "T   "))+
      scale_x_discrete(limits = 1:nc,
                       labels=parse(text=lab.Plot))+
      ylim(-.01, 1.01)+
      theme_bw()+
      ggtitle(title)+
      xlab(NULL)+ylab(expression(paste("Relative Affinity")))
  } else {
    if (is.null(lengths)) {
      stop('If regs is specified, lengths for the regions must also be included.')
    } else if (is.null(regs)) {
      numReg = length(lengths[[1]])
      regs = list(as.character(1:numReg))
    }
    rects = getRegions(regs, lengths)
    psamPlot = ggplot2::ggplot(values.Plot)+
      ggplot2::geom_rect(data = rects, 
                aes(xmin = xstart,
                    xmax = xend,
                    ymin = -Inf,
                    ymax = Inf, 
                    fill = col), 
                alpha = 0.4)+
      ggplot2::geom_point(data = values.Plot, aes(x=PosID, y=Affinity, colour=N), size = 1, position = position_dodge(width = .2)) +
      ggplot2::geom_errorbar(data = values.Plot, aes(x=PosID, ymax = PlotErrorMax, ymin=PlotErrorMin, colour = N), width = 1, position = position_dodge(width = .2))+
      ggplot2::scale_colour_manual(name = "", values=c("green","blue","#FFC000","red"), breaks=c("A", "C", "G", "T"),
                          labels=c("A           ", "C   ", "G   ", "T   "))+
      ggplot2::scale_x_discrete(limits = 1:nc,
                       labels=parse(text=lab.Plot))+
      ggplot2::ylim(-.01, 1.01)+
      ggplot2::theme_bw()+
      ggplot2::ggtitle(title)+
      ggplot2::xlab(NULL)+ggplot2::ylab(expression(paste("Relative Affinity")))
  }
  
  
  return(psamPlot)
}


#' Get sequence labels for psam ratio plot.
#' 
#' Gets list of nucleotides maximizing ratio of PSAMs from 2 objects of class 'model' across the models' footprints.
#' 
#' @param model1 First model used to extract mono-nucleotide feature top ratio labels.
#' @param model2 Second model used to extract mono-nucleotide feature top ratio labels.
#' @return Vector of nucleotides, with nucleotide at index i being the nucleotide that maximizes the ratio of model1 nucleotide PSAM to model2 nucleotide PSAM at the i'th position.
getTopRatioLabels = function(model1, model2, iter = NULL) {
  PSAM1 = getPSAM(model1, iter)
  PSAM2 = getPSAM(model2, iter)
  PSAM.Ratio = PSAM1/PSAM2
  lab.plot = c()
  Ns = c("A", "C", "G", "T")
  for (i in 1:model1@fpLen) {
    lab.plot = c(lab.plot, Ns[(PSAM.Ratio[,i] == max(PSAM.Ratio[,i]))][1])
  }
  return (lab.plot)
}

#' Dot plot of ratios between nucleotide psams from 2 objects of class 'model'.
#' 
#' Returns standardized PSAM plot of ratios between two model objects' nucleotide affinities.
#' 
#' @param model1 Model1 from which to take first set of nucleotide feature values.
#' @param model2 Model2 from which to take first set of nucleotide feature values.
#' @param title Title for plot.
#' @param maxErr Maximum height of error bars to use in plotting.
#' @param iter Iteration of model to use for mono nucleotide values.
#' @param regs Optional parameter giving labels to use for region shading. If used, lengths parameter must also be specified.
#' @param lengths Length of region parameters to be included.
#' @export
plotPSAMComp.Ratios = function(model1, model2, title="", yTitle = "", maxErr=1, iter=NULL, regs=NULL, lengths = NULL) {
  values.Plot1 = getPlotValues(model1@features@N, maxErr, iter)
  values.Plot2 = getPlotValues(model2@features@N, maxErr, iter)
  values.CompPlot = merge(values.Plot1, values.Plot2, by = c("PosID", "N"), all = TRUE, sort = FALSE)
  values.CompPlot$Affinity.Ratio = values.CompPlot$Affinity.x/values.CompPlot$Affinity.y
  values.CompPlot$SE.Ratio = sqrt(values.CompPlot$SE.x^2+values.CompPlot$SE.y^2)
  values.CompPlot$PlotErrorMax = values.CompPlot$Affinity.Ratio*(1+values.CompPlot$SE.Ratio)
  values.CompPlot$PlotErrorMin = values.CompPlot$Affinity.Ratio*(1-values.CompPlot$SE.Ratio)
  values.CompPlot$PlotErrorMax[values.CompPlot$PlotErrorMax>1+maxErr] = 1+maxErr
  values.CompPlot$PlotErrorMin[values.CompPlot$PlotErrorMin<1-maxErr] = 1-maxErr
  lab.Plot = getTopRatioLabels(model1, model2, iter)
  nc = ncol(getPSAM(model1))
  minPlotPoint = min(values.CompPlot$PlotErrorMin)
  maxPlotPoint = max(values.CompPlot$PlotErrorMax)
  ymin = ((minPlotPoint %/% .25))*.25
  ymax = ((maxPlotPoint %/% .25)+1)*.25

  if ((is.null(regs)) & (is.null(lengths))) {
    psamRatioPlot = ggplot2::ggplot(values.CompPlot)+
      ggplot2::geom_point(data = values.CompPlot, aes(x=PosID, y=Affinity.Ratio, colour=N), size = 1, position = position_dodge(width = .2)) +
      ggplot2::geom_errorbar(data = values.CompPlot, aes(x=PosID, ymax = PlotErrorMax, ymin=PlotErrorMin, colour = N), width = 1, position = position_dodge(width = .2))+
      ggplot2::scale_colour_manual(name = "", values=c("green","blue","#FFC000","red"), breaks=c("A", "C", "G", "T"),
                          labels=c("A           ", "C   ", "G   ", "T   "))+
      ggplot2::scale_x_discrete(limits = 1:nc,
                       labels=parse(text=lab.Plot))+
      ggplot2::ylim(ymin, ymax)+
      ggplot2::theme_bw()+
      ggplot2::ggtitle(title)+
      ggplot2::xlab(NULL)+
      ggplot2::ylab(yTitle)
  } else {
    if (is.null(lengths)) {
      stop('If regs is specified, lengths for the regions must also be included.')
    } else if (is.null(regs)) {
      numReg = length(lengths[[1]])
      regs = list(as.character(1:numReg))
    }
    rects = getRegions(regs, lengths)
    psamRatioPlot = ggplot2::ggplot(values.CompPlot)+
      ggplot2::geom_rect(data = rects, 
                         aes(xmin = xstart,
                             xmax = xend,
                             ymin = -Inf,
                             ymax = Inf, 
                             fill = col), 
                         alpha = 0.4)+
      ggplot2::geom_point(data = values.CompPlot, aes(x=PosID, y=Affinity.Ratio, colour=N), size = 1, position = position_dodge(width = .2)) +
      ggplot2::geom_errorbar(data = values.CompPlot, aes(x=PosID, ymax = PlotErrorMax, ymin=PlotErrorMin, colour = N), width = 1, position = position_dodge(width = .2))+
      ggplot2::scale_colour_manual(name = "", values=c("green","blue","#FFC000","red"), breaks=c("A", "C", "G", "T"),
                                   labels=c("A           ", "C   ", "G   ", "T   "))+
      ggplot2::scale_x_discrete(limits = 1:nc,
                                labels=parse(text=lab.Plot))+
      ggplot2::ylim(ymin, ymax)+
      ggplot2::theme_bw()+
      ggplot2::ggtitle(title)+
      ggplot2::xlab(NULL)+
      ggplot2::ylab(yTitle)
  }
  
  
  return(psamRatioPlot)
}

#' Scatter plot showing nucleotide ddG values from 'model1' vs. nucleotide ddG values from 'model2'.
#' 
#' Returns standardized scatter plot comparing two sets of mononucleotide features.
#' 
#' @param model1 Model1 from which to take first set of nucleotide feature values.
#' @param model2 Model2 from which to take first set of nucleotide feature values.
#' @param title Title for plot.
#' @param xTitle Xlabel for plot.
#' @param yTitle Ylabel for plot.
#' @param iter Iteration of model to use for mono nucleotide values.
#' @param regs Optional parameter giving labels to use for region shading. If used, lengths parameter must also be specified.
#' @param lengths Length of region parameters to be included.
#' @export
plotPSAMComp.Scatter = function(model1, model2, title="", xTitle="", yTitle = "", iter=NULL) {
  values.Plot1 = getPlotValues(model1@features@N, iteration = iter)
  values.Plot2 = getPlotValues(model2@features@N, iteration = iter)
  values.CompPlot = merge(values.Plot1, values.Plot2, by = c("PosID", "N"), all = TRUE, sort = FALSE)
  values.CompPlot$ddG.x = -log(values.CompPlot$Affinity.x)
  values.CompPlot$ddG.y = -log(values.CompPlot$Affinity.y)
  values.CompPlot$PlotErrorMax.x = values.CompPlot$ddG.x+values.CompPlot$SE.x
  values.CompPlot$PlotErrorMax.y = values.CompPlot$ddG.y+values.CompPlot$SE.y
  values.CompPlot$PlotErrorMin.x = values.CompPlot$ddG.x-values.CompPlot$SE.x
  values.CompPlot$PlotErrorMin.y = values.CompPlot$ddG.y-values.CompPlot$SE.y
  
  if (model1@rcSymmetric == TRUE) {
    colMax = (ncol(model1@features@N@N.values) %/% 2)+1
    values.CompPlot = values.CompPlot[values.CompPlot$PosID <= colMax,]
  }
  
  xMinPlotPoint = min(values.CompPlot$PlotErrorMin.x)
  xMaxPlotPoint = max(values.CompPlot$PlotErrorMax.x)
  yMinPlotPoint = min(values.CompPlot$PlotErrorMin.y)
  yMaxPlotPoint = max(values.CompPlot$PlotErrorMax.y)
  xmin = ((xMinPlotPoint %/% .01)-2)*.01
  xmax = ((xMaxPlotPoint %/% .01)+3)*.01
  ymin = ((yMinPlotPoint %/% .01)-2)*.01
  ymax = ((yMaxPlotPoint %/% .01)+3)*.01
  #maxPlotPoint = max(values.CompPlot$PlotErrorMax)
  
  psamScatterPlot = ggplot2::ggplot(values.CompPlot)+
    ggplot2::geom_abline(slope = 1, intercept = 0, col = "black")+
    ggplot2::geom_point(data = values.CompPlot, aes(x=ddG.x, y=ddG.y, colour=N), size = 1, pch=16) +
    ggplot2::geom_errorbar(data = values.CompPlot, aes(x=ddG.x, ymax = PlotErrorMax.y, ymin=PlotErrorMin.y, colour = N), width = .01)+
    ggplot2::geom_errorbarh(data = values.CompPlot, aes(x=ddG.x, y = ddG.y, xmax = PlotErrorMax.x, xmin=PlotErrorMin.x, colour = N), height = .01)+
    ggplot2::scale_colour_manual(name = "", values=c("green","blue","#FFC000","red"), breaks=c("A", "C", "G", "T"),
                                 labels=c("A           ", "C   ", "G   ", "T   "))+
    ggplot2::coord_cartesian(xlim=c(xmin,xmax),
                             ylim=c(ymin,ymax),
                             expand=FALSE)+
    ggplot2::theme_bw()+
    ggplot2::ggtitle(title)+
    ggplot2::xlab(xTitle)+
    ggplot2::ylab(yTitle)
  
  return(psamScatterPlot)
}


#' Add shape profiles for k-mers. 
#' 
#' Adds Shape parameter values at all offsets of k-mers contained in data table input. 
#' 
#' @param data Data table containing sequences to be profiled. Must include a "Kmer" variable column with character objects.
#' @param shapeTable Table with shape parameters to be used for scoring.
#' @param shape Shape parameter to be used in scoring.
#' @export
addShapeValues = function(data, shapeTable, shape) {
  klen = nchar(data$Kmer[1])
  shapeParamLen = nchar(row.names(shapeTable)[1])
  for (i in 1:(klen-shapeParamLen+1)) {
    sLabel = paste("Shape.", shape, i, sep = "")
    data[,sLabel] = shapeTable[stringi::stri_sub(data$Kmer, i, (i+shapeParamLen-1)), shape]
  }
  return(data)
}


#' Set ggplot panel size
#' 
#' Sets absolute panel size of ggplot objects (not including labels).
#' 
#' @param p Ggplot object.
#' @param width Output width of ggplotGrob object.
#' @param height Output height of ggplotGrob object.
#' @export
set_panel_size <- function(p=NULL, width=ggplot2::unit(3, "cm"), height=ggplot2::unit(3, "cm")){
  g = ggplot2::ggplotGrob(p)
  panelIndex.w = g$layout$l[g$layout$name=="panel"]
  panelIndex.h = g$layout$t[g$layout$name=="panel"]
  g$widths[[panelIndex.w]] <- width
  g$heights[[panelIndex.h]] <- height
  class(g) <- c("fixed", class(g), "ggplot")
  g
}

#' Get Vertical Plot Height for object of class model
#' 
#' Gets vertical plot height to match output of plot(model)
#' 
#' @param model Object of class \linkS4class{model}.
#' @param iter Optional parameter indicating iteration to use for plotting. Default is latest iteration value in model.
#' @export
verticalPlot_height <- function(model, iter) {
  vPh = 0
  if (missing(iter)) {
    iter = model@iteration
  }
  if ((model@useFixedValuesOffset.N == TRUE) | (!all(model@features@N@N.set == 0))) {
    vPh = vPh + 4
  } else if ((model@useFixedValuesOffset.N == FALSE) & (iter == 0)) {
    vPh = vPh + 4
  }
  if (model@includeView == TRUE) {
    vPh = vPh +3
  } else if (model@includeDNAstrand == TRUE) {
    vPh = vPh + 2
  }
  numR = length(model@rounds[[1]])
  if (numR > 1) { 
    vPh = vPh + 2
  } 
  if (model@includeShape == TRUE) {
    if ((!all(model@features@Shape@Shape.set == 0)) | (model@useFixedValuesOffset.Shape == TRUE)) {
      vPh = vPh+3*(nrow(model@features@Shape@Shape.values)) 
    } else if ((model@useFixedValuesOffset.Shape == FALSE) & (iter == 0)) {
      vPh = vPh+3*(nrow(model@features@Shape@Shape.values))
    }
  }
  if ((numR <= 1) & (model@includeView == FALSE) & (model@includeDNAstrand == FALSE)) {
    vPh = vPh + 2
  }
  return (vPh)
}

