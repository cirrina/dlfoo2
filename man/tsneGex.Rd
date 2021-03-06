% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_dlfoo2_AnalysisWrap.R
\name{tsneGex}
\alias{tsneGex}
\title{function tsneGex wrapper for RtSNE}
\usage{
tsneGex(
  eset,
  perplex = seq(5, 50, by = 5),
  runName,
  runFolder,
  doPlotTaxonomy = F,
  theta = 0.5,
  max_iter = 5000,
  floor_value = NULL,
  z_score = F,
  log = F,
  variance = NULL
)
}
\arguments{
\item{eset}{expressionSet object}

\item{perplex}{vector of integers. Suggested max is n-1  3}

\item{runName}{name of plot/analyses - subfolder}

\item{runFolder}{higher level directory where to plot/save objects}

\item{doPlotTaxonomy}{if to load rcc phenotype from file & plot}

\item{theta}{Rtsne theta value}

\item{max_iter}{Rtsne max_iter value}

\item{floor_value}{Rtsne floor_value}

\item{z_score}{if to use z-scores (False default)}

\item{log}{if to log matrix(False default)}

\item{variance}{numeric between 0 and 1 describing the ratio of Features to keep (0.8 keeps the 80 percent highest vatying Features)}
}
\value{
write tables and tSNEres to file
}
\description{
Version 1.0 20190418

tSNE Aanalysis using \code{\link[Rtsne]{Rtsne}} function in stats package.
}
\seealso{
Other tsne: 
\code{\link{plotViolinScatterRankBarGeneric}()},
\code{\link{scatterPoint_annotationPlot_SHINY}()},
\code{\link{scatterPoint_annotationPlot}()},
\code{\link{tsne2df}()},
\code{\link{tsneMe}()},
\code{\link{tsnePlot}()},
\code{\link{tsneShinyTaxEvaluation}()}

Other analysis: 
\code{\link{ConsensusClusterPlus.dl}()},
\code{\link{ConsensusCluster_process}()},
\code{\link{corrCcle}()},
\code{\link{corrGtex}()},
\code{\link{meAnalysisCorrGexRandom}()},
\code{\link{meAnalysisCorrGex}()},
\code{\link{tsneMe}()}
}
\concept{analysis}
\concept{tsne}
