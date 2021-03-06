% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_dlfoo2_PlotWrap.R
\name{coOccurence_bubble}
\alias{coOccurence_bubble}
\title{Generic plot for plottng co-occurrence freuency matrix using Bubble blot.}
\usage{
coOccurence_bubble(
  my.tab = NULL,
  my.columns = NULL,
  color.key = "size",
  title.text = NULL,
  order.by.clustering.x = T,
  order.by.clustering.y = T,
  border.col = NA,
  reverse.x = F,
  reverse.y = F
)
}
\arguments{
\item{my.tab, }{data frame with two grouping values to compare overlap. Should be factors to set order}

\item{my.columns}{two columns with group assignments that should be compared in frequency matrix. Y-axis and X-axis, respectively}

\item{color.key}{color key to match to primary column (y-axis)}

\item{order.by.clustering.x, }{let order of terms/groups be decided by clustering. Otherwise defined by levels}

\item{order.by.clustering.y, }{let order of terms/groups be decided by clustering. Otherwise defined by levels}

\item{reverse.x}{if to reverse x-axis}

\item{reverse.y}{if to reverse y-axis}
}
\value{
ggplot object
}
\description{
Framction of y-axis group members (defined by my.columns[1]) that fall within x-axis groups, i.e. Sum of each plot row fraction equals 1.

v1.0 20191025.
}
\seealso{
Other plotwrap: 
\code{\link{enrichmentPlot_bubble}()}
}
\concept{Co-occurrence}
\concept{plotwrap}
