% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_dlfoo2_AnalysisWrap.R
\name{enrichmentAnalysis_fgsea}
\alias{enrichmentAnalysis_fgsea}
\title{perform Enrichment analysis fgsea using rank values}
\usage{
enrichmentAnalysis_fgsea(
  rank.list = NULL,
  run.folder = NULL,
  run.name = NULL,
  pathway.list.gmt = NULL,
  my_ontologies = c("h.all", "c2.cp.biocarta", "c2.cp.kegg", "c2.cp.reactome",
    "c3.tft", "c5.bp", "c5.cc", "c5.mf"),
  entrez.or.symbol = c("symbol"),
  minSize = 10,
  maxSize = 500,
  nperm = 1e+06,
  minLeadingEdge = 3,
  adjp.filter = 0.05,
  leadingEdgeMin = 4
)
}
\arguments{
\item{rank.list, }{list of named numeric vectors. One for each entity (sample/group etc). Each numeric vector is named with gene identifiers. no duplicated (duplicated names allowed)}

\item{run.folder, }{name of higher level directory}

\item{run.name, }{name of run (directory created)}

\item{entrez.or.symbol}{what identifiers to use}

\item{minSize}{fgsea param min size}

\item{maxSize}{fgsea param max size}

\item{nperm}{fgsea param nperm}

\item{my.ontologies, }{what mSig or other signatues to load}
}
\value{
write results to utput directory
\itemize{
  \item _fgsea_allOntologies: tab delimited text file with all significant ontologies.
}
}
\description{
https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

https://bioconductor.org/packages/release/bioc/manuals/fgsea/man/fgsea.pdf

updated 20190916. Saveing complete list as well as significant ontologies
}
\seealso{
Other enrichment: 
\code{\link{enrichmentAnalysis_ClusterProfiler}()},
\code{\link{enrichmentPlot_bubble}()},
\code{\link{gsea_RankList}()}

Other gene ontology: 
\code{\link{enrichmentAnalysis_ClusterProfiler}()},
\code{\link{gsea_RankList}()}

Other msig: 
\code{\link{enrichmentPlot_bubble}()},
\code{\link{gsea_RankList}()}

Other gsea: 
\code{\link{enrichmentPlot_bubble}()},
\code{\link{geneRankPlot}()},
\code{\link{gsea_RankList}()}

Other fgsea: 
\code{\link{enrichmentPlot_bubble}()},
\code{\link{fgseaCoocPlotrWrap}()}
}
\concept{enrichment}
\concept{fgsea}
\concept{gene ontology}
\concept{gsea}
\concept{msig}
