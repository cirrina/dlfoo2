% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_dlfoo2_FeatureAnnotation.R
\name{features2genomicRanges}
\alias{features2genomicRanges}
\title{Function to cerate Genomic Ranges fdata object of gene coordinates using Biomart
previously genes2genomicRanges and probes2genomicRanges}
\usage{
features2genomicRanges(
  ids = "all",
  id_type = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  collapse.to.unique = T,
  primary_id = NULL,
  genome.version = c("hg19", "hg38"),
  keep.chromosomes = c(1:22, "X", "Y", "MT")
)
}
\arguments{
\item{ids}{Input identfier - A vector of >1 gene ensembl ids. if 'all' then all genes in mart are retrieved}

\item{id_type}{what identifier type for biomart lookup - one of c("ensembl_gene_id","hgnc_symbol", "entrezgene")}

\item{collapse.to.unique}{If output should contain uniqe ids (per row) or if multiple rows are allowed for out. The non-collapsed table will be stored in metadata slot}

\item{primary_id}{If collapse.to.unique = T. This id type will be used as primary (per row). Other Ids will be piped. If multiople regions are present per primary region - The first instance will be acquired.}

\item{genome.version}{what genome should be used - one of c('hg19','grch38')}

\item{keep.chromosomes}{what chrosomosomes to keep. set to NULL if to keep all chromosomes}
}
\value{
GRanges object with genomic start/stop for supplied identifiers (or all if 'all' option used)
}
\description{
Function to cerate Genomic Ranges fdata object of gene coordinates using Biomart
previously genes2genomicRanges and probes2genomicRanges
}
\seealso{
Other granges: 
\code{\link{chip_genePeaks2occupancyTab}()},
\code{\link{chip_peak2gene_maxValue}()},
\code{\link{chip_peak2gene}()},
\code{\link{expandRanges}()},
\code{\link{liftoverGRanges}()},
\code{\link{loadRefgene}()},
\code{\link{meAnalysisCorTableStats}()},
\code{\link{meAnalysisCorrGexRandom}()},
\code{\link{meAnalysisCorrGex}()},
\code{\link{meMergeCpGs}()},
\code{\link{segmentPlotter}()}

Other feature annotation: 
\code{\link{entrez2chr}()},
\code{\link{entrez2symbol}()},
\code{\link{featureGetBM}()},
\code{\link{homoloGene}()},
\code{\link{illuminav4_2fdata}()}
}
\concept{feature annotation}
\concept{granges}
