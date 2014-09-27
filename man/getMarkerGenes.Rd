\name{getMarkerGenes}
\alias{getMarkerGenes}
\title{Marker Gene Detection

}
\description{
    Function to detect marker genes using microarray gene expression data sets
}
\usage{
    getMarkerGenes(data.mat, samples2compare="all", annotate=TRUE, chip=NULL,
score.cutoff=1)
}

\arguments{
    \item{data.mat}{
    The microarray data matrix with probe sets corresponding to rows
and samples corresponding to columns.
}
    \item{samples2compare}{
    A character vector with the sample names to be compared (e.g.
c("liver", "lung", "brain")). By default all samples are used.
}
    \item{annotate}{
    A boolean value. If TRUE the gene symbol and the entrez gene id are shown.
}

    \item{chip}{
    Chip name.
}
    \item{score.cutoff}{
    A value in the interval [0,1] to filter the markers according to the
specificity score. The default value is 1 (no filtering).}

}
\details{
    For each marker in the output list, the probe set and
the corresponding score are shown. If annotate is TRUE, the gene symbol and the
entrez gene id are shown. The score is used to
rank the markers according to their specificity. The score values range from 0
to 1. 
Values near 0 would indicate high specificity and large values closer to 1 would
indicate low specificity. 
}
\value{
    A list with marker genes associated with each sample type.
}

\author{
    Khadija El Amrani <khadija.el-amrani@charite.de>
}





\examples{
data("ds2.mat")
res.list <- getMarkerGenes(ds2.mat, samples2compare="all", annotate=TRUE,
chip="hgu133a", score.cutoff=1)
names(res.list)
## show the first 20 markers of liver
res.list[["liver_markers"]][1:20]
}
\keyword{ marker genes }
\keyword{ microarray data }
