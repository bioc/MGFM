\name{getHtmlpage}
\alias{getHtmlpage}
\title{Function to build HTML pages to show marker genes

}
\description{
    Function to build HTML pages to show marker genes
}
\usage{
    getHtmlpage(markers.list,chip, directory=getwd())
}

\arguments{
  \item{markers.list}{
    List of marker genes, returned by the function \code{\link{getMarkerGenes}}.
}

  \item{chip}{
    Chip name.
}
  \item{directory}{
    Path to the directory where to save the html pages.}

}
\details{
    This function is based on the function \link[annotate]{htmlpage} from the R-package
'annotate'.
}
\value{
    This function is used only for the side effect of creating HTML tables.
}

\author{
Khadija El Amrani <khadija.el-amrani@charite.de>
}

\examples{
data("ds2.mat")
res.list <- getMarkerGenes(ds2.mat, samples2compare="all", annotate=TRUE, chip="hgu133a",
score.cutoff=1)
getHtmlpage(res.list, chip="hgu133a",directory=getwd())
}

\keyword{ marker genes }
\keyword{ microarray data }
