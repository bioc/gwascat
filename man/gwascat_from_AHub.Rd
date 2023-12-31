% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/updater.R
\name{gwascat_from_AHub}
\alias{gwascat_from_AHub}
\title{grab an image of EBI GWAS catalog from AnnotationHub}
\usage{
gwascat_from_AHub(tag = "AH91571", simple = FALSE, fixNonASCII = TRUE)
}
\arguments{
\item{tag}{character(1) defaults to "AH91571" which is the 3.30.2021 image}

\item{simple}{logical(1) if TRUE, just returns data.frame as retrieved from EBI; defaults to FALSE}

\item{fixNonASCII}{logical(1) if TRUE, use iconv to identify and eliminate non-ASCII content}
}
\value{
If `simple`, a data.frame is returned based on TSV data produced by EBI.  Otherwise, non-ASCII
content is processed according to the value of
`fixNonASCII` and a `gwaswloc` instance is returned, which has a concise
show method.  This can be coerced to a simple
GRanges instance with as(..., "GRanges").  The reference build is GRCh38.
}
\description{
grab an image of EBI GWAS catalog from AnnotationHub
}
\examples{
gwcat = gwascat_from_AHub()
gwcat
}
