% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ldtagr.R
\name{ldtagr}
\alias{ldtagr}
\title{expand a list of variants by including those in a VCF with LD exceeding some
threshold; uses snpStats ld()}
\usage{
ldtagr(
  snprng,
  tf,
  samples,
  genome = "hg19",
  lbmaf = 0.05,
  lbR2 = 0.8,
  radius = 1e+05
)
}
\arguments{
\item{snprng}{a named GRanges for a single SNP.  The name must correspond to
the name that will be assigned by genotypeToSnpMatrix (from VariantTools) to
the corresponding column of a SnpMatrix.}

\item{tf}{TabixFile instance pointing to a bgzipped tabix-indexed VCF file}

\item{samples}{a vector of sample identifiers, if excluded, all samples used}

\item{genome}{tag like 'hg19'}

\item{lbmaf}{lower bound on variant MAF to allow consideration}

\item{lbR2}{lower bound on R squared for regarding SNP to be incorporated}

\item{radius}{radius of search in bp around the input range}
}
\value{
a GRanges with names corresponding to 'new' variants and mcols
fields 'paramRangeID' (base variant input) and 'R2'
}
\description{
expand a list of variants by including those in a VCF with LD exceeding some
threshold; uses snpStats ld()
}
\note{
slow but safe approach.  probably a matrix method could be substituted
using the nice sparse approach already in snpStats
}
\examples{
 cand = GenomicRanges::GRanges("1", IRanges::IRanges(113038694, width=1))
 names(cand) = "rs883593"
 requireNamespace("VariantAnnotation")
 expath = dir(system.file("vcf", package="gwascat"), patt=".*exon.*gz$", full=TRUE)
 tf = Rsamtools::TabixFile(expath)
 ldtagr( cand, tf, lbR2 = .8)
}
\author{
VJ Carey
}
\keyword{models}
