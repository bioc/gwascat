#' expand a list of variants by including those in a VCF with LD exceeding some
#' threshold
#' 
#' expand a list of variants by including those in a VCF with LD exceeding some
#' threshold
#' 
#' uses snpStats ld()
#' 
#' @param snprng a named GRanges for a single SNP.  The name must correspond to
#' the name that will be assigned by genotypeToSnpMatrix (from VariantTools) to
#' the corresponding column of a SnpMatrix.
#' @param tf TabixFile instance pointing to a bgzipped tabix-indexed VCF file
#' @param samples a vector of sample identifiers, if excluded, all samples used
#' @param genome tag like 'hg19'
#' @param lbmaf lower bound on variant MAF to allow consideration
#' @param lbR2 lower bound on R squared for regarding SNP to be incorporated
#' @param radius radius of search in bp around the input range
#' @return a GRanges with names corresponding to 'new' variants and mcols
#' fields 'paramRangeID' (base variant input) and 'R2'
#' @note slow but safe approach.  probably a matrix method could be substituted
#' using the nice sparse approach already in snpStats
#' @author VJ Carey
#' @keywords models
#' @examples
#' 
#' require(GenomicRanges)
#' if (requireNamespace("gQTLstats")) {
#'    # install gQTLstats to test this function
#'  cand = GRanges("1", IRanges(113038694, width=1))
#'  names(cand) = "rs883593"
#'  require(VariantAnnotation)
#'  expath = dir(system.file("vcf", package="gwascat"), patt=".*exon.*gz$", full=TRUE)
#'  tf = TabixFile(expath)
#'  ldtagr( cand, tf, lbR2 = .8)
#' }
#' # should do with 1000 genomes in S3 bucket and gwascat
#' 
#' @export ldtagr
ldtagr = function( snprng, tf, samples, genome="hg19",
   lbmaf=.05, lbR2=.8, radius=100000 ) {
#
# given a GRanges with name and address of a SNP with representation in a
# VCF referenced by tf, give the names and addresses of all
# SNP in LD with the original locus at R2 at least lbR2
# 
  stopifnot(length(snprng)==1)
  snpid = names(snprng)
  stopifnot(length(snpid)==1)
  if (!requireNamespace("gQTLstats")) stop("install gQTLstats to use this function")
  if (!requireNamespace("DelayedArray")) stop("install DelayedArray to use this function")
  if (!requireNamespace("snpStats")) stop("install snpStats to use this function")
  quer = gQTLstats::queryVCF( gr=snprng+radius, vcf.tf=tf, 
         samps=samples, genome=genome, getSM=TRUE )
  empty = GRanges()
  mcols(empty) = DataFrame(paramRangeID = factor(), R2 = numeric())
  vcfrng = DelayedArray::rowRanges(quer$readout)
  gt = quer$sm$genotypes
  if (!(snpid %in% colnames(gt))) {
    message(paste0("NOTE: ", snpid, " not in VCF at given radius, returning empty GRanges"))
    return(empty)
    }
  map = quer$sm$map
  cs = snpStats::col.summary(gt)
  mafok = which(cs[,"MAF"] >= lbmaf)
  stopifnot(length(mafok)>0)
  gt = gt[,mafok]
  if (!(snpid %in% colnames(gt))) {
    message(paste0("NOTE: ", snpid, " has MAF too low, returning empty GRanges."))
    return(empty)
    }
  ldvec = snpStats::ld(gt[,snpid], gt, stats="R.squared")
  extrng = vcfrng[ thec <- colnames(theld <- ldvec[, which(as.numeric(ldvec) > lbR2),drop=FALSE ]) ]
  extrng$R2 = theld
  #list(gt=gt, map=map, vcfrng=vcfrng, ldvec=ldvec,
  #   extrng = vcfrng[ colnames(ldvec[, 
  #   which(as.numeric(ldvec) > lbR2),drop=FALSE ]) ])
  names(extrng) = make.names(thec)
  extrng[, c("paramRangeID", "R2")]
}
