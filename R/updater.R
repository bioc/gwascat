
#' read NHGRI GWAS catalog table and construct associated GRanges instance
#' records for which clear genomic position cannot be determined are dropped
#' from the ranges instance
#' an effort is made to use reasonable data types for GRanges metadata, so some
#' qualifying characters such as (EA) in Risk allele frequency field will
#' simply be omitted during coercion of contents of that field to numeric.
#' @importFrom GenomeInfoDb genome genome<- seqnames seqlevelsStyle seqlevelsStyle<- seqlevels seqlevels<- seqinfo seqinfo<-
#' @import readr
#' 
#' @param table.url string identifying the .txt file curated at EBI/EMBL
#' @param fixNonASCII logical, if TRUE, non-ASCII characters as identified by
#' iconv will be replaced by asterisk
#' @param genome character string: 'GRCh38' is default and yields current image
#' as provided by EMBL/EBI; 'GRCh37' yields a realtime liftOver to hg19
#' coordinates, via AnnotationHub storage of the chain files. Any other value
#' yields an error.
#' @param withOnt logical indicating whether 'alternative' (ontology-present,
#' includes repetition of loci with one:many ontological mapping) or 'full'
#' (ontology-absent, one record per locus report) version of distributed table
#' @return a GRanges instance
#' @author VJ Carey
#' @keywords models
#' @examples
#' 
#' # if you have good internet access
#'   if (interactive()) {
#'      newcatr = makeCurrentGwascat()
#'      newcatr
#'      }
#' 
#' @export makeCurrentGwascat
makeCurrentGwascat = function(table.url=
  "http://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
   fixNonASCII=TRUE, genome="GRCh38", withOnt=TRUE) {
 stopifnot(genome %in% c("GRCh37", "GRCh38"))
 tf = tempfile()
 if (!withOnt) table.url = sub("alternative", "full", table.url)
 tst = try(download.file(table.url, destfile=tf))
 if (inherits(tst, "try-error")) stop("could not complete download")

 ct = readr::cols(
  .default = col_character(),
  `DATE ADDED TO CATALOG` = col_date(format = ""),
  PUBMEDID = col_double(),
  DATE = col_date(format = ""),
  CHR_ID = col_character(),
  CHR_POS = col_double(),
  UPSTREAM_GENE_DISTANCE = col_double(),
  DOWNSTREAM_GENE_DISTANCE = col_double(),
  MERGED = col_double(),
  SNP_ID_CURRENT = col_double(),
  INTERGENIC = col_double(),
  `P-VALUE` = col_double(),
  PVALUE_MLOG = col_double(),
  `OR or BETA` = col_double()
)

 tab = readr::read_tsv(tf, col_types=ct) #, sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
 message(paste0("formatting gwaswloc instance..."))
 tab = as.data.frame(tab)
 if (fixNonASCII) tab = fixNonASCII(tab)
 cur_plus = gwdf2GRanges(tab, extractDate=as.character(Sys.Date()))
 cur = cur_plus$okrngs
 nogr = cur_plus$nogr
 seqlevelsStyle(cur) = "NCBI"
 cursn = seqlevels(cur)
 data(si.hs.38)  
 seqinfo(cur) = si.hs.38[cursn]
 if (genome == "GRCh37") cur = lo38to19(cur)
 metadata(cur) = list(
    date.created = date(),
    creation = match.call(),
    badpos = nogr,
    sessInfo.creation = sessionInfo()
    )
 message("done.")
 cur
}

lo38to19 = function(gwwl) {
 if (!requireNamespace("AnnotationHub")) stop("install AnnotationHub to use this function")
 message("starting liftover from GRCh38 to GRCh37")
 stopifnot(genome(gwwl)[1] == "GRCh38")
 ah = AnnotationHub::AnnotationHub()
 ii = AnnotationHub::query(ah, "UCSC liftOver chain file from hg38 to hg19")
 ch = ah[[names(ii)]]
 seqlevelsStyle(gwwl) = "UCSC"
 g19 = liftOver( as(gwwl, "GRanges"), ch )
 message("liftover complete.")
 e = elementNROWS(g19)
 dr = which(e==0)
 if (length(dr)>0) g19 = g19[-dr]
 g19 = unlist(g19)
 metadata(g19)$conversion = "liftOver"
 metadata(g19)$sessInfo.creation.liftOver = sessionInfo()
 seqlevelsStyle(g19) = "NCBI"
 genome(g19) = "GRCh37"
 cursn = seqlevels(g19)
 data(si.hs.37)
 seqinfo(g19) = si.hs.37[cursn]
 new("gwaswloc", extractDate=date(), g19)
}

