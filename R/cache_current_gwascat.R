
#' use BiocFileCache to retrieve and keep an image of the tsv file distributed by EBI
#' @import BiocFileCache
#' @importFrom readr read_tsv problems
#' @param url character(1) url to use
#' @param cache BiocFileCache::BiocFileCache instance
#' @param refresh logical(1) force download and recaching
#' @param \dots passed to bfcadd
#' @note will If query of cache with 'ebi.ac.uk/gwas' returns 0-row tibble,
#' will populate cache with bfcadd.  Uses readr::read_tsv on cache content to return tibble.
#' The etag field does not seem to be used at EBI, thus user must check for updates.
#' @return a tibble as produced by readr::read_tsv, with attributes extractDate (as
#' recorded in cache as `access_time`, and problems (a tibble returned by read_tsv).
#' @export
get_cached_gwascat = function(url="http://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
                        cache=BiocFileCache::BiocFileCache(), refresh=FALSE, ...) {
  chk = BiocFileCache::bfcquery(cache, "ebi.ac.uk/gwas")
  if (nrow(chk)==0 | refresh) {
      nca = BiocFileCache::bfcadd(cache, url, ...)
      chk = BiocFileCache::bfcquery(cache, "ebi.ac.uk/gwas")
      }
  ans = suppressMessages({ suppressWarnings({
      readr::read_tsv(BiocFileCache::bfcrpath(cache)[[rev(chk$rid)[1]]]) # use rev to get latest addition
      }) })
  pb = readr::problems(ans)
  attr(ans, "extractDate") = chk$access_time
  attr(ans, "problems") = pb
  ans
}
  
#' produce a GRanges from gwascat tibble
#' @param x a tibble from `get_cached_gwascat()`
#' @param short logical(1) if TRUE only keep selected columns in mcols
#' @param for_short character() column names to keep in mcols
#' @param genome_tag character(1) defaults to "GRCh38"
#' @export
as_GRanges = function(x, short=TRUE, for_short=c("PUBMEDID", "DATE", "DISEASE/TRAIT",
   "SNPS"), genome_tag = "GRCh38") {
  if (!requireNamespace("GenomeInfoDb")) stop("install GenomeInfoDb to use this package")
  bad = which(is.na(x$CHR_POS))
  lbad = length(bad)
  if (lbad>0) {
    message(sprintf("dropping %d records that have NA for CHR_POS\n", lbad))
    x = x[-bad,]
    }
  has_semi = grep(";", x$CHR_POS)
  lsemi = length(has_semi)
  if (lsemi>0) {
    message(sprintf("%d records have semicolon in CHR_POS; splitting and using first entry.\n", lsemi))
    npos = strsplit(x$CHR_POS, ";")
    x$CHR_POS = vapply(npos, function(z) z[1], character(1))
#    print(summary(as.numeric(x$CHR_POS)))
#    bad = which(is.na(as.numeric(x$CHR_POS)))[1:4]
#    print(as.data.frame(x[bad,]))
    }
  has_x = grep(" x ", x$CHR_POS)
  lhasx = length(has_x)
  if (lhasx > 0) {
    message(sprintf("%d records have ' x ' in CHR_POS indicating multiple SNP effects, using first.\n", lhasx))
    nx = strsplit(x$CHR_POS, " x ")
    x$CHR_POS = vapply(nx, function(z) z[1], character(1))
    }
  ans = GenomicRanges::GRanges(x$CHR_ID, IRanges(as.numeric(x$CHR_POS), width=1))
  if (short) mcols(ans) = x[,for_short]
  else mcols(ans) = x
  GenomeInfoDb::genome(ans) = genome_tag
  ans
}
#
#
#
#> colnames(x)
# [1] "DATE ADDED TO CATALOG"      "PUBMEDID"                  
# [3] "FIRST AUTHOR"               "DATE"                      
# [5] "JOURNAL"                    "LINK"                      
# [7] "STUDY"                      "DISEASE/TRAIT"             
# [9] "INITIAL SAMPLE SIZE"        "REPLICATION SAMPLE SIZE"   
#[11] "REGION"                     "CHR_ID"                    
#[13] "CHR_POS"                    "REPORTED GENE(S)"          
#[15] "MAPPED_GENE"                "UPSTREAM_GENE_ID"          
#[17] "DOWNSTREAM_GENE_ID"         "SNP_GENE_IDS"              
#[19] "UPSTREAM_GENE_DISTANCE"     "DOWNSTREAM_GENE_DISTANCE"  
#[21] "STRONGEST SNP-RISK ALLELE"  "SNPS"                      
#[23] "MERGED"                     "SNP_ID_CURRENT"            
#[25] "CONTEXT"                    "INTERGENIC"                
#[27] "RISK ALLELE FREQUENCY"      "P-VALUE"                   
#[29] "PVALUE_MLOG"                "P-VALUE (TEXT)"            
#[31] "OR or BETA"                 "95% CI (TEXT)"             
#[33] "PLATFORM [SNPS PASSING QC]" "CNV"                       
#[35] "MAPPED_TRAIT"               "MAPPED_TRAIT_URI"          
#[37] "STUDY ACCESSION"            "GENOTYPING TECHNOLOGY"     
#
