#' use BiocFileCache to retrieve and keep an image of the tsv file
#' @import BiocFileCache
#' @param url character(1) url to use
#' @param cache BiocFileCache::BiocFileCache instance
#' @param \dots passed to bfcadd
#' @note will If query of cache with 'ebi.ac.uk/gwas' returns 0-row tibble,
#' will populate cache with bfcadd.  Uses readr::read_tsv on cache content to return tibble.
#' The etag field does not seem to be used at EBI, thus user must check for updates.
#' @export
get_cached_gwascat = function(url="http://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
                        cache=BiocFileCache::BiocFileCache(), ...) {
  chk = BiocFileCache::bfcquery(cache, "ebi.ac.uk/gwas")
  if (nrow(chk)==0) {
      nca = BiocFileCache::bfcadd(url, ...)
      chk = BiocFileCache::bfcquery(cache, "ebi.ac.uk/gwas")
      }
  ans = readr::read_tsv(BiocFileCache::bfcrpath(cache)[[chk$rid]])
  attr(ans, "extractDate") = chk$access_time
  ans
}
  
