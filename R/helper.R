#' Check if rowname is duplicated and add binary is_duplicated
#' Helper function
#'
#' @import magrittr
#' @export
#' @param dat data frame
#' @param rowname rowname to check as character
add_is_duplicated <- function(dat, rowname) {
  dat[paste0(rowname, '_is_duplicated')] <- duplicated(dat[rowname]) | duplicated(dat[rowname], fromLast=TRUE)
  dat %>% return()
}
