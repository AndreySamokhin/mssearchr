#==============================================================================#
#' Perform the library search within R
#'
#' @description
#'   Perform library search using a custom implementation of the Identity (EI
#'   Normal) or Similarity (EI Simple) algorithm. Pairwise comparison of two
#'   mass spectra is implemented in C.
#'
#' @param msp_objs_u,msp_objs_l
#'   A list of nested lists. Each nested list is a mass spectrum. Each nested
#'   list must contain at least three elements: (1) \code{name} (a string) -
#'   compound name (or short description); (2) \code{mz} (a numeric/integer
#'   vector) - m/z values of mass spectral peaks; (3) \code{intst} (a
#'   numeric/integer vector) - intensities of mass spectral peaks. Letters 'u'
#'   and 'l' stand for unknown and library respectively). Mass spectra should be
#'   pre-processed using the \code{\link{PreprocessMassSpectra}} function.
#' @param algorithm
#' A string. Library search algorithm. Either the Identity EI Normal
#' (\code{identity_normal}) or Similarity EI Simple (\code{similarity_simple})
#' algorithm.
#' @param n_hits
#'   An integer value. The maximum number of hits (i.e., candidates) to display.
#' @param hitlist_columns
#'   A character vector. Three columns are always present in the returned
#'   hitlist: \code{name}, \code{mf} (i.e., the match factor), and \code{idx}
#'   (i.e., the index of the respective library mass spectrum in the
#'   \code{msp_objs_l} list). Some additional columns can be added using the
#'   \code{hitlist_columns} argument (e.g., \code{cas_no}, \code{formula},
#'   \code{inchikey}, etc.). Only scalar values (i.e., an atomic vector of unit
#'   length) are allowed.
#' @param mz_min,mz_max
#'   An integer value. Boundaries of the m/z range (all m/z values out of this
#'   range are not taken into account when the match factor is calculated).
#' @param comments
#'   Any R object. Some additional information. It is saved as the 'comments'
#'   attribute of the returned list.
#'
#' @return
#'   Return a list of data frames. Each data frame is a hitlist (i.e., list of
#'   possible candidates). Each hitlist always contains three columns:
#'   \code{name}, \code{mf} (i.e., the match factor), and \code{idx} (i.e., the
#'   index of the respective library mass spectrum in the \code{msp_objs_l}
#'   list). Additional columns can be extracted using the \code{hitlist_columns}
#'   argument. Library search options are saved as the
#'   \code{library_search_options} attribute.
#'
#' @examples
#' # Reading the 'alkanes.msp' file
#' msp_file <- system.file("extdata", "alkanes.msp", package = "mssearchr")
#'
#' # Pre-processing
#' msp_objs_u <- PreprocessMassSpectra(ReadMsp(msp_file)) # unknown mass spectra
#' msp_objs_l <- PreprocessMassSpectra(massbank_alkanes)  # library mass spectra
#'
#' # Searching using the Identity algorithm
#' hitlists <- LibrarySearch(msp_objs_u, msp_objs_l,
#'                           algorithm = "identity_normal", n_hits = 10L,
#'                           hitlist_columns = c("formula", "smiles", "db_no"))
#'
#' # Printing a hitlist for the first compound from the 'alkanes.msp' file
#' print(hitlists[[1]][1:5, ])
#'
#' #>        name       mf idx formula        smiles                db_no
#' #> 1  UNDECANE 950.5551  11  C11H24   CCCCCCCCCCC MSBNK-{...}-JP006877
#' #> 2  UNDECANE 928.4884  72  C11H24   CCCCCCCCCCC MSBNK-{...}-JP005760
#' #> 3  DODECANE 905.7546  74  C12H26  CCCCCCCCCCCC MSBNK-{...}-JP006878
#' #> 4 TRIDECANE 891.7862  41  C13H28 CCCCCCCCCCCCC MSBNK-{...}-JP006879
#' #> 5  DODECANE 885.6247  42  C12H26  CCCCCCCCCCCC MSBNK-{...}-JP005756
#'
#' @export
#'
#==============================================================================#
LibrarySearch <- function(msp_objs_u,
                          msp_objs_l,
                          algorithm = c("identity_normal", "similarity_simple"),
                          n_hits = 100L,
                          hitlist_columns = c("formula", "mw", "smiles"),
                          mz_min = NULL,
                          mz_max = NULL,
                          comments = NULL) {

  #--[ Input check ]------------------------------------------------------------

  # 'msp_objs_u'
  if (!is.list(msp_objs_u) || !is.list(msp_objs_u[[1]]) ||
      is.null(msp_objs_u[[1]]$mz) || is.null(msp_objs_u[[1]]$intst)
      || is.null(msp_objs_u[[1]]$name)) {
    stop("'msp_objs_u' is invalid.")
  }

  # 'msp_objs_l'
  if (!is.list(msp_objs_l) || !is.list(msp_objs_l[[1]]) ||
      is.null(msp_objs_l[[1]]$mz) || is.null(msp_objs_l[[1]]$intst)
      || is.null(msp_objs_l[[1]]$name)) {
    stop("'msp_objs_l' is invalid.")
  }

  # 'algorithm'
  # The 'match.arg()' function is used.

  # 'n_hits'
  if (!is.null(n_hits)) {
    if (!is.numeric(n_hits) || length(n_hits) != 1L ||
        abs(n_hits - as.integer(n_hits)) > 1e-10 || as.integer(n_hits) < 1L) {
      stop("'n_hits' must be a non-negative integer value.")
    }
  }

  # 'hitlist_columns'
  if (!is.null(hitlist_columns)) {
    if (!is.character(hitlist_columns) || length(hitlist_columns) < 1L) {
      stop("'hitlist_columns' must be a character vector.")
    }
  }

  # 'mz_min'
  if (!is.null(mz_min)) {
    if (!is.numeric(mz_min) || length(mz_min) != 1L ||
        abs(mz_min - as.integer(mz_min)) > 1e-10 || as.integer(mz_min) < 0L) {
      stop("'mz_min' must be a non-negative integer value.")
    }
  }

  # 'mz_max'
  if (!is.null(mz_max)) {
    if (!is.numeric(mz_max) || length(mz_max) != 1L ||
        abs(mz_max - as.integer(mz_max)) > 1e-10 || as.integer(mz_max) < 0L) {
      stop("'mz_max' must be a non-negative integer value.")
    }
    if (!is.null(mz_min)) {
      if (as.integer(mz_min) >= as.integer(mz_max)) {
        stop("'mz_max' must be greater than 'mz_min'.")
      }
    }
  }

  # 'comments'
  # Any R object can be passed.



  #--[ Modifying input arguments ]----------------------------------------------

  algorithm <- match.arg(algorithm)
  is_identity <- (algorithm == "identity_normal")
  if (is.null(mz_min)) {
    mz_min <- -1L
  }
  if (is.null(mz_max)) {
    mz_max <- -1L
  }
  hitlist_columns <- unique(hitlist_columns[hitlist_columns != "name"])



  #--[ Processing ]-------------------------------------------------------------

  wms_u <- .CalculateWeightedIntensities(msp_objs_u, algorithm)
  wms_l <- .CalculateWeightedIntensities(msp_objs_l, algorithm)

  out <- lapply(seq_along(msp_objs_u), function(idx_u) {
    mf <- vapply(seq_along(msp_objs_l), function(idx_l) {
      .CalcMatchFactorC(wms_u[[idx_u]], wms_l[[idx_l]], is_identity,
                        mz_min, mz_max)
    }, numeric(1))
    new_order <- order(mf, decreasing = TRUE)
    # 'hl' is an abbreviation for 'hitlist'
    hl_idxs <- new_order[seq(1L, min(n_hits, length(new_order)))]
    hl_names <- vapply(msp_objs_l[hl_idxs], function(x) {
      x$name
    }, character(1))
    hl_mfs <- mf[hl_idxs]
    addl_columns <- lapply(hitlist_columns, function(field_name) {
      unlist(lapply(msp_objs_l[hl_idxs], function(x) {
        if (!is.atomic(x[[field_name]]) || length(x[[field_name]]) > 1L) {
          stop(paste0("'", field_name, "'", " is not a scalar value."))
        }
        if (is.null(x[[field_name]])) {
          return(NA)
        } else {
          return(x[[field_name]])
        }
      }), recursive = FALSE)
    })
    hl <- as.data.frame(
      c(list(hl_names), list(hl_mfs), list(hl_idxs), addl_columns),
      col.names = c("name", "mf", "idx", hitlist_columns)
    )
    return(hl)
  })



  #--[ Attributes ]-------------------------------------------------------------

  library_search_options <- list(algorithm = algorithm,
                                 mz_min = mz_min,
                                 mz_max = mz_max,
                                 n_hits = n_hits)
  attr(out, "library_search_options") <- library_search_options
  attr(out, "comments") <- comments



  #--[ Output ]-----------------------------------------------------------------

  return(out)
}





#==============================================================================#
#' Calculate weighted intensities
#'
#' @description
#'   Calculate weighted intensities.
#'
#' @param msp_objs
#'   A list of nested lists. See the \code{\link{LibrarySearch}} function.
#' @inheritParams LibrarySearch
#'
#' @return
#'   Return a matrix with three or four columns:
#'   \describe{
#'     \item{\code{mz}}{Mass values.}
#'     \item{\code{intst}}{Original intensities.}
#'     \item{\code{w1}}{Weighted intensities which are used to calculate the
#'     cosine similarity (for the Identity algorithm: \code{w1 = sqrt(mz *
#'     intst)}; for the Similarity algorithm: \code{w1 = sqrt(intst)}).}
#'     \item{\code{w2}}{Weighted intensities which are used to calculate ratio
#'     of adjacent peaks (it is calculated only for the Identity algorithm:
#'     \code{w2 = sqrt(intst)}).}
#'   }
#'
#' @noRd
#'
#==============================================================================#
.CalculateWeightedIntensities <- function(msp_objs,
                                          algorithm) {

  #--[ Input check ]------------------------------------------------------------

  # It is a hidden function. Any check is not performed.



  #--[ Processing ]-------------------------------------------------------------

  wms <- lapply(msp_objs, function(msp) {

    if (is.null(attr(msp, "preprocessed")) || !attr(msp, "preprocessed")) {
      stop("'msp_objs' must be preprocessed using 'PreprocessMassSpectra()'.")
    }

    if (algorithm == "identity_normal") {
      temp <- c(msp$mz, msp$intst, sqrt(msp$intst * msp$mz), sqrt(msp$intst))
      return(matrix(temp, ncol = 4L))

    } else if (algorithm == "similarity_simple") {
      return(matrix(c(msp$mz, msp$intst, sqrt(msp$intst)), ncol = 3L))
    } else {
      stop("'algorithm' takes an unsupported value.")
    }
  })



  #--[ Output ]-----------------------------------------------------------------

  return(wms)
}




