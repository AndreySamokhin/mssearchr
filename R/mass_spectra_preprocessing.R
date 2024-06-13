#==============================================================================#
#' Pre-process mass spectra
#'
#' @description
#'   Pre-process mass spectra. Pre-processing includes rounding/binning,
#'   sorting, and normalization.
#'
#' @param msp_objs
#'   A list of nested lists. Each nested list is a mass spectrum. Each nested
#'   list must contain at least three elements: (1) the \code{name} element
#'   (a string) - compound name (or short description); (2) the \code{mz}
#'   element (a numeric/integer vector) - m/z values of mass spectral peaks; (3)
#'   the \code{intst} (a numeric/integer vector) - intensities of mass
#'   spectral peaks.
#' @param bin_boundary
#'   A numeric value. The position of a bin boundary (it can be considered as a
#'   'rounding point'). The \code{bin_boundary} argument must be in the
#'   following range: \code{0.01 <= bin_boundary <= 0.99}. The default value is
#'   \code{0.649}. This value is used in the AMDIS software and it is close to
#'   the optimal rounding rule proposed in our research (Khrisanfov, M.;
#'   Samokhin, A. A General Procedure for Rounding m/z Values in Low‐resolution
#'   Mass Spectra. Rapid Comm Mass Spectrometry 2022, 36 (11), e9294.
#'   https://doi.org/10.1002/rcm.9294).
#' @param remove_zeros
#'   An integer value. If \code{TRUE}, all m/z values with zero intensity are
#'   excluded from mass spectra. It should be taken into account that all
#'   zero-intensity peaks presented in a mass spectrum are considered as 'trace
#'   peaks' in the case of MS Search software. As a result, the presence/absence
#'   of such peaks can influence the value of the match factor.
#'
#' @details
#'   Pre-processing includes the following steps:
#'   \itemize{
#'     \item Calculating a nominal mass spectrum. All floating point m/z values
#'     are rounded to the nearest integer using the value of the
#'     \code{bin_boundary} argument. Intensities of peaks with identical m/z
#'     values are summed.
#'     \item Intensities of mass spectral peaks are normalized to 999 (as it is
#'     done in the MS Search software).
#'     \item Intensities of mass spectral peaks are rounded to the nearest
#'     integer.
#'     \item If the \code{remove_zeros} argument is \code{TRUE}, all
#'     zero-intensity peaks are removed from the mass spectrum.
#'     \item The \code{preprocessed} attribute is added and set to \code{TRUE}
#'     for the respective mass spectrum.
#'   }
#'
#' @return
#'   A list of nested lists. Each nested list is a mass spectrum. Only the
#'   \code{mz} and \code{intst} elements of each nested list are
#'   modified during the pre-processing step.
#'
#' @examples
#' # Original mass spectra of chlorine and methane
#' msp_objs <- list(
#'   list(name = "Chlorine",
#'        mz = c(34.96885, 36.96590, 69.93771, 71.93476, 73.93181),
#'        intst = c(0.83 * c(100, 32), c(100, 63.99, 10.24))),
#'   list(name = "Methane",
#'        mz = c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19),
#'        intst = c(0, 0, 25, 75, 155, 830, 999, 10, 0, 0))
#' )
#' matrix(c(msp_objs[[1]]$mz, msp_objs[[1]]$intst), ncol = 2) # Chlorine
#' matrix(c(msp_objs[[2]]$mz, msp_objs[[2]]$intst), ncol = 2) # Methane
#'
#' # Pre-processed mass spectra of chlorine and methane
#' pp_msp_objs <- PreprocessMassSpectra(msp_objs, remove_zeros = TRUE)
#' matrix(c(pp_msp_objs[[1]]$mz, pp_msp_objs[[1]]$intst), ncol = 2) # Chlorine
#' matrix(c(pp_msp_objs[[2]]$mz, pp_msp_objs[[2]]$intst), ncol = 2) # Methane
#'
#' @export
#'
#==============================================================================#
PreprocessMassSpectra <- function(msp_objs,
                                  bin_boundary = 0.649,
                                  remove_zeros = TRUE) {

  #--[ Input check ]------------------------------------------------------------

  # 'msp_objs'
  if (!is.list(msp_objs) || !is.list(msp_objs[[1]]) ||
      is.null(msp_objs[[1]]$mz) || is.null(msp_objs[[1]]$intst) ||
      is.null(msp_objs[[1]]$name)) {
    stop("'msp_objs' is invalid.")
  }

  # 'bin_boundary'
  if (!is.numeric(bin_boundary) || length(bin_boundary) != 1) {
    stop("'bin_boundary' must be a numeric value.")
  }
  if (bin_boundary < 0.01 || bin_boundary > 0.99) {
    stop("'bin_boundary' must be in the range from 0.01 to 0.99.")
  }

  # 'remove_zeros'
  if (!is.logical(remove_zeros) || length(remove_zeros) != 1) {
    stop("'remove_zeros' must be a logical value.")
  }



  #--[ Pre-processing ]---------------------------------------------------------

  for (i in seq_along(msp_objs)) {

    if (is.null(attr(msp_objs[[i]], "preprocessed")) ||
        !attr(msp_objs[[i]], "preprocessed")) {

      if (length(msp_objs[[i]]$mz) != length(msp_objs[[i]]$intst)) {
        stop("'length(msp_objs[[", i, "]]$mz) != ",
             "length(msp_objs[[", i, "]]$intst)'")
      }

      # Both 'mz' and 'intst' can be overwritten after initializing.
      mz_int <- as.integer(ceiling(msp_objs[[i]]$mz - bin_boundary))
      mz <- sort(unique(mz_int))
      if (length(mz) == length(mz_int)) {
        intst <- msp_objs[[i]]$intst[order(mz_int)]
      } else {
        intst <- vapply(mz, function(x) {
          sum(msp_objs[[i]]$intst[x == mz_int])
        }, numeric(1L))
      }

      # Mass spectra are normalized to 999 in the MS Search (NIST) software.
      intst <- as.integer(999 * (intst / max(intst)) + 0.5)

      if (remove_zeros) {
        mask <- (intst > 0L)
        intst <- intst[mask]
        mz <- mz[mask]
      }

      msp_objs[[i]]$mz <- mz
      msp_objs[[i]]$intst <- intst
      attr(msp_objs[[i]], "preprocessed") <- TRUE
    }
  }



  #--[ Output ]-----------------------------------------------------------------

  return(invisible(msp_objs))
}




