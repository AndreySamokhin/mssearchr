#==============================================================================#
#' Parse a single mass spectrum read from an msp-file
#'
#' @description
#'   Parse a single mass spectrum read from an msp-file.
#'
#' @param input_lines
#'   A character vector. Each element of the character vector is a line from an
#'   msp-file.
#' @param line_no
#'   An integer. The line in the original file containing the 'Name' field. It
#'   is used only to identify the position of a problematic mass spectrum in
#'   the original msp-file.
#'
#' @details
#'   It is a hidden helper, which is called from the \code{\link{ReadMsp}}
#'   function.
#'
#' @return
#'   Return a named list. Most elements of the list are strings containing some
#'   metadata (e.g., \code{name}, \code{cas_no}, \code{formula}, etc.). Mass
#'   spectral data are stored in two numeric vectors (\code{mz} and
#'   \code{intst}) containing  m/z values and intensities of mass spectral
#'   peaks. Synonyms (i.e., \code{synon}) are stored as a character vector.
#'
#' @noRd
#'
#==============================================================================#
.ParseSingleSpectrum <- function(input_lines,
                                 line_no) {

  #--[ Input check ]------------------------------------------------------------

  # It is a hidden function. Any check is not performed.



  #--[ Metadata and mass spectral data ]----------------------------------------

  mask_metadata <- grepl("^[a-zA-Z].*?:", input_lines)
  mask_ms <- (input_lines != "" & !mask_metadata)
  if (max(which(mask_metadata)) > min(which(mask_ms))) {
    stop("Incorrect msp format.")
  }



  #--[ Parsing metadata ]-------------------------------------------------------

  ptrn <- "^([a-zA-Z][^:]*):\\s(.*)$"
  field_values <- trimws(gsub(ptrn, "\\2", input_lines[mask_metadata]))
  temp <- trimws(gsub(ptrn, "\\1", input_lines[mask_metadata]))
  field_names <- gsub("\\W+", "_", gsub("#", "_no", tolower(temp)))

  # Checking that the 'Name' field is present.
  if (!is.element("name", field_names)) {
    stop("The 'name' field is absent (line_no: ", line_no, ").")
  }

  # Checking that all field names (except 'Synon') are unique.
  if (length(unique(field_names[field_names != "synon"])) !=
      sum(field_names != "synon")) {
    stop("Field names are not unique (line_no: ", line_no, ").")
  }

  # Checking that 'MZ' and 'Intst' fields are absent.
  if (is.element("mz", field_names)) {
    stop("The 'mz' field is present (line_no: ", line_no, ").")
  }
  if (is.element("intst", field_names)) {
    stop("The 'intst' field is present (line_no: ", line_no, ").")
  }

  idxs <- which(!(field_names %in% c("synon", "num_peaks")))
  msp_obj <- as.list(field_values[idxs])
  names(msp_obj) <- field_names[idxs]

  # Processing the 'NIST#' field.
  if (is.null(msp_obj$nist_no) &&
      (!is.null(msp_obj$cas_no) || !is.null(msp_obj$compound_rep))) {
    if (!is.null(msp_obj$cas_no)) {
      fied_preceding_nist_no <- "cas_no"
    } else { # i.e., '!is.null(msp_obj$compound_rep)'
      fied_preceding_nist_no <- "compound_rep"
    }
    if (grepl("NIST#:", msp_obj[[fied_preceding_nist_no]])) {
      # According to the format specification, the 'NIST#' filed is present on
      # the same line with the 'CAS#' field:
      #   CAS#: 00-00-0;  NIST#: 000000
      temp <- sub(";\\s*nist#:\\s+", ";", msp_obj[[fied_preceding_nist_no]],
                  ignore.case = TRUE)
      temp <- unlist(strsplit(temp, ";", fixed = TRUE))
      # The 'ignore.case' argument is not available for 'strsplit()'.
      msp_obj[[fied_preceding_nist_no]] <- temp[[1]]
      msp_obj$nist_no <- temp[[2]]
    }
  }

  # Processing the 'Synon' fields.
  if (is.element("synon", field_names)) {
    msp_obj$synon <- field_values[field_names == "synon"]
  }



  #--[ Parsing mass spectral data ]---------------------------------------------

  temp <- paste(input_lines[mask_ms], collapse = " ")
  temp <- trimws(gsub("[ \\t,;:\\[\\]{}()]+", " ", temp, perl = TRUE))
  raw_ms <- unlist(strsplit(temp, " ", fixed = TRUE))
  msp_obj$mz <- as.numeric(raw_ms[c(TRUE, FALSE)])
  # if (max(abs(msp_obj$mz - as.integer(msp_obj$mz))) < 1e-10) {
  #   msp_obj$mz <- as.integer(msp_obj$mz)
  # }
  msp_obj$intst <- as.numeric(raw_ms[c(FALSE, TRUE)])

  # Checking that the 'Num Peaks' field is present and the specified number of
  # peaks is correct.
  if (!is.element("num_peaks", field_names)) {
    warning("The 'Num Peaks' field is absent (line_no: ", line_no, ").")
  } else {
    num_peaks <- as.numeric(field_values[[match("num_peaks", field_names)]])
    if (length(msp_obj$mz) != num_peaks) {
      warning("The 'Num Peaks' field is incorrect (line_no: ", line_no, ").")
    }
  }



  #--[ Output ]-----------------------------------------------------------------

  return(msp_obj)
}





#==============================================================================#
#' Read mass spectra from an msp-file (NIST format)
#'
#' @description
#'   Read an msp-file containing mass spectra in the NIST format. The complete
#'   description of the format can be found in the NIST Mass Spectral Search
#'   Program manual. A summary is presented below in the "Description of the
#'   NIST format" section.
#'
#' @param input_file
#'   A string. The name of a file.
#'
#' @details
#'   Data from an msp-file are read without any modification (e.g., the order of
#'   mass values is not changed, zero-intensity peaks are preserved, etc.).
#'
#' @return
#'   Return a list of nested lists. Each nested list is a mass spectrum. Almost
#'   all metadata fields (e.g., "Name", "CAS#", "Formula", "MW", etc.) are
#'   represented as strings. All "Synon" fields are merged into a single
#'   character vector. Mass values and intensities are represented as numeric
#'   vectors (\code{mz} and \code{intst}). Names of fields are slightly
#'   modified:
#'   \itemize{
#'     \item names are converted to lowercase;
#'     \item hash symbols are replaced with \code{_no};
#'     \item any other special character is replaced with an underscore
#'     character.
#'   }
#'
#' @section Description of the NIST format:
#'   The summary was prepared using the NIST Mass Spectral Search Program manual
#'   v.2.4 (2020).
#'   \itemize{
#'     \item An msp-file can contain as many spectra as wanted.
#'     \item Each spectrum must start with the "Name" field. There must be
#'     something in this field.
#'     \item The "Num Peaks" field is also required. It must contain the number
#'     of mass/intensity pairs.
#'     \item Some optional fields (e.g. "Comments", "Formula", "MW") can be
#'     between the "Name" and "Num Peaks" fields.
#'     \item When a spectrum is exported from the NIST library it also
#'     contains the "NIST#" and "DB#" fields. The "NIST#" field is on the same
#'     line as the "CAS#" field and separated by a semicolon.
#'     \item Each field should be on a separate line (the "NIST#" field is an
#'     exception from this rule)
#'     \item The mass/intensity list begins on the line following the "Num
#'     Peaks" field. The peaks need not be normalized, and the masses need not
#'     be ordered. The exact spacing and delimiters used for the mass/intensity
#'     pairs are unimportant. The following characters are accepted as
#'     delimiters: '\code{space}', '\code{tab}', '\code{,}', '\code{;}',
#'     '\code{:}'. Parentheses, square brackets and curly braces  ('\code{(}',
#'     '\code{(}', '\code{[}', '\code{]}', '\code{\{}', and '\code{\}}') are
#'     also allowed.
#'     \item The "Name" field can be up to 511 characters.
#'     \item The "Comments" field can be up to 1023 characters.
#'     \item The "Formula" field can be up to 23 characters.
#'     \item The "Synon" field may be repeated.
#'   }
#'
#' @examples
#' # Reading the 'alkanes.msp' file
#' msp_file <- system.file("extdata", "alkanes.msp", package = "mssearchr")
#' msp_objs <- ReadMsp(msp_file)
#'
#' # Plotting the first mass spectrum from the 'msp_objs' list
#' par_old <- par(yaxs = "i")
#' plot(msp_objs[[1]]$mz, msp_objs[[1]]$intst,
#'      ylim = c(0, 1000), main = msp_objs[[1]]$name,
#'      type = "h", xlab = "m/z", ylab = "Intensity", bty = "l")
#' par(par_old)
#'
#' @export
#'
#==============================================================================#
ReadMsp <- function(input_file) {

  #--[ Input check ]------------------------------------------------------------

  # 'input_file'
  if (!is.character(input_file) || length(input_file) != 1) {
    stop("'input_file' must be a string.")
  }
  if (!file.exists(input_file)) {
    stop("The file ", "'", basename(input_file), "'", " does not exist.")
  }



  #--[ Reading and pre-processing ]---------------------------------------------

  all_lines <- readLines(input_file)
  first_line_idxs <- grep("^name: ", all_lines, ignore.case = TRUE)
  if (length(first_line_idxs) == 0L) {
    stop("The msp-file is invalid (the 'Name' field is absent).")
  }
  last_line_idxs <- c(first_line_idxs[-1] - 1L, length(all_lines))



  #--[ Parsing ]----------------------------------------------------------------

  msp_objs <- lapply(seq_along(first_line_idxs), function(ms_no) {
    idxs <- seq(first_line_idxs[[ms_no]], last_line_idxs[[ms_no]])
    out <- .ParseSingleSpectrum(all_lines[idxs], first_line_idxs[[ms_no]])
    return(out)
  })



  #--[ Output ]----------------------------------------------------------------

  return(invisible(msp_objs))
}





#==============================================================================#
#' Write mass spectra in an msp-file (NIST format)
#'
#' @description
#'   Write mass spectra in an msp-file (NIST format).
#'
#' @param msp_objs
#'   A list of nested lists. Each nested list is a mass spectrum. Each nested
#'   list must contain at least three elements: (1) \code{name} (a string) -
#'   compound name (or short description); (2) \code{mz} (a numeric/integer
#'   vector) - m/z values of mass spectral peaks; (3) \code{intst} (a
#'   numeric/integer vector) - intensities of mass spectral peaks.
#' @param output_file
#'   A string. The name of a file.
#' @param fields
#'   A character vector. Names of elements in an R list (not the original field
#'   names from an msp-file) to be exported. For example, if only CAS number is
#'   needed to be exported, the \code{'cas_no'} (not \code{'cas#'}) should be
#'   passed. If \code{NULL}, all fields are exported. The output file always
#'   contains the 'Name' field, the 'Num Peaks' field, and the mass/intensity
#'   list.
#'
#' @details
#'   Names of all fields are exported in lower case. It does not cause any
#'   problem in the case of the MS Search (NIST) software (however correct
#'   operation with other software products has not been tested). Only in a few
#'   cases hash symbols and spaces are restored:
#'   \itemize{
#'     \item the \code{cas_no} element is exported as the 'cas#' field;
#'     \item the \code{nist_no} element is exported as the 'nist#' field;
#'     \item the \code{num_peaks} element is exported as the 'num peaks' field.
#'   }
#'
#' @return
#'   \code{NULL} is returned.
#'
#' @examples
#'
#' \dontshow{
#' .old_wd <- setwd(tempdir())
#' }
#' # Exporting mass spectra
#' # Only 'Name', 'SMILES', 'Formula', and 'Num Peaks' fields are exported.
#' WriteMsp(massbank_alkanes[1:3], "test.msp", fields = c("smiles", "formula"))
#' \dontshow{
#' setwd(.old_wd)
#' }
#'
#' @export
#'
#==============================================================================#
WriteMsp <- function(msp_objs,
                     output_file,
                     fields = NULL) {

  #--[ Input check ]------------------------------------------------------------

  # 'msp_objs'
  if (!is.list(msp_objs) || !is.list(msp_objs[[1]]) ||
      is.null(msp_objs[[1]]$mz) || is.null(msp_objs[[1]]$intst) ||
      is.null(msp_objs[[1]]$name)) {
    stop("'msp_objs' is invalid.")
  }

  # 'output_file'
  if (!is.character(output_file) || length(output_file) != 1) {
    stop("'output_file' must be a string.")
  }

  # 'fields'
  if (!is.null(fields)) {
    if (!is.character(fields) || length(fields) == 0L) {
      stop("'temp' must be a character vector.")
    }
  }



  #--[ 'fields' ]---------------------------------------------------------------

  if (!is.null(fields)) {
    mask <- !(fields %in% c("name", "mz", "intst"))
    fields <- unique(fields[mask])
  }



  #--[ Writing data to msp-file ]-----------------------------------------------

  all_lines <- unlist(lapply(msp_objs, function(msp) {

    if (is.null(msp$name)) {
      stop("The '...$name' element is absent.")
    }

    # 'fields_to_export'
    all_fields <- names(msp)
    if (is.null(fields)) {
      mask <- !(all_fields %in% c("name", "mz", "intst"))
    } else {
      mask <- all_fields %in% fields
    }
    fields_to_export <- all_fields[mask]

    # Treating 'irregular' fields ('synon', 'db#', 'cas#', and 'nist#')
    irreg_lines <- list()
    if ("synon" %in% fields_to_export) {
      irreg_lines[[length(irreg_lines) + 1L]] <- paste0("synon: ", msp$synon)
    }
    if ("db_no" %in% fields_to_export) {
      irreg_lines[[length(irreg_lines) + 1L]] <- paste0("db#: ", msp$db_no)
    }
    if (("cas_no" %in% fields_to_export) && ("nist_no" %in% fields_to_export)) {
      irreg_lines[[length(irreg_lines) + 1L]] <-
        paste0("cas#: ", msp$cas_no, ";  ", "nist#: ", msp$nist_no)
    } else if ("cas_no" %in% fields_to_export) {
      irreg_lines[[length(irreg_lines) + 1L]] <- paste0("cas#: ", msp$cas_no)
    } else if ("nist_no" %in% fields_to_export) {
      irreg_lines[[length(irreg_lines) + 1L]] <- paste0("nist#: ", msp$nist_no)
    }

    # Treating 'regular' fields
    reg_fields_to_export <-
      setdiff(fields_to_export, c("synon", "db_no", "cas_no", "nist_no"))
    if (length(reg_fields_to_export) > 0L) {
      reg_lines <- paste0(reg_fields_to_export, ": ", msp[reg_fields_to_export])
    } else {
      reg_lines <- NULL
    }

    # Preparing text lines for exporting
    out <- c(paste0("name: ", msp$name), # 'name'
             unlist(irreg_lines), # 'synon', 'db#', 'cas#', and 'nist#'
             reg_lines, # all regular fields (e.g., 'formula', 'inchikey', 'mw')
             paste0("num peaks: ", length(msp$mz)), # 'num peaks'
             paste0(msp$mz, " ", msp$intst), # mass spectral data
             "")
    return(out)
  }))

  writeLines(all_lines, output_file)



  #--[ Output ]-----------------------------------------------------------------

  return(invisible(NULL))
}




