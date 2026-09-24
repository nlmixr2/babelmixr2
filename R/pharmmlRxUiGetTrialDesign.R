#' Classify each model covariate as continuous or categorical
#'
#' rxode2 carries no covariate type information at the UI level, so this needs
#' the dataset.  A character or factor column is categorical; anything else is
#' continuous, which is what a numeric column means.
#'
#' `bblDatToNonmem()` converts a character covariate to integer codes (`"F"`,
#' `"M"` become `1`, `2`), so the level labels survive only in the *original*
#' data.  The codes are derived here the same way -- `factor()` ordering -- and
#' then checked against the converted column, so a mismatch is an error rather
#' than a silently wrong mapping.
#'
#' @param ui rxode2 UI
#'
#' @param data The original (pre-conversion) dataset, or `NULL`.  With `NULL`
#'   every covariate is reported continuous, which is all that can be known
#'   from the model alone.
#'
#' @param nmData The `bblDatToNonmem()` conversion of `data`, used to check the
#'   derived category codes, or `NULL` to skip the check.
#'
#' @return named list, one entry per covariate, each with `type` and -- for a
#'   categorical covariate -- `levels` and `codes`
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.pharmmlCovariateInfo <- function(ui, data = NULL, nmData = NULL) {
  .covs <- ui$allCovs
  .ret <- list()
  for (.c in .covs) {
    .col <- .pharmmlDataColumn(data, .c)
    if (is.null(.col) || !(is.character(.col) || is.factor(.col))) {
      .ret[[.c]] <- list(type = "continuous")
      next
    }
    .lev <- levels(factor(.col))
    .codes <- seq_along(.lev)
    .pharmmlCheckCategoryCodes(.c, .col, .lev, nmData)
    .ret[[.c]] <- list(type = "categorical", levels = .lev, codes = .codes)
  }
  .ret
}

#' Verify the derived category codes against the converted dataset
#'
#' The level -> code mapping is *derived* here (`factor()` ordering) rather
#' than read out of `bblDatToNonmem()`, which only returns the codes.  If that
#' conversion ever coded categories differently -- by first appearance, say --
#' the derived mapping would be silently wrong, and a PharmML consumer would
#' read the covariate with its levels transposed.  So it is checked row by row
#' against the converted column and an error is raised on any mismatch.
#'
#' @param name covariate name
#' @param orig the original (character or factor) column
#' @param levels the derived levels
#' @param nmData the converted NONMEM dataset, or `NULL` to skip the check
#' @return Nothing, called for the side effect of erroring
#' @noRd
.pharmmlCheckCategoryCodes <- function(name, orig, levels, nmData) {
  if (is.null(nmData)) {
    return(invisible())
  }
  .got <- .pharmmlDataColumn(nmData, name)
  if (is.null(.got)) {
    return(invisible())
  }
  .expect <- match(as.character(orig), levels)
  if (
    length(.got) != length(.expect) ||
      !isTRUE(all.equal(as.numeric(.got), as.numeric(.expect)))
  ) {
    stop(
      "cannot map the categories of covariate '",
      name,
      "' to PharmML: the codes `bblDatToNonmem()` produced do not match the ",
      "levels derived from the data (expected ",
      paste(levels, "=", seq_along(levels), collapse = ", "),
      ")",
      call. = FALSE
    )
  }
  invisible()
}

#' Look a column up in a dataset, case-insensitively
#'
#' The converted NONMEM dataset upper-cases covariate names, so a lookup by the
#' model's own spelling has to ignore case.
#'
#' @param data data frame or NULL
#' @param name column name
#' @return the column, or NULL when absent
#' @noRd
.pharmmlDataColumn <- function(data, name) {
  if (is.null(data)) {
    return(NULL)
  }
  .w <- which(tolower(names(data)) == tolower(name))
  if (length(.w) != 1L) {
    return(NULL)
  }
  data[[.w]]
}

#' The NONMEM-format dataset for a model
#'
#' @param ui rxode2 UI
#' @param data original dataset
#' @return data frame, with the internal row-number column dropped
#' @noRd
.pharmmlNonmemData <- function(ui, data) {
  .ret <- bblDatToNonmem(ui, data)
  .ret[, names(.ret) != "nlmixrRowNums", drop = FALSE]
}

# NONMEM column name -> PharmML ds:Column/@columnType.  The names are the
# standard slots `getStandardColNames()` reports.
.pharmmlColumnType <- c(
  id = "id",
  time = "idv",
  amt = "dose",
  rate = "rate",
  dur = "duration",
  evid = "evid",
  cmt = "cmt",
  ss = "ss",
  ii = "ii",
  addl = "addl",
  dv = "dv",
  mdv = "mdv",
  dvid = "dvid",
  cens = "censoring",
  limit = "limit"
)

#' The PharmML valueType for a data column
#'
#' @param x the column
#' @return "int" or "real"
#' @noRd
.pharmmlValueType <- function(x) {
  if (is.integer(x)) {
    return("int")
  }
  if (is.numeric(x) && all(is.na(x) | x == trunc(x))) {
    return("int")
  }
  "real"
}

#' Map a NONMEM dataset column to its PharmML columnType
#'
#' @param col column name as it appears in the converted dataset
#' @param std the `getStandardColNames()` result for that dataset
#' @param covs model covariate names
#' @return a PharmML columnType
#' @noRd
.pharmmlColumnTypeOf <- function(col, std, covs) {
  .slot <- names(std)[which(!is.na(std) & std == col)]
  if (length(.slot) == 1L && .slot %in% names(.pharmmlColumnType)) {
    return(.pharmmlColumnType[[.slot]])
  }
  if (tolower(col) %in% tolower(covs)) {
    return("covariate")
  }
  "undefined"
}

#' The model symbol a dataset column maps to
#'
#' @param col column name in the converted dataset
#' @param std the `getStandardColNames()` result
#' @param ui rxode2 UI
#' @return an emitted `ct:SymbRef`, or `NULL` when the column maps to nothing
#' @noRd
.pharmmlColumnSymbRef <- function(col, std, ui) {
  .slot <- names(std)[which(!is.na(std) & std == col)]
  if (length(.slot) == 1L) {
    if (.slot == "time") {
      return(.pmlNode("ct:SymbRef", attrs = c(symbIdRef = "t")))
    }
    if (.slot == "dv") {
      # single-endpoint models map DV straight onto the observation
      .predDf <- ui$predDf
      if (nrow(.predDf) == 1L) {
        return(.pmlNode(
          "ct:SymbRef",
          attrs = c(
            blkIdRef = "om1",
            symbIdRef = paste0(paste(.predDf$cond[1]), "_obs")
          )
        ))
      }
      return(NULL)
    }
    # id/evid/cmt/amt and friends are structural data columns; PharmML knows
    # what they mean from columnType alone, so they need no symbol mapping
    return(NULL)
  }
  .cov <- ui$allCovs[tolower(ui$allCovs) == tolower(col)]
  if (length(.cov) == 1L) {
    return(.pmlNode(
      "ct:SymbRef",
      attrs = c(blkIdRef = .pmlBlk[["covariate"]], symbIdRef = .cov)
    ))
  }
  NULL
}

#' PharmML TrialDesign
#'
#' The dataset is exported in NONMEM format, which PharmML maps natively
#' through `ds:Column/@columnType` -- the same route the upstream
#' `example1_NONMEM.xml` takes.  Reusing `bblDatToNonmem()` means the exported
#' data is byte-identical to what `est="nonmem"` would write.
#'
#' @param ui rxode2 UI
#'
#' @param data The original dataset
#'
#' @param dataFile Name the emitted document should use to refer to the CSV
#'
#' @param indent Indent depth
#'
#' @return character(1) holding the `TrialDesign` element
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.pharmmlTrialDesign <- function(ui, data, dataFile = "data.csv", indent = 0L) {
  .nm <- .pharmmlNonmemData(ui, data)
  .std <- getStandardColNames(.nm)
  .covInfo <- .pharmmlCovariateInfo(ui, data, .nm)

  .mappings <- character(0)
  .columns <- character(0)

  for (.i in seq_along(names(.nm))) {
    .col <- names(.nm)[.i]
    .sym <- .pharmmlColumnSymbRef(.col, .std, ui)
    if (!is.null(.sym)) {
      .children <- c(
        .pmlNode("ds:ColumnRef", attrs = c(columnIdRef = .col)),
        .sym
      )
      .cov <- ui$allCovs[tolower(ui$allCovs) == tolower(.col)]
      if (
        length(.cov) == 1L && identical(.covInfo[[.cov]]$type, "categorical")
      ) {
        .info <- .covInfo[[.cov]]
        .maps <- vapply(
          seq_along(.info$levels),
          function(.j) {
            .pmlNode(
              "ds:Map",
              attrs = c(
                dataSymbol = as.character(.info$codes[.j]),
                modelSymbol = .info$levels[.j]
              )
            )
          },
          character(1),
          USE.NAMES = FALSE
        )
        .children <- c(
          .children,
          .pmlNode("ds:CategoryMapping", children = .maps)
        )
      }
      .mappings <- c(
        .mappings,
        .pmlNode("design:ColumnMapping", children = .children)
      )
    }
    .columns <- c(
      .columns,
      .pmlNode(
        "ds:Column",
        attrs = c(
          columnId = .col,
          columnType = .pharmmlColumnTypeOf(.col, .std, ui$allCovs),
          valueType = .pharmmlValueType(.nm[[.i]]),
          columnNum = as.character(.i)
        )
      )
    )
  }

  .dataSet <- .pmlNode(
    "ds:DataSet",
    children = c(
      .pmlNode("ds:Definition", children = .columns),
      .pmlNode(
        "ds:ExternalFile",
        attrs = c(oid = "dataOid"),
        children = c(
          .pmlText("ds:path", dataFile),
          .pmlText("ds:format", "CSV"),
          .pmlText("ds:delimiter", "COMMA")
        )
      )
    )
  )

  .pmlNode(
    "design:TrialDesign",
    children = .pmlNode(
      "design:ExternalDataSet",
      attrs = c(toolName = "NONMEM", oid = "nmOid"),
      children = c(.mappings, .dataSet)
    ),
    indent = indent
  )
}

#' Write the NONMEM-format dataset a PharmML document refers to
#'
#' @param ui rxode2 UI
#' @param data original dataset
#' @param file path to write to
#' @return `file`, invisibly
#' @noRd
.pharmmlWriteData <- function(ui, data, file) {
  utils::write.csv(
    .pharmmlNonmemData(ui, data),
    file,
    row.names = FALSE,
    quote = FALSE,
    na = "."
  )
  invisible(file)
}
