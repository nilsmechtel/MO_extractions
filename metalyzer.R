library(xlsx)
library(dplyr)

# Define functions --------------------------------------------------------

# read the full Excel sheet ##
openFile <- function(object) {
  full_sheet <- openxlsx::read.xlsx(object@file_path,
                                    sheet = object@sheet,
                                    colNames = FALSE,
                                    skipEmptyRows=FALSE,
                                    skipEmptyCols = FALSE)
  full_sheet <- as.matrix(full_sheet)
  full_sheet[full_sheet == "NA"] <- NA # all entries are strings
  return(full_sheet)
}

## data range ##
getDataRange <- function(object) {
  row_class <- which(object@full_sheet == "Class") %% nrow(object@full_sheet) # row of header "Class"
  col_class <- ceiling(which(object@full_sheet == "Class") / nrow(object@full_sheet)) # column of header "Class"
  
  row_sample_type <- which(object@full_sheet == "Sample Type") %% nrow(object@full_sheet) # row of "Sample Type"
  col_sample_type <- ceiling(which(object@full_sheet == "Sample Type") / nrow(object@full_sheet)) # column of "Sample Type"
  
  rows_data <- which(object@full_sheet[, col_sample_type] == "Sample") # rows
  names(rows_data) <- NULL
  cols_data <- (col_class+1):ncol(object@full_sheet) # columns
  
  data_ranges <- list("class_row"=row_class, "class_col"=col_class,
                      "sample_type_row"=row_sample_type, "sample_type_col"=col_sample_type,
                      "data_rows"=rows_data, "data_cols"=cols_data)
  return(data_ranges)
}

## metabolites with classes ##
readMetabolties <- function(object) {
  metabolites <- object@full_sheet[object@data_ranges[["class_row"]]-1, object@data_ranges[["data_cols"]]] # metabolites are a row above classes
  classes <- object@full_sheet[object@data_ranges[["class_row"]], object@data_ranges[["data_cols"]]]
  names(metabolites) <- classes
  return(metabolites)
}

## raw data ##
readRawData <- function(object) {
  raw_data <- as.data.frame(object@full_sheet[object@data_ranges[["data_rows"]], object@data_ranges[["data_cols"]]])
  colnames(raw_data) <- object@metabolites
  raw_data[] <- sapply(raw_data[], as.numeric)
  return(raw_data)
}

## meta data ##
readMetaData <- function(object) {
  meta_header <- paste(object@full_sheet[object@data_ranges[["sample_type_row"]], 1:object@data_ranges[["class_col"]]])
  meta_data <- as.data.frame(object@full_sheet[object@data_ranges[["data_rows"]], 1:object@data_ranges[["class_col"]]])
  colnames(meta_data) <- meta_header
  meta_data <- meta_data[colSums(is.na(meta_data)) != nrow(meta_data)] # remove full NA columns
  nas <- grep("NA", colnames(meta_data))
  colnames(meta_data)[nas] <- paste0("X", seq_along(nas)) # replace NAs in colnames with X and consecutive numbers
  return(meta_data)
}

## read background color ##
readBGColor <- function(object) {
  wb <- loadWorkbook(object@file_path)
  sheet1 <- getSheets(wb)[[object@sheet]]
  rows <- getRows(sheet1)
  cells <- getCells(rows)
  styles <- sapply(cells, getCellStyle) # get style of each cell
  bg <- sapply(styles, function(style) { # get background color of each cell
    fg  <- style$getFillForegroundXSSFColor()
    rgb <- tryCatch(fg$getRgb(), error = function(e) "")
    rgb <- paste(rgb, collapse = "")
    return(rgb)
  })
  row_index <- as.numeric(unlist(lapply(strsplit(names(bg),split = "\\."),
                                        "[", 1))) # row indexes are at first position
  row_index <- row_index - min(row_index) + 1 # start at 1
  col_index <- as.numeric(unlist(lapply(strsplit(names(bg), split = "\\."),
                                        "[", 2))) # column indexes are at second position
  col_index <- col_index - min(col_index) + 1 # start at 1
  mat_BG <- matrix("", ncol = max(col_index), nrow = max(row_index))
  for (i in 1:length(bg)) { # fill background color matrix
    mat_BG[row_index[i], col_index[i]] <- bg[i]
  }
  df_BG <- as.data.frame(mat_BG)[object@data_ranges[["data_rows"]], object@data_ranges[["data_cols"]]]
  colnames(df_BG) <- object@metabolites
  df_BG[is.na(object@raw_data)] <- NA
  return(df_BG)
}

## assign back ground color to quantification status ##
assignQuantStatus <- function(df_BG) {
  df_BG[df_BG == "6a5acd"] <- "LOD"
  df_BG[df_BG == "87ceeb"] <- "LOQ"
  df_BG[df_BG == "00cd66"] <- "Valid"
  df_BG[df_BG == "ffff33"] <- "Out of calibration range"
  return(df_BG)
}

## filter matabolism indicators ##
filterClasses_f <- function(object, class_name="Metabolism Indicators") {
  if (class_name %in% names(object@metabolites)) {
    not_indicators <- which(names(object@metabolites) != class_name)
    metabolites <- object@metabolites[not_indicators]
    if (nrow(object@raw_data) > 0) {
      object@raw_data <- object@raw_data[, metabolites]
    }
    if (nrow(object@quant_status) > 0) {
      object@quant_status <- object@quant_status[, metabolites]
    }
    object@metabolites <- metabolites
  } else {
    print("-------------------------------------")
    print(paste("No", class_name, "to filter! Returning original object."))
  }
  return(object)
}


# Define metalyzer-class --------------------------------------------------

setClass("metalyzer",
         slots=list(
           ## input ##
           file_path="character",
           sheet="numeric",
           ## output ##
           metabolites="character",
           raw_data="data.frame",
           quant_status="data.frame",
           meta_data="data.frame",
           ## intern ##
           full_sheet="matrix",
           data_ranges="list"
         )
)

## open file and read data ##
setGeneric("readData", function(object) standardGeneric("readData"))
setMethod("readData",
          "metalyzer",
          function(object) {
            object@full_sheet <- openFile(object)
            object@data_ranges <- getDataRange(object)
            object@metabolites <- readMetabolties(object)
            object@raw_data <- readRawData(object)
            object@meta_data <- readMetaData(object)
            return(object)
          }
)

## read quantification status ##
setGeneric("readQuantStatus", function(object) standardGeneric("readQuantStatus"))
setMethod("readQuantStatus",
          "metalyzer",
          function(object) {
            df_BG <- readBGColor(object)
            object@quant_status <- assignQuantStatus(df_BG)
            return(object)
          }
)

## filter metabolism indicators ##
setGeneric("filterClasses", function(object, class_name="Metabolism Indicators") standardGeneric("filterClasses"))
setMethod("filterClasses",
          "metalyzer",
          function(object, class_name) {
            return(filterClasses_f(object, class_name))
          }
)