library(MetAlyzer)
library(stringr)

# Colors ------------------------------------------------------------------

classes <- c("Acylcarnitines", "Alkaloids", "Amine Oxides", "Aminoacids",
             "Aminoacids Related", "Bile Acids", "Biogenic Amines", "Carboxylic Acids",
             "Ceramides", "Cholesterol Esters", "Cresols", "Diacylglycerols",
             "Dihydroceramides", "Fatty Acids", "Glycerophospholipids",
             "Glycosylceramides", "Hormones", "Indoles Derivatives",
             "Nucleobases Related", "Sphingolipids", "Sugars", "Triacylglycerols",
             "Vitamins & Cofactors")

col_vec_classes <- c("#a6cee3", "#1f78b4", "#ffff99", "#b2df8a", "#33a02c", "#fb9a99",
                        "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#8dd3c7",
                        "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69",
                        "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
names(col_vec_classes) <- classes

polarity <- c("FIA", "LC", "LC", "LC", "LC", "LC", "LC", "LC", "FIA", "FIA",
              "LC", "FIA","FIA", "LC", "FIA",  "FIA", "LC", "LC", "LC", "FIA",
              "FIA", "FIA", "LC")
names(polarity) <- classes
polarity <- sort(polarity, decreasing = TRUE)

classes <- names(polarity)
lc_classes <- names(polarity[which(polarity == "LC")])
fia_classes <- names(polarity[which(polarity == "FIA")])

lc_col_vec <- col_vec_classes[lc_classes]
fia_col_vec <- col_vec_classes[fia_classes]

breaks <- c("LC:", lc_classes,
            "FIA:", fia_classes, "", " ", "  ")
values <- c("white", lc_col_vec,
            "white", fia_col_vec, rep("white", 3))
names(values) <- NULL

col_vec_methods <- c("#B9DBF4", "#8785BA", "#B56478", "#F3C2A2", "#17506C",
                        "#72BEB7", "#B8E2DE")
methods <- c(
  "1" = "100 IPA",
  "2" = "IPA/ACN/H2O (3/2/2)",
  "4" = "MeOH/ACN/H2O+FA (2/2/1)",
  "6" = "MeOH/CHCl3/H2O", #
  "3" = "75 EtOH/MTBE", #
  "5" = "MeOH/MTBE/MeOH", #
  "7" = "2x MeOH/MTBE/MeOH" #
)
names(col_vec_methods) <- methods

tissues <- c("Z" = "Zebrafish Liver",
             "D" = "Drosophila",
             "ML" = "Mouse Liver",
             "MN" = "Mouse Kidney")
col_vec_tissues <- c("#FD6A02", "#3C89D0", "#ACDF87", "#68BB59")
names(col_vec_tissues) <- tissues

# Loading data ------------------------------------------------------------

file_path <- "/Users/nilsm/OneDrive - bwedu/Uni/HiWi COS/4_SFB_Ectraction_Tests/MO_extractions/data/Extraction_test_model_organisms.xlsx"
obj_organisms <- MetAlyzerDataset(file_path = file_path)

obj_organisms <- filterMetabolites(obj_organisms, class_name = "Metabolism Indicators")

# Adjust meta data --------------------------------------------------------

# extract methods from meta data
method_vec <- sapply(metaData(obj_organisms)$X1, function(m) {
  factor(methods[m], levels = methods)
})

# extract tissue from meta data
pattern <- paste0("(", paste(names(tissues), collapse = "|"), ")")
tissue_vec <- sapply(metaData(obj_organisms)$`Sample Identification`, function(id) {
  x <- str_extract_all(id, pattern)[[1]][1]
  if (length(x) > 0) {
    for (i in 1:4) {
      x <- stringr::str_replace(x, names(tissues)[i], tissues[i])
    }
  }
  return(factor(x, levels = tissues))
})

obj_organisms <- updateMetaData(obj_organisms, Method, method_vec)
obj_organisms <- updateMetaData(obj_organisms, Tissue, tissue_vec)

obj_organisms <- filterMetaData(obj_organisms, Tissue, keep = tissues)

## plotting data
obj_organisms <- createPlottingData(obj_organisms, Tissue, Method)
obj_organisms <- imputePlottingData(obj_organisms, Tissue, Metabolite)
obj_organisms <- transformPlottingData(obj_organisms)
obj_organisms <- performANOVA(obj_organisms, Method)

## update class levels
plotting_data <- plottingData(obj_organisms)
plotting_data$Class <- factor(plotting_data$Class, levels = breaks)
obj_organisms <- setPlottingData(obj_organisms, plotting_data)

## save R Data ##
save(breaks, values, methods, col_vec_methods, col_vec_tissues, plotting_data,
     file = "data/plotting.RDATA")
