library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(agricolae)

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

source("metalyzer.R")

obj_organisms <- new("metalyzer",
                     file_path = "data/Extraction_test_model_organisms.xlsx",
                     sheet = 1)
obj_organisms <- readData(obj_organisms)
obj_organisms <- readQuantStatus(obj_organisms)
obj_organisms <- filterClasses(obj_organisms)


# Prepare data for plotting -----------------------------------------------

## extract tissue from meta data
pattern <- paste0("(", paste(names(tissues), collapse = "|"), ")")
tissue_vec <- sapply(obj_organisms@meta_data$`Sample Identification`, function(id) {
  x <- str_extract_all(id, pattern)[[1]][1]
  if (length(x) > 0) {
    for (i in 1:4) {
      x <- str_replace(x, names(tissues)[i], tissues[i])
    }
  }
  return(x)
})

## extract method from meta data
method_vec <- obj_organisms@meta_data$X1

## add tissue and method as columns
tmp_data <- obj_organisms@raw_data[!is.na(method_vec),]
tmp_data$Tissue <- tissue_vec[!is.na(method_vec)]
tmp_data$Method <- method_vec[!is.na(method_vec)]

tmp_status <- obj_organisms@quant_status[!is.na(method_vec),]
tmp_status$Tissue <- tissue_vec[!is.na(method_vec)]
tmp_status$Method <- method_vec[!is.na(method_vec)]

## reshape data
gathered_data <- gather(tmp_data, key = Metabolite, value = Concentration, -Tissue, -Method)
gathered_status <- gather(tmp_status, key = Metabolite, value = Status, -Tissue, -Method)
gathered_data$Status <- gathered_status$Status

assign_class <- function(metabolite) {
  metabolites <- obj_organisms@metabolites
  classes <- names(metabolites)
  class <- classes[which(metabolites == metabolite)]
  class <- factor(class, levels = breaks)
  return(class)
}

assign_method <- function(number) {
  factor(methods[number], levels = methods)
}
  
plotting_data <- gathered_data %>%
  mutate(Tissue = factor(Tissue, levels = tissues),
         Class = sapply(Metabolite, assign_class),
         Metabolite = factor(Metabolite, levels = obj_organisms@metabolites),
         Method = sapply(Method, assign_method)) %>%
  relocate(Class, Metabolite, .after = Tissue) %>%
  arrange(Tissue, Metabolite, Method)

## calculate coeficcient of variation (CV) ##
set_treshold <- function(x) {
  if (is.na(x)) {
    v <- NA
  } else if (x < 0.1) {
    v <- "max10"
  } else if (x < 0.2) {
    v <- "max20"
  } else if (x < 0.3) {
    v <- "max30"
  } else {
    v <- "more30"
  }
  return(v)
}

# per tissue, metabolite and method
plotting_data <- plotting_data %>%
  group_by(Tissue, Metabolite, Method) %>%
  mutate(CV = sd(Concentration) / mean(Concentration), # CV = SD / Mean
         CV_tresh = sapply(CV, set_treshold),
         .before = Status)

## filter based on quantification status ##
filter_satatus <- function(vec) {
  sum(vec %in% c("Valid", "LOQ")) / length(vec) > 0.5 # TRUE if > 50% valid/LOQ
}

# per tissue, metabolite and method
plotting_data <- plotting_data %>%
  group_by(Tissue, Metabolite, Method) %>%
  mutate(Valid = filter_satatus(Status),
         .after = Status)

## perform zero imputation with 20% of the minimal positive value ##
zero_imputation <- function(vec) {
  non_zero <- vec[vec > 0]
  imp_v <- ifelse(length(non_zero) > 0, min(non_zero) / 5, NA)
  vec[vec == 0] <- imp_v
  return(vec)
}

# per tissue and metabolite
plotting_data <- plotting_data %>%
  group_by(Tissue, Metabolite) %>%
  mutate(imp_Conc = zero_imputation(Concentration))

## log2 transformation
log2_transform <- function(vec) {
  vec[vec > 0 & !is.na(vec)] <- log2(vec[vec > 0 & !is.na(vec)])
  return(vec)
}
plotting_data$log2_Conc <- log2_transform(plotting_data$imp_Conc)

## calculate ANOVA ##
calc_anova <- function(m_vec, v_vec, valid_vec) {
  # if all concentration values equal 0 (log2 -> NA) or no method achieves valid
  # concentrations no ANOVA is calculated (output: NA)
  if (all(is.na(v_vec)) | sum(valid_vec) == 0) {
    group_vec <- as.character(rep(NA, length(m_vec)))
  } else {
    tmp_df <- data.frame(Method = as.character(m_vec),
                         Values = as.numeric(v_vec))
    # ANOVA
    anova <- aov(Values ~ Method, data = tmp_df)
    # Tukey post-hoc; each methods gets assigned to a group
    groups <- HSD.test(anova, 'Method', group=TRUE)$groups %>%
      select(-Values) %>%
      rownames_to_column('Method') %>%
      mutate(Method = factor(Method, levels = levels(m_vec))) %>%
      arrange(Method) %>%
      deframe() %>%
      toupper()
    
    group_vec <- sapply(m_vec, function(m) groups[m])
  }
  return(group_vec)
}

# per tissue and metabolite
plotting_data <- plotting_data %>%
  group_by(Tissue, Metabolite) %>%
  mutate(Group = calc_anova(Method, log2_Conc, Valid))

## save R Data ##
save(breaks, values, methods, col_vec_methods, col_vec_tissues, plotting_data,
     file = "data/plotting.RDATA")
