library(dplyr)
library(ggplot2)
library(stringr)
library(venn)
library(xlsx)

load("data/plotting.RDATA")

valid_data <- filter(plotting_data, valid_replicates)

# Bar plot ----------------------------------------------------------------

data <- valid_data %>%
  group_by(Tissue, Method, Class) %>%
  summarise(Number = n())

ggplot(data = data,
       aes(x = Method,
           y = Number,
           fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual("Classes",
                    breaks = breaks,
                    values = values,
                    drop = FALSE,
                    guide = guide_legend(order=2, ncol = 2)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 8)) +
  xlab("") +
  facet_grid(~Tissue)


# Venn diagrams -----------------------------------------------------------

intersection <- function(set_lists, indexes) {
  incl_tissue <- names(set_lists)[indexes]
  subset_set_lists <- set_lists[incl_tissue]
  not_incl_tissue <- names(set_lists)[-indexes]
  not_incl_metab <- unique(unlist(set_lists[not_incl_tissue]))
  
  y <- not_incl_metab
  
  if (length(subset_set_lists) == 1) {
    x <- subset_set_lists[[1]]
    out_vec <- setdiff(x, y)
  } else if (length(subset_set_lists) %in% c(2,3)) {
    x <- intersect(subset_set_lists[[1]], subset_set_lists[[2]])
    if (length(subset_set_lists) == 3) {
      x <- intersect(x, subset_set_lists[[3]])
    }
    out_vec <- setdiff(x, y)
  } else {
    out_vec <- intersect(intersect(subset_set_lists[[1]], subset_set_lists[[2]]),
                         intersect(subset_set_lists[[3]], subset_set_lists[[4]]))
  }
  return(out_vec)
}

write_venn_metabos <- function(method, set_lists, file, append = TRUE) {
  combis <- expand.grid(c(1,2,3,4), c(1,2,3,4), c(1,2,3,4), c(1,2,3,4))
  unique_combis <- unique(lapply(1:nrow(combis), function(i) {
    sort(unique(unlist(combis[i,])))
  }))
  venn_lists <- lapply(unique_combis, function(indexes) {
    intersection(set_lists, indexes)
  })
  names(venn_lists) <- sapply(unique_combis, function(indexes) {
    paste(names(set_lists)[indexes], collapse = "_")
  })
  
  max_len <- max(sapply(venn_lists, length))
  venn_lists <- lapply(venn_lists, function(vec) {
    out_vec <- c(vec, rep(NA, max_len-length(vec)))
    return(out_vec)
  })
  out_df <- data.frame(bind_cols(venn_lists))
  write.xlsx(out_df, file, showNA = FALSE, row.names = FALSE,
             sheetName = str_replace_all(method, "/", "-"), append = append)
}

for (method in levels(valid_data$Method)) {
  set_lists <- lapply(levels(valid_data$Tissue), function(tissue) {
    valid_data %>%
      filter(Tissue == tissue,
             Method == method) %>%
      ungroup() %>%
      select(Metabolite) %>%
      unlist() %>%
      as.character()
  })
  names(set_lists) <- levels(valid_data$Tissue)
  
  venn(set_lists, ilab=TRUE, zcolor = col_vec_tissues,
       ilcs = 2, sncs = 1.5, box = FALSE)
  
  append <- ifelse(method == levels(valid_data$Method)[1], FALSE, TRUE)
  write_venn_metabos(method, set_lists, "data/metabolite_coverage.xlsx",
  append = append)
}
