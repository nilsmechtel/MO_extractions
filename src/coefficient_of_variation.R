library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggradar)
library(ggpubr)
library(gridExtra)

load("data/plotting.RDATA")

valid_data <- plotting_data %>%
  filter(valid_replicates) %>%
  group_by(Tissue, Class, Metabolite, Method, CV, CV_thresh) %>%
  summarise()

# Bar plots ---------------------------------------------------------------

data1 <- valid_data %>%
  group_by(Tissue, Method, CV_thresh) %>%
  summarise(Number = n())

ggplot(data = data1,
       aes(x = Method,
           y = Number,
           fill = CV_thresh)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual("CV",
                    values = c("#5AAA46", "#A0AA46", "#D37538", "#C84D4C"),
                    labels = c("0-10 %", "11-20 %", "21-30 %", ">30%")) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 8)) +
  xlab("") +
  facet_grid(~Tissue)

## Median CV + MAD
data2 <- valid_data %>%
  group_by(Tissue, Method) %>%
  summarise(Median_CV = median(CV, na.rm = TRUE), # median CV
            MAD_CV = mad(CV, na.rm = TRUE)) # MAD CV

ggplot(data2, aes(Method, Median_CV)) +
  geom_bar(aes(fill = Tissue), 
           position = "dodge", 
           stat = "identity") +
  geom_errorbar(aes(ymin = Median_CV, 
                    ymax = Median_CV+MAD_CV, 
                    group = Tissue),
                width = .2,
                position = position_dodge(.9), 
                col = "black") +
  scale_fill_manual(values = col_vec_tissues) +
  xlab("") +
  ylab("Median CV + MAD") +
  theme_light() +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Radar plots -------------------------------------------------------------
gg_radarchart <- function(df, col_vector_methods, title, l.pos="bottom") {
  ggradar(rownames_to_column(df[3:nrow(df),], "group"),
          values.radar = c("0 %", paste(df["Max", 1]/2, "%"), paste(df["Max", 1], "%")),
          grid.min = 0, grid.mid = df["Max", 1]/2, grid.max = df["Max", 1],
          axis.labels = c(
            "metabolites\nCV 0-10%",
            "metabolites\nCV 11-20%",
            "metabolites\nCV 21-30%",
            "metabolites\nCV > 30%"
          ),
          group.line.width = 1,
          group.point.size = 3,
          group.colours = col_vector_methods,
          background.circle.colour = "white",
          gridline.mid.colour = "grey",
          plot.title = title,
          legend.position = l.pos) + 
    theme(plot.title = element_text(hjust = 0.5))
}

p_list <- lapply(unique(plotting_data$Tissue), function(tissue) {
  # percentage of metabolites below each threshold
  tmp_df <- valid_data  %>%
    filter(Tissue == tissue) %>%
    group_by(Method, CV_thresh) %>%
    summarise(Percentage = n() / 630 * 100) %>%
    spread(CV_thresh, Percentage) %>%
    data.frame() %>%
    column_to_rownames("Method")
  
  # first 2 rows determine the limits for each axis
  max_min <- data.frame(matrix(c(rep(50, 4), rep(0, 4)), ncol = 4, byrow = TRUE))
  colnames(max_min) <- colnames(tmp_df)
  rownames(max_min) <- c("Max", "Min")
  df <- rbind(max_min, tmp_df)
  
  gg_radarchart(df, col_vec_methods, tissue)
})
legend <- get_legend(p_list[[4]]) # save legend
p_list <- lapply(p_list, function(p) {
  p + theme(legend.position="none") # remove legends
})
p_list[[5]] <- legend

grid.arrange(
  grobs = p_list,
  widths = c(1,3,3,1),
  heights = c(4,4,1),
  layout_matrix = rbind(c(1,1,2,2),
                        c(3,3,4,4),
                        c(NA,5,5,NA))
)
