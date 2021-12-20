library(dplyr)
library(ggplot2)

load("data/plotting.RDATA")

data <- plotting_data %>%
  group_by(Tissue, Method, Status) %>%
  summarise(Number = n()) %>%
  group_by(Tissue, Method) %>%
  mutate(Percentage = Number / sum(Number) * 100,
         Status = factor(Status, levels = c("Valid", "LOQ", "LOD")))

ggplot(data = data,
       aes(x = Method,
           y = Percentage,
           fill = Status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#00cd66", "#87ceeb", "#6a5acd")) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 8)) +
  xlab("") +
  ylab("%") +
  facet_grid(~Tissue)
