library(dplyr)
library(ggplot2)

load("data/plotting.RDATA")

max_metabolites <- plotting_data %>%
  filter(!is.na(ANOVA_group)) %>%
  group_by(Tissue, Metabolite) %>%
  summarise() %>%
  group_by(Tissue) %>%
  summarise(max = n())

group_A_data <- plotting_data %>%
  filter(sapply(ANOVA_group, function(x) grepl("A", x))) %>%
  group_by(Tissue, Method, Class, Metabolite, ANOVA_group) %>%
  summarise()

data1 <- group_A_data %>%
  group_by(Tissue, Method, Class) %>%
  summarise(Number = n())

ggplot(data1,
       aes(x=Method,
           y=Number,
           fill=Class)) +
  geom_hline(data=max_metabolites,
             aes(yintercept=max), linetype = 2) +
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
  ylab("Number") + 
  xlab("") +
  facet_grid(~Tissue)


data2 <- group_A_data %>%
  group_by(Tissue, Metabolite) %>%
  mutate(n_opt = sum(sapply(ANOVA_group, function(x) grepl("A", x)))) %>%
  filter(n_opt == 1,
         sapply(ANOVA_group, function(x) grepl("A", x))) %>%
  group_by(Tissue, Method, Class) %>%
  summarise(Number = n())

ggplot(data2,
       aes(x=Method,
           y=Number,
           fill=Class)) +
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
  ylab("Number") + 
  xlab("") +
  facet_grid(~Tissue)
