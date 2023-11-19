#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

dat <- read_csv(file = args[1], col_types = "f") 
dat %>%
  distinct(across(everything())) %>%
  mutate(across(
    !starts_with(c("tag")),
    .fns = function(x)
      str_split(x, "_", simplify = T)[, 1]
  )) %>%
  group_by(across(everything())) %>%
  summarise(n = n()) %>%
  unite(remove = F,
        !starts_with(c("tag", "n")),
        col = "combination",
        sep = "_") %>%
  mutate(combination = str_remove_all(combination, "NA_|_NA")) %>%
  pivot_longer(values_to = "set",!starts_with(c("tag", "n", "combination"))) %>%
  select(-name) %>%
  filter(!is.na(set)) %>%
  group_by(set, tag) %>%
  mutate(total = sum(n), r = n / total) %>%
  mutate(combination = if_else(set == combination, ".", combination)) %>%
  #
  ggplot(aes(
    x = tag,
    y = r,
    fill = combination,
    group = combination
  )) +
  geom_bar(
    width = 1,
    stat = "identity",
    colour = "black",
    linetype = "dashed",
    alpha = .8,
    size = 0.1
  ) +
  geom_label(
    aes(
      y = r,
      label = if_else(
        r > 0.1,
        if_else(
          str_ends(combination, set),
          str_remove_all(combination, paste0("_", set)),
          str_remove_all(combination, paste0(set, "_")),
        ) %>%
          str_replace_all("_", "\n"),
        as.character(NA)
      ),
      group = combination
    ),
    fontface = "bold",
    color = "black",
    fill = "white",
    lineheight = .9,
    alpha = .5,
    position = position_stack(vjust = .5),
    angle = 0,
    size = 3.5
  ) +
  facet_grid(set ~ ., switch = "y") +
  scale_fill_viridis_d(option = "D") +
  scale_x_discrete(name = NULL, position = "top") +
  scale_y_continuous(
    name = NULL,
    labels = NULL,
    breaks = NULL,
    limits = c(0, 1)
  ) +
  guides(fill = guide_legend(
    title = "Combination",
    nrow = 3,
    byrow = T
  )) +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(colour = "navy", fill = "lightskyblue"),
    legend.key.size = unit(2, "mm"),
    legend.text = element_text(face = "bold"),
    panel.background = element_blank(),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 10),
    axis.text.x.top = element_text(angle = 0, face = "bold", size = 8),
    axis.ticks = element_blank(),
    axis.line.y.left = element_blank(),
    strip.background = element_blank(),
  ) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) -> p
ggsave(
  plot = p,
  filename = "collected_barplot.png",
  device = "png",
  width = 1.8 * nlevels(dat$tag)
)
