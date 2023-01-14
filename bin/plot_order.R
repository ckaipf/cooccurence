#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

read_delim(file = args[1], delim = ",") %>%
  mutate(across(everything(), .fns = ~ str_split(., "_", simplify = T)[, 1])) %>%
  group_by(across(everything())) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(f  = n / sum(n)) %>%
  arrange(n) %>%
  rowid_to_column() %>%
  pivot_longer(cols = !starts_with(c("n", "rowid", "f")), names_to = "x") %>%
  filter(!is.na(value)) %>% 
  mutate(max_y = max(rowid),
         max_x = max(as.numeric(x))) %>%
  { max_y <<- unlist(distinct(select(., max_y))); .} %>%
  ggplot(aes(x = as.numeric(x), y = rowid, group = rowid, label = value)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_x_continuous(expand = c(0.5, 0)) +
  annotate("rect", xmin=-0.7, xmax=-0.3, ymin=0, ymax=max_y+1, alpha=0.8, fill="lightblue", col = "darkblue") +
  geom_point(aes(x = -0.5, col = f), size = 5) +
  scale_color_viridis_c(name = "Frequency") + 
  geom_text(hjust=.5, vjust=-1) +
  theme(
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.background = element_rect(fill="lightblue",
                                     size=0.5, 
                                     linetype="solid", 
                                     colour ="darkblue"),
    legend.position = "left"
  ) +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  ) -> p 
ggsave(plot = p, filename = paste(args[2], "_freq_of_orders.png", sep = ""), device = "png", width = 10)
