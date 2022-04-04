#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.files = [
"1.gff",
"2.gff",
"3.gff"
]

params.default_bedtools_parameters = "-s -k 1"
params.default_distance = 150
params.config = "example.config"
params.tag = "example_run_0"
params.min_comb_freq = 0.05

workflow {
  main:
  pairedIntersections(params.files, params.config,)
}

workflow pairedIntersections {
  take: files
        config
  main:

  files_ch = Channel.fromPath(files) | \
          sortGff | \
          map { it -> [it.getSimpleName(), it] }

  config = Channel.fromPath(config) | \
      splitCsv() | \
      map { row -> (row[0] == params.tag) ? row : null } | \
      map { row -> (files.contains(row[1] + ".gff") && files.contains(row[2] + ".gff")) ? row : null}    

  
  parameters = config | \
    map { it -> [ it[1], it[2], it[3] ] }

  distances = config | \
    map { it -> [ it[1], it[2], it[4] ] }

  x  = files_ch.combine(files_ch) | \
    filter(it -> it[0] != it[2]) | \
    map {it -> [it[0], it[2], it[1], it[3]]} |  \
    join(parameters, by: [0,1], remainder: true) |  \
    map { it -> (it[4] == null) ? it[0..3] + params.default_bedtools_parameters : it } | \
    closestBed | \
    join(distances, by: [0,1], remainder: true) |  \
    map { it -> (it[3] == null) ? it[0..2] + params.default_distance : it } | \
    map { it -> ['"'+it[0]+'"', '"'+it[1]+'"', '"'+it[3]+'"', '"'+it[2]+'"'] } | \
    toList | \
    fullJoin 
    
    x | (plotVenn & barPlot) 
}


process sortGff {
  input:
  path(a)
  output:
  path("${a}.sorted")

  """
  sort -k1,1 -k4,4n -k5,5n ${a} > ${a}.sorted
  """
}

process closestBed {
//  publishDir ".", mode: "copy"
  input:
  tuple val(i), val(j), path(a), path(b), val(parameters)
  output:
  tuple val(i), val(j), path("*.closest")

"""
bedtools closest -D a ${parameters} -a ${a} -b ${b} > ${i}${j}.closest
"""
}

process fullJoin {
  publishDir ".", mode: "copy"
  input:
  val closestOutput
  output:
  path "*.RData"
  
  script:
  //  Groovy list to R list string
  input_list = closestOutput
      .collect{f -> f.toString().replace("[", "c(").replace("]", ")")}
      .toString()
      .minus("[")
      .minus("]")

  """
  #!/usr/bin/Rscript
  library(tidyverse)

  na_match = "na"

  x <- list(${input_list}) 
  permutations <- map(x, ~ paste0(chuck(., 1), chuck(., 2))) %>% unlist
  sets <- map(x, ~ .[1]) %>% unlist(use.names = F) %>% unique
  names(x) <- permutations
  names(permutations) <- permutations

  pairs <- x %>%
  map(function(y) {
    read_delim(file = chuck(y, 4), delim = "\t", col_names = F, col_types = "ccccccccccccccccccc") %>%
    mutate("{chuck(y, 1)}" := paste0(chuck(y, 1), X4, X5, X7)) %>%
    mutate("{chuck(y, 2)}" := paste0(chuck(y, 2), X13, X14, X16)) %>%
    mutate(distance = abs(as.numeric(X19))) %>%
    select(!starts_with("X")) %>%
    distinct(across(everything()))
    }
  )

  pairs_filtered <- 
    map2(.x = pairs, .y = x, .f = ~ filter(.x, distance < as.numeric(chuck(.y, 3)))) %>%
    map(~ select(., -distance))

  sets_total <- sets %>%
    map(~ permutations[str_starts(permutations, .)]) %>%
    map(~ first(.)) %>%
    map(~ map(., .f = function(i) chuck(pairs, i))) %>% 
    flatten %>%
    map(~ .[1]) 

  sets_full <- sets %>%
    map(~ permutations[str_starts(permutations, .)]) %>%
    map(~ map(., .f = function(i) chuck(pairs_filtered, i))) %>%
    map(~ reduce(., function(acc, y) full_join(acc, y, na_matches = na_match))) 

  joined_filtered <- sets_full %>%
    reduce(function(acc, y) full_join(acc, y, na_matches = na_match))

  venn_complete <- sets_total %>%
    reduce(.f = function(acc, y) full_join(acc, y, na_matches = na_match), .init = joined_filtered) %>%
    rowid_to_column()
    
  save(venn_complete, sets, file = "${params.tag}.RData")
  """
}

process plotVenn{
  publishDir ".", mode: "copy"
  input:
  path Rdata
  output:
  path "*.png"

  """
  #!/usr/bin/Rscript
  library(tidyverse)
  load("${Rdata}")
  
  names(sets) <- sets
  ids <- sets %>%
    map(~ filter(venn_complete, !is.na(.data[[.]]))) %>%
    map(~ select(., rowid) %>% unlist %>% unique) 

  p <- ggVennDiagram::ggVennDiagram(ids)
  ggsave("${params.tag}_venn.png", plot = p, device = "png")
  save(venn_complete, sets, file = "${params.tag}.RData")
  """
}

process barPlot {
  publishDir ".", mode: "copy"
  input:
  path Rdata
  output:
  path "*.png"

"""
 #!/usr/bin/Rscript
  library(tidyverse)
  load("${Rdata}")

map2(.x = sets, .y = sets, .f = function(x, y) {
  seq(sets) %>%
    map(
      ~ {combn(sets, .) %>% t %>% data.frame() %>%filter(if_any(everything(), .fns = ~ . %in% c(x)))})%>%
    map(.f = ~ pmap(., c)) %>% 
    flatten %>% 
    map(as.vector) %>%
    set_names(map_chr(., .f = function(x) paste(x, collapse =  "_" ) %>% as.vector)) %>% 
    map(~ filter(venn_complete, 
                 if_all(all_of(.), .fns = ~ !is.na(.)) & if_all(sets[!sets %in% .], .fns = is.na)
    )
    ) %>%
    map(~ distinct(., .data[[x]])) %>%
    map(count) %>%
    map(unlist) %>%
    bind_cols() %>%
  pivot_longer(cols = everything(), names_to = "combination") %>%
  bind_cols(., set = y)
}) %>%
  bind_rows %>%
  group_by(set) %>%
  mutate(total = sum(value)) %>%
  ungroup %>%  
  group_by(set) %>% 
  mutate(r = value / total) %>%
  mutate(combination = if_else(set == combination, "disjoint", combination )) %>%
  group_by(set) %>% 
  ggplot(aes(
    x = 1, 
    y = r, 
    fill = combination, 
    group = combination)
    ) +
  geom_bar(stat = "identity", colour = "black", linetype = "dashed", alpha = .8, size = 0.1) +
  geom_label(aes(y = r, 
                label = if_else(r > ${params.min_comb_freq}, 
                                  if_else(str_ends(combination, set), 
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
            position = position_stack(vjust = .5), angle = 0, size = 3.5) + 
  facet_grid(. ~ set) +
  scale_fill_viridis_d(option = "B") +
  scale_y_continuous(expand = c(0,0), name = "Percentage", breaks = c(0, 1), limits = c(0, 1), labels = c("0%", "100%")) +
 guides(fill = guide_legend(title = "Combination", ncol = length(sets), byrow = T)) +
  theme(
    legend.position = "bottom",
    panel.background = element_blank(),
    
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    
    legend.key.size = unit(4, "mm"),
    legend.text = element_text(face = "bold"),

    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 10),
    axis.line.y = element_line()
    ) -> p
     
 ggsave("${params.tag}_barPlot.png", plot = p, device = "png", width = 8, height = 8)
  """
}