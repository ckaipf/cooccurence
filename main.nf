#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow {
  main:
  cooccurrence(params.files, params.config)
}

workflow cooccurrence {
  take: 
    files
    config
  main:
  files_ch = Channel.fromPath(files) | \
          sortGff | \
          map { it -> [it.getSimpleName(), it] }
  
  config = Channel.fromPath(config) | \
      splitCsv() | \
      map { row -> (row[0] == params.tag) ? row : null } |  \
      map { row -> (files.any { it.contains(row[1] + ".gff")} && files.any { it.contains(row[2] + ".gff")}) ? row : null} 

  parameters = config | \
    map { it -> [ it[1], it[2], it[3] ] }

  distances = config | \
    map { it -> [ it[1], it[2], it[4] ] }

  files_ch.combine(files_ch) | \
    filter(it -> it[0] != it[2]) | \
    map {it -> [it[0], it[2], it[1], it[3]]} |  \
    join(parameters, by: [0,1], remainder: true) |  \
    map { it -> (it[4] == null) ? it[0..3] + params.bedtools_default_parameters : it } |  \
    closestBed | \
    join(distances, by: [0,1], remainder: true) | \
    map { it -> (it[3] == null) ? it[0..2] + params.default_distance : it } | \
    monitorParams

    monitorParams.out.data | \
    rearrange | \
    catFiles | \
    collectFile(name: "collected_file.txt")  | \
    buildCompleteGraphs

    monitorParams.out.warnings.toList().forEach { if(it) log.warn(it.trim()) }
    buildCompleteGraphs.out.warnings.view()

    emit: 
      set_csv = buildCompleteGraphs.out.set_csv
      order_csv = buildCompleteGraphs.out.order_csv
}

workflow plot {
  main:
    cooccurrence(params.files, params.config)
    cooccurrence.out.set_csv | (plotVenn & barPlot)
    cooccurrence.out.order_csv | plotOrder
}

/*
    k is not large enough if
      exists max(observed_distances) < parameter(distance)
*/
process monitorParams {
  container 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'
  errorStrategy "ignore"
  input:
  tuple val(i), val(j), path(a), val(parameters)
  output:
  tuple val(i), val(j), path(a), val(parameters), emit: data
  stdout emit: warnings
  
  """
  sort -k1,1 -k4,4n -k5,5n -k7,7 ${a} | \
  awk 'function abs(v) {return v < 0 ? -v : v}
  { 
    for(i=1;i<=NF-1;i++) printf \$i"\t"; print abs(\$19)
  }' | \
  bedtools groupby -g 1,4,5,7 -c 19 -o max,count | \
  awk '{
    if(\$5<${parameters} && \$6==${params.bedtools_k})
      {print "Pairs may have been skipped for \\033[4m${i}\\033[24m and \\033[4m${j}\\033[24m. Consider to increase parameter k.";exit 0}
  }' 
  """
}

process catFiles {
  container 'debian:stable'
 // publishDir ".", mode: "copy"
  input:
  tuple val(i), val(j), path(csv), val(parameters)
  output:
  path "*.gz"
 """
 awk 'BEGIN{OFS=","}{print "${i}","${j}","${parameters}",\$0}' ${csv} | \
 gzip \
 > ${csv}.gz
 """
}

process sortGff {
  container 'debian:stable'

  input:
  path(a)
  output:
  path("${a}.sorted")

  """
  sort -k1,1 -k4,4n -k5,5n ${a} > ${a}.sorted
  """
}

process closestBed {
  container 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'
  input:
  tuple val(i), val(j), path(a), path(b), val(parameters)
  output:
  tuple val(i), val(j), path("*.closest")


"""
bedtools closest -D a ${parameters} -k ${params.bedtools_k}  ${params.bedtools_global} -a ${a} -b ${b} > ${i}${j}.closest
"""
}

process rearrange {
  container 'debian:stable-slim'
  input:
  tuple val(i), val(j), path(file), val(parameters)
  output:
  tuple val(i), val(j), path("*.csv"), val(parameters)

"""
awk  'BEGIN{OFS=","} {print "${i}_"\$4"_"\$5"_"\$7, "${j}_"\$13"_"\$14"_"\$16, \$19}' ${file} > ${file}.csv
"""
}


process buildCompleteGraphs {
  container 'python:slim'
  publishDir "${params.dir}/${params.tag}", mode: "copy"
  input:
  path csv
  output:
  path "*_set.csv.gz", emit: set_csv
  path "*_order.csv.gz", emit: order_csv
    stdout emit: warnings

  script:
  """
  build_complete_graphs.py ${csv} ${params.tag}
  """
}

process plotVenn{
  publishDir "${params.dir}/${params.tag}", mode: "copy"
  input:
  path csv
  output:
  path "*.png"

  """
  #!/usr/bin/Rscript
  library(tidyverse)
  library(ggVennDiagram)

  venn_complete <- read_delim(file = "${csv}", delim = ",")  
  sets <- colnames(venn_complete)
  names(sets) <- sets
  names(venn_complete) <- sets
  venn_complete <- venn_complete %>%
    rowid_to_column()

  ids <- sets %>%
    map(~ filter(venn_complete, !is.na(.data[[.]]))) %>%
    map(~ select(., rowid) %>% unlist %>% unique) 

  venn <- Venn(ids)
  data <- process_data(venn)
  region_data <- venn_region(data)
  ggplot() +
    geom_sf(aes(fill = count), alpha = .7, data = venn_region(data), show.legend = F) +
    geom_sf(size = 2, color = "grey", data = venn_setedge(data), show.legend = F) +
    geom_sf_text(aes(label = name), data = venn_setlabel(data), fontface = "bold") +
    geom_sf_label(aes(label=count), alpha = .7, fontface = "bold", data = venn_region(data)) +
    scale_fill_viridis_c() + 
    theme_void() -> p

  ggsave("${params.tag}_venn.png", plot = p, device = "png", width = 6, height = 6)
  """
}

process barPlot {
  publishDir "${params.dir}/${params.tag}", mode: "copy"
  input:
  path csv
  output:
  path "*.png"

"""
 #!/usr/bin/Rscript
  library(tidyverse)

  venn_complete <- read_delim(file = "${csv}", delim = ",")  
  sets <- colnames(venn_complete)
  names(sets) <- sets
  names(venn_complete) <- sets
  venn_complete <- venn_complete %>%
    rowid_to_column()

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
 guides(fill = guide_legend(title = "Combination", nrow = length(sets), byrow = T)) +
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

 ggsave("${params.tag}_barPlot.png", plot = p, device = "png", width = 6, height = 6)
  """
}

process plotOrder {
  publishDir "${params.dir}/${params.tag}", mode: "copy"
  input:
  path csv
  output:
  path "*.png"
 
  """
   #!/usr/bin/Rscript
    library(tidyverse)

read_delim(file = "${csv}", delim = ",") %>%
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
  ) -> p 
ggsave(plot = p, filename = "${params.tag}_freq_of_orders.png", device = "png", width = 10)
"""
}


workflow collectRuns {
  main:
  Channel.fromPath( params.dir + "/**_set.csv.gz") | \
    zcatFilesAppendFileName | \
    collect | \
    appendFiles | \
    collectedBarPlot
}

process zcatFilesAppendFileName {
//  publishDir "${params.dir}", mode: "copy"
  input:
  path f
  output:
  path "*t.csv"
 """
 zcat ${f} | \
 awk -F ',' -v OFS=',' '{if(NR==1) {print "tag",\$0} else {print "${f.simpleName}",\$0}}' \
 > ${f.simpleName}.t.csv
 """
}

process appendFiles {
  publishDir "${params.dir}", mode: "copy"
  input:
  path f
  output:
  path "collected.csv"
 """
 head -1 ${f[0]} > collected.csv
 tail -q -n +2 ${f} >> collected.csv
 """
}

process collectedBarPlot {
  publishDir "${params.dir}/", mode: "copy"
  input:
  path csv
  output:
  path "*.png"

  """
  #!/usr/bin/Rscript
  library(tidyverse)

  read_delim("${csv}", delim = ",")  %>%
  distinct(across(everything())) %>%
  mutate(across(!starts_with(c("tag")), .fns = function(x) str_split(x, "_", simplify = T)[,1])) %>%
  group_by(across(everything())) %>%
  summarise(n = n()) %>%
  unite(remove = F, !starts_with(c("tag", "n")), col = "combination", sep = "_") %>%
  mutate(combination = str_remove_all(combination, "NA_|_NA")) %>% 
  pivot_longer(values_to = "set", !starts_with(c("tag", "n", "combination"))) %>% 
  select(-name) %>%
  filter(!is.na(set)) %>% 
  group_by(set, tag) %>%
  mutate(total = sum(n), r = n / total) %>% 
  mutate(combination = if_else(set == combination, ".", combination )) %>%
  #
 ggplot(aes(
    x = tag, 
    y = r, 
    fill = combination, 
    group = combination)
  ) +
  geom_bar(width=1, stat = "identity", colour = "black", linetype = "dashed", alpha = .8, size = 0.1) +
  geom_label(aes(y = r, 
                 label = if_else(r > 0.1, 
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
  facet_grid(set ~ ., switch = "y") +
  scale_fill_viridis_d(option = "D") +
  scale_x_discrete(name = NULL, position = "top") +
  scale_y_continuous(name = NULL, labels = NULL, breaks = NULL, limits = c(0, 1)) +
  guides(fill = guide_legend(title = "Combination", nrow = 3, byrow = T)) +
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
  ) -> p
  ggsave(plot = p, filename = "collected_barplot.png", device = "png", width = 10)
  """
}