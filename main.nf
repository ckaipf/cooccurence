#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
params.files = [
"example/1.gff",
"example/2.gff",
"example/3.gff"
]
*/

params.files = [
"example/genes.gff",
"example/promoters.gff",
"example/terminators.gff"
]

params.default_bedtools_parameters = "-s -k 1"
params.default_distance = 50
params.config = "example.config"
params.tag = "example_run_RegulonDB"
params.min_comb_freq = 0.05
params.k = 6

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
    map { it -> (it[4] == null) ? it[0..3] + params.default_bedtools_parameters : it } |  \
    closestBed | \
    join(distances, by: [0,1], remainder: true) | \
    map { it -> (it[3] == null) ? it[0..2] + params.default_distance : it } | \
    monitorParams

    monitorParams.out.data | \
    rearrange | \
    map { it -> it.collect {x -> '"' + x +'"'}} | \
    map { it -> [it[0], it[1], it[3], it[2]] }  | \
    toList | \
    buildCompleteGraphs

    
    monitorParams.out.warnings.toList().forEach { if(it) print it.trim() }
    //buildCompleteGraphs.out.warnings.view()

    emit: 
      csv = buildCompleteGraphs.out.csv
}

workflow plot {
  main:
    cooccurrence(params.files, params.config)
    cooccurrence.out.csv | (plotVenn & barPlot)
}

/*
    k is not large enough if
      exists max(observed_distances) < parameter(distance)
*/
process monitorParams {
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
    if(\$5<${parameters} && \$6==${params.k})
      {print "[Warning][monitorParams] Pairs may have been skipped for ${i} and ${j}. Consider to increase the parameter k.";exit 0}
  }' 
  """
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
  input:
  tuple val(i), val(j), path(a), path(b), val(parameters)
  output:
  tuple val(i), val(j), path("*.closest")


"""
bedtools closest -D a ${parameters} -k ${params.k} -a ${a} -b ${b} > ${i}${j}.closest
"""
}

process rearrange {
  input:
  tuple val(i), val(j), path(file), val(parameters)
  output:
  tuple val(i), val(j), path("*.csv"), val(parameters)

"""
awk  'BEGIN{OFS=","} {print "${i}_"\$4"_"\$5"_"\$7, "${j}_"\$13"_"\$14"_"\$16, \$19}' ${file} > ${file}.csv
"""
}


process buildCompleteGraphs {
  publishDir ".", mode: "copy"
  input:
  val file
  output:
  path "*.csv", emit: csv
    stdout emit: warnings

  """
  #!/usr/bin/env python3
import collections
import os
import csv

# Disjoint-set data structure, without ranks
# Group 2-tuples if they share a parent node
# Partition of edges in connected components
class DisjointSet:
    def __init__(self):
        self.parent_pointer = collections.defaultdict(lambda: None)
        self.groups = collections.defaultdict(set)
    
    def find(self, x):
        curr_parent = self.parent_pointer[x]
        if curr_parent:
            updated_parent = self.find(curr_parent)
            self.parent_pointer[x] = updated_parent
            return updated_parent
        return x

    def union(self, x, y):
        parent_x, parent_y = self.find(x), self.find(y)
        if parent_x != parent_y:
            self.parent_pointer[parent_x] = parent_y
    
    def group(self):
        for x in self.parent_pointer:
            self.groups[self.find(x)].add(x)

# Fold: build complete graphs (k_n) from (k_n-1)
#   For each k_n
#    if exists v not in k_n: k_n U v = n and
#    for all u in k_n exists (u,v) then k_n+1 exists
def complete_graphs(vs: list, edges: set) -> set:
    # Init
    complete_graphs = set([frozenset([v]) for v in vs])
    acc = set()
    def f(complete_graphs: set, edges: set, vs: list, n: int, acc: list) -> set:
      # Complete graph of max |V| vertices
      if n == len(vs) + 1: 
        return acc
      k_n_1 = set()
      for v in vs:
        for graph in complete_graphs:
            # K_n graph has n vertices  and for all u in graph, exists edges (u,v)
            if len(graph.union({v})) == n and all(any((i,j) in edges for i,j in [(v,u),(u,v)]) for u in graph):
                k_n_1.add(graph.union({v}))
                if graph in acc:
                  acc.remove(graph)
                if {v} in acc:
                  acc.remove({v})
                acc.add(graph.union({v}))
      return f(k_n_1, edges, vs, n+1, acc)
    return f(complete_graphs, edges, vs, 1, acc)

# Connect DisjointSet and complete graphs computations 
def complete_graphs_in_components(edges: list) -> list:
    disjoint_edges = DisjointSet()
    # Build tree
    for u,v in edges:
        disjoint_edges.union(str(u), str(v))
    # Build components
    disjoint_edges.group()
    components = list()
    for v in disjoint_edges.groups.values():
        v.discard(None)
        components.append(v)
    return [complete_graphs(component, edges) for component in components]

def to_data_frame(complete_graphs: list, prefixes: list) -> dict:
    acc = []
    col = {prefix: [] for prefix in prefixes}
    for component in complete_graphs:
        for graph in component:
            for k in col.keys():
                # Add exception here if > 1, there should not be graphs with edges from the same set
                 l = [v for v in graph if v.startswith(k)]
                 if l:
                    col[k].append(l[0]) 
                 else:
                    col[k].append(None)
    return col


# IO
ids = set()
edges = set()
vs = set()
for file in ${file}:
  ids.add(file[1])
  with open(file[3]) as f:
      for line in f.readlines():
          v,u,distance = line.strip().split(",")
          vs.add(v)
          vs.add(u)
          if abs(int(distance)) < abs(int(file[2])):
            edges.add((v,u))

# To return also k_0 graphs, loops are added
# Should slow down the computation and should be replaced by something more efficient
edges.update([(vertice,vertice) for vertice in vs])

ids = list(ids)
gs = complete_graphs_in_components(edges)
cols = to_data_frame(gs, ids)

with open("groups.csv", "w") as the_file:
     wr = csv.writer(the_file)
     wr.writerow(ids)
     for x in zip(*[v for v in cols.values()]): 
        wr.writerow(x)
  """
}

process plotVenn{
  publishDir ".", mode: "copy"
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

  ggsave("${params.tag}_venn.png", plot = p, device = "png")
  save(venn_complete, sets, file = "${params.tag}.RData")
  """
}

process barPlot {
  publishDir ".", mode: "copy"
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