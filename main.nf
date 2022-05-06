#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


params.dir = "example"
params.files = [
  params.dir + "/genes.gff",
  params.dir + "/promoters.gff",
  params.dir + "/terminators.gff"
]


params.bedtools_default_parameters = "-s"
params.bedtools_k = 6
params.bedtools_global = ""
params.default_distance = 50

params.config = params.dir + "/example.config"
params.tag = "example_run_RegulonDB"
params.min_comb_freq = 0.05


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

    
    monitorParams.out.warnings.toList().forEach { if(it) print it.trim() }
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
      {print "[Warning][monitorParams] Pairs may have been skipped for ${i} and ${j}. Consider to increase the parameter k.";exit 0}
  }' 
  """
}

process catFiles {
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
bedtools closest -D a ${parameters} -k ${params.bedtools_k}  ${params.bedtools_global} -a ${a} -b ${b} > ${i}${j}.closest
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
  publishDir "${params.dir}/${params.tag}", mode: "copy"
  input:
  path csv
  output:
  path "*_set.csv.gz", emit: set_csv
  path "*_order.csv.gz", emit: order_csv
    stdout emit: warnings

  """
  #!/usr/bin/env python3
import collections, os, csv, gzip

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

# Fold: build k_n+1 from k_n
#   For each k_n
#    if exists v not in k_n: k_n U v = n and
#    for all u in k_n exists (u,v) then k_n+1 exists
def complete_graphs(vs: list, es: set) -> set:
    def build_kn(complete_graphs: set, n: int, acc: set) -> set:
      if n == len(vs) + 1: 
        return acc
      ks_step = set()
      for v in vs:
        for graph in complete_graphs:
            if len(graph | {v}) == n and all((v,u) in es or (u,v) in es for u in graph):
                ks_step |= {graph | {v}}
                acc |= {graph | {v}}
                acc -= {graph} - {frozenset({v})}
      return build_kn(ks_step, n+1, acc)
    return build_kn({frozenset([v]) for v in vs}, 1, set())

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

def cols_by_set(complete_graphs: list, prefixes: list) -> dict:
    acc, col = [], {prefix: [prefix] for prefix in prefixes}
    for component in complete_graphs:
        for graph in component:
            for k in col.keys():
                 l = [v for v in graph if v.startswith(k)]
                 assert len(l) <= 1, "Multiple vertices shared the same prefix"
                 if l:
                    col[k].append(l[0]) 
                 else:
                    col[k].append(None)
    return col

def cols_by_order(complete_graphs: list, prefixes: list) -> dict:
    n = len(prefixes)
    acc, col = [], {i: [i] for i in range(0,n)}
    for component in complete_graphs:
        for graph in component:
           l = [tuple(v.split("_")) for v in graph]
           # Definition of ordering, (+) i < j if i.start < j.start, (-) i < j if i.stop*-1 < j.stop*-1 
           l = sorted(l, key = lambda x: -1 * int(x[2]) if x[3] == "-" else int(x[1]))
           for i in range(0, n):
            try:
              col[i].append("_".join((str(s) for s in l[i])))
            except IndexError:
              col[i].append(None) 
    return col

# IO
ids, es, vs = set(), set(), set()
with gzip.open('${csv}', mode='rt') as csvfile:
  rows = csv.reader(csvfile, delimiter=',', quotechar='|')
  for row in rows:
    i,j,max_distance,v,u,distance = row
    vs |= {v,u}
    ids |= {i,j}
    if abs(int(distance)) < abs(int(max_distance)):
      es |= {(v,u)}

# To return also k_1 graphs, loops are added
# Should slow down the computation and should be replaced by something more efficient
es.update([(vertice,vertice) for vertice in vs])
ids, gs = list(ids), complete_graphs_in_components(es)

for f in [cols_by_set, cols_by_order]:
   cols = f(gs, ids)
   with gzip.open(".".join([f.__name__, "csv", "gz"]), "wt") as out:
     wr = csv.writer(out)
     for x in zip(*[cols[k] for k in sorted(cols.keys())]): 
        wr.writerow(x) 
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
  Channel.fromPath( params.dir + "/**_set.csv.gz") | view | \
    zcatFilesAppendFileName | \
    collect | \
    appendFiles
}

process zcatFilesAppendFileName {
//  publishDir "${params.dir}", mode: "copy"
  input:
  path f
  output:
  path "*t.csv"
 """
 zcat ${f} | \
 awk -F ',' -v OFS=',' '{if(NR==1) {print "file",\$0} else {print "${f.simpleName}",\$0}}' \
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