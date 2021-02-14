### Configuramos el entorno de los resultados ###
args <- commandArgs(TRUE)
DATA_DIR <- paste(args[1], "/results", sep = "")
setwd(DATA_DIR)

### Configuramos el entorno de las librerÃ?as ###
print(getwd())
.libPaths(c(.libPaths(), paste(getwd(),"software/deps_r",sep="/")))
library("STRINGdb")
library("ggplot2")
library("org.Hs.eg.db")
library("enrichplot")
library("ggnewscale")
print(.libPaths())

set.seed(365)

##################### MAPEO ADRIAN #############################

string_db <- STRINGdb$new(version = "11", 
                          species = 9606, 
                          score_threshold = 400, 
                          input_directory = "")

proteins.table <- read.delim(file = 'files/string_protein_annotations.tsv', sep = '\t', header = TRUE)
proteins.names <- data.frame(proteins.table[1])
colnames(proteins.names)[1] <- "genes"
proteins.mapped.table <- string_db$map(proteins.names, "genes", removeUnmappedRows = TRUE)
proteins.mapped.string_ids <- proteins.mapped.table$STRING_id

string.network <- string_db$get_graph()
proteins.mapped.network <- string_db$get_subnetwork(proteins.mapped.string_ids)

proteins.mapped.network.comp <- igraph::components(proteins.mapped.network)
nodes_to_remove <- names(proteins.mapped.network.comp$membership[proteins.mapped.network.comp$membership!=1]) 

proteins.mapped.network <- igraph::delete_vertices(proteins.mapped.network, nodes_to_remove)

nombres_cortos <- gsub("9606.ENSP00000", "", igraph::V(proteins.mapped.network)$name)
igraph::V(proteins.mapped.network)$name <- nombres_cortos

# Graficamos la red actual
dev.new(width = 750, height = 530, unit = "px")
png("proteins_mapped_network.png")
plot(proteins.mapped.network,
     vertex.color = "tomato",
     vertex.size = igraph::degree(proteins.mapped.network)/10,
     vertex.label.color = "black",
     vertex.label.family = "sans",
     vertex.label.cex = 0.5#,
     # layout="layout.kamada.kawai"
)
dev.off()



######################## ENRIQUECIMIENTO ADRIAN #############################

print("##################################################################")
print("##################################################################")
print("Enriquecimiento Funcional con STRINGdb")

# El enriquecimiento se aplicarÃ¡ a aquellas comunidades que tengan un mayor grado de conectividad, 
# ya que resultan ser las mÃ¡s interesantes. Por ello hacemos lo siguiente:
best.communities <- sort(community.connectedness, decreasing = TRUE)[1:5]

# DiseÃ±amos la funciÃ³n de enriquecimiento:
enriquecimiento <- function(cluster) {
  # Observamos la representaciÃ³n de los genes del cluster
  png(paste("function_red_cluster_plot_", cluster, ".png", sep = ""))
  plot(proteins.mapped.network.lc, type = "graph", clusterids = cluster)
  dev.off()
  
  # Extraemos los ids de las proteÃ?nas del cluster
  proteins.cluster.string_ids <- paste0("9606.ENSP00000", linkcomm::getNodesIn(proteins.mapped.network.lc, clusterids = cluster, type = "names"))
  
  # Enriquecimiento funcional con GO
  enrichmentGO <- string_db$get_enrichment(proteins.cluster.string_ids, category = "Process", methodMT = "fdr", iea = TRUE)
  enrichmentGO$ontology <- rep("GO")
  
  # Enriquecimiento funcional con KEGG
  enrichmentKEGG <- string_db$get_enrichment(proteins.cluster.string_ids, category = "KEGG", methodMT = "fdr", iea = TRUE)
  enrichmentKEGG$ontology <- rep("KEGG")
  
  # Enriquecimiento funcional con Pfam
  enrichmentPfam <- string_db$get_enrichment(proteins.cluster.string_ids, category = "Pfam", methodMT = "fdr", iea = TRUE)
  enrichmentPfam$ontology <- rep("Pfam")
  
  # Juntamos las tres
  enrichment <- rbind(enrichmentGO, enrichmentKEGG, enrichmentPfam)
  data.table::setcolorder(enrichment, c(11, c(1:10)))
  
  # Agrupamos filas repetidas indicando las ontologÃ?as de origen
  for (i in enrichment$term[!duplicated(enrichment$term)]) {
    num_filas <- which(enrichment$term == i)
    ontologias <- paste(enrichment[num_filas, 1], collapse = ", ")
    nueva_fila <- enrichment[num_filas[1],]
    nueva_fila$ontology <- ontologias
    enrichment <- enrichment[-num_filas, ]
    enrichment <- rbind(enrichment, nueva_fila)
  }
  
  return(enrichment)
}

res_79 <- enriquecimiento(79)
# View(res_79)
write.csv(res_79, "enriquecimiento_funcional_c079.csv")

res_84 <- enriquecimiento(84)
# View(res_84)
write.csv(res_84, "enriquecimiento_funcional_c084.csv")

res_58 <- enriquecimiento(58)
# View(res_58)
write.csv(res_58, "enriquecimiento_funcional_c058.csv")

# res_104 <- enriquecimiento(104)
# View(res_104)
# write.csv(res_104, "enriquecimiento_funcional_c104.csv")

res_83 <- enriquecimiento(83)
# View(res_83)
write.csv(res_83, "enriquecimiento_funcional_c083.csv")

