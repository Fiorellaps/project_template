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
library("readxl")
library("igraph")
library("linkcomm")
library("data.table")
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

######################## LINK COMMUNITIES FIORELLA #############################
proteins.mapped.network.df <- igraph::as_data_frame(proteins.mapped.network, what="edges")
proteins.mapped.network.df<-proteins.mapped.network.df[c("from","to","combined_score")]
proteins.mapped.network.lc <- linkcomm::getLinkCommunities(proteins.mapped.network.df, hcmethod = "single",plot=FALSE)
print(proteins.mapped.network.lc)

# link communities summary plot
png(file="linkcomm_summary.png")
plot(proteins.mapped.network.lc, type = "summary")
dev.off()
# link communities dendogram
# png(file="linkcomm_dend.png")
# plot(proteins.mapped.network.lc, type = "dend")
# dev.off()

# link communities Fruchterman Reingold layout
# png(file="linkcomm_layout.fruchterman.reingold.png")
# png(file="linkcomm_graph.png")
# plot(proteins.mapped.network.lc, type = "graph", vlabel=FALSE)#Le he quitado , layout = "layout.fruchterman.reingold"
# dev.off()

# link communities Spencer Circle layout
# png(file="linkcomm_spencer.circle.png")
# plot(proteins.mapped.network.lc, type = "graph", layout = "spencer.circle", vlabel=FALSE)
# dev.off()

# link communities Spencer Circle layout
# displays only the nodes that belong to 5 or more communities
png(file="linkcomm_common_nodes_5.png")
plot(proteins.mapped.network.lc, type = "graph", shownodesin = 5, vlabel=FALSE) # Le he quitado  layout = "spencer.circle",
dev.off()

# link communities Node pies
png(file="linkcomm_node.pies.png")
plot(proteins.mapped.network.lc, type = "graph", shownodesin = 3, node.pies = TRUE, vlabel=FALSE)
dev.off()

# see link communities members
png(file="linkcomm_communities_members.png")
plot(proteins.mapped.network.lc, type = "members")
dev.off()

# obtener comunidades completamente anidadas dentro de una comunidad más grande de nodos
print("##################################################################")
print("##################################################################")
print("Comunidades completamente anidadas dentro de una comunidad más grande de nodos")
linkcomm::getAllNestedComm(proteins.mapped.network.lc)

# comprobar comunidades anidadas
print("##################################################################")
print("##################################################################")
print("Comprobando comunidades anidadas")
linkcomm::getAllNestedComm(proteins.mapped.network.lc)[1]

# indica que la comunidad 11 esta completamente anidad por los nodos de la comunidad 74
# png(file="linkcomm_communities_11.74.png")
# plot(proteins.mapped.network.lc, type = "graph", clusterids = c(11,74))
# dev.off()

# Visualizar relaciones entre comunidades
print("##################################################################")
print("##################################################################")
print("Visualizando relaciones entre comunidades")
lc.communities.relations <- linkcomm::getClusterRelatedness(proteins.mapped.network.lc, hcmethod = "average",cutat = 0.5) #average link method cutting at 0.5
linkcomm::getClusterRelatedness(proteins.mapped.network.lc, hcmethod = "ward.D2") #ward.D2 link method

# get smaller number of communities (meta-communities)

lc.meta.communities <- linkcomm::meta.communities(proteins.mapped.network.lc, hcmethod = "average", deepSplit = 0)

png(file="linkcomm_metacommunities_summary.png")
plot(lc.meta.communities, type = "summary")
dev.off()
png(file="linkcomm_metacommunities_graph.png")
dev.off()
plot(proteins.mapped.network.lc, type = "graph", vlabel=FALSE)#Le he quitado , layout = layout.fruchterman.reingold
png(file="linkcomm_metacommunities_members.png")
plot(lc.meta.communities, type = "members")
dev.off()


# Community centrality
print("##################################################################")
print("##################################################################")
print("Visualizando community centrality")
community.centrality <- linkcomm::getCommunityCentrality(proteins.mapped.network.lc)

# modularity of the communities
community.connectedness <- linkcomm::getCommunityConnectedness(proteins.mapped.network.lc, conn = "modularity") 
png(file="linkcomm_communities_modularity.png")
plot(proteins.mapped.network.lc, type = "commsumm", summary = "modularity")
dev.off()

# Focus on one linkcomm
# plot one cluster that I have chosen randomly
png(file="cluster12_graph.png")
plot(proteins.mapped.network.lc, type = "graph", clusterids = 12, vlabel=FALSE)
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

