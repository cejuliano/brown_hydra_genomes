library(phytools)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#get list of filenames for orthofinder trees
treeFiles <- list.files("Results_Sep15_1/Resolved_Gene_Trees/", full.names = F)

unlink("geneTreePDF", recursive = T)
dir.create("geneTreePDF",showWarnings = F)

#create a pdf for each orthofinder tree
for (treeFile in treeFiles) {
  #extract orthogroup name from tree filename
  ogName <- gsub("_tree.txt","",treeFile)
  
  #import tree
  tree <- read.newick(file = paste0("Results_Sep15_1/Resolved_Gene_Trees/",treeFile))
  
  #fix certain problematically formated gene IDs to make the tree more readable
  tree$tip.label <- gsub("_Parent_.*","",tree$tip.label)
  tree$tip.label <- gsub("Sc4wPfr_.*ID_","",tree$tip.label)
  tree$tip.label <- gsub("mediterranea_.*gene","mediterranea_gene",tree$tip.label)
  tree$tip.label <- gsub(" protein AED.*","",tree$tip.label)
  tree$tip.label <- gsub("_annot","",tree$tip.label)
  
  #export pdf of tree plot
  pdf(file = paste0("geneTreePDF/",ogName,"_tree.pdf"), width = 35, height = (0.2 * length(tree$tip.label)))
  plotTree(tree)
  add.scale.bar(length = 0.1, lwd = 2)
  #nodelabels(text=tree$node.label, adj = c(1.25,-1.25), frame = 'none', cex = 0.8)
  dev.off()
}
