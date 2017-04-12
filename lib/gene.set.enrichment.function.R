# x.filtered <- readRDS("./cache/filtered.counts.DGEList.rds")
# go.list <- readRDS("./data/GO.gene.sets.rds")
# contrast <- contrast.levels[,1]
# go.enrichment(counts.obj = x.filtered, go.list = go.list, contrast = contrast)
# counts.obj <- x.filtered

go.enrichment <- function(counts.obj, design, go.list=go.list, contrast) {
	require(limma)
	print(paste0("Running gene set enrichment analysis for ", paste(contrast, collapse = " vs ")))
	## load DGEList
	## load gene set object
	
	## Convert the genes ids in the gene set object to indices
	go.list.idx <- ids2indices(gene.sets = go.list, identifiers = counts.obj$genes$GeneID, remove.empty = T)
	head(go.list.idx)
	
	## set up the contrasts
	contr <- paste(contrast, collapse = " - ")
	contr <- do.call(makeContrasts, list(contr, levels = des))
	
	camera.x <- camera(y = counts.obj, index = go.list.idx, design = design, contrast = contr)
	# head(camera.x)
	
	## Annotate results
	# head(go.anno)
	camera.x <- cbind(camera.x, go.anno[match(rownames(camera.x), go.anno$go.term), -1])
	# head(camera.x)
	
	
	write.csv(camera.x, file = paste0("./reports/", contrast[1], "vs", contrast[2], ".GO.enrichment.camera.csv"), quote = T)
}

# summary(camera.x$FDR < 0.05)
# roast.x <- mroast(x.filtered, index = go.list.idx, design = des, contrast = contr, nrot = 10000)
# ?roast
# head(roast.x, 30)
# 
# roast.x <- romer(x.filtered, index = go.list.idx, design = des, contrast = contr)
# head(roast.x)

# de.table <- read.csv("./reports/IMAGOvsME.voom.logFC=2.results.table.csv")
# head(de.table)
# go.list.idx.de <- ids2indices(gene.sets = go.list, identifiers = de.table$genes, remove.empty = T)
# geneSetTest(index = go.list.idx.de[[1]], statistics = de.table$t, type = "t", alternative = "either")
# ?geneSetTest
