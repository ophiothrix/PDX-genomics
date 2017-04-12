# contrast.levels <- combn(levels(x.filtered$samples$group), 2)
# 
# contrast <- contrast.levels[,3]
# y <- x.filtered
# 
# compare.conditions.edgeR.QLF(x.filtered, contrast.levels[,4], target.logFC = 2)
# 
# apply(contrast.levels, 2, function(v) compare.conditions.edgeR.QLF(x.filtered, v, target.logFC = 0))


# Differential expression analysis with edgeR -----------------------------

## The function takes a filtered DGEList object and a vector of two group names and preforms differential expression analysis of Group1 vs Group2 using edgeR glm quasi-likelihood method.
# The function writes summary plots (MA plot and Volcano plot) and results for all the genes in the DGEList object (DE and non-DE).
## The group names in the vector supplied to the function should match the group names in the x$samples$group, where "x" is the DGEList object used.
## There is an option to specify a target logFC threshold to limit the results to the genes with the greater effect size. Rather that doing a combination of FDR and logFC cut off - the function actually tests for significant deviation from the target logFC. The default is logFC 0.
compare.conditions.edgeR.QLF <- function(y, contrast, target.logFC = 0, prefix = "") {
	require(edgeR)
	require(limma)
	print(paste0("Performing differential expression analysis with edgeR for ", contrast[1], " VS ", contrast[2], " against target logFC=", target.logFC))
	
	## Fit the quasi-likelihood glm
	fit <- glmQLFit(y, des, robust=TRUE)
	## Select the conditions to be compared
	contr <- paste(contrast, collapse = " - ")
	contr <- do.call(makeContrasts, list(contr, levels = des))
	
	## Compute quasi-likelihood test statistics for the given contrast against the target logFC
	results <- glmTreat(fit, contrast = contr, lfc = target.logFC)
	
	## Correct the p values
	top.table <- topTags(results, n = nrow(results), sort.by = "none")$table
	top.table$status <- 0
	top.table$status[top.table$FDR < 0.05 & top.table$logFC > 0] <- 1
	top.table$status[top.table$FDR < 0.05 & top.table$logFC < 0] <- -1
	status <- top.table$status
	table(status)
	
	## Plot summary plots
	pdf(paste0("./graphs/", prefix, contrast[1], "vs", contrast[2], ".edgeR.logFC=", target.logFC, ".summary.plots.pdf"), 16, 8)
	par(mfrow = c(1, 2))
	
	## Plot MA plot marking genes with FDR < 0.05
	plot(top.table$logCPM, top.table$logFC, type = "n", main = paste0(contrast[1], " vs ", contrast[2], "\n Up - ", sum(status== 1), ", Down - ", sum(status == -1), "\ntarget logFC = ", target.logFC), xlab = "Average expression (log2 CPM)", ylab = "logFC")
	points(top.table$logCPM[status == 0], top.table$logFC[status == 0], pch = 16, cex = 0.3, col = "grey")
	points(top.table$logCPM[status == 1], top.table$logFC[status == 1], col = "#8B000099", pch = 16, cex = 0.7)
	points(top.table$logCPM[status == -1], top.table$logFC[status == -1], col = "#00008B99", pch = 16, cex = 0.7)
	
	## Plot Volcano plot - -log10 raw p values vs logFC, mark genes with FDR < 0.05
	plot(top.table$logFC, -log10(top.table$PValue), type = "n", ylab = "-log10 uncorrected P value", xlab = "logFC", main = paste0(contrast[1], " vs ", contrast[2], "\n Up - ", sum(status== 1), ", Down - ", sum(status == -1), "\ntarget logFC = ", target.logFC))
	points(top.table$logFC[status == 0], -log10(top.table$PValue)[status == 0], cex = 0.3, col = "grey", pch = 16)
	points(top.table$logFC[status == 1], -log10(top.table$PValue)[status == 1], col = "#8B000099", pch = 16, cex = 0.7)
	points(top.table$logFC[status == -1], -log10(top.table$PValue)[status == -1], col = "#00008B99", pch = 16, cex = 0.7)
	dev.off()
	# Write out the results table
	write.csv(top.table[order(top.table$PValue),], file = paste0("./reports/", prefix, contrast[1], "vs", contrast[2], ".edgeR.logFC=", target.logFC, ".results.table.csv"), quote = T)
	return(top.table[order(top.table$PValue),])
}

# Differential expression analysis with limma voom -----------------------------

## The function takes a filtered DGEList object and a vector of two group names and preforms differential expression analysis of Group1 vs Group2 using limma voom pipeline with sample weights method.
# The function writes summary plots (MA plot and Volcano plot) and results for all the genes in the DGEList object (DE and non-DE).
## The group names in the vector supplied to the function should match the group names in the x$samples$group, where "x" is the DGEList object used.
## There is an option to specify a target logFC threshold to limit the results to the genes with the greater effect size. Rather that doing a combination of FDR and logFC cut off - the function actually tests for significant deviation from the target logFC. The default is logFC 0.
compare.conditions.limma.voom <- function(y, contrast, target.logFC = 0, use.sample.weights = T, prefix = "") {
	require(limma)
	print(paste0("Performing differential expression analysis with limma voom for ", contrast[1], " VS ", contrast[2], " against target logFC=", target.logFC))
	
	# 	des <- model.matrix(~0 + y$samples$group)
	# 	colnames(des) <- levels(y$samples$group)
	
	## Compute gene-wise precision weights and sample weights
	if (use.sample.weights) {
		v <- voomWithQualityWeights(y, des, plot=F)
	} else {
		v <- voom(y, des, plot=F)
	}
	
	## Fit linear model using precision weights and sample weights (if selected) from voom
	vfit <- lmFit(v, des)
	## Select the conditions to be compared
	contr <- paste(contrast, collapse = " - ")
	contr <- do.call(makeContrasts, list(contr, levels = des))
	
	## Compute the empirical Bayes statistic using the selected target logFC
	vfit <- contrasts.fit(vfit, contrasts = contr)
	results <- treat(vfit, lfc = target.logFC)
	
	## Correct the p values
	top.table <- topTreat(results, n = nrow(results), sort.by = "none")
	head(top.table)
	top.table$status <- 0
	top.table$status[top.table$adj.P.Val < 0.05 & top.table$logFC > 0] <- 1
	top.table$status[top.table$adj.P.Val < 0.05 & top.table$logFC < 0] <- -1
	status <- top.table$status
	table(status)
	
	## Plot summary plots
	pdf(paste0("./graphs/", prefix, contrast[1], "vs", contrast[2], ".voom.logFC=", target.logFC, ".summary.plots.pdf"), 16, 8)
	par(mfrow = c(1, 2))
	
	## Plot MA plot marking genes with FDR < 0.05
	plot(top.table$AveExpr, top.table$logFC, type = "n", main = paste0(contrast[1], " vs ", contrast[2], "\n Up - ", sum(status== 1), ", Down - ", sum(status == -1), "\ntarget logFC = ", target.logFC), xlab = "Average expression", ylab = "logFC")
	points(top.table$AveExpr[status == 0], top.table$logFC[status == 0], pch = 16, cex = 0.3, col = "grey")
	points(top.table$AveExpr[status == 1], top.table$logFC[status == 1], col = "#8B000099", pch = 16, cex = 0.7)
	points(top.table$AveExpr[status == -1], top.table$logFC[status == -1], col = "#00008B99", pch = 16, cex = 0.7)
	
	## Plot Volcano plot - -log10 raw p values vs logFC, mark genes with FDR < 0.05
	plot(top.table$logFC, -log10(top.table$P.Value), type = "n", ylab = "-log10 uncorrected P value", xlab = "logFC", main = paste0(contrast[1], " vs ", contrast[2], "\n Up - ", sum(status== 1), ", Down - ", sum(status == -1), "\ntarget logFC = ", target.logFC))
	points(top.table$logFC[status == 0], -log10(top.table$P.Value)[status == 0], cex = 0.3, col = "grey", pch = 16)
	points(top.table$logFC[status == 1], -log10(top.table$P.Value)[status == 1], col = "#8B000099", pch = 16, cex = 0.7)
	points(top.table$logFC[status == -1], -log10(top.table$P.Value)[status == -1], col = "#00008B99", pch = 16, cex = 0.7)
	dev.off()
	# Write out the results table
	write.csv(top.table[order(top.table$P.Value),], file = paste0("./reports/", prefix, contrast[1], "vs", contrast[2], ".voom.logFC=", target.logFC, ".results.table.csv"), quote = T)
	return(top.table[order(top.table$P.Value),])
}


# Differential expression analysis with removal of unwanted variance (RUV) -----------------------------

## The function takes a filtered DGEList object and a vector of two group names and preforms differential expression analysis of Group1 vs Group2 using RUV to remove the unwanted variation (e.g. batch effects) and quasi-likelihood edgeR pipeline to compute the statistics.
# The function writes summary plots (MA plot and Volcano plot) and results for all the genes in the DGEList object (DE and non-DE).
## The group names in the vector supplied to the function should match the group names in the x$samples$group, where "x" is the DGEList object used.
## There is an option to specify a target logFC threshold to limit the results to the genes with the greater effect size. Rather that doing a combination of FDR and logFC cut off - the function actually tests for significant deviation from the target logFC. The default is logFC 0.
compare.conditions.RUV <- function(y, contrast, target.logFC = 0, prefix = "") {
	require(edgeR)
	require(RUVSeq)
	print(paste0("Performing differential expression analysis with RUV for ", contrast[1], " VS ", contrast[2], " against target logFC=", target.logFC))
	
	## First we run a standard edgeR analysis to detect least DE genes. We do not consider the target logFC for that step.
	print("Performing first pass analysis to derive empirical negative gene set")
	## Fit the quasi-likelihood glm
	fit <- glmQLFit(y, des, robust=TRUE)
	## Select the conditions to be compared
	contr <- paste(contrast, collapse = " - ")
	contr <- do.call(makeContrasts, list(contr, levels = des))
	
	## Compute quasi-likelihood test statistics for the given contrast against the target logFC
	results <- glmQLFTest(fit, contrast = contr)
	
	## Correct the p values
	top.table <- topTags(results, n = nrow(results), sort.by = "none")$table
	## Get the proportion of DE genes at FDR 0.05
	
	## Halve the number of all non-DE genes and take that number of genes with the highest p values as the empirical negative control set for RUV
	neg.set <- round(sum(top.table$FDR > 0.05)/2, 0)
	empirical <- rownames(top.table)[order(top.table$PValue, decreasing = T)][1:neg.set]
	## Check that the order of genes in the counts object and in the top.table is the same
	if (!all(rownames(y$counts) == rownames(top.table))) {
		stop("The order of the genes in the DGEList does not match that of the DE results table")
	}
	
	## Make an expressionSet object for RUV
	set <- newSeqExpressionSet(as.matrix(y$counts), phenoData = data.frame(y$samples$group, row.names=colnames(y$counts)))
	
	# colors <- brewer.pal(3, "Set2")
	# plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[groups])
	# plotPCA(set, col=colors[groups], cex=1.2)
	
	# Normalise the expression values
	set <- betweenLaneNormalization(set, which="upper")
	# plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[groups]) 
	# plotPCA(set, col=colors[groups], cex=1.2)
	
	## Compute the "batch factors" for the dataset using RUVg
	set2 <- RUVg(x = set, cIdx = empirical, k=1)
	print(pData(set2))
	
	print("Performing second pass analysis using RUV generated normalisation factors as a covariate.")
	## Make a new design matrix and include the RUV normalisation factor as a covariate
	des1 <- model.matrix(~0 + x.filtered$samples$group + pData(set2)$W_1)
	colnames(des1) <- c(levels(x.filtered$samples$group), "W_1")
	
	## Fit the quasi-likelihood glm
	fit <- glmQLFit(y, des1, robust=TRUE)
	
	## Select the conditions to be compared
	contr <- paste(contrast, collapse = " - ")
	contr <- do.call(makeContrasts, list(contr, levels = des1))
	
	## Compute quasi-likelihood test statistics for the given contrast against the target logFC
	results <- glmTreat(fit, contrast = contr, lfc = target.logFC)
	
	## Correct the p values
	top.table <- topTags(results, n = nrow(results), sort.by = "none")$table
	top.table$status <- 0
	top.table$status[top.table$FDR < 0.05 & top.table$logFC > 0] <- 1
	top.table$status[top.table$FDR < 0.05 & top.table$logFC < 0] <- -1
	status <- top.table$status
	table(status)
	
	## Plot summary plots
	pdf(paste0("./graphs/", prefix, contrast[1], "vs", contrast[2], ".RUV.logFC=", target.logFC, ".summary.plots.pdf"), 16, 8)
	par(mfrow = c(1, 2))
	
	## Plot MA plot marking genes with FDR < 0.05
	plot(top.table$logCPM, top.table$logFC, type = "n", main = paste0(contrast[1], " vs ", contrast[2], "\n Up - ", sum(status== 1), ", Down - ", sum(status == -1), "\ntarget logFC = ", target.logFC), xlab = "Average expression (log2 CPM)", ylab = "logFC")
	points(top.table$logCPM[status == 0], top.table$logFC[status == 0], pch = 16, cex = 0.3, col = "grey")
	points(top.table$logCPM[status == 1], top.table$logFC[status == 1], col = "#8B000099", pch = 16, cex = 0.7)
	points(top.table$logCPM[status == -1], top.table$logFC[status == -1], col = "#00008B99", pch = 16, cex = 0.7)
	
	## Plot Volcano plot - -log10 raw p values vs logFC, mark genes with FDR < 0.05
	plot(top.table$logFC, -log10(top.table$PValue), type = "n", ylab = "-log10 uncorrected P value", xlab = "logFC", main = paste0(contrast[1], " vs ", contrast[2], "\n Up - ", sum(status== 1), ", Down - ", sum(status == -1), "\ntarget logFC = ", target.logFC))
	points(top.table$logFC[status == 0], -log10(top.table$PValue)[status == 0], cex = 0.3, col = "grey", pch = 16)
	points(top.table$logFC[status == 1], -log10(top.table$PValue)[status == 1], col = "#8B000099", pch = 16, cex = 0.7)
	points(top.table$logFC[status == -1], -log10(top.table$PValue)[status == -1], col = "#00008B99", pch = 16, cex = 0.7)
	dev.off()
	# Write out the results table
	write.csv(top.table[order(top.table$PValue),], file = paste0("./reports/", prefix, contrast[1], "vs", contrast[2], ".RUV.logFC=", target.logFC, ".results.table.csv"), quote = T)
	return(top.table[order(top.table$PValue),])
}

## Given the contrast (a vector of 2 group names) and target logFC, plot a venn diagram of the overlap between the genes
## This function relies on the summary tables having been produced for limma.voom, edgeR and RUV methods
plot.gene.overlaps <- function(contrast, target.logFC) {
	require(Vennerable)
	limma.voom <- read.csv(paste0("./reports/", contrast[1], "vs", contrast[2], ".voom.logFC=", target.logFC, ".results.table.csv"))
	head(limma.voom)
	
	edger <- read.csv(paste0("./reports/", contrast[1], "vs", contrast[2], ".edgeR.logFC=", target.logFC, ".results.table.csv"))
	head(edger)
	
	ruv <- read.csv(paste0("./reports/", contrast[1], "vs", contrast[2], ".RUV.logFC=", target.logFC, ".results.table.csv"))
	head(ruv)
	
	venn.obj <- Venn(list(edgeR = edger$genes[edger$FDR < 0.05], limma.voom = limma.voom$genes[limma.voom$adj.P.Val < 0.05], RUV = ruv$genes[ruv$FDR < 0.05]))
	pdf(paste0("./graphs/", "gene.set.overlap.", contrast[1], "vs", contrast[2], ".logFC=", target.logFC, ".pdf"))
	plot(venn.obj)
	dev.off()
}
