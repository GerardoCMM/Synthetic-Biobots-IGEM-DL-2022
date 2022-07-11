setwd("/home/gerardo/Documentos/Gera/SyntheticBiobots/2022/Metods/Exon_usage")

library("DEXSeq")

countFiles <- c('./ERR4099794.tabular',
               './ERR4099795.tabular',
               './ERR4099796.tabular',
               './ERR4099797.tabular',
               './ERR4099798.tabular',
               './ERR4099799.tabular',
               './ERR4099800.tabular',
               './ERR4099801.tabular',
               './ERR4099802.tabular',
               './ERR4099803.tabular',
               './ERR4099804.tabular',
               './ERR4099805.tabular')

#flattenedFile = "genome_flattened.gff"

flattenedFile = "./Transcripts/genome_flattened_2.gff"

sampleTable <- data.frame(condition = factor(c(rep("fruit_20d",3),rep("fruit_40d",3),rep("leaf",3),rep("panicle",3))))

sampleTable$condition <- relevel(sampleTable$condition,ref="fruit_40d")

dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

genesForSubset = c(rownames(difInf40Counts),"Pn3.4770")

dxd = dxd[geneIDs( dxd ) %in% genesForSubset,]

dxd = estimateSizeFactors( dxd )

dxd = estimateDispersions( dxd )

plotDispEsts( dxd )

dxd = testForDEU( dxd )

dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

dxr1 = DEXSeqResults( dxd )

dfdxr1 = as.data.frame(dxr1)

table ( dxr1$padj < 0.05 )

table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )

plotMA( dxr1, cex=0.8 )

plotDEXSeq( dxr1, "Pn16.1237", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "Pn16.1237", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "Pn16.1237", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "Pn16.1237", splicing = TRUE, norCounts=FALSE,expression=FALSE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )



plotDEXSeq( dxr1, "Pn2.2377", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "Pn2.2377", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "Pn2.2377", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "Pn2.2377", splicing = TRUE, norCounts=FALSE,expression=FALSE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )




plotDEXSeq( dxr1, "Pn3.4770", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "Pn3.4770", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "Pn3.4770", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "Pn3.4770", splicing = TRUE, norCounts=FALSE,expression=FALSE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "Pn3.4770", splicing = TRUE, displayTranscripts=TRUE, norCounts=FALSE,expression=FALSE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "Pn3.4770", splicing = TRUE, displayTranscripts=TRUE, norCounts=FALSE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
