library(Seurat)
library(edgeR)

seurat.data <- Seurat::Read10X("Data")
seurat <- CreateSeuratObject(seurat.data)

er <- edgeR::read10X("matrix.mtx.gz",
                     "features.tsv.gz",
                     "barcodes.tsv.gz",
                     "Data",
                     TRUE
                     )

s <- c(object.size(er), object.size(seurat))
labs <- c("EdgeR", "Seurat")

barplot(s, main="Size of Objects",
        xlab="bytes",
        names.arg = labs)

perc.mito <- PercentageFeatureSet(seurat, pattern = "^mt-")

plot(perc.mito, seurat$nCount_RNA)
