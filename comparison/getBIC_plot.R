library(clustNet)
library(ggplot2)
stats <- readRDS("./amlIC.rds")

# Create an empty matrix with 7 rows and 5 columns
result_matrix <- matrix(nrow = 7, ncol = 5)

# Populate the matrix with elements from stats$output
for (x in 1:35) {
  # Calculate the row and column indices for the current x
  row_index <- ((x - 1) %% 7) + 1
  col_index <- ((x - 1) %/% 7) + 1
  
  # Extract the value from stats$output and store it in the matrix
  result_matrix[row_index, col_index] <- stats$output[[x]]$testBIC
}

# Print the result_matrix
print(result_matrix)
bics <- result_matrix

chiVec <- c(0.001,0.5,1,2,3)
ncol <- length(chiVec)
minK=3
maxK=9
nrow <- maxK-minK+1
BICrange = 100

ks <- rep(0, ncol)
minbics <- rep(0, ncol)
minbics <- apply(bics, 2, min, na.rm = TRUE)
bics <- t(bics)
bics <- bics - minbics
topbics <- BICrange
bics[which(bics > topbics)] <- topbics
bicsscaled <- t(t(bics)/colSums(bics))
divergy <- bics
rownames(divergy) <- chiVec
colnames(divergy) <- c(minK:maxK)
meltdivergy <- reshape2::melt(divergy)


ggplot2::ggplot(data = meltdivergy, ggplot2::aes(x = Var1, 
                                                 y = Var2, fill = value)) + ggplot2::geom_tile()
middycol <- c(0.8, 0.2, 0)
ggheatmap <- ggplot2::ggplot(data = meltdivergy, ggplot2::aes(Var1, 
                                                              Var2, fill = value)) + ggplot2::geom_tile() + ggplot2::xlab(expression(chi)) + 
  ggplot2::ylab("k") + ggplot2::scale_fill_gradient2(high = "#FAFAFF", 
                                                     low = "#117777", mid = "#88BBBB", space = "Lab", 
                                                     na.value = "grey75", midpoint = topbics/2, limit = c(0, 
                                                                                                          topbics), name = "BIC\nchange\n") + ggplot2::scale_y_continuous(breaks = c(minK:maxK)) + 
  ggplot2::theme_minimal() + ggplot2::theme(axis.title.x = ggplot2::element_text(vjust = -1), 
                                            axis.title.y = ggplot2::element_text(angle = 0, hjust = -0.5, 
                                                                                 vjust = 0.505)) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, 
                                                                                                                                                      vjust = 0.5, size = 20, hjust = 0.6), axis.text.y = ggplot2::element_text(angle = 0, 
                                                                                                                                                                                                                                vjust = 0.5, size = 20, hjust = 1), legend.text = ggplot2::element_text(size = 20), 
                                                                                                                  axis.title = ggplot2::element_text(size = 30), legend.title = ggplot2::element_text(size = 24)) + 
  ggplot2::theme(legend.key.size = ggplot2::unit(2, 
                                                 "line")) + ggplot2::theme(plot.margin = ggplot2::unit(c(0, 
                                                                                                         0, 0.4, 0.4), "cm"))
print(ggheatmap)

# pdf(file = "./AML_BIC.pdf",width = 12, height = 6)
BIC_aml <- ggheatmap + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                           panel.border = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                           axis.ticks = ggplot2::element_blank(), legend.title=element_text(size=22), 
                           plot.title = element_text(color="black", size=24, face="bold")) +
  ggtitle("AML only")
# dev.off()

### MDS
library(clustNet)
library(ggplot2)
stats <- readRDS("./mdsIC.rds")

# Create an empty matrix with 7 rows and 5 columns
result_matrix <- matrix(nrow = 7, ncol = 5)

# Populate the matrix with elements from stats$output
for (x in 1:35) {
  # Calculate the row and column indices for the current x
  row_index <- ((x - 1) %% 7) + 1
  col_index <- ((x - 1) %/% 7) + 1
  
  # Extract the value from stats$output and store it in the matrix
  result_matrix[row_index, col_index] <- stats$output[[x]]$testBIC
}

# Print the result_matrix
print(result_matrix)
bics <- result_matrix

chiVec <- c(0.001,0.5,1,2,3)
ncol <- length(chiVec)
minK=3
maxK=9
nrow <- maxK-minK+1
BICrange = 100

ks <- rep(0, ncol)
minbics <- rep(0, ncol)
minbics <- apply(bics, 2, min, na.rm = TRUE)
bics <- t(bics)
bics <- bics - minbics
topbics <- BICrange
bics[which(bics > topbics)] <- topbics
bicsscaled <- t(t(bics)/colSums(bics))
divergy <- bics
rownames(divergy) <- chiVec
colnames(divergy) <- c(minK:maxK)
meltdivergy <- reshape2::melt(divergy)


ggplot2::ggplot(data = meltdivergy, ggplot2::aes(x = Var1, 
                                                 y = Var2, fill = value)) + ggplot2::geom_tile()
middycol <- c(0.8, 0.2, 0)
ggheatmap <- ggplot2::ggplot(data = meltdivergy, ggplot2::aes(Var1, 
                                                              Var2, fill = value)) + ggplot2::geom_tile() + ggplot2::xlab(expression(chi)) + 
  ggplot2::ylab("k") + ggplot2::scale_fill_gradient2(high = "#FAFAFF", 
                                                     low = "#117777", mid = "#88BBBB", space = "Lab", 
                                                     na.value = "grey75", midpoint = topbics/2, limit = c(0, 
                                                                                                          topbics), name = "BIC\nchange\n") + ggplot2::scale_y_continuous(breaks = c(minK:maxK)) + 
  ggplot2::theme_minimal() + ggplot2::theme(axis.title.x = ggplot2::element_text(vjust = -1), 
                                            axis.title.y = ggplot2::element_text(angle = 0, hjust = -0.5, 
                                                                                 vjust = 0.505)) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, 
                                                                                                                                                      vjust = 0.5, size = 20, hjust = 0.6), axis.text.y = ggplot2::element_text(angle = 0, 
                                                                                                                                                                                                                                vjust = 0.5, size = 20, hjust = 1), legend.text = ggplot2::element_text(size = 20), 
                                                                                                                  axis.title = ggplot2::element_text(size = 30), legend.title = ggplot2::element_text(size = 24)) + 
  ggplot2::theme(legend.key.size = ggplot2::unit(2, 
                                                 "line")) + ggplot2::theme(plot.margin = ggplot2::unit(c(0, 
                                                                                                         0, 0.4, 0.4), "cm"))
print(ggheatmap)

# pdf(file = "./MDS_BIC.pdf",width = 12, height = 6)
BIC_mds <- ggheatmap + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                           panel.border = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                           axis.ticks = ggplot2::element_blank(), legend.title=element_text(size=22), 
                           plot.title = element_text(color="black", size=24, face="bold")) +
  ggtitle("MDS only")



pdf(file = "./MDS_AML_BIC.pdf",width = 12, height = 8)
ggarrange(BIC_aml,BIC_mds, ncol = 1, nrow = 2,labels = letters[1:4])
dev.off()
