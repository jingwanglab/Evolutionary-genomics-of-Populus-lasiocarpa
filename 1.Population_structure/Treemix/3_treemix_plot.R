#install.packages('OptM')
library(OptM)
setwd('./treemix/treemix_out')
library(RColorBrewer)
library(R.utils)
prefix='TreeMix'
source("./treemix/plotting_funcs.R") # here you need to add the path
pdf('./treemix/treemix.tree.pdf',width=15,height=10)
par(mfrow=c(2,2))
for(edge in 0:3){
  plot_tree(cex=2,paste0(prefix,".",edge))
  title(paste(edge,"edges"))
}
dev.off()
pdf('./treemix/treemix.heatmap.pdf',width=1,height=10)
par(mfrow=c(2,2))
for(edge in 0:3){
  plot_resid(stem=paste0(prefix,".",edge),pop_order='./treemix/pop_order.txt')
}
dev.off()
