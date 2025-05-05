#! /export/software/R/R4.1.0/bin/Rscript --no-save --no-restore
.libPaths("/export/software/miniconda3.scale.envs/R4.1.0/lib/R/library")
library(LEA)
project = NULL
setwd("/export2/home/longzq/lzq/dyy/lfmm/lfmm_input")
###default burnin iter
project = lfmm("variant.lfmm", "bioclimate.txt", K = 2, repetitions = 5, CPU = 5, project = "new")
save(project,file="lzq_project")
z.pdry = z.scores(project, K = 2, d = 1)
z.pdry <- apply(z.pdry, 1, median)
lambda.pdry = median(z.pdry^2)/qchisq(0.5, df = 1)
p.pdry.adj = pchisq(z.pdry^2/lambda.pdry, df = 1, lower = FALSE)
lfmm.results <- cbind(z.pdry,p.pdry.adj)
write.table(lfmm.results, "lzq_result", sep="\t", quote=F, row.names=F)
write.table(lambda.pdry, "lzq_scale", sep="\t", quote=F, row.names=F)
