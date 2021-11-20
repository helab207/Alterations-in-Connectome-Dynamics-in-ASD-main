library(clusterProfiler)
library(enrichplot)
library(ggplot2)

d = read.csv(file="PLS1_OrdGeneList.csv",header=TRUE)
geneList=d[,2]
names(geneList)=as.character(d[,1])
geneList = sort(geneList,decreasing = TRUE)
geneSet = clusterProfiler::read.gmt("CustomGeneSets.gmt")
res<-GSEA(geneList,exponent = 1,
          nPerm = 10000,
          minGSSize = 10,
          maxGSSize = 750,
          pvalueCutoff = 1,
          pAdjustMethod = "BH",
          TERM2GENE=geneSet,
          TERM2NAME = NA,
          verbose = TRUE,
          seed = FALSE,
          by = "fgsea"
)
summary(res)
GeneSetOrder<-res$ID
NES_value<-res$NES
Padjust<-res$p.adjust
Results<-data.frame(GeneSetOrder,NES_value,Padjust)

save(res, file="res.RData")
save(Results,file="Results.RData")
write.csv(Results, file = "Results.csv",row.names = FALSE)


######################################################################

# visualizing the results of GSEA

#load("res.RData")

gseaplot2(res, 
         geneSetID="Parikshak_Down",
         title="", 
         color = "#2F5C85",
         base_size = 11,
         rel_heights = c(1.5, 0.5, 1),
         subplots = 1:3,
         pvalue_table = FALSE,
         ES_geom = "line")
ggsave('Parikshak_Down.tiff', width =20,height = 14.5,dpi = 300, limitsize = FALSE,units = "cm")

#####
gseaplot2(res, 
          geneSetID="krishnanGenelevel_all",
          title="", 
          color = "#2F5C85",
          base_size = 11,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")
ggsave('krishnanGenelevel_all.tiff', width =20,height = 14.5,dpi = 300, limitsize = FALSE,units = "cm")

#####
gseaplot2(res, 
          geneSetID="Parikshak_Up",
          title="", 
          color = "#2F5C85",
          base_size = 11,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")
ggsave('Parikshak_Up.tiff', width =20,height = 14.5,dpi = 300, limitsize = FALSE,units = "cm")

#####
gseaplot2(res, 
          geneSetID="GWAS",
          title="", 
          color = "#2F5C85",
          base_size = 11,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")
ggsave('GWAS.tiff', width =20,height = 14.5,dpi = 300, limitsize = FALSE,units = "cm")

#####
gseaplot2(res, 
          geneSetID="Sanders",
          title="", 
          color = "#2F5C85",
          base_size = 11,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")
ggsave('Sanders.tiff', width =20,height = 14.5,dpi = 300, limitsize = FALSE,units = "cm")

#####
gseaplot2(res, 
          geneSetID="Satterstrom",
          title="", 
          color = "#2F5C85",
          base_size = 11,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")
ggsave('Satterstrom.tiff', width =20,height = 14.5,dpi = 300, limitsize = FALSE,units = "cm")

#####
gseaplot2(res, 
          geneSetID="krishnanNegtive",
          title="", 
          color = "#2F5C85",
          base_size = 11,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")
ggsave('krishnanNegtive.tiff', width =20,height = 14.5,dpi = 300, limitsize = FALSE,units = "cm")








