#set the path to the gsea gmt files (need to be downloaded first from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
gsea_dir <- ""

pathways.hallmark <- gmtPathways(file.path(gsea_dir, "h.all.v7.1.symbols.gmt"))
pathways.bio_processes <- gmtPathways(file.path(gsea_dir, "c5.bp.v7.1.symbols.gmt"))
pathways.mol_funs <- gmtPathways(file.path(gsea_dir, "c5.mf.v7.1.symbols.gmt"))
pathways.cell_comp <- gmtPathways(file.path(gsea_dir, "c5.cc.v7.1.symbols.gmt"))
pathways.kegg <- gmtPathways(file.path(gsea_dir, "c2.cp.kegg.v7.1.symbols.gmt"))

