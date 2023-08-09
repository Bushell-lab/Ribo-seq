#set the path to the gsea gmt files (need to be downloaded first from https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp)
gsea_dir <- ""

pathways.hallmark <- gmtPathways(file.path(gsea_dir, "mh.all.v2022.1.Mm.symbols.gmt"))
pathways.bio_processes <- gmtPathways(file.path(gsea_dir, "m5.go.bp.v2022.1.Mm.symbols.gmt"))
pathways.mol_funs <- gmtPathways(file.path(gsea_dir, "m5.go.mf.v2022.1.Mm.symbols.gmt"))
pathways.cell_comp <- gmtPathways(file.path(gsea_dir, "m5.go.cc.v2022.1.Mm.symbols.gmt"))
pathways.curated <- gmtPathways(file.path(gsea_dir, "m2.all.v2022.1.Mm.symbols.gmt"))
