#Check this script is in the same directory as the .gmt files below, or add path accordingly (need to be downloaded first from https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp)

pathways.hallmark <- gmtPathways("h.all.v7.1.symbols.gmt")
pathways.bio_processes <- gmtPathways("c5.bp.v7.1.symbols.gmt")
pathways.mol_funs <- gmtPathways("c5.mf.v7.1.symbols.gmt")
pathways.cell_comp <- gmtPathways("c5.cc.v7.1.symbols.gmt")
pathways.kegg <- gmtPathways("c2.cp.kegg.v7.1.symbols.gmt")



