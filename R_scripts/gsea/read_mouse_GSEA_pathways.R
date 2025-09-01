#Check this script is in the same directory as the .gmt files below, or add path accordingly (need to be downloaded first from https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp)

pathways.hallmark <- gmtPathways("mh.all.v2022.1.Mm.symbols.gmt")
pathways.bio_processes <- gmtPathways("m5.go.bp.v2022.1.Mm.symbols.gmt")
pathways.mol_funs <- gmtPathways("m5.go.mf.v2022.1.Mm.symbols.gmt")
pathways.cell_comp <- gmtPathways("m5.go.cc.v2022.1.Mm.symbols.gmt")
pathways.curated <- gmtPathways("m2.all.v2022.1.Mm.symbols.gmt")

