import pandas as pd

gene_info = pd.read_csv(
    "/data/chanjlab/Xenium.lung.NE_plasticity.010124/ref/Final_Gene_Panel_Design.081523/Xenium_hLung_v1_metadata.csv"
)
gene_info_2 = pd.read_table(
    "/data/chanjlab/Xenium.lung.NE_plasticity.010124/ref/Final_Gene_Panel_Design.081523/Xenium.custom_gene_panel.v5.txt"
)
genes = pd.concat([gene_info.Gene, gene_info_2.Gene])
annotations = pd.concat([gene_info.Annotation, gene_info_2.Info])
out = pd.DataFrame({"gene": genes, "annotation": annotations})
out.to_csv("../data/gene_information.csv")
