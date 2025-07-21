Análise de Gene ID e Anotação:
  # Carregar pacotes necessários para anotação
  # Se você ainda não os tem, instale-os primeiro:
  # if (!requireNamespace("BiocManager", quietly = TRUE))
  #    install.packages("BiocManager")
  # BiocManager::install("AnnotationDbi")
  # BiocManager::install("org.Mm.eg.db") # Para Mus musculus
  library(AnnotationDbi)
# 1. Instalar BiocManager (se ainda não estiver instalado)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# 2. Instalar os pacotes do Bioconductor necessários
BiocManager::install("AnnotationDbi")
# Certifique-se de que os pacotes necessários estão carregados
library(AnnotationDbi)
library(org.Mm.eg.db) # Banco de dados de anotação para Mus musculus
setwd("D:/bioinformatica/GEO/results.GSE145840")
# Carregar pacotes necessários
library(AnnotationDbi)
library(org.Mm.eg.db) # Banco de dados de anotação para Mus musculus
library(writexl) # Para salvar em formato .xlsx
print("Iniciando a anotação dos DEGs com Símbolo e Nome do Gene (CORREÇÃO FINAL do keytype) e salvando em .xlsx...")
# Carregar os resultados DESeq2 (use esta linha se você começou uma nova sessão do R)
deseq_results_per_tissue <- readRDS("D:/bioinformatica/deseq2_results_per_tissue_GSE145840.rds")
# Criar uma nova lista para armazenar os resultados anotados
deseq_results_annotated_per_tissue <- list()
for (current_tissue in names(deseq_results_per_tissue)) {
  print(paste("Anotando resultados para o tecido:", current_tissue))
  # Acessar os resultados brutos do DESeq2 para o tecido atual
  res_df <- deseq_results_per_tissue[[current_tissue]]
  # Os GeneIDs são os rownames; precisamos convertê-los para uma coluna para o mapeamento
  res_df$GeneID <- rownames(res_df)
  # Mapear GeneIDs (que são RIKEN IDs ou aliases) para ENTREZID
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = res_df$GeneID,
                       column = "ENTREZID",
                       keytype = "ALIAS", # *** AQUI ESTÁ A CORREÇÃO FINAL ***
                       multiVals = "first")
  # Adicionar os ENTREZID ao dataframe
  res_df$ENTREZID <- entrez_ids
  # Mapear ENTREZID para SYMBOL e GENENAME
  symbols <- mapIds(org.Mm.eg.db,
                    keys = res_df$ENTREZID,
                    column = "SYMBOL",
                    keytype = "ENTREZID",
                    multiVals = "first")
  gene_names <- mapIds(org.Mm.eg.db,
                       keys = res_df$ENTREZID,
                       column = "GENENAME",
                       keytype = "ENTREZID",
                       multiVals = "first")
  # Adicionar as novas colunas ao dataframe
  res_df$Symbol <- symbols
  res_df$GeneName <- gene_names
  # Reorganizar as colunas para ter GeneID, Symbol, GeneName no início
  res_df <- res_df[, c("GeneID", "Symbol", "GeneName", setdiff(names(res_df), c("GeneID", "Symbol", "GeneName", "ENTREZID")))]
  # Armazenar o dataframe anotado na nova lista
  deseq_results_annotated_per_tissue[[current_tissue]] <- res_df
  print(paste("   Anotação concluída para", current_tissue, ". Total de genes anotados:", sum(!is.na(res_df$Symbol))))
}
print("Anotação de todos os resultados DESeq2 concluída!")
# Visualizar o cabeçalho dos resultados anotados para um tecido (ex: LIVER)
if ("LIVER" %in% names(deseq_results_annotated_per_tissue)) {
  print("Cabeçalho dos resultados anotados para LIVER:")
  print(head(deseq_results_annotated_per_tissue$LIVER))
}
# --- SALVANDO EM ARQUIVO .XLSX COM MÚLTIPLAS ABAS ---
# O writexl::write_xlsx recebe uma lista de dataframes e salva cada um em uma aba.
output_excel_path_xlsx <- file.path(getwd(), results_folder, "DEGs_all_tissue_completo.xlsx")
