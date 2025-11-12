#install.packages(c("arrow","dplyr"))
library(arrow)
library(dplyr) 

caminho_base <- file.path("3W", "dataset")

caminhos_dos_arquivos <- list.files(
  path = caminho_base,
  pattern = "\\.parquet$",
  recursive = TRUE, 
  full.names = TRUE 
)

dataset_dhsv_united <- open_dataset(
  caminhos_dos_arquivos, 
  format = "parquet"
)

df_dhsv_united <- dataset_dhsv_united %>%
  collect()

print(paste("N° de linhas: ", nrow(df_dhsv_united)))
print(paste("N° de colunas",ncol(df_dhsv_united)))
#print(colnames(df_sample)) #nome das colunas
print(head(df_dhsv_united))