import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import seaborn as sns
import os

# Pasta onde os graficos (png e pdf) serao salvos
pasta_saida = "C:/Users/Viviane Jotta/Desktop/UFMG/Doutorado/Projeto - painel solar/Metodologia/MAGs/DRAM/DRAM_spades"

##################################################################
# PASSO 0 - CRIAR DICIONARIO COM CATEGORIAS E KO_IDs DE INTERESSE
##################################################################

print('######################################################')
print('PASSO 0 - DICIONARIO CATEGORIAS')
print('######################################################')

ko_list = {
 "Transporte_Fe": {
   "K02014": "TonB-dep. receptor",
   "K04758": "feoA",
   "K04759": "feoB",
   "K02016": "ABC.FEV.S",
   "K02012": "afuA/fbpA",
   "K22736": "VIT",
   "K25283": "yclO/ceuC",
   "K09826": "irr",
 },
 "Biofilme_EPS": {
   "K21023": "DGC",
   "K01768": "adenilato ciclase",
   "K11935": "pgaA",
   "K11936": "pgaC",
   "K11937": "pgaD",
   "K21006": "pelA",
   "K20997": "pslA",
   "K04334": "csgA",
   "K04333": "csgD",
   "K19336": "bdcA",
   "K23778": "rcdA",
   "K19334": "tabA",
   "K19449": "XRE",
   "K13654": "GntR",
   "K21976": "ariR",
   "K19688": "bssR",
   "K13650": "mqsR",
   "K25884": "mucR",
   "K16566": "exoY",
   "K16567": "exoQ",
   "K16568": "exoZ",
   "K25885": "exoR",
 },
 "Quorum_sensing": {
   "K07173": "luxS",
   "K22954": "luxI",
   "K18098": "rpaR",
   "K13060": "lasI",
   "K18304": "lasR",
   "K18099": "rhlR",
 },
 "Quimiotaxia_motilidade": {
   "K03406": "mcp",
   "K03407": "cheA",
   "K03408": "cheW",
   "K03412": "cheB",
   "K00575": "cheR",
   "K13924": "cheB/cheR",
   "K03414": "cheZ",
   "K02556": "motA",
   "K02557": "motB",
 },
 "Estresse_oxidativo": {
   "K03781": "katE",
   "K03782": "katG",
   "K04564": "sodA/sodB",
   "K04565": "sodC",
   "K04761": "oxyR",
 },
 "Choque_termico": {
   "K04077": "GroEL",
   "K04078": "GroES",
   "K04043": "dnaK",
   "K03687": "grpE",
   "K03686": "dnaJ",
   "K05801": "DnaJ-like",
 },
 "Reparo_UV": {
   "K03701": "uvrA",
   "K03702": "uvrB",
   "K03703": "uvrC",
   "K03657": "uvrD",
 },
 "Carotenoides_pigmentos": {
   "K02291": "crtB",
   "K10027": "crtI",
   "K06443": "crtY/lcyB",
   "K22502": "crtY",
   "K15746": "crtZ",
   "K09836": "crtW",
   "K02292": "crtO",
 },
 "Regulacao_sinalizacao": {
   "K03088": "sigma-70 ECF",
   "K20374": "sistema 2 componentes",
 },
}

ko_categoria = {}
ko_gene = {}

for i, j in ko_list.items():
  for a,b in j.items():
    if a in ko_categoria:
      print("duplicado:", a)
    ko_categoria[a] = i
    ko_gene[a] = b

#Confere se tem KO repetido
tamanho = [len(v) for v in ko_list.values()] # cria uma lista com todos os tamanhos dos valores

if len(ko_categoria) == sum(tamanho):
  print("Sem KO repetido")
else:
  print("KO repetido")

##################################################################
# PASSO 1 - LER ARQUIVO, CONFERIR E FILTRAR
##################################################################

print('######################################################')
print('PASSO 1 - LER ARQUIVO, CONFERIR E FILTRAR')
print('######################################################')

#Lendo o arquivo annotations.tsv
df = pd.read_csv("C:/Users/Viviane Jotta/Desktop/UFMG/Doutorado/Projeto - painel solar/Metodologia/MAGs/DRAM/DRAM_spades/annotations.tsv", sep="\t", usecols=["fasta", "ko_id"])

# Conferência
print(f"O número de linhas e colunas é: {df.shape}")
print("As primeiras 10 linhas:")
print(df.head(10))
print(f"O numero de NaN em ko_id: {df['ko_id'].isna().sum()}")

print(f"O numero de hits com ko_id:{df.shape[0]-df['ko_id'].isna().sum()}")
df["ko_id"].str.contains(",", na=False).sum() #Confere se tem mais de um KO na mesma celula

# Filtrando df para aparecer apenas quem tem ko_id
df_2 = df[df["ko_id"].notna()]

#########################################################
# PASSO 2 - CONFERIR E FILTRAR KO DE INTERESSE
#########################################################

print('######################################################')
print('PASSO 2 - FILTRAR KO DE INTERESSE')
print('######################################################')

# Conferir se tem mais de um ko na coluna ko_id
df_2[df_2["ko_id"].str.contains(",")]
#Acrescentar no codigo depois o que acontece se tiver

# Filtra ko de interesse
df_filtrado = df_2[df_2["ko_id"].isin(ko_gene.keys())]
df_filtrado

# Conta quantas ocorrencias tem daquele gene
df_filtrado["ko_id"].value_counts() # Serve mais pra validação

# vamos fazer uma segunda tabela de frequencias
tabela = pd.crosstab(df_filtrado["fasta"], df_filtrado["ko_id"])
print(f'O formato da tabela é: {tabela.shape} linhas x colunas')
print(tabela.head(10))

##############################################################
# PASSO 3 - COLOCAR GENEROS NA TABELA
##############################################################

print('######################################################')
print('PASSO 3 - COLOCAR GENEROS NA TABELA')
print('######################################################')

#Lendo o arquivo master_table.tsv
genus_table = pd.read_csv("C:/Users/Viviane Jotta/Desktop/UFMG/Doutorado/Projeto - painel solar/Metodologia/MAGs/GTDB-tk/GTDB-Tk - spades/master_table_168MAGs.tsv", sep="\t", usecols=["MAG_ID", "genus"])
print(f"O número de linhas e colunas é: {genus_table.shape}")

# Conferindo se existe duplicata na master_table
if genus_table["MAG_ID"].duplicated().any() == False:
  print('Não existem IDs de MAGs duplicados')
else:
  print('Existem MAGs IDs duplicados:')
  print(genus_table["MAG_ID"][genus_table["MAG_ID"].duplicated()])

#Comparando se os Nomes das MAGs batem entre as colunas

ids_genus = set(genus_table["MAG_ID"])
ids_tabela = set(tabela.index)

print(f"IDs na genus_table: {len(ids_genus)}")
print(f"IDs na tabela (crosstab): {len(ids_tabela)}")

if ids_genus == ids_tabela:
    print("Conferência concluída: Todos os IDs batem")
else:
    print("Conferência concluída: Alguns IDs NÃO batem")

    so_na_genus  = ids_genus  - ids_tabela
    so_na_tabela = ids_tabela - ids_genus

    if so_na_genus:
        print(f"Só na genus_table ({len(so_na_genus)}): {sorted(so_na_genus)}")
    if so_na_tabela:
        print(f"Só na tabela/crosstab ({len(so_na_tabela)}): {sorted(so_na_tabela)}")

# Colocando o genero na tabela
tabela_final = pd.merge(tabela, genus_table, left_on=tabela.index, right_on="MAG_ID")
if tabela_final['genus'].isna().sum() == 0: #Conferir de tem NaN
  print("Nenhum NaN na tabela final")
else:
  print("Existe um NaN na tabela final")

# Transformando o MAG_id em index novamente
tabela_final = tabela_final.set_index("MAG_ID")
#print(tabela_final)

#Filtrando ko_genes para aparecer apenas o que existe na tabela final
ko_gene_set = set(ko_gene.keys())
tabela_colunas_set = set(tabela_final.columns)

ko_faltante = ko_gene_set - tabela_colunas_set
print(f"KO que nao estao nos dados: {ko_faltante}")

ko_list_2 = {}
for categoria,ko_da_categoria in ko_list.items():
  dicionario = {}
  for ko,gene in ko_da_categoria.items():
    if ko not in ko_faltante:
      dicionario[ko] = gene
  ko_list_2[categoria] = dicionario

soma_ko_list_2 = 0
for i,j in ko_list_2.items():
  soma_ko_list_2 += len(j)
soma_ko_list = 0
for i,j in ko_list.items():
  soma_ko_list += len(j)

assert soma_ko_list - soma_ko_list_2 == len(ko_faltante), (
    f"Conferência falhou: esperava remover {len(ko_faltante)} KOs, "
    f"mas a diferença foi {soma_ko_list - soma_ko_list_2}. "
    "Verifique o laço que cria o ko_list_2."
)
print('KOs vazios removidos, se existirem')
# Atualizando ko_categoria e ko_gene para os proximos passos
ko_categoria = {}
ko_gene = {}

for categoria, ko_da_categoria in ko_list_2.items():
  for ko, gene in ko_da_categoria.items():
    if ko in ko_categoria:
      print("duplicado:", a)
    ko_categoria[ko] = categoria
    ko_gene[ko] = gene

#Confere se tem KO repetido
tamanho = [len(v) for v in ko_list_2.values()] # cria uma lista com todos os tamanhos dos valores

if len(ko_categoria) == sum(tamanho):
  print("Sem KO repetido")
else:
  print("KO repetido")

#############################################################
# PASSO 4 - FAZENDO TABELA AGREGADA DE GENERO
##############################################################

print('######################################################')
print('PASSO 4 - TABELA AGREGADA DE GENERO')
print('######################################################')

print(tabela_final["genus"].nunique()) #Quantos generos tem no total
tabela_final["genus"].value_counts() # lista de frequencia de cada genero

#Fazendo uma tabela agregada para a imagem fica menor e conseguir colocar no artigo/dissertacao
tabela_final_agregada = tabela_final.groupby("genus")[list(ko_gene.keys())].mean()
#Essa linha esta agrupando por genero as linhas e para uma lista de KO (colunas) esta calculando a média da copia de genes
print(tabela_final_agregada.shape)

#Modificando a tabela_final pra colocar o nome dos Generos junto com as MAG ids
#Criando uma nova coluna para a tabela final
tabela_final["MAG_genus"] = tabela_final.index + "_" + tabela_final["genus"]

#Ordenando as linhas para generos iguais ficares proximos
tabela_final = tabela_final.sort_values(by='genus')

# Remover a coluna genus
tabela_final = tabela_final.drop(columns="genus")

# Colocando a nova coluna como indice
tabela_final = tabela_final.set_index("MAG_genus")

print("As primeiras 10 linhas:")
print(tabela_final.head(10))

#############################################################
# PASSO 5 - PREPARANDO OS DADOS PARA PLOTAR GRÁFICO
##############################################################

print('######################################################')
print('PASSO 5 - PREPARANDO OS DADOS PARA PLOTAR GRAFICOS')
print('######################################################')

# montando uma lista ordenada de KOs por categoria

ordem_ko = []
for i in ko_list_2.keys():
  for j in ko_list_2[i]:
    ordem_ko.append(j)

#Montando uma lista com as categorias na mesma ordem

lista_categorias = []
for k in ordem_ko:
  lista_categorias.append(ko_categoria[k])

# Tranformar isso em series para o pandas e para fazer o heatmap
serie_categorias = pd.Series(lista_categorias, index = ordem_ko)

# Mapear categorias em cores escolhidas
mapa_cores = {"Transporte_Fe":"#564138",
              "Biofilme_EPS":"#2E86AB",
              "Quorum_sensing":"#CB9CF2",
              "Quimiotaxia_motilidade": "#F5F749",
              "Estresse_oxidativo": "#F24236",
              "Choque_termico": "#0C7C59",
              "Reparo_UV": "#6200B3",
              "Carotenoides_pigmentos": "#FF8C42",
              "Regulacao_sinalizacao": "#8EA604"}

serie_cores = serie_categorias.map(mapa_cores)

# Reordenar as colunas da tabela agregada para aparecer na ordem da ordem_ko
tabela_agregada_plot = tabela_final_agregada[ordem_ko]

##############################################################
# PASSO 6 - PLOTAGEM DE GRPAFICOS
##############################################################

print('######################################################')
print('PASSO 6 - PLOTANDO OS GRAFICOS')
print('######################################################')

print('Grafico 1. Heatmap agregado por gênero (Numero de copias)')


# 1. Desenha o heatmap
g1 = sns.clustermap(tabela_agregada_plot,
                    col_colors = serie_cores,
                    cmap="YlOrRd",
                    vmin=1, vmax=52,
                    annot=True,
                    fmt=".1f",
                    row_cluster=False,
                    col_cluster=False,
                    mask=(tabela_agregada_plot == 0),
                    figsize=(25, 12),
                    cbar_pos=None)

# 2. Cria a legenda das categorias (aqueles quadradinhos coloridos com o nome do lado)
handles = [Patch(color=cor, label=cat) for cat, cor in mapa_cores.items()]

# Adiciona a legenda no gráfico.
g1.ax_heatmap.legend(handles=handles,
                     title="Categorias",
                     bbox_to_anchor=(-0.30, 0.5),
                     loc="center left",
                     frameon=False,
                     fontsize=14)

# 3. Adiciona barra de cores (logo abaixo da legenda de categorias)
# As coordenadas são [esquerda, baixo, largura, altura] na figura inteira (0 a 1)
ax_cb = g1.fig.add_axes([0.02, 0.60, 0.015, 0.15]) # <-- MUDANÇA AQUI (ajuste a posição conforme a imagem)
ColorbarBase(ax_cb,
             cmap=plt.get_cmap("YlOrRd"),
             norm=Normalize(vmin=1, vmax=52),
             orientation="vertical",
             label="Média de cópias")
#Cola o grafico a esquerda ignorando o dendograma
g1.fig.subplots_adjust(left=0.0)

# Salvando o grafico 1
plt.rcParams["svg.fonttype"] = "none"   # mantém o texto editável no SVG
g1.savefig(os.path.join(pasta_saida, "heatmap_genero.png"), dpi=300, bbox_inches="tight")
g1.savefig(os.path.join(pasta_saida, "heatmap_genero.pdf"), bbox_inches="tight")

#plt.show()


# GRAFICO 2 - COMPLETO

print('Grafico 2. Heatmap completo (Numero de copias)')
# Reordenar as colunas da tabela agregada para aparecer na ordem da ordem_ko
tabela_final_plot = tabela_final[ordem_ko]
print(f'Shape da tabela final (MAGs vcs KO_id): {tabela_final_plot.shape}')

# 1. Desenha o heatmap
g2 = sns.clustermap(tabela_final_plot,
                    col_colors = serie_cores,
                    cmap="YlOrRd",
                    vmin=1, vmax=52,
                    annot=True,
                    fmt=".1f",
                    row_cluster=False,
                    col_cluster=False,
                    mask=(tabela_final_plot == 0),
                    figsize=(25, 35),
                    cbar_pos=None)

# Força o matplotlib a mostrar TODOS os nomes do eixo Y
g2.ax_heatmap.set_yticks(range(len(tabela_final_plot.index)))
g2.ax_heatmap.set_yticklabels(tabela_final_plot.index, rotation=0, fontsize=8)

# 2. Cria a legenda das categorias (aqueles quadradinhos coloridos com o nome do lado)
handles = [Patch(color=cor, label=cat) for cat, cor in mapa_cores.items()]

# Adiciona a legenda no gráfico.
g2.ax_heatmap.legend(handles=handles,
                     title="Categorias",
                     bbox_to_anchor=(-0.30, 0.5),
                     loc="center left",
                     frameon=False,
                     fontsize=14)

# 3. Adiciona barra de cores (logo abaixo da legenda de categorias)
# As coordenadas são [esquerda, baixo, largura, altura] na figura inteira (0 a 1)
ax_cb = g2.fig.add_axes([0.02, 0.60, 0.015, 0.15]) # <-- MUDANÇA AQUI (ajuste a posição conforme a imagem)
ColorbarBase(ax_cb,
             cmap=plt.get_cmap("YlOrRd"),
             norm=Normalize(vmin=1, vmax=52),
             orientation="vertical",
             label="Copia no genoma")
#Cola o grafico a esquerda ignorando o dendograma
g2.fig.subplots_adjust(left=0.0)

# Salvando o grafico 1
plt.rcParams["svg.fonttype"] = "none"   # mantém o texto editável no SVG
g2.savefig(os.path.join(pasta_saida, "heatmap_completo.png"), dpi=300, bbox_inches="tight")
g2.savefig(os.path.join(pasta_saida, "heatmap_completo.pdf"), bbox_inches="tight")

#plt.show()

#############################################################
# Grafico presença/ausencia AGREGADO
#############################################################

print('Grafico 3. Heatmap agregado por gênero (Presenca/Ausencia)')

# Transformar em presença e ausencia
tabela_presenca_agregada = (tabela_agregada_plot > 0).astype(int)

# 1. Desenha o heatmap
g3 = sns.clustermap(tabela_presenca_agregada,
                    col_colors = serie_cores,
                    cmap="YlOrRd",
                    vmin=0, vmax=1,
                    row_cluster=False,
                    col_cluster=False,
                    figsize=(22, 12),
                    cbar_pos=None)

# Força o matplotlib a mostrar TODOS os nomes do eixo Y
g3.ax_heatmap.set_yticks(range(len(tabela_presenca_agregada.index)))
g3.ax_heatmap.set_yticklabels(tabela_presenca_agregada.index, rotation=0, fontsize=8)

# 2. Legenda das categorias
handles = [Patch(color=cor, label=cat) for cat, cor in mapa_cores.items()]
g3.fig.legend(handles=handles,
              title="Categorias",
              bbox_to_anchor=(0.02, 0.50), # Mais para cima
              loc="center left",
              frameon=False,
              fontsize=14,
              title_fontsize=16)

# 3. Legenda de presença ou ausência
cor_vinho = plt.get_cmap("YlOrRd")(1.0) 

handles_presenca = [Patch(facecolor="white", label="Ausente", edgecolor="black"),
                    Patch(color=cor_vinho, label="Presente")] 
g3.fig.legend(handles=handles_presenca, 
              title="Presença", 
              bbox_to_anchor=(0.02, 0.72), # Mais para baixo
              loc="center left", 
              frameon=False,
              fontsize=14,
              title_fontsize=16)

# Salvando o grafico 3
plt.rcParams["svg.fonttype"] = "none"   # mantém o texto editável no SVG
g3.savefig(os.path.join(pasta_saida, "heatmap_agregado_PA.png"), dpi=300, bbox_inches="tight")
g3.savefig(os.path.join(pasta_saida, "heatmap_agregado_PA.pdf"), bbox_inches="tight")

#plt.show()

#############################################################
# Grafico presença/ausencia COMPLETO
#############################################################

print('Grafico 4. Heatmap completo (Presenca/Ausencia)')

# Transformar em presença e ausencia
tabela_presenca_completa = (tabela_final_plot > 0).astype(int)

# 1. Desenha o heatmap
g4 = sns.clustermap(tabela_presenca_completa,
                    col_colors = serie_cores,
                    cmap="YlOrRd",
                    vmin=0, vmax=1,
                    row_cluster=False,
                    col_cluster=False,
                    #mask=(tabela_final_plot == 0),
                    figsize=(20, 32),
                    cbar_pos=None)

# Força o matplotlib a mostrar TODOS os nomes do eixo Y
g4.ax_heatmap.set_yticks(range(len(tabela_presenca_completa.index)))
g4.ax_heatmap.set_yticklabels(tabela_presenca_completa.index, rotation=0, fontsize=8)

# 2. Legenda das categorias
handles = [Patch(color=cor, label=cat) for cat, cor in mapa_cores.items()]
g4.fig.legend(handles=handles,
              title="Categorias",
              bbox_to_anchor=(0.02, 0.63), # Mais para cima
              loc="center left",
              frameon=False,
              fontsize=14,
              title_fontsize=16)

# 3. Legenda de presença ou ausência
cor_vinho = plt.get_cmap("YlOrRd")(1.0) 

handles_presenca = [Patch(facecolor="white", label="Ausente", edgecolor="black"),
                    Patch(color=cor_vinho, label="Presente")] 
g4.fig.legend(handles=handles_presenca, 
              title="Presença", 
              bbox_to_anchor=(0.02, 0.72), # Mais para baixo
              loc="center left", 
              frameon=False,
              fontsize=14,
              title_fontsize=16)

# Salvando o grafico 3
plt.rcParams["svg.fonttype"] = "none"   # mantém o texto editável no SVG
g4.savefig(os.path.join(pasta_saida, "heatmap_completo_PA.png"), dpi=300, bbox_inches="tight")
g4.savefig(os.path.join(pasta_saida, "heatmap_completo_PA.pdf"), bbox_inches="tight")

#plt.show()