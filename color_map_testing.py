import pandas as pd
from pycirclize import Circos
import matplotlib.pyplot as plt
from ncbi.datasets.openapi.model.v1_orientation import V1Orientation
from sklearn import preprocessing
import numpy as np

gene_exp_data = pd.read_excel('C:/Users/cheif/Downloads/Expression data CR10.xlsx', header=1)

walnut_gene_meta = pd.read_csv('https://raw.githubusercontent.com/nitink23/WGKB/main/ncbi_dataset.tsv', delimiter='\t')

full_data = pd.merge(gene_exp_data, walnut_gene_meta, on='Gene ID', how='inner')

data = pd.DataFrame({
    'Chromosome': ['chr01', 'chr02', 'chr03', 'chr04', 'chr05', 'chr06', 'chr07', 'chr08', 'chr09', 'chr10', 
                   'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'Pltd'],  
    'Size (bp)': [45207397, 37821870, 35064427, 34823025, 22562875, 39020271, 52418484, 30564197, 24263475,
                  37707155, 37114715, 31492331, 39757759, 28841373, 20407330, 28711772, 160537]
})

sectors = {str(row['Chromosome']): row['Size (bp)'] for index, row in data.iterrows()}
circos = Circos(sectors, space=5)

for index, sector_obj in enumerate(circos.sectors):

    # Plot sector name
    sector_obj.text(f"{sector_obj.name}", r=110, size=15)

    # Add scatter plot points based on genomic ranges and orientation
    scatter_track = sector_obj.add_track((95, 100), r_pad_ratio=0.1)
    scatter_track.axis()
    
    # Scatter values in given chromosome
    chr_data = full_data[full_data['Chromosome'] == str(index+1)].reset_index(drop=True)
    
    # Normalize column and get color map
    vectorized_color_col = np.array(chr_data['log2FC CR10 pini / mock']).reshape(-1, 1)

    if vectorized_color_col.shape[0] == 0:
        continue
    mean = vectorized_color_col.mean()
    sd = vectorized_color_col.std()

    normalized_arr = (vectorized_color_col - mean) / sd
    cmap = plt.get_cmap('coolwarm')
    colors = [cmap(value) for value in chr_data['log2FC CR10 pini / mock']]

    for row_num, row in chr_data.iterrows():
        y = 1 if row["Orientation"] == 'plus' else 0
        scatter_track.scatter([row['Begin']], [y], color=colors[row_num])
        
fig = circos.plotfig()
plt.show()
