import streamlit as st
import pandas as pd
from pycirclize import Circos
import matplotlib.pyplot as plt
from ncbi.datasets import GeneApi
from ncbi.datasets.openapi import ApiClient
from ncbi.datasets.openapi.rest import ApiException
from ncbi.datasets.openapi.model.v1_orientation import V1Orientation
import numpy as np


def main() -> None:
    data = pd.DataFrame({
    'Chromosome': ['chr01', 'chr02', 'chr03', 'chr04', 'chr05', 'chr06', 'chr07', 'chr08', 'chr09', 'chr10', 
                   'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'Pltd'],  
    'Size (bp)': [45207397, 37821870, 35064427, 34823025, 22562875, 39020271, 52418484, 30564197, 24263475,
                  37707155, 37114715, 31492331, 39757759, 28841373, 20407330, 28711772, 160537]})
    gene_colors = {}
    chrom_dict = {
        'NC_049901.1': 'chr01',
        'NC_049902.1': 'chr02',
        'NC_049903.1': 'chr03',
        'NC_049904.1': 'chr04',
        'NC_049905.1': 'chr05',
        'NC_049906.1': 'chr06',
        'NC_049907.1': 'chr07',
        'NC_049908.1': 'chr08',
        'NC_049909.1': 'chr09',
        'NC_049910.1': 'chr10',
        'NC_049911.1': 'chr11',
        'NC_049912.1': 'chr12',
        'NC_049913.1': 'chr13',
        'NC_049914.1': 'chr14',
        'NC_049915.1': 'chr15',
        'NC_049916.1': 'chr16',
        'NC_028617.1': 'Pltd'}

    st.title('Custom Circos Plot Generator with Gene Metadata Integration')

    # Allow user to upload optional gene expression file
    gene_exp_file = st.file_uploader('Upload a gene expression file (optional, must be .csv, .xls or .xlsx)', type=["csv", "xls", "xlsx"])

    # Read walnut gene metadata file straight from GitHub
    url = 'https://raw.githubusercontent.com/nitink23/WGKB/main/ncbi_dataset.tsv'
    walnut_gene_meta = pd.read_csv(url, delimiter='\t')

    all_genes = None
    full_data = None
    gene_exp_gene_ids, gene_exp_avg_exp, gene_exp_pini_mock, gene_exp_capsici_mock, gene_exp_pini_capsici = None, None, None, None, None

    
    if gene_exp_file:
        try:
            gene_exp_df, gene_exp_gene_ids, gene_exp_avg_exp, gene_exp_pini_mock, gene_exp_capsici_mock, gene_exp_pini_capsici = read_gene_exp_file(gene_exp_file)
            full_data = pd.merge(gene_exp_df, walnut_gene_meta, on='Gene ID', how='inner')

            all_genes = st.checkbox('Display genes in gene expression file')

            # Check to see if there are string values in each of the columns provided
            for val_1, val_2, val_3, val_4, val_5 in zip(full_data[gene_exp_gene_ids], full_data[gene_exp_avg_exp], full_data[gene_exp_pini_mock], \
                                                         full_data[gene_exp_capsici_mock], full_data[gene_exp_pini_capsici]):
                values = (val_1, val_2, val_3, val_4, val_5)

                if all(isinstance(value, str) for value in values):

                    st.warning('At least one of the columns selected contains string values. The app can only read integers.')
                    return
        except KeyError:
            st.warning('WARNING: The columns entered cannot be parsed correctly.')

    # User input for multiple gene IDs
    gene_id_input = st.text_input('Enter Gene IDs (space-separated)')
    genomic_ranges = []

    if gene_id_input:
        
        gene_ids = [gene_id.strip() for gene_id in gene_id_input.split() if gene_id.strip().isdigit()]
        
        # Ensure gene_ids are in correct format
        for gene_id in gene_ids:
            for char in gene_id:
                if not char.isdigit():
                    st.error(f'ERROR: {gene_id} is not a valid gene id')
                    break

        # Allow user to select a color for each gene ID
        for index, gene_id in enumerate(gene_ids):

            # Assign a unique key to each color picker using the gene ID and index
            color = st.color_picker(f"Pick a color for Gene ID {gene_id}", '#ff0000', key=f"color_picker_{gene_id}_{index}")
            gene_colors[int(gene_id)] = color

        gene_metadata = fetch_gene_metadata(gene_ids)
        genomic_ranges = display_gene_meta(gene_metadata, gene_ids)

        # Check to see if genomic ranges are valid, if invalid print warning and stop function
        if len(genomic_ranges) <= 0:
            gene_id_input = []
            genomic_ranges = []
            st.warning(f"No valid genomic ranges found for Gene IDs: {', '.join(gene_ids)}")
            return

    if all_genes or gene_id_input:
        
        try:
            display_circos_plot(data, genomic_ranges, chrom_dict, gene_colors, full_data, all_genes, \
                                gene_exp_gene_ids, gene_exp_avg_exp, gene_exp_pini_mock, gene_exp_capsici_mock, gene_exp_pini_capsici)
        except KeyError:
            st.warning('WARNING: The columns entered cannot be parsed correctly.')


def fetch_gene_metadata(gene_ids):
    # Fetch gene metadata for all gene IDs
    api_client = ApiClient()  # Properly initialize the API client
    gene_api = GeneApi(api_client)
    try:
        gene_metadata = gene_api.gene_metadata_by_id([int(gene_id) for gene_id in gene_ids])
        return gene_metadata
    except ApiException as e:
        st.error(f"An error occurred while fetching gene metadata: {e}")
        return None


def display_formatted_response(response):
    # Function to format and display the full API response as a useful table
    if response and 'genes' in response:
        gene_info_list = []

        # Iterate through each gene in the response
        for gene_data in response['genes']:
            gene_info = gene_data.get('gene', {})
            annotations = gene_info.get('annotations', [])
            genomic_ranges = gene_info.get('genomic_ranges', [])
            chromosomes = ', '.join(gene_info.get('chromosomes', []))
            symbol = gene_info.get('symbol', 'N/A')
            gene_id = gene_info.get('gene_id', 'N/A')
            description = gene_info.get('description', 'N/A')
            taxname = gene_info.get('taxname', 'N/A')
            common_name = gene_info.get('common_name', 'N/A')

            for annotation in annotations:
                for genomic_range in genomic_ranges:
                    for loc in genomic_range.get('range', []):
                        accession_version = genomic_range.get('accession_version', 'N/A')
                        begin = loc.get('begin', 'N/A')
                        end = loc.get('end', 'N/A')
                        orientation = loc.get('orientation', 'N/A')

                        # Append gene metadata to the list
                        gene_info_list.append({
                            'Gene ID': gene_id,
                            'Symbol': symbol,
                            'Description': description,
                            'Taxname': taxname,
                            'Common Name': common_name,
                            'Chromosomes': chromosomes,
                            'Accession Version': accession_version,
                            'Range Start': begin,
                            'Range End': end,
                            'Orientation': orientation,
                            'Release Date': annotation.get('release_date', 'N/A'),
                            'Release Name': annotation.get('release_name', 'N/A')
                        })

        # Convert the list of gene info to a DataFrame
        gene_df = pd.DataFrame(gene_info_list)

        # Adjust the index to start from 1 instead of the default 0
        gene_df.index = range(1, len(gene_df) + 1)  # This sets the index starting from 1

        # Display the DataFrame with the adjusted index
        st.dataframe(gene_df)

    else:
        st.warning("No gene information available in the API response.")


def read_gene_exp_file(gene_exp_file):
    if str(gene_exp_file).endswith('.csv'):
        gene_exp_df = pd.read_csv(gene_exp_file)
        csv = True
    else:
        gene_exp_df = pd.read_excel(gene_exp_file)
        csv = False

    # See if column column headers are in row 1 or row 2
    col_row = st.selectbox('Are the column headers in row 1 or row 2?', ['Row 1', 'Row 2'])

    if col_row == 'Row 1':
        header = 0
    else:
        header = 1
    
    # Re-read gene_exp_df with correct headers
    if csv:
        gene_exp_df = pd.read_csv(gene_exp_file, header = header)
    else:
        gene_exp_df = pd.read_excel(gene_exp_file, header = header)
    
    # Allow user to specify which columns denote which values
    gene_exp_gene_ids = st.selectbox("Select which column has the GeneIDs", gene_exp_df.columns)
    gene_exp_avg_exp = st.selectbox("Select which column has the average expression level", gene_exp_df.columns)
    gene_exp_pini_mock = st.selectbox("Select which column has the log2FC CR10 pini / mock", gene_exp_df.columns)
    gene_exp_capsici_mock = st.selectbox("Select which column has the log2FC CR10 capsici / mock", gene_exp_df.columns)
    gene_exp_pini_capsici = st.selectbox("Select which column has the log2FC CR10 pini / capsici", gene_exp_df.columns)
    
    # Rename col to 'Gene ID' and create merged dataframe
    gene_exp_df = gene_exp_df.rename(columns={gene_exp_gene_ids: 'Gene ID'})
    gene_exp_gene_ids = 'Gene ID'

    return gene_exp_df, gene_exp_gene_ids, gene_exp_avg_exp, gene_exp_pini_mock, gene_exp_capsici_mock, gene_exp_pini_capsici
    

def get_colormap(data_col):
    # input: column to plot that needs colors
    # output: colormap
    vectorized_color_col = np.array(data_col).reshape(-1, 1)

    # Default to see if column is empty
    if vectorized_color_col.shape[0] == 0:
        return ['#ff0000'] * len(data_col)
    
    min_val = vectorized_color_col.min()
    max_val = vectorized_color_col.max()
    
    normalized_arr = (vectorized_color_col - min_val) / (max_val - min_val)

    cmap = plt.get_cmap('coolwarm')
    colors = [cmap(value) for value in normalized_arr]
    return colors


def display_gene_meta(gene_metadata, gene_ids: list) -> list:
    # Log the entire API response as a table
    st.write("Full API Response:")
    display_formatted_response(gene_metadata.to_dict())  # Convert to dict for tabular format

    # Check if genes exist in the response
    if hasattr(gene_metadata, 'genes') and len(gene_metadata.genes) > 0:
        st.write(f"Fetched metadata for Gene IDs: {', '.join(gene_ids)}")

        # Parse the genomic ranges including orientation
        genomic_ranges = []
        for gene_data in gene_metadata.genes:
            gene_info = gene_data.gene

            if 'genomic_ranges' in gene_info and gene_info['genomic_ranges']:
                for genomic_range in gene_info['genomic_ranges']:
                    chrom = genomic_range['accession_version']
                    for loc in genomic_range['range']:
                        start = int(loc['begin'])
                        end = int(loc['end'])
                        orientation = loc.get('orientation', 'plus')  # Get 'orientation' from the correct location
                        genomic_ranges.append((chrom, start, end, orientation, gene_info['gene_id']))
            else:
                st.warning(f"No genomic ranges available for gene {gene_info['gene_id']}")
        return genomic_ranges
    else:
        st.warning(f'No gene ID metadata found for {gene_ids}')


def add_point(chrom: str, chrom_dict: dict, inputted_point: bool, sector_obj, orientation, scatter_track, scatter2_track, scatter3_track, \
              scatter4_track, bar_track, x: int, gene_colors, gene_id: int, full_data, gene_exp_gene_ids: str, gene_exp_pini_mock: str, \
              gene_exp_capsici_mock: str, gene_exp_pini_capsici: str, gene_exp_avg_exp: str, color_1, color_2, color_3, bar_color):
    
    if chrom in chrom_dict:
        chrom_name = chrom_dict[chrom]
    else:
        st.warning(f"Chromosome {chrom} not found in mapping. Skipping.")
        return  # Skip if chromosome is not found

    # Ensure the chromosome name matches the sector (e.g., 'Pltd', 'chr01')
    if chrom_name == sector_obj.name:

        # Extract the orientation value
        if isinstance(orientation, V1Orientation):
            orientation_value = orientation.value
        elif isinstance(orientation, str):
            orientation_value = orientation.strip().lower()
        else:
            st.warning(f"Unexpected type for orientation at start {x}. Skipping this entry.")
            return  
            
        # Assign y value based on orientation
        y = 1 if orientation_value == 'plus' else 0

        if full_data is not None:

            if str(gene_id) in full_data[gene_exp_gene_ids].astype(str).values:
                
                # Get the log2FC values from the user's gene expression file
                log2fc_PM = full_data[full_data[gene_exp_gene_ids] == int(gene_id)][gene_exp_pini_mock].values[0]
                log2fc_CM = full_data[full_data[gene_exp_gene_ids] == int(gene_id)][gene_exp_capsici_mock].values[0]
                log2fc_PC = full_data[full_data[gene_exp_gene_ids] == int(gene_id)][gene_exp_pini_capsici].values[0]
                avg_exp_lvl = full_data[full_data[gene_exp_gene_ids] == int(gene_id)][gene_exp_avg_exp].values[0]

                # Calculate vmin and vmax from the range of log2FC values in the expression file
                vmin_PM = full_data[gene_exp_pini_mock].min()
                vmax_PM = full_data[gene_exp_pini_mock].max()

                vmin_CM = full_data[gene_exp_capsici_mock].min()
                vmax_CM = full_data[gene_exp_capsici_mock].max()

                vmin_PC = full_data[gene_exp_pini_capsici].min()
                vmax_PC = full_data[gene_exp_pini_capsici].max()

                vmax_ael = full_data[gene_exp_avg_exp].max()

                # Plot point with correct color and size specifications
                if inputted_point:
                    scatter_track.scatter([x], [y], color=gene_colors[gene_id], s=70)
                    scatter2_track.scatter([x], [log2fc_PM], color=gene_colors[gene_id], vmin=vmin_PM, vmax=vmax_PM, s=70)
                    scatter3_track.scatter([x], [log2fc_CM], color=gene_colors[gene_id], vmin=vmin_CM, vmax=vmax_CM, s=70)
                    scatter4_track.scatter([x], [log2fc_PC], color=gene_colors[gene_id], vmin=vmin_PC, vmax=vmax_PC, s=70)
                    bar_track.bar([x], [avg_exp_lvl], ec=gene_colors[gene_id], lw=0.9, vmin=0, vmax=vmax_ael)
                else:
                    scatter2_track.scatter([x], [log2fc_PM], color=color_1, vmin=vmin_PM, vmax=vmax_PM)
                    scatter3_track.scatter([x], [log2fc_CM], color=color_2, vmin=vmin_CM, vmax=vmax_CM)
                    scatter4_track.scatter([x], [log2fc_PC], color=color_3, vmin=vmin_PC, vmax=vmax_PC)
                    bar_track.bar([x], [avg_exp_lvl], ec=bar_color, lw=0.9, vmin=0, vmax=vmax_ael)

            else:
                # Display a warning message if the gene ID is not found in the uploaded file
                st.warning(f"No metadata available for Gene ID {gene_id} in the uploaded file.")
        else:
            scatter_track.scatter([x], [y], color=gene_colors[gene_id], s=70)


def display_circos_plot(data: dict, genomic_ranges: list, chrom_dict: dict, gene_colors: dict, full_data, all_genes, \
                        gene_exp_gene_ids, gene_exp_avg_exp, gene_exp_pini_mock, gene_exp_capsici_mock, gene_exp_pini_capsici) -> None:
    # Prepare Circos plot with sectors (using chromosome sizes)
    sectors = {str(row[1]['Chromosome']): row[1]['Size (bp)'] for row in data.iterrows()}
    circos = Circos(sectors, space=5)

    color_1, color_2, color_3, bar_color = None, None, None, None

    if full_data is not None:
        colors_1 = get_colormap(full_data[gene_exp_pini_mock])
        colors_2 = get_colormap(full_data[gene_exp_capsici_mock])
        colors_3 = get_colormap(full_data[gene_exp_pini_capsici])

        bar_color = st.color_picker(f"Pick a color for bar plot on track 5", '#ff0000', key=f"color_picker")

    for index, sector_obj in enumerate(circos.sectors):
        # Plot sector name
        sector_obj.text(f"{sector_obj.name}", r=110, size=15)

        # Add scatter plot points based on genomic ranges and orientation
        scatter_track = sector_obj.add_track((95, 100), r_pad_ratio=0.1)
        scatter_track.axis()

        # Scatter track for log2FC CR10 pini / mock values from gene expression
        scatter2_track = sector_obj.add_track((75, 90), r_pad_ratio=0.1)
        scatter2_track.axis()

        # Scatter track for log2FC CR10 capsici / mock
        scatter3_track = sector_obj.add_track((55, 70), r_pad_ratio=0.1)
        scatter3_track.axis()

        # Scatter track for log2FC CR10 pini / capsici 
        scatter4_track = sector_obj.add_track((35, 50), r_pad_ratio=0.1)
        scatter4_track.axis()

        # Add bar track for Expression level
        bar_track = sector_obj.add_track((10, 30), r_pad_ratio=0.1)
        bar_track.axis()

        # Add y-ticks only on the left side of the chr01 sector
        if sector_obj.name == 'chr01':  

            scatter_track.yticks([0, 1], list("-+"), vmin=0, vmax=1, side="left")

            scatter2_track.yticks(y=np.linspace(full_data[gene_exp_pini_mock].min(), full_data[gene_exp_pini_mock].max(), num=5), \
                                  labels=[f"{round(tick)}" for tick in np.linspace(full_data[gene_exp_pini_mock].min(), full_data[gene_exp_pini_mock].max(), num=5)], \
                                    vmin=full_data[gene_exp_pini_mock].min(), vmax=full_data[gene_exp_pini_mock].max(), side="left" )
            
            scatter3_track.yticks(y=np.linspace(full_data[gene_exp_capsici_mock].min(), full_data[gene_exp_capsici_mock].max(), num=5), \
                                  labels=[f"{round(tick)}" for tick in np.linspace(full_data[gene_exp_capsici_mock].min(), full_data[gene_exp_capsici_mock].max(), num=5)], \
                                     vmin=full_data[gene_exp_capsici_mock].min(), vmax=full_data[gene_exp_capsici_mock].max(), side="left" )
            
            scatter4_track.yticks(y=np.linspace(full_data[gene_exp_pini_capsici].min(), full_data[gene_exp_pini_capsici].max(), num=5), \
                                  labels=[f"{round(tick)}" for tick in np.linspace(full_data[gene_exp_pini_capsici].min(), full_data[gene_exp_pini_capsici].max(), num=5)], \
                                     vmin=full_data[gene_exp_pini_capsici].min(), vmax=full_data[gene_exp_pini_capsici].max(), side="left")
        
        inputted_point = False

        # If use indicated they want all genes in gene expression file to be included
        if all_genes:

            # Get subset from full_data of all the rows with a specific chromosome
            chr_data = full_data[full_data['Chromosome'] == str(index+1)].reset_index(drop=True)

            # For each row in the chromosome subset, example: all chromosome 15 rows
            for row_num, row in chr_data.iterrows():

                x = row['Begin']
                color_1 = colors_1[row_num]
                color_2 = colors_2[row_num]
                color_3 = colors_3[row_num]
                chrom = chr_data['Accession'][row_num]
                orientation = chr_data['Orientation'][row_num]
                gene_id = int(chr_data['Gene ID'][row_num])

                add_point(chrom, chrom_dict, inputted_point, sector_obj, orientation, scatter_track, scatter2_track, scatter3_track, scatter4_track, \
                          bar_track, x, gene_colors, gene_id, full_data, gene_exp_gene_ids, gene_exp_pini_mock, \
                          gene_exp_capsici_mock, gene_exp_pini_capsici, gene_exp_avg_exp, color_1, color_2, color_3, bar_color)

        # If user manually included gene IDs
        if genomic_ranges:

            # Indicate that the point being added is from the geneIDs entered, thus larger and different color
            inputted_point = True

            for chrom, x, _, orientation, gene_id in genomic_ranges:
                add_point(chrom, chrom_dict, inputted_point, sector_obj, orientation, scatter_track, scatter2_track, scatter3_track, scatter4_track, \
                          bar_track, x, gene_colors, int(gene_id), full_data, gene_exp_gene_ids, gene_exp_pini_mock, \
                          gene_exp_capsici_mock, gene_exp_pini_capsici, gene_exp_avg_exp, color_1, color_2, color_3, bar_color)

    # Render the plot using Matplotlib
    fig = circos.plotfig()

    # Display the plot in Streamlit
    st.pyplot(fig)


main()
