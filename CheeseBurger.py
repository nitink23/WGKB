import streamlit as st
import pandas as pd
from pycirclize import Circos
import matplotlib.pyplot as plt
from ncbi.datasets import GeneApi
from ncbi.datasets.openapi import ApiClient
from ncbi.datasets.openapi.rest import ApiException
import numpy as np


def main() -> None:
    data = pd.DataFrame({
    'Chromosome': ['NC_049901.1', 'NC_049902.1', 'NC_049903.1', 'NC_049904.1', 'NC_049905.1', 'NC_049906.1', 
                   'NC_049907.1', 'NC_049908.1', 'NC_049909.1', 'NC_049910.1', 'NC_049911.1', 'NC_049912.1',
                   'NC_049913.1', 'NC_049914.1', 'NC_049915.1', 'NC_049916.1', 'NC_028617.1'],  
    'Size (bp)': [45207397, 37821870, 35064427, 34823025, 22562875, 39020271, 52418484, 30564197, 24263475,
                  37707155, 37114715, 31492331, 39757759, 28841373, 20407330, 28711772, 160537]})

    st.title('Custom Circos Plot Generator with Gene Metadata Integration')
    st.markdown("###### Created by Paulo Zaini, Adam Hetherwick, Hibiki Ono, and Nitin Kanchi")
    st.markdown("### Enter GeneIDs and/or a gene expression file to get started")

    # Allow user to upload optional gene expression file
    gene_exp_file = st.file_uploader('Upload a gene expression file (optional, must be .csv, .xls or .xlsx)', type=["csv", "xls", "xlsx"])

    # Read walnut gene metadata file straight from GitHub
    url = 'https://raw.githubusercontent.com/nitink23/WGKB/main/ncbi_dataset.tsv'
    walnut_gene_meta = pd.read_csv(url, delimiter='\t')

    full_data, full_data_cols = None, []

    if gene_exp_file:
        try:
            gene_exp_df = read_gene_exp_file(gene_exp_file)
            full_data = pd.merge(gene_exp_df, walnut_gene_meta, on='Gene ID', how='inner')
            full_data_cols = full_data.columns

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
                    st.error(f'ERROR: {gene_id} is not a valid Gene ID')
                    return

        gene_metadata = fetch_gene_metadata(gene_ids)
        genomic_ranges = display_gene_meta(gene_metadata, gene_ids)

        # Allow user to select a color for each gene ID
        for index, row in enumerate(genomic_ranges):

            # Assign a unique key to each color picker using the gene ID and index
            color = st.color_picker(f"Pick a color for Gene ID {row[4]}", '#ff0000', key=f"color_picker_{row[4]}_{index}")
            row.append(color)

        # Check to see if genomic ranges are valid, if invalid print warning and stop function
        if len(genomic_ranges) <= 0:
            gene_id_input = []
            genomic_ranges = []
            st.error(f"No valid genomic ranges found for Gene IDs: {', '.join(gene_ids)}")
            return

    if full_data is not None or gene_id_input:
        
        # Allow user to specify how many tracks they want to visualize
        available_tracks = [1, 2, 3, 4, 5] if full_data is not None else [1]
        num_tracks = st.selectbox("How many tracks would you like to visualize? (Upload gene expression file to view more than 1 track)", available_tracks)
        bar_correction = 0
        bar = False

        if full_data is not None and num_tracks > 1:
            bar = st.checkbox(f"Would you like to include a bar plot? (Can only be visualized on bottom track)")
            bar_correction = 1 if bar else 0
        
        # For each track, allow user to specify what data goes on which track
        track_cols = {}
        bar_color = None

        # Add 'Gene Location' column if user entered Gene IDs
        available_cols = list(full_data_cols)
        if genomic_ranges:
            available_cols = ['Gene Location'] + available_cols

        # Allow user to select which data they want in which tracks
        for track in range(num_tracks-bar_correction):

            desired_col = st.selectbox(f"Select which data you would like to visualize in track {track+1}:", available_cols)

            if desired_col == 'Gene Location':
                desired_data = genomic_ranges
            else:
                desired_data = full_data[desired_col]

                if invalid_col(full_data, desired_col):
                    return

            track_cols[desired_col] = ['dot', desired_data]
        
        # Check to see if the user wants to add a bar track
        if bar:
            desired_col = st.selectbox(f"Select which data you would like to visualize in bar track:", available_cols)
            bar_color = st.color_picker(f"Pick a color for bar plot", '#ff0000', key=f"color_picker_bar_track")
            if invalid_col(full_data, desired_col):
                return
            if desired_col != 'Gene Location':
                desired_data = full_data[desired_col]
                track_cols[desired_col] = ['bar', desired_data]

        try:
            if st.button('Click to plot!'):
                st.write('Plotting!')
                display_circos_plot(data, full_data, track_cols, bar_color)
        except (KeyError):
            st.error('WARNING: There was an error displaying the plot.')
            return


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


def invalid_col(full_data, desired_col) -> bool:
    # Checks for string values in given column
    if desired_col != 'Gene Location':
        if full_data[desired_col].apply(lambda x: isinstance(x, str)).all():
            st.error(f'The \'{desired_col}\' column contains string values and thus cannot be visualized on the Circos plot. The app can only read integers.')
            return True
    return False


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
    
    # See if column column headers are in row 1 or row 2
    col_row = st.selectbox('Select which row the column headers are in:', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    
    # Re-read gene_exp_df with correct headers
    if str(gene_exp_file).endswith('.csv'):
        gene_exp_df = pd.read_csv(gene_exp_file, header = col_row-1)
    else:
        gene_exp_df = pd.read_excel(gene_exp_file, header = col_row-1)

    return gene_exp_df
    

def get_colormap(data_col):
    # Convert data_col to a numpy array
    vectorized_color_col = np.array(data_col)

    # Handle case where data_col is empty or has NaNs
    if vectorized_color_col.size == 0 or np.isnan(vectorized_color_col).all():
        return ['#ff0000'] * len(data_col)
    
    # Filter out NaNs for normalization
    finite_data = vectorized_color_col[np.isfinite(vectorized_color_col)]
    min_val = finite_data.min() if finite_data.size > 0 else 0
    max_val = finite_data.max() if finite_data.size > 0 else 1
    
    # Normalize the array, ensuring the range is 0-1
    if min_val == max_val:
        normalized_arr = np.zeros_like(vectorized_color_col)  # Set all to 0 if min and max are equal
    else:
        normalized_arr = (vectorized_color_col - min_val) / (max_val - min_val)

    # Apply colormap
    cmap = plt.get_cmap('coolwarm')
    colors = [cmap(value) if np.isfinite(value) else (0.5, 0.5, 0.5, 1) for value in normalized_arr]

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
                        genomic_ranges.append([chrom, start, end, orientation, gene_info['gene_id']])
            else:
                st.warning(f"No genomic ranges available for gene {gene_info['gene_id']}")
        return genomic_ranges
    else:
        st.warning(f'No gene ID metadata found for {gene_ids}')


def get_chrom_num(key: str) -> str:
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
    if key in chrom_dict.keys():
        return chrom_dict[key]
    else:
        return None


def display_circos_plot(data: dict, full_data, track_cols: dict, bar_color) -> None:
    
    # Prepare Circos plot with sectors (using chromosome sizes)
    sectors = {str(row[1]['Chromosome']): row[1]['Size (bp)'] for row in data.iterrows()}
    circos = Circos(sectors, space=5)

    lower = 105

    # Format: track_cols[desired_col] = [plot_type, desired_data]
    for index, (key, val) in enumerate(track_cols.items()):

        upper = lower - 5
        if key == 'Gene Location':
            width = 5
        elif val[0] == 'bar':
            width = 75/len(track_cols)
        else:
            width = 50/len(track_cols)
        lower = upper - width

        # Get the colors for given column
        full_data['color_' + str(index)] = get_colormap(val[1])

        for chrom_num, sector_obj in enumerate(circos.sectors):
            
            # Plot sector name
            sector_obj.text(f"{get_chrom_num(sector_obj.name)}", r=110, size=10)

            # Add sector with correct sizing per number of tracks
            track = sector_obj.add_track((lower, upper), r_pad_ratio=0.1)
            track.axis()

            # Add y-ticks only on the left side of the chr01 sector
            if get_chrom_num(sector_obj.name) == 'chr01':  
                if key == 'Gene Location':
                    track.yticks([0, 1], list("-+"), vmin=0, vmax=1, side="left")
                else:
                    track.yticks(y=np.linspace(val[1].min(), val[1].max(), num=5), \
                                 labels=[f"{round(tick)}" for tick in np.linspace(val[1].min(), val[1].max(), num=5)], \
                                 vmin=val[1].min(), vmax=val[1].max(), side="left", label_size=4)
        
            # If given track is not Gene Location
            if key != 'Gene Location':
                
                if full_data is not None:

                    if val[0] == 'dot':

                        # Get subset from full_data of all the rows with a specific chromosome
                        chr_data = full_data[full_data['Chromosome'] == str(chrom_num+1)].reset_index(drop=True)

                        track.scatter(chr_data['Begin'].tolist(), chr_data[key].tolist(), color=chr_data['color_' + str(index)], cmap='coolwarm', vmin=val[1].min(), vmax=val[1].max(), s=5)

                    else:
                        track.bar(chr_data['Begin'].tolist(), chr_data[key].tolist(), ec=bar_color, lw=0.9, vmin=0, vmax=val[1].max())

            # If user manually included gene IDs
            else:
                for chrom, x, _, orientation, _, color_to_use in val[1]:
                    if chrom == sector_obj.name:
                        y = 1 if orientation.value == 'plus' else 0
                        track.scatter([x], [y], color=color_to_use, s=20)

    # Add colorbar for expression level columns
    exp_cols = [col for col, val in track_cols.items() if col != 'Gene Location' and val[0] == 'dot']
    for index, col in enumerate(exp_cols):
        circos.colorbar(
            bounds=(1-0.25*index, 1.1, 0.02, 0.3),
            vmin=full_data[col].min(),
            vmax=full_data[col].max(),
            cmap="coolwarm",
            orientation="vertical",
            label=col,
            label_kws=dict(size=10, color="black"),
            tick_kws=dict(labelsize=8, colors="black")
        )

    plt.tight_layout()

    # Render the plot using Matplotlib
    fig = circos.plotfig()

    # Display the plot in Streamlit
    st.pyplot(fig)


main()
