import streamlit as st
import pandas as pd
from pycirclize import Circos
import matplotlib.pyplot as plt
from ncbi.datasets import GeneApi
from ncbi.datasets.openapi import ApiClient
from ncbi.datasets.openapi.rest import ApiException
from ncbi.datasets.openapi.model.v1_orientation import V1Orientation

# Load data (Chromosome sizes)
data = pd.DataFrame({
    'Chromosome': ['chr01', 'chr02', 'chr03', 'chr04', 'chr05', 'chr06', 'chr07', 'chr08', 'chr09', 'chr10', 
                   'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'Pltd'],  
    'Size (bp)': [45207397, 37821870, 35064427, 34823025, 22562875, 39020271, 52418484, 30564197, 24263475,
                  37707155, 37114715, 31492331, 39757759, 28841373, 20407330, 28711772, 160537]
})

# Non-numbered chromosome mapping (for plastid, mitochondrion, etc.)
non_numbered_chromosome = {
    'NC_028617.1': 'Pltd'  # Plastid for Juglans regia (English walnut)
}

numbered_chromosomes = {
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
    'NC_049916.1': 'chr16'
}

# Function to fetch gene metadata using GeneApi
def fetch_gene_metadata(gene_ids):
    api_client = ApiClient()  # Properly initialize the API client
    gene_api = GeneApi(api_client)
    try:
        # Fetch gene metadata for all gene IDs
        gene_metadata = gene_api.gene_metadata_by_id([int(gene_id) for gene_id in gene_ids])
        return gene_metadata
    except ApiException as e:
        st.error(f"An error occurred while fetching gene metadata: {e}")
        return None

# Function to format and display the full API response as a useful table
def display_formatted_response(response):
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

# Streamlit app
st.title('Custom Circos Plot Generator with Gene Metadata Integration')

# Allow user to upload optional gene expression file
gene_exp_file = st.file_uploader('Upload a gene expression file (optional, must be .csv, .xls or .xlsx)', type=["csv", "xls", "xlsx"])

# Read gene expression columns based on file type
if gene_exp_file is not None:
    
    if str(gene_exp_file).endswith('.csv'):
        gene_exp_df = pd.read_csv(gene_exp_file)
        csv = True
    else:
        gene_exp_df = pd.read_excel(gene_exp_file)
        csv = False

    # See if column headers are in row 1 or row 2
    col_row = st.selectbox('Are the headers in row 1 or row 2?', ['Row 1', 'Row 2'])

    if col_row == 'Row 1':
        header = 0
    else:
        header = 1
    
    # Re-read gene_exp_df with correct headers
    if csv:
        gene_exp_df = pd.read_csv(gene_exp_file, header = header)
    else:
        gene_exp_df = pd.read_excel(gene_exp_file, header = header)
    
    gene_exp_cols = gene_exp_df.columns

    # Allow user to specify which columns denote which values
    gene_exp_gene_ids = st.selectbox("Select which column has the GeneIDs", gene_exp_cols)
    gene_exp_avg_exp = st.selectbox("Select which column has the average expression level", gene_exp_cols)
    gene_exp_pini_mock = st.selectbox("Select which column has the log2FC CR10 pini / mock", gene_exp_cols)
    gene_exp_capsici_mock = st.selectbox("Select which column has the log2FC CR10 capsici / mock", gene_exp_cols)
    gene_exp_pini_capsici = st.selectbox("Select which column has the log2FC CR10 pini / capsici", gene_exp_cols)

# User input for multiple gene IDs
gene_id_input = st.text_input('Enter Gene IDs (space-separated)')

if gene_id_input:
    gene_ids = [gene_id.strip() for gene_id in gene_id_input.split() if gene_id.strip().isdigit()]

    if gene_ids:
        # Allow user to select a color for each gene ID
        gene_colors = {}
        for index, gene_id in enumerate(gene_ids):
            # Assign a unique key to each color picker using the gene ID and index
            color = st.color_picker(f"Pick a color for Gene ID {gene_id}", '#ff0000', key=f"color_picker_{gene_id}_{index}")
            gene_colors[gene_id] = color

        gene_metadata = fetch_gene_metadata(gene_ids)

        if gene_metadata:
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

                if len(genomic_ranges) > 0:
                    # Prepare Circos plot with sectors (using chromosome sizes)
                    sectors = {str(row['Chromosome']): row['Size (bp)'] for index, row in data.iterrows()}
                    circos = Circos(sectors, space=5)

                    for sector_obj in circos.sectors:
                        # Plot sector name
                        sector_obj.text(f"{sector_obj.name}", r=110, size=15)

                        # Add scatter plot points based on genomic ranges and orientation
                        scatter_track = sector_obj.add_track((95, 100), r_pad_ratio=0.1)
                        scatter_track.axis()

                        for chrom, start, end, orientation, gene_id in genomic_ranges:
                            # Check if the accession_version belongs to non-numbered chromosomes (e.g., plastid)
                            if chrom in non_numbered_chromosome:
                                chrom_name = non_numbered_chromosome[chrom]
                            else:
                                chrom_name = numbered_chromosomes.get(chrom, None)

                            if chrom_name is None:
                                st.warning(f"Chromosome {chrom} not found in mapping. Skipping.")
                                continue  # Skip if chromosome is not found

                            # Ensure the chromosome name matches the sector (e.g., 'Pltd', 'chr01')
                            if chrom_name == sector_obj.name:
                                x = start

                                # Extract the orientation value
                                if isinstance(orientation, V1Orientation):
                                    orientation_value = orientation.value
                                elif isinstance(orientation, str):
                                    orientation_value = orientation.strip().lower()
                                else:
                                    st.warning(f"Unexpected type for orientation at start {start}. Skipping this entry.")
                                    continue

                                # Assign y value based on orientation
                                y = 1 if orientation_value == 'plus' else 0 if orientation_value == 'minus' else None
                                if y is not None:
                                    # Use the color selected for the gene ID
                                    scatter_track.scatter([x], [y], color=gene_colors[str(gene_id)], label=f'Gene Start: {start} (Orientation: {orientation_value})')

                    # Render the plot using Matplotlib
                    fig = circos.plotfig()

                    # Display the plot in Streamlit
                    st.pyplot(fig)
                else:
                    st.warning(f"No valid genomic ranges found for Gene IDs: {', '.join(gene_ids)}")
            else:
                st.warning(f"No genes found for Gene IDs: {', '.join(gene_ids)}")
    else:
        st.error("Please enter valid numeric Gene IDs")
else:
    st.error("Please enter Gene IDs")

