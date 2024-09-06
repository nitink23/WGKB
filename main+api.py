import streamlit as st
import pandas as pd
from pycirclize import Circos
import matplotlib.pyplot as plt
from ncbi.datasets import GeneApi
from ncbi.datasets.openapi import ApiClient
from ncbi.datasets.openapi.rest import ApiException

# Load data (Chromosome sizes)
data = pd.DataFrame({
    'Chromosome': ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', 'Pltd'],
    'Size (bp)': [45207397, 37821870, 35064427, 34823025, 22562875, 39020271, 52418484, 30564197, 24263475,
                  37707155, 37114715, 31492331, 39757759, 28841373, 20407330, 28711772, 160537]
})

# Function to fetch gene metadata using GeneApi
def fetch_gene_metadata(gene_id):
    api_client = ApiClient()  # Properly initialize the API client
    gene_api = GeneApi(api_client)
    try:
        # Fetch gene metadata by gene ID (int conversion needed)
        gene_metadata = gene_api.gene_metadata_by_id([int(gene_id)])
        return gene_metadata
    except ApiException as e:
        st.error(f"An error occurred while fetching gene metadata: {e}")
        return None

# Streamlit app
st.title('Custom Circos Plot Generator with Gene Metadata Integration')

# User input for gene ID
gene_id_input = st.text_input('Enter Gene ID')

# Fetch gene metadata based on input
if gene_id_input.isdigit():  # Ensure valid integer input
    gene_metadata = fetch_gene_metadata(gene_id_input)

    if gene_metadata:
        # Log the entire API response for debugging
        st.write(f"Full API Response: {gene_metadata}")
        
        # Check if genes exist in the response
        if hasattr(gene_metadata, 'genes') and len(gene_metadata.genes) > 0:
            st.write(f"Fetched metadata for Gene ID: {gene_id_input}")
            
            # Parse the genomic ranges
            genomic_ranges = []
            for gene_data in gene_metadata.genes:
                gene_info = gene_data.gene
                
                if 'genomic_ranges' in gene_info and gene_info['genomic_ranges']:
                    for genomic_range in gene_info['genomic_ranges']:
                        for loc in genomic_range['range']:  # Access the range list within genomic_ranges
                            genomic_ranges.append((
                                genomic_range['accession_version'],  # Chromosome (accession)
                                int(loc['begin']),  # Start
                                int(loc['end'])     # End
                            ))
                else:
                    st.warning(f"No genomic ranges available for gene {gene_info['gene_id']}")
            
            if len(genomic_ranges) > 0:
                # Prepare Circos plot with sectors
                sectors = {row['Chromosome']: row['Size (bp)'] for index, row in data.iterrows()}
                circos = Circos(sectors, space=5)

                for sector_obj in circos.sectors:
                    # Plot sector name
                    sector_obj.text(f"Sector: {sector_obj.name}", r=110, size=15)
                    
                    # Add scatter plot points based on genomic ranges
                    points_track = sector_obj.add_track((15, 100), r_pad_ratio=0.1)
                    points_track.axis()
                    
                    for chrom, start, end in genomic_ranges:
                        if chrom.endswith(sector_obj.name):
                            # Plot the gene's start position as a point on the chromosome sector
                            points_track.scatter([start], [0], color='red', label=f'Gene Start: {start}')

                # Render the plot using Matplotlib
                fig = circos.plotfig()

                # Display the plot in Streamlit
                st.pyplot(fig)
            else:
                st.warning(f"No valid genomic ranges found for Gene ID: {gene_id_input}")
        else:
            st.warning(f"No genes found for Gene ID: {gene_id_input}")
else:
    st.error("Please enter a valid numeric Gene ID")
