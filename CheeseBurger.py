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
        track_correction = 0
        line, bar = False, False

        if full_data is not None and num_tracks >= 1:
            bar = st.checkbox("Would you like one of these tracks to be a bar plot? (Can only be visualized on bottom track)")
            line = st.checkbox("Would you like one of these tracks to be a line plot?")
            track_correction +=1 if bar else 0
            track_correction += 1 if line else 0
        
        # For each track, allow user to specify what data goes on which track
        track_cols = []
        bar_color = None

        # Add 'Gene Location' column if user entered Gene IDs
        available_cols = list(full_data_cols)
        if genomic_ranges:
            available_cols = ['Gene Location'] + available_cols

        # Allow user to select which data they want in which tracks
        for track in range(num_tracks-track_correction):

            desired_col = st.selectbox(f"Select which data you would like to visualize in track {track+1}:", available_cols)

            if desired_col == 'Gene Location':
                desired_data = genomic_ranges
                omit_pct = 0
            else:
                desired_data = full_data[desired_col]

                if invalid_col(full_data, desired_col):
                    return
                
                # Add slider for user to omit certain values close to mean
                omit_pct = st.slider(
                    label="Choose a cutoff to omit points close to the mean",
                    min_value=0,
                    max_value=100,
                    value=0,
                    step=1,
                    key='track_slider_' + str(track)
                )

            track_cols.append([desired_col, 'dot', desired_data, omit_pct])
        
        # Add data to track_cols dict when user wants line or bar track
        for track_type in ['line' if line else None, 'bar' if bar else None]:

            if track_type:

                desired_col = st.selectbox(f"Select which data you would like to visualize in {track_type} track:", available_cols)

                if track_type == 'bar':
                    bar_color = st.color_picker(f"Pick a color for bar plot", '#ff0000', key=f"color_picker_bar_track")

                if invalid_col(full_data, desired_col):
                    return
                
                if desired_col != 'Gene Location':
                    desired_data = full_data[desired_col]
                    track_cols.append([desired_col, track_type, desired_data, 0])

        # Plot genes not associated with any chromosome
        if st.checkbox('Would you like to plot the genes not associated with a chromosome?'):

            desired_col = st.selectbox(f"Select which data you would like to visualize in remaining genes track:", available_cols)
            rem_genes_color = st.color_picker(f"Pick a color for remaining genes plot", '#ff0000', key=f"rem_genes_color")
            omit_pct = st.slider(
                    label="Choose a cutoff to omit points close to the mean",
                    min_value=0,
                    max_value=100,
                    value=0,
                    step=1,
                    key='track_slider_rem_genes'
                    )

            if st.button('Click to plot genes without chromosome'):
                if desired_col != 'Gene Location':
                    plot_rem_genes(full_data, desired_col, omit_pct, rem_genes_color)
                else:
                    st.warning('Cannot visualize Gene Location on remaining genes track.')

        # Plot main Circos plot
        try:
            if st.button('Click to plot!'):
                st.write('Plotting!')
                display_circos_plot(data, full_data, track_cols, bar_color, genomic_ranges)
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
            st.error(f'The \'{desired_col}\' column contains string values and thus cannot be visualized on the Circos plot. The app can only read integers or floats.')
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


def display_circos_plot(data: dict, full_data, track_cols: list, bar_color, genomic_ranges) -> None:
    
    # Prepare Circos plot with sectors (using chromosome sizes)
    sectors = {str(row[1]['Chromosome']): row[1]['Size (bp)'] for row in data.iterrows()}
    circos = Circos(sectors, space=5)

    lower = 105

    # Format: track_cols[0] = [desired_col, plot_type, desired_data, omit_pct]
    for index, [desired_col, plot_type, desired_data, omit_pct] in enumerate(track_cols):

        upper = lower - 5
        if desired_col == 'Gene Location':
            width = 5
        elif plot_type == 'bar':
            width = 75/len(track_cols)
        else:
            width = 50/len(track_cols)
        lower = upper - width

        # Get the colors for given column if it is not Gene Location
        if desired_col != 'Gene Location':
            full_data['color_' + str(index)] = get_colormap(desired_data)

        for chrom_num, sector_obj in enumerate(circos.sectors):
            
            # Plot sector name
            sector_obj.text(f"{get_chrom_num(sector_obj.name)}", r=110, size=10)

            # Add sector with correct sizing per number of tracks
            track = sector_obj.add_track((lower, upper), r_pad_ratio=0.1)
            track.axis()

            # Add x-ticks for gene location
            if index == 0:  # Apply ticks only to the outermost track

                # Major ticks (every 20 Mb)
                track.xticks_by_interval(
                    interval=20 * 1_000_000,
                    label_orientation="vertical",
                    tick_length=2,
                    label_formatter=lambda v: f"{v / 1_000_000:.0f}Mb",
                    label_size=5
                )

                # Minor ticks (every 5 Mb)
                track.xticks_by_interval(
                    interval=5 * 1_000_000,
                    tick_length=1,
                    show_label=False,  # No labels for minor ticks
                )
                
            # Add y-ticks only on the left side of the chr01 sector
            if get_chrom_num(sector_obj.name) == 'chr01':  

                if desired_col == 'Gene Location':
                    track.yticks([0, 1], list("-+"), vmin=0, vmax=1, side="left")

                # Only add tick marks to tracks that are not bar or line    
                elif plot_type not in ['bar', 'line']:
                    track.yticks(y=np.linspace(desired_data.min(), desired_data.max(), num=5), \
                                 labels=[f"{round(tick)}" for tick in np.linspace(desired_data.min(), desired_data.max(), num=5)], \
                                 vmin=desired_data.min(), vmax=desired_data.max(), side="left", label_size=7-index)
        
            # If given track is not Gene Location
            if desired_col != 'Gene Location':

                # Get subset from full_data of all the rows with a specific chromosome
                chr_data = full_data[full_data['Chromosome'] == str(chrom_num+1)].reset_index(drop=True)
                
                if full_data is not None:

                    if plot_type == 'dot':

                        # Remove points close to the mean if user selected:
                        if omit_pct > 0:

                            lower_bound = desired_data.quantile((100 - omit_pct) / 200)
                            upper_bound = desired_data.quantile(1 - (100 - omit_pct) / 200)

                            chr_data = chr_data[(chr_data[desired_col] < lower_bound) | (chr_data[desired_col] > upper_bound)]

                        track.scatter(chr_data['Begin'].tolist(), chr_data[desired_col].tolist(), color=chr_data['color_' + str(index)], 
                                      cmap='coolwarm', vmin=desired_data.min(), vmax=desired_data.max(), s=5)

                    if not pd.isna(chr_data[desired_col].max()):

                        if plot_type == 'line':

                            # Add solid line connecting data points
                            chr_data_sorted = chr_data.sort_values(by='Begin')
                            track.line(chr_data_sorted['Begin'].tolist(), chr_data_sorted[desired_col].tolist(), vmin=desired_data.min(), vmax=full_data[desired_col].max())

                        elif plot_type == 'bar':

                            # Ensure the column has values in it
                            vmin = 0 if full_data[desired_col].min() > 0 else full_data[desired_col].min()
                            track.bar(chr_data['Begin'].tolist(), chr_data[desired_col].tolist(), ec=bar_color, lw=0.9, vmin=vmin, vmax=full_data[desired_col].max())

            # If user manually included gene IDs
            else:
                for chrom, x, _, orientation, _, color_to_use in desired_data:
                    if chrom == sector_obj.name:
                        y = 1 if orientation.value == 'plus' else 0
                        track.scatter([x], [y], color=color_to_use, s=20)

    # Add colorbar for expression level columns
    exp_cols = [[desired_col, plot_type, index] for index, [desired_col, plot_type, _, _] in enumerate(track_cols) if desired_col != 'Gene Location']

    for index, (label, plot_type, track_index) in enumerate(exp_cols):

        if plot_type == 'bar' or plot_type == 'line':

            # Shorten title if we are examining Expression level
            if label == 'Expression level (average normalized FPKM of all 36 samples)':
                label = 'Average normalized FPKM expression'

            circos.colorbar(
                bounds=(0.92, 1+index*0.1, 0.2, 0),
                vmin=full_data[label].min(),
                vmax=full_data[label].max(),
                orientation="horizontal",
                label=f'Track {track_index+1}: ' + label,
                label_kws=dict(size=8, color="black"),
                tick_kws=dict(labelsize=8, colors="black")
            )

        else:
            circos.colorbar(
                bounds=(0.92, 1+index*0.1, 0.2, 0.02),
                vmin=full_data[label].min(),
                vmax=full_data[label].max(),
                cmap="coolwarm",
                orientation="horizontal",
                label=f'Track {track_index+1}: ' + label,
                label_kws=dict(size=8, color="black"),
                tick_kws=dict(labelsize=8, colors="black")
            )

    plt.tight_layout()

    # Render the plot using Matplotlib
    fig = circos.plotfig()

    # Add legend for selected gene IDs
    if 'Gene Location' in [col for col, _, _, _ in track_cols]:
         
        scatter_legend = circos.ax.legend(
            handles=[plt.Line2D([0], [0], color=row[-1], marker='o', ls='None', ms=8) for row in genomic_ranges],  
            labels=[row[4] for row in genomic_ranges],  
            bbox_to_anchor=(-0.14, 1.1),
            loc='upper left',
            fontsize=8,
            title="Gene ID",
            handlelength=2
        )

        scatter_legend._legend_title_box._text_pad = 5
        circos.ax.add_artist(scatter_legend)

    # Add legend for sliders
    removal_legend_items = [
        (f"Track {track_index + 1}: Omit {omit_pct}% around median", track_index)
        for track_index, (_, _, _, omit_pct) in enumerate(track_cols)
        if omit_pct > 0
    ]
    # Add removal legend if there are items
    if removal_legend_items:

        # Create handles and labels for all tracks in one go
        legend_handles = [
            plt.Line2D([0], [0], color="black", linestyle='None')
            for _ in removal_legend_items
        ]
        legend_labels = [label for label, _ in removal_legend_items]

        removal_legend = circos.ax.legend(
            handles=legend_handles,  
            labels=legend_labels,  
            bbox_to_anchor=(-0.17, 0),
            loc='lower left',
            title="Median Percentile Removal",
            fontsize=6,
            title_fontsize=7,
            handlelength=1.5
        )

        removal_legend._legend_title_box._text_pad = 5
        circos.ax.add_artist(removal_legend)
        
    # Display the plot in Streamlit
    st.pyplot(fig)


def plot_rem_genes(full_data, desired_col, omit_pct, color) -> None:

    rem_genes_df = full_data[pd.isna(full_data['Chromosome'])].reset_index()
    desired_data = rem_genes_df[desired_col]

    sectors = {rem_genes_df['Accession'].iloc[0]: max(rem_genes_df.index.tolist())}
    circos = Circos(sectors, space=10)

    for sector_obj in circos.sectors:

        track = sector_obj.add_track((70, 100), r_pad_ratio=0.1)
        track.axis()

        # Add y-ticks only on the left side
        track.yticks(y=np.linspace(min(desired_data), max(desired_data), num=5), \
                     labels=[f"{round(tick)}" for tick in np.linspace(min(desired_data), max(desired_data), num=5)], \
                     vmin=min(desired_data), vmax=max(desired_data), side="left", label_size=7)

        if omit_pct > 0:

            lower_bound = desired_data.quantile((100 - omit_pct) / 200)
            upper_bound = desired_data.quantile(1 - (100 - omit_pct) / 200)

            desired_data = desired_data[(desired_data < lower_bound) | (desired_data > upper_bound)]

        rem_genes_x = desired_data.index.tolist()
        desired_data = desired_data.tolist()

        track.scatter(x=rem_genes_x, y=desired_data, color=color, 
                      vmin=min(desired_data), vmax=max(desired_data), s=12)
        
        # Add y-axis gridlines manually
        track.grid()

        # Add xticks for index
        indices = list(range(len(desired_data)))

        # Major ticks
        track.xticks_by_interval(
            interval=max(indices) // 10,
            show_label=True,
            label_orientation="vertical",
            tick_length=2,
            label_size=10,
            outer=False,
        )

        # Minor ticks
        minor_interval = 1 if max(indices) // 50 == 0 else max(indices) // 50
        track.xticks_by_interval(
            interval=minor_interval,  # Minor ticks at a smaller interval
            outer=False,
            tick_length=1,  # Shorter tick marks for minor ticks
            show_label=False,  # No labels for minor ticks
        )
        
    fig_2 = circos.plotfig()
        
    # Add sector name as a legend
    scatter_legend = circos.ax.legend(
        bbox_to_anchor=(1.05, 1.05),
        loc='upper right',
        fontsize=8,
        title=sector_obj,
        handlelength=2
    )

    scatter_legend._legend_title_box._text_pad = 5
    circos.ax.add_artist(scatter_legend)

    if omit_pct > 0:

        legend_handles = [plt.Line2D([0], [0], color="black", linestyle='None')]
        legend_labels = [f"Remaining Genes Track: Omit {omit_pct}% around median"]

        removal_legend = circos.ax.legend(
            handles=legend_handles,  
            labels=legend_labels,  
            bbox_to_anchor=(-0.1, 0),
            loc='lower left',
            title="Median Percentile Removal",
            fontsize=9,
            title_fontsize=10,
            handlelength=1.5
        )

        removal_legend._legend_title_box._text_pad = 5
        circos.ax.add_artist(removal_legend)

    st.pyplot(fig_2)

main()
