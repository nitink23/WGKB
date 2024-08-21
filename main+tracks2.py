import streamlit as st
import pandas as pd
import numpy as np
from pycirclize import Circos
import io

# Load data
data = pd.read_csv('internship Genome//all_transcripts.txt', delimiter='\t', index_col='GeneID')
metadata = pd.read_csv('internship Genome//all_metadata.txt', delimiter='\t')

# Streamlit app
st.title('Custom Circos Plot Generator')

# Select a condition to compare
condition_choice = st.selectbox(
    'Choose comparison type',
    ['Variability among replicates (single condition)', 'Comparison of conditions (averaged)']
)

# Filter columns based on selection
if condition_choice == 'Variability among replicates (single condition)':
    condition = st.selectbox('Select condition', ['Chandler Mock'])
    columns = metadata[(metadata['Cultivar'] == 'Chandler') & (metadata['Treatment'] == 'Mock')]['Sample'].tolist()
elif condition_choice == 'Comparison of conditions (averaged)':
    conditions = ['Chandler Mock', 'Chandler P.cap', 'Chandler P.pini']
    selected_conditions = st.multiselect('Select conditions to compare', conditions)
    
    columns = []
    if 'Chandler Mock' in selected_conditions:
        mock_samples = metadata[(metadata['Cultivar'] == 'Chandler') & (metadata['Treatment'] == 'Mock')]['Sample']
        data['Chandler_Mock_Avg'] = data[mock_samples].mean(axis=1)
        columns.append('Chandler_Mock_Avg')

    if 'Chandler P.cap' in selected_conditions:
        cap_samples = metadata[(metadata['Cultivar'] == 'Chandler') & (metadata['Treatment'] == 'P. cap')]['Sample']
        data['Chandler_Pcap_Avg'] = data[cap_samples].mean(axis=1)
        columns.append('Chandler_Pcap_Avg')
    
    if 'Chandler P.pini' in selected_conditions:
        pini_samples = metadata[(metadata['Cultivar'] == 'Chandler') & (metadata['Treatment'] == 'P. pini')]['Sample']
        data['Chandler_Ppini_Avg'] = data[pini_samples].mean(axis=1)
        columns.append('Chandler_Ppini_Avg')

# Select rows (GeneIDs)
rows = st.multiselect('Select rows (GeneIDs)', data.index)

# Select plot types for each track
plot_types = st.multiselect('Select plot types for tracks', ['Line', 'Scatter', 'Bar'], format_func=lambda x: x)

# Filter data and prepare the plot
if columns and rows:
    filtered_data = data.loc[rows, columns]

    # Debugging: Show filtered data
    st.write("Filtered Data:")
    st.write(filtered_data)

    try:
        # Prepare Circos plot
        sectors = {col: len(filtered_data) for col in columns}  # Use length of filtered data for sector sizes
        circos = Circos(sectors, space=10, start=90, end=360, endspace=False)

        for i, col in enumerate(columns):
            sector = circos.get_sector(col)  # Use get_sector to retrieve sector by name
            x = np.arange(len(filtered_data))
            y = filtered_data[col].values
            
            # Outer Track - for sector labels
            outer_track = sector.add_track((95, 100))
            outer_track.text(col, color="white")
            outer_track.axis(fc="grey")
            outer_track.xticks_by_interval(interval=10, label_orientation="vertical")

            # Add tracks based on selected plot types
            for plot_type in plot_types:
                if plot_type == 'Line':
                    line_track = sector.add_track((75, 90), r_pad_ratio=0.1)
                    line_track.axis()
                    line_track.line(x, y, color="blue")
                elif plot_type == 'Scatter':
                    scatter_track = sector.add_track((45, 70), r_pad_ratio=0.1)
                    scatter_track.axis()
                    scatter_track.scatter(x, y, color="green", s=3)
                elif plot_type == 'Bar':
                    bar_track = sector.add_track((15, 40), r_pad_ratio=0.1)
                    bar_track.axis()
                    bar_track.bar(x, y, color="orange")
            

        # Plot figure
        fig = circos.plotfig()

        # Debugging: Ensure figure is created
        st.write("Plot created successfully.")

        # Display plot
        st.pyplot(fig)

        # Allow user to download the plot
        buf = io.BytesIO()
        fig.savefig(buf, format="png")
        st.download_button(
            label="Download Plot",
            data=buf.getvalue(),
            file_name="circos_plot.png",
            mime="image/png"
        )

    except Exception as e:
        st.error(f"An error occurred: {e}")

    # Option to create another plot
    if st.button('Create another plot'):
        st.experimental_rerun()
else:
    st.write("Please select columns and rows to generate the plot.")
