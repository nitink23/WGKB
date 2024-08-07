import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pycirclize import Circos
import io

# Load data
data = pd.read_csv('internship Genome//all_transcripts.txt', delimiter='\t', index_col='GeneID')
metadata = pd.read_csv('internship Genome//all_metadata.txt', delimiter='\t')

# Ensure each column has the corresponding cultivar and treatment
column_info = {col: (metadata[metadata['Sample'] == col]['Cultivar'].values[0], metadata[metadata['Sample'] == col]['Treatment'].values[0]) for col in data.columns}

# Streamlit app
st.title('Enhanced Circos Plot Generator')

# Select columns and rows
columns = st.multiselect('Select columns', data.columns, format_func=lambda x: f"{x} ({column_info[x][0]}, {column_info[x][1]})")
rows = st.multiselect('Select rows', data.index)

# Select plot types for each track
plot_types = st.multiselect('Select plot types for tracks', ['Line', 'Scatter', 'Bar'], format_func=lambda x: x)

# Filter data
if columns and rows:
    filtered_data = data.loc[rows, columns]

    # Debugging: Show filtered data
    st.write("Filtered Data:")
    st.write(filtered_data)
    
    try:
        # Prepare Circos plot
        sectors = {col: len(filtered_data) for col in columns}  # Use length of filtered data for sector sizes
        circos = Circos(sectors, space=5)
        
        for col in columns:
            sector = circos.get_sector(col)  # Use get_sector to retrieve sector by name
            x = np.arange(len(filtered_data))
            y = filtered_data[col].values
            
            # Debugging: Show sector name and data
            st.write(f"Sector: {col}")
            st.write(f"x: {x}")
            st.write(f"y: {y}")
            
            # Add tracks based on selected plot types
            for plot_type in plot_types:
                if plot_type == 'Line':
                    line_track = sector.add_track((75, 100), r_pad_ratio=0.1)
                    line_track.axis()
                    line_track.xticks_by_interval(1)
                    line_track.line(x, y)
                elif plot_type == 'Scatter':
                    points_track = sector.add_track((45, 70), r_pad_ratio=0.1)
                    points_track.axis()
                    points_track.scatter(x, y)
                elif plot_type == 'Bar':
                    bar_track = sector.add_track((15, 40), r_pad_ratio=0.1)
                    bar_track.axis()
                    bar_track.bar(x, y)
        
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


