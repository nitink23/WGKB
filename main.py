import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from pycirclize import Circos
import io

# Load data
data = pd.read_csv('internship Genome//all_transcripts.txt', delimiter='\t', index_col='GeneID')
metadata = pd.read_csv('internship Genome//all_metadata.txt', delimiter='\t')

# Ensure each column has the corresponding cultivar and treatment
column_info = {col: (metadata[metadata['Sample'] == col]['Cultivar'].values[0], metadata[metadata['Sample'] == col]['Treatment'].values[0]) for col in data.columns}

# Streamlit app
st.title('Circos Plot Generator')

# Select columns and rows
columns = st.multiselect('Select columns', data.columns, format_func=lambda x: f"{x} ({column_info[x][0]}, {column_info[x][1]})")
rows = st.multiselect('Select rows', data.index)

# Filter data
if columns and rows:
    filtered_data = data.loc[rows, columns]

    circos = Circos.initialize_from_matrix(
    filtered_data,
    cmap="tab10",
    start=-265,
    end=95,
    space=5,
    r_lim=(93, 100),
    label_kws=dict(r=94, size=12, color="white"),
    link_kws=dict(ec="black", lw=0.5),
    )
    fig = circos.plotfig()

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

    # Option to create another plot
    if st.button('Create another plot'):
        st.experimental_rerun()
    
else:
    st.write("Please select columns and rows to generate the plot.")
