"""
Created on Fri Jan 17 13:53:51 2025

@author: Alex Fernandes
@GitHub: alexfernandes8a
"""

import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objs as go
import pandas as pd
import scanpy as sc
from scipy.stats import zscore
import numpy as np
import argparse
import os
import requests
import base64
import gdown
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run the Dash app.")
parser.add_argument("--host", default="127.0.0.1", help="Host IP address")
parser.add_argument("--port", type=int, default=8050, help="Port number")
args = parser.parse_args()

# Define the file path and URL
file_path = os.getenv("GDRIVE_FILE_PATH")
url = os.getenv("GDRIVE_URL")
md5 = os.getenv("GDRIVE_MD5")

if not url:
    raise ValueError("GDRIVE_URL environment variable is not set.")

# Function to download the file if it doesn't exist
def download_file(url, file_path):
    if not os.path.exists(file_path):
        print(f"File not found at {file_path}. Downloading from {url}...")
        gdown.download(url, file_path, fuzzy=True, quiet=False)
        gdown.cached_download(url, file_path, hash=md5, quiet=False)
        print(f"File downloaded and saved to {file_path}.")
    else:
        print(f"File already exists at {file_path}.")

# Download the file if necessary
download_file(url, file_path)

# Load your preprocessed data
adata = sc.read_h5ad(file_path)

# Define variables for easier modification
CELLTYPE_COLUMN = "scDeepSort"  # Column in adata.obs for cell types
AGE_COLUMN = "age"  # Column in adata.obs for age
AGE_ORDER = ["E11", "E12", "E14", "E16", "E18", "P0", "P2", "P5", "P8", "P14"]  # Age order

# Ensure the age column is categorical and ordered
adata.obs[AGE_COLUMN] = pd.Categorical(adata.obs[AGE_COLUMN], categories=AGE_ORDER, ordered=True)

# Create a Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Load GitHub icon
github_icon_base64 = base64.b64encode(open("figures/github-icon.png", "rb").read()).decode('ascii')
github_icon = f"data:image/png;base64,{github_icon_base64}"

# Layout of the app
app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H1("scChronoScope", style={"font-size": "3rem"}),
            html.H2("Single-Cell RNA-Seq Data Explorer", style={"font-size": "1.5rem", "margin-top": "-10px"}),
            html.Div([
                html.Span("Developed by Alex Fernandes "),
                html.A(
                    html.Img(src=github_icon, style={"width": "20px", "height": "20px", "vertical-align": "middle"}),
                    href="https://github.com/alexfernandes8a",
                    target="_blank"
                )
            ], style={"margin-top": "10px"})
        ], className="mb-4")
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("Select Gene to Visualize"),
            dcc.Dropdown(
                id="gene-selector",
                options=[{"label": gene, "value": gene} for gene in adata.var_names],
                placeholder="Select a gene..."
            )
        ], width=6),
        dbc.Col([
            html.Label("Select Cell Type to Visualize"),
            dcc.Dropdown(
                id="cell-type-selector",
                options=[{"label": ct, "value": ct} for ct in adata.obs[CELLTYPE_COLUMN].unique()],
                placeholder="Select a cell type..."
            )
        ], width=6)
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("Select UMAP Visualization Type"),
            dbc.RadioItems(
                id="umap-type-selector",
                options=[
                    {"label": "2D UMAP", "value": "2D"},
                    {"label": "3D UMAP", "value": "3D"}
                ],
                value="2D",  # Default to 2D
                inline=True
            )
        ], width=6)
    ]),
    dbc.Row([
        dbc.Col(dcc.Graph(id="umap-gene-expression"), width=6),  # Plot 1
        dbc.Col(dcc.Graph(id="umap-cell-type"), width=6)  # Plot 2
    ]),
    dbc.Row([
        dbc.Col(dcc.Graph(id="expression-by-cell-type-boxplot"), width=12)  # Plot 3
    ]),
    dbc.Row([
        dbc.Col(dcc.Graph(id="expression-by-age-lineplot"), width=12)  # Plot 4
    ]),
    dbc.Row([
        dbc.Col(dcc.Graph(id="expression-heatmap"), width=12)  # Heatmap
    ])
])

# Callback to update plots based on gene, cell type, and UMAP type selection
@app.callback(
    [Output("umap-gene-expression", "figure"),
     Output("umap-cell-type", "figure"),
     Output("expression-by-cell-type-boxplot", "figure"),
     Output("expression-by-age-lineplot", "figure"),
     Output("expression-heatmap", "figure")],  # Add output for heatmap
    [Input("gene-selector", "value"),
     Input("cell-type-selector", "value"),
     Input("umap-type-selector", "value")]
)
def update_plots(selected_gene, selected_cell_type, umap_type):
    # Base UMAP coordinates
    if umap_type == "2D":
        umap_coords = pd.DataFrame({
            "UMAP1": adata.obsm["X_umap"][:, 0],
            "UMAP2": adata.obsm["X_umap"][:, 1],
            "Cell Type": adata.obs[CELLTYPE_COLUMN],
            "Age": adata.obs[AGE_COLUMN]
        })
    else:  # 3D UMAP
        umap_coords = pd.DataFrame({
            "UMAP1": adata.obs["umap_coord1"],
            "UMAP2": adata.obs["umap_coord2"],
            "UMAP3": adata.obs["umap_coord3"],
            "Cell Type": adata.obs[CELLTYPE_COLUMN],
            "Age": adata.obs[AGE_COLUMN]
        })

    # Calculate dynamic dot size based on the number of cells
    num_cells = umap_coords.shape[0]
    base_size = 5  # Base size for dots
    dot_size = max(1, base_size - 0.01 * num_cells)  # Adjust size based on number of cells

    # Plot 1: Gene Expression (UMAP)
    if selected_gene:
        umap_coords["Expression"] = adata[:, selected_gene].X.toarray().flatten()
        if umap_type == "2D":
            gene_fig = px.scatter(
                umap_coords, x="UMAP1", y="UMAP2", color="Expression",
                title=f"Gene Expression: {selected_gene}",
                color_continuous_scale="Viridis",
                size_max=dot_size  # Set dynamic dot size
            )
        else:  # 3D UMAP
            gene_fig = px.scatter_3d(
                umap_coords, x="UMAP1", y="UMAP2", z="UMAP3", color="Expression",
                title=f"Gene Expression: {selected_gene}",
                color_continuous_scale="Viridis"
            )
            gene_fig.update_traces(marker=dict(size=dot_size))  # Set dynamic dot size for 3D
    else:
        if umap_type == "2D":
            gene_fig = px.scatter(
                umap_coords, x="UMAP1", y="UMAP2",
                title="Gene Expression (No Gene Selected)",
                size_max=dot_size  # Set dynamic dot size
            )
        else:  # 3D UMAP
            gene_fig = px.scatter_3d(
                umap_coords, x="UMAP1", y="UMAP2", z="UMAP3",
                title="Gene Expression (No Gene Selected)"
            )
            gene_fig.update_traces(marker=dict(size=dot_size))  # Set dynamic dot size for 3D

    # Plot 2: Cell Type (UMAP)
    if selected_cell_type:
        umap_coords["Highlight"] = umap_coords["Cell Type"] == selected_cell_type
        if umap_type == "2D":
            cell_type_fig = px.scatter(
                umap_coords, x="UMAP1", y="UMAP2", color="Highlight",
                title=f"Cell Type: {selected_cell_type}",
                color_discrete_map={True: "red", False: "lightgray"},
                size_max=dot_size  # Set dynamic dot size
            )
        else:  # 3D UMAP
            cell_type_fig = px.scatter_3d(
                umap_coords, x="UMAP1", y="UMAP2", z="UMAP3", color="Highlight",
                title=f"Cell Type: {selected_cell_type}",
                color_discrete_map={True: "red", False: "lightgray"}
            )
            cell_type_fig.update_traces(marker=dict(size=dot_size))  # Set dynamic dot size for 3D
    else:
        if umap_type == "2D":
            cell_type_fig = px.scatter(
                umap_coords, x="UMAP1", y="UMAP2", color="Cell Type",
                title="Cell Type Annotations",
                size_max=dot_size  # Set dynamic dot size
            )
        else:  # 3D UMAP
            cell_type_fig = px.scatter_3d(
                umap_coords, x="UMAP1", y="UMAP2", z="UMAP3", color="Cell Type",
                title="Cell Type Annotations"
            )
            cell_type_fig.update_traces(marker=dict(size=dot_size))  # Set dynamic dot size for 3D

    # Plot 3: Gene Expression by Cell Type (Box Plot)
    if selected_gene:
        expression_data = pd.DataFrame({
            "Cell Type": adata.obs[CELLTYPE_COLUMN],
            "Expression": adata[:, selected_gene].X.toarray().flatten()
        })
        boxplot_fig = px.box(
            expression_data, x="Cell Type", y="Expression", color="Cell Type",
            title=f"Expression of {selected_gene} by Cell Type",
            color_discrete_sequence=px.colors.qualitative.Plotly  # Use consistent colors
        )
    else:
        boxplot_fig = px.box(
            title="Select a gene to view expression by cell type"
        )

    # Plot 4: Gene Expression by Age (Line Plot)
    if selected_gene:
        lineplot_data = pd.DataFrame({
            "Cell Type": adata.obs[CELLTYPE_COLUMN],
            "Age": adata.obs[AGE_COLUMN],
            "Expression": adata[:, selected_gene].X.toarray().flatten()
        })
        # Aggregate expression by cell type and age
        lineplot_data = lineplot_data.groupby(["Cell Type", "Age"], observed=False).mean().reset_index()
        lineplot_fig = px.line(
            lineplot_data, x="Age", y="Expression", color="Cell Type",
            title=f"Expression of {selected_gene} by Age",
            labels={"Expression": "Mean Expression", "Age": "Age"},
            color_discrete_sequence=px.colors.qualitative.Plotly  # Use consistent colors
        )
    else:
        lineplot_fig = px.line(
            title="Select a gene to view expression by age"
        )

    # Heatmap: Relative Gene Expression by Cell Type and Age
    if selected_gene:
        # Aggregate expression by cell type and age
        heatmap_data = pd.DataFrame({
            "Cell Type": adata.obs[CELLTYPE_COLUMN],
            "Age": adata.obs[AGE_COLUMN],
            "Expression": adata[:, selected_gene].X.toarray().flatten()
        })
        heatmap_data = heatmap_data.groupby(["Cell Type", "Age"], observed=False).mean().reset_index()

        # Drop rows with NaN values in the Expression column
        heatmap_data = heatmap_data.dropna(subset=["Expression"])

        # Calculate relative expression (z-score normalization)
        heatmap_data["Relative Expression"] = zscore(heatmap_data["Expression"])

        # Pivot the data for the heatmap
        heatmap_pivot = heatmap_data.pivot(index="Cell Type", columns="Age", values="Relative Expression")

        # Replace NaN values with 0 (or another appropriate value)
        heatmap_pivot = heatmap_pivot.fillna(0)

        # Create the heatmap
        heatmap_fig = go.Figure(data=go.Heatmap(
            z=heatmap_pivot.values,
            x=heatmap_pivot.columns,
            y=heatmap_pivot.index,
            colorscale="Viridis",
            colorbar=dict(title="Relative Expression (Z-Score)")
        ))
        heatmap_fig.update_layout(
            title=f"Relative Expression of {selected_gene} by Cell Type and Age",
            xaxis_title="Age",
            yaxis_title="Cell Type",
            height=600  # Adjust height as needed
        )
    else:
        heatmap_fig = go.Figure()  # Empty heatmap if no gene is selected
        heatmap_fig.update_layout(
            title="Select a gene to view relative expression by cell type and age"
        )

    return gene_fig, cell_type_fig, boxplot_fig, lineplot_fig, heatmap_fig

# Run the app
if __name__ == "__main__":
    app.run_server(host=args.host, port=args.port, debug=True)
