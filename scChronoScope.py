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

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run the Dash app.")
parser.add_argument("--host", default="127.0.0.1", help="Host IP address")
parser.add_argument("--port", type=int, default=8050, help="Port number")
args = parser.parse_args()

# Load your preprocessed data
adata = sc.read_h5ad("./BrianClark_logp1.h5ad")

# Define variables for easier modification
CELLTYPE_COLUMN = "scDeepSort"  # Column in adata.obs for cell types
AGE_COLUMN = "age"  # Column in adata.obs for age
AGE_ORDER = ["E11", "E12", "E14", "E16", "E18", "P0", "P2", "P5", "P8", "P14"]  # Age order

# Ensure the age column is categorical and ordered
adata.obs[AGE_COLUMN] = pd.Categorical(adata.obs[AGE_COLUMN], categories=AGE_ORDER, ordered=True)

# Create a Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Layout of the app
app.layout = dbc.Container([
    dbc.Row([
        dbc.Col(html.H1("Single-Cell RNA-Seq Data Explorer"), className="mb-4")
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
            dcc.Dropdown(
                id="umap-type-selector",
                options=[
                    {"label": "2D UMAP", "value": "2D"},
                    {"label": "3D UMAP", "value": "3D"}
                ],
                value="2D",  # Default to 2D
                clearable=False
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

# Callback to synchronize zoom/pan between Plot 1 and Plot 2
@app.callback(
    [Output("umap-gene-expression", "figure", allow_duplicate=True),
     Output("umap-cell-type", "figure", allow_duplicate=True)],
    [Input("umap-gene-expression", "relayoutData"),
     Input("umap-cell-type", "relayoutData")],
    [State("umap-gene-expression", "figure"),
     State("umap-cell-type", "figure")],
    prevent_initial_call=True  # Prevent this callback from running on page load
)
def sync_zoom(relayout_data1, relayout_data2, gene_fig, cell_type_fig):
    # Determine which plot triggered the callback
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update

    triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]

    # Update the zoom/pan state for both plots
    if triggered_id == "umap-gene-expression" and relayout_data1:
        if "xaxis.range" in relayout_data1:
            cell_type_fig["layout"]["xaxis"]["range"] = relayout_data1["xaxis.range"]
        if "yaxis.range" in relayout_data1:
            cell_type_fig["layout"]["yaxis"]["range"] = relayout_data1["yaxis.range"]
    elif triggered_id == "umap-cell-type" and relayout_data2:
        if "xaxis.range" in relayout_data2:
            gene_fig["layout"]["xaxis"]["range"] = relayout_data2["xaxis.range"]
        if "yaxis.range" in relayout_data2:
            gene_fig["layout"]["yaxis"]["range"] = relayout_data2["yaxis.range"]

    return gene_fig, cell_type_fig

# Run the app
if __name__ == "__main__":
    app.run_server(host=args.host, port=args.port, debug=True)
