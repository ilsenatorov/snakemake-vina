import json

import dash_bio as dashbio
import dash_bio.utils.ngl_parser as ngl_parser
import pandas as pd
import plotly.express as px

import dash
from dash import html
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate

receptor_path = "file:/home/ilya/Projects/docking/datasets/top10x10/resources/receptors/"
ligand_path = "file:/home/ilya/Projects/docking/datasets/top10x10/results/run_qvina/"

app = dash.Dash(__name__)

summary = pd.read_csv("../datasets/top10x10/summary.csv")
fig = px.scatter(summary, "score", "runtime", hover_name="pair", width=1000, height=600)

app.layout = html.Div(
    [
        dash.dcc.Graph(figure=fig, id="fig"),
        html.Pre(id="click-data"),
        dashbio.NglMoleculeViewer(id="molecule", height=1000, width=1000),
    ]
)


@app.callback(
    Output("click-data", "children"),
    Input("fig", "clickData"),
)
def show_molecule(clickData):
    return json.dumps(clickData, indent=2)


@app.callback(
    Output("molecule", "data"),
    Output("molecule", "molStyles"),
    Input("fig", "clickData"),
)
def return_molecule(clickData):
    if clickData is None:
        raise PreventUpdate

    pair = clickData["points"][0]["hovertext"]
    mol_id = pair.split("_")[0]

    molstyles_dict = {
        "representations": ["ball+stick"],
        "chosenAtomsColor": "white",
        "chosenAtomsRadius": 1,
        "molSpacingXaxis": 100,
    }

    data_list = [
        ngl_parser.get_data(data_path=receptor_path, pdb_id=mol_id, color="red", reset_view=True, local=False),
        ngl_parser.get_data(data_path=ligand_path, pdb_id=pair, color="blue", reset_view=True, local=False),
    ]

    return data_list, molstyles_dict


if __name__ == "__main__":
    app.run_server(debug=True)
