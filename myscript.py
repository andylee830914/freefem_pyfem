import freefem
import traceback
import numpy as np
import plotly.graph_objects as go


mesh = freefem.mesh
tri = freefem.tri

fig = go.Figure(data=[
    go.Mesh3d(
        x=mesh[:,0],
        y=mesh[:,1],
        z=mesh[:,2],
        colorbar_title='z',
        colorscale=[[0, 'gold'],
                    [0.5, 'mediumturquoise'],
                    [1, 'magenta']],
        # Intensity of each vertex, which will be interpolated and color-coded
        intensity=mesh[:,2],
        # intensitymode='cell',
        # i, j and k give the vertices of triangles
        # here we represent the 4 triangles of the tetrahedron surface
        i = tri[:,0],
        j = tri[:,1],
        k = tri[:,2],
        name='y',
        showscale=True
    )
])
from dash import Dash, dcc, html
app = Dash()
graph_style = {"flex": 1, "minwidth": 700,'height': '100vh'}
app.layout = html.Div([
    dcc.Graph(figure=fig, responsive=True, style=graph_style)
],style={"display": "flex", "flexwrap": "wrap"})

app.run_server(debug=True, use_reloader=False)  # Turn off reloader if inside Jupyter
