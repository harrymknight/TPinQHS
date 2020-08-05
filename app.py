import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State, ClientsideFunction
from dash.exceptions import PreventUpdate
from flask_caching import Cache
import plotly.express as px
import plotly.graph_objects as go

import numpy as np
import pandas as pd
import wavefunc as wf
import electron as et

import time

dpoints = 400
xpoints = 800
domain = np.linspace(-5, 5, dpoints)
x_k = np.linspace(-7.5, 7.5, xpoints)
grads = [i*0.01 for i in range(11)]

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

server = app.server

app.layout = html.Div(
    [
        dbc.Row(
            [
                dbc.Col(html.Div(
                    [
                        dcc.Graph(
                            id='graph-with-slider', 
                        ),
                        dcc.Dropdown(
                            id='slider-dropdown',
                            options=[
                                {'label': 'Temperature Gradient', 'value': 'tempgrad'},
                                {'label': 'Average Chemical Potential', 'value': 'chempot'}
                            ],
                            value='tempgrad'
                        ),
                        html.Div(id='grad-slider-container', children=[
                            dcc.Slider(
                                id='grad-slider',
                                min=0,
                                max=0.1,
                                value=0,
                                marks={grad : "{:.2f}".format(grad) for grad in grads},
                                step=None
                            ),
                        ], style={'display': 'block'}),
                        html.Div(id='chem-slider-container', children=[
                            dcc.Slider(
                                id='chem-slider',
                                min=0,
                                step=0.01,
                                value=2,
                                updatemode='drag'
                            ),
                        ], style={'display': 'none'}),
                    ]
                )),
                dbc.Col(html.Div(
                    [
                        dcc.Graph(
                            id='wavefunc-graph',
                        ),
                        dcc.Dropdown(
                            id='wavefunc-dropdown',
                            options=[
                                {'label' : 'Wavefunction', 'value' : 'wavefunc'},
                                {'label' : 'Probability Density Function', 'value' : 'pdf'}
                            ],
                            value='wavefunc'
                        ),
                        html.Div(id='expectation-value-container', children=[
                            dcc.Checklist(
                                id='expectation-value',
                                options=[
                                    {'label' : 'Average Electron Position', 'value' : 'ticked'}
                                ]
                            )
                        ], style={'display' : 'none'}),
                    ]
                )),
            ],
            #no_gutters=True,
        ),
        dcc.Store(id='graph-with-slider-data', storage_type='memory'),
        dcc.Store(id='wavefunc-graph-data', storage_type='memory'),
        dcc.Store(id='wavefunc-data', storage_type='memory')
    ]
)

@app.callback(
    [Output('grad-slider-container', 'style'),
     Output('chem-slider-container', 'style'),
     Output('chem-slider', 'max')],
    [Input('slider-dropdown', 'value')],
    [State('graph-with-slider', 'figure')]
)
def choose_slider(value, figure):
    if figure is None:
        raise PreventUpdate
    maximum = max(figure['data'][1]['y'])
    if value == 'chempot':
        gsc_style = {'display': 'none'}
        csc_style = {'display': 'block'}
    else:
        gsc_style = {'display': 'block'}
        csc_style = {'display': 'none'}
    return gsc_style, csc_style, maximum

@app.callback(
    [Output('graph-with-slider-data', 'data'),
     Output('wavefunc-graph-data', 'data'),
     Output('graph-with-slider', 'hoverData')],
    [Input('grad-slider', 'value')],
    [State('graph-with-slider', 'figure'),
     State('wavefunc-graph', 'figure')])
def update_figure(selected_grad, landau_level_fig, wavefunc_fig):
    
    if landau_level_fig != None:
        chempot = landau_level_fig['layout']['shapes'][0]['y0']
    else:
        chempot = 2

    ll1 = wf.get_eigenvalues(x_k, selected_grad, 0)
    ll2 = wf.get_eigenvalues(x_k, selected_grad, 1)
    lls = pd.DataFrame(zip(x_k, ll1, ll2))
    lls.columns=['x_k', 'LL1', 'LL2']
    landau_level_fig = px.line(lls, x='x_k', y=['LL1', 'LL2'])

    landau_level_fig.for_each_trace(
        lambda trace: trace.update(
            overwrite=True,
            line=dict(color='rgba(149, 165, 166, 0.4)'),
            hovertemplate=None,
            hoverinfo='skip',
            showlegend=False,
            ) if trace.name == "LL1" or trace.name == "LL2" else (),
    )

    landau_level_fig.add_shape(
        dict(
            type="line",
            x0=x_k.min(),
            y0=chempot,
            x1=x_k.max(),
            y1=chempot,
            line=dict(
                color="Black",
                width=3,
                dash="dash"
            )
        )
    )

    landau_level_fig.add_trace(go.Scatter(
        x=None, 
        y=None, 
        name='Filled LL1', 
        mode='lines',
    ))

    landau_level_fig.add_trace(go.Scatter(
        x=None, 
        y=None, 
        name='Filled LL2', 
        mode='lines'
    ))

    landau_level_fig.add_trace(go.Scatter(
        x=None, 
        y=None, 
        name='Electrons', 
        mode='markers', 
        marker={'size' : 8}, 
        hoverinfo='skip'
    ))

    landau_level_fig.add_trace(go.Scatter(
        x=None, 
        y=None, 
        name='Electrons', mode='markers', 
        marker={'size' : 8}, 
        hoverinfo='skip'
    ))

    landau_level_fig.update_layout(margin={'l': 20, 'b': 30, 'r': 10, 't': 10}, xaxis = dict(nticks=4, range = [x_k.min(),x_k.max()],),)

    initial_wavefunc = wf.get_wavefunc(0, lls['LL1'][xpoints//2], selected_grad, 0, xpoints)
    wavefunc_fig = px.line(x=initial_wavefunc['x'], y=initial_wavefunc['y'])

    wavefunc_fig.add_shape(
        dict(
            type="line",
            x0=0,
            y0=0,
            x1=0,
            y1=0,
            line=dict(
                color="Black",
                width=3,
                dash="dash"
            )
        )
    )

    hoverData=dict(
        points=[dict(
            curveNumber=2,
            x=lls['x_k'][xpoints//2],
            y=lls['LL1'][xpoints//2],
            pointNumber=200,
            pointIndex=200,
        )]
    )

    return landau_level_fig.to_dict(), wavefunc_fig.to_dict(), hoverData

app.clientside_callback(
    ClientsideFunction(
        namespace='clientside',
        function_name='landau_level'
    ),
    Output("graph-with-slider", "figure"),
    [Input("graph-with-slider-data", "data"),
     Input("chem-slider", "value")],
)

@app.callback(
    Output('wavefunc-data', 'data'),
    [Input('graph-with-slider', 'hoverData')],
    [State('grad-slider', 'value')])
def pass_wavefunc_data(hoverData, selected_grad):
    x_0 = hoverData['points'][0]['x']
    eigenvalue = hoverData['points'][0]['y']
    level = hoverData['points'][0]['curveNumber']-2
    wavefuncData = wf.get_wavefunc(x_0, eigenvalue, selected_grad, level, xpoints)
    return dict(
            x=wavefuncData['x'].tolist(),
            y=wavefuncData['y'].tolist(),
    )

@app.callback(
    Output('expectation-value-container', 'style'),
    [Input('wavefunc-dropdown', 'value')])
def update_dropdown(dropdown):
    if dropdown == 'pdf':
        return {'display': 'block'}
    else:
        return {'display': 'none'}

app.clientside_callback(
    ClientsideFunction(
        namespace='clientside',
        function_name='wavefunc'
    ),
    Output("wavefunc-graph", "figure"),
    [Input("wavefunc-graph-data", "data"),
     Input("wavefunc-data", "data"),
     Input("wavefunc-dropdown", "value"),
     Input("expectation-value", "value")],
)

if __name__ == '__main__':
    app.run_server(debug=True)