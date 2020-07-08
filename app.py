import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import plotly.express as px

import array
import numpy as np
import pandas as pd
import wavefunc as wf

import time

dpoints = 400
xpoints = 400
domain = np.array([i for i in range(dpoints)])
x_k = np.linspace(-7.5, 7.5, xpoints)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = html.Div(
    [
        dbc.Row(
            [
                dbc.Col(html.Div([dcc.Graph(
                    id='graph-with-slider', 
                    hoverData={'points': [{'curveNumber': 0, 'pointNumber': 200, 'pointIndex': 200, 'x': 200, 'y': 0.5}]}
                )])),
                dbc.Col(html.Div([
                    dcc.Dropdown(
                        id='func-dropdown',
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
                            {'label' : 'Average Wavefunction Position', 'value' : 'x-expectation'}
                        ])
                    ], style={'display' : 'none'}),
                    dcc.Graph(
                        id='func',
                        figure=px.line(wf.get_wavefunc(0, 0.5, 0, 0, xpoints), width=900, height=900),
                    )])),
            ]
        ),
        dbc.Row(dbc.Col(html.Div([dcc.Slider(
            id='grad-slider',
            min=0,
            max=0.1,
            value=0,
            marks={i*0.01 : "{:.2f}".format(i*0.01) for i in range(11)},
            step=None,
        )])))
    ]
)


@app.callback(
    Output('graph-with-slider', 'figure'),
    [Input('grad-slider', 'value')])
def update_figure(selected_grad):
    lls = zip(wf.get_eigenvalues(x_k, selected_grad, 0), wf.get_eigenvalues(x_k, selected_grad, 1))

    fig = px.line(lls, width=900, height=900)

    fig.update_layout(transition_duration=500)

    return fig


@app.callback(
    [Output('func', 'figure'), Output('expectation-value-container', 'style')],
    [Input('graph-with-slider', 'hoverData'),
     Input('func-dropdown', 'value'), 
     Input('expectation-value', 'value')],
    [State('grad-slider', 'value')])
def plot_func(hoverData, func, value, grad):
    x_0 = x_k[hoverData['points'][0]['pointIndex']]
    eigenvalue = hoverData['points'][0]['y']
    level = hoverData['points'][0]['curveNumber']
    wavefunc = wf.get_wavefunc(x_0, eigenvalue, grad, level, xpoints)
    toggle = {'display' : 'none'}
    if func == 'pdf':
        wavefunc = [i**2 for i in wavefunc]
        toggle = {'display' : 'block'}
    fig = px.line(wavefunc, width=800, height=800)
    if value and toggle['display'] == 'block':
        spacing = (10+0.7*x_0)/dpoints
        x_expectation = (sum(i*j*spacing for i, j in zip(wavefunc, x_k)))
        print(x_expectation)
        fig = fig.add_shape(
        dict(
            type="line",
            x0=x_expectation,
            y0=0,
            x1=x_expectation,
            y1=1,
            line=dict(
                color="RoyalBlue",
                width=3
            )))
    return fig, toggle


if __name__ == '__main__':
    app.run_server(debug=True)