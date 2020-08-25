import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State, ClientsideFunction
from dash.exceptions import PreventUpdate
import plotly.express as px
import plotly.graph_objects as go
import dash_defer_js_import as dji
import os

import numpy as np
import pandas as pd
import wavefunc as wf

dpoints = 400
xpoints = 1000
domain = np.linspace(-5, 5, dpoints)
x_k = np.linspace(-5.5, 5.5, xpoints)

text_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "text")
introduction_text = open(os.path.join(text_path, "introduction.md"), "r").read()
edgevsbulk_text = open(os.path.join(text_path, "edgevsbulk.md"), "r").read()
thermopower_text = open(os.path.join(text_path, "thermopower.md"), "r").read()
cases_text = open(os.path.join(text_path, "cases.md"), "r").read()

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
mathjax_script = dji.Import(src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/latest.js?config=TeX-MML-AM_CHTML")

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

server = app.server

app.index_string = """
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            <script type="text/x-mathjax-config">
            MathJax.Hub.Config({
                tex2jax: {
                inlineMath: [ ['$','$'],],
                processEscapes: true
                }
            });
            </script>
            {%renderer%}
        </footer>
    </body>
</html>
"""

app.layout = html.Div(
    [
        dcc.Markdown(introduction_text),
        dcc.Markdown(edgevsbulk_text),
        dcc.Markdown(thermopower_text),
        dcc.Markdown(cases_text),
        dbc.Row(
            [
                dbc.Col(html.Div(
                    [
                        dcc.Graph(
                            id='system-graph',
                                config={
                                    'displayModeBar': False
                                }
                        )
                    ]
                ), style = {'width': '100%', 'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})
            ]
        ),
        dcc.Store(id='system-data', storage_type='session'),
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
                                {'label': 'Temperature Gradient', 'value':'grad'},
                                {'label': 'Average Chemical Potential', 'value': 'chempot'},
                                {'label': 'System Temperature', 'value':'temp'}
                            ],
                            value='grad'
                        ),
                        html.Div(id='grad-slider-container', children=[
                            dcc.Slider(
                                id='grad-slider',
                                min=0,
                                max=0.1,
                                value=0,
                                marks={
                                    0: {'label': 'Minimum', 'style': {'color': '#77b0b1'}},
                                    0.1: {'label': 'Maximum', 'style': {'color': '#f50'}}
                                },
                                step=0.01
                            ),
                        ], style={'display': 'block'}),
                        html.Div(id='chem-slider-container', children=[
                            dcc.Slider(
                                id='chem-slider',
                                min=0,
                                marks={int(i*0.1) if not i%10 else i*0.1 : "{:.1f}".format(i*0.1) for i in range(0, 46, 5)},
                                step=0.01,
                                value=1,
                                updatemode='drag'
                            ),
                        ], style={'display': 'none'}),
                        html.Div(id='temp-slider-container', children=[
                            dcc.Slider(
                                id='temp-slider',
                                min=0.2,
                                max=0.8,
                                step=0.01,
                                marks={i*0.1 : "{:.1f}".format(i*0.1) + 'K' for i in range(11)},
                                updatemode='drag',
                                value=0.5
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
        ),
        dcc.Store(id='graph-with-slider-data', storage_type='memory'),
        dcc.Store(id='wavefunc-graph-data', storage_type='memory'),
        dcc.Store(id='wavefunc-data', storage_type='memory'),
        dbc.Row(
            [
                dbc.Col(html.Div(
                    [
                        dcc.Graph(
                            id='diffusion-thermopower-graph',
                        )
                    ]
                ))
            ]
        ),
        dbc.Row(
            [
                dbc.Col(html.Div(
                    [
                        dcc.Graph(
                            id='diffusion-thermopower-2-graph',
                        )
                    ]
                ))
            ]
        ),
        mathjax_script
    ]
)

@app.callback(
    [Output('grad-slider-container', 'style'),
     Output('chem-slider-container', 'style'),
     Output('temp-slider-container', 'style'),
     Output('chem-slider', 'max')],
    [Input('slider-dropdown', 'value')],
    [State('graph-with-slider', 'figure')]
)
def choose_slider(value, figure):
    if figure is None:
        raise PreventUpdate
    maximum = max(figure['data'][1]['y'])
    gsc_style = {'display': 'none'}
    csc_style = {'display': 'none'}
    tsc_style = {'display': 'none'}
    if value == 'chempot':
        csc_style = {'display': 'block'}
    elif value == 'grad':
        gsc_style = {'display': 'block'}
    elif value == 'temp':
        tsc_style = {'display': 'block'}
    return gsc_style, csc_style, tsc_style, maximum

@app.callback(
    [Output('graph-with-slider-data', 'data'),
     Output('graph-with-slider', 'hoverData')],
    [Input('grad-slider', 'value')],
    [State('graph-with-slider', 'figure')])
def update_landau_level_graph(selected_grad, landau_level_fig):
    
    if landau_level_fig != None:
        chempot = landau_level_fig['layout']['shapes'][0]['y0']
    else:
        chempot = 1

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
                width=2,
                dash="dash"
            )
        )
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
                width=2,
                dash="dash"
            )
        )
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
                width=2,
                dash="dash"
            )
        )
    )

    landau_level_fig.add_trace(go.Scatter(
        x=None, 
        y=None, 
        name='LL1', 
        mode='lines',
        opacity=1,
        line=dict(
            color="IndianRed",
            width=4,
        )
    ))

    landau_level_fig.add_trace(go.Scatter(
        x=None, 
        y=None, 
        name='LL2', 
        mode='lines',
        opacity=1,
        line=dict(
            color="RoyalBlue",
            width=4,
        )
    ))

    landau_level_fig.add_trace(go.Scatter(
        x=None, 
        y=None, 
        name='LL1 Electrons', 
        mode='markers', 
        marker=dict(
            size=8,
            opacity=1
        ), 
        hoverinfo='skip'
    ))

    landau_level_fig.add_trace(go.Scatter(
        x=None, 
        y=None, 
        name='LL2 Electrons', mode='markers', 
        marker=dict(
            size=8,
            opacity=1
        ),
        hoverinfo='skip'
    ))

    landau_level_fig.update_layout(
        margin={'l': 20, 'b': 30, 'r': 10, 't': 10}, 
        xaxis = dict(
            nticks=16, 
            title=r'$x_k \mathrm{ } \left[ l_B \right]$',
            range = [x_k.min(),x_k.max()],         
            gridcolor="rgb(238, 238, 238)", 
            gridwidth=1, 
            zerolinecolor= "rgb(238, 238, 238)",
            zerolinewidth=1
        ),
        yaxis = dict(
            nticks=16,
            title=r'$E_{n{,}k} \mathrm{ } \left[ \hbar \omega_c \right]$',
            range = [0,lls['LL2'].max()+0.5],
            gridcolor="rgb(238, 238, 238)", 
            gridwidth=1, 
            zerolinecolor= "rgb(238, 238, 238)",
            zerolinewidth=1
        ),
        showlegend=True,
        legend_title_text=None,
        legend=dict(
            orientation="h",
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.8,
            bgcolor='rgba(0,0,0,0)',
            font=dict(
                size=14,
                color="Black"
            ),
        ),
        plot_bgcolor='rgba(0,0,0,0)'
    )

    hoverData=dict(
        points=[dict(
            curveNumber=2,
            x=lls['x_k'][xpoints//2],
            y=lls['LL1'][xpoints//2],
            pointNumber=xpoints//2,
            pointIndex=xpoints//2,
        )]
    )

    return landau_level_fig.to_dict(), hoverData

app.clientside_callback(
    ClientsideFunction(
        namespace='clientside',
        function_name='landau_level'
    ),
    Output("graph-with-slider", "figure"),
    [Input("graph-with-slider-data", "data"),
     Input("chem-slider", "value"),
     Input("temp-slider", "value")],
    [State("grad-slider", "value")],
)

@app.callback(
    Output('wavefunc-graph-data', 'data'),
    [Input('graph-with-slider-data', 'data')],
    [State('grad-slider', 'value')],
)
def update_wavefunc_graph(fig_dict, selected_grad):
    eigenvalue = fig_dict['data'][0]['y'][xpoints//2]
    x_0 = fig_dict['data'][0]['x'][xpoints//2]

    initial_wavefunc = wf.get_wavefunc(x_0, eigenvalue, selected_grad, 0, xpoints)
    wavefunc_fig = go.Figure(go.Scatter(
        x=initial_wavefunc['x'],
        y=initial_wavefunc['y'],
        line_color='DarkSlateGrey'
    ))

    wavefunc_fig.add_shape(
        dict(
            type="line",
            x0=0,
            y0=0,
            x1=0,
            y1=0,
            line=dict(
                color="Black",
                width=2,
                dash="dash"
            )
        )
    )

    wavefunc_fig.update_layout(
        margin={'l': 20, 'b': 30, 'r': 10, 't': 10},
        xaxis=dict(
            title=r'$\textrm{Displacement from } x_k \mathrm{ } \left[ l_B \right]$',
            gridcolor="rgb(238, 238, 238)", 
            gridwidth=1, 
            zerolinecolor= "rgb(238, 238, 238)",
            zerolinewidth=1
        ),
        yaxis=dict(
            title=r'$\psi \mathrm{ } \left[ {l_B}^{-\frac{1}{2}} \right]$',
            gridcolor="rgb(238, 238, 238)", 
            gridwidth=1, 
            zerolinecolor= "rgb(238, 238, 238)",
            zerolinewidth=1
        ),
        plot_bgcolor='rgba(0,0,0,0)'
    )

    return wavefunc_fig.to_dict()

@app.callback(
    Output('wavefunc-data', 'data'),
    [Input('graph-with-slider', 'hoverData')],
    [State('grad-slider', 'value')])
def pass_wavefunc_hoverData(hoverData, selected_grad):
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

@app.callback(
    [Output("system-data", "data")],
    [Input("graph-with-slider-data", "data")])
def pass_system_data(fig_dict):
    return dict(
        min=dict(zip(range(1,3), [min(fig_dict['data'][0]['y']), min(fig_dict['data'][1]['y'])])),
        max=dict(zip(range(1,3), [max(fig_dict['data'][0]['y']), max(fig_dict['data'][1]['y'])])),
    ),

app.clientside_callback(
    ClientsideFunction(
        namespace='clientside',
        function_name='system'
    ),
    Output("system-graph", "figure"),
    [Input("system-data", "data"),
     Input("temp-slider", "value"),
     Input("chem-slider", "value")],
    [State("grad-slider", "value")]    
)

app.clientside_callback(
    ClientsideFunction(
        namespace='clientside',
        function_name='diffusion_thermopower'
    ),
    Output("diffusion-thermopower-graph", "figure"),
    [Input("temp-slider", "value"),
     Input("chem-slider", "value")],
    [State("diffusion-thermopower-graph", "figure")]
)

app.clientside_callback(
    ClientsideFunction(
        namespace='clientside',
        function_name='diffusion_thermopower_2'
    ),
    Output("diffusion-thermopower-2-graph", "figure"),
    [Input("temp-slider", "value"),
     Input("chem-slider", "value")],
    [State("diffusion-thermopower-2-graph", "figure")]
)

if __name__ == '__main__':
    app.run_server(debug=True)