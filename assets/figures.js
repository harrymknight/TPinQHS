function get_offset_coords(x_k, eigenvalues, selected_grad) {
    let cutoff = [];
    let found = false;
    let sigdigits = 3;
    let index = 0;
    for (var [[i, j], [k, l]] of zip(pairwise(x_k), pairwise(eigenvalues))) {
        current_grad = Math.round(get_gradient(i, j, k, l) * 10**sigdigits)/(10**sigdigits);
        if (cutoff.length != 2) {
            if (!found && (selected_grad - current_grad < 0.02)) {
                cutoff.push({x:j, y:l, point:index});
                found = true;
            }
            if (found && (current_grad - selected_grad > 0.02)) {
                cutoff.push({x:i, y:k, point:index});
            }
        }
        index++;
    }
    return cutoff;
}

function fill_landau_level(landau_level, offset_coords, filled_landau_level, selected_chem, left_selected_chem, right_selected_chem) {
    landau_level_left = landau_level['y'].slice(0, offset_coords[0]['point']);
    landau_level_centre = landau_level['y'].slice(offset_coords[0]['point'], offset_coords[1]['point']);
    landau_level_right = landau_level['y'].slice(offset_coords[1]['point'], landau_level.length);
    let landau_level_min_index = get_min_index_of_array(landau_level_left, left_selected_chem);
    let landau_level_max_index = get_max_index_of_array(landau_level_right, right_selected_chem);
    if (Number.isInteger(landau_level_min_index) && Number.isInteger(landau_level_max_index)) {
        landau_level_max_index += offset_coords[1]['point'];
        filled_landau_level['x'] = landau_level['x'].slice(landau_level_min_index, landau_level_max_index);
        filled_landau_level['y'] = landau_level['y'].slice(landau_level_min_index, landau_level_max_index);
    }
    else if (landau_level_min_index === false || landau_level_max_index === false) {
        if (landau_level_min_index === false) {
            landau_level_min_index = get_min_index_of_array(landau_level_centre, selected_chem);
            landau_level_min_index += offset_coords[0]['point'];
        }
        if (landau_level_max_index === false) {
            landau_level_max_index = get_max_index_of_array(landau_level_centre, selected_chem);
            landau_level_max_index += offset_coords[0]['point'];
        }
        filled_landau_level['x'] = landau_level['x'].slice(landau_level_min_index, landau_level_max_index);
        filled_landau_level['y'] = landau_level['y'].slice(landau_level_min_index, landau_level_max_index);
        if (landau_level_min_index === false || landau_level_max_index === false) {
            filled_landau_level['x'] = null;
            filled_landau_level['y'] = null;
        }
    }
}

function get_electron_coords(x_k, eigenvalues, chempot, chem_offset, selected_temp) {
    let dist = 0;
    let dists = [];
    const _spacing = 400;
    let closest = [];
    let coords = {x: [], y:[]};
    for (var [[i, j], [k, l]] of zip(pairwise(x_k), pairwise(eigenvalues))) {
        dists.push(dist);
        dist += (j - i)**2 + (l - k)**2;
    }
    const npoints = Math.floor(dist*_spacing)+1
    points = [];
    for (var m=1;m<npoints;m++) {
        points.push(m*dist/(npoints-1))
    }
    const selected_chem = chempot['y0'];
    let offset;
    for (var [[[i, j], [k, l]], m] of zip(zip(pairwise(dists), pairwise(x_k)), eigenvalues)) {
        if (i <= points[0] && points[0] < j) {
            if (k < chempot['x0']) {
                offset = -chem_offset;
            } else if (k > chempot['x1']) {
                offset = chem_offset;
            } else
                offset = 0;
            if (Math.random() < FD(m+i-points[0], selected_chem + offset, selected_temp)) {
                closest.push([k, l, points.shift()-i, m]);
            }
            else {
                points.shift();
            }
        }
    }
    for (var [j, k, l, m] of closest) {
        coords['x'].push(k+(k-j)*0.5);
        coords['y'].push(m - l);
    }
    return coords
}

function get_expectation_value(wavefunc) {
    let sum = 0;
    for (var [[i, j], k] of zip(pairwise(wavefunc['x']), wavefunc['y'])) {
        sum += k * i * (j-i);
    }
    return sum;
}

window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        landau_level: function (fig_dict, selected_chem, selected_temp, selected_grad) {

            if (!fig_dict) {
                throw "landau level data not loaded, aborting update.";
            }
            
            let fig_dict_copy = {...fig_dict};
            
            let lls = {
                '1': {
                    level:fig_dict_copy['data'][0],
                    fill:fig_dict_copy['data'][2],
                    electrons:fig_dict_copy['data'][4],
                    min:Math.min.apply(Math, fig_dict_copy['data'][0]['y']),
                    max:Math.max.apply(Math, fig_dict_copy['data'][0]['y'])
                },
                '2':{
                    level:fig_dict_copy['data'][1],
                    fill:fig_dict_copy['data'][3],
                    electrons:fig_dict_copy['data'][5],
                    min:Math.min.apply(Math, fig_dict_copy['data'][1]['y']),
                    max:Math.max.apply(Math, fig_dict_copy['data'][1]['y'])
                }
            }

            const chempot = fig_dict_copy['layout']['shapes'][0];
            const left_chempot = fig_dict_copy['layout']['shapes'][1];
            const right_chempot = fig_dict_copy['layout']['shapes'][2];
            const offset_coords = get_offset_coords(lls['1']['level']['x'], lls['1']['level']['y'], selected_grad);
            const chem_offset = (offset_coords[1]['y'] - offset_coords[0]['y'])/2;
            const left_selected_chem = selected_chem - chem_offset;
            const right_selected_chem = selected_chem + chem_offset;
            chempot['y0'] = selected_chem;
            chempot['y1'] = selected_chem;
            chempot['x0'] = offset_coords[0]['x'];
            chempot['x1'] = offset_coords[1]['x'];
            left_chempot['y0'] = left_selected_chem;
            left_chempot['y1'] = left_selected_chem;
            left_chempot['x1'] = offset_coords[0]['x'];
            right_chempot['y0'] = right_selected_chem;
            right_chempot['y1'] = right_selected_chem;
            right_chempot['x0'] = offset_coords[1]['x'];
            fig_dict_copy['layout']['annotations'] = [];
            fig_dict_copy['layout']['annotations'].push(
                    {
                        x:(left_chempot.x1 - left_chempot.x0)/2 + left_chempot.x0,
                        y:left_chempot.y0+0.2,
                        showarrow:false,
                        text:"$ \\LARGE{\\mu_l} $",
                        xref:"x",
                        yref:"y",
                        font: {
                            color: "Black"
                        }
                    },
                    {
                        x:(right_chempot.x1-right_chempot.x0)/2 + right_chempot.x0,
                        y:right_chempot.y0+0.2,
                        showarrow:false,
                        text:"$ \\LARGE{\\mu_r} $",
                        xref:"x",
                        yref:"y",
                        font: {
                            color: "Black"
                        }
                    },
                    {
                        x:(chempot.x1-chempot.x0)/2 + chempot.x0,
                        y:chempot.y0+0.2,
                        showarrow:false,
                        text:"$ \\LARGE{\\mu_{av}} $",
                        xref:"x",
                        yref:"y",
                        font: {
                            color: "Black"
                        }
                    }
            )

            let temp;
            fill_landau_level(lls['1']['level'], offset_coords, lls['1']['fill'], selected_chem, left_selected_chem, right_selected_chem);
            //Destructuring assignment gives identical output but does not update electron trace (Mutating an inplace object does not force an update)
            temp = get_electron_coords(lls['1']['level']['x'], lls['1']['level']['y'], chempot, chem_offset, selected_temp);
            lls['1']['electrons']['x'] = temp.x;
            lls['1']['electrons']['y'] = temp.y;
            fill_landau_level(lls['2']['level'], offset_coords, lls['2']['fill'], selected_chem, left_selected_chem, right_selected_chem);
            temp = get_electron_coords(lls['2']['level']['x'], lls['2']['level']['y'], chempot, chem_offset, selected_temp);
            lls['2']['electrons']['x'] = temp.x;
            lls['2']['electrons']['y'] = temp.y;

            return fig_dict_copy;
        },

        wavefunc: function (fig_dict, wavefuncData, dropdown, expectation) {

            if (!fig_dict) {
                throw "wavefunction data not loaded, aborting update.";
            }
            
            let fig_dict_copy = {...fig_dict};
            let wavefunc = fig_dict_copy['data'][0];
            wavefunc['x'] = wavefuncData['x'];
            wavefunc['y'] = wavefuncData['y'];
            let x_expectation = fig_dict_copy['layout']['shapes'][0];
            x_expectation['x0'] = 0;
            x_expectation['x1'] = 0;
            x_expectation['y1'] = 0;
            fig_dict_copy.layout.yaxis.title = '$\\psi \\mathrm{ } \\left[ {l_B}^{-\\frac{1}{2}} \\right]$';

            if (wavefunc['x'][0] > -5) {
                wavefunc['x'].unshift(-5);
                wavefunc['y'].unshift(0);
            } else {
                wavefunc['x'].push(5);
                wavefunc['y'].push(0);
            }
            
            if (dropdown == "pdf") {
                let wavefunc = fig_dict_copy['data'][0];
                wavefunc['y'] = Array.from(wavefunc['y'], i => i**2);
                fig_dict_copy.layout.yaxis.title = '$|\\psi|^2 \\mathrm{ } \\left[ {l_B}^{-1} \\right]$';
                if (expectation == 'ticked') {
                    expectation_value = get_expectation_value(wavefunc)
                    x_expectation['x0'] = expectation_value;
                    x_expectation['x1'] = expectation_value;
                    x_expectation['y1'] = Math.max.apply(Math, wavefunc['y']);
                }
            }

            return fig_dict_copy;
        },

        system: function(chempot, selected_temp, selected_chem, selected_grad) {
            
            let coords = [];
            coords.push(squircle(4.5, 800, 0, 0));
            coords.push(squircle(4, 800, 0, 0));
            fig_dict = {
                data: [
                    {
                        type: "scatter",
                        name: "LL1",
                        x: coords[0]['x'],
                        y: coords[0]['y'],
                        line: {
                            color:"IndianRed",
                            width:3
                        }
                    },
                    {
                        type: "scatter",
                        name: "LL2",
                        x: coords[1]['x'],
                        y: coords[1]['y'],
                        line: {
                            color:"RoyalBlue",
                            width:3
                        }
                    },
                    ],
                layout: {
                    shapes: [
                        {
                            type:"rect",
                            x0:-5,
                            y0:-5,
                            x1:5,
                            y1:5,
                            line:{
                                color: "DarkSlateGrey",
                            },
                        },
                        {
                            type:"line",
                            x0:-5,
                            y0:0,
                            x1:5,
                            y1:0,
                            line:{
                                color: "Black",
                                dash: "dash"
                            },
                        },
                        {
                            type:"circle",
                            x0:-0.8,
                            y0:-3.0,
                            x1:0.8,
                            y1:-1.4,
                            line:{
                                color: "Black",
                            },
                        },
                        {
                            type:"circle",
                            x0:-0.2,
                            y0:-2.4,
                            x1:0.2,
                            y1:-2.0,
                            line:{
                                color: "Black",
                            },
                            fillcolor: "Black",
                        },
                    ],
                    annotations: [
                        {
                            x:0.61,
                            y:0.266,
                            showarrow:false,
                            text:"$ \\LARGE{\\mathrm{B}} $",
                            xref:"paper",
                            yref:"paper",
                            font: {
                                color: "Black"
                            }
                        },
                        {
                            x:0.5,
                            y:0.766,
                            showarrow:false,
                            text:"$ \\LARGE{\\nabla \\mathrm{T}} $",
                            xref:"paper",
                            yref:"paper",
                            font: {
                                color: "Black"
                            }
                        },
                        {
                            x: -1,
                            y: 2,
                            ax: 1,
                            ay: 2,
                            xref: 'x',
                            yref: 'y',
                            axref: 'x',
                            ayref: 'y',
                            text: '',
                            showarrow: true,
                            arrowhead: 7,
                            arrowhead: 2,
                            arrowcolor: 'black'
                        }

                    ],
                    margin: {'l': 20, 'b': 30, 'r': 10, 't': 10},
                    xaxis: {nticks: 16, range: [-5.1,5.1], showgrid: false, zeroline: false, showticklabels: false},
                    yaxis: {nticks: 16, range: [-5.1,5.1], showgrid: false, zeroline: false, showticklabels: false},
                    paper_bgcolor: 'rgba(0,0,0,0)',
                    plot_bgcolor: 'rgba(0,0,0,0)',
                    showlegend: true,
                    legend_title_text: false,
                    legend: {
                        orientation: "v",
                        yanchor: "top",
                        y: 0.98,
                        xanchor: "left",
                        x: 0.86,
                        bgcolor: 'rgba(0,0,0,0)',
                        font: {
                            size: 14,
                            color: "Black"
                        },
                    },
                    hovermode: 'closest',
                    height: 600,
                    width: 600,
                }
            }

            triggered = dash_clientside.callback_context.triggered.map(t => t["prop_id"].split('.')[0]);
            if (triggered != selected_chem) {
                integers = []
                for (let i=0;i<100;i++) {
                    integers.push(i)
                }
                for (let [i,j] of pairwise(integers)) {
                    fig_dict['layout']['shapes'].push(
                        {
                            type:"rect",
                            x0:-5+i*0.1,
                            y0:-5,
                            x1:-5+j*0.1,
                            y1:5,
                            line:{
                                color: "Tomato",
                                width: 0.8,
                            },
                            layer:'below',
                            fillcolor:'Tomato',
                            opacity:(selected_grad*10)*((100-i)*0.01)
                        }
                    )
                }
            }
            let current = {
                max: {
                    '1': I(chempot['max']['1'], chempot['min']['1'], selected_temp),
                    '2': I(chempot['max']['2'], chempot['min']['2'], selected_temp)
                }
            }
            fig_dict['data'][0]['opacity'] = I(selected_chem, chempot['min']['1'], selected_temp)/current['max']['1'];
            fig_dict['data'][1]['opacity'] = I(selected_chem, chempot['min']['2'], selected_temp)/current['max']['2'];

            return fig_dict;
        },

        diffusion_thermopower: function (selected_temp, selected_chem, current_fig_dict) {
            const xaxis = [];
            const points = 400;
            const start = 0;
            const end = 2;
            const interval = (end-start)/(points-1);
            for (let i=0;i<points;i++) {
                xaxis.push(start + i*interval)
            }
            thermopower = {};
            if (current_fig_dict) {
                thermopower['1'] = [...current_fig_dict['data'][0]['y']];
                thermopower['2'] = [...current_fig_dict['data'][1]['y']];
                thermopower['1,2'] = [...current_fig_dict['data'][2]['y']];
                thermopower['1+2'] = [...current_fig_dict['data'][3]['y']];
                thermopower['limit'] = {
                    lower: current_fig_dict['layout']['shapes'][0]['y0'],
                    upper: current_fig_dict['layout']['shapes'][0]['y1']
                }
            } else {
                thermopower['1'] = Array(points);
                thermopower['2'] = Array(points);
                thermopower['1,2'] = Array(points);
                thermopower['1+2'] = Array(points);
                thermopower['limit'] = {
                    lower: 0,
                    upper: 0
                }
            }
            fig_dict = {
                data: [
                    {
                        type: "scatter",
                        name: "1",
                        x: xaxis,
                        y: thermopower['1'],
                        line: {
                            color:"IndianRed",
                            width:2
                        }
                    },
                    {
                        type: "scatter",
                        name: "2",
                        x: xaxis,
                        y: thermopower['2'],
                        line: {
                            color:"RoyalBlue",
                            width:2
                        }
                    },
                    {
                        type: "scatter",
                        name: "1,2",
                        x: xaxis,
                        y: thermopower['1,2'],
                        line: {
                            color:"ForestGreen",
                            width:2
                        }
                    },
                    {
                        type: "scatter",
                        name: "1+2",
                        x: xaxis,
                        y: thermopower['1+2'],
                        line: {
                            color: "DarkMagenta",
                            width: 2,
                            dash: "dash"
                        }
                    },
                    ],
                layout: {
                    shapes: [
                        {
                            type:"line",
                            x0:selected_chem,
                            y0:thermopower.limit.lower,
                            x1:selected_chem,
                            y1:thermopower.limit.upper,
                            line:{
                                color: "Black",
                                width: 2,
                                dash: "dash"
                            }
                        }
                    ],
                    margin: {'t': 10},
                    xaxis: {nticks:16, range:[0, 2], title: '$\\mu_{av} \\mathrm{ } \\left[ \\hbar\\omega_{c} \\right]$'},
                    yaxis: {nticks:16, title: '$S_{xx } \\left[ k_{B} \\right]$'},
                    showlegend: true,
                    legend_title_text: false,
                    legend: {
                        orientation: "h",
                        yanchor: "top",
                        y: 0.99,
                        xanchor: "left",
                        x: 0.8,
                        bgcolor: 'rgba(0,0,0,0)',
                        font: {
                            size: 14,
                            color: "Black"
                        }
                    }
                }
            }

            triggered = dash_clientside.callback_context.triggered.map(t => t["prop_id"].split('.')[0]);
            if (triggered != 'chem-slider') {
                for (let i=0;i<points;i++) {
                    fig_dict['data'][0]['y'][i] = dmudT(xaxis[i], ll1_min, selected_temp)*_k;
                    fig_dict['data'][1]['y'][i] = dmudT(xaxis[i], ll2_min, selected_temp)*_k;
                    fig_dict['data'][2]['y'][i] = dmu2dT(xaxis[i], ll1_min, ll2_min, selected_temp)*_k;
                    fig_dict['data'][3]['y'][i] = dmu3dT(xaxis[i], ll1_min, ll2_min, selected_temp)*_k;
                }
                fig_dict['layout']['shapes'][0]['y1'] = Math.max.apply(Math, fig_dict['data'][2]['y']) + 0.1;
                fig_dict['layout']['shapes'][0]['y0'] = Math.min.apply(Math, fig_dict['data'][0]['y']) - 0.1;
            }

            return fig_dict;
        },

        diffusion_thermopower_2: function (selected_temp, selected_chem, current_fig_dict) {
            const xaxis = [];
            const points = 400;
            const start = 0.2;
            const end = 0.8;
            const interval = (end-start)/(points-1);
            for (let i=0;i<points;i++) {
                xaxis.unshift(start + i*interval);
            }
            thermopower = {};
            if (current_fig_dict) {
                thermopower['1'] = [...current_fig_dict['data'][0]['y']];
                thermopower['2'] = [...current_fig_dict['data'][1]['y']];
                thermopower['1,2'] = [...current_fig_dict['data'][2]['y']];
                thermopower['1+2'] = [...current_fig_dict['data'][3]['y']];
                thermopower['limit'] = {
                    lower: current_fig_dict['layout']['shapes'][0]['y0'],
                    upper: current_fig_dict['layout']['shapes'][0]['y1']
                }
            } else {
                thermopower['1'] = Array(points);
                thermopower['2'] = Array(points);
                thermopower['1,2'] = Array(points);
                thermopower['1+2'] = Array(points);
                thermopower['limit'] = {
                    lower: 0,
                    upper: 0
                }
            }
            fig_dict = {
                data: [
                    {
                        type: "scatter",
                        name: "1",
                        x: xaxis,
                        y: thermopower['1'],
                        line: {
                            color:"IndianRed",
                            width:2
                        }
                    },
                    {
                        type: "scatter",
                        name: "2",
                        x: xaxis,
                        y: thermopower['2'],
                        line: {
                            color:"RoyalBlue",
                            width:2
                        }
                    },
                    {
                        type: "scatter",
                        name: "1,2",
                        x: xaxis,
                        y: thermopower['1,2'],
                        line: {
                            color:"ForestGreen",
                            width:2
                        }
                    },
                    {
                        type: "scatter",
                        name: "1+2",
                        x: xaxis,
                        y: thermopower['1+2'],
                        line: {
                            color: "DarkMagenta",
                            width: 2,
                            dash: "dash"
                        }
                    },
                    ],
                layout: {
                    shapes: [
                        {
                            type:"line",
                            x0:E*_k/selected_temp,
                            y0:thermopower.limit.lower,
                            x1:E*_k/selected_temp,
                            y1:thermopower.limit.upper,
                            line:{
                                color: "Black",
                                width: 2,
                                dash: "dash"
                            }
                        }
                    ],
                    margin: {'t': 10},
                    xaxis: {nticks:16, title: '$\\frac{\\hbar \\omega_c}{k_{B} \\mathrm{T}}$'},
                    yaxis: {nticks:16, title: '$S_{xx } \\left[ k_{B} \\right]$'},
                    showlegend: true,
                    legend_title_text: false,
                    legend: {
                        orientation: "h",
                        yanchor: "top",
                        y: 0.99,
                        xanchor: "left",
                        x: 0.8,
                        bgcolor: 'rgba(0,0,0,0)',
                        font: {
                            size: 14,
                            color: "Black"
                        }
                    }
                }
            }

            triggered = dash_clientside.callback_context.triggered.map(t => t["prop_id"].split('.')[0]);
            if (triggered != 'temp-slider') {
                if (selected_chem <= 2) {
                    for (let i=0;i<points;i++) {
                        fig_dict['data'][0]['y'][i] = dmudT(selected_chem, ll1_min, xaxis[i])*_k;
                        fig_dict['data'][1]['y'][i] = dmudT(selected_chem, ll2_min, xaxis[i])*_k;
                        fig_dict['data'][2]['y'][i] = dmu2dT(selected_chem, ll1_min, ll2_min, xaxis[i])*_k;
                        fig_dict['data'][3]['y'][i] = dmu3dT(selected_chem, ll1_min, ll2_min, xaxis[i])*_k;
                    }
                    maxandmin = []
                    for (let i=0;i<4;i++) {
                        maxandmin.push(fig_dict['data'][i]['y'][0])
                    }
                    fig_dict['layout']['shapes'][0]['y1'] = Math.max.apply(Math, maxandmin);
                    fig_dict['layout']['shapes'][0]['y0'] = Math.min.apply(Math, maxandmin);
                }
            }

            for (let i=0;i<points;i++) {
                xaxis[i] = E*_k/xaxis[i]
            }

            return fig_dict;
        }
    }
});