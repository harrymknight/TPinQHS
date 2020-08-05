function get_max_index_of_array(array, selected_chem) {
    for (let i=array.length;i>-1;i--) {
        if (array[i] < selected_chem) {
            return i;
        }
    }
}

function get_min_index_of_array(array, selected_chem) {
    for (let i=0;i<array.length;i++) {
        if (array[i] < selected_chem) {
            return i;
        }
    }
}

function zip(a, b) {
    let c = [];
    for (let i=0;i<a.length;i++) {
        c.push([a[i], b[i]]);
    }
    return c;
}

function pairwise(a) {
    let b = a.slice(0, a.length-2);
    let c = a.slice(1, a.length);
    return zip(b,c);
}

function get_electron_coords(x_k, eigenvalues) {
    let dist = 0;
    let dists = [];
    const _spacing = 50;
    let closest = [];
    let coords = {x: [], y:[]};
    for (var [[i, j], [k, l]] of zip(pairwise(x_k), pairwise(eigenvalues))) {
        dists.push(dist)
        dist += (j - i)**2 + (l - k)**2
    }
    const npoints = Math.floor(dist*_spacing)+1
    const integers = [];
    for (var m=1;m<npoints;m++) {
        integers.push(m);
    }
    const points = integers.map(integer => integer*dist/(npoints-1));
    for (var [[[i, j], [k, l]], m] of zip(zip(pairwise(dists), pairwise(x_k)), eigenvalues)) {
        if (i <= points[0] && points[0] < j) {
            closest.push([k, l, points.shift()-i, m]);
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
        landau_level: function (fig_dict, selected_chem) {

            if (!fig_dict) {
                throw "Figure data not loaded, aborting update.";
            }
            
            let fig_dict_copy = {...fig_dict};
            const chempot = fig_dict_copy['layout']['shapes'][0];
            chempot['y0'] = selected_chem;
            chempot['y1'] = selected_chem;

            let filled_ll1 = fig_dict_copy['data'][2];
            let filled_ll2 = fig_dict_copy['data'][3];
            let ll1_electrons = fig_dict_copy['data'][4];
            let ll2_electrons = fig_dict_copy['data'][5];
            let ll1_electron_coords;
            let ll2_electron_coords;
            //assume below chemical potential
            filled_ll1['x'] = null;
            filled_ll1['y'] = null;
            filled_ll2['x'] = null;
            filled_ll2['y'] = null;
            ll1_electrons['x'] = null;
            ll1_electrons['y'] = null;
            ll2_electrons['x'] = null;
            ll2_electrons['y'] = null;
            //if not then update for each case
            if (selected_chem > 0.5) {
                let ll1 = fig_dict['data'][0];
                let ll1_min_index = get_min_index_of_array(ll1['y'], selected_chem);
                let ll1_max_index = get_max_index_of_array(ll1['y'], selected_chem);
                filled_ll1['x'] = ll1['x'].slice(ll1_min_index, ll1_max_index);
                filled_ll1['y'] = ll1['y'].slice(ll1_min_index, ll1_max_index);
                ll1_electron_coords = get_electron_coords(filled_ll1['x'], filled_ll1['y'])
                ll1_electrons['x'] = ll1_electron_coords['x'];
                ll1_electrons['y'] = ll1_electron_coords['y'];
                if (selected_chem > 1.5) {
                    let ll2 = fig_dict['data'][1];
                    let ll2_min_index = get_min_index_of_array(ll2['y'], selected_chem);
                    let ll2_max_index = get_max_index_of_array(ll2['y'], selected_chem);
                    filled_ll2['x'] = ll2['x'].slice(ll2_min_index, ll2_max_index);
                    filled_ll2['y'] = ll2['y'].slice(ll2_min_index, ll2_max_index);
                    ll2_electron_coords = get_electron_coords(filled_ll2['x'], filled_ll2['y']);
                    ll2_electrons['x'] = ll2_electron_coords['x'];
                    ll2_electrons['y'] = ll2_electron_coords['y'];
                }
            }
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
                if (expectation == 'ticked') {
                    expectation_value = get_expectation_value(wavefunc)
                    x_expectation['x0'] = expectation_value;
                    x_expectation['x1'] = expectation_value;
                    x_expectation['y1'] = Math.max.apply(Math, wavefunc['y']);
                }
            }

            return fig_dict_copy;
        }
    }
});