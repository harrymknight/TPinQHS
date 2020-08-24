function get_max_index_of_array(array, selected_chem) {
    for (let i=array.length;i>-1;i--) {
        if (array[i] < selected_chem) {
            return i;
        }
    }
    return false;
}

function get_min_index_of_array(array, selected_chem) {
    for (let i=0;i<array.length;i++) {
        if (array[i] < selected_chem) {
            return i;
        }
    }
    return false;
}

function get_gradient(x0, x1, y0, y1) {
    return (y1 - y0)/(x1 - x0)
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
