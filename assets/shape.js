function squircle (radius, npoints, a, b) {
    const interval = (2*radius)/(npoints-1);
    const xpoints = [];
    for (let i=-npoints/2;i<npoints/2 + 1;i++) {
        xpoints.push(interval*i + a);
    }
    //taking the positive root
    ypoints = xpoints.slice(1, xpoints.length - 1).map(i => (radius**4 - (i - a)**4)**(1/4) + b);
    ypoints.unshift(0);
    ypoints.push(0);
    const coords = {
        x:[],
        y:[]
    }
    //exploiting symmetry under reflection
    for (var [x, y] of zip(xpoints, ypoints)) {
        coords.x.push(x);
        coords.y.push(y);
    }
    for (var [x, y] of zip(xpoints, ypoints)) {
        coords.x.push(-x);
        coords.y.push(-y);
    }
    return coords;
}