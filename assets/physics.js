const h = 6.62607004e-34;
const _h = 1/h;
const hbar = 1.0545718e-34;
const k = 1.38064852e-23;
const _k = 1/k;
const e = 1.60217662e-19;
const _e = 1/e;
const B = 15;
const m_e = 9.10938356e-31;
const w = e*B/m_e;
const E = hbar*w;
const ll1_min = 0.5;
const ll2_min = 1.5;

function exp(u, T) {
    return Math.exp(u*E/(k*T));
}

/* Fermi-Dirac distribution */
function FD(mu_av, E_0, T) {
    return 1/(exp(mu_av - E_0, T) + 1);
}

/* equation for current at edge */
function I(mu_av, E_0, T) {
    u = mu_av - E_0;
    return e*k*T/h*Math.log(exp(u, T) + 1);
}

/* Landau level conductivity, sigma_H, is sigma * e */
function sigma(u, T, offset=0.1, width=0.16) {
    return 0.5*Math.tanh(width*u*E/(k*(T+offset))) + 0.5;
}

/* note missing prefactor of e/h */
function dIdT(mu_av, E_0, T) {
    u = mu_av - E_0;
    return (k*Math.log(exp(u, T) + 1) - (u*E/T)*(exp(u, T)/(exp(u, T) + 1)));
}

/* note missing prefactor of e/h */
function dIdU(mu_av, E_0, T) {
    u = mu_av - E_0;
    return exp(u, T)/(exp(u, T) + 1);
}

/* Ratio of delta mu and delta T given by equations expressing independent conservation of current in edge */
function dmudT(mu_av, E_0, T) {
    u = mu_av - E_0;
    return -dIdT(mu_av, E_0, T)/(dIdU(mu_av, E_0, T) + sigma(u, T));
}

/* Ratio of delta mu and delta T given by case 2 */
function dmu2dT(mu_av, E_01, E_02, T) {
    u1 = mu_av - E_01;
    if (u1 > 0.4) {
        return (-dIdT(mu_av, E_01, T) - sigma(u1, T)*dmudT(mu_av, E_02, T))/dIdU(mu_av, E_01, T);
    } else {
        return 0;
    }
}

/* Ratio of delta mu and delta T given by case 3 */
function dmu3dT(mu_av, E_01, E_02, T) {
    u1 = mu_av - E_01;
    u2 = mu_av - E_02;
    return (-dIdT(mu_av, E_01, T) - dIdT(mu_av, E_02, T))/(dIdU(mu_av, E_01, T) + dIdU(mu_av, E_02, T) + sigma(u1, T) + sigma(u2, T));
}