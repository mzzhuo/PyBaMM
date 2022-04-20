
import pybamm
import numpy as np
import matplotlib.pyplot as plt


timer = pybamm.Timer()
    
#%
# parameters
# constants
# gas constant R = 8.314 [J.K-1.mol-1]
R = pybamm.constants.R
#%
# Faraday constant F = 96485.33212 [C·mol-1]
F = pybamm.constants.F

# physical parameters
c_p_max = pybamm.Parameter("Maximum lithium concentration in positive electrode [mol.m-3]")
c_n_max = pybamm.Parameter("Maximum lithium concentration in negative electrode [mol.m-3]")

# the constant coefficients in the temperature dependent diffusivities
D_p_ref = pybamm.Parameter("Reference lithium diffusivity in positive particle core [m2.s-1]")
D_s_ref = pybamm.Parameter("Reference lithium diffusivity in positive particle shell [m2.s-1]")
D_o_ref = pybamm.Parameter("Reference oxygen diffusivity in positive particle shell [m2.s-1]")
D_n_ref = pybamm.Parameter("Reference lithium diffusivity in negative electrode [m2.s-1]")

# geometry parameters
A_p = pybamm.Parameter("Positive electroactive surface area in total [m2]")
A_n = pybamm.Parameter("Negative electroactive surface area in total [m2]")
# radii of positive and negative particles
R_p = pybamm.Parameter("Positive electrode particle radius [m]")
R_n = pybamm.Parameter("Negative electrode particle radius [m]")

# chracteristic time scale as reference
t_ref = R_p ** 2 / D_p_ref

# ds / dt = k_1 - k_2 *c_o
k_1_ref = pybamm.Parameter("Forward chemical reaction (degradation) coefficient [m.s-1]")
k_2_ref = pybamm.Parameter("Reverse chemical reaction (degradation) coefficient [m4.mol-1.s-1]")

# electrochemical reaction parameters
alpha = pybamm.Parameter("Charge transfer coefficient [-]")
k_p_ref = pybamm.Parameter("Reference exchange-current density in positive particle [m2.5.mol-0.5.s-1]")
k_n_ref = pybamm.Parameter("Reference exchange-current density in negative particle [m2.5.mol-0.5.s-1]")

# 
aEne_diff_p = pybamm.Parameter("Activation energy for lithium diffusivity in positive particle core [J.mol-1]")
aEne_diff_s = pybamm.Parameter("Activation energy for lithium diffusivity in positive particle shell [J.mol-1]")
aEne_diff_o = pybamm.Parameter("Activation energy for oxygen diffusivity in positive particle shell [J.mol-1]")
aEne_diff_n = pybamm.Parameter("Activation energy for lithium diffusivity in negative electrode [J.mol-1]")
# # 
aEne_chem_reac_for = pybamm.Parameter("Activation energy for forward degradation reaction [J.mol-1]")
aEne_chem_reac_rev = pybamm.Parameter("Activation energy for reverse degradation reaction [J.mol-1]")
# 
aEne_surf_reac_p = pybamm.Parameter("Activation energy for exchange-current density in positive particle [J.mol-1]")
aEne_surf_reac_n = pybamm.Parameter("Activation energy for exchange-current density in negative particle [J.mol-1]")

# non-changing parameters
c_e  = pybamm.Parameter("Electrolyte concentration [mol.m-3]")
theta_1  = pybamm.Parameter("Electrolyte resistance coefficient 1 [Ohm]")
theta_2  = pybamm.Parameter("Electrolyte resistance coefficient 2 [Ohm.K-1]")

# cell parameters
rho = pybamm.Parameter("Cell density [kg.m-3]")
vol = pybamm.Parameter("Cell volume [m3]")
c_rho  = pybamm.Parameter("Specific heat capacity of the cell [J.kg-1.K-1]")
hA  = pybamm.Parameter("Heat transfer coefficient x cell surface area [J.s-1.K-1]")

T_ref = pybamm.Parameter("Reference temperature [K]")
T_amb = pybamm.Parameter("Ambient temperature [K]")

# simulation related parameters
t_end = pybamm.Parameter("Total simulation time (half C) [s]")
t_end_scaled = t_end / t_ref

# discharge -> positive
I_app = pybamm.Parameter("Current function [A]")

# passivation shell layer Parameters (Positive Electrode)
# 

phi = pybamm.Parameter("Ratio of maximum lithium concentration between core and shell [-]")

c_o_core = pybamm.Parameter("Constant oxygen concentration in particle core [mol.m-3]")
c_o_core_scaled = c_o_core / c_p_max 

# initial conditions
c_p_0 = pybamm.Parameter("Initial normalized lithium concentration in positive particle core [-]")
c_s_0 = pybamm.Parameter("Initial normalized lithium concentration in positive particle shell [-]")
c_o_0 = pybamm.Parameter("Initial normalized oxygen concentration in positive particle shell [-]")
c_n_0 = pybamm.Parameter("Initial normalized lithium concentration in negative electrode [-]")
s_0 = pybamm.Parameter("Initial normalized phase boundary location [-]")

# normalized critical lithium concentration in core at which phase transition may occur
c_p_thrd = pybamm.Parameter("Normalized threshold lithium concentration in positive particle core [-]")

#%
model = pybamm.BaseModel(name='DegradationModel')


# normalized variables by c_p_max
c_p = pybamm.Variable("Normalized lithium concentration in positive particle core [-]", domain="positive particle core")
c_s = pybamm.Variable("Normalized lithium concentration in positive particle shell [-]", domain="positive particle shell")
c_o = pybamm.Variable("Normalized oxygen concentration in positive particle shell [-]", domain="positive particle shell")
# normalized by c_n_max
c_n = pybamm.Variable("Normalized lithium concentration in negative particle [-]", domain="negative particle")


# Variables that depend on time only are created without a domain
Q = pybamm.Variable("Charge capacity [A.h]")
        
        
        
# surface value of concentration variables via extrapolation
c_p_surf = pybamm.boundary_value(c_p, "right")
c_s_cent = pybamm.boundary_value(c_s, "left")
c_s_surf = pybamm.boundary_value(c_s, "right")
c_o_cent = pybamm.boundary_value(c_o, "left")
c_n_surf = pybamm.boundary_value(c_n, "right")

c_p_N = pybamm.boundary_cell_value(c_p, "right")
c_s_1 = pybamm.boundary_cell_value(c_s, "left")
c_o_1 = pybamm.boundary_cell_value(c_o, "left")

# boundary cell length: node to the edge
dx_cp = pybamm.boundary_cell_length(c_p, "right")
dx_cs = pybamm.boundary_cell_length(c_s, "left")
dx_co = pybamm.boundary_cell_length(c_o, "left")

# normalized moving phase boundary location by R_p
s = pybamm.Variable("Normalized moving phase boundary location [-]")
# time derivative of phase boundary location
# c_o_cent need revisit

#%
# normalized temp. by T_ref
T = pybamm.Variable("Normalized temperature [-]")
# T = T_ref


#%
# arrhenius takes non-scaled temp. with unit [K]
D_p = D_p_ref * pybamm.arrhenius(aEne_diff_p, T_ref, T * T_ref)
D_s = D_s_ref * pybamm.arrhenius(aEne_diff_s, T_ref, T * T_ref)
D_o = D_o_ref * pybamm.arrhenius(aEne_diff_o, T_ref, T * T_ref)
D_n = D_n_ref * pybamm.arrhenius(aEne_diff_n, T_ref, T * T_ref)

# we use D_p_ref to normalize other diffusivities
D_p_scaled = D_p / D_p_ref
D_s_scaled = D_s / D_p_ref
D_o_scaled = D_o / D_p_ref
D_n_scaled = D_n / D_p_ref

k_1 = k_1_ref * pybamm.arrhenius(aEne_chem_reac_for, T_ref, T * T_ref)
k_2 = k_2_ref * pybamm.arrhenius(aEne_chem_reac_rev, T_ref, T * T_ref)

k_1_scaled = k_1 * t_ref / R_p 
k_2_scaled = k_2 * t_ref / R_p * c_p_max

s_dot = -(k_1_scaled - k_2_scaled * c_o_cent) / c_o_core_scaled * pybamm.EqualHeaviside(c_p_surf, c_p_thrd)
# pybamm.sigmoid(c_p_surf, c_p_thrd, 10)
# .evaluate().astype(float)
# s_dot = pybamm.Scalar(0) 
# matlab code
# x = -10:0.01:10;
# y = (1 + tanh(1*x)) ./ 2;
# figure
# plot(x,y,'-')

k_p = k_p_ref * pybamm.arrhenius(aEne_surf_reac_p, T_ref, T * T_ref)
k_n = k_n_ref * pybamm.arrhenius(aEne_surf_reac_n, T_ref, T * T_ref)


#%
# boundary flux
# molar flux of lithium at the boundary [mol.m-2.s-1]
J_p = -I_app / (F * A_p) 
J_n = I_app / (F * A_n) 

# scaled flux by nondimensionalization and variable change
J_p_scaled = (1 - s) * J_p * R_p / (D_p_ref * c_p_max)
J_n_scaled = J_n * R_n / (D_p_ref * c_n_max)



# spatial variables
eta = pybamm.SpatialVariable("eta", domain=["positive particle core"], coord_sys="cartesian")
chi = pybamm.SpatialVariable("chi", domain=["positive particle shell"], coord_sys="cartesian")
r_n = pybamm.SpatialVariable("r_n", domain=["negative particle"], coord_sys="spherical polar")
# r_n = pybamm.SpatialVariable("r_n", domain=["negative particle"], coord_sys="cartesian")

# rhs of equations
# u = pybamm.PrimaryBroadcastToEdges(1, ["positive particle core"])
# 
# model.rhs[c_p] = eta * s_dot / s * pybamm.div(c_p * u) + D_p_scaled / ((eta * s) ** 2) * pybamm.div(eta ** 2 * pybamm.grad(c_p))

model.rhs[c_p] = pybamm.inner(eta * s_dot / s, pybamm.grad(c_p)) + D_p_scaled / ((eta * s) ** 2) * pybamm.div(eta ** 2 * pybamm.grad(c_p))
# model.rhs[c_p] = pybamm.div(c_p * u)
# model.rhs[c_p] = eta * s_dot / s * pybamm.grad(c_p) + D_p_scaled / ((eta * s) ** 2) * pybamm.div(eta ** 2 * pybamm.grad(c_p))
# model.rhs[c_p] = D_p_scaled / ((eta * s) ** 2) * pybamm.div(eta * pybamm.grad(c_p))
# model.rhs[c_p] = pybamm.div(eta * pybamm.grad(c_p))
#%
v = pybamm.PrimaryBroadcastToEdges(1, ["positive particle shell"])
# 
model.rhs[c_s] = (
    (1 - chi) * s_dot / (1 - s) * pybamm.div(c_s * v)
    + D_s_scaled / (phi * (1 - s) ** 2 * (chi * (1 - s) + s) ** 2) * pybamm.div((chi * (1 - s) + s) ** 2 * pybamm.grad(c_s))
    )

model.rhs[c_o] = (
    (1 - chi) * s_dot / (1 - s) * pybamm.div(c_o * v)
    + D_o_scaled / ((1 - s) ** 2 * (chi * (1 - s) + s) ** 2) * pybamm.div((chi * (1 - s) + s) ** 2 * pybamm.grad(c_o))
    )

# model.rhs[c_n] = (R_p / R_n) ** 2 * (-pybamm.div(-D_n_scaled * pybamm.grad(c_n)))
# D_n_scaled has to go out
model.rhs[c_n] = (R_p / R_n) ** 2 * D_n_scaled * pybamm.div(pybamm.grad(c_n))
# model.rhs[c_n] = (R_p / R_n) ** 2 * D_n_scaled / (r_n ** 2) * pybamm.div(r_n ** 2 * pybamm.grad(c_n))

model.rhs[s] = s_dot


model.rhs[Q] = -I_app * t_ref / 3600

#%
inputs_p = {"Positive particle stoichiometry": c_s_surf}
inputs_n = {"Negative particle stoichiometry": c_n_surf}

# dUdT (derivative of OCV wit respect to Temp), entropy at reference temp.
dUdT_p = 1e-3 * pybamm.FunctionParameter("Positive electrode OCP entropy coefficient [V.T-1]", inputs_p)
dUdT_n = 1e-3 * pybamm.FunctionParameter("Negative electrode OCP entropy coefficient [V.T-1]", inputs_n)
# OCP (open-circuit potential) as a function of concentration at reference temp.
U_p_ref = pybamm.FunctionParameter("Positive electrode OCP at reference temp [V]", inputs_p)
U_n_ref = pybamm.FunctionParameter("Negative electrode OCP at reference temp [V]", inputs_n)

# electrode open-circuit potential as function of stochiometry and temp 
U_p = U_p_ref + dUdT_p * T_ref * (T - 1)
U_n = U_n_ref + dUdT_n * T_ref * (T - 1)


m_p = -I_app / (F * A_p * k_p * c_e ** 0.5 * c_p_max * (1 - c_s_surf) ** 0.5 * c_s_surf ** 0.5)
m_n = -I_app / (F * A_n * k_n * c_e ** 0.5 * c_n_max * (1 - c_n_surf) ** 0.5 * c_n_surf ** 0.5)
# overpotential [V], temp back to dimensional one
eta_p = 2 * R * (T * T_ref) / F * pybamm.log(((m_p ** 2 + 4) ** 0.5 + m_p) / 2)
eta_n = 2 * R * (T * T_ref) / F * pybamm.log(((m_n ** 2 + 4) ** 0.5 + m_n) / 2)
# temp back to dimensional one
R_cell = theta_1 + theta_2 * (T * T_ref - T_amb);

V_cell = U_p - U_n + eta_p + eta_n - I_app * R_cell

# model.rhs[T] = pybamm.Scalar(0)
model.rhs[T] = (
    R_p ** 2 / (rho * vol * c_rho * D_p_ref)
    * (-I_app * (dUdT_p - dUdT_n) * T - I_app / T_ref * (V_cell - U_p + U_n) - hA * (T - T_amb / T_ref))
    )
# model.rhs[T] = (
#     R_p ** 2 / (rho * vol * c_rho * D_p_ref)
#     * (- hA * (T - T_amb / T_ref))
#     )


#%
# boundary conditions 
# -----------------------------------------------------------------------------------
# c = pybamm.VarCombination(c_p, c_s, s, s_dot)
# model.boundary_conditions[c] = {
#     "left":  (pybamm.Scalar(0), "Neumann"),
#     "right": (-J_p_scaled / D_s_scaled, "Neumann")
#     # "right": (pybamm.Scalar(0), "Neumann")
#     }

c_shared = (
    (D_p_scaled * (1 - s) / dx_cp * c_p_N + D_s_scaled * s / dx_cs * c_s_1 
        + s_dot * s * (1 - s) * c_o_core_scaled) 
    / (D_p_scaled * (1 - s) / dx_cp + D_s_scaled * s / dx_cs 
        + s_dot * s * (1 - s) * (1 - phi))
    )

# c_shared = (
#     (D_p_scaled * (1 - s) / dx_cp * c_p_N + D_s_scaled * s / dx_cs * c_s_1)
#     / (D_p_scaled * (1 - s) / dx_cp + D_s_scaled * s / dx_cs)
#     )
# 
# c_shared = (c_p_N + c_s_1) / 2

rbc_cp = (c_shared - c_p_N) / dx_cp
lbc_cs = (c_s_1 - c_shared) / dx_cs

model.boundary_conditions[c_p] = {
    "left":  (pybamm.Scalar(0), "Neumann"),
    "right": (rbc_cp, "Neumann")
    }

model.boundary_conditions[c_s] = {
    "left":  (lbc_cs, "Neumann"),
    "right": (-J_p_scaled / D_s_scaled, "Neumann")
    # "right": (pybamm.Scalar(0), "Neumann")
    }
# -----------------------------------------------------------------------------------
c_o_b =  ((1 - s) * s_dot * dx_co * c_o_core_scaled - D_o_scaled * c_o_1) / ((1 - s) * s_dot * dx_co - D_o_scaled)
lbc_co = (c_o_1 - c_o_b) / dx_co
# or
# lbc_co = (1 - s) * s_dot * (c_o_1 - c_o_core_scaled) / (dx_co * (1 - s) * s_dot - D_o_scaled)

model.boundary_conditions[c_o] = {
    # "left":  (s_dot * (1 - s) * (c_o_core_scaled - c_o_cent) / D_o_scaled, "Neumann"),
    "left":  (lbc_co, "Neumann"), 
    "right": (pybamm.Scalar(0), "Dirichlet")
    }


model.boundary_conditions[c_n] = {
    "left":  (pybamm.Scalar(0), "Neumann"), 
    "right": (-J_n_scaled / D_n_scaled, "Neumann")
    }

#%
# initial conditions 
model.initial_conditions[c_p] = c_p_0
model.initial_conditions[c_s] = c_s_0
model.initial_conditions[c_o] = c_o_0 * (1 - chi ** 2)
model.initial_conditions[c_n] = c_n_0
model.initial_conditions[s] = s_0
model.initial_conditions[T] = T_amb / T_ref
model.initial_conditions[Q] = pybamm.Scalar(0)

model.variables = {
    "Normalized lithium concentration in positive particle core [-]": c_p,
    "Normalized lithium concentration in positive core-shell interface - core [-]":  c_p_surf, 
    "Normalized lithium concentration in positive particle shell [-]": c_s,
    "Normalized lithium concentration in positive core-shell interface - shell [-]": c_s_cent,
    "Normalized lithium concentration in positive shell surface [-]": c_s_surf, 
    "Normalized oxygen concentration in positive particle shell [-]": c_o,
    "Normalized oxygen concentration in positive core-shell interface [-]": c_o_cent, 
    "Normalized lithium concentration in negative particle [-]": c_n,
    "Normalized lithium concentration in negative particle surface [-]": c_n_surf,
    "Normalized moving phase boundary location [-]": s,
    "Normalized temperature [-]": T,
    "Normalized interface concentration [-]": c_shared,
    "Cell voltage [V]": V_cell,
    "Charge capacity [A.h]": Q
}


# model.events = [pybamm.Event("Minimum positive shell surface concentration", c_s_surf - 0.015)]
model.events = [
    pybamm.Event(
        "Maximum voltage cut-off", 
        V_cell - 4.2, 
        pybamm.EventType.TERMINATION
        )]

geometry = {"negative particle": {r_n: {"min": pybamm.Scalar(0.0), "max": pybamm.Scalar(1)}}, 
            "positive particle core": {eta: {"min": pybamm.Scalar(0.0), "max": pybamm.Scalar(1)}},
            "positive particle shell": {chi: {"min": pybamm.Scalar(0.0), "max": pybamm.Scalar(1)}}
            }

#%
# param = pybamm.ParameterValues("pybamm\myDevelop\data\parameters_abir.csv")
param = pybamm.ParameterValues("pybamm/myDevelop/data/parameters_abir.csv")
# change string type to symbol
[param.update({k: eval(v)}) for k, v in param.items() if isinstance(v, str)]

#%
param.process_model(model)
param.process_geometry(geometry)

#%
submesh_types = {"positive particle core": pybamm.MeshGenerator(pybamm.Uniform1DSubMesh),
                 "positive particle shell": pybamm.MeshGenerator(pybamm.Uniform1DSubMesh),
                 "negative particle": pybamm.MeshGenerator(pybamm.Uniform1DSubMesh)}

var_pts = {eta: 20, chi: 5, r_n: 20}
mesh = pybamm.Mesh(geometry, submesh_types, var_pts)

#%
spatial_methods = {"positive particle core": pybamm.FiniteVolume(),
                   "positive particle shell": pybamm.FiniteVolume(),
                   "negative particle": pybamm.FiniteVolume()}

#%

disc = pybamm.Discretisation(mesh, spatial_methods)

#%
# params = {
#     D_p_scaled: param.process_symbol(D_p_scaled),
#     D_s_scaled: param.process_symbol(D_s_scaled),
#     c_o_core_scaled: param.process_symbol(c_o_core_scaled),
#     phi: param.process_symbol(phi),
#     }
vars_mix = [D_p_scaled, D_s_scaled, c_o_core_scaled, phi]

#%
model_disc = disc.process_model(model, vars_mix, param, inplace=False);


#%
solver = pybamm.ScipySolver()
# solver = pybamm.CasadiSolver(mode="safe")

#%%
# t_end_scaled_value = param.process_symbol(t_end_scaled).value
t_end_scaled_value = param.evaluate(t_end_scaled)

t_eval = np.linspace(0, t_end_scaled_value, 101)
# t_eval = t_eval[0:99]
#%
solution = solver.solve(model_disc, t_eval)

print(timer.time())

#%%

c_p_out = solution["Normalized lithium concentration in positive particle core [-]"]
c_s_out = solution["Normalized lithium concentration in positive particle shell [-]"]
c_o_out = solution["Normalized oxygen concentration in positive particle shell [-]"]
c_n_out = solution["Normalized lithium concentration in negative particle [-]"]
s_out = solution["Normalized moving phase boundary location [-]"]
T_out = solution["Normalized temperature [-]"]
c_shared_out = solution["Normalized interface concentration [-]"]
c_s_surf_out = solution["Normalized lithium concentration in positive shell surface [-]"]
V_cell_out = solution["Cell voltage [V]"]
Q_out = solution["Charge capacity [A.h]"]

# f = open("c_p.txt", "w")
# f.write("Woops! I have deleted the content!")
# f.close()

#%

# times = np.linspace(0, t_end_scaled_value, 11)
times = solution.t
eta_plot = mesh["positive particle core"].edges  # numerical gridpoints
chi_plot = mesh["positive particle shell"].edges 
r_n_plot = mesh["negative particle"].edges 

s_out_value = s_out(t=times)
T_out_value = T_out(t=times)
c_shared_out_value = c_shared_out(t=times)
V_cell_out_value = V_cell_out(t=times)
Q_out_value = Q_out(t=times)
c_s_surf_out_value = c_s_surf_out(t=times)


c_p_out_value = c_p_out(t=times, x=eta_plot)
c_s_out_value = c_s_out(t=times, x=chi_plot)
c_o_out_value = c_o_out(t=times, x=chi_plot)
c_n_out_value = c_n_out(t=times, r=r_n_plot)


#%% plot
plt.figure(figsize=(8, 6))
plt.plot(Q_out_value, V_cell_out_value, 'b-')
plt.ylim(2.8,4.2)

plt.xlabel("Charge capacity [A.h]")
plt.ylabel("Voltage [V]")

#%% plot
plt.figure(figsize=(8, 6))

plt.plot(solution.t, s_out(solution.t), '-o')
plt.xlabel("Time [s]")
plt.ylabel("Moving interface location [m]")


#%%
# f = open("pybamm\myDevelop\output\c_p.txt", "w")
# np.savetxt(f, c_p_out_value.transpose(), fmt='%1.6e', delimiter=', ')
# f.close()

# f = open("pybamm\myDevelop\output\c_s.txt", "w")
# np.savetxt(f, c_s_out_value.transpose(), fmt='%1.6e', delimiter=', ')
# f.close()

# f = open("pybamm\myDevelop\output\c_o.txt", "w")
# np.savetxt(f, c_o_out_value.transpose(), fmt='%1.6e', delimiter=', ')
# f.close()

# f = open("pybamm\myDevelop\output\c_n.txt", "w")
# np.savetxt(f, c_n_out_value.transpose(), fmt='%1.6e', delimiter=', ')
# f.close()

# #%
# np.savetxt("pybamm\myDevelop\output\s.txt", np.column_stack((times, s_out_value)), fmt='%10.6f, %10.6f')
# np.savetxt("pybamm\myDevelop\output\Temp.txt", np.column_stack((times, T_out_value)), fmt='%10.6f, %10.6f')

# np.savetxt("pybamm\myDevelop\output\V_cell.txt", np.column_stack((times, V_cell_out_value)), fmt='%10.6f, %10.6f')

#%%
# Plot ------------------------------------------------------------------------
# eta_plot = mesh["positive particle core"].nodes  # numerical gridpoints
# chi_plot = mesh["positive particle shell"].nodes 
eta_plot = mesh["positive particle core"].edges  # numerical gridpoints
chi_plot = mesh["positive particle shell"].edges 
# 
# eta_plot = np.linspace(0.0, 1, 20)
# chi_plot = np.linspace(0.0, 1, 5)

plot_times = np.linspace(0, t_end_scaled_value, 11)
# plot_times = t[0::22]

plt.figure(figsize=(12, 8))
cmap = plt.get_cmap("inferno")
for i, tp in enumerate(plot_times):
    color = cmap(float(i) / len(plot_times))
    plt.plot(
        eta_plot * s_out(t=tp),
        c_p_out(tp, x=eta_plot),
        "-o",
        color=color,
        label="core" if i == 0 else ""
    )
    plt.plot(
        chi_plot * (1 - s_out(t=tp)) + s_out(t=tp),
        c_s_out(tp, x=chi_plot),
        "x",
        color=color,
        label="shell" if i == 0 else ""
    )
    
plt.ylim(0, 1.0)

plt.xlabel("r", fontsize=16)
plt.ylabel("c", fontsize=16)
plt.legend()
plt.show()




#%%
# post-process, so that the solution can be called at any time t or space r
# 
# Plot ------------------------------------------------------------------------
# eta_plot = mesh["positive particle core"].nodes  # numerical gridpoints
# chi_plot = mesh["positive particle shell"].nodes 
# r_n_plot = mesh["negative particle"].nodes 

# eta_plot = mesh["positive particle core"].edges  # numerical gridpoints
# chi_plot = mesh["positive particle shell"].edges 
# r_n_plot = mesh["negative particle"].edges 

plot_times = np.linspace(0, t_end_scaled_value, 11)

f,  ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2 ,figsize=(15,10))
cmap = plt.get_cmap("inferno")
for i, t in enumerate(plot_times):
    color = cmap(float(i) / len(plot_times))

    ax1.plot(
        eta_plot,
        c_p_out(t, x=eta_plot),
        "-o",
        color=color,
        label="core" if i == 0 else ""
    )
    ax2.plot(
        chi_plot,
        c_s_out(t, x=chi_plot),
        "-x",
        color=color,
        label="shell" if i == 0 else ""
    )
    
    ax3.plot(
        chi_plot,
        c_o_out(t, x=chi_plot),
        "-o",
        color=color,
        label="oxygen" if i == 0 else ""
    )
    
    ax4.plot(
        r_n_plot,
        c_n_out(t, r=r_n_plot),
        "-o",
        color=color,
        label="negative" if i == 0 else ""
    )
    
    
# ax1.set_xlabel("r_n", fontsize=16)
# ax1.set_ylabel("c_n", fontsize=16)
ax1.set_xlabel(r"$\eta$", fontsize=16)
ax1.set_ylabel(r"$c_p$", fontsize=16)
ax2.set_xlabel(r"$\chi$", fontsize=16)
ax2.set_ylabel(r"$c_s$", fontsize=16)
ax3.set_xlabel(r"$\chi$", fontsize=16)
ax3.set_ylabel(r"$c_o$", fontsize=16)
ax4.set_xlabel(r"$r_n$", fontsize=16)
ax4.set_ylabel(r"$c_n$", fontsize=16)

# ax1.legend()
# plt.show()


#%% plot

rr_n = np.linspace(0.001, 1, 20)
plot_times = np.linspace(0, t_end_scaled_value, 11)

plt.figure(figsize=(8, 6))
cmap = plt.get_cmap("inferno")
for i, t in enumerate(plot_times):
    color = cmap(float(i) / len(plot_times))
    
    plt.plot(rr_n, c_n_out(t, r=rr_n), "-o", color=color)
    
plt.xlabel("$r_n$")
plt.ylabel(r"$c_n$")




