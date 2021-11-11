import pybamm
import numpy as np
# import os
import matplotlib.pyplot as plt

#
# %%

# model = pybamm.lithium_ion.BasicSPM()
model = pybamm.lithium_ion.BasicSPMTest({"PE degradation": "yes"})


#%%
geometry = model.default_geometry

param = model.default_parameter_values
param.process_model(model)
param.process_geometry(geometry)


#%
submesh_types = model.default_submesh_types
var_pts = model.default_var_pts

mesh = pybamm.Mesh(geometry, submesh_types, var_pts)

spatial_methods = model.default_spatial_methods

disc = pybamm.Discretisation(mesh, spatial_methods)

#%
model_disc = disc.process_model(model, inplace=False)

#%%
# solver = pybamm.ScipySolver()
solver = pybamm.CasadiSolver()

t_eval = np.linspace(0, 7300, 100)

solution = solver.solve(model_disc, t_eval)



#%%

# sim = pybamm.Simulation(model)
# sim.solve([0, 3600])

#%%
output_variables = [
    "X-averaged negative particle concentration",
    #---------------------------------------------------------
    # "Current [A]",
    "Terminal voltage",
    # "X-averaged positive electrode temperature [K]",
    # "X-averaged time derivative of moving phase boundary location",
    # "c_c_N",
    # "c_s_1",
    # "D_c",
    # "D_s",
    # "shared concentration"
    #---------------------------------------------------------
    "X-averaged positive core concentration",
    "X-averaged positive shell concentration",
    # "X-averaged positive particle concentration",
]

plot = pybamm.QuickPlot(
    solution,
    output_variables,
    # time_unit="seconds",
    # spatial_unit="um",
)

plot.dynamic_plot()

# sim.plot(output_variables)





#%%

cn = model.variables["X-averaged negative particle concentration"]
cn_rav = pybamm.r_average(cn)
disc_cnrav = disc.process_symbol(cn_rav)




#%%
from scipy.sparse import csr_matrix, issparse

# eqn = model.rhs[list(model.rhs.keys())[2]]

eqn = model.tesym

#%%

zzz = disc.process_symbol(eqn)

#%%
left, right = eqn.children

disc_left = disc.process_symbol(left)
disc_right = disc.process_symbol(right)



#%%
# zzz.children[1].children[1].children[0].entries.toarray()
zzz.children[1].children[1].children[0].evaluate()













#%%
solution = sim.solution
c_c_out = solution["X-averaged positive core concentration"]
c_s_out = solution["X-averaged positive shell concentration"]

#%
timescale = solution.timescale_eval
lengthscale_core = solution.length_scales_eval["positive core"]
lengthscale_shell = solution.length_scales_eval["positive shell"]

#%
time_in_sec = solution.t * timescale

eta_plot = np.linspace(0.0, 1, 30) * lengthscale_core
chi_plot = np.linspace(0.0, 1, 20) * lengthscale_shell


c_c_out_value = c_c_out(t=time_in_sec, x=eta_plot)
c_s_out_value = c_s_out(t=time_in_sec, x=chi_plot)

#%%
plot_times = np.linspace(0, max(time_in_sec), 11)

plt.figure(figsize=(12, 8))
cmap = plt.get_cmap("inferno")
for i, tp in enumerate(plot_times):
    color = cmap(float(i) / len(plot_times))
    plt.plot(
        eta_plot * 0.5,
        c_c_out(tp, x=eta_plot),
        "-o",
        color=color,
        label="core" if i == 0 else ""
    )
    plt.plot(
        ( chi_plot / lengthscale_shell * (1 - 0.5) + 0.5) * lengthscale_shell,
        c_s_out(tp, x=chi_plot),
        "x",
        color=color,
        label="shell" if i == 0 else ""
    )
    
# plt.ylim(0, 1.0)

plt.xlabel("r", fontsize=16)
plt.ylabel("c", fontsize=16)
plt.legend()
plt.show()

