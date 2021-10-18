import pybamm
import numpy as np
# import os
import matplotlib.pyplot as plt



#%%
# ---------------------------------------------------------
#
model = pybamm.lithium_ion.SPM({"PE phase transition": "yes"})

geometry = model.default_geometry

param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Zhuo2021)
# param = model.default_parameter_values
param["Current function [A]"] = -1.675
#
param.process_model(model)
#
param.process_geometry(geometry)

#%

submesh_types = model.default_submesh_types
var_pts = model.default_var_pts

mesh = pybamm.Mesh(geometry, submesh_types, var_pts)

spatial_methods = model.default_spatial_methods

disc = pybamm.Discretisation(mesh, spatial_methods)

model_disc = disc.process_model(model, inplace=False);

#%%
solver = pybamm.ScipySolver()
# solver = pybamm.CasadiSolver(mode="safe")

# timescale = param.evaluate(model.timescale)
t_eval = np.linspace(0, 7300, 100)

solution = solver.solve(model_disc, t_eval)

#%%
output_variables = [
    "X-averaged negative particle concentration [mol.m-3]",
    # "X-averaged positive particle concentration [mol.m-3]",
    #---------------------------------------------------------
    "X-averaged positive core concentration [mol.m-3]",
    #---------------------------------------------------------
    # "X-averaged positive shell concentration [mol.m-3]",
    # #---------------------------------------------------------
    "X-averaged positive shell concentration of oxygen [mol.m-3]",
    #---------------------------------------------------------
    # "X-averaged moving phase boundary location",
    # #---------------------------------------------------------
    # "X-averaged shared concentration at core-shell interface [mol.m-3]",
    #---------------------------------------------------------
    # "Current [A]",
    "Terminal voltage [V]",
    # "X-averaged positive electrode temperature [K]",
    # "X-averaged time derivative of moving phase boundary location",
]
#%
output_variables.extend(
    [
        #---------------------------------------------------------
        [
            "Total lithium in positive electrode [mol]",
            "Total lithium in negative electrode [mol]",
            "Total lithium [mol]",
        ],
        "Loss of lithium inventory [%]",
        # "Electrolyte concentration [mol.m-3]",
        # "Negative electrode potential [V]",
        # "Electrolyte potential [V]",
        # "Positive electrode potential [V]",
    ]
)

#%%
# solution = sim.solution
plot = pybamm.QuickPlot(
    solution,
    output_variables,
    # time_unit="seconds",
    # spatial_unit="um",
)
plot.dynamic_plot()














#%%
model = pybamm.lithium_ion.SPM({"PE phase transition": "yes"})
param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Zhuo2021)
#%
experiment = pybamm.Experiment(
    [
        (
            "Charge at 1.0 C until 4.2 V",
            # "Hold at 4.2 V until C/50",
            # # "Rest for 30 minutes",
            # "Discharge at 1.0 C until 2.9 V",
            # "Hold at 2.9 V until C/10",
            # "Rest for 30 minutes",
        )
    ] * 1
)

sim = pybamm.Simulation(
    model, experiment=experiment,
    parameter_values=param,
)
solution = sim.solve()

#%%
# param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Zhuo2021)
# param["Current function [A]"] = -1.675

# sim = pybamm.Simulation(model, parameter_values=param)
# solution = sim.solve([0, 7300])

#%%
sim.plot(output_variables)








