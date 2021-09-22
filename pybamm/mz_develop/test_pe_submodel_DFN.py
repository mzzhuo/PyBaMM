import pybamm
import numpy as np
# import os
import matplotlib.pyplot as plt


# %%
model = pybamm.lithium_ion.DFN({"PE phase transition": "yes"})

#%
# param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Zhuo2021)
# param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Chen2020)

#%%
# 1 C = 3.35
# param["Current function [A]"] = -3.35 * 0.5
# sim = pybamm.Simulation(model, parameter_values=param)
# solution = sim.solve([0, 7500])

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
    # parameter_values=param,
    # solver=pybamm.CasadiSolver("fast with events"),
)
solution = sim.solve(calc_esoh=False)

# #%%
# output_variables = [
#     # "Negative particle surface concentration [mol.m-3]",
#     "Negative electrode volume-averaged concentration [mol.m-3]",
#     # "Electrolyte concentration [mol.m-3]",
#     # "Positive particle surface concentration [mol.m-3]",
#     "Positive electrode volume-averaged concentration [mol.m-3]",
#     "Current [A]",
#     "Negative electrode surface potential difference [V]",
#     # "Electrolyte potential [V]",
#     "Positive electrode surface potential difference [V]",
#     "Terminal voltage [V]",
#     "X-averaged negative electrode interfacial current density [A.m-2]",
#     "X-averaged positive electrode interfacial current density [A.m-2]",   
#     # "Exchange current density",
#     "Total current density [A.m-2]",
# ]

# sim.plot(output_variables)


#%%    
output_variables = [
    # "Moving phase boundary location",
    "X-averaged moving phase boundary location",
    #---------------------------------------------------------
    # "Shared concentration at core-shell interface [mol.m-3]",     
    "X-averaged shared concentration at core-shell interface [mol.m-3]",
    #---------------------------------------------------------
    # "Electrolyte concentration [mol.m-3]",
    "Current [A]",
    # "Negative electrode interfacial current density [A.m-2]",
    # "Positive electrode interfacial current density [A.m-2]",  
    # "X-averaged negative electrode interfacial current density [A.m-2]",
    # "X-averaged positive electrode interfacial current density [A.m-2]",  
    "Terminal voltage [V]",
    # "X-averaged positive electrode temperature [K]",
    # "Time derivative of moving phase boundary location",
    # "X-averaged time derivative of moving phase boundary location",
    [
    "Total lithium in positive electrode [mol]",
    "Total lithium in negative electrode [mol]",
    "Total lithium [mol]",
    ],
    # "LAM_ne [%]",
    # "LAM_pe [%]",
    "Loss of lithium inventory [%]"
]

sim.plot(output_variables)



#%%
which_average = "X"

if which_average == "X":
    output_variables = [
        "X-averaged negative particle concentration [mol.m-3]",
        "X-averaged positive core concentration [mol.m-3]",
        "X-averaged positive shell concentration [mol.m-3]",
        "X-averaged positive shell concentration of oxygen [mol.m-3]",
    ]
elif which_average == "R":
    output_variables = [
        "R-averaged negative particle concentration",
        "R-averaged positive core concentration [mol.m-3]",
        "R-averaged positive shell concentration [mol.m-3]",
        "R-averaged positive shell concentration of oxygen [mol.m-3]",
    ]
else:
    output_variables = [
        "Negative particle concentration [mol.m-3]",
        "Positive core concentration [mol.m-3]",
        "Positive shell concentration [mol.m-3]",
        "Positive shell concentration of oxygen [mol.m-3]"
    ] 
    
    
output_variables.extend(
    [
        # "Moving phase boundary location",
        "X-averaged moving phase boundary location",
        #---------------------------------------------------------
        # "Shared concentration at core-shell interface [mol.m-3]",     
        "X-averaged shared concentration at core-shell interface [mol.m-3]",
        "Loss of lithium inventory [%]",
        #---------------------------------------------------------
        # "Electrolyte concentration [mol.m-3]",
        "Current [A]",
        "Terminal voltage [V]",
        # "X-averaged positive electrode temperature [K]",
        # "Time derivative of moving phase boundary location",
        # "X-averaged time derivative of moving phase boundary location",
        [
        "Total lithium in positive electrode [mol]",
        "Total lithium in negative electrode [mol]",
        "Total lithium [mol]",
        ],
        # "Electrolyte concentration [mol.m-3]",
        # "Negative electrode potential [V]",
        # "Electrolyte potential [V]",
        # "Positive electrode potential [V]",
    ]
)


sim.plot(output_variables)

#%%

output_variables = [
        "Current [A]",
        "Terminal voltage [V]",
        # "X-averaged moving phase boundary location",
        # "X-averaged shared concentration at core-shell interface [mol.m-3]",
        # "Moving phase boundary location",
        # "Shared concentration at core-shell interface [mol.m-3]",     
        "Loss of lithium inventory [%]",
        [
        "Total lithium in positive electrode [mol]",
        "Total lithium in negative electrode [mol]",
        "Total lithium [mol]",
        ],
        "X-averaged moving phase boundary location",
        #---------------------------------------------------------
        # "Shared concentration at core-shell interface [mol.m-3]",     
        # "X-averaged shared concentration at core-shell interface [mol.m-3]",
        "X-averaged positive shell surface concentration [mol.m-3]",
        # "Negative particle concentration [mol.m-3]",
        # "Positive shell concentration of oxygen [mol.m-3]",
]



sim.plot(output_variables,n_rows=2)






#%%
model = pybamm.lithium_ion.DFN({"PE phase transition": "yes"})
geometry = model.default_geometry

# param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Chen2020)
param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Zhuo2021)
# param = model.default_parameter_values
param["Current function [A]"] = -1.675


#%
param.process_model(model)
#%
param.process_geometry(geometry)

#%%

submesh_types = model.default_submesh_types
var_pts = model.default_var_pts

mesh = pybamm.Mesh(geometry, submesh_types, var_pts)

spatial_methods = model.default_spatial_methods

disc = pybamm.Discretisation(mesh, spatial_methods)
#%%
model_disc = disc.process_model(model, inplace=False);

#%%
# solver = pybamm.ScipySolver()
solver = pybamm.CasadiSolver()

timescale = param.evaluate(model.timescale)
t_eval = np.linspace(0, 7300, 100)

solution = solver.solve(model_disc, t_eval)

#%%
which_average = "X"

if which_average == "X":
    output_variables = [
        "X-averaged negative particle concentration [mol.m-3]",
        "X-averaged positive core concentration [mol.m-3]",
        "X-averaged positive shell concentration [mol.m-3]",
        "X-averaged positive shell concentration of oxygen [mol.m-3]",
    ]
elif which_average == "R":
    output_variables = [
        "R-averaged negative particle concentration",
        "R-averaged positive core concentration [mol.m-3]",
        "R-averaged positive shell concentration [mol.m-3]",
        "R-averaged positive shell concentration of oxygen [mol.m-3]",
    ]
else:
    output_variables = [
        "Negative particle concentration [mol.m-3]",
        "Positive core concentration [mol.m-3]",
        "Positive shell concentration [mol.m-3]",
        "Positive shell concentration of oxygen [mol.m-3]"
    ] 
    
    
output_variables.extend(
    [
        # "Moving phase boundary location",
        "X-averaged moving phase boundary location",
        #---------------------------------------------------------
        # "Shared concentration at core-shell interface [mol.m-3]",     
        "X-averaged shared concentration at core-shell interface [mol.m-3]",
        #---------------------------------------------------------
        # "Current [A]",
        "Terminal voltage [V]",
        # "X-averaged positive electrode temperature [K]",
        # "Time derivative of moving phase boundary location",
        # "X-averaged time derivative of moving phase boundary location",
        [
        "Total lithium in positive electrode [mol]",
        "Total lithium in negative electrode [mol]",
        "Total lithium [mol]",
        ],
        "LLI [%]",
        # "Electrolyte concentration [mol.m-3]",
        # "Negative electrode potential [V]",
        # "Electrolyte potential [V]",
        # "Positive electrode potential [V]",
    ]
)
#%%
# output_variables = [
#     # "Negative particle surface concentration [mol.m-3]",
#     "Negative electrode volume-averaged concentration [mol.m-3]",
#     # "Electrolyte concentration [mol.m-3]",
#     # "Positive particle surface concentration [mol.m-3]",
#     "Positive electrode volume-averaged concentration [mol.m-3]",
#     "Current [A]",
#     "Negative electrode surface potential difference [V]",
#     # "Electrolyte potential [V]",
#     "Positive electrode surface potential difference [V]",
#     "Terminal voltage [V]",
#     "X-averaged negative electrode interfacial current density [A.m-2]",
#     "X-averaged positive electrode interfacial current density [A.m-2]",   
#     # "Exchange current density",
#     "Total current density [A.m-2]"
# ]

plot = pybamm.QuickPlot(
    solution,
    output_variables,
    # time_unit="seconds",
    # spatial_unit="um",
)

plot.dynamic_plot()



