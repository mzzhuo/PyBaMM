import pybamm
import numpy as np
# import os
import matplotlib.pyplot as plt
import pybamm.mz_develop.output_module as outmod

#%%
# ---------------------------------------------------------
#
model = pybamm.lithium_ion.SPM({"PE degradation": "yes"})
#%%
geometry = model.default_geometry

param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Zhuo2021)
# param = model.default_parameter_values
param["Current function [A]"] = -1.675
# param["Current function [A]"] = 0
#
param.process_model(model)
#
param.process_geometry(geometry)

#%%

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
# t_eval = np.linspace(0, 1800, 100)

solution = solver.solve(model_disc, t_eval)

#%%
# solution = sim.solution
output_variables = outmod.output_variables_spm
plot = pybamm.QuickPlot(
    solution,
    output_variables,
    # time_unit="seconds",
    # spatial_unit="um",
)
plot.dynamic_plot()








#%%
import pybamm
import numpy as np
# import os
import matplotlib.pyplot as plt
import pybamm.mz_develop.output_module as outmod

#%%
model = pybamm.lithium_ion.SPM({"PE degradation": "yes"})
#%
param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Zhuo2021)
#%%
experiment = pybamm.Experiment(
    [
        (
            "Charge at 0.5 C until 4.2 V",
            "Hold at 4.2 V until C/50",
            "Rest for 60 minutes",
            "Discharge at 0.1 C until 2.8 V",
            "Hold at 2.8 V until C/50",
            "Rest for 60 minutes",
        )
    ] * 20
)

sim = pybamm.Simulation(
    model, experiment=experiment,
    parameter_values=param,
)
#%%
solution = sim.solve()


#%%
output_variables = outmod.output_variables_spm
sim.plot(output_variables)  

#%%
total_cycles = len(solution.cycles)
Q_dis_cycles = []
cycle_numbers = np.arange(1, total_cycles+1)
for i in range(total_cycles):
    Q_dis_cyc = solution.cycles[i].steps[3]["Discharge capacity [A.h]"].entries
    Q_dis = Q_dis_cyc[-1] - Q_dis_cyc[0]
    Q_dis_cycles.append(Q_dis)
    
plt.figure(figsize=(8, 6))
# markerfacecolor
plt.plot(cycle_numbers, Q_dis_cycles, 'o', mfc='none')
 






