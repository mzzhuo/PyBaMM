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
            "Discharge at 0.5 C until 2.8 V",
            "Hold at 2.8 V until C/50",
            "Rest for 60 minutes",
        )
    ] * 20,
    # need for case I
    period="0.5 minute",  
)

sim = pybamm.Simulation(
    model, experiment=experiment,
    parameter_values=param,
)


#%%
solution = sim.solve(calc_esoh=False)


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

#%%
f = open("pybamm\mz_develop\output\capacity.txt", "a")
np.savetxt(f, np.c_[cycle_numbers, Q_dis_cycles], fmt='%12.6f', delimiter=', ')
f.close()
#%%
plt.figure(figsize=(8, 6))
# markerfacecolor
plt.plot(cycle_numbers, Q_dis_cycles, 'o', mfc='none', label="LAM and LLI")
plt.legend()

#%%
# np.save('Q_dis_cycles_2', Q_dis_cycles_2)
Q_dis_cycles_0 = np.load('Q_dis_cycles_0.npy')
Q_dis_cycles_1 = np.load('Q_dis_cycles_1.npy')

#%%
plt.figure(figsize=(8, 6))
# markerfacecolor
plt.plot(cycle_numbers, Q_dis_cycles_0, '-o', mfc='none', label="only LAM")
plt.plot(cycle_numbers, Q_dis_cycles_1, '-s', mfc='none', label="LAM and LLI")
plt.plot(cycle_numbers, Q_dis_cycles_2, '-*', mfc='none', label="LAM, LLI, and shell")
plt.legend()

#%%
# experiment2 = pybamm.Experiment(
#     [
#         (
#             "Charge at 1.0 C until 4.2 V",
#             # "Hold at 4.2 V until C/50",
#             # "Rest for 60 minutes",
#             # "Discharge at 0.5 C until 3.0 V",
#             # "Hold at 3.0 V until C/50",
#             # "Rest for 60 minutes",
#         )
#     ] * 1
# )

# #%
# new_model = model.set_initial_conditions_from(solution, inplace=False)
# #%
# new_sim = pybamm.Simulation(
#     new_model, 
#     experiment=experiment2,
#     parameter_values=param,
# )

# solution2 = new_sim.solve(calc_esoh=False)

# %%

# >> 49340*0.72*0.745*4.7227e-6

# ans =

#     0.1250

# >> 34257*0.85*0.694*6.1852e-6

# ans =

#     0.1250

# >> 49340*0.72*0.745*4.7227e-6*96487/3600

# ans =

#     3.3500











#%%
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
