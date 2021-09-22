
import pybamm
import numpy as np
import matplotlib.pyplot as plt

#%% simulate the experiments
timer = pybamm.Timer()

experiment = pybamm.Experiment(
    [
        (
            "Charge at 1.0 C until 4.2 V",
            "Hold at 4.2 V until C/50",
            "Rest for 30 minutes",
            "Discharge at 1.0 C until 2.8 V",
            # "Hold at 2.8 V until C/50",
            "Rest for 30 minutes"
        )
    ] * 1
)

#%
model = pybamm.lithium_ion.PeDegradationModel(name='DegradationModel')

#%%
sim = pybamm.Simulation(model, experiment=experiment)


#
solution = sim.solve()

print(timer.time())

#%%
sim.plot(
    [
        "Discharge capacity [A.h]",
        "Current [A]",
        "Terminal voltage [V]",
        # "Temperature [Kelvin]",
        # 
        "Normalized lithium concentration in positive particle core [-]",
        "Normalized moving phase boundary location [-]", 
        "Normalized lithium concentration in positive particle shell [-]",
        # 
        # "Normalized oxygen concentration in positive particle shell [-]",
        # "Normalized lithium concentration in negative particle [-]",
        # "J_p_scaled",
        "Normalized interface concentration [-]",
        "Normalized lithium concentration in positive shell surface [-]"
    ],
    time_unit = "minutes",
)








#%% run simulation user-specified

model = PeDegradationModel(name='DegradationModel')

sim = pybamm.Simulation(model)

t_eval = np.linspace(0, model.t_end, 101)

solution = sim.solve(t_eval = t_eval)




#%%

# eta = pybamm.SpatialVariable(
#     "eta", domain=["positive particle core"], coord_sys="cartesian")
# chi = pybamm.SpatialVariable(
#     "chi", domain=["positive particle shell"], coord_sys="cartesian")
# r_n = pybamm.SpatialVariable(
#     "r_n", domain=["negative particle"], coord_sys="spherical polar")
     
# #%
# geometry = {"negative particle": {r_n: {"min": pybamm.Scalar(0.0), "max": pybamm.Scalar(1)}}, 
#             "positive particle core": {eta: {"min": pybamm.Scalar(0.0), "max": pybamm.Scalar(1)}},
#             "positive particle shell": {chi: {"min": pybamm.Scalar(0.0), "max": pybamm.Scalar(1)}}
#             }

# #%
# param = pybamm.ParameterValues("pybamm/myDevelop/data/parameters_abir.csv")
# # change string type to symbol
# # [param.update({k: eval(v)}) for k, v in param.items() if isinstance(v, str)]

# #%
# param.process_model(model)
# param.process_geometry(geometry)

# #%
# submesh_types = {"positive particle core": pybamm.MeshGenerator(pybamm.Uniform1DSubMesh),
#                   "positive particle shell": pybamm.MeshGenerator(pybamm.Uniform1DSubMesh),
#                   "negative particle": pybamm.MeshGenerator(pybamm.Uniform1DSubMesh)}

# var_pts = {eta: 20, chi: 5, r_n: 20}
# mesh = pybamm.Mesh(geometry, submesh_types, var_pts)

# #%
# spatial_methods = {"positive particle core": pybamm.FiniteVolume(),
#                     "positive particle shell": pybamm.FiniteVolume(),
#                     "negative particle": pybamm.FiniteVolume()}

# #%

# disc = pybamm.Discretisation(mesh, spatial_methods)

# #%

# model_disc = disc.process_model(model, inplace=False);


# #%%
# solver = pybamm.ScipySolver()
# # solver = pybamm.CasadiSolver(mode="safe")

# # simulation related parameters
        
# t_end_scaled_value = param.process_symbol(model.t_end_scaled).value

# t_eval = np.linspace(0, t_end_scaled_value, 101)
# # t_eval = t_eval[0:99]
# #%%
# solution = solver.solve(model_disc, t_eval)








