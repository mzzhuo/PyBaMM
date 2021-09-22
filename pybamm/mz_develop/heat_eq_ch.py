#
# Solve the transient heat equation with a spatially-dependent source term
#

import pybamm
import numpy as np
import matplotlib.pyplot as plt

#%%
# Numerical solution ----------------------------------------------------------

# Start with a base model
model = pybamm.BaseModel()

#%
# Define the variables and parameters
# Note: we need to define the spatial variable x here too, so we can use it
# to write down the source term
x_le = pybamm.SpatialVariable("x_le", domain="rod_left", coord_sys="cartesian")
T_le = pybamm.Variable("Temperature_le", domain="rod_left")
k_le = pybamm.Parameter("Thermal diffusivity in the left")

x_ri = pybamm.SpatialVariable("x_ri", domain="rod_right", coord_sys="cartesian")
T_ri = pybamm.Variable("Temperature_ri", domain="rod_right")
k_ri = pybamm.Parameter("Thermal diffusivity in the right")

# Write the governing equations
N_le = -k_le * pybamm.grad(T_le)  # Heat flux
Q_le = 1 - pybamm.Function(np.abs, x_le - 1)  # Source term
dTdt_le = -pybamm.div(N_le) + Q_le
model.rhs[T_le] = dTdt_le  # add to model

N_ri = -k_ri * pybamm.grad(T_ri)  # Heat flux
Q_ri = 1 - pybamm.Function(np.abs, x_ri - 1)  # Source term
dTdt_ri = -pybamm.div(N_ri) + Q_ri
model.rhs[T_ri] = dTdt_ri  # add to model

# T_infc = pybamm.Variable("Temperature at interface")

T = pybamm.Concatenation(T_le, T_ri)

#%

# Add the boundary and initial conditions
# model.boundary_conditions[T_le] = {"left": (0, "Dirichlet"), "right": (pybamm.Scalar(1.0), "Dirichlet")}
# model.boundary_conditions[T_ri] = {"left": (pybamm.Scalar(1.0), "Dirichlet"), "right": (pybamm.Scalar(0), "Dirichlet")}

#%
# model.boundary_conditions = {T_le: {"left": (pybamm.Scalar(0), "Dirichlet"), "right": (pybamm.Scalar(1.0), "Dirichlet")},
#                               T_ri: {"left": (pybamm.Scalar(1.0), "Dirichlet"), "right": (pybamm.Scalar(0), "Dirichlet")}}

# model.boundary_conditions = {T_le: {"left": (pybamm.Scalar(0), "Dirichlet"), "right": (T_ri_infc, "Dirichlet")},
#                               T_ri: {"left": (T_le_infc, "Dirichlet"), "right": (pybamm.Scalar(0), "Dirichlet")}}
                             
# model.boundary_conditions = {T_le: {"left": (pybamm.Scalar(0), "Dirichlet"), "right": (T_le_infc, "Dirichlet")},
#                              T_ri: {"left": (T_le_infc, "Dirichlet"), "right": (pybamm.Scalar(0), "Dirichlet")}}



# -------------------------------------------------------------------------------------------------------------------
# model.boundary_conditions = {T: {"left": (pybamm.Scalar(0), "Dirichlet"), "right": (pybamm.Scalar(0), "Dirichlet")}}

# or
T_le_N = pybamm.boundary_cell_value(T_le, "right")
T_ri_1 = pybamm.boundary_cell_value(T_ri, "left")

# boundary cell length: node to the edge
dx_T_le = pybamm.boundary_cell_length(T_le, "right")
dx_T_ri = pybamm.boundary_cell_length(T_ri, "left")

# T_le_infc = pybamm.boundary_value(T_le, "right")
# T_ri_infc = pybamm.boundary_value(T_ri, "left")


T_shared = ((k_le / dx_T_le * T_le_N + k_ri / dx_T_ri * T_ri_1 ) / 
            (k_le / dx_T_le + k_ri / dx_T_ri ))
# T_shared = (T_le_N + T_ri_1) / 2


T_le_infc  = (T_shared - T_le_N) / dx_T_le
T_ri_infc = (T_ri_1 - T_shared) / dx_T_ri




model.boundary_conditions = {T_le: {"left": (pybamm.Scalar(0), "Dirichlet"), "right": (T_le_infc, "Neumann")},
                             T_ri: {"left": (T_ri_infc, "Neumann"), "right": (pybamm.Scalar(0), "Dirichlet")}}


# -------------------------------------------------------------------------------------------------------------------


#%
# {2 * x_le - x_le ** 2} is a set object so no curly
model.initial_conditions[T_le] = 2 * x_le - x_le ** 2
model.initial_conditions[T_ri] = 2 * x_ri - x_ri ** 2
#%
# zzz = model.initial_conditions

# model.initial_conditions = {T_le: 2 * x_le - x_le ** 2,
#                             T_ri: 2 * x_ri - x_ri ** 2}


#%


# zzz = model.boundary_conditions
# zzz = model.initial_conditions



#%
# Add desired output variables
model.variables = {
    "Temperature in the left": T_le, 
    "Heat flux in the left": N_le, 
    # "Heat source in the left": Q_le,
    "Temperature in the right": T_ri, 
    "Heat flux in the right": N_ri, 
    # "Heat source in the right": Q_ri
}

#%
# Define geometry
geometry = {"rod_left":  {x_le: {"min": pybamm.Scalar(0), "max": pybamm.Scalar(1)}}, 
            "rod_right": {x_ri: {"min": pybamm.Scalar(1), "max": pybamm.Scalar(2)}}}

# Set parameter values
param = pybamm.ParameterValues({"Thermal diffusivity in the left": 0.75, 
                                "Thermal diffusivity in the right": 1.5})

#%
# Process model and geometry
param.process_model(model)
#
param.process_geometry(geometry)


#%
# Pick mesh, spatial method, and discretise
submesh_types = {"rod_left": pybamm.Uniform1DSubMesh,
                 "rod_right": pybamm.Uniform1DSubMesh}
var_pts = {x_le: 10, x_ri: 5}

mesh = pybamm.Mesh(geometry, submesh_types, var_pts)
spatial_methods = {"rod_left": pybamm.FiniteVolume(),
                   "rod_right": pybamm.FiniteVolume()}

#%
disc = pybamm.Discretisation(mesh, spatial_methods)

#%
model_disc = disc.process_model(model,inplace=False)



#%%
# Solve
solver = pybamm.ScipySolver()
t = np.linspace(0, 1, 100)
solution = solver.solve(model_disc, t)

#%
# Extract output variables
T_out_le = solution["Temperature in the left"]
T_out_ri = solution["Temperature in the right"]

N_out_le = solution["Heat flux in the left"]
N_out_ri = solution["Heat flux in the right"]

# #%%
# # Exact solution -------------------------------------------------------
# N = 100  # number of Fourier modes to sum
# # k_val = param["Thermal diffusivity"]  # extract value of diffusivity
# k_val = 0.75

# # Fourier coefficients
# def q(n):
#     return (8 / (n ** 2 * np.pi ** 2)) * np.sin(n * np.pi / 2)


# def c(n):
#     return (16 / (n ** 3 * np.pi ** 3)) * (1 - np.cos(n * np.pi))


# def b(n):
#     return c(n) - 4 * q(n) / (k_val * n ** 2 * np.pi ** 2)


# def T_n(t, n):
#     return (4 * q(n) / (k_val * n ** 2 * np.pi ** 2)) + b(n) * np.exp(
#         -k_val * (n * np.pi / 2) ** 2 * t
#     )


# # Sum series to get the source term
# def Q_exact(n):
#     out = 0
#     for n in np.arange(1, N):
#         out += q(n) * np.sin(n * np.pi * x / 2)
#     return out


# # Sum series to get the temperature
# def T_exact(x, t):
#     out = 0
#     for n in np.arange(1, N):
#         out += T_n(t, n) * np.sin(n * np.pi * x / 2)
#     return out


#%%
# Plot ------------------------------------------------------------------------
# x_le_nodes = mesh["rod_left"].nodes  # numerical gridpoints
# x_ri_nodes = mesh["rod_right"].nodes  # numerical gridpoints

x_le_nodes = np.linspace(0, 1, 20)
x_ri_nodes = np.linspace(1, 2, 20)

xx = np.linspace(0, 2, 101)  # fine mesh to plot exact solution
plot_times = np.linspace(0, 1, 2)

plt.figure(figsize=(15, 8))
cmap = plt.get_cmap("inferno")
for i, t in enumerate(plot_times):
    color = cmap(float(i) / len(plot_times))
    plt.plot(
        x_le_nodes,
        T_out_le(t, x=x_le_nodes),
        "o",
        color=color,
        label="left temp" if i == 0 else ""
    )
    plt.plot(
        x_ri_nodes,
        T_out_ri(t, x=x_ri_nodes),
        "x",
        color=color,
        label="right temp" if i == 0 else ""
    )
    # plt.plot(
    #     xx, T_exact(xx, t), "-", color=color, label="Exact (t={})".format(plot_times[i])
    # )
plt.xlabel("x", fontsize=16)
plt.ylabel("T", fontsize=16)
plt.legend()
plt.show()


#%%
# Plot ------------------------------------------------------------------------
x_le = np.linspace(0, 1, 20)
x_ri = np.linspace(1, 2, 20)

plot_times = np.linspace(0, 1, 2)

plt.figure(figsize=(15, 8))
cmap = plt.get_cmap("inferno")
for i, t in enumerate(plot_times):
    color = cmap(float(i) / len(plot_times))
    plt.plot(
        x_le,
        N_out_le(t, x=x_le),
        "o",
        color=color,
        label="left flux" if i == 0 else ""
    )
    plt.plot(
        x_ri,
        N_out_ri(t, x=x_ri),
        "x",
        color=color,
        label="right flux" if i == 0 else ""
    )
plt.xlabel("x", fontsize=16)
plt.ylabel("N", fontsize=16)
plt.legend()
plt.show()












