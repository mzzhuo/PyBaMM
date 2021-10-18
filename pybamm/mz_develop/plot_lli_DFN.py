
c_c_out = solution["Positive core concentration [mol.m-3]"]
c_o_out = solution["Positive shell concentration of oxygen [mol.m-3]"]

s_out = solution["X-averaged moving phase boundary location"]
s_out_full = solution["Moving phase boundary location"]

# c_tot = solution["Total lithium [mol]"]
lli = solution["Loss of lithium inventory [%]"]
lam = solution["X-averaged loss of active material in positive electrode (MZ)"]
lam_full = solution["Loss of active material in positive electrode (MZ)"]

V_cell = solution["Terminal voltage [V]"]
I = solution["Current [A]"]

#%%
timescale = solution.timescale_eval
lengthscale_core = solution.length_scales_eval["positive core"]
lengthscale_shell = solution.length_scales_eval["positive shell"]

#%
time_in_sec = solution.t * timescale

xed = 86.7E-6 + 12E-6 + np.linspace(0.0, 1, 20) * 66.2E-6

s_out_value = s_out(t=time_in_sec)
s_out_full_value = s_out_full(t=time_in_sec, x=xed)

lli_value = lli(t=time_in_sec)
lam_value = lam(t=time_in_sec)
lam_full_value = lam_full(t=time_in_sec, x=xed)

V_cell_value = V_cell(t=time_in_sec)
I_value = I(t=time_in_sec)
#%%
eta_plot = np.linspace(0.0, 1, 30) * lengthscale_core

# c_c_out_value = c_c_out(t=time_in_sec, r=eta_plot)

#%%
f = open("pybamm\mz_develop\output\DFN_phase_locationXave.txt", "w")
np.savetxt(f, np.c_[time_in_sec, s_out_value], fmt='%1.6e', delimiter=', ')
f.close()
#%%
f = open("pybamm\mz_develop\output\DFN_phase_location_60mins.txt", "w")
np.savetxt(f, np.c_[xed, s_out_full_value[:,60]], fmt='%1.6e', delimiter=', ')
f.close()
#%%
f = open("pybamm\mz_develop\output\DFN_phase_location_70mins.txt", "w")
np.savetxt(f, np.c_[xed, s_out_full_value[:,70]], fmt='%1.6e', delimiter=', ')
f.close()

#%%
#%%
f = open("pybamm\mz_develop\output\DFN_V_cell.txt", "w")
np.savetxt(f, np.c_[time_in_sec, V_cell_value], fmt='%1.6e', delimiter=', ')
f.close()

f = open("pybamm\mz_develop\output\DFN_I.txt", "w")
np.savetxt(f, np.c_[time_in_sec, I_value], fmt='%1.6e', delimiter=', ')
f.close()

#%%
f = open("pybamm\mz_develop\output\DFN_lli.txt", "w")
np.savetxt(f, np.c_[time_in_sec, lli_value], fmt='%1.6e', delimiter=', ')
f.close()

#%%
f = open("pybamm\mz_develop\output\DFN_lam.txt", "w")
np.savetxt(f, np.c_[time_in_sec, lam_value], fmt='%1.6e', delimiter=', ')
f.close()

#%%
f = open("pybamm\mz_develop\output\DFN_lam_110mins.txt", "w")
np.savetxt(f, np.c_[xed, lam_full_value[:,110]], fmt='%1.6e', delimiter=', ')
f.close()