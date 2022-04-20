
# c_c_out = solution["Positive core concentration [mol.m-3]"]
# c_o_out = solution["Positive shell concentration of oxygen [mol.m-3]"]

# s_out = solution["X-averaged moving phase boundary location"]


# # c_tot = solution["Total lithium [mol]"]
# lli = solution["Loss of lithium inventory [%]"]
# lam = solution["X-averaged loss of active material in positive electrode (MZ)"]
# lam_full = solution["Loss of active material in positive electrode (MZ)"]

# V_cell = solution["Terminal voltage [V]"]
# I = solution["Current [A]"]

#%%
timescale = solution.timescale_eval
lengthscale_core = solution.length_scales_eval["positive core"]
lengthscale_shell = solution.length_scales_eval["positive shell"]

#%
time_in_sec = solution.t * timescale
time_in_sec = solution["Time [s]"].entries
time_output = np.linspace(0, max(time_in_sec), 1001) 



# s_out_value = s_out(t=time_in_sec)
# 

# lli_value = lli(t=time_in_sec)
# lam_value = lam(t=time_in_sec)
# lam_full_value = lam_full(t=time_in_sec, x=xed)

# V_cell_value = V_cell(t=time_in_sec)
# I_value = I(t=time_in_sec)

#%%

chi = np.linspace(0.0, 1, 20) * lengthscale_shell
xed = 86.7E-6 + 12E-6 + np.linspace(0.0, 1, 20) * 66.2E-6

time_plot = 55*60

c_o = solution["Positive shell concentration of oxygen [mol.m-3]"]

c_o_55 = []
x = []
y = []
for i in xed:
    for j in chi:
        x.append(i)
        y.append(j)
        c_o_55.append(c_o(t=time_plot, r=j, x=i))

x = np.asarray(x)
y = np.asarray(y)
c_o_55 = np.asarray(c_o_55)
#%%
        

f = open("pybamm\mz_develop\output\DFN_c_o_55.csv", "w")
np.savetxt(f, np.c_[x, y, c_o_55], fmt='%1.6e', delimiter=', ')
f.close()



#%%
I = solution["Current [A]"].entries

f = open("pybamm\mz_develop\output\DFN_current.txt", "w")
np.savetxt(f, np.c_[time_in_sec, I], fmt='%1.6e', delimiter=', ')
f.close()



#%%
eta_plot = np.linspace(0.0, 1, 30) * lengthscale_core

# c_c_out_value = c_c_out(t=time_in_sec, r=eta_plot)

#%%

xed = 86.7E-6 + 12E-6 + np.linspace(0.0, 1, 20) * 66.2E-6
s_out = solution["Moving phase boundary location"]
s_out_value = s_out_full(t=95*60, x=xed)

f = open("pybamm\mz_develop\output\DFN_phase_location_95mins.txt", "w")
np.savetxt(f, np.c_[xed, s_out_value], fmt='%1.6e', delimiter=', ')
f.close()


#%%

xed = 86.7E-6 + 12E-6 + np.linspace(0.0, 1, 20) * 66.2E-6
lam_out = solution["Loss of active material in positive electrode (MZ)"]
lam_out_value = lam_out(t=95*60, x=xed)

f = open("pybamm\mz_develop\output\DFN_lam_95mins.txt", "w")
np.savetxt(f, np.c_[xed, lam_out_value], fmt='%1.6e', delimiter=', ')
f.close()

#%%

xed = 86.7E-6 + 12E-6 + np.linspace(0.0, 1, 20) * 66.2E-6
c_surf = solution["Positive core surface concentration [mol.m-3]"]
c_surf_value = c_surf(t=0*60, x=xed)

f = open("pybamm\mz_develop\output\DFN_c_surf_0mins.txt", "w")
np.savetxt(f, np.c_[xed, c_surf_value], fmt='%1.6e', delimiter=', ')
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
lli_cyc = solution["LLI_cyc"](t=time_output)
f = open("pybamm\mz_develop\output\DFN_lli_cyc.txt", "w")
np.savetxt(f, np.c_[time_output, lli_cyc], fmt='%1.6e', delimiter=', ')
f.close()


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
f = open("pybamm\mz_develop\output\DFN_lam_110mins.txt", "w")
np.savetxt(f, np.c_[xed, lam_full_value[:,110]], fmt='%1.6e', delimiter=', ')
f.close()