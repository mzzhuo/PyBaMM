#%%

# c_c_out = solution.cycles[0].steps[2]["X-averaged positive core concentration [mol.m-3]"]
# s_out = solution["X-averaged moving phase boundary location"]
# c_shared_out = solution["X-averaged lithium concentration at core-shell interface [mol.m-3]"]
# c_tot = solution["Total lithium [mol]"]
# lli = solution["LLI"]
# lli_cyc = solution["LLI_cyc"]


# V_cell = solution["Terminal voltage [V]"]

# Q = solution.cycles[0]["Discharge capacity [A.h]"]
#%%
timescale = solution.timescale_eval
lengthscale_core = solution.length_scales_eval["positive core"]
lengthscale_shell = solution.length_scales_eval["positive shell"]

#%
# time_in_sec = solution.t * timescale
time_in_sec = solution["Time [s]"].entries
time_output = np.linspace(0, max(time_in_sec), 1001) 

# s_out_value = s_out(t=time_in_sec)
# c_shared_out_value = c_shared_out(t=time_in_sec)

# lli_value = lli(t=time_in_sec)
# lli_cyc_value = lli_cyc(t=time_in_sec)
# lam_value = lam(t=time_in_sec)
# V_cell_value = V_cell(t=time_in_sec)

#%
eta_plot = np.linspace(0.0, 1, 20) * lengthscale_core
# c_c_out_value = c_c_out(t=time_in_sec, x=eta_plot)

#%%
I = solution["Current [A]"].entries

f = open("pybamm\mz_develop\output\eg_current.txt", "w")
np.savetxt(f, np.c_[time_in_sec, I], fmt='%1.6e', delimiter=', ')
f.close()

#%%
s_out = solution["X-averaged moving phase boundary location"](t=time_output)
f = open("pybamm\mz_develop\output\sce3_phase_location.txt", "w")
np.savetxt(f, np.c_[time_output, s_out], fmt='%1.6e', delimiter=', ')
f.close()

#%%
lam = solution["X-averaged loss of active material in positive electrode (MZ)"](t=time_output)
f = open("pybamm\mz_develop\output\sce3_lam.txt", "w")
np.savetxt(f, np.c_[time_output, lam], fmt='%1.6e', delimiter=', ')
f.close()


#%%
s_out = solution.cycles[0].steps[2]["X-averaged moving phase boundary location"].entries

c_c_out = solution.cycles[0].steps[2]["X-averaged positive core concentration [mol.m-3]"].entries


tps = [1, 180, 360]
for tp in tps:
    th = tp/60
    path = "pybamm\mz_develop\output\eg_c_c_" + str(int(th)) + "h.txt"
    f = open(path, "w")
    np.savetxt(f, np.c_[eta_plot * s_out[tp], c_c_out[:,tp]], fmt='%1.6e', delimiter=', ')
    f.close()
    
# tps = [0, 10800, 21600]
# for tp in tps:
#     th = tp/3600
#     path = "pybamm\mz_develop\output\c_c_" + str(int(th)) + "h.txt"
#     f = open(path, "w")
#     np.savetxt(f, np.c_[eta_plot * s_out(t=tp), c_c_out(tp, x=eta_plot)], fmt='%1.6e', delimiter=', ')
#     f.close()

#%%
lli = solution["LLI"].entries
f = open("pybamm\mz_develop\output\eg_lli.txt", "w") 
np.savetxt(f, np.c_[time_in_sec, lli], fmt='%1.6e', delimiter=', ')
f.close()

#%%
lli_cyc = solution["LLI_cyc"](t=time_output)
f = open("pybamm\mz_develop\output\sce3_lli_cyc.txt", "w")
np.savetxt(f, np.c_[time_output, lli_cyc], fmt='%1.6e', delimiter=', ')
f.close()

#%%
soc = solution["SoC"](t=time_output)
f = open("pybamm\mz_develop\output\sce3_soc.txt", "w")
np.savetxt(f, np.c_[time_output, soc], fmt='%1.6e', delimiter=', ')
f.close()

#%%
shelleta = solution["X-averaged PE shell layer overpotential [V]"](t=time_output)
f = open("pybamm\mz_develop\output\sce3_shelloverpoten.txt", "w")
np.savetxt(f, np.c_[time_output, shelleta], fmt='%1.6e', delimiter=', ')
f.close()

#%%
c_p_surf = solution["X-averaged positive core surface concentration"](t=time_output)
c_n_surf = solution["X-averaged negative particle surface concentration"](t=time_output)

f = open("pybamm\mz_develop\output\sce3_c_surf.txt", "w")
np.savetxt(f, np.c_[time_output, c_n_surf, c_p_surf], fmt='%1.6e', delimiter=', ')
f.close()

#%%
M_p = solution["Total cyclable lithium in positive electrode [mol]"](t=time_output)
M_n = solution["Total cyclable lithium in negative electrode [mol]"](t=time_output)
M = solution["Total cyclable lithium in particles [mol]"](t=time_output)

f = open("pybamm\mz_develop\output\sce3_M.txt", "w")
np.savetxt(f, np.c_[time_output, M_p, M_n, M], fmt='%1.6e', delimiter=', ')
f.close()

#%%
f = open("pybamm\mz_develop\output\V_cell.txt", "w")
np.savetxt(f, np.c_[time_in_sec, V_cell_value], fmt='%1.6e', delimiter=', ')
f.close()

#%%
plot_times = np.linspace(0, max(time_in_sec), 11)

plt.figure(figsize=(10, 5))
cmap = plt.get_cmap("inferno")
for i, tp in enumerate(plot_times):
    color = cmap(float(i) / len(plot_times))
    plt.plot(
        eta_plot * s_out(t=tp),
        c_c_out(tp, x=eta_plot),
        "-",
        color=color,
        label="core" if i == 0 else ""
    )
    
# plt.ylim(0, 1.0)

plt.xlabel("r", fontsize=16)
plt.ylabel("c", fontsize=16)
plt.legend()
plt.show()