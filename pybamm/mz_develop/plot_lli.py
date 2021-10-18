#%%

c_c_out = solution["X-averaged positive core concentration [mol.m-3]"]
s_out = solution["X-averaged moving phase boundary location"]
c_shared_out = solution["X-averaged shared concentration at core-shell interface [mol.m-3]"]
# c_tot = solution["Total lithium [mol]"]
lli = solution["Loss of lithium inventory [%]"]
lam = solution["X-averaged loss of active material in positive electrode (MZ)"]

V_cell = solution["Terminal voltage [V]"]

Q = solution.cycles[0]["Discharge capacity [A.h]"]
#%
timescale = solution.timescale_eval
lengthscale_core = solution.length_scales_eval["positive core"]
lengthscale_shell = solution.length_scales_eval["positive shell"]

#%
time_in_sec = solution.t * timescale



s_out_value = s_out(t=time_in_sec)
c_shared_out_value = c_shared_out(t=time_in_sec)

lli_value = lli(t=time_in_sec)
lam_value = lam(t=time_in_sec)
V_cell_value = V_cell(t=time_in_sec)

eta_plot = np.linspace(0.0, 1, 30) * lengthscale_core
c_c_out_value = c_c_out(t=time_in_sec, x=eta_plot)

#%%
# f = open("pybamm\mz_develop\output\c_c.txt", "w")
# np.savetxt(f, np.c_[eta_plot/lengthscale_core, c_c_out_value], fmt='%1.6e', delimiter=', ')
# f.close()


#%%

tps = [0, 10800, 21600]
for tp in tps:
    th = tp/3600
    path = "pybamm\mz_develop\output\case3_c_c_" + str(int(th)) + "h.txt"
    f = open(path, "w")
    np.savetxt(f, np.c_[eta_plot * s_out(t=tp), c_c_out(tp, x=eta_plot)], fmt='%1.6e', delimiter=', ')
    f.close()


#%%
f = open("pybamm\mz_develop\output\case3_phase_location.txt", "w")
np.savetxt(f, np.c_[time_in_sec, s_out_value], fmt='%1.6e', delimiter=', ')
f.close()

#%%
f = open("pybamm\mz_develop\output\case3_lli.txt", "w")
np.savetxt(f, np.c_[time_in_sec, lli_value], fmt='%1.6e', delimiter=', ')
f.close()

#%%
f = open("pybamm\mz_develop\output\case3_lam.txt", "w")
np.savetxt(f, np.c_[time_in_sec, lam_value], fmt='%1.6e', delimiter=', ')
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