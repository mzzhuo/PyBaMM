#%%

c_c_out = solution["X-averaged positive core concentration [mol.m-3]"]
c_s_out = solution["X-averaged positive shell concentration [mol.m-3]"]
s_out = solution["X-averaged moving phase boundary location"]

c_tot = solution["Total lithium [mol]"]

c_shared_out = solution["X-averaged shared concentration at core-shell interface [mol.m-3]"]

#%
timescale = solution.timescale_eval
lengthscale_core = solution.length_scales_eval["positive core"]
lengthscale_shell = solution.length_scales_eval["positive shell"]

#%
time_in_sec = solution.t * timescale

eta_plot = np.linspace(0.0, 1, 30) * lengthscale_core
chi_plot = np.linspace(0.0, 1, 20) * lengthscale_shell

s_out_value = s_out(t=time_in_sec)
c_shared_out_value = c_shared_out(t=time_in_sec)

c_tot_value = c_tot(t=time_in_sec)

c_c_out_value = c_c_out(t=time_in_sec, x=eta_plot)
c_s_out_value = c_s_out(t=time_in_sec, x=chi_plot)

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
    plt.plot(
        ( chi_plot / lengthscale_shell * (1 - s_out(t=tp)) + s_out(t=tp)) * lengthscale_shell,
        c_s_out(tp, x=chi_plot),
        "--",
        color=color,
        label="shell" if i == 0 else ""
    )
    
# plt.ylim(0, 1.0)

plt.xlabel("r", fontsize=16)
plt.ylabel("c", fontsize=16)
plt.legend()
plt.show()