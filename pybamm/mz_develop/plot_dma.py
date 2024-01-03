time_in_sec = solution["Time [s]"].entries
L_sei = solution["SEI thickness [m]"].entries

#%%
plt.figure(figsize=(8, 6))
# markerfacecolor
plt.plot(time_in_sec, L_sei, 'o', mfc='none')