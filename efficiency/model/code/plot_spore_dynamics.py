import numpy
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp

import simulation_utils
import config



d_max = simulation_utils.d_max
n_s_0 = simulation_utils.n_s_0
r_max_spore = simulation_utils.r_max_spore
y_v = simulation_utils.y_v
y_s = simulation_utils.y_s


r_0 = 10**3
n_v_0 = 10**5

state_spores_0 = [r_0, n_v_0, n_s_0]
state_0 = [r_0, n_v_0]

r_max = 2
k = 3.2*r_0
s = 0.1
#sigma = 1000.2*r_0
sigma = 5

k_spore = 0.01*r_0



#t_span_48 = (0.0, 24.0)
t_span_48 = (0.0, 48.0)

#t_eval = numpy.linspace(0, 48, 1000)

params = (y_v, r_max, k)


params_spores = (y_v, y_s, r_max, k, d_max, s, sigma, r_max_spore, k_spore)
#ivp_48_spires = solve_ivp(simulation_utils.growth_cycle_spores, t_span_48, state_spores_0, args=params_spores, dense_output=True)


t = numpy.linspace(0, 48, 1000)


def plot_cr():

    ivp_48 = solve_ivp(simulation_utils.growth_cycle, t_span_48, state_0, args=params, dense_output=True)
    z = ivp_48.sol(t)

    fig, ax = plt.subplots(figsize=(4,4))

    ax.plot(t, z[0,:], ls='-', c='k', lw=2, label='Resource concentration')
    ax.plot(t, z[1,:], ls='-', c='b', lw=2, label='Cell density')

    ax.set_yscale('log', basey=10)

    ax.set_xlabel("Time (hrs)", fontsize = 12)
    ax.set_ylabel("Variable", fontsize = 12)

    ax.set_ylim([1e-3, 1e10])

    ax.legend(loc='lower left')


    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig_name = "%stest_cr.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




def plot_cr_spores():

    ivp_48 = solve_ivp(simulation_utils.growth_cycle_spores, t_span_48, state_spores_0, args=params_spores, dense_output=True)
    z = ivp_48.sol(t)

    fig, ax = plt.subplots(figsize=(4,4))

    ax.plot(t, z[0,:], ls='-', c='k', lw=2, label='Resource concentration')
    ax.plot(t, z[1,:], ls='-', c='b', lw=2, label='Cell density')
    ax.plot(t, z[2,:], ls='-', c='r', lw=2, label='Spore density')

    ax.set_yscale('log', basey=10)

    ax.set_xlabel("Time (hrs)", fontsize = 12)
    ax.set_ylabel("Variable", fontsize = 12)

    ax.set_ylim([1e-3, 1e10])

    ax.legend(loc='lower left')

    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig_name = "%scr_spores.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



    # save data for Canan as flat file
    data_file_path = '%sconsumer_resource_spore_model.tsv' % config.data_directory
    data_file = open(data_file_path, 'w')
    
    time_str =  'hours\t' + ",".join(str(x) for x in t)
    resource_str = 'resource_concentration\t' + ",".join(str(x) for x in z[0,:].tolist())
    cell_str = 'cell_concentration\t' + ",".join(str(x) for x in z[1,:].tolist())
    spore_str = 'spore_concentration\t' + ",".join(str(x) for x in z[2,:].tolist())

    data_file.write(time_str)
    data_file.write('\n')
    
    data_file.write(resource_str)
    data_file.write('\n')

    data_file.write(cell_str)
    data_file.write('\n')

    data_file.write(spore_str)
    data_file.write('\n')
    
    data_file.close()






def plot_cost_vs_efficiency():

    y_s_range = numpy.logspace(-3, 3, base=10, num=50)
    y_s_range *= y_s

    x_axis = y_v/y_s_range

    efficiency_24_all = []
    efficiency_48_all = []

    for y_s_i in y_s_range:

        params_spores_i = (y_v, y_s_i, r_max, k, d_max, s, sigma, r_max_spore, k_spore)


        ivp_48 = solve_ivp(simulation_utils.growth_cycle_spores, t_span_48, state_spores_0, args=params_spores_i, dense_output=True)
        z_24 = ivp_48.sol([24])
        z_48 = ivp_48.sol([48])

        efficiency_24 = z_24[2,:][0] / (z_24[2,:][0] + z_24[1,:][0])
        efficiency_48 = z_48[2,:][0] / (z_48[2,:][0] + z_48[1,:][0])

        efficiency_24_all.append(efficiency_24)
        efficiency_48_all.append(efficiency_48)



    fig, ax = plt.subplots(figsize=(4,4))

    #ax.plot(x_axis, efficiency_24_all, ls='--', c='k', lw=2, label='24 hrs.')
    ax.plot(x_axis, efficiency_48_all, ls='-', c='k', lw=2)

    ax.set_xscale('log', basex =10)
    ax.set_yscale('log', basey=10)

    ax.set_xlabel("Ratio of energetic costs," + r'$ \frac{\mathrm{Spore}}{\mathrm{Cell}}$', fontsize = 12)
    ax.set_ylabel("Sporulation efficiency", fontsize = 12)

    ax.axvline(x=1/20, ls=':', c='k', lw=1.5, label='Empirical estimate')
    ax.legend(loc='lower left')

    #ax.set_ylim([1e-3, 1e10])
    #ax.legend(loc='lower left')


    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig_name = "%scr_spores_efficiency.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


    # save data for Canan as flat file
    data_file_path = '%sefficiency_data.tsv' % config.data_directory
    data_file = open(data_file_path, 'w')
    time_str =  'spore_vs_cell_atp_ratio\t' + ",".join(str(x) for x in x_axis)
    efficiency_str = 'sporulation_efficiency\t' + ",".join(str(x) for x in efficiency_48_all)
    data_file.write(time_str)
    data_file.write('\n')
    data_file.write(efficiency_str)
    data_file.close()




plot_cr_spores()
plot_cost_vs_efficiency()




#result = odeint(growth_cycle_spores, state_spores_0, t, args=(y_s,))