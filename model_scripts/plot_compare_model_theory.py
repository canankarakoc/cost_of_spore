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

r_max = 2
#k = 3.2*r_0
k = 3.2*r_0
s = 0.1
#sigma = 1000.2*r_0
sigma = 5
k_spore = 1*r_0



state_spores_0 = [r_0, n_v_0, n_s_0]


def plot_example():

    t_range = numpy.linspace(0, 48, num=1000, endpoint=True) 

    # example
    params_spores = (y_v, y_s, r_max, k, d_max, s, sigma, r_max_spore, k_spore)
    ivp_48 = solve_ivp(simulation_utils.growth_cycle_spores, (0.0, 48.0), state_spores_0, args=params_spores, dense_output=True)

    z_48 = ivp_48.sol(t_range)
    r_48 = z_48[0,:]
    n_v_48 = z_48[1,:]
    n_s_48 = z_48[2,:]

    # set small negative values to zero.
    n_v_48[(n_v_48 <0)] = 0
    n_s_48[(n_s_48 <0)] = 0

    efficiency_48 = n_s_48/(n_s_48 + n_v_48)
    stationry_s = n_s_48[-1]

    batch_prediction = (1 + (y_v/stationry_s) * r_0 + (n_v_0/y_v) - (stationry_s/y_s) )**-1
    chemostat_prediction = (r_max_spore*k)/(r_max*k_spore)

    fig, ax = plt.subplots(figsize=(4,4))

    last_zero_idx = numpy.where(efficiency_48==0)[0][-1]
    t_range_to_plot = t_range[last_zero_idx+1:]
    efficiency_48_to_plot = efficiency_48[last_zero_idx+1:]


    ax.plot(t_range, efficiency_48, lw=2, c='dodgerblue', label='Numerical solution')
    ax.axhline(y=batch_prediction, ls=':', c='k', label='Batch culture prediction')
    ax.axhline(y=chemostat_prediction, ls='--', c='k', label='Chemostat prediction')

    ax.set_xlim([0, 48])
    ax.set_ylim([1e-6, 2*chemostat_prediction])
    ax.set_yscale('log', basey=10)
    ax.legend(loc='upper left')
    ax.set_xlabel("Time (hrs.)", fontsize = 10)
    ax.set_ylabel("Sporulation efficiency", fontsize = 10)

    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig_name = "%scompare_model_theory.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





def plot_prediction_accuracy():

    fig = plt.figure(figsize = (8, 8))
    fig.subplots_adjust(bottom= 0.15)

    ax_monod_compare = plt.subplot2grid((2, 2), (0, 0), colspan=1)
    ax_monod_accuracy = plt.subplot2grid((2, 2), (0, 1), colspan=1)

    ax_signal_compare = plt.subplot2grid((2, 2), (1, 0), colspan=1)
    ax_signal_accuracy = plt.subplot2grid((2, 2), (1, 1), colspan=1)

    ax_monod_compare.text(-0.1, 1.05, "a)", fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_monod_compare.transAxes)
    ax_monod_accuracy.text(-0.1, 1.05, "b)", fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_monod_accuracy.transAxes)
    ax_signal_compare.text(-0.1, 1.05, "c)", fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_signal_compare.transAxes)
    ax_signal_accuracy.text(-0.1, 1.05, "d)", fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_signal_accuracy.transAxes)

    # predict efficiency for various monod constants.
    # simulation_utils.r_max_spore = 1/8
    r_max_spore_all = [1/16, simulation_utils.r_max_spore, 1/4]
    k_spore_all = numpy.logspace(-2.5, 2.5, num=50, endpoint=True, base=10)
    all_monod_data = []
    for r_max_spore_i_idx, r_max_spore_i in enumerate(r_max_spore_all):
        
        efficiency_48_monod_all = []
        batch_prediction_monod_all = []
        
        for k_spore_i in k_spore_all:
            
            # multiply by initial resources, by definition.
            params_spores_i = (y_v, y_s, r_max, k, d_max, s, sigma, r_max_spore_i, k_spore_i*r_0)
            ivp_48_i = solve_ivp(simulation_utils.growth_cycle_spores, (0.0, 48.0), state_spores_0, args=params_spores_i, dense_output=True)
            z_48_i = ivp_48_i.sol([48])

            r_48_i = z_48_i[0,:][0]
            n_v_48_i = z_48_i[1,:][0]
            n_s_48_i = z_48_i[2,:][0]
            efficiency_48_i = n_s_48_i/(n_s_48_i + n_v_48_i)
            batch_prediction_i = (1 + (y_v/n_s_48_i) * r_0 + (n_v_0/y_v) - (n_s_48_i/y_s) )**-1

            efficiency_48_monod_all.append(efficiency_48_i)
            batch_prediction_monod_all.append(batch_prediction_i)


        efficiency_48_monod_all = numpy.asarray(efficiency_48_monod_all)
        batch_prediction_monod_all = numpy.asarray(batch_prediction_monod_all)

        error = numpy.absolute((efficiency_48_monod_all-batch_prediction_monod_all)/efficiency_48_monod_all)
        #accuracy = 1 - error
        idx_only_positive = (batch_prediction_monod_all>0)

        all_monod_data.append(efficiency_48_monod_all[idx_only_positive])
        all_monod_data.append(batch_prediction_monod_all[idx_only_positive])

        ax_monod_compare.scatter(efficiency_48_monod_all[idx_only_positive], batch_prediction_monod_all[idx_only_positive], s=8, alpha=0.5, c=simulation_utils.colors_dict[r_max_spore_i_idx], label=r'$r_{\mathrm{max}}^{(s)} = $' + str(r_max_spore_i) + r'$\,\mathrm{hrs.^{-1}}$')
        ax_monod_accuracy.scatter(k_spore_all, error, s=8, c=simulation_utils.colors_dict[r_max_spore_i_idx], label=r'$r_{\mathrm{max}}^{(s)} = $' + str(r_max_spore_i) + r'$\,\mathrm{hrs.^{-1}}$')


    #### ax_monod_compare
    all_monod_data = numpy.concatenate(all_monod_data).ravel()
    ax_monod_compare.set_xlim([min(all_monod_data), max(all_monod_data)])
    ax_monod_compare.set_ylim([min(all_monod_data), max(all_monod_data)])
    ax_monod_compare.plot([min(all_monod_data), max(all_monod_data)], [min(all_monod_data), max(all_monod_data)], ls=':', lw=2, c='k', label='1:1')

    ax_monod_compare.set_xscale('log', basex=10)
    ax_monod_compare.set_yscale('log', basey=10)
    ax_monod_compare.legend(loc='upper left', fontsize=8)
    ax_monod_compare.set_xlabel("Simulated sporulation efficiency", fontsize = 10)
    ax_monod_compare.set_ylabel("Predicted sporulation efficiency", fontsize = 10)


    ### ax_monod_accuracy
    #ax_monod_accuracy.axhline(y=1, lw=2, ls=':', c='k')
    ax_monod_accuracy.set_xscale('log', basex=10)
    ax_monod_accuracy.set_yscale('log', basey=10)
    ax_monod_accuracy.legend(loc='upper right',  fontsize=8)
    ax_monod_accuracy.set_xlabel("Rescaled sporulation Monod constant, " + r'$K_{s}/R_{0}$' , fontsize = 10)
    ax_monod_accuracy.set_ylabel("Relative error of prediction", fontsize = 10)


    # manipulations of sporulation signal
    sigma_all = [10**-2, 10**0, 10**2]
    sigma_label_all = [r'$10^{-5}$', r'$10^{-3}$', r'$10^{-1}$']
    s_all = numpy.logspace(-2.5, 2.5, num=50, endpoint=True, base=10)


    all_signal_data = []
    for sigma_i_idx, sigma_i in enumerate(sigma_all):
        
        efficiency_48_signal_all = []
        batch_prediction_signal_all = []
        
        for s_i in s_all:
            
            # multiply by initial resources, by definition.
            params_spores_i = (y_v, y_s, r_max, k, d_max, s_i, sigma_i, r_max_spore, k_spore)
            ivp_48_i = solve_ivp(simulation_utils.growth_cycle_spores, (0.0, 48.0), state_spores_0, args=params_spores_i, dense_output=True)
            z_48_i = ivp_48_i.sol([48])

            r_48_i = z_48_i[0,:][0]
            n_v_48_i = z_48_i[1,:][0]
            n_s_48_i = z_48_i[2,:][0]
            
            efficiency_48_i = n_s_48_i/(n_s_48_i + n_v_48_i)
            batch_prediction_i = (1 + (y_v/n_s_48_i) * r_0 + (n_v_0/y_v) - (n_s_48_i/y_s) )**-1

            efficiency_48_signal_all.append(efficiency_48_i)
            batch_prediction_signal_all.append(batch_prediction_i)


        efficiency_48_signal_all = numpy.asarray(efficiency_48_signal_all)
        batch_prediction_signal_all = numpy.asarray(batch_prediction_signal_all)

        error = numpy.absolute((efficiency_48_signal_all-batch_prediction_signal_all)/efficiency_48_signal_all)
        #accuracy = 1 - error
        idx_only_positive = (batch_prediction_signal_all>0)

        all_signal_data.append(efficiency_48_signal_all[idx_only_positive])
        all_signal_data.append(batch_prediction_signal_all[idx_only_positive])

        signal_label = r'$R_{\mathrm{min}}/R_{0} = $' + sigma_label_all[sigma_i_idx]
        ax_signal_compare.scatter(efficiency_48_signal_all[idx_only_positive], batch_prediction_signal_all[idx_only_positive], s=8, alpha=0.5, c=simulation_utils.colors_dict[sigma_i_idx], label=signal_label)
        ax_signal_accuracy.scatter(s_all, error, s=8, c=simulation_utils.colors_dict[sigma_i_idx], label=signal_label)


    #### ax_signal_compare
    all_signal_data = numpy.concatenate(all_signal_data).ravel()
    ax_signal_compare.set_xlim([min(all_signal_data), max(all_signal_data)])
    ax_signal_compare.set_ylim([min(all_signal_data), max(all_signal_data)])
    ax_signal_compare.plot([min(all_signal_data), max(all_signal_data)], [min(all_signal_data), max(all_signal_data)], ls=':', lw=2, c='k', label='1:1')

    ax_signal_compare.set_xscale('log', basex=10)
    ax_signal_compare.set_yscale('log', basey=10)
    ax_signal_compare.legend(loc='upper left', fontsize=8)
    ax_signal_compare.set_xlabel("Simulated sporulation efficiency", fontsize = 10)
    ax_signal_compare.set_ylabel("Predicted sporulation efficiency", fontsize = 10)



    ### ax_monod_accuracy
    #ax_monod_accuracy.axhline(y=1, lw=2, ls=':', c='k')
    ax_signal_accuracy.set_xscale('log', basex=10)
    ax_signal_accuracy.set_yscale('log', basey=10)
    ax_signal_accuracy.legend(loc='upper right',  fontsize=8)
    ax_signal_accuracy.set_xlabel('"Sharpness" of spore formation initiation, ' + r'$\sigma$', fontsize = 10)
    ax_signal_accuracy.set_ylabel("Relative error of prediction", fontsize = 10)



    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig_name = "%sprediction_accuracy.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




# just hard code this....
#fig = plt.figure(figsize = (, 4))
#fig.subplots_adjust(bottom= 0.15)

#ax = plt.subplot2grid((1, 3), (0, 0), colspan=1)



plot_prediction_accuracy()


