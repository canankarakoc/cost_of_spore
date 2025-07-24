import config
import pickle
#from scipy.stats import uniform
import numpy
from scipy.integrate import odeint, solve_ivp
import simulation_utils

import signal
import time


import matplotlib.pyplot as plt

numpy.random.seed(123456789)


def handler(signum, frame):
    raise Exception("Time's up!")

signal.signal(signal.SIGALRM, handler)


efficiency_sim_path = '%ssim.pickle' % config.data_directory



d_max = simulation_utils.d_max
n_s_0 = simulation_utils.n_s_0
r_max_spore = simulation_utils.r_max_spore
y_v = simulation_utils.y_v
y_s = simulation_utils.y_s


n_iter = 1000

t_span_24 = (0.0, 24.0)
t_span_48 = (0.0, 48.0)

# initial conditions
r_0_uniform = 10**numpy.random.uniform(1, 4, size=n_iter)
n_v_0_uniform = 10**numpy.random.uniform(4, 8, size=n_iter)

# growth parameters
# r_max = 2 (hr^-1) ==> cellular division every 30 minutes
r_max_uniform = 10**numpy.random.uniform(-1, 2, size=n_iter)
k_uniform = 10**numpy.random.uniform(1, 7, size=n_iter)

# spore initiation parameters
#s_uniform = 10**numpy.random.uniform(-9, -2, size=n_iter)
s_uniform = 10**numpy.random.uniform(-1, 2, size=n_iter)
sigma_uniform = 10**numpy.random.uniform(1, 5, size=n_iter)

# spore growth parameters
k_spore_uniform = 10**numpy.random.uniform(1, 7, size=n_iter)


# range of spore costs

y_s_range = numpy.logspace(-2, 2, base=10, num=20)
y_s_range *= y_v



def plot_spore_response():

    r_0_range = numpy.logspace(-2, 4, base=10, num=1000)

    fig, ax = plt.subplots(figsize=(4,4))

    for n in range(40):

        s = s_uniform[n]
        sigma = sigma_uniform[n]

        spore_response_n = simulation_utils.d_max/ (1 + numpy.exp(s*(r_0_range-sigma)))

        print(spore_response_n)

        ax.plot(r_0_range, spore_response_n, ls='-', c='dodgerblue', lw=2, alpha=0.6)


    ax.set_xscale('log', basex=10)

    ax.set_xlabel("Supplied resource concentration, " + r'$R_{0}$', fontsize = 12)
    ax.set_ylabel("Spore formation initiation rate, " + r'$f(R_{0})$', fontsize = 12)

    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig_name = "%sspore_response.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def run_simulation():

    efficiency_dict = {}

    for y_s_i in y_s_range:

        print(y_s_i/y_v)

        spore_efficiency_24_all = []
        spore_efficiency_48_all = []
        for n in range(n_iter):

            if (n+1)%100 == 0:
                print("%d simulations complete....." % (n+1))

            r_0 = r_0_uniform[n]
            n_v_0 = n_v_0_uniform[n]
            r_max = r_max_uniform[n]
            s = s_uniform[n]
            sigma = sigma_uniform[n]
            k = k_uniform[n]
            k_spore = k_spore_uniform[n]

            state_0 = [r_0, n_v_0, n_s_0]

            # 10 seconds per integration
            signal.alarm(10)

            try:

                params = (y_v, y_s_i, r_max, k, d_max, s, sigma, r_max_spore, k_spore)
                ivp_24 = solve_ivp(simulation_utils.growth_cycle_spores, t_span_24, state_0, args=params)
                ivp_48 = solve_ivp(simulation_utils.growth_cycle_spores, t_span_48, state_0, args=params)

                spore_efficiency_24 = ivp_24.y[:,-1][2] / (ivp_24.y[:,-1][1] + ivp_24.y[:,-1][2])
                spore_efficiency_48 = ivp_48.y[:,-1][2] / (ivp_48.y[:,-1][1] + ivp_48.y[:,-1][2])

                #print(spore_efficiency_24)
                spore_efficiency_24_all.append(spore_efficiency_24)
                spore_efficiency_48_all.append(spore_efficiency_48)

                #time.sleep(11)


            except Exception:
                #print('out of time')
                continue 

            #print(s, sigma, spore_efficiency_24)

        
        efficiency_dict[y_s_i] = {}
        efficiency_dict[y_s_i][24] = spore_efficiency_24_all
        efficiency_dict[y_s_i][48] = spore_efficiency_48_all


    with open(efficiency_sim_path, 'wb') as outfile:
        pickle.dump(efficiency_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)





#print(spore_efficiency_24_all)
    
#print("Median 24 hrs. %f" % round(numpy.median(spore_efficiency_24_all), 3) )
#print("Median 48 hrs. %f" % round(numpy.median(spore_efficiency_48_all), 3) )



def plot_spore_costs():

    sim_dict = pickle.load(open(efficiency_sim_path, "rb"))

    fig = plt.figure(figsize = (8, 4))
    fig.subplots_adjust(bottom= 0.1, wspace=0.15)

    for t_idx, t in enumerate([24, 48]):

        ax = plt.subplot2grid((1, 2), (0, t_idx))

        #efficiency_all = [sim_dict[y_s_i]  for y_s_i in y_s_range]

        lower_quantile_all = []
        upper_quantile_all = []
        mean_all = []
        for y_s_i in y_s_range:

            efficiency_all = numpy.asarray(sim_dict[y_s_i][t])
            efficiency_all = numpy.sort(efficiency_all)

            #print(efficiency_all[int(n_iter*0.025)])

            mean_all.append(numpy.median(efficiency_all))

            #print(numpy.quantile(efficiency_all, 0.975))
            lower_quantile_all.append(numpy.quantile(efficiency_all, 0.25))
            upper_quantile_all.append(numpy.quantile(efficiency_all, 0.75))
        

        mean_all = numpy.asarray(mean_all)
        upper_quantile_all = numpy.asarray(upper_quantile_all)
        lower_quantile_all = numpy.asarray(lower_quantile_all)

        cost_range_all = y_s_range

        ax.fill_between(y_v/cost_range_all, lower_quantile_all, upper_quantile_all, alpha=0.2, color='dodgerblue', zorder=1)

        ax.plot(y_v/cost_range_all, mean_all, c='k', ls='-', lw=3, zorder=2)

        ax.set_ylim([0,1])
        ax.set_xscale('log', basex=10)
        ax.set_title('%d hours' % t)


        ax.set_xlabel("Ratio of spore and cell costs", fontsize = 12)
        ax.set_ylabel("Sporulation efficiency", fontsize = 12)


    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig_name = "%stest_efficiency_range.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()






def plot_function():


    sim_dict = pickle.load(open(efficiency_sim_path, "rb"))

    spore_efficiency_24_all = sim_dict[23357214.690901212][24]
    spore_efficiency_48_all = sim_dict[23357214.690901212][48]



    fig = plt.figure(figsize = (8, 4))
    fig.subplots_adjust(bottom= 0.1, wspace=0.15)

    ax_24 = plt.subplot2grid((1, 2), (0, 0))
    ax_48 = plt.subplot2grid((1, 2), (0, 1))


    ax_24.hist(spore_efficiency_24_all, 50, density=True, histtype='step', facecolor='darkblue', alpha=0.75)
    ax_24.set_xlabel("Sporulation efficiency", fontsize = 12)
    ax_24.set_ylabel("Probability", fontsize = 12)
    ax_24.set_title('24 hours')


    ax_48.hist(spore_efficiency_48_all, 50, density=True, histtype='step', facecolor='darkblue', alpha=0.75)
    ax_48.set_xlabel("Sporulation efficiency", fontsize = 12)
    ax_48.set_ylabel("Probability", fontsize = 12)
    ax_48.set_title('48 hours')


    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig_name = "%stest_efficiency.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



#run_simulation()

#plot_spore_response()


plot_spore_costs()
plot_function()
