
"""
Calculates the VAF for a silver run, assuming all atoms are silver atoms, and that velocities
are given in Ang/ps, as standard in GULP.
By Leonid Kahle (Spring 2018)
"""

import json

def get_vaf(jsonfile, output, stepsize_t=1, stepsize_tau=1, t_end=5.0, init_time=0):

    #~ gs = gridspec.GridSpec(1,1, left=0.15, bottom=0.15, right=0.95, top=0.95)
    #~ fig = plt.figure()
    #~ ax = fig.add_subplot(gs[0])

    print 'Reading velocities...'
    with open(jsonfile) as f:
        traj_dict = json.load(f)
        velocities = traj_dict['velocities']
        times = traj_dict['times']

    print 'Calculating VAF(t)...'
    # I estimate the timestep, assuming sampling is equidistant in time:
    timestep = (times[-1]  - times[0]) / len(times)
    tau_index_start = int(init_time/timestep) # If I don't start averaging from beginning
    t_index_end = int(float(t_end) / stepsize_t / timestep )
    tau_index_end = len(velocities) - stepsize_t*t_index_end # avoid index out of bond
    # Getting number of atoms
    nat = len(velocities[0])
    # instantiating list storing VAF
    vaf = [0.0]*(t_index_end)
    # I quickly check that I actually have something to average over:
    if tau_index_end < tau_index_start:
        raise RuntimeError("My starting index is larger than final index for time average\n"
                "Too large init time?")
    # Normalization
    norm_factor = 1.0 / len(range(tau_index_start, tau_index_end, stepsize_tau)) / nat / 3
    # I'm looping over the time t
    for t_index in range(t_index_end):
        # Here I perform the ensemble/time average, which means that I integrate also over tau
        for tau_index in range(tau_index_start, tau_index_end, stepsize_tau):
            # I directly take mean over atoms
            for at in range(nat):
                # And directions:
                for d in range(3):
                    # Also, I directly normalize with norm_factor, this gives me the results in a temperature K
                    vaf[t_index] += norm_factor*velocities[tau_index+stepsize_t*t_index][at][d]*velocities[tau_index][at][d]
    print 'Writing results to {}...'.format(output)
    with open(output, 'w') as f:
        f.write('# Time   VAF(t)\n')
        for idx, vaf_of_t in enumerate(vaf):
            f.write('{:.6f}  {:.6f}\n'.format(stepsize_t*timestep*idx, vaf_of_t))
    print "Done"
    #~ ax.plot(stepsize_t*timestep*np.arange(t_index_end)+t_index_start, vaf)
    #~ plt.xlabel('Time [ps]')
    #~ plt.ylabel(r'VAF $\quad M_{Ag} k_B^{-1}\langle v(t)v(0)\rangle \quad [K]$')
    #~ plt.legend(loc=1)
    #~ plt.show()






if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('jsonfile', help='The json output, produces with parser.py'
        ' from a gulp trajectory')
    parser.add_argument('-o', '--output', help='The output file to write, defaults to vaf.dat', default='vaf.dat')
    parser.add_argument('--init-time', type=float, help='The time (in ps) to start sampling from, by default that is 0', default=0.0)
    parser.add_argument('--stepsize-t', type=int, help='The stepsize on the time index, by default the sampling time is used',
        default=1)
    parser.add_argument('--stepsize-tau', type=int, help='The stepsize over the time averaging of the trajectory. This can be safely set to a higher values'
            ' than 1 (default), if this function is too slow', default=1)
    parser.add_argument('-e', '--t-end', help='Until that time in ps the VAF is calculated, defaults to 5.', type=float, default=5.0)
    # Getting command line arguments:
    parsed_args = parser.parse_args()
    # Calling main function
    get_vaf(**vars(parsed_args))
