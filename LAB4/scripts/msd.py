
"""
Calculates the MSD for a silver run, assuming all atoms are silver atoms, and that velocities
are given in Ang as standard in GULP.
By Leonid Kahle (Spring 2018)
"""

import json
# some constants
# MASS of silver
MASS = 107.87
# Boltzmann constant in this unit system
kbT = 0.8314457e0 # in units of A, ps and atomic mass u
# I now set the starting uindex for time as a constant
#~ t_index_start = 0


def get_msd(jsonfile, output, stepsize_t=1, stepsize_tau=1, t_end=5.0, init_time=0):
    print 'Reading positions...'
    with open(jsonfile) as f:
        traj_dict = json.load(f)
        positions = traj_dict['positions']
        times = traj_dict['times']
    # I estimate the timestep, assuming sampling is equidistant in time:
    timestep = (times[-1]  - times[0]) / len(times)
    dt = timestep * stepsize_t
    tau_index_start = int(init_time / timestep) # If I don't start averaging from beginning
    t_index_end = int(float(t_end) / dt )
    tau_index_end = len(positions) - stepsize_t*t_index_end # avoid index out of bond
    nat = len(positions[0])
    msd = [0.0]*(t_index_end)
    # I quickly check that I actually have something to average over:
    if tau_index_end < tau_index_start:
        raise RuntimeError("My starting index is larger than final index for time average\n"
                "Too large init time?")
    print 'Calculating  MSD(t)...'
    # Normalization 
    norm_factor = 1.0 / nat / 3 / len(range(tau_index_start, tau_index_end, stepsize_tau))
    #~ norm_factor = MASS / kbT / len(range(tau_index_start, tau_index_end, stepsize_tau))
    # I'm looping over the time t
    for t_index in range(t_index_end):
        # Here I perform the ensemble mean, which means that I integrate also over tau
        for tau_index in range(tau_index_start, tau_index_end, stepsize_tau):
            for at in range(nat):
                for d in range(3):
                    # Also, I directly normalize with norm_factor, this gives me the results in a temperature K
                    msd[t_index] += norm_factor*(positions[tau_index+stepsize_t*t_index][at][d] - positions[tau_index][at][d])**2

    #~ times = stepsize_t*timestep*np.arange(t_index_end)+t_index_start

    # Using least squares to get slope m and y-intercepts c
    #~ A = np.vstack([times[fit_index_start:], np.ones(len(times[fit_index_start:]))]).T
    #~ m, c = np.linalg.lstsq(A, msd[fit_index_start:])[0]
    #~ fit = m*times[fit_index_start:] + c
    # For convenience, I print some things
    #~ print '@ {}\n   slope: {:.4f}\n   residual: {:.4f}'.format(jsonfile, m, np.abs(msd[fit_index_start:] - fit).mean())

    print 'Calculating d/dt MSD(t)...'
    # differentiating in time:
    dmsd_dt = [0.0] * t_index_end
    for t_index in range(1, t_index_end-1):
        dmsd_dt[t_index] = ( msd[t_index+1] - msd[t_index-1] ) / (2*dt)
    dmsd_dt[t_index_end-1] = ( msd[t_index_end-1] - msd[t_index_end-2] ) / (dt)
    print 'Writing results to {} ...'.format(output)
    with open(output, 'w') as f:
        f.write('# Time t    MSD(t)    d/dt MSD(t)\n')
        for t_index, (msd, dmsd_dt) in enumerate(zip(msd, dmsd_dt)):
            f.write('{:.6f}   {:.6f}   {:.6f}\n'.format(dt*t_index, msd, dmsd_dt))
    print 'Done'






if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('jsonfile', help='The json output, produced with parser.py from a gulp trajectory')
    parser.add_argument('-o', '--output', help='write MSD to this output file. Defaults to msd.dat', default='msd.dat')

    parser.add_argument('--init-time', type=float, help='The time (in ps) to start sampling from, by default that is 0', default=0.0)
    parser.add_argument('--stepsize-t', type=int, help='The stepsize on the time index, by default the sampling time is used',
        default=1)
    parser.add_argument('--stepsize-tau', type=int, help='The stepsize over the time averaging of the trajectory. This can be safely set to a higher values'
            ' than 1 (default), if this function is too slow', default=1)
    parser.add_argument('-e', '--t-end', help='Until that time in ps the MSD is calculated, defaults to 5.', type=float, default=5.0)
    # Getting command line arguments:
    parsed_args = parser.parse_args()
    # Calling main function
    get_msd(**vars(parsed_args))
