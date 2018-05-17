#!/usr/bin/env python


"""
BACKGROUND:
This script will calculate the 3D radial distribution function from the custom json
WARNING: ONLY works in orthorhombic system
AUTHOR: Leonid Kahle (2018)
"""

import json, math

def calc_rdf(json_file, n_bins, init_time, stepsize_t, output='rdf.dat', plot=False):

    print 'Reading positions...'
    with open(json_file,'r') as f:
        all_dict = json.load(f)
    positions = all_dict['positions']
    times = all_dict['times']
    cell = all_dict['cell']
    
    if any([cell[i][j] for i in range(3) for j in range(3) if i!=j]):
        raise ValueError("Non-orthorhombic cell was given")

    lattice_vecs = [cell[i][i] for i in range(3)]
    inv_lattice_vecs = [1.0/cell[i][i] for i in range(3)]
    
    # The max_distance we need to consider is from the center of the rhombohedron to an edge!
    max_distance = min([cell[i][i]/2. for i in range(3)])
    nat = len(positions[0])

    dr = float(max_distance/n_bins)              # Bin width
    hist = [0]*(n_bins)

    if init_time >= times[-1]:
        raise Exception('The init_time is too large: no points are found!')

    for start_index, t in enumerate(times):
        if t > init_time:
            break

    print "Calculating g(r)..."

    # density of ideal gas:
    density_ideal = float(nat) / (lattice_vecs[0]*lattice_vecs[1]*lattice_vecs[2])
    
    # I want the entire histogram to integrate to NAT-1
    factor = 2. / float(nat) / len(range(start_index, len(positions), stepsize_t))

    for t_index in range(start_index, len(positions), stepsize_t):
    #~ for frac_pos_of_t in wrapped_fractional_coords:
        # loop over particle pairs
        positions_t = positions[t_index]
        for i in range(nat):
            for j in range(i+1,nat):
                distance_real_sq = 0.0
                for d in range(3):
                    dist_frac = inv_lattice_vecs[d]*(positions_t[i][d] - positions_t[j][d]) % 1.0
                    
                    if dist_frac < -0.5:
                        dist_frac += 1.
                    elif dist_frac > 0.5:
                        dist_frac = 1. - dist_frac
                    distance_real_sq += (dist_frac*lattice_vecs[d])**2
                bin = int(math.sqrt(distance_real_sq)/dr)
                if bin < n_bins:
                    hist[bin] += factor

    # Normalize to ideal gas density, and integrate
    integral = [0.0]*n_bins
    print 'Integrating g(r)...'
    for i in range(n_bins):
        vb = 4./3. * math.pi *( ((i+1)*dr)**3 - (i*dr)**3 ) # volume between r(i) and r(i+1)
        number_ideal =  vb * density_ideal # number of particles in a real gas
        # integral:
        if i:
            integral[i] = integral[i-1] + hist[i]
        hist[i] /= number_ideal
    print 'Writing results to {} ...'.format(output)
    with open(output, 'w') as f:
        f.write('# distance     RDF     integral\n')
        for idx, (r, i) in enumerate(zip(hist, integral)):
            f.write('{:.6f}  {:.6f}  {:.6f}\n'.format(dr*(idx+0.5), r, i))
    print 'Done'

#-------------------------------------------------------------------#

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    
    parser.add_argument('json_file', help='The json output, produces with parser.py'
        ' from a gulp trajectory')
    parser.add_argument('-n', '--n-bins', help='number of bins', type=int, default=100)
    parser.add_argument('-o', '--output', help='The output, defaults to rdft.dat', default='rdf.dat')
    parser.add_argument('--init-time', type=float, help='The time (in ps) to start sampling from, by default that is 0', default=0.0)
    parser.add_argument('--stepsize-t', type=int, help='The stepsize on the time index, by default the sampling time is used',
        default=1)

    parsed_args = parser.parse_args()

    calc_rdf(**vars(parsed_args))

