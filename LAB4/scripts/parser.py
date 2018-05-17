"""
read the 'file.trg' and convert the data in a json
"""

import sys, json

def read_cell(output_file):
    with open(output_file,'r') as f:
        lines = f.readlines()
        cell = []
        for count,line in enumerate(lines, start=2):
            if "Cartesian lattice vectors (Angstroms)" in line:
                cell = [map(float, lines[count+i].split()) for i in range(3)]
                return cell



def parse_gulp_output(trajectory_fname, output_fname, json_fname, e_of_t_fname):
    """
    :param trajectory_fname: The name of the trajectory produced by GULP
    :param json_fname: The filename I will dump my result to, as a json
    :param e_of_t_fname: The filename I will write extensive properties of the system to
    """
    
    if not json_fname.endswith('.json'):
        json_fname+='.json'

    with open(trajectory_fname,'r') as f:
        lines = f.readlines()

    potential_energies = []
    kinetic_energies = []
    times = []
    temperatures = []
    positions = []
    velocities = []

    nat = int(lines[1].split()[0])

    for count,line in enumerate(lines, start=1):
        if 'Time/KE/E/T' in line:
            time, k, p, temp = map(float, lines[count].split())
            times.append(time)
            kinetic_energies.append(k)
            potential_energies.append(p)
            temperatures.append(temp)

        if 'Coordinates' in line:
            positions.append([
                map(float, lines[count+i].split()) for i in range(nat)])

            
        if 'Velocities' in line:
            velocities.append([
                map(float, lines[count+i].split()) for i in range(nat)])

    all_dict = {}
    all_dict['velocities'] = velocities
    all_dict['positions'] = positions
    all_dict['times'] = times
    all_dict['temperatures'] = temperatures
    all_dict['kinetic_energies'] = kinetic_energies
    all_dict['potential_energies'] = potential_energies
    all_dict['cell'] = read_cell(output_fname)

    with open(json_fname,'w') as f:
        json.dump(all_dict,f)

    with open(e_of_t_fname,'w') as f:
        f.write('# \t t \t T\tE\tK\n')
        for i in zip(times, temperatures, potential_energies, kinetic_energies):
            f.write( '%s\t%s\t%s\t%s\n' % (i[0],i[1],i[2],i[3]) )

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('trajectory', help='The trajectory file produced by GULP')
    parser.add_argument('output', help='The output file produced by GULP')
    parser.add_argument('json_output', help='The name of the JSON-output')
    parser.add_argument('-e', '--e-of-t', default='E_of_t.dat',
        help='The filename to write extensive quantities as a function of time, default: E_of_t.dat')

    parsed_args = parser.parse_args()
    parse_gulp_output(parsed_args.trajectory, parsed_args.output, parsed_args.json_output, parsed_args.e_of_t)
