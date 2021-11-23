from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

import argparse

def make_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--coords')
	parser.add_argument('--top')
	parser.add_argument('--xml')

	return parser

parser = make_parser()


if __name__ == '__main__':
	args = parser.parse_args()

	coordinate_file = args.coords
	topology_file = args.top

	output_file = args.xml

	crd = AmberInpcrdFile(coordinate_file)
	prmtop = AmberPrmtopFile(topology_file)

	system = prmtop.createSystem(nonbondedMethod = PME, nonbondedCutoff = 1 * nanometer,
		constraints = HBonds)

	system_xml = XmlSerializer.serialize(system)

	with open(output_file, 'w') as f:
		f.write(system_xml)

	exit()

	
