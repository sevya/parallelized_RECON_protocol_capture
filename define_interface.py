#!/usr/bin/env python2.7

from Bio.PDB import PDBParser
from Bio.PDB import PDBExceptions
import warnings
import sys
import math
from optparse import OptionParser

standard_aa_names={"ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "PHE":"F", "GLY":"G", "HIS":"H", "ILE":"I", "LYS":"K", 
                   "LEU":"L", "MET":"M", "ASN":"N", "PRO":"P", "GLN":"Q", "ARG":"R", "SER":"S", "THR":"T", "VAL":"V",
                   "TRP":"W", "TYR":"Y", "TYS":"Y"}

'''
Calculates the minimal distance between any pair of atoms in two residues

@param pose BioPDB Model object to operate on
@param residue1 Residue 1 to measure distance. Residue is a tuple of (position, chain, aminoacid_id)
@param residue2 Residue 2 to measure distance. Residue is a tuple of (position, chain, aminoacid_id)

@return int value of the minimum pairwise distance of atoms in two residues
'''
def min_distance( pose, residue1, residue2 ):
	min_dist = -1
	seqpos1, chain1, resid1 = residue1
	seqpos2, chain2, resid2 = residue2
	for atom1 in pose[ chain1 ][ int(seqpos1) ]:
		if 'H' in atom1.get_name():
			continue
		for atom2 in pose[ chain2 ][ int(seqpos2) ]:
			if 'H' in atom2.get_name():
				continue
			if atom1 - atom2 < min_dist or min_dist == -1:
				min_dist = atom1 - atom2

	return min_dist

'''
Splits all the residues in a complex into two sides, divided by
which side of the protein-protein interface they are located on

@param pose BioPDB Model object to operate on
@param side1_chains List of chains that constitute side1
@param side2_chains List of chains that constitute side2

@return Two lists of residues in a tuple, where the first list is all
the residues on side1, and the second all the residues on side2
'''
def get_residue_sides( pose, side1_chains, side2_chains ):
	side1_residues = []
	for chain in side1_chains:
		for res in pose[ chain ]:
			seqpos = res.get_id()[1]
			resid  = standard_aa_names[ res.get_resname() ]
			side1_residues.append( ( str(seqpos), chain, resid ) )

	side2_residues = []
	for chain in side2_chains:
		for res in pose[ chain ]:
			seqpos = res.get_id()[1]
			resid  = standard_aa_names[ res.get_resname() ]
			side2_residues.append( ( str(seqpos), chain, resid ) )
	return side1_residues, side2_residues

'''
Calculates all residues with a heavy atom within X A of an atom
on the opposing chain, where X A is defined by nearby_atom_cutoff

@param pose BioPDB Model object to operate on
@param side1_chains List of chains that constitute side1
@param side2_chains List of chains that constitute side2
@param nearby_atom_cutoff distance cutoff

@return Two lists of residues in a tuple, where the first list is all
the residues on side1 that fall in the distance cutoff, 
and the second all the residues on side2 that fall in the distance cutoff.
'''
def sc_distance_cutoff( pose, side1_chains, side2_chains, nearby_atom_cutoff ):
	side1_residues, side2_residues = get_residue_sides( pose, side1_chains, side2_chains )
	
	side1_interface = []
	side2_interface = []
	for res1 in side1_residues:
		for res2 in side2_residues:
			if min_distance( pose, res1, res2 ) < nearby_atom_cutoff:
				side1_interface.append( res1 )
				side2_interface.append( res2 )
	return [ set( side1_interface ), set( side2_interface ) ]

'''
Calculates all residues with a Calpha within X A of an atom
on the opposing chain, where X A is defined by nearby_atom_cutoff

@param pose BioPDB Model object to operate on
@param side1_chains List of chains that constitute side1
@param side2_chains List of chains that constitute side2
@param ca_cutoff distance cutoff

@return Two lists of residues in a tuple, where the first list is all
the residues on side1 that fall in the distance cutoff, 
and the second all the residues on side2 that fall in the distance cutoff.
'''
def ca_distance_cutoff( pose, side1_chains, side2_chains, ca_cutoff ):
	side1_residues, side2_residues = get_residue_sides( pose, side1_chains, side2_chains )
	
	side1_interface = []
	side2_interface = []
	for res1 in side1_residues:
		for res2 in side2_residues:
			atom1 = pose[ res1[1] ][ int(res1[0]) ][ 'CA' ]
			atom2 = pose[ res2[1] ][ int(res2[0]) ][ 'CA' ]
			if atom1 - atom2 < ca_cutoff:
				side1_interface.append( res1 )
				side2_interface.append( res2 )
	return [ set( side1_interface ), set( side2_interface ) ]

'''
Calculates all residues with a Cbeta within X A of an atom
on the opposing chain, where X A is defined by nearby_atom_cutoff

@param pose BioPDB Model object to operate on
@param side1_chains List of chains that constitute side1
@param side2_chains List of chains that constitute side2
@param ca_cutoff distance cutoff

@return Two lists of residues in a tuple, where the first list is all
the residues on side1 that fall in the distance cutoff, 
and the second all the residues on side2 that fall in the distance cutoff.
'''
def cb_distance_cutoff( pose, side1_chains, side2_chains, cb_cutoff ):

	side1_residues, side2_residues = get_residue_sides( pose, side1_chains, side2_chains )
	side1_interface = []
	side2_interface = []
	for res1 in side1_residues:
		for res2 in side2_residues:
			type1 = 'CB' if res1[2] != 'G' else 'CA'
			type2 = 'CB' if res2[2] != 'G' else 'CA'
			atom1 = pose[ res1[1] ][ int(res1[0]) ][ type1 ]
			atom2 = pose[ res2[1] ][ int(res2[0]) ][ type2 ]
			if atom1 - atom2 < cb_cutoff:
				side1_interface.append( res1 )
				side2_interface.append( res2 )
	return [ set( side1_interface ), set( side2_interface ) ]



'''
@file define_interface.py
@brief Generates a resfile defining a protein-protein interface.
Calculates the interface based on a side chain distance cutoff.
Can also be set to repack the residues on the opposing side of
the interface.
@author Alex Sevy (alex.sevy@gmail.com)
'''
if __name__ == '__main__':
	usage = "%prog [options] <pdb_file>"
	parser=OptionParser(usage)
	parser.add_option("--side1",dest="side1",help="the chains that make up one side of the interface (as a string, e.g. 'AB')", default="")
	parser.add_option("--side2",dest="side2",help="the chains that make up the other side of the interface (as a string, e.g. 'CD')", default="")
	parser.add_option("--nearby_atom_cutoff",dest="nearby_atom_cutoff",help="SC distance cutoff to define a residue as part of the interface. If any SC atom from a residue on one side is within this cutoff of a residue on the other side it's considered to be in the interface. Default=7.0",default=7.0)
	parser.add_option("--output",dest="output",default="out", help="Output name for resfile")
	parser.add_option("--native",dest="native",action='store_true',help='Just repack the residues on the side flagged "design side"', default=False)
	parser.add_option("--design-side",dest="design_side",default='1',help="Side of interface to design - either 1 or 2. Defaults to 1.")
	parser.add_option("--repack",dest="repack",action='store_true',default=False, help='Repack side of the interface not being designed')
	(options,args)= parser.parse_args()
	
	if len(args) < 1:
		parser.error('specify a pdb file or use -h')
	elif len(args) > 1:
		print 'Warning: only the first pdb is considered by this script. The rest will be ignored'

	warnings.simplefilter('ignore',PDBExceptions.PDBConstructionWarning)	

	print 'Processing',args[0]
	parser = PDBParser()
	structure = parser.get_structure( 'X', args[ 0 ] )

	## As implemented thid script only considers side chain distances. Can easily
	## be reconfigured for Ca or Cb distances though
	side1, side2 = sc_distance_cutoff( structure[ 0 ], 
									   list( options.side1 ),
									   list( options.side2 ),
									   float(options.nearby_atom_cutoff) )
	
	side_dict = {
					'1' : sorted( side1, key=lambda res: ( res[1], int(res[0]) ) ),
					'2' : sorted( side2, key=lambda res: ( res[1], int(res[0]) ) )
				}

	opposing_side = '2' if options.design_side == '1' else '1'
	packer_aas = 'NATAA' if options.native else 'ALLAA'

	with open(options.output+'.resfile', 'w')  as out:
		out.write('NATRO\nEX 1 EX 2\nstart\n')
		
		## This iterator lets you put in multiple design sides
		for side in list(options.design_side):
			for seqpos, chain, resid in side_dict[ side ]:
				out.write( ' '.join([ seqpos, chain, packer_aas ])+'\n' )

		if options.repack:
			for seqpos, chain, resid in side_dict[ opposing_side ]:
				out.write( ' '.join([ seqpos, chain, 'NATAA' ])+'\n' )

