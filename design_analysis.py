#! /usr/bin/env python2.7

## Sevy's version of DeepAnalysis
from optparse import OptionParser
from Bio.PDB import PDBParser
import subprocess, gzip
import multiprocessing as mp
from sevy_utils import Utils

if __name__ == '__main__':
	usage = "%prog [options] <pdb_files>"
	parser=OptionParser(usage)
	parser.add_option("--prefix",dest="prefix", help="Prefix for output files", default="out")
	parser.add_option("--native", '-n', dest="native", help="Native PDB for comparison", default="")
	parser.add_option("--format",dest="format", help="Sequence logo format", default="png")
	parser.add_option('--title',dest='title', help='Title for sequence logo', default='')
	parser.add_option('--res', '--resfile', dest='resfile',help='Resfile',default='')
	parser.add_option('--multiproc', '-m', dest='multiproc',help='Should I run this job on multiple processors?',default=False, action='store_true')
	parser.add_option('--debug',dest='debug',help='Should I print additional information', default=False, action='store_true')
	parser.add_option('--units',dest='units',help='What units for the weblogo? bits or probability, default bits',default='bits')
	parser.add_option('--weblogo-path',dest='weblogo_path',help='Path to the weblogo executable',default='weblogo')
	(options,args)= parser.parse_args()
	## Step 1: parse resfile to see which chains and residues I am designing
	if options.resfile == '':
		designable_residues = Utils.parse_designable_residues( args[0], [('*', '*')] )
	else:
		designable_residues = Utils.parse_resfile( options.resfile )

	## Step 2: find out the sequence of all of my models

	def partial_sequence_from_pdb( arg ):
		return Utils.sequence_from_pdb( arg, designable_residues )
		
	if options.multiproc:
		## Use all available CPUs
		pool = mp.Pool( processes=mp.cpu_count() )

		sequences = pool.map( partial_sequence_from_pdb, args )
	else:
		sequences = [ Utils.sequence_from_pdb( arg, designable_residues) for arg in args ]


	## Step 3: find out the sequence of my native
	if options.native == '':
		native_sequence = ''
	else:
		native_sequence = Utils.sequence_from_pdb( options.native, designable_residues )

	## Step 4: make fasta from my sequences
	with open(options.prefix+'.fasta', 'w') as out:
		modified_array = [ (name.split('/')[-1].split('.')[0], sequence) for name, sequence in zip( args, sequences ) ]
		fasta_str = Utils.fasta_from_sequences( modified_array )
		out.write( fasta_str )

	## Step 4: make .tab frequency file
	position_dict = {}
	## Changed by AMS 9.28.16 - need to get unique values from Utils.amino_acids
	## When I added sulfated tyrosine as Y this made Y a duplicate entry
	aa_list = sorted( list(set(Utils.amino_acids.values())) )
	for resno, chain in designable_residues:
		key = ','.join([chain, str(resno)])
		position_dict[ key ] = {aa : 0 for aa in aa_list}

	for sequence in sequences:
		for currentAA, des in zip(sequence, designable_residues):
			resno, chain = des
			key = ','.join([chain, str(resno)])
			position_dict[ key ][ currentAA ] += 1

	## Normalize position dictionary 
	for resno, chain in designable_residues:
		key = ','.join([chain, str(resno)])
		total = sum(position_dict[key].values())
		if total == 0:
			continue
		for aa in aa_list:
			position_dict[key][aa] /= float(total)

	## Write out to a file
	with open(options.prefix+'.tab', 'w') as out:
		out.write( '\t'.join( ['Chain,Position'] + aa_list ) + '\n' )
		for resno, chain in designable_residues:
			key = ','.join([chain, str(resno)])
			out.write( '\t'.join( [key] + [str(position_dict[key][aa]) for aa in aa_list] ) + '\n' )

	## Step 5: make sequence logo from my fasta
	weblogo_path = options.weblogo_path
	annotation = []
	if options.native == '':
		annotation = [chain+str(resno) for resno, chain in designable_residues ]
	else:
		for currentAA, des in zip(native_sequence, designable_residues):
			resno, chain = des
			annotation.append( chain+str(resno)+':'+currentAA )

	if options.units == 'bits':
		weblogo_command = [weblogo_path, '-f', options.prefix+'.fasta', '-o', options.prefix+'_seq_log.'+options.format, '-F', options.format, '-A', 'protein', '-U', 'bits', '-i', '1', '-s', 'large', '-t', options.title, '--annotate', ','.join(annotation), '-S', '4.32', '--errorbars', 'NO', '-c', 'chemistry', '--fineprint', 'sevy_analysis', '--scale-width', 'Yes', '--composition', 'equiprobable', '-n', '40', '-W', '30', '-y', 'bits']

	elif options.units == 'probability':
		weblogo_command = [weblogo_path, '-f', options.prefix+'.fasta', '-o', options.prefix+'_seq_log.'+options.format, '-F', options.format, '-A', 'protein', '-U', 'probability', '-i', '1', '-s', 'large', '-t', options.title, '--annotate', ','.join(annotation), '-S', '1.0', '--errorbars', 'NO', '-c', 'chemistry', '--fineprint', 'sevy_analysis', '--scale-width', 'Yes', '--composition', 'equiprobable', '-n', '40', '-W', '30', '-y', 'Frequency']

	if options.debug:
		print ' '.join( weblogo_command )

	output = subprocess.check_output(weblogo_command)
	# print output





