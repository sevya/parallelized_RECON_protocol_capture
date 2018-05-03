#!/usr/bin/env python2.7
import json
from operator import itemgetter
from find_interface_residues import *

def define_second_shell(structure, iface_res, second_shell_dist_cutoff):
	second_shell = set()
	all_atom_list =  Selection.unfold_entities(structure, 'A') # A for atoms
	residue_list = Selection.unfold_entities(structure, 'R')
	ns = NeighborSearch(all_atom_list)

	for residue in iface_res:
		atom_list = Selection.unfold_entities(residue, 'A') # A for atoms
		for atom in atom_list:
			distance = 0
			center = atom.get_coord()
			neighbors = ns.search(center, second_shell_dist_cutoff,level='A')
			neighbor_list = Selection.unfold_entities(neighbors, 'R') # R for residues
			for neighbor in neighbor_list:
				second_shell.add(neighbor)
	return second_shell

def print_secondaryresfile(sec_resfile_name, structure_residues, epitope, paratope, secondary_paratope, cdrs):
	sec_resfile = open(sec_resfile_name, 'w')
	sec_resfile.write("NATRO\n")
	sec_resfile.write("start\n")
	residue_chain_dict = {}
	for res in structure_residues:
		residue_chain_dict[res.get_id()[1]] = res.get_parent().get_id()
	for res in epitope:
		sec_resfile.write(" ".join([str(res), residue_chain_dict[res], "NATAA EX 1 EX 2", "\n"]))
	for res in cdrs:
		if res in paratope:
			continue
		sec_resfile.write(" ".join([str(res), residue_chain_dict[res], "NATAA EX 1 EX 2","\n"]))
	for res in secondary_paratope:
		if res in paratope or res in cdrs:
			continue
		sec_resfile.write(" ".join([str(res), residue_chain_dict[res], "NATAA EX 1 EX 2","\n"]))

def print_correspondence_entity_resfile(corr_filename, entity_filename, analysis_corr_filename, structure_residues, paratope, cdrs):
	corr_file = open(corr_filename,'w')
	analysis_corr_file = open(analysis_corr_filename,'w')
	entity_file = open(entity_filename, 'w')
	entity_file.write(str(len(paratope))+"\n")
	entity_file.write("ALLAA EX 1 EX 2\n")
	entity_file.write("start\n")
	residue_chain_dict = {}
	for res in structure_residues:
		residue_chain_dict[res.get_id()[1]] = res.get_parent().get_id()
	for i,res in enumerate(paratope):
		corr_file.write(" ".join([str(i+1),str(res),residue_chain_dict[res],"\n"]))
		entity_file.write(" ".join([str(i+1), residue_chain_dict[res], 'ALLAA EX 1 EX 2',"\n"]))
	analysis_residues = set()
	for res in paratope:
		analysis_residues.add(int(res)+1)
		analysis_residues.add(int(res))
		analysis_residues.add(int(res)-1)
	#for res in cdrs:
	#	analysis_residues.add(res)
	for i,res in enumerate(sorted(analysis_residues)):
		analysis_corr_file.write(" ".join([str(i+1),str(res),residue_chain_dict[res],"\n"]))
	
	
def print_out_pymol(pml_filename, structure_residues, epitope, paratope, secondary_paratope):
	pml_file = open(pml_filename,'w')
	pml_file.write("as ribbon\n")
	pml_file.write("util.cbc\n")
	residue_chain_dict = {}
	for res in structure_residues:
		residue_chain_dict[res.get_id()[1]] = res.get_parent().get_id()
	epitope_chain_res_dict = {}
	paratope_chain_res_dict = {}
	sec_paratope_chain_res_dict = {}
	for resid in epitope:
		try:
			epitope_chain_res_dict[residue_chain_dict[resid]].append(str(resid))
		except KeyError:
			epitope_chain_res_dict[residue_chain_dict[resid]] = [str(resid)]
	for resid in paratope:
		try:
			paratope_chain_res_dict[residue_chain_dict[resid]].append(str(resid))
		except KeyError:
			paratope_chain_res_dict[residue_chain_dict[resid]] = [str(resid)]
	for resid in secondary_paratope:
		try:
			sec_paratope_chain_res_dict[residue_chain_dict[resid]].append(str(resid))
		except KeyError:
			sec_paratope_chain_res_dict[residue_chain_dict[resid]] = [str(resid)]
	for chain in sorted(epitope_chain_res_dict):
		pml_file.write("show sticks, chain "+chain+" and resi "+"+".join(epitope_chain_res_dict[chain])+" and not name N+C+O\n")
	for chain in sorted(paratope_chain_res_dict):
		pml_file.write( "show sticks, chain "+chain+" and resi "+"+".join(paratope_chain_res_dict[chain])+" and not name N+C+O\n")
	for chain in sorted(sec_paratope_chain_res_dict):
		pml_file.write( "show lines, chain "+chain+" and resi "+"+".join(sec_paratope_chain_res_dict[chain])+" and not name N+C+O\n")
	pml_file.write("hide everything, elem H\n")
	pml_file.write("color atomic, not name c*\n")
	return



if __name__ == '__main__':
	usage = "%prog [options] <pdb_file>"
	parser=OptionParser(usage)
	parser.add_option("--side1",dest="pat1", default="",
								help="the chains that make up one side of the interface (as a \
								string, e.g. 'AB')")
	parser.add_option("--side2",dest="pat2", default="",
								help="the chains that make up the other side of the interface (as \
								a string, e.g. 'CD')")
	parser.add_option("--CB_dist_cutoff",dest="CB_dist_cutoff",
								help="Primary cutoff distance. All residues within this big cutoff \
								(CB_dist_cutoff) of residues on the other chain are evaluated \
								further. Default=10.0",default="10.0")
	parser.add_option("--nearby_atom_cutoff",dest="nearby_atom_cutoff",
								help="Secondary cutoff 1: iterate through all the side chain atoms \
								in the residue of interest and check to see if the distance to a \
								residue across the interface is less than the nearby atom cutoff \
								(nearby_atom_cutoff), if so then they are an interface residue. \
								Default=6",default="5.5")
	parser.add_option("--vector_angle_cutoff",dest="vector_angle_cutoff",
								help="If a residue does not pass Secondary cutoff 1, then two \
								vectors are drawn, a CA-CB vector and a vector from CB to a CB atom\
								 on the neighboring chain. The dot product between these two \
								vectors is then found and if the angle between them \
								(vector_angle_cutoff) is less than some cutoff then they are \
								classified as interface. The vector cannot be longer \
								than some other distance (vector_dist_cutoff). \
								Default=75 degrees",default="75")
	parser.add_option("--vector_dist_cutoff",dest="vector_dist_cutoff",default="9",
								help="maximal vector length. Default=9")
	parser.add_option("--second_shell_dist_cutoff","--sec_dist", dest="second_shell_dist_cutoff",
								help="Cutoff for residues near interface residues, measured by Cb-Cb \
								distance. Default=2",default="2")
	parser.add_option("--output",dest="output",default="corr", help="Output name for resfile and json \
								correspondence file")
	parser.add_option("--native",dest="native",action='store_true',default=False)
	parser.add_option("--design-side",dest="design",default='1',help="Side of interface to design - either 1 or 2. Defaults to 1.")
	parser.add_option("--repack",dest="repack",action='store_true',default=False, help='Repack side of the interface not being designed')
	(options,args)= parser.parse_args()
	
	warnings.simplefilter('ignore',PDBExceptions.PDBConstructionWarning)	
	
	print "Analyzing",len(args),"pdbs."

	group1_chains = []
	for chain in options.pat1:
		group1_chains.append(chain)
	group2_chains = []
	for chain in options.pat2:
		group2_chains.append(chain)

	side1_dict = {}
	side2_dict = {}
	side1_list = []
	side2_list = []
	for pdb in args:
		print 'Processing',pdb
		parser = PDBParser()
		structure = parser.get_structure('X',pdb)
		structure_residues = [residue for residue in structure.get_residues()]
		dist_shell = distance_cutoff(structure_residues, group1_chains, group2_chains, \
			options.nearby_atom_cutoff, options.CB_dist_cutoff)
		iface_shell = vector_cutoff(structure_residues, dist_shell, options.CB_dist_cutoff, \
			options.vector_dist_cutoff, options.vector_angle_cutoff, group1_chains, group2_chains)
		second_shell = define_second_shell(structure, iface_shell, float(options.second_shell_dist_cutoff ))
		name = pdb[:-4]
		for res in sorted(iface_shell):
			if res.get_parent().get_id() in group1_chains:
				letter = res.get_parent().get_id()
				number = res.get_id()[1]
				key = str(number)+' '+letter
				if key not in side1_list:
					side1_list.append(key)
				if name not in side1_dict:
					side1_dict[name] = []
				side1_dict[name].append([res.get_parent().get_id(), res.get_id()[1]])
			if res.get_parent().get_id() in group2_chains:
				letter = res.get_parent().get_id()
				number = res.get_id()[1]
				key = str(number)+' '+letter
				if name not in side2_dict:
					side2_dict[name] = []
				if key not in side2_list:
					side2_list.append(key)
				side2_dict[name].append([res.get_parent().get_id(), res.get_id()[1]])

		# for res in sorted(second_shell):
		# 	if res.get_parent().get_id() in group1_chains:
		# 		try:
		# 			epitope_dict[HA].add(res.get_id()[1])
		# 		except KeyError:
		# 			epitope_dict[HA] = set(res.get_id()[1])
		# 	if res.get_parent().get_id() in group2_chains:
		# 		secondary_paratope_set.add( res.get_id()[1])
	side1_list.sort(key=lambda a: (a.split()[1], int(a.split()[0])))
	side2_list.sort(key=lambda a: (a.split()[1], int(a.split()[0])))
	if options.design=='1':
		#json.dump(side1_dict, open(options.output+".json", 'w'))
		pml_string = ''
		with open(options.output+'.resfile', 'w')  as out:
			out.write('NATRO\nEX 1 EX 2\nstart\n')
			pml_string+='as cartoon\nselect designable, '
			# for key in side1_dict:
			# 	for letter, number in side1_dict[key]:
			# 		out.write(str(number)+' '+letter+' ALLAA EX 1 EX 2\n')
			# 		pml_string+='(resi '+str(number)+' and chain '+letter+') or '
			for key in side1_list:
				if options.native:
					out.write(key+' NATAA\n')
				else:
					out.write(key+' ALLAA\n')
			if options.repack:
				for key in side2_list:
					out.write(key+' NATAA\n')
				pml_string+='(resi '+key.split()[0]+' and chain '+key.split()[1]+') or '
		#with open(options.output+'.pml', 'w') as out:
		#	out.write(pml_string[:-3])
		#	out.write('\ncolor blue, designable\nshow sticks,designable')
			
			
	else:
		#json.dump(side2_dict, open(options.output+".json", 'w'))
		pml_string=''
		with open(options.output+'.resfile', 'w') as out:
			out.write('NATRO\nEX 1 EX 2\nstart\n')
			pml_string+='as cartoon\nselect designable, '
			# for key in side2_dict:
			# 	for letter, number in side2_dict[key]:
			# 		out.write(str(number)+' '+letter+' ALLAA EX 1 EX 2\n')
			# 		pml_string+='(resi '+str(number)+' and chain '+letter+') or '
			for key in side2_list:
				if options.native:
					out.write(key+' NATAA\n')
				else:
					out.write(key+' ALLAA\n')
				pml_string+='(resi '+key.split()[0]+' and chain '+key.split()[1]+') or '
			if options.repack:
				for key in side1_list:
					out.write(key+' NATAA\n')
		#with open(options.output+'.pml', 'w') as out:
		#	out.write(pml_string[:-3])
		#	out.write('\ncolor blue, designable\nshow sticks,designable')
	exit()

	paratope_set = set()
	secondary_paratope_set = set()
	epitope_set = set()
	cdr_res_set = set()
	
	epitope_dict = {}

	Ab_name = "tmp"
	for pdb in args:
		print pdb
		parser = PDBParser()
		structure = parser.get_structure('X',pdb)
		structure_residues = [residue for residue in structure.get_residues()]

		name_bits = pdb.strip().split('/')[-1].split('_')
		HA = name_bits[0]
		Ab_name = name_bits[3]
		Ab_pdbid = name_bits[4]

		cdr_res = [ int(line.strip().split()[0]) for line in open("/home/nannemdp/antibody-antigen_interactions/single_state_maturation/automated_workflow/template_prep/native_xtal_pieces/"+Ab_pdbid+"_Ab_cdr_res.list",'r').readlines() ]
		HA_length = 0
		for chain in structure[0]:
			for res in chain:
				if chain.get_id() in group1_chains and int(res.get_id()[1]) > HA_length:
					HA_length = res.get_id()[1]
		for res in cdr_res:
			cdr_res_set.add(int(res))
		
		iface_shell = distance_cutoff(structure_residues, group1_chains, group2_chains, options.nearby_atom_cutoff, options.CB_dist_cutoff)
		iface_shell = vector_cutoff(structure_residues, iface_shell, options.CB_dist_cutoff, options.vector_dist_cutoff, options.vector_angle_cutoff, group1_chains, group2_chains)
		
		second_shell = define_second_shell(structure, iface_shell, float(options.second_shell_dist_cutoff ))
		

		for res in sorted(iface_shell):
			if res.get_parent().get_id() in group1_chains:
				try:
					epitope_dict[HA].add(res.get_id()[1])
				except KeyError:
					epitope_dict[HA] = set([res.get_id()[1]])
			if res.get_parent().get_id() in group2_chains:
				paratope_set.add( res.get_id()[1])
		for res in sorted(second_shell):
			if res.get_parent().get_id() in group1_chains:
				try:
					epitope_dict[HA].add(res.get_id()[1])
				except KeyError:
					epitope_dict[HA] = set(res.get_id()[1])
			if res.get_parent().get_id() in group2_chains:
				secondary_paratope_set.add( res.get_id()[1])
	
	entity_name = Ab_name+'.entres'
	corr_filename = Ab_name+'.corr'
	analysis_corr_filename = Ab_name+'_analysis.corr'
	
	parser = PDBParser()
	structure = parser.get_structure('X',args[0])
	structure_residues = [residue for residue in structure.get_residues()]

	print_correspondence_entity_resfile( corr_filename, entity_name, analysis_corr_filename, structure_residues, sorted(paratope_set), sorted(cdr_res_set))
	print epitope_dict
	for HA in epitope_dict:
		print HA
		sec_resfile_name = HA+".2res"
		epitope_set = epitope_dict[HA]
		print_secondaryresfile(sec_resfile_name, structure_residues, sorted(epitope_set), sorted(paratope_set), sorted(secondary_paratope_set), sorted(cdr_res_set))

		pml_filename = Ab_name+"_"+HA+".pml"
		print_out_pymol(pml_filename, structure_residues, sorted(epitope_set), sorted(paratope_set), sorted(secondary_paratope_set)+sorted(cdr_res_set))
