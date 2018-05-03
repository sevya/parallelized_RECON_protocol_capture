#!/usr/bin/env python2.7

import sys
from optparse import OptionParser, OptionGroup
from Bio.PDB import *
import warnings

usage = "%prog [options] --new_chain_order <comma-delimited list> <input_pdb> <output_pdb>"
parser=OptionParser(usage)
parser.add_option("--new_chain_order",dest="new_chain_order",help="the new order of chains. All chains must be specified")
parser.add_option("--new_chain_ids",dest="new_chain_ids",help="the new ids of the reordered chains. Type 'rename' for 'A,B,C...'. Default is conservation of chain id.",default="")

numbergroup = OptionGroup(parser,"Numbering Options","options which control the renumbering of the new pdb.")
numbergroup.add_option("-n",dest="start",help="residue number to start with, default is 1",default=1)
numbergroup.add_option("--norestart",dest="norestart",help="continue numbering through entire pdb instead of restarting at each chain",default=False, action="store_true")
numbergroup.add_option("--preserve",dest="preserve",help="preserve insertion code and heteroflags",default=False, action="store_true")
numbergroup.add_option("--keep_table",dest="keep_table",help="append non-ATOM lines, such as the Rosetta table, to the end of the pdb text file",default=False, action="store_true")
parser.add_option_group(numbergroup)

(options,args) = parser.parse_args()

if len(args) < 2:
	parser.error("specify input and output pdb_names")
	exit()

new_chain_order = options.new_chain_order.split(',')
if options.new_chain_ids == "rename":
        new_chain_ids = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','0','1','2','3','4','5','6','7','8','9']
elif options.new_chain_ids == "":
        new_chain_ids = new_chain_order
else:
        new_chain_ids = options.new_chain_ids.split(',')
if options.keep_table:
	other_lines = []
	pdb_text_file = open(args[0],'r')
	for line in pdb_text_file.readlines():
		if line.strip() == "":
			continue
		elif line.strip().split()[0] == "ATOM":
			continue
		elif line.strip().split()[0] == "TER":
			continue
		else:
			other_lines.append(line.strip())
	pdb_text_file.close()

parser = PDBParser(PERMISSIVE=1)
warnings.simplefilter('ignore',PDBExceptions.PDBConstructionWarning)

input_pdb = parser.get_structure('X',args[0])

output_structure_builder = StructureBuilder.StructureBuilder()
output_structure_builder.init_structure(args[1])
output_structure_builder.init_model(0)
output_struct = output_structure_builder.get_structure()

for model in output_struct:
	for new_order_id in new_chain_order:
		model.add(input_pdb[0][new_order_id])
		
id_iterator = 0
for model in output_struct:
	for chain in model:
		chain.id = new_chain_ids[id_iterator]
		id_iterator += 1

residue_id = int(options.start)

chain_id = ""
for residue in output_struct.get_residues():
	chain = residue.get_parent()
	if(chain_id != chain.get_id() and not options.norestart):
		chain_id = chain.get_id()
		residue_id=int(options.start)
	if(options.preserve):
		hetero = residue.id[0]
		insert = residue.id[2]
		residue.id=(hetero,residue_id,insert)
	else:
		residue.id=(' ',residue_id,' ')
	residue_id +=1
	
io = PDBIO()
io.set_structure(output_struct)
io.save(args[1])
print args[1],'successful'
if options.keep_table:
	pdb_text_file = open(args[1],'a')
	for line in other_lines:
		pdb_text_file.write(line)
		pdb_text_file.write("\n")

##os.system("/sb/meiler/scripts/capture_command.sh " + ' '.join([pipes.quote(x) for x in sys.argv]))
