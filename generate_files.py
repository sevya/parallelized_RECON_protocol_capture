#! /usr/bin/env python2.7

import sys, os, shutil
sys.path.append('/dors/meilerlab/apps/Linux2/x86_64/python/2.7.6/lib/python2.7/site-packages/')
from Bio import pairwise2
from Bio import SeqIO
import multiprocessing as mp
import time
import numpy as np

def hamming( a, b ):
	return np.mean([aa==bb for aa,bb in zip(a,b)])

def make_grishin( target_id, target_aln, template_id, template_aln, outfile ):
	with open(outfile, 'w') as out:
		out.write('## '+target_id+' '+template_id+'.pdb\n')
		out.write('#\n')
		out.write('score from program: 0\n')
		out.write('0 '+target_aln+'\n')
		out.write('0 '+template_aln+'\n')

def seq_to_aln_numbering( seqpos, seq, aln ):
	counter = -1
	for ii, aa in enumerate(aln):
		if aa != '-':
			counter += 1
		if counter == seqpos:
			return ii

def aln_to_seq_numbering( alnpos, seq, aln ) :
	num_gaps = aln[:alnpos].count('-')
	return alnpos - num_gaps

def define_disulfides( target_id, target_aln, template_aln_4hkx ):
	disulfides_4hkx = ((42,87), (7,19))

	target_disulfides = []
	template_seq_4hkx = template_aln_4hkx.replace('-','')
	target_seq = target_aln.replace('-','')
	for disulf in disulfides_4hkx:
		disulf_in_aln = (seq_to_aln_numbering(disulf[0], template_seq_4hkx, template_aln_4hkx),
						 seq_to_aln_numbering(disulf[1], template_seq_4hkx, template_aln_4hkx))
		target_disulfides.append( 
			(aln_to_seq_numbering( disulf_in_aln[0], target_seq, target_aln ),
			 aln_to_seq_numbering( disulf_in_aln[1], target_seq, target_aln )) 
			)

	with open( target_id+'.disulf', 'w') as out:
		for a,b in target_disulfides:
			## Adjust index to be one-based
			out.write(str(a+1)+' '+str(b+1)+'\n')

if __name__ == '__main__':
	input_file = sys.argv[1]
	
	base_name = input_file.split('.')[0]

	target_sequences = [(str(item.id),str(item.seq)) 
					for item in SeqIO.parse(input_file+'.fasta','fasta')]
	template_sequences = {str(item.id): str(item.seq)
					for item in SeqIO.parse('H1_templates.fasta','fasta')}

	for target_id, target_seq in target_sequences:

		similarities = []
		for template_id, template_seq in template_sequences.items():
			alignment = pairwise2.align.localxs(target_seq, template_seq,-2, -5)
			target_aln = alignment[0][0]
			template_aln = alignment[0][1]
			similarities.append( (template_id, target_aln, template_aln, hamming(target_aln, template_aln)) )

		similarities.sort(key=lambda a: a[-1], reverse=True)

		with open('template_identity.txt','w') as out:
			for template_id, target_aln, template_aln, similarity in similarities[:5]:
				out.write(template_id+' '+str(similarity)+'\n')
		for template_id, target_aln, template_aln, similarity in similarities[:5]:
			make_grishin( target_id, target_aln, template_id, template_aln, target_id+'_'+template_id+'.grishin')
			
		## Make disulfides based on 4hkx
		alignment = pairwise2.align.localxs(target_seq, template_sequences['4hkx'],-2, -5)
		target_aln = alignment[0][0]
		aln_4hkx = alignment[0][1]
		define_disulfides( target_id, target_aln, aln_4hkx )
