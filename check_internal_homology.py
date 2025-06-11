import sys
import getopt
from Bio import SeqIO
from collections import defaultdict
from itertools import combinations
import re

def tryp(proseq, miss_cleavage):
    peptides = []
    cut_sites = [0]

    for i in range(len(proseq) - 1):
        if (proseq[i] == 'K' or proseq[i] == 'R') and proseq[i + 1] != 'P':
            cut_sites.append(i + 1)

    if cut_sites[-1] != len(proseq):
        cut_sites.append(len(proseq))

    if len(cut_sites) > 2:
        for j in range(len(cut_sites) - 1):
            for k in range(j + 1, min(j + 2 + miss_cleavage, len(cut_sites))):
                peptides.append(proseq[cut_sites[j]:cut_sites[k]])
    else:
        peptides.append(proseq)
    return peptides

def count_mismatches(peptide1, peptide2):
	return sum(c1 != c2 for c1, c2 in zip(peptide1, peptide2))

def find_match(protein_seq, peptide):
	try:
		start_index = protein_seq.index(peptide)
		end_index = start_index + len(peptide)
		return start_index, end_index
	except ValueError:
		return -1, -1

def check_internal_homology(proteome_input, gene_input):
	AA_dict = {'A': 71.04, 'R': 156.10, 'N': 114.04, 'D': 115.03, 'C': 103.00, 'E': 129.04, 'Q': 128.06, 'G': 57.02, 'H': 137.06, 'I': 113.08, 
			   'L': 113.08, 'K': 128.09, 'M': 131.04, 'F': 147.07, 'P': 97.05, 'S': 87.03, 'T': 101.04, 'W': 186.08, 'Y': 163.06, 'V': 99.07}

	count = 0
	handle = SeqIO.parse(proteome_input,'fasta')
	seq_dic = defaultdict(dict)
	pro_dic = dict()
	for protein in handle:
		pro_id = str(protein.id)
		pro_seq = str(protein.seq)
		pro_dic[pro_id] = pro_seq
		peptide_list = tryp(pro_seq,2)
		count = 0
		for peptide in peptide_list:
			peptide_mass = 0.0
			for i in range(len(peptide)):
				try:
					peptide_mass = peptide_mass + AA_dict[peptide[i]]
				except KeyError:
					peptide_mass = 0
			if peptide_mass > 500 and len(peptide) >= 6:
				count += 1
				key = (f"{protein.id}-{count}")
				pep_len = len(peptide)
				seq_dic[pep_len][key] = peptide

	genes = SeqIO.parse(gene_input,'fasta')
	gene_dic = dict()
	for rec in genes:
		gen_seq = str(rec.seq)
		key = str(rec.id)
		gene_dic[key] = gen_seq

	peptides_for_removal = set()
	for length in seq_dic:
		for key1, key2 in combinations(seq_dic[length].keys(), 2):
			test_pep = seq_dic[length][key1]
			comp_pep = seq_dic[length][key2]
			mismatches = count_mismatches(test_pep,comp_pep)
			if mismatches == 1:
				if test_pep not in peptides_for_removal:
					peptides_for_removal.add(test_pep)
				if comp_pep not in peptides_for_removal:
					peptides_for_removal.add(comp_pep)

			elif mismatches == 0:
				test_key = re.split('[-]',key1)
				test_key = test_key[0]
				comp_key = re.split('[-]',key2)
				comp_key = comp_key[0]
				index_pep1 = find_match(pro_dic[test_key],test_pep)
				index_pep2 = find_match(pro_dic[comp_key],comp_pep)
				if index_pep1[0] > -1 and index_pep2[0] > -1:
					start_gene1 = (index_pep1[0]*3)
					end_gene1 = (index_pep1[1]*3)+2
					start_gene2 = (index_pep2[0]*3)
					end_gene2 = (index_pep2[1]*3)+2
					gene_seq1 = gene_dic[test_key][start_gene1:end_gene1]
					gene_seq2 = gene_dic[comp_key][start_gene2:end_gene2]
					if gene_seq1 != gene_seq2:
						peptides_for_removal.add(test_pep)

	return peptides_for_removal

if __name__ == '__main__':
	try:
		options, remainder = getopt.getopt(sys.argv[1:], '', ['proteome_input=','gene_input='])
	except getopt.GetoptError:
		print("Example: python check_internal_homology.py --proteome_input custom_proteome.fasta, --gene_input genes.fasta")
	
	for opt, arg in options:
		if opt == '--proteome_input': 
			proteome_input = arg
		elif opt == '--gene_input': 
			gene_input = arg
		else: 
			print(f"Warning! Command-line argument: {opt} not recognized. Exiting...")
			sys.exit(2)

	if remainder:
		print(f"Unrecognized arguments: {remainder}. These were ignored.")
	
	check_internal_homology(proteome_input, gene_input)
