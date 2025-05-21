import sys
import getopt
import re
from Bio import SeqIO
from check_internal_homology import check_internal_homology
import pickle

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

def mutate(proteome_input,miss,muts,gene_input=None):
	AA_dict = {
	'A':71.04,
	'R':156.10,
	'N':114.04,
	'D':115.03,
	'C':103.00,
	'E':129.04,
	'Q':128.06,
	'G':57.02,
	'H':137.06,
	'I':113.08,
	'L':113.08,
	'K':128.09,
	'M':131.04,
	'F':147.07,
	'P':97.05,
	'S':87.03,
	'T':101.04,
	'W':186.08,
	'Y':163.06,
	'V':99.07}

	print("Checking for Internal Sequence Homologies...")
	homologous_peptides = check_internal_homology(proteome_input,gene_input)


	handle=SeqIO.parse(proteome_input,'fasta')
	if muts >= 1 and muts <=3 and gene_input == False:
		split_proteome_input = re.split('[.]', proteome_input)
		output_file = 'mut_'+proteome_input
		homologous_file = f'homology_mut_{split_proteome_input[0]}.pkl'
		# Only outputs xle seqs for single sub mutants
		xle_seqs_out = f'xle_seqs_{split_proteome_input[0]}.pkl'
	elif muts >= 1 and muts <=3 and gene_input:
		split_proteome_input = re.split('[.]', proteome_input)
		output_file = 'ih_mut_'+proteome_input
		homologous_file = f'homology_ih_mut_{split_proteome_input[0]}.pkl'
		# Only outputs xle seqs for single sub mutants
		xle_seqs_out = f'xle_seqs_{split_proteome_input[0]}.pkl'
	# Trypsinzing without mutating doesn't currecntly check for homologous peptides
	elif muts == 0:
		output_file = 'tryp_'+proteome_input
		split_proteome_input = re.split('[.]', proteome_input)
		homologous_file = f'homology_ih_mut_{split_proteome_input[0]}.pkl'

	xle_seqs = dict()
	with open(output_file, 'w') as output:
		print("Trypsinizing and Mutagenizing...")
		if gene_input:
			for record in handle:
				proseq=str(record.seq)
				peptide_list=tryp(proseq,miss)
				count = 0
				for peptide in peptide_list:
					if peptide not in homologous_peptides:
						count2 = 1
						peptide_mass = 0.0
						for i in range(len(peptide)):
							try:
								peptide_mass = peptide_mass + AA_dict[peptide[i]]
							except KeyError:
								peptide_mass = 0
						if peptide_mass > 500 and len(peptide) >= 6:
							count += 1
							output.write(f">wild-{record.id}-{count}\n")
							output.write(f"{peptide}\n")
							#Make every single mutation in peptide
							if muts >= 1:
								for k in range(len(peptide)):
									for j in AA_dict:
										if j != peptide[k]:
											new_pep = peptide[:k] + j + peptide[(k+1):]
											#mut_count += 1
											dic_key = f"mutant-{record.id}-{count}.{count2}.{peptide[k]}{k+1}{j}"
											output.write(f">{dic_key}\n")
											output.write(f"{new_pep}\n")
											if (peptide[k] == 'L' and j == 'I') or (peptide[k] == 'I' and j == 'L'):
												xle_seqs[dic_key] = new_pep
											count2 += 1
											#Make every double mutation in peptide
											if muts >= 2:
												for kk in range(k+1,len(new_pep)):
													for jj in AA_dict:
														if jj != peptide[kk]:
															new_pep2 = new_pep[:(kk)] + jj + new_pep[(kk+1):]
															if j == new_pep[kk]:
																if peptide[k] != jj:
																	output.write(f">double-mutant-{record.id}-{count}.{count2}.{peptide[k]}{k+1}{j}.{new_pep[kk]}{kk+1}{jj}\n")
																	output.write(f"{new_pep2}\n")
																	count2 += 1
															if j != new_pep[kk]:
																output.write(f">double-mutant-{record.id}-{count}.{count2}.{peptide[k]}{k+1}{j}.{new_pep[kk]}{kk+1}{jj}\n")
																output.write(f"{new_pep2}\n")
																count2 += 1
															#Make every triple mutation in peptide
															if muts == 3:
																for kkk in range(k+2, len(new_pep2)):
																	for jjj in AA_dict:
																		if jjj != peptide[kkk]:
																			new_pep3 = new_pep2[:(kkk)] + jjj + new_pep2[(kkk+1):]
																			output.write(f">triple-mutant-{record.id}-{count}.{count2}.{peptide[k]}{k+1}{j}.{new_pep[kk]}{kk+1}{jj}.{new_pep2[kkk]}{kkk+1}{jjj}\n")
																			output.write(f"{new_pep3}\n")
																			count2 += 1

	if xle_seqs:
		with open(xle_seqs_out, 'wb') as xle:
			pickle.dump(xle_seqs, xle)
	if "mut" in output_file:
		with open(homologous_file, 'wb') as f:
			pickle.dump(homologous_peptides, f)
		return output_file

if __name__ == '__main__':
	try:
		options, remainder = getopt.getopt(sys.argv[1:],'', ['proteome_input=',
															'miss=',
															'muts=',
															'gene_input='])
	except getopt.GetoptError:
		#Gene input and therefore internal homology check is optional
		print("Example: python tryp_and_mutate.py --proteome_input custom_proteome.fasta --miss 2 --muts 1 --gene_input genes.fasta")
		sys.exit(2)

	for opt, arg in options:
		if opt == '--proteome_input':
			proteome_input=arg
		elif opt == '--miss':
			miss=int(arg)
		elif opt == '--muts':
			muts=int(arg)
		elif opt == '--gene_input':
			gene_input= arg
		else:
			print(f"Warning! Command-line argument: {opt} not recognized. Exiting..."); sys.exit(2)

	if muts > 3:
		print("Warning! Muts cannot exceed 3. Exiting..."); sys.exit(2)

	if remainder:
		print(f"Unrecognized arguments: {remainder}. These were ignored.")

	mutate(proteome_input,miss,muts,gene_input)