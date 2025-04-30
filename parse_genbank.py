import getopt
import sys
import re
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.Data import CodonTable

def translate_gene(gene, tt):
	aa_codons = CodonTable.unambiguous_dna_by_id[tt]
	stop_codons = aa_codons.stop_codons
	aa_codons_dict = dict(aa_codons.forward_table.items())
	for stop in stop_codons:
		aa_codons_dict[stop] = ''

	start_codons = aa_codons.start_codons

	gene_length = len(gene)
	translated_gene = ''
	codon = gene[0:3]

	if codon in start_codons:
		translated_gene += 'M'
	elif codon == 'GTG':
		translated_gene += 'M'
	elif codon == 'TTG':
		translated_gene += 'M'
	elif codon == 'CTG':
		translated_gene += 'M'
	else:
		translated_gene += aa_codons_dict[codon]

	for nt_codon in range(3, gene_length, 3):
		codon = gene[nt_codon:nt_codon + 3]
					
		if nt_codon + 3 < gene_length and codon in stop_codons:
			translated_gene += 'U'
		else:
			translated_gene += aa_codons_dict[codon]
	return translated_gene

def extract_ids(feat, parse, default_counter=[0], raw=False):
	qualifier_dic = {"lt": "locus_tag",
					 "g": "gene",
					 "p": "protein_id"}

	key = qualifier_dic.get(parse)
	if key and key in feat.qualifiers:
		locus_tag = str(feat.qualifiers[key][0])
		# Debug a record
		# if locus_tag == "FUN_003051":
		# 	print(f"Record Details for locus_tag {locus_tag}: {feat.qualifiers}")
		# 	print(f'Translation qualifier for locus tag {locus_tag}\n{feat.qualifiers["translation"][0]}')
		# 	sys.exit()
		if not raw:
			for ch in ['.', '-']:
				locus_tag = locus_tag.replace(ch, '')
		return locus_tag
	else:
		default_counter[0] += 1
		return f"no_locus_tag_{default_counter[0]}"


def parse_genbank(gb, tt, parse):
	recs = [rec for rec in SeqIO.parse(gb, "genbank")]

	split_gb = re.split('[.]', gb)
	file_name = split_gb[0]
	nt_output_file = f"{file_name}_genes.fasta"
	aa_output_file = f"{file_name}_proteome.fasta"
	mismatch_set = set()
	nt_dict = dict()
	aa_dict = dict()
	count = 0
	partial_codon_count = 0
	protein_count = 0
	no_locus_count = 0
	pseudo_count = 0
	cds_count = 0
	locus_tag_set = set()

	with open(nt_output_file, 'w') as ntf, open(aa_output_file, 'w') as aaf:
		for rec in recs:
			feats = [feat for feat in rec.features if feat.type == "CDS"]
			for feat in feats:
				cds_count += 1
				if "pseudo" in feat.qualifiers or "pseudogene" in feat.qualifiers:
					pseudo_count += 1
					continue

				location = str(feat.location)
				if "(+)" in location:
					location = "location=(+)"
				elif "(-)" in location:
					location = "location=(-)"

				locus_tag = extract_ids(feat, parse)
				if locus_tag in locus_tag_set and locus_tag is not None:
					print(f"{locus_tag} is a dulplicate ID; ensure parsing argument used corresponds to unique IDs; exiting")
					sys.exit(-1)
				else:
					locus_tag_set.add(locus_tag)

				try:
					nt_sequence = feat.extract(rec.seq)
				except ValueError as e:
					print(e)
					continue
				try:
					codon_start = int(feat.qualifiers["codon_start"][0])
					nt_sequence = nt_sequence[(codon_start-1):]
				except KeyError:
					pass

				partial_codon = len(nt_sequence) % 3

				if partial_codon != 0:
					partial_codon_count += 1
					#print(f"Warning: Sequence length for {locus_tag} is not a multiple of 3. Adjusting by removing {partial_codon} nucleotide(s).")
					nt_sequence = nt_sequence[:-partial_codon]


				try:
					gb_protein = feat.qualifiers["translation"][0]
				except KeyError:
					print("Error, missing translation; likely pseudogene; exiting; resolve and try again")
					print(feat)
					sys.exit(-1)

				if not re.match(r"^[ATGC]+$", str(nt_sequence)):
					print(f"Warning: Non-standard characters found in sequence for {locus_tag}; exiting")
					print(nt_sequence)
					sys.exit(1)

				try:
					# print(f'Translating locus_tag: {locus_tag}\n')
					# print(f'Translating nuceleotides: {nt_sequence}\n')
					translated_seq = translate_gene(nt_sequence, tt)
					# print(f'Translated sequence is: {translated_seq}\n')
				except TranslationError:
					if "transl_except" in feat.qualifiers:
						translated_seq = translate_gene(nt_sequence, tt)
						if translated_seq == gb_protein:
							print(f"{locus_tag} Translation exception successfully handled")
						else:
							translated_seq = translated_seq.replace('U', 'O')
							if translated_seq == gb_protein:
								print(f"{locus_tag} Translation exception successfully handled")
							else:
								print(f"{locus_tag} Translation exception failed; continuing")
								print(feat)
								continue

				# print(f"Processing locus_tag: {locus_tag}\n")
				# print(f"Nucleotide sequence: {nt_sequence}\n")
				# print(f"Translated sequence: {translated_seq}\n")
				# print(f"GB Protein: {gb_protein}\n")

				if translated_seq == gb_protein:
					if locus_tag:
						protein_count += 1
						ntf.write(f">{locus_tag} {location}\n")
						aaf.write(f">{locus_tag}\n")
						ntf.write(f"{nt_sequence}\n")
						aaf.write(f"{gb_protein}\n")
						nt_dict[locus_tag] = nt_sequence
						aa_dict[locus_tag] = gb_protein
					else:
						count += 1
						protein_count += 1
						locus_tag = f"no_locus_tag_{count}"
						no_locus_count += 1
						ntf.write(f">{locus_tag} {location}\n")
						aaf.write(f">{locus_tag}\n")
						ntf.write(f"{nt_sequence}\n")
						aaf.write(f"{translated_seq}\n")
						nt_dict[locus_tag] = nt_sequence
						aa_dict[locus_tag] = gb_protein
				else:
					print(f"Mismatch for locus_tag: {locus_tag}\n")
					print(f"Nucleotide sequence: {nt_sequence}\n")
					print(f"Translated sequence: {translated_seq}\n")
					print(f"GB Protein: {gb_protein}\n")
					if locus_tag:
						mismatch_set.add(locus_tag)
						nt_dict[locus_tag] = nt_sequence
						aa_dict[locus_tag] = gb_protein
					else:
						count += 1
						locus_tag = f"no_locus_tag_{count}"
						mismatch_set.add(locus_tag)
						nt_dict[locus_tag] = nt_sequence
						aa_dict[locus_tag] = gb_protein

		if mismatch_set:
			print(f"The following {len(mismatch_set)} of {cds_count} IDs have sequences that do not match between DNA and Protein:")
			# for item in mismatch_set:
			# 	print(item)

		if partial_codon_count != 0:
			print(f"{partial_codon_count} sequences were not divisble by 3")

		if no_locus_count != 0:
			print(f"{no_locus_count} genes did not have a locus tag")

		print(f"{protein_count} genes and proteins were parsed of {cds_count} CDSs")
		print(f"{len(set(aa_dict.values()))} unique sequences")
		print(f"{pseudo_count} pseudogenes were skipped")
		print("Next, please run a mass-spec search with the generated proteome output, then execute the following:")
		print(f"python generate_custom_database.py --gen_xlsx {file_name}_generic.xlsx --proteome_file {file_name}_proteome.fasta --gene_file {file_name}_genes.fasta")


def print_usage():
	print("Usage:")
	print("Example: python parse_genbank.py --gb species.gbff --tt 11 --parse lt")
	print("\nArguments:")
	print("--gb   :    GenBank file to be parsed (.gb, .gbk, or.gbff) [gbff seems to work best]")
	print("--tt   :    Translation table number (e.g., 11)")
	print("--parse:    Parse on either locus tags with option \"lt\" or gene names with option \"g\", or proteins IDs with option \"p\"")

if __name__ == '__main__':
	gb = None
	tt = None
	parse = None
	try:
		options, remainder = getopt.getopt(sys.argv[1:], 'h', ['gb=', 'tt=', 'parse=', 'help'])
	except getopt.GetoptError:
		print_usage()
		sys.exit(2)

	for opt, arg in options:
		if opt in ('-h', '--help'):
			print_usage()
			sys.exit(0)
		elif opt == '--gb':
			gb = arg
		elif opt == '--tt':
			tt = int(arg)
		elif opt == '--parse':
			parse = arg
		else:
			print(f"Warning! Command-line argument: {opt} not recognized. Exiting...")
			sys.exit(2)

	if remainder:
		print(f"Unrecognized arguments: {remainder}. These were ignored.")

	if not gb:
		print("Error: Missing --gb argument (GenBank file)")
		print_usage()
		sys.exit(2)
	if not tt:
		print("Error: Missing --tt argument (translation table number)")
		print_usage()
		sys.exit(2)
	if not parse:
		print("Error: Missing --tt argument (translation table number)")
		print_usage()
		sys.exit(2)

	parse_genbank(gb, tt, parse)
