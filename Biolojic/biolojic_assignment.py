'''
Biolojic Assignment
Submitted by Nomi Hadar
March 2019
'''
import argparse
import re
import math
from random import randrange

# files names
CODONS_FILE = "codon_frequence_table.txt"
FASTA_FILE = "fasta_file.txt"
RECOGNITION_SITES_FILE = "recognition_sites.txt"
OUTPUT_NAME = "backtranslated_sequences.txt"

CODON_SIZE = 3 # codon size
NUM_CODONS = 64 # total number of codons
STOP_CODON = '*' #sign for stop codon

MAX_ITERATIONS = 30 #max iterations

# regular expressions
CODON = "([AUCG]{3})" #DNA codon
AMINO_ACID = "([A-Z\*])" # amino-acid. asterisk is stop codon
FLOAT = "(\d+\.\d+)" #float number
FASTA_SEQ_NAME = "(>.*)" # sequence name in fasta format
FASTA_SEQ = "([A-Z]*\*?)" # sequence in fasta format
RECOGNITION_SITE = '5\'-(.*)-3\'' # recognition site of restriction enzyme

# error message
ERROR_MSG = "Faile to remove recognition sites from sequence: {}"

###############################################################################

def read_codons_table(frequence_file):
    '''
    Read codons table into a dictionary
    :param file of codon usage frequence table
    :return: amino_to_codons dictionary in the format: {amino_acid: {codon: fraction}}
    and codon_to_amino dictionary in the format: {codon: amino_acid}
    '''

    # read file
    with open(frequence_file, 'r') as f:
        file_content = f.read()

    # parse frequencies table
    pattern = "{}\s+{}\s+{}".format(CODON, AMINO_ACID, FLOAT)
    table_entries = re.findall(pattern, file_content)

    #assert that the number of codons is 64
    assert len(table_entries) == NUM_CODONS

    # convert table to a dictionary
    amino_to_codons = {}
    codon_to_amino = {}
    for entry in table_entries:
        codon = entry[0].replace("U", "T") # convert U to T
        amino_acid = entry[1]
        fraction = float(entry[2])

        # add to amino_to_codons dict
        if amino_acid not in amino_to_codons:
            amino_to_codons[amino_acid] = {}
        amino_to_codons[amino_acid][codon] = fraction

        # add to codons_to_amino dict
        codon_to_amino[codon] = amino_acid

    return (amino_to_codons, codon_to_amino)


def create_random_choice_arrays(amino_to_codons):
    '''
    For each amino acid create an array of all its codons,
    where the proportion of each codon in the array
    is relative to its fraction in the codons table.
    The arrays will be used to do a random choice by frequency.
    :param amino_to_codons (dictionary)
    :return: dictionary in the format {amino_acid: [X, X, X, X, ..., Y, Y]}
    where X, Y are codons (array length is between 99 to 101)
    '''

    dic = {}
    #for each amino acid in codons table
    for amino_acid, codons_dic in amino_to_codons.items():
        # create an empty array
        dic[amino_acid] = []
        # for each codon of current amino acid
        for codon, fraction in codons_dic.items():
            # create an array in size fraction*100
            # where all elements are current codon
            array = [codon] * round(fraction*100)
            #extend array of current amino acid
            dic[amino_acid].extend(array)

    return dic


def read_fasta_file(fasta_file, add_stop_codon=True):
    '''
    :param fasta file of protein sequences
    :param flag whether or not to add a stop codon in the end of the sequence
    :return: dictionary in the format: {sequence_name: sequence}
    and list of sequences names (to keep order in output)
    '''

    # read fasta file
    with open(fasta_file, 'r') as f:
        file_content = f.read()

    # parse fasta file
    pattern = "{}\s+{}".format(FASTA_SEQ_NAME, FASTA_SEQ)
    rows = re.findall(pattern,file_content)

    # create dictionary
    seqs_dic = {}
    for seq_name, seq in rows:
        if add_stop_codon and seq[-1] != STOP_CODON:
                seq += STOP_CODON
        seqs_dic[seq_name] = seq

    #save order of
    seqs_order = [seq_name for seq_name, seq in rows]

    return (seqs_dic, seqs_order)


def read_recognition_sites(recognition_file):
    '''
    :param file of recognition sites of restriction enzymes
    :return: dictionary in the format: {restriction_enzyme: recognition_site}
    '''

    # read file
    with open(recognition_file, 'r') as f:
        file_content = f.read()

    # parse file
    pattern = "R[^-\s]*\s+(.*)\s+{}".format(RECOGNITION_SITE)
    rows = re.findall(pattern, file_content)

    # convert table to a dictionary
    dic = {}
    for res_enzyme, rec_site in rows:
        # remove backslash which symbolizes the cut position
        # replace N with all nucleotides
        # replace round brackets with square brackets
        rec_site = rec_site.replace("/", "") \
                            .replace("N", "[ACTG]") \
                            .replace("(", "[") \
                            .replace(")", "]")
        dic[res_enzyme] = rec_site

    return dic


def split_to_codons(dna_seq):
    '''
    :param DNA sequence (string)
    :return: sequence splitted into a list of codons ['TTA', 'ACG', ...]
    '''
    return [dna_seq[i:i+CODON_SIZE] for i in range(0, len(dna_seq), CODON_SIZE)]


def translate(dna_seq, codon_to_amino):
    '''
    convert a DNA sequence into a protein sequence
    :param DNA sequence (string)
    :param codon_to_amino (dic)
    :return: sequence converted to amino acids (string)
    '''

    protein_seq = ""
    # for each codon in the given DNA sequence
    for codon in split_to_codons(dna_seq):
        protein_seq += codon_to_amino[codon]
    return protein_seq


def backtranslate(prot_seq, codons_arrays, exclude_codon=''):
    '''
    convert a protein sequence into a DNA sequence
    :param protein sequence (string)
    :param codons arrays (dic)
    :param exclude codon (string)
    :return: sequence converted to DNA (string)
    '''

    dna_seq = ""
    # for each amino acid in the given protein sequence
    for amino_acid in prot_seq:
        # get array of codons corresponding to current amino acid
        codons_array = codons_arrays[amino_acid]

        #if codon specified - do not choose this codon
        if exclude_codon:
            codons_array = [codon for codon in codons_array
                            if codon != exclude_codon]

        # draw an index in range of the array size
        random_index = randrange(0, len(codons_array))
        # get the codon in that index
        codon = codons_array[random_index]
        # add codon to the DNA sequence
        dna_seq += codon

    return dna_seq


def contains_recognition_site(recognition_sites, dna_sequence):
    '''
    :param recognition_sites (dict)
    :param dna sequence (string)
    :return: True if sequence contains one of the recognition sites, False o.w.
    '''
    for res_enzyme, rec_site in recognition_sites.items():
        if re.search(rec_site, dna_sequence):
            return True
    return False


def replace_recognition_sites(recognition_sites, dna_seq,
                              codon_to_amino, codons_arrays):
    '''
    Replace all recognition sites found in a given DNA sequence
    by changing one of the codons contained in the recognition site
    :param recognition sites (dict)
    :param DNA sequence (string)
    :param codons_arrays (dict)
    :param codon_to_amino (dict)
    :return: a transformed DNA sequence
    '''

    # split sequence into a list of codons
    seq_by_codons = split_to_codons(dna_seq)

    # for each recognition site
    for res_enzyme, rec_site in recognition_sites.items():
        # for each match of recognition site in dna sequence
        for match in re.compile(rec_site).finditer(dna_seq):

            # draw an index in the range of match
            index = randrange(match.start(), match.end())

            # get start index of codon which contains the drawn codon
            codon_index = math.floor(index/CODON_SIZE)

            #get codon and its amino acid
            codon = seq_by_codons[codon_index]
            amino_acid = codon_to_amino[codon]

            # backtranslate codon to another codon
            new_codon = backtranslate(amino_acid, codons_arrays, codon)

            # replace codon by the new one
            seq_by_codons[codon_index] = new_codon

    new_seq = "".join(seq_by_codons)
    return new_seq


def write_output(seqs_dic, seqs_order, output_name):
    '''
    :param sequences (dict)
    :param list of sequences names in order of input file (list)
    :param output name
    '''
    output = "\n".join("{}\n{}".format(seq_name, seqs_dic[seq_name])
                       for seq_name in seqs_order)
    with open(output_name, 'w') as f:
        f.write(output)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-fasta_file', type=str,
                        required=False, default=FASTA_FILE,
                        help='fasta file of protein sequences')
    parser.add_argument('-codons_file', type=str,
                        required=False, default=CODONS_FILE,
                        help='file of codons table')
    parser.add_argument('-recognition_file', type=str,
                        required=False, default=RECOGNITION_SITES_FILE,
                        help='file with recognition sites of restriction enzymes')
    parser.add_argument('-output', type=str,
                        required=False, default=OUTPUT_NAME,
                        help='file output name of dna sequences')
    args = parser.parse_args()

    # read fasta file of proteins sequences
    (prot_seqs, seqs_order) = read_fasta_file(args.fasta_file)

    # read codon usage frequence table into dictionaries
    (amino_to_codons, codon_to_amino) = read_codons_table(args.codons_file)

    # create random choice arrays of codons per amino acid
    codons_arrays = create_random_choice_arrays(amino_to_codons)

    # convert each protein sequence into a DNA sequence
    dna_seqs = {}
    for seq_name, prot_seq in prot_seqs.items():
        dna_seq = backtranslate(prot_seq, codons_arrays)
        dna_seqs[seq_name] = dna_seq

    # read recognition sites of restriction enzymes
    rec_sites = read_recognition_sites(args.recognition_file)

    # for each back-translated sequence,
    # remove recognition sites of restriction enzymes
    seqs_no_rec_sites = {}
    for seq_name, dna_seq in dna_seqs.items():
        new_seq = dna_seq
        # limit number of attempts to replace recognition sites
        for iteration in range(MAX_ITERATIONS):
            # if sequence does not contain a recognition site - save sequence
            if (not contains_recognition_site(rec_sites, new_seq)):
                seqs_no_rec_sites[seq_name] = new_seq
                break
            # if sequence contains recognition sites - replace these sites
            new_seq = replace_recognition_sites(rec_sites, new_seq, codon_to_amino, codons_arrays)
        else:
            print (ERROR_MSG.format(seq_name))

    # verify that the DNA sequence still translates
    # into the original protein sequence
    for seq_name, dna_seq in seqs_no_rec_sites.items():
        original_protein_seq = prot_seqs[seq_name]
        translated_protein_seq = translate(dna_seq, codon_to_amino)
        assert translated_protein_seq == original_protein_seq

    #write output
    write_output(seqs_no_rec_sites, seqs_order, args.output)

    '''
    prints:
    
    print ("Amino acid to codons mapper:")
    print (amino_to_codons)

    print("\nCodon to amino acid mapper:")
    print(codon_to_amino)

    print ("\nProteins sequences:")
    print (prot_seqs)

    print("\nBacktranslated sequences (before replacing recognition sites):")
    print(dna_seqs)

    print("\nRecognition sites:")
    print(rec_sites)

    print("\nBacktranslated sequences (no recognition sites):")
    print(seqs_no_rec_sites)

    '''
