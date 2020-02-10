'''
EMENDO Assignment
February 2020
Submitted by Nomi Hadar
'''
import argparse
import pandas as pd
import numpy as np
import requests
import json

pd.set_option('display.max_columns', 20)
pd.options.display.max_colwidth = 100

#define human DNA hg38
GENOME = 'hg38'

# define PAM sequence (5' to 3')
PAM = 'NGG'

#define name of output file
OUTPUT_FILE = "output.fasta"

#define upstream & downstream lengths
UPSTREAM_LEN = 200
DOWNSTREAM_LEN = 200

#define minimal and maximal lengths of amplicons
MIN_LEN = 420
MAX_LEN = 423


def get_amplicon(row):
    '''
    Request sequence from UCSC server
    :param attributes: dictionary of attributes for request
    :return: DNA sequence
    '''

    # define server
    SERVER = "https://api.genome.ucsc.edu"

    # define request format
    REQUEST = "/getData/sequence?" \
              "genome={genome};" \
              "chrom={chrom};" \
              "start={start};" \
              "end={end};"

    #make a request with required attributes
    req = REQUEST.format(genome=GENOME,
                         chrom=row['chromosome'],
                         start=row['start'],
                         end=row['end'])

    #send request
    r = requests.get(SERVER + req)  # , headers={"Content-Type": "text/x-fasta"}

    #if request not o.k.
    if not r.ok:
        r.raise_for_status()
        return np.nan

    #convert request result to a dictionary
    d = json.loads(r.text)
    #exract sequence
    sequence = d['dna']

    return sequence

def create_amplicon_id(row):
    '''
    Create an amplicon id
    :param row: row of data frame
    :return:
    '''
    id = "{guide_name}_{length}bp_{position}|{genome}_{chrom}_{start}:{end}" \
        .format(genome=GENOME,
                guide_name=row['guide_name'],
                chrom=row['chromosome'],
                start=row['start'],
                end=row['end'],
                position=row['position'],
                length=row['new_length'])
    if row['length'] in ['ref', 'alt']:
        id += ", " + row['length']

    return id

def get_start_end_positions(row):
    '''
    Get start and end positions of desired sequence
    :param row: row in data frame
    :return: three new columns:
        'new_length': guide RNA length, not including PAM sequence.
        'start' and 'end': start and end positions of sequence, including additional upstream and downstream BPs.
    '''

    #if guide length is given
    if "bp" in row['length']:
        #get length of guide RNA
        new_len = int(row['length'].replace("bp", ""))
    else:  # assume PAM is NGG - length of 3
        new_len = len(row['sequence'][:-3])

    #compute start and end positions
    start = row['position'] - UPSTREAM_LEN
    end = row['position'] + new_len + DOWNSTREAM_LEN

    return pd.Series((new_len, start, end))

def write_fasta_file(df):
    '''
    Write output as a FASTA file
    :param df: data frame with aplicons and aplicons ids
    :return: none
    '''

    with open(OUTPUT_FILE, 'w') as f:
        for id, amplicon in zip(df['amplicon_id'], df['amplicon']):
            f.write(">" + id + "\n")
            f.write(amplicon + "\n")


def main(input_file):

    # read input file into a pandas data frame
    df = pd.read_csv(input_file)

    # replace names of two first columns
    df.rename(columns={df.columns[0]: 'guide_name', df.columns[1]: 'length'}, inplace=True)

    # fill missing values in by propagating last valid observation forward to next valid
    df.fillna(method='ffill', inplace=True)

    # remove space from chromosome column (e.g., 'chr 15' --> 'chr15')
    df['chromosome'] = df['chromosome'].str.replace(" ", "")

    # get start and end positions of desired sequence
    df[['new_length', 'start', 'end']] = df.apply(get_start_end_positions, axis=1)

    # create amplicons ids
    df['amplicon_id'] = df.apply(create_amplicon_id, axis=1)

    # get amplicons
    df['amplicon'] = df.apply(get_amplicon, axis=1)

    ############# validation #############

    # get amplicon length for validation
    df['amplicon_len'] = df['amplicon'].str.len()

    # assert amplicons lengths are correct
    assert(df['amplicon_len'].min() >= MIN_LEN)
    assert(df['amplicon_len'].max() <= MAX_LEN)

    # print amplicon ids of duplicae rows
    print("Duplicate guide RNAs:")
    print(df[df.duplicated()]['amplicon_id'])

    # remove duplicates rows, and pring their ids
    df.drop_duplicates(inplace=True)

    # if amplicons ids are NOT unique
    if df['amplicon_id'].nunique() != df['amplicon_id'].shape[0]:

        # print duplicate ids
        print("\nDuplicate amplicons IDs:")
        dup_ids = df[df['amplicon_id'].duplicated()]['amplicon_id']
        print(dup_ids)

        # enumerate duplicates
        dup = df['amplicon_id'].duplicated(keep=False)
        index = df.groupby('amplicon_id').cumcount() + 1
        df.loc[dup, 'amplicon_id'] = df['amplicon_id'].map(str) + ", " + index.astype(str)

    ############# write output #############

    # write output as a FASTA file
    write_fasta_file(df)


if __name__ == "__main__":

    #get path to input file
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help='input file')
    args = parser.parse_args()

    #call main function
    main(args.input_file)
