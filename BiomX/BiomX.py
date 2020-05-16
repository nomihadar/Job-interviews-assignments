'''
BiomX Assignment
February 2020
Submitted by Nomi Hadar
'''


import numpy as np
import re
import wikipedia

def getCeleberiy(celeb_name):
    '''
    :param celeb_name: string of celebrity name
    :return: A link to the Wikipedia page of the input celebrity,
            The year in which the input celebrity was born

    '''

    titles = [celeb_name.title()]

    #get suggestions for spelling correction
    suggests = wikipedia.suggest(celeb_name)
    if suggests and type(suggests) == type([]):
        titles.extend(suggests)
    elif suggests:
        titles.append(suggests)

    #for each possible article title
    for title in titles:
        try:
            page = wikipedia.page(title)
            break
        except wikipedia.exceptions.PageError as e:
            continue
        except wikipedia.exceptions.DisambiguationError as e:
            pass
            print(e.options)
    else:
        #raise Exception("All items weren't successful")
        print("in PageError")
        url = 'Error: Page {} Not Found'.format(celeb_name)
        year = None
        return (url, year)

    #get URL
    url = page.url

    #get page content
    text = page.content

    #search for born year
    r1 = 'born [A-Za-z]+ [1-9]{1,2}, ([0-9]{2,4})'
    m = re.search(r1, text)
    year = m.group(1) if m else 'Not Found'

    return (url, year)

def simulate_paired_end_sequencing(genome, avg_chunk_len=400, avg_read_len=150, n_reads=100):
    '''
    function that simulates paired-end sequencing data of a given genome
    :param genome: string of A/C/G/T, genome is 10,000 bp long
    :param avg_chunk_len: Average length of DNA chunk "broken" out of the genome.
    :param avg_read_len: Average read length.
    :param n_reads: Total number of paired-end reads.
    :return: NGS (crude) data files (FASQ)
    '''

    # --------------------------------------------------------------------------
    # definitions
    # --------------------------------------------------------------------------
    OUTPUT_FILE = "BiomX_output_R{}.fastq"  # define output files
    PRECENT = 0.2 # percent, used to determine standard-deviation
    Q = lambda p: -10 * np.log10(p)  # Phred or Q score (not used)

    # --------------------------------------------------------------------------

    # draw DNA chunk lengths
    chunks_lens = np.random.normal(avg_chunk_len,
                                   avg_chunk_len * PRECENT,
                                   n_reads).round().astype(int)
    # draw read lengths
    reads_lens = np.random.normal(avg_read_len,
                                  avg_read_len * PRECENT,
                                  n_reads).round().astype(int)

    # draw indexes of "breaking points" of genome.
    # range of drawing: [0, length of genome - minimal length of read)
    indexes = np.random.randint(0, len(genome)-np.min(reads_lens), n_reads)

    #number of zeros to fill in seq id
    n_zfill = len(str(n_reads)) - 1

    read1 = {}
    read2 = {}
    for i in range(n_reads):

        #get chunk of DNA
        chunk = genome[indexes[i]:indexes[i]+int(chunks_lens[i])]

        #get left and right reads
        left_read = chunk[:reads_lens[i]]
        right_read = chunk[-reads_lens[i]:]

        # assume sequencing quality is perfect
        left_quality = "~" * len(left_read)
        right_quality = "~" * len(left_read)

        #create seq id
        seq_id = "@SEQ_{}, length:{},".format(str(i+1).zfill(n_zfill), len(left_read))

        #insert to dict
        read1[seq_id + " /1"] = (left_read, left_quality)
        read2[seq_id + " /2"] = (right_read, right_quality)

    #--------------------------------------------------------------------------
    # write FASTQ files
    # --------------------------------------------------------------------------
    def write_fastq(output_file, reads):
        with open(output_file, 'w') as f:
            for id, val in reads.items():
                read, quality = val
                f.write(id + "\n")
                f.write(read + "\n")
                f.write("+\n")
                f.write(quality + "\n")

    write_fastq(OUTPUT_FILE.format(1), read1)
    write_fastq(OUTPUT_FILE.format(2), read2)

if __name__ == "__main__":

    res = getCeleberiy("Barak Obama")
    print(res)

    simulate_paired_end_sequencing("A" * 10000, avg_chunk_len=400, avg_read_len=500, n_reads=100)
