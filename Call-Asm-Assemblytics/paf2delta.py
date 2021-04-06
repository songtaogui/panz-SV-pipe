#!/usr/bin/env python
# https://gist.githubusercontent.com/malonge/d347a0061eab77d9ebeb9181b7a9ff31/raw/44c832760cd43bf17eca37bfeae52cf7ec3299f4/paf2delta.py
import gzip
import argparse
from collections import defaultdict


class Alignment:

    def __init__(self, in_r_start, in_r_end, in_q_start, in_q_end, in_cigar, in_strand, in_num_matches, in_aln_len):
        self.r_start = int(in_r_start) + 1
        self.r_end = int(in_r_end)
        self.q_start = int(in_q_start) + 1
        self.q_end = int(in_q_end)
        self.cigar = in_cigar
        self.strand = in_strand
        self.num_matches = int(in_num_matches)
        self.aln_len = int(in_aln_len)

        self.parsed_cigar = []
        self._parse_cigar()

        if self.strand == "-":
            self.parsed_cigar = self.parsed_cigar[::-1]
            self.q_start, self.q_end = self.q_end, self.q_start

    def _parse_cigar(self):
        """ Given a CIGAR string, create a list of with each CIGAR operation as its own element. """
        cigar_chars = {
            'M',
            'I',
            'D',
            'N',
            'S',
            'H',
            'P',
            '=',
            'X'
        }

        this_field = ''
        for char in self.cigar:
            this_field += char
            if char in cigar_chars:
                self.parsed_cigar.append(this_field)
                this_field = ''


def write_delta(in_alns, in_r_lens, in_q_lens, in_file_name):
    with open(in_file_name, 'w') as f:
        f.write('file1 file2\n')
        f.write('NUCMER\n')

        # Iterate over each reference-query header pair
        for r_header in in_alns.keys():
            for q_header in in_alns[r_header].keys():
                f.write('>%s %s %r %r\n' % (r_header, q_header, in_r_lens[r_header], in_q_lens[q_header]))
                for z in in_alns[r_header][q_header]:
                    f.write('%r %r %r %r %r %r %r\n' % (
                        z.r_start,
                        z.r_end,
                        z.q_start,
                        z.q_end,
                        z.aln_len - z.num_matches,
                        z.aln_len - z.num_matches,
                        0
                    ))
                    # Continue with the cigar string
                    offsets = []
                    cigar = z.parsed_cigar
                    if cigar[0][-1] == 'S' or cigar[0][-1] == 'H':
                        cigar = cigar[1:-1]
                    else:
                        cigar = cigar[:-1]

                    counter = 1
                    for code in cigar:
                        if code[-1] == 'M':
                            counter += int(code[:-1])
                        elif code[-1] == 'D':
                            offsets.append(counter)
                            num_I = int(code[:-1])
                            for i in range(1, num_I):
                                offsets.append(1)
                            counter = 1
                        elif code[-1] == 'I':
                            offsets.append(-1*counter)
                            num_I = int(code[:-1])
                            for i in range(1, num_I):
                                offsets.append(-1)
                            counter = 1
                        else:
                            raise ValueError('Unexpected CIGAR code')
                    offsets.append(0)
                    offsets = [str(a) for a in offsets]
                    f.write('\n'.join(offsets) + '\n')


def main():
    parser = argparse.ArgumentParser(description="Convert a PAF file to a nucmer delta file.\nPAF file must be created with a CIGAR string '-c' Minimap2 parameter ")
    parser.add_argument("paf_file", metavar="<alns.paf>", type=str, help="PAF file to convert (gziped file allowed).")

    args = parser.parse_args()
    paf_file = args.paf_file
    alns = dict()

    # Dictionaries to store reference and query sequence lengths
    r_chr_lens = dict()
    q_chr_lens = dict()

    if paf_file[-3:] == ".gz":
        f = gzip.open(paf_file)
    else:
        f = open(paf_file, 'r')

    for line in f:
        if not isinstance(line, str):
            fields = line.decode("utf-8").split('\t')
        else:
            fields = line.split('\t')

        # Get the reference/query sequence lengths
        r_header = fields[5]
        q_header = fields[0]
        if r_header not in r_chr_lens:
            r_chr_lens[r_header] = int(fields[6])
        if q_header not in q_chr_lens:
            q_chr_lens[q_header] = int(fields[1])

        # Get the rest of the info and instantiate the Alignment object
        cigar_string = ''
        for i in fields[12:]:
            if i.startswith('cg:'):
                cigar_string = i.split(":")[2]
        if not cigar_string:
            raise ValueError("PAF file must contain a CIGAR string. Use 'minimap2 -c'")

        x = Alignment(
            fields[7],
            fields[8],
            fields[2],
            fields[3],
            cigar_string,
            fields[4],
            fields[9],
            fields[10]
        )

        # Add the alignments to the nested dictionary (first key=reference header, second key=query header)
        if r_header not in alns:
            alns[r_header] = defaultdict(list)
        alns[r_header][q_header].append(x)
    f.close()

    write_delta(alns, r_chr_lens, q_chr_lens, paf_file + '.delta')


if __name__ == "__main__":
    main()

