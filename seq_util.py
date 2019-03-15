'''
Created on May 16, 2018

@author: vdeshpande
'''
import re
from collections import defaultdict


# Basepair properties
basepair_color = defaultdict(lambda: 'k', {'A': 'g', 'C': 'b', 'G': 'orange', 'T': 'r', 'N': 'k', 'R': 'k'})
complementary_nucleotide = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                            'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                            'N': 'N', 'n': 'n'}
def reverse_complement(seq):
    """Returns reverse complemented sequence of input string.
    """
    return ''.join([complementary_nucleotide[c] for c in seq[::-1]])


def parse_cigar(cigar):
    """Parses cigar string into list of tuples (len, type).
    """
    return [(int(t[0]), t[1]) for t in re.findall('(\d+)([A-Z]{1})', cigar)]


def get_cigar_insertions(cigar):
    """Parses cigar to get list of gap (insertion) tuples (0-based gap_position, gap_size).
    0-based gap_position is index of position of base-pair before which gap is inserted.
    e.g. if gap precedes any sequence, then gap_position = 0.
    TODO: Handles cigar flags S, H, P
    """
    cigar_tuples = parse_cigar(cigar)
    gap_list = []
    position = 0
    for t in cigar_tuples:
        if t[1] == 'I':
            gap_list.append((position, t[0]))
        elif t[1] in 'MDN=X':
            position += t[0]
    return gap_list


def get_cigar_query_align_len(cigar):
    """Parses cigar to get number of bases in query covered by cigar
    TODO: Handles cigar flags S, H, P
    """
    cigar_tuples = parse_cigar(cigar)
    num_query_bases = 0
    for t in cigar_tuples:
        if t[1] in 'MIN=XS':
            num_query_bases += t[0]
    return num_query_bases
            

def parse_graph_cigar(graph_cigar):
    """Parses graph cigar string into list of tuples (node_id, cigar).
    """
    return re.findall('([^\[]+)\[([^\]]*)\]', graph_cigar)


def get_node_cigar(read_seq, ref_seq):
    """ Compares read sequence and reference sequence position by position and returns cigar string
    Only supports cigar flags M, X, I and D
    Args:
        read_seq: read sequence
        ref_seq: reference sequence
    """
    match_status = ''
    match_count = 0
    cigar = ''
    for i in range(len(read_seq)):
        if read_seq[i] == ref_seq[i] and read_seq[i] == '-':
            continue
        elif read_seq[i] == '-':
            new_match_status = 'D'
        elif ref_seq[i] == '-':
            new_match_status = 'I'
        elif read_seq[i].lower() == ref_seq[i].lower():
            new_match_status = 'M'
        else:
            new_match_status = 'X'
        
        if new_match_status == match_status:
            match_count += 1
        else:
            if match_count > 0:
                cigar += str(match_count) + match_status
            match_count = 1
            match_status = new_match_status
    if match_count > 0:
        cigar += str(match_count) + match_status
    return cigar


def get_EH_graph_cigar(read_seq, ref_seq, repeat_unit, offset=0):
    """Obtain graph_cigar corresponding to EH v2.5 read alignment to reference sequence.
    Assumes that the alignment always intersects the repeat unit and indels are
    marked by '-' symbol.
    Args:
        read_seq: Read sequence from EH log file
        ref_seq: Reference sequence with repeat units replaced by 'R'
        repeat_unit: Sequence of the repeat unit
        offset: First position in the node to which the read is aligned
    Returns:
        graph_cigar string
    """
    graph_cigar = ''
    if len(ref_seq) > 0:
        left_flank_size = ref_seq.find('R')
        right_flank_size = len(ref_seq) - ref_seq.rfind('R') - 1
    if left_flank_size > 0:
        cigar = get_node_cigar(read_seq[:left_flank_size], ref_seq[:left_flank_size])
        graph_cigar += '0[%s]' % cigar
        offset = 0
    ref_node_seq = ''
    read_node_seq = ''
    for seq_index in range(left_flank_size, len(read_seq) - right_flank_size):
        read_node_seq += read_seq[seq_index] 
        if ref_seq[seq_index] != '-':
            ref_node_seq += repeat_unit[offset % len(repeat_unit)]
            offset += 1
        else:
            ref_node_seq += '-'
        if offset % len(repeat_unit) == 0:
            cigar = get_node_cigar(read_node_seq, ref_node_seq)
            graph_cigar += '1[%s]' % cigar
            read_node_seq = ''
    if right_flank_size > 0:
        cigar = get_node_cigar(read_seq[-1 * right_flank_size:], ref_seq[-1 * right_flank_size:])
        graph_cigar += '2[%s]' % cigar
    return graph_cigar
