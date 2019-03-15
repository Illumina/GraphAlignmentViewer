'''
Created on May 13, 2018

@author: vdeshpande
'''
from functools import total_ordering


@total_ordering
class ReadAlignment(object):
    """Alignment of a single read in the read_align file.
    Attributes:
        query_sequence: Sequence of aligned read
        reference_sequence: Sequence of reference to which read is aligned
        graph_cigar: Graph CIGAR string associated with alignment
        offset: Position within the first node of the first base-pair in the graph CIGAR
    Ordered by:
        read_seq.upper()
    """


    def __init__(self, query_sequence, reference_sequence, graph_cigar, offset, read_qual='', query_name=''):
        """Initializes read alignment class object.
        Args:
            query_sequence: Sequence of aligned read
            reference_sequence: Sequence of reference to which read is aligned
            graph_cigar: Graph CIGAR string associated with alignment
            offset: Position within the first node of the first base-pair in the graph CIGAR
            read_qual: Binary base call qualities (40: high, <40: low)
            query_name: Name of the read            
        """
        if read_qual == '':
            self.query_sequence = query_sequence
        else:
            self.query_sequence = ''
            for c in zip(query_sequence, read_qual):
                if c[1] == 40:
                    self.query_sequence += c[0]
                else:
                    self.query_sequence += c[0].lower()
        self.reference_sequence = reference_sequence
        self.graph_cigar = graph_cigar
        self.offset = offset
        self.query_name= query_name


    def __eq__(self, y):
        """Equality of two ReadAlignment objects based on equality of read_seq attribute
        """
        return self.query_sequence.upper() == y.query_sequence.upper()


    def __lt__(self, y):
        """Less than comparison of two read_align objects bases on comparison of read_seq attribute
        """
        return self.query_sequence.upper() < y.query_sequence.upper()