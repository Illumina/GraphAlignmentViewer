'''
Created on May 13, 2018

@author: vdeshpande
'''
from seq_util import get_cigar_insertions, parse_cigar

class NodeInsertionList(list):
    '''
    List of insertions in a path node sorted by position.
    Each insertion is a tuple (position, size).
    Each position is unique and the max size of insertion is stored for each position.
    '''


    def __init__(self, node_id=None):
        '''
        Create empty list and set self.node_id
        '''
        list.__init__(self)
        self.node_id = node_id
    
    
    def update(self, cigar, offset=0):
        """ Update insertion_list using cigar string
        Args:
            cigar: linear cigar string
            offset: offset of first base of cigar with respect to node 
        """
        index = 0
        for gap in get_cigar_insertions(cigar):
            while index < len(self) and self[index][0] < gap[0] + offset:
                index += 1
            if index < len(self) and gap[0] + offset == self[index][0]:
                self[index] = (gap[0] + offset, max(gap[1], self[index][1]))
            else:
                self.insert(index, (gap[0] + offset, gap[1]))
                
                
    def get_total_insertion_size(self):
        """ Get sum of all insertions in node
        """
        return sum([i[1] for i in self])
        
        
    def get_alignment_coordinates(self, cigar, query_pos, query_seq, node_coordinate, node_offset, reference_graph, show_insertions=False):
        """Report list of path coordinates of each base-pair and gap in the node alignment.
        TODO: Handles cigar flags S, H, P
        Args:
            cigar: linear cigar string
            query_pos: index of first position in query sequence corresponding to the cigar
            query_seq: entire sequence of query
            node_coordinate: coordinate of the node in the path
            node_offset: index of first position in node sequence corresponding to the cigar
            reference_graph: ReferenceGraph for locus
            show_insertions: Boolean to show full sequence of insertions
        Returns:
            alignment_coordinates: a list of 3-tuples
                Each 3-tuple consists of (coordinate, query_seq, ref_seq)
                Gaps are represented by "-" 
        """
        alignment_coordinates = []
        cigar_tuples = parse_cigar(cigar)
        reference_sequence = reference_graph[self.node_id].seq
        query_index = query_pos
        reference_index = node_offset
        insertion_index = 0
        coordinate = node_coordinate + node_offset
        while insertion_index < len(self) and self[insertion_index][0] < node_offset:
            if show_insertions:
                coordinate += self[insertion_index][1]
            insertion_index += 1
        for t in cigar_tuples:
            if t[1] == 'S':
                query_index += t[0]
            elif t[1] == 'I':
                if insertion_index >= len(self) or self[insertion_index][0] != reference_index or self[insertion_index][1] < t[0]:
                    raise Exception("no reference insertion corresponding to query insertion")
                if show_insertions:
                    for ti in range(t[0]):
                        alignment_coordinates.append((coordinate,
                                                      query_seq[query_index],
                                                      '-'))
                        query_index += 1
                        coordinate += 1
                    for gap_i in range(t[0], self[insertion_index][1]):
                        alignment_coordinates.append((coordinate, '-', '-'))
                        coordinate += 1
                else:
                    alignment_coordinates.append((coordinate - 0.5,
                                                      r"$^{^{\mathbf{\bigvee}}}$",
                                                      '-'))
                    alignment_coordinates.append((coordinate - 0.5,
                                                      r"$^{^{^{\mathbf{%d}}}}$" % t[0],
                                                      '-'))
                    alignment_coordinates.append((coordinate - 0.5,
                                                      "|",
                                                      '-'))
                    query_index += t[0]
                insertion_index += 1
            elif t[1] == 'D':
                for ti in range(t[0]):
                    if insertion_index < len(self) and reference_index == self[insertion_index][0]:
                        if show_insertions:
                            for gap_i in range(self[insertion_index][1]):
                                alignment_coordinates.append((coordinate, '-', '-'))
                                coordinate += 1
                        insertion_index += 1
                    alignment_coordinates.append((coordinate,
                                                  '-',
                                                  reference_sequence[reference_index]))
                    reference_index += 1
                    coordinate += 1
            elif t[1] in 'MX=N':
                for ti in range(t[0]):
                    if insertion_index < len(self) and  reference_index == self[insertion_index][0]:
                        if show_insertions:
                            for gap_i in range(self[insertion_index][1]):
                                alignment_coordinates.append((coordinate, '-', '-'))
                                coordinate += 1
                        insertion_index += 1
                    if query_index >= len(query_seq):
                        raise Exception("query_index out of range", self.node_id, query_index, len(query_seq), cigar, query_pos, query_seq, node_coordinate, node_offset)
                    if reference_index >= len(reference_sequence):
                        print ("reference_index out of range", self.node_id, reference_index, len(reference_sequence), cigar, query_pos, query_seq, node_coordinate, node_offset)
                    alignment_coordinates.append((coordinate,
                                                  query_seq[query_index],
                                                  reference_sequence[reference_index]))
                    query_index += 1
                    reference_index += 1
                    coordinate += 1
        return alignment_coordinates
                            
            
