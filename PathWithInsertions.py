'''
Created on May 17, 2018

@author: vdeshpande
'''
from NodeInsertionList import NodeInsertionList
from seq_util import parse_graph_cigar, get_cigar_query_align_len

class PathWithInsertions(list):
    '''PathWithInsertions is a list of path nodes with a record
    of all insertions(gaps) in each node (NodeInsertionList).
    '''


    def __init__(self, node_id_list=[]):
        '''Constructor
        Aruguments:
            node_id_list: Ordered list of node id's in the path
        '''
        list.__init__(self, [NodeInsertionList(node_id) for node_id in node_id_list])
        
        
    def get_alignment_start_index(self, graph_cigar):
        """Find index of first node in alignment from graph_cigar string
        Arguments:
            graph_cigar: graph cigar string
        Returns:
            Ordered list of node indices in PathWithInsertions.
        """
        graph_cigar_tuples = parse_graph_cigar(graph_cigar)
        if len(graph_cigar_tuples) == 0:
            return 0
        first_id = int(graph_cigar_tuples[0][0])
        
        # left aligned if single node or in-repeat
        if first_id == int(graph_cigar_tuples[-1][0]):
            first_start_index = 0
            while first_start_index < len(self):
                current_node = self[first_start_index]
                if current_node.node_id == first_id:
                    return first_start_index
                first_start_index += 1
        cigar_index = 0
        while cigar_index < len(graph_cigar_tuples):
            if int(graph_cigar_tuples[cigar_index][0]) != int(graph_cigar_tuples[0][0]):
                break
            cigar_index += 1
        first_count = cigar_index
        first_end_index = 0
        current_node_id = self[first_end_index].node_id
        while first_end_index < len(self):
            if first_end_index == len(self) - 1:
                return first_end_index - first_count + 1
            next_node_id = self[first_end_index + 1].node_id
            if current_node_id == first_id and next_node_id != first_id:
                return first_end_index - first_count + 1
            first_end_index += 1
            current_node_id = next_node_id
        alignment_start_index = first_end_index - first_count + 1
        return alignment_start_index


    def update(self, graph_cigar, offset=0):
        """Update nodes with insertions in the read alignment.
        Arguments:
            graph_cigar: graph cigar string of read alignment.
            offset: 0-based position of first base in the first node in the alignment 
        """
        graph_cigar_tuples = parse_graph_cigar(graph_cigar)
        path_index = self.get_alignment_start_index(graph_cigar)
        cigar_tuple_index = 0
        for node_alignment in graph_cigar_tuples:
            while path_index < len(self) and self[path_index].node_id != int(node_alignment[0]):
                path_index += 1 
            if cigar_tuple_index == 0:
                self[path_index].update(node_alignment[1], offset)
            else:
                self[path_index].update(node_alignment[1])
            cigar_tuple_index += 1
            path_index += 1
            
            
    def get_alignment_coordinates(self, query_sequence, reference_graph, graph_cigar, offset=0, show_insertions=False):
        """Update nodes with insertions in the read alignment.
        Arguments:
            query_sequence: sequence of aligned read
            reference_graph: ReferenceGraph for locus
            graph_cigar: graph cigar string of read alignment.
            offset: 0-based position of first base in the first node in the alignment
            show_insertions: Boolean to show full sequence of insertions
        Returns:
            alignment_coordinates: a list of 3-tuples
                Each 3-tuple consists of (coordinate, query_bp, ref_bp)
                Gaps are represented by "-"
        """
        INTERNODE_GAP = 1
        alignment_coordinates = []
        graph_cigar_tuples = parse_graph_cigar(graph_cigar)
        path_index = self.get_alignment_start_index(graph_cigar)
        node_coordinate = 0
        cigar_tuple_index = 0
        query_pos = 0
        for i in range(path_index):
            node_coordinate += len(reference_graph[self[i].node_id].seq)
            if show_insertions:
                node_coordinate += self[i].get_total_insertion_size()
            node_coordinate += INTERNODE_GAP
        for node_alignment in graph_cigar_tuples:
            while path_index < len(self) and self[path_index].node_id != int(node_alignment[0]):
                node_coordinate += len(reference_graph[self[path_index].node_id].seq)
                if show_insertions:
                    node_coordinate += self[path_index].get_total_insertion_size()
                path_index += 1
                node_coordinate += INTERNODE_GAP
            cigar = node_alignment[1]
            #TODO: can we rely on graph cigar to always include bases from node extremes? Then we do not need offset
            if cigar_tuple_index == 0:
                node_offset = offset
            else:
                node_offset = 0
            node_align = self[path_index].get_alignment_coordinates(cigar,
                                                                query_pos,
                                                                query_sequence,
                                                                node_coordinate,
                                                                node_offset,
                                                                reference_graph,
                                                                show_insertions)                
            alignment_coordinates += node_align
            query_pos += get_cigar_query_align_len(cigar)
            cigar_tuple_index += 1
            node_coordinate += len(reference_graph[self[path_index].node_id].seq)
            if show_insertions:
                node_coordinate += self[path_index].get_total_insertion_size()
            node_coordinate += INTERNODE_GAP
            path_index += 1
            
        return alignment_coordinates


    def get_node_separators(self, reference_graph, show_insertions=False):
        """Return coordinate of node separators.
        Args:
            reference_graph: ReferenceGraph for locus
            show_insertions: Boolean to show full sequence of insertions
        Returns:
            alignment_coordinates: a list of 2-tuples
                Each 3-tuple consists of (coordinate, preceding node_id, following node_id)
                Gaps are represented by "-"
        """
        INTERNODE_GAP = 1
        node_coordinate = 0
        node_separators = []
        previous_node = reference_graph[self[0].node_id] if len(self) > 0 else None
        for i in range(len(self)):
            node_separators.append((node_coordinate - INTERNODE_GAP/2.0, reference_graph[previous_node.node_id], reference_graph[self[i].node_id])) 
            node_coordinate += len(reference_graph[self[i].node_id].seq)
            if show_insertions:
                node_coordinate += self[i].get_total_insertion_size()
            node_coordinate += INTERNODE_GAP
            previous_node = self[i]
        return node_separators
