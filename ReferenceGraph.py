'''
Created on Jun 1, 2018

@author: vdeshpande
'''
import re
from ReferenceNode import ReferenceNode
import seq_util

class ReferenceGraph(list):
    '''
    Ordered list of ReferenceNode objects
    '''


    def __init__(self, repeat_spec, reference_fasta=None, flank_size=1000):
        '''
        Create reference graph from input specs
        Args:
           repeat_spec: Repeat specification for locus (dict inferred from json)
           reference_fasta: Reference FASTA file
           flank_size: Size of flanking sequences to append to locus. 
        '''
        list.__init__(self)
        position_list = []
        self.flank_size = flank_size
        if 'TargetRegion' in repeat_spec:
            if type(repeat_spec["TargetRegion"]) == str:
                self.target_region = repeat_spec["TargetRegion"]
                position_list = [repeat_spec["TargetRegion"]]
            else:
                self.target_region = (repeat_spec["TargetRegion"][0].split('-')[0]
                                      + '-' + repeat_spec["TargetRegion"][-1].split('-')[1])
                position_list = repeat_spec["TargetRegion"]
        elif 'ReferenceRegion' in repeat_spec:
            if type(repeat_spec["ReferenceRegion"]) == str:
                self.target_region = repeat_spec["ReferenceRegion"]
                position_list = [repeat_spec["ReferenceRegion"]]
            else:
                self.target_region = (repeat_spec["ReferenceRegion"][0].split('-')[0] 
                                      + '-' +repeat_spec["ReferenceRegion"][-1].split('-')[1])
                position_list = repeat_spec["ReferenceRegion"]
        else:
            self.target_region = (repeat_spec["ReferenceLoci"][0].split('-')[0]
                + '-' +repeat_spec["ReferenceLoci"][-1].split('-')[1])
            position_list = repeat_spec["ReferenceLoci"]
        if 'RepeatUnit' in repeat_spec:
            self.repeat_unit = repeat_spec['RepeatUnit']
            self.locus_structure = repeat_spec['RepeatUnit']
        elif 'LocusStructure' in repeat_spec:
            self.repeat_unit = re.sub(r'\*|\+', '', repeat_spec['LocusStructure'])
            self.locus_structure = repeat_spec['LocusStructure']
        else:
            self.repeat_unit = repeat_spec['RegionStructure']
            self.locus_structure = repeat_spec['RegionStructure']
        if 'RepeatId' in repeat_spec:
            self.repeat_id = repeat_spec['RepeatId']
            node_names = [self.repeat_id]
        elif 'LocusId' in repeat_spec:
            self.repeat_id = repeat_spec['LocusId']
            if type(repeat_spec["ReferenceRegion"]) == str:
                repeat_loci = [repeat_spec["ReferenceRegion"]]
            else:
                repeat_loci = repeat_spec["ReferenceRegion"]
            if 'VariantId' in repeat_spec:
                if type(repeat_spec['VariantId']) == str:
                    node_names = [repeat_spec["VariantId"]]
                else:
                    node_names = [vid for vid in repeat_spec["VariantId"]]
            else:
                if len(repeat_loci) > 1:
                    node_names = [repeat_spec["LocusId"] + '_' + rr for rr in repeat_loci]
                else:
                    node_names = [repeat_spec["LocusId"]]
        else:
            self.repeat_id = [s for s in repeat_spec['RepeatIds'] if 'ignore' not in s][0]
            node_names = [s for s in repeat_spec['RepeatIds']]
        self.chrom = self.target_region.split(':')[0]
        self.start = int(self.target_region.split(':')[1].split('-')[0])
        self.end = int(self.target_region.split(':')[1].split('-')[1])
        node_id = 0
        node_name_index = 0
        self.node_dict = {}
        self.repeat_unit = self.repeat_unit
        self.repeat_id = self.repeat_id
        if reference_fasta is None:
            self.append(ReferenceNode(node_id, 'N' * self.flank_size, False, "%s:%s-%s" % (self.chrom, self.start - self.flank_size, self.start)))
        else:
            self.append(ReferenceNode(node_id, reference_fasta.fetch(self.chrom, self.start - self.flank_size, self.start), False, "%s:%s-%s" % (self.chrom, self.start - flank_size, self.start)))
        node_id += 1
        if '(' not in self.repeat_unit:
            # Then entire string is a single loop
            self.append(ReferenceNode(1, self.repeat_unit, True, node_names[node_name_index], "%s:%s-%s" % (self.chrom, self.start, self.end)))
            node_id += 1
            node_name_index += 1
        else:
            repeat_unit_index = 0
            pos_i = 0
            while repeat_unit_index < len(self.repeat_unit):
                index2 = repeat_unit_index + 1
                while index2 < len(self.repeat_unit):
                    if self.repeat_unit[index2] in '()':
                        break
                    index2 += 1
                if self.repeat_unit[repeat_unit_index] == '(':
                    self.append(ReferenceNode(node_id, str(self.repeat_unit[repeat_unit_index + 1:index2]), True,
                                              node_names[node_name_index], position_list[pos_i]))
                    node_name_index += 1
                    pos_i += 1
                else:
                    start = self.start if (pos_i == 0) else position_list[pos_i - 1].split('-')[1]
                    end = self.end if (pos_i == len(position_list)) else position_list[pos_i].split('-')[0].split(':')[1]
                    self.append(ReferenceNode(node_id, str(self.repeat_unit[repeat_unit_index:index2]), False, "%s:%s-%s" % (self.chrom, start, end)))
                node_id += 1
                if index2 < len(self.repeat_unit) and self.repeat_unit[index2] == ')':
                    repeat_unit_index = index2 + 1
                else:
                    repeat_unit_index = index2
        if reference_fasta is None:
            self.append(ReferenceNode(node_id, 'N' * flank_size, False, "%s:%s-%s" % (self.chrom, self.end, flank_size)))
        else:
            self.append(ReferenceNode(node_id, reference_fasta.fetch(self.chrom, self.end,
                                                                     self.end + flank_size), False, "%s:%s-%s" % (self.chrom, self.end, flank_size)))
        


    def get_graph_cigar_seq(self, graph_cigar, offset=0):
        '''
        Get reference sequence corresponding to graph_cigar alignment
        Args:
            graph_cigar: graph_cigar string
            offset: first matched position in the first node (0-based)
        Returns:
            reference sequence
        '''
        graph_cigar_seq = ''
        graph_cigar_tuples = seq_util.parse_graph_cigar(graph_cigar)
        for graph_cigar_tuple in enumerate(graph_cigar_tuples):
            if graph_cigar_tuple[0] == 0:
                graph_cigar_seq += self[int(graph_cigar_tuple[1][0])].get_cigar_seq(graph_cigar_tuple[1][1], offset)
            else:
                graph_cigar_seq += self[int(graph_cigar_tuple[1][0])].get_cigar_seq(graph_cigar_tuple[1][1])
        return graph_cigar_seq
        

    def __str__(self):
        ''' Returns concatenation of str representation of all nodes, except flanking nodes
        '''
        return ''.join([str(node) for node in self[1:-1]])
