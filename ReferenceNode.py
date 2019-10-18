'''
Created on Jun 1, 2018

@author: vdeshpande
'''
import seq_util


class ReferenceNode(object):
    '''
    A node in the reference sequence graph
    Attributes:
        node_id: int NodeId
        seq: sequence of reference
        is_repeat: true if node is a self loop      
    '''


    def __init__(self, node_id, seq, is_repeat, node_name='', reference_region=""):
        '''
        Initialize ReferenceNode
        Arguments:
            node_id: id of node
            seq: sequence of the node
            is_repeat: True if node is a self-loop else False
        '''
        self.node_id = node_id
        self.seq = seq
        self.is_repeat = is_repeat
        self.node_name = node_name
        self.reference_region = reference_region


    def get_cigar_seq(self, cigar, offset=0):
        '''
        Get sequence of node clipped by offset and cigar length
        Args:
            cigar: cigar string
            offset: first matched position in the node
        
        '''
        end = offset
        for cigar_tuple in seq_util.parse_cigar(cigar):
            if cigar_tuple[1] in 'MDX=':
                end += cigar_tuple[0]
        return self.seq[offset:offset + end]

    
    def __str__(self):
        '''Returns sequence representation of node. If node is a self loop then adds parenthesis
        '''
        if self.is_repeat:
            return '(%s)' % self.seq
        else:
            return self.seq