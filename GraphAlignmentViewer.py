#!/usr/bin/env python 
#
# STR_pileup
# Copyright(c) 2018 Illumina, Inc.
#
# Author: Viraj Deshpande < vdeshpande@illumina.com >
# Concept: Egor Dolzhenko < edolzhenko@illumina.com >, Michael Eberle < meberle@illumina.com >
#
# This program is free software: you can redistribute it and / or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see < http: // www.gnu.org / licenses / >.
#

import argparse
import json
import os
import bisect
from time import time
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pysam

matplotlib.mathtext.SHRINK_FACTOR = 0.8
matplotlib.mathtext.SCRIPT_SPACE = 0.4

from seq_util import basepair_color, reverse_complement, parse_graph_cigar,\
    get_EH_graph_cigar
from ReadAlignment import ReadAlignment
from ReferenceGraph import ReferenceGraph
from PathWithInsertions import PathWithInsertions

TSTART = time()


def get_EH_ReadAlignment(read, variant, repeat_unit, flank_size=1000):
    """Get ReadAlignment object from read alignment record from EHv2.5 json log files
    Arguments:
        read: YAML read alignment record produced by EH
        variant: Variant string in log file produced by EH
        repeat_unit: Sequence of repeat unit
        flank_size: Size of flank used by EH
    Returns:
        ReadAlignment corresponding to read alignment.
    """
    if 'INREPEAT' in variant:
        read_seq = read['irr'] if 'irr' in read else read['un_ir']
        ref_seq = read['anc'] if 'anc' in read else read['un_ma']
        if ref_seq.count(repeat_unit) < reverse_complement(ref_seq).count(repeat_unit):
            read_seq = reverse_complement(read_seq)
            ref_seq = reverse_complement(ref_seq)
        if ref_seq.upper().find(repeat_unit) != -1:
            offset = (len(repeat_unit) - ref_seq.upper().find(repeat_unit)) % len(repeat_unit)
        else:
            offset = (len(repeat_unit) - read_seq.upper().find(repeat_unit)) % len(repeat_unit)
        graph_cigar = get_EH_graph_cigar(read_seq, ''.join(['R' if b != '-' else '-' for b in ref_seq]), repeat_unit, offset)
    else:
        [read_seq, _, mapped_ref_seq] = read['align'].split('\n')[:3]
        if len(mapped_ref_seq) > 0 and mapped_ref_seq[0] == 'R':
            left_flank_size = 0
            offset = len(repeat_unit) - (mapped_ref_seq.rfind('R') % len(repeat_unit)) - 1
            repeat_seq_index = offset
        else:
            left_flank_size = mapped_ref_seq.find('R')            
            offset = flank_size - len([c for c in mapped_ref_seq[:left_flank_size] if c != '-'])
            repeat_seq_index = 0                        
        graph_cigar = get_EH_graph_cigar(read_seq, mapped_ref_seq, repeat_unit, offset)
        ref_seq = ''
        mapped_seq_index = left_flank_size
        while mapped_seq_index < len(mapped_ref_seq):
            if mapped_ref_seq[mapped_seq_index] == 'R':
                ref_seq += repeat_unit[repeat_seq_index]
                repeat_seq_index  = (repeat_seq_index + 1) % len(repeat_unit)
            else:
                ref_seq += mapped_ref_seq[mapped_seq_index]
            mapped_seq_index += 1
    return ReadAlignment(read_seq, ref_seq, graph_cigar, offset, query_name=read['name'])


def parse_read_align_EH(file_name, repeat_id_list=[], repeat_graphs={}, flank_size=1000):
    """Parses YAML-formatted read alignment file (log file) in EHv2.5 format into ReadAlignment objects.
    Args:
        file_name: Name of log file (include path)
        repeat_id_list (optional): List of repeat_ids to be parsed (Default=all)
    Returns:
        A dict mapping repeat IDs to lists of ReadAlignment objects.
    """
    import yaml
    y = yaml.load(open(file_name))
    read_aligns = defaultdict(lambda: [], {})
    for repeat_id in y:
        if len(repeat_id_list) > 0 and repeat_id not in repeat_id_list:
            continue
        site_reads = []
        repeat_unit = repeat_graphs[repeat_id].repeat_unit
        for variant in y[repeat_id]:
            for read in y[repeat_id][variant]:
                site_reads.append(get_EH_ReadAlignment(read, variant, repeat_unit, flank_size))
        read_aligns[repeat_id] = site_reads
    return read_aligns


def parse_read_align_graphEH(file_name, repeat_id_list=[]):
    """Parses YAML-formatted read alignment log file EHv2.7 format into ReadAlignment objects.
    Args:
        file_name: Name of log file (include path)
        repeat_id_list (optional): List of repeat_ids to be parsed (Default=all)
    Returns:
        A dict mapping repeat IDs to lists of ReadAlignment objects.
    """
    import yaml
    y = yaml.load(open(file_name))
    read_aligns = defaultdict(lambda: [], {})
    for repeat_id in y:
        if repeat_id_list != [] and repeat_id not in repeat_id_list:
            continue
        site_reads = []
        i = 0
        for read in y[repeat_id]:
            i += 1
            if 'alignment' not in read or read['alignment'] is None:
                continue
            read_seq = read['alignment'].split('\n')[2].replace(':', '').replace('-', '')
            ref_seq = read['alignment'].split('\n')[0].replace(':', '').replace('-', '')
            offset = int(read['path'].split('-')[0].strip('()').split('@')[1])
            site_reads.append(ReadAlignment(read_seq, ref_seq, read['graph_cigar'], offset=offset, query_name=read['name']))
        read_aligns[repeat_id] = site_reads
        print (file_name, repeat_id, len(site_reads))
    return read_aligns


def parse_realigned_bam_graphEH(file_name, repeat_id_list=[], repeat_graphs=None):
    """Parses bam read alignment file in EHv3 format into ReadAlignment objects.
    Args:
        file_name: Name of realigned bam file (include path)
        repeat_id_list (optional): List of repeat_ids to be parsed (Default=all)
    Returns:
        A dict mapping repeat IDs to lists of ReadAlignment objects.
    """
    y = pysam.AlignmentFile(file_name, 'rb')
    read_aligns = defaultdict(lambda: [], {})
    for read in y.fetch(until_eof=True):
        repeat_id = read.get_tag('XG').split(',')[0]
        if repeat_id not in repeat_graphs:
            continue
        if len(repeat_id_list) > 0 and repeat_id not in repeat_id_list:
            continue
        offset = int(read.get_tag('XG').split(',')[1])
        graph_cigar = read.get_tag('XG').split(',')[2]
        read_seq = read.query_sequence
        read_qual = read.query_qualities
        ref_seq = repeat_graphs[repeat_id].get_graph_cigar_seq(graph_cigar, offset)
        read_aligns[repeat_id].append(ReadAlignment(read_seq, ref_seq, graph_cigar, offset=offset, read_qual=read_qual, query_name=read.query_name))
    return read_aligns
        
            
def parse_read_align(file_name, repeat_id_list=[], file_format='v3', repeat_graphs={}, flank_size=1000):
    """Parses read alignment file for given EH version into ReadAlignment objects.
    Args:
        file_name: Name of log file (include path)
        repeat_id_list (optional): List of repeat_ids to be parsed (Default=all)
    Returns:
        A dict mapping repeat IDs to lists of ReadAlignment objects.
    """
    if file_format == 'v2.5':
        return parse_read_align_EH(file_name, repeat_id_list, repeat_graphs=repeat_graphs, flank_size=flank_size)
    elif file_format == 'v3.0.0-rc1':
        return parse_read_align_graphEH(file_name, repeat_id_list)
    elif file_format == 'v3':
        return parse_realigned_bam_graphEH(file_name, repeat_id_list, repeat_graphs=repeat_graphs)


def get_repeat_path(read_alignments, reference_graph):
    """Get a PathWithInsertions corresponding to highest count of each repeat_unit
    Arguments:
        read_alignments: dict mapping sample name to list of ReadAlignments
        reference_graph: ReferenceGraph for repeat locus
    """
    #TODO: Infer repeat_graph from read_alignments if repeat_specs not provided.
    max_node_count = defaultdict(lambda: 0, {})
    for sample in read_alignments:
        for read_alignment in read_alignments[sample]:
            node_count = defaultdict(lambda: 0, {})
            for node in parse_graph_cigar(read_alignment.graph_cigar):
                node_count[int(node[0])] += 1
            for node_id in node_count:
                max_node_count[node_id] = max(max_node_count[node_id], node_count[node_id])
    if max_node_count[reference_graph[0].node_id] == 0:
        max_node_count[reference_graph[0].node_id] = 1
    if max_node_count[reference_graph[-1].node_id] == 0:
        max_node_count[reference_graph[-1].node_id] = 1
    node_id_list = []
    for node in reference_graph:
        node_id_list += [node.node_id] * max_node_count[node.node_id]
    repeat_path = PathWithInsertions(node_id_list)
    print ("Reference path = %s" % [n.node_id for n in repeat_path])
    for sample in read_alignments:
        for read_alignment in read_alignments[sample]:
            repeat_path.update(read_alignment.graph_cigar, read_alignment.offset)
    return repeat_path


def sort_alignments_by_genotype(read_align_list, reference_graph, repeat_path, node_grouping_list=None, show_insertions=False):
    """Sort ReadAlignments by genotype length of prioritized nodes, spanning status and alignment position.
    Arguments:
        read_align_list: List of ReadAlignment to the locus in the sample
        reference_graph: ReferenceGraph for STR locus
        repeat_path: PathWithInsertions which to which alignments can be mapped.
        node_grouping_list: Ordered list of ReferenceNodes in decreasing order of priority. (Default: left to right, repeat units first)
        show_insertions: Boolean to show full sequence of insertions
    Returns:
        Sorted list of ReadAlignment
    """
    sorted_read_align_list = []
    spanning_genotypes = get_spanning_genotypes(read_align_list, reference_graph)
    for read_alignment in read_align_list:
        gt_labels = get_genotype_classification(read_alignment, reference_graph, spanning_genotypes, node_grouping_list)
        gt_status = []
        for gt in gt_labels:
            if gt[1] == '\leq':
                gt_status.append((gt[2], -1))
            elif gt[1] == '=':
                gt_status.append((gt[2], 0))
            else:
                gt_status.append((gt[2], 1))
        a = read_alignment
        coords = repeat_path.get_alignment_coordinates(a.query_sequence,
                                                    reference_graph,
                                                    a.graph_cigar, a.offset, show_insertions=show_insertions)[0][0]
        gt_status.append(coords)
        sorted_read_align_list.append((gt_status, read_alignment))
    sorted_read_align_list.sort(reverse=True)
    return [ra[1] for ra in sorted_read_align_list]


def get_spanning_genotypes(read_align_list, reference_graph):
    '''From a list of ReadAlignments, return a dict mapping each repetitive node_id to ordered list of genotypes spanned by a read.
    Arguments:
        read_align_list: List of ReadAlignment to the locus in the sample
        reference_graph: ReferenceGraph for STR locus
        node_grouping_list: Ordered list of ReferenceNodes in decreasing order of priority. (Default: left to right, repeat units first)
    Returns:
        dict mapping node_id to list of spanned genotypes
    '''
    spanning_genotypes = defaultdict(lambda: [], {})
    for read_alignment in read_align_list:
        node_count = np.zeros(len(reference_graph))
        left_flank_uncovered = np.ones(len(reference_graph))
        right_flank_uncovered = np.ones(len(reference_graph))
        previous_node_id = -1
        for node in parse_graph_cigar(read_alignment.graph_cigar):
            node_count[int(node[0])] += 1
            if previous_node_id != -1 and previous_node_id != int(node[0]):
                for i in range(previous_node_id, int(node[0])):
                    right_flank_uncovered[i] = 0
                    left_flank_uncovered[i + 1] = 0
            previous_node_id = int(node[0])
        for node in reference_graph:
            if node.is_repeat and right_flank_uncovered[node.node_id] == 0 and left_flank_uncovered[node.node_id] == 0:
                if node_count[node.node_id] not in spanning_genotypes[node.node_id]:
                    spanning_genotypes[node.node_id].append(node_count[node.node_id])
    for node_id in spanning_genotypes:
        spanning_genotypes[node_id].sort()
    return spanning_genotypes


def get_genotype_classification(read_alignment, reference_graph, spanning_genotypes, node_grouping_list=None):
    '''From a list of ReadAlignments, return a dict mapping each repetitive node_id to ordered list of genotypes spanned by a read.
    Arguments:
        read_align_list: List of ReadAlignment to the locus in the sample
        reference_graph: ReferenceGraph for STR locus
        node_grouping_list: Ordered list of ReferenceNodes in decreasing order of priority. (Default: left to right, repeat units first)
    Returns:
        list of 3 tuples for each node in node_grouping_list: (ReferenceNode, CompareStatus ("=", ">" or "\leq"), Genotype String)
    '''
    node_count = np.zeros(len(reference_graph))
    left_flank_uncovered = np.ones(len(reference_graph))
    right_flank_uncovered = np.ones(len(reference_graph))
    previous_node_id = -1
    for cigar_node in parse_graph_cigar(read_alignment.graph_cigar):
        node_count[int(cigar_node[0])] += 1
        if previous_node_id != -1 and previous_node_id != int(cigar_node[0]):
            for i in range(previous_node_id, int(cigar_node[0])):
                right_flank_uncovered[i] = 0
                left_flank_uncovered[i + 1] = 0
        previous_node_id = int(cigar_node[0])
    gt_class = []
    if node_grouping_list is None:
        node_grouping_list = ([n.node_id for n in reference_graph if (n.is_repeat and 'ignore' not in n.node_name)]
                                  + [n.node_id for n in reference_graph if (n.is_repeat and 'ignore' in n.node_name)]
                                  + [n.node_id for n in reference_graph if not n.is_repeat])
    for node_id in node_grouping_list:
        graph_node = reference_graph[node_id]
        if not graph_node.is_repeat:
            continue
        gt_compare_status = "="
        if (node_count[graph_node.node_id], left_flank_uncovered[graph_node.node_id], right_flank_uncovered[graph_node.node_id]) == (0, 1, 1):
            gt = -1
        elif ((left_flank_uncovered[graph_node.node_id] != right_flank_uncovered[graph_node.node_id]) or
              node_count[graph_node.node_id] > 0 and (left_flank_uncovered[graph_node.node_id], right_flank_uncovered[graph_node.node_id]) == (1, 1)):
            spanning_index = bisect.bisect_left(spanning_genotypes[graph_node.node_id], node_count[graph_node.node_id])
            if spanning_index >= len(spanning_genotypes[graph_node.node_id]):
                gt_compare_status = '>'
                if len(spanning_genotypes[graph_node.node_id]) == 0:
                    gt = 0
                else:
                    gt = spanning_genotypes[graph_node.node_id][-1]
            else:
                gt = spanning_genotypes[graph_node.node_id][spanning_index]
                gt_compare_status = '\leq'
        else:
            gt = node_count[graph_node.node_id]
        gt_string = r'$gt%s %s %s$' % (graph_node.node_id, gt_compare_status, int(gt))
        if gt > 0 and gt_compare_status == '\leq':
            spanning_index = bisect.bisect_left(spanning_genotypes[graph_node.node_id], node_count[graph_node.node_id])
            if spanning_index > 0:
                gt_string = r"$%s < gt%s \leq %s$" % (int(spanning_genotypes[graph_node.node_id][spanning_index - 1]), graph_node.node_id, int(gt))
        gt_class.append((graph_node, gt_compare_status, gt, gt_string))
    return gt_class


def get_pileup(read_align_list, reference_graph, repeat_path=None, node_grouping_list=None, show_insertions=False):
    """ Get genotype_pileup_list for given sample
    Args:
        read_align_list: List of read_alignments to repeat_id (see parse_read_align)
        reference_graph: ReferenceGraph for STR locus
        repeat_path: PathWithInsertions for read alignments (inferred if not provided)
        node_grouping_list: Ordered list of ReferenceNodes in decreasing order of priority. (Default: left to right, repeat units first)
        show_insertions: Boolean to show full sequence of insertions       
    Returns:
        Pileup: (PathWithInsertions, List of 3-tuples for each alignment (ReadAlignment, AlignmentCoordinates, GenotypeString))
        where:
            AlignmentCoordinates is a list of (Position, QueryBasePair, RefBasePair)
    """
    if repeat_path is None:
        repeat_path = get_repeat_path({None:read_align_list}, reference_graph)

    ra = sort_alignments_by_genotype(read_align_list, reference_graph, repeat_path, node_grouping_list, show_insertions=show_insertions)
    coords = [repeat_path.get_alignment_coordinates(a.query_sequence,
                                                    reference_graph,
                                                    a.graph_cigar, a.offset, show_insertions=show_insertions)
                for a in ra]
    spanning_genotypes = get_spanning_genotypes(ra, reference_graph)
    genotype_labels = [get_genotype_classification(a, reference_graph, spanning_genotypes, node_grouping_list) for a in ra]
    return (repeat_path, [z for z in zip(ra, coords, genotype_labels)])


def get_pileup_list(all_sample_read_aligns, sample_names, reference_graph, node_grouping_list=None, show_insertions=False):
    """ Generate genotype_pileups for all sample
    Args:
        all_sample_read_aligns: dict: {sample_name : {repeat_id: [list of read aligns]}}
        sample_names: ordered list of sample names in pileup
        reference_graph: ReferenceGraph for STR locus
        node_grouping_list: Ordered list of ReferenceNodes in decreasing order of priority. (Default: left to right, repeat units first)
        show_insertions: Boolean to show full sequence of insertions   
    Returns:
        List of 3-tuples with pileups for each sample
            e.g. [(sample_name, PathWithInsertions, Pileup (see get_pileup))]
    """
    for s in all_sample_read_aligns:
        repeat_aligns = [a for a in all_sample_read_aligns[s] if len([n for n in parse_graph_cigar(a.graph_cigar) if reference_graph[int(n[0])].is_repeat]) > 0]
        all_sample_read_aligns[s] = repeat_aligns
    repeat_path = get_repeat_path(all_sample_read_aligns, reference_graph)

    pileup_list = []
    for sample_name in sample_names:
        if sample_name in all_sample_read_aligns:
            genotype_pileup_list = get_pileup(all_sample_read_aligns[sample_name], reference_graph, repeat_path, node_grouping_list, show_insertions=show_insertions)
        else:
            genotype_pileup_list = get_pileup([], reference_graph, repeat_path, node_grouping_list, show_insertions=show_insertions)
        pileup_list.append((sample_name, repeat_path, genotype_pileup_list[1]))
    return pileup_list


def plot_pileup(repeat_id, pileup_list, genotypes, reference_graph, greyscale=False, output_name='', title_prefix='', show_read_names=False, show_insertions=False, flank_clip_size=-1, dpi=100, pdf=False):
    """Plot read pileups for all samples for give repeat_id
    Args:
        repeat_id: Repeat ID for which to plot pileup
        pileup_list: See get_pile_list()
        genotypes: Genotypes reported by EH for each sample
        reference_graph: ReferenceGraph for locus
        greyscale: (optional): Set True for greyscale scheme. Default: color(IGV scheme).
        output_prefix: (optional) Default ''
        title_prefix: String to be added as a prefix to the title
        show_read_names: Show names of reads alongside the read
        show_insertions: Boolean to show full sequence of insertions
    Outputs:
        PNG and PDF files: output_prefix + "_" + repeat_id
    """
    repeat_path = pileup_list[0][1]
    node_separators = repeat_path.get_node_separators(reference_graph, show_insertions=show_insertions)
    ref_path_seq = ''.join([reference_graph[node.node_id].seq for node in repeat_path])
    ref_graph_cigar = ''.join(['%s[%sM]' % (str(node.node_id), len(reference_graph[node.node_id].seq)) for node in repeat_path])
    ref_coords = repeat_path.get_alignment_coordinates(ref_path_seq, reference_graph, ref_graph_cigar, show_insertions=show_insertions)

    if flank_clip_size == -1:
        xoffset = min([position[0] for sample in pileup_list for align in sample[2] for position in align[1]])
        xmax = max([position[0] for sample in pileup_list for align in sample[2] for position in align[1]])
    else:
        xoffset = repeat_path.get_alignment_coordinates('N', reference_graph, '%s[1M]' % reference_graph[0].node_id, offset=len(reference_graph[0].seq) - flank_clip_size, show_insertions=show_insertions)[0][0]
        xmax = repeat_path.get_alignment_coordinates('N', reference_graph, '%s[1M]' % reference_graph[-1].node_id, offset=flank_clip_size, show_insertions=show_insertions)[0][0]
    if show_read_names:
        read_name_max = max([len(alignment_pileup[0].query_name) for site_pileup in pileup_list for alignment_pileup in site_pileup[2]] + [len("Reference")])
        xlen = xmax + read_name_max - xoffset
    else:
        xlen = xmax - xoffset
    ylen = 0 
    for site_pileup in pileup_list:
        ylen += 3.5 + 3.5 * len(set([str(align[2]) for align in site_pileup[2]])) + len(site_pileup[2])
    
    
    xscale = 1.0 / 10
    yscale = 2.0 / 10
    fontsize = 8
    margin = 0.5
    fig = plt.figure(figsize=(xlen * xscale + 2 *
                                margin,  ylen * yscale + 2 * margin))
    if title_prefix == "":
        fig.suptitle("%s:%s" % (repeat_id, reference_graph.locus_structure), fontsize=4 * fontsize, y=1 - margin / (ylen * yscale + 2 * margin), va='bottom')
    else:
        fig.suptitle("%s %s:%s" % (title_prefix, repeat_id, reference_graph.locus_structure), fontsize=4 * fontsize, y=1 - margin / (ylen * yscale + 2 * margin), va='bottom')
    ax = fig.add_axes([margin / (xlen * xscale + 1), margin / (ylen * yscale + 1),
                       xlen * xscale / (xlen * xscale + 2 * margin), ylen * yscale / (ylen * yscale + 2 * margin)])
    ax.set_xlim(xoffset, xoffset + xlen)
    ax.set_ylim(ylen, 0)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ypos = 0
    ax.axhline(ypos + 1, linestyle='--', c='0.5')

    for site_pileup in pileup_list:
        gt_list = []
        gtnone = True
        for node_id in range(len(reference_graph)):
            if not reference_graph[node_id].is_repeat:
                continue
            if site_pileup[0] in genotypes and genotypes[site_pileup[0]] is None:
                continue
            if site_pileup[0] in genotypes and reference_graph[node_id].node_name in genotypes[site_pileup[0]]:
                node_name = reference_graph[node_id].node_name
                gt = '/'.join([str(s) for s in genotypes[site_pileup[0]][node_name]])
            else:
                gt = '??'
            gtnone = False
            gt_list.append("gt%d:%s:%s" % (node_id, reference_graph[node_id].seq, gt))
        if gtnone:
            sample_str = site_pileup[0]
        else:
            sample_str = site_pileup[0] + ' (Prediction - ' + '; '.join(gt_list) + ')'
        ax.text((xoffset + xoffset + xlen) / 2, ypos + 3, sample_str, ha='center',
                va='bottom', fontsize=3 * fontsize)
        ypos += 3
        sample_start = ypos
        previous_genotype = None
        for alignment_pileup in site_pileup[2]:
            if alignment_pileup[2] != previous_genotype:
                ax.axhline(ypos + 0.5, linestyle=':', c='0.25')
                left_coord = xoffset
                for node_separator in node_separators:
                    if node_separator[1].node_id == node_separator[2].node_id:
                        continue
                    for gt in alignment_pileup[2]:
                        if gt[0].node_id == node_separator[1].node_id:
                            ax.text((left_coord + node_separator[0]) / 2, ypos + 2, gt[3], ha='center', va='bottom', fontsize=2 * fontsize)
                            break
                    left_coord = node_separator[0]
                for gt in alignment_pileup[2]:
                    if gt[0].node_id == repeat_path[-1].node_id:
                        ax.text((left_coord + xmax) / 2, ypos + 2, gt[3], ha='center', va='bottom', fontsize=2 * fontsize)
                        break
                ypos += 2.5
                for position in ref_coords:
                    if position[0] < xoffset or position[0] >= xmax:
                        continue
                    if greyscale:
                        ax.text(position[0] + 0.5, ypos + 1, position[1], ha='center',
                                va='center', fontsize=fontsize, weight='bold')
                    else:
                        c = basepair_color[position[1].upper()]
                        ax.text(position[0]+ 0.5, ypos + 1, position[1], ha='center',
                                va='center', fontsize=fontsize, weight='bold', color=c)
                if show_read_names:
                    ax.text(xmax + 2, ypos + 1, "Reference", color='k', ha='left', va='center')
                previous_genotype = alignment_pileup[2]
                
                
                for node_separator in node_separators:
                    if node_separator[1].node_id != node_separator[2].node_id:
                        ax.plot((node_separator[0], node_separator[0]), (ypos - 2, ypos + 1.5), color="k", linewidth=0.5)
                
                current_node_count = 1
                ax.axhline(ypos + 1.5, zorder=2)
                ax.axhline(ypos + 0.5, zorder=2)
                for node_separator in node_separators:
                    if node_separator[1].is_repeat:
                        if len(node_separator[1].seq) > 3 or current_node_count % 5 == 0:
                            ax.plot((node_separator[0], node_separator[0]), (ypos, ypos + 1.5), color='k')
                            ax.text(node_separator[0] - 0.5, ypos, str(current_node_count), ha='right', va='center', fontsize=2 * fontsize)
                    if node_separator[1].node_id == node_separator[2].node_id:
                        current_node_count += 1
                    else:
                        current_node_count = 1
                read_parity = 0
                read_count = 0
                ypos += 1
            ppos = -1
            if read_parity == 0:
                ax.axhline(ypos + 0.95, color="0.9", linewidth=15)
                read_color = 'magenta'
            else:
                read_color = 'cyan'
            if read_count == 0:
                ax.axhline(ypos + 0.5, zorder=2)
                read_count += 1
                
            for position in alignment_pileup[1]:
                if position[0] < xoffset or position[0] >= xmax:
                    continue
                if ppos != -1 and ppos + 1 < position[0]:
                    ax.plot((ppos + 1, position[0]), (ypos + 1, ypos + 1), color=read_color)
                if greyscale:
                    if position[1] != position[2]:
                        c = 'r'
                    elif position[1].islower():
                        c = '0.5'
                    else:
                        c = 'k'
                else:
                    c = basepair_color[position[1].upper()]
                ax.text(position[0] + 0.5, ypos + 1, position[1], ha='center',
                        va='center', fontsize=fontsize, color=c, weight='bold')
                ppos = position[0]
            for node_separator in node_separators:
                if node_separator[1].node_id != node_separator[2].node_id:
                    ax.plot((node_separator[0], node_separator[0]), (ypos + 0.5, ypos + 1.5), color="k", linewidth=0.5)
            if show_read_names:
                ax.text(xmax + 2, ypos + 1, alignment_pileup[0].query_name, color='k', ha='left', va='center')
            read_parity = 1 - read_parity
            ypos += 1
        ax.axhline(ypos + 0.5, linestyle='--', c='0.5')
        if show_read_names:
            ax.plot((xmax + 1, xmax + 1), (sample_start + 0.5, ypos + 0.5), color='k')
            # ax.axvline(xmax + 1, color='k')
    if pdf: 
        fig.savefig(output_name + '.pdf')
    else:
        fig.savefig(output_name + '.png', dpi=dpi)
    plt.close(fig)


def get_args():
    """ Get arguments
    """
    parser = argparse.ArgumentParser(
        description="Create a pileup of read alignments from ExpansionHunter")
    parser.add_argument("--variant_catalog", dest='variant_catalog', type=str,
                        help="Path to variant catalog json file used to run EH (default: \"variant_catalog.json\")", required=True)
    read_align_group = parser.add_mutually_exclusive_group(required=True)
    read_align_group.add_argument(
        "--read_align", dest='read_align_file', type=str, help="Read alignment log file from EH output")
    read_align_group.add_argument(
        "--read_align_list", dest='read_align_file_list', type=str, help="CSV file with paths to EH output for multiple samples. Column1: Sample name, Column2: Read alignment file from EH, Column 3(optional) VCF file")
    parser.add_argument("--gt_file", dest='gt_file', type=str,
                                help=" (Optional) Output file (VCF or JSON) from EH to display GT predicted in output image")
    parser.add_argument("--file_format", dest='file_format', type=str,
                        help="Format of read alignments from EH. [\"v3\": BAM(default), \"v2.5\": YAML]", default="v3")
    parser.add_argument("--greyscale", action="store_true",
                        help="Show nucleotides in greyscale: high quality match - black, low quality match - grey, mismatch - red (Default IGV color scheme)")
    parser.add_argument("--locus_id", dest='locus_id',
                        type=str, help="Comma-separated locus IDs for which to plot pileup. Default: All repeats", default="")
    parser.add_argument("--output_prefix", dest='output_prefix',
                        type=str, help="Prefix of output file. Output filename(s) \"<OUTPUT_PREFIX>_<CHROM>-<START>-<REPEATUNIT>.alignment.png\"(\".pdf\") corresponding to the position of the first repeat unit in the node grouping. If node grouping is \"NONE\" or \"ALL\", then position corresponds to the first repeat unit in the locus. (default: No prefix. Output filename(s) \"<CHROM>-<START>-<REPEATUNIT>.alignment.png\"(\".pdf\"))", default="")
    parser.add_argument("--output_dir", dest='output_dir',
                        type=str, help="Output directory (default=Current working directory)", default="")
    parser.add_argument("--title_prefix", dest='title_prefix',
                        type=str, help="Prefix text to be appended to title of the plot (default:\"\")", default="")
    parser.add_argument("--reference_fasta", dest='reference_fasta',
                        type=str, help="Fasta file with .fai index for reference sequence. If not provided, flanks are set to 'N'. Default: None", default=None)
    parser.add_argument("--node_grouping", dest='node_grouping',
                        type=str, help="Comma-separated list of node indices (left flank=0) to group and sort reads by genotype. \"NONE\": sort reads only by position, \"ALL\": group by all repeat nodes from left to right. (default: create a separate image for each repeat unit)", default="")
    parser.add_argument("--show_read_names", action="store_true",
                        help="Show names of reads")
    parser.add_argument("--show_insertions", action="store_true",
                        help="Show full sequence of insertions")
    parser.add_argument("--region_extension_length", dest="flank_size", type=int,
                        help="Size of nodes flanking the locus used for generating the read alignments. Default: 1000", default=1000)
    parser.add_argument("--region_extension_clip_length", dest="flank_clip_size", type=int,
                        help="Number of basepairs of flanking regions to display. `-1`: Infer from maximum span of reads overlapping the locus. (default 20)", default=20)
    graphics_format_group = parser.add_mutually_exclusive_group()
    graphics_format_group.add_argument("--dpi", dest='dpi', type=int,
                                help=" (Optional) Resolution of output PNG image. Default: 100", default=100)
    graphics_format_group.add_argument("--pdf", dest='pdf', action='store_true',
                                help=" (Optional) Output PDF image instead of PNG")
    
    return parser.parse_args()


def get_sample_filepaths(args):
    """Parse arguments to get sample paths and output prefix
    Args:
        Argparse arguments
    Returns:
        List of 3-tuples for each samples: (sample_name, read_align_file, vcf_file)
    """
    sample_list = []
    # output_prefix: prefix of final output images
    # Populate sample_list and output_prefix
    if args.read_align_file is not None:
        sample_name = os.path.splitext(
            os.path.basename(args.read_align_file))[0]
        if args.gt_file is not None:
            sample_list = [(sample_name, args.read_align_file, args.gt_file)]
        else:
            sample_list = [(sample_name, args.read_align_file, None)]
        # output_prefix = sample_name
    else:
        in_root = os.path.dirname(args.read_align_file_list)
        for l in open(args.read_align_file_list):
            if len(l.strip()) == 0 or l[0] == '#':
                continue
            ll = l.strip().split(',')
            if len(ll) >= 3:
                sample_list.append((ll[0], os.path.join(
                    in_root, ll[1]), os.path.join(in_root, ll[2])))
            else:
                sample_list.append((ll[0], os.path.join(
                    in_root, ll[1]), None))
        # output_prefix = os.path.splitext(
        #     os.path.basename(args.read_align_file_list))[0]
    # if args.output_prefix != "":
    #     output_prefix = args.output_prefix
    #     prefix_flag = True
    # else
    #     prefix_flag = False
    return (sample_list)#, output_prefix , prefix_flag)


def populate_repeat_graphs(specs_path=None, specs_format='v3', reference_fasta=None, flank_size=1000, sample_list=None):
    """Parse EH repeat-specs to create a ReferenceGraph for each Locus ID
    Args:
        specs_path: Path to repeats-specs (JSON file for v3, directory of JSON files for v2.5)
        specs_format: Version of EH v3 or v2.5
        reference_fasta: Path to reference fasta file.
        flank_size: Size of flanks set when aligning reads in EH
        sample_list: Sample list (sample_name, read_align_file) if specs not provided and graph to be inferred from v2.5 alignments.
    Returns:
        dict mapping locus_id to ReferenceGraphs
    """
    repeat_graphs = {}
    if specs_path is None:
        if specs_format == 'v2.5':
            for sample in sample_list:
                eh_json = json.load(open(sample[1]))
                for repeat_id in eh_json:
                    if repeat_id != 'BamStats':
                        repeat_graphs[repeat_id] = ReferenceGraph(eh_json[repeat_id],
                                                                  reference_fasta=reference_fasta,
                                                                  flank_size=flank_size)
    elif specs_format == 'v2.5':
        for spec in os.listdir(specs_path):
            specs_json = json.load(open(os.path.join(specs_path, spec)))
            repeat_graphs[specs_json['RepeatId']] = ReferenceGraph(specs_json,
                                                                   reference_fasta=reference_fasta,
                                                                   flank_size=flank_size)
    elif specs_format in ['v3.0.0-rc1', 'v3']:
        repeat_specs = json.load(open(specs_path))
        if type(repeat_specs) == dict:
            repeat_specs = [repeat_specs]
        for specs_json in repeat_specs:
            if 'RepeatId' in specs_json:
                repeat_id = specs_json['RepeatId']
            elif 'LocusId' in specs_json:
                repeat_id = specs_json['LocusId']
            else:
                repeat_id = [s for s in specs_json['RepeatIds'] if 'ignore' not in s][0]
            refgraph = ReferenceGraph(specs_json,
                                      reference_fasta=reference_fasta,
                                      flank_size=flank_size)
            if ('|') in refgraph.repeat_unit:
                continue
            repeat_graphs[repeat_id] = refgraph
    return repeat_graphs


def get_EH_genotypes(sample_list, file_format='vcf'):
    """Parse EH JSON outputs from list of samples to obtain reported genotypes for each sample.
    """
    # List of genotypes for each sample for each repeat_id
    # genotype[sample][repeat_id] = "genotype1/genotype2"
    genotypes = {}
    for sample in sample_list:
        if sample[2] == None:
            sample_genotypes = None
        elif file_format == 'json':
            eh_json = json.load(open(sample[2]))
            sample_genotypes = {repeat_id: eh_json[repeat_id]['Genotype']
                                for repeat_id in eh_json if repeat_id != 'BamStats'}
            sample_genotypes = {repeat_id: [int(g) for g in sample_genotypes[repeat_id].split(
                '/') if g != ''] for repeat_id in sample_genotypes}
        elif file_format == 'tsv':
            graphEH_out = [l.strip().split() for l in open(sample[2])]
            if len(graphEH_out) == 0 or len(graphEH_out[0]) == 3:
                sample_genotypes = {l[1]:l[2] for l in graphEH_out}
                sample_genotypes = {repeat_id: [int(g) for g in sample_genotypes[repeat_id].split(
                    '/') if g != ''] for repeat_id in sample_genotypes}
            else:
                sample_genotypes = defaultdict(lambda:[], {})
                for l in graphEH_out:
                    sample_genotypes[l[1]].append((l[2], [int(gt) for gt in l[3].split('/')]))
        elif file_format in ['vcf', 'v3', 'v2.5', 'v3.0.0-rc1']:
            sample_genotypes = defaultdict(lambda: [-1, -1], {})
            if os.path.splitext(sample[2])[1] == '.vcf':
                for l in open(sample[2]):
                    ll = l.strip().split()
                    if len(ll) == 0 or ll[0][0] == '#':
                        continue
                    info_dict = {f.split('=')[0]:f.split('=')[1] for f in ll[7].split(';')}
                    if 'VARID' in info_dict:
                        site = info_dict['VARID']
                    elif 'REPID' in info_dict:
                        site = info_dict['REPID']
                    else:
                        continue
                    if 'REF' in info_dict:
                        try:
                            refgt = int(info_dict['REF'])
                        except:
                            continue
                    else:
                        continue
                    gtlist = [refgt] + [int(f.strip('<>STR')) for f in ll[4].split(',') if f != '.']
                    feature_dict = {f[0]:f[1] for f in zip(ll[8].split(':'), ll[9].split(':'))}
                    if 'GT' in feature_dict:
                        genotype = [gtlist[0] if g == '.' else gtlist[int(g)] for g in feature_dict['GT'].split('/')]
                    else:
                        continue
                    genotype.sort()
                    sample_genotypes[site] = genotype
            else:
                sample_genotypes = defaultdict(lambda: [-1, -1], {})
                if file_format in ['v3', 'v3.0.0-rc1']:
                    eh_json = json.load(open(sample[2]))["LocusResults"]
                    for locus_id in eh_json.keys():
                        if 'Variants' not in eh_json[locus_id]:
                            continue
                        for variant_id in eh_json[locus_id]["Variants"].keys():
                            if 'Genotype' not in eh_json[locus_id]["Variants"][variant_id]:
                                continue
                            sample_genotypes[variant_id] = eh_json[locus_id]["Variants"][variant_id]['Genotype']
                else:
                    eh_json = json.load(open(sample[2]))
                    sample_genotypes = {repeat_id: eh_json[repeat_id]['Genotype']
                                        for repeat_id in eh_json if repeat_id != 'BamStats'}   
                sample_genotypes = {repeat_id: [int(g) for g in sample_genotypes[repeat_id].split(
                    '/') if g != ''] for repeat_id in sample_genotypes}    
        genotypes[sample[0]] = sample_genotypes
    return genotypes


def main():
    args = get_args()
    sample_list = get_sample_filepaths(args)
    if args.output_prefix == "":
        output_prefix = ""
    else:
        output_prefix = args.output_prefix + '_'
    if args.reference_fasta is None:
        reference_fasta = None
    else:
        reference_fasta = pysam.FastaFile(args.reference_fasta)
    repeat_graphs = populate_repeat_graphs(args.variant_catalog, args.file_format, reference_fasta=reference_fasta, flank_size=args.flank_size)
    genotypes = get_EH_genotypes(sample_list, args.file_format)
    
    if args.locus_id != '':
        for repid in args.locus_id.split(','):
            if repid not in repeat_graphs:
                raise ValueError('repeat_id "%s" not in repeat-specs' % repid)

    # dict mapping from sample to pileup_coordinates for each (selected) STR region in the sample
    all_sample_read_aligns = {s[0]: parse_read_align(s[1], args.locus_id.split(',') if args.locus_id != '' else [], args.file_format, repeat_graphs, flank_size=args.flank_size) for s in sample_list}
    
    # set of all (selected) sites reported in all samples
    site_list = set([])
    for sample in all_sample_read_aligns:
        for repeat_id in all_sample_read_aligns[sample]:
            site_list.add(repeat_id)

    for repeat_id in site_list:
        if args.node_grouping != "":
            if args.node_grouping == 'NONE':
                node_grouping_list = []
                output_name = ['%s%s-%s.alignment' % (output_prefix, n.reference_region.split('-')[0].replace(':', '-'), n.seq)
                                 for n in repeat_graphs[repeat_id] if n.is_repeat][0]
            elif args.node_grouping == 'ALL':
                node_grouping_list = None
                output_name = ['%s%s-%s.alignment' % (output_prefix, n.reference_region.split('-')[0].replace(':', '-'), n.seq)
                                 for n in repeat_graphs[repeat_id] if n.is_repeat][0]
            else:
                node_grouping_list = [int(i) for i in args.node_grouping.split(',')]
                for node_id in node_grouping_list:
                    if not repeat_graphs[repeat_id][node_id].is_repeat:
                        raise Exception("node id %s in --node_grouping is not a repeat unit" % node_id)
                output_name = ['%s%s-%s.alignment' % (output_prefix, n.reference_region.split('-')[0].replace(':', '-'), n.seq)
                                 for n in repeat_graphs[repeat_id] if n.node_id in node_grouping_list][0]
            
            if args.output_dir != "":
                if not os.path.exists(args.output_dir):
                    os.makedirs(args.output_dir)
                output_name = os.path.join(args.output_dir, output_name)
            pileup_list = get_pileup_list({s: all_sample_read_aligns[s][repeat_id] for s in all_sample_read_aligns},
                                        [s[0] for s in sample_list], repeat_graphs[repeat_id],
                                        node_grouping_list=node_grouping_list, show_insertions=args.show_insertions)
            # if prefix_flag and len(site_list) == 1:
            plot_pileup(repeat_id, pileup_list, genotypes, repeat_graphs[repeat_id], args.greyscale, output_name, args.title_prefix, args.show_read_names, show_insertions=args.show_insertions, flank_clip_size=args.flank_clip_size, dpi=args.dpi, pdf=args.pdf)
            # else:
            #     plot_pileup(repeat_id, pileup_list, genotypes, repeat_graphs[repeat_id],
            #                 args.greyscale, output_prefix + "_" + repeat_id, args.title_prefix, args.show_read_names, show_insertions=args.show_insertions, flank_clip_size=args.flank_clip_size, dpi=args.dpi, pdf=args.pdf)
        else:
            for n in repeat_graphs[repeat_id]:
                if not n.is_repeat:
                    continue
                node_grouping_list = [n.node_id]
                output_name = '%s%s-%s.alignment' % (output_prefix, n.reference_region.split('-')[0].replace(':', '-'), n.seq)
                if args.output_dir != "":
                    if not os.path.exists(args.output_dir):
                        os.makedirs(args.output_dir)
                    output_name = os.path.join(args.output_dir, output_name)
                pileup_list = get_pileup_list({s: all_sample_read_aligns[s][repeat_id] for s in all_sample_read_aligns},
                                            [s[0] for s in sample_list], repeat_graphs[repeat_id],
                                            node_grouping_list=node_grouping_list, show_insertions=args.show_insertions)
                plot_pileup(repeat_id, pileup_list, genotypes, repeat_graphs[repeat_id], args.greyscale, output_name, args.title_prefix, args.show_read_names, show_insertions=args.show_insertions, flank_clip_size=args.flank_clip_size, dpi=args.dpi, pdf=args.pdf)
                


if __name__ == '__main__':
    main()
