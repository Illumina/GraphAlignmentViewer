# GraphAlignmentViewer:
Visualize the pileup of read alignments to a STR locus (or any path in the sequence graph) for one or more samples.
## Description:
  Genotyping STR loci is a difficult problem due to multiple reasons including flanks with similar genomic sequences, difficulty in aligning reads sequencing errors and sample contamination. ExpansionHunter is a tool that can automatically predict the genotypes of STR loci based on read alignments. However, due to the difficult nature of the problem, there may be errors in the predicted genotype. GraphAlignmentViewer, is a standalone Python3 script that creates a visualization of all read alignments to an STR locus to enable visual inspection of the genotype predictions made by ExpansionHunter. The script can visualize one or more samples, including trios to compare genotype calls in multiple samples. For each sample the script requires genotypes predicted by ExpansionHunter in BAM format for EH version 3 or YAML format for EH version 2.5.

The image generated by GraphAlignmentViewer displays alignments of reads to a modified reference genome. The modified reference corresponds to sequence of the STR genotype with the largest count for each repeat unit in the realigned BAM file from ExpansionHunter. As a result, the modified reference and aligned reads are (vertically) sectioned into consecutive STR and non-STR regions (graph nodes). Next, GraphAlignmentViewer groups all aligned reads from a sample into (horizontal) blocks. Each block consists of all reads with a unique STR count for the selected repeat unit. Within each block, GraphAlignmentViewer displays the sequence of modified reference in the top row. This is followed by the aligned sequences of reads with the corresponding STR count. This grouping allows a user to quickly compare and count all similar reads to each other as well as identify difference with the modified reference. Insertions and deletions are denoted by "^" and "-" symbols respectively. If reads from multiple samples are provided, then the reads from each sample are simply stacked vertically from top to bottom. Finally, GraphAlignmentViewer displays the genotyped reported by ExpansionHunter and the locus structure in the title of the plot.

Of note, if the STR is longer than read length, then the sample would contains reads aligned completely inside the STR region (in-repeat reads). In such a case, the modified reference corresponds to the smallest genotype with whole repeat units larger than the read length. In case of a large repeat expansion, the user will be able to observe a block with large number of in-repeat reads aligned to the section corresponding to the expanded repeat unit.

**Note:** The current version only supports graph specifications consisting of linear paths and self-loops, but not branching nodes.

**Example of a visualization generated by GraphAlignmentViewer:**
![Sample image](/images/GB18_mother.png)

## Usage:
`python3 GraphAlignmentViewer.py --help`
### Use case 1: Single sample
`python3 GraphAlignmentViewer.py --variant_catalog VARIANT_CATALOG --read_align READ_ALIGN_FILE [--gt_file GT_FILE]`
### Use case 2: Comparison of genotypes from multiple samples
`python3 GraphAlignmentViewer.py --variant_catalog VARIANT_CATALOG --read_align_list READ_ALIGN_FILE_LIST`
## Requirements:
* Python3
* Matplotlib
* Pysam
* Numpy
* PyYAML (if visualizing output from versions older than v3)
## Inputs:
### Common inputs:
1. `VARIANT_CATALOG`: Path to variant catalog JSON file used to run ExpansionHunter
### Use case 1: Single sample
1. `READ_ALIGN_FILE`: Read alignment file generated be ExpansionHunter. BAM for v3, YAML for v2.5.
2. `GT_FILE` (optional): VCF or JSON output of genotype calls from ExpansionHunter
### Use case 2: Comparison of genotypes from multiple samples
1. `READ_ALIGN_FILE_LIST`: A 2/3 column text CSV file containing the list of EH output for all samples:
* Column1: Sample name,
* Column2: Read alignment from EH output similar to `READ_ALIGN_FILE`,
* Column3 (optional): VCF or JSON file from EH output

The pileups for the samples are shown in the order reported in the input file.  
**Note:** If the paths to the read alignment BAM/YAML and VCF files are not absolute, they are chosen relative to the location of the sample list file.


## Optional Arguments:
| Option | Argument | Default | Description |
|:--:|:--:|:--:|:--|
|`--file_format` | `FILE_FORMAT` | `v3` | Format of read alignments from EH. [`v3`: BAM, `v2.5`: YAML] |
|`--locus_id` | `LOCUS_ID` | Plot pileups for all loci | Comma-separated list of locus IDs for which to plot pileup |
|`--greyscale` | `-`      | Nucleotides colored in IGV color scheme | Show nucleotides in greyscale: high quality match - black, low quality match - grey, mismatch - red |
|`--show_read_names` | `-` | Do not display read names | Display read names next to the read alignment |
|`--show_insertions` | `-` | Do not display inserted sequences | Display full sequences of insertions |
|`--output_prefix` | `OUTPUT_PREFIX` | No prefix. Output filename(s): `<CHROM>-<START>-<REPEATUNIT>.alignment.png`(`.pdf`) | Prefix of output file. Output filename(s): `<OUTPUT_PREFIX>_<CHROM>-<START>-<REPEATUNIT>.alignment.png`(`.pdf`) corresponding to the position of the first repeat unit in the node grouping. If node grouping is `NONE` or `ALL`, then position corresponds to the first repeat unit in the locus. |
|`--output_dir` | `OUTPUT_DIR` | `Current working directory` | Output directory |
|`--title_suffix` | `TITLE_PREFIX` | "" | Prefix text to be appended to title of the plot |
|`--reference_fasta` | `REFERENCE_FASTA` | Represent flanks with 'N's | Indexed FASTA file for reference sequence |
|`--node_grouping` | `NODE_GROUPING` | Create a separate image for each repeat unit | Comma-separated list of node indices (left flank=`0`) to group and sort reads by genotype. `NONE`: sort reads only by position, `ALL`: group by all repeat nodes from left to right. |
|`--region_extension_length` | `REGION_EXTENSION_LENGTH` (`INT`) | `1000` | Size of nodes flanking the region structure used for generating the read alignments |
|`--region_extension_clip_length` | `REGION_EXTENSION_CLIP_LENGTH` (`INT`) | `20` | Number of basepairs of flanking regions to display. `-1`: Infer from maximum span of reads overlapping the locus. |
| `--dpi` | `DPI` (`INT`) | `100` | Resolution of output PNG image |
| `--pdf` | `-` | Output PNG image | Output PDF vector graphics image instead of PNG |




## Outputs:
The script produces 1 file per repeat unit or `NODE_GROUPING`: `<OUTPUT_DIR>/<CHROM>-<START>-<REPEATUNIT>.alignment.png`. If the `--pdf` flag is set then it produces the file `<OUTPUT_DIR>/<CHROM>-<START>-<REPEATUNIT>.alignment.pdf` instead.

If the flag `--output_prefix` is set then the output file name is `<OUTPUT_DIR>/<OUTPUT_PREFIX>_<CHROM>-<START>-<REPEATUNIT>.alignment.png` or `<OUTPUT_DIR>/<OUTPUT_PREFIX>_<CHROM>-<START>-<REPEATUNIT>.alignment.pdf` 
