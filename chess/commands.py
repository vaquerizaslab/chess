import argparse
import textwrap
import sys


def chess_parser():
    usage = '''\
        chess <command> [options]

        Commands:
            sim             Calculate structural similarity
            oe              Transform a Hi-C matrix to an observed/expected matrix
            pairs           Make window pairs for chess genome scan
            background   Generate background region BED files
            filter          Filter output of chess sim and save as BED file
            extract         Extract specific regions that are significantly different
            crosscorrelate  Get structural clusters from the extracted submatrices

        Run chess <command> -h for help on a specific command.
        '''
    parser = argparse.ArgumentParser(
        description='CHESS: Compare Hi-C Experiments using Structural Similarity',
        usage=textwrap.dedent(usage)
    )

    parser.add_argument(
        '--version', dest='print_version',
        action='store_true',
        help='Print version information'
    )
    parser.set_defaults(print_version=False)

    parser.add_argument(
        '--verbose', '-v', dest='verbosity',
        action='count',
        default=0,
        help='Set verbosity level: Can be chained like '
             '"-vvv" to increase verbosity. Default is to show '
             'errors, warnings, and info messages (same as "-vv"). '
             '"-v" shows only errors and warnings, "-vvv" shows errors '
             'warnings, info, and debug messages in addition.')

    parser.add_argument(
        '-l', '--log-file', dest='log_file',
        help='Path to file in which to save log.'
    )

    parser.add_argument('command', nargs='?', help='Subcommand to run')

    return parser


def sim_parser():
    parser = MyParser(
            description='Compare structures between pairs of chromatin '
                        'contact matrices using the structural similarity '
                        'index. If --background-query or --background-regions '
                        'are specified, compute p-value and z-value '
                        'after obtaining a background distribution '
                        'of similarity values of the reference in the pair '
                        'to specfied background regions. '
                        'The input matrices are expected to be balanced, '
                        'e.g. by Knight-Ruiz matrix balancing. '
                        'The contacts are automatically converted internally '
                        'to observed / expected values. In your input files '
                        'are already observed / expected transformed, '
                        'set the --oe-input flag.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'reference_contacts',
        type=str,
        help='Balanced contact matrix for the reference sample '
             'in one of the following formats: fanc .hic, '
             'juicer .hic@<resolution>,'
             'cooler .cool@<resolution> or .mcool@<resolution>, '
             'sparse format '
             '(each line: '
             '<row region index> <column region index> <matrix value>). '
             'If the file is in sparse format, the corresponding regions '
             ' BED file needs to be passed via --reference-regions.')

    parser.add_argument(
        '--reference-regions', dest="reference_regions",
        type=str,
        help='BED file (no header) with regions corresponding to '
             'the number of rows in the provided reference matrix.')

    parser.add_argument(
        'query_contacts',
        type=str,
        help='Balanced contact matrix for the query sample '
             'in one of the following formats: fanc .hic, '
             'juicer .hic@<resolution>,'
             'cooler .cool@<resolution> or .mcool@<resolution>, '
             'sparse format '
             '(each line: '
             '<row region index> <column region index> <matrix value>). '
             'If the file is in sparse format, the corresponding regions '
             ' BED file needs to be passed via --query-regions.')

    parser.add_argument(
        '--query-regions', dest="query_regions",
        type=str,
        help='BED file (no header) with regions corresponding to '
             'the number of rows in the provided query matrix.')

    parser.add_argument(
        'pairs',
        type=str,
        help='Region pairs to compare. '
             'Expected to be in bedpe format with chrom1, start1, ... '
             'corresponding to reference and chrom2, start2, ... to query.')

    parser.add_argument(
        'out',
        type=str,
        help='Path to outfile.')

    parser.add_argument(
        '--background-regions', dest='background_regions',
        type=str,
        help='BED file with regions to be used for background calculations. '
             'If provided, CHESS will generate Z-scores and P-values for '
             'similarities.')

    parser.add_argument(
        '--background-query', dest='background_query',
        default=False,
        action='store_true',
        help='Use every region of the same size as the '
             'reference from the query genome as background. '
             'Useful, for example, as background for inter-species '
             'comparisons.')

    parser.add_argument(
        '-p', dest='threads',
        type=int,
        default=1,
        help='Number of cores to use.')

    parser.add_argument(
        '--keep-unmappable-bins',
        action='store_true',
        default=False,
        help='Disable deletion of unmappable bins.')

    parser.add_argument(
        '--mappability-cutoff',
        type=float,
        default=0.1,
        help='Low pass threshold for fraction of unmappable bins. '
             'Matrices with a higher content of unmappable bins '
             'will not be considered. '
             'Unmappable bins will be deleted from matrices '
             'passing the filter.')

    parser.add_argument(
        '-r', '--relative-windowsize',
        type=float,
        default=1,
        help='Relative window size value '
             'for the win_size param in the ssim function. '
             'Fraction of the matrix size.')

    parser.add_argument(
        '-a', '--absolute-windowsize',
        type=int,
        default=None,
        help='Absolute window size value in bins '
             'for the win_size param in the ssim function. '
             'Overwrites -r.')

    parser.add_argument(
        '--oe-input',
        action='store_true',
        default=False,
        help='Use if input contacts are already observed / expected '
             'transformed.')

    parser.add_argument(
        '--limit-background',
        action='store_true',
        default=False,
        help='Restrict background computation to the syntenic / paired '
             'chromosome as indicated in the pairs file.')

    return parser


def oe_parser():
    parser = MyParser(
        description='Convert a sparse Hi-C matrix to '
                    'observed/expected format.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'input_matrix',
        type=str,
        help='Balanced Hi-C matrix in sparse format. '
             '(each line: '
             '<row region index> <column region index> <matrix value>)')

    parser.add_argument(
        'regions',
        type=str,
        help='BED file (no header) with regions corresponding to '
             'the number of rows in the provided reference matrix.')

    parser.add_argument(
        'output_matrix',
        type=str,
        help='Obs/exp transformed matrix (same as input matrix format)')

    return parser


def pairs_parser():
    parser = MyParser(
        description='Make window pairs for CHESS genome scan. '
                    'Write all positions of a sliding window of specified '
                    'size with specified step in the specified genome to '
                    'the outfile which can be directly used to run CHESS sim.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'genome',
        type=str,
        help='UCSC genome identifier (as recognized by pybedtools), '
             'or path to tab-separated chrom sizes file with columns '
             '<chromosome name> <chromosome size>. '
             'Will use the path only if no USCS entry with that name is found'
             ', or --file-input is specified.')

    parser.add_argument(
        'window',
        type=int,
        help='''Window size in base pairs''')

    parser.add_argument(
        'step',
        type=int,
        help='''Step size in base pairs''')

    parser.add_argument(
        'output',
        type=str,
        help='''Path to output file''')

    parser.add_argument(
        '--file-input',
        action='store_true',
        default=False,
        help='Will not check for USCS entry of genome'
             'input with pybedtools if set'
        )

    parser.add_argument(
        '--chromosome',
        type=str,
        help='Produce window pairs only for the specified chromosome')

    return parser


def background_parser():
    parser = MyParser(
        description='Generate BED file with regions to be used'
                    'in CHESS background calculations.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'genome_or_region',
        type=str,
        help='UCSC genome identifier (as recognized by pybedtools), '
             'OR path to tab-separated chrom sizes file with columns '
             '<chromosome name> <chromosome size> OR region identifier in '
             'the format <chromosome>:<start>-<end>. '
             'Will try options in the order listed.')

    parser.add_argument(
        'window',
        type=int,
        help='Window size in base pairs')

    parser.add_argument(
        'step',
        type=int,
        help='Step size in base pairs')

    parser.add_argument(
        'output',
        type=str,
        help='Path to output file')

    parser.add_argument(
        '--strand',
        type=str,
        help='[+/-] .Generate regions on this strand only. '
             'Default: both strands')

    return parser


def extract_parser():
    parser = MyParser(
        description='Extract the specific regions that are different '
                    'between the regions identified by CHESS.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'pairs',
        type=str,
        help='Region pairs that have been identified to '
             'contain structural differences. Expected to be in '
             'bedpe format with chrom1, start1, ... '
             'corresponding to reference and chrom2, start2, ... to query.')

    parser.add_argument(
        'reference_contacts',
        type=str,
        help='Balanced contact matrix for the reference sample '
             'in one of the following formats: fanc .hic, '
             'juicer .hic@<resolution>,'
             'cooler .cool@<resolution> or .mcool@<resolution>, '
             'sparse format '
             '(each line: '
             '<row region index> <column region index> <matrix value>). '
             'If the file is in sparse format, the corresponding regions '
             ' BED file needs to be passed via --reference-regions.')

    parser.add_argument(
        '--reference-regions', dest="reference_regions",
        type=str,
        help='BED file (no header) with regions corresponding to '
             'the number of rows in the provided reference matrix.')

    parser.add_argument(
        'query_contacts',
        type=str,
        help='Balanced contact matrix for the query sample '
             'in one of the following formats: fanc .hic, '
             'juicer .hic@<resolution>,'
             'cooler .cool@<resolution> or .mcool@<resolution>, '
             'sparse format '
             '(each line: '
             '<row region index> <column region index> <matrix value>). '
             'If the file is in sparse format, the corresponding regions '
             ' BED file needs to be passed via --query-regions.')

    parser.add_argument(
        '--query-regions', dest="query_regions",
        type=str,
        help='BED file (no header) with regions corresponding to '
             'the number of rows in the provided query matrix.')

    parser.add_argument(
        'out',
        type=str,
        help='Path to output directory.')

    parser.add_argument(
        '--windowsize',
        type=int,
        default=3,
        help='Window size to average the bins according to their spatial '
             'closeness and their radiometric similarity, '
             'by default the windows size is the 3 x 3 bins. '
             'Larger values will average bins with larger differences.')

    parser.add_argument(
        '--sigma-spatial', dest="sigma_spatial",
        type=int,
        default=3,
        help='Gaussian function of the Euclidean distance between '
             'two bins and its standard deviation. '
             'Larger values will average bins with larger differences.')

    parser.add_argument(
        '--size-medianfilter', dest="size_medianfilter",
        type=int,
        default=9,
        help='Windows size used to scan and smooth the contained bins. '
             'Higher values will smooth larger figures, '
             'while smaller values will consider subtle signals (i.e. loops).')

    parser.add_argument(
        '--closing-square', dest="closing_square",
        type=int,
        default=8,
        help='Side length of the square used to remove noise, '
             'and fill structures. '
             'Larger values will enclose larger structures and remove '
             'punctuate or looping structures.')

    return parser


def crosscorrelate_parser():
    parser = MyParser(
        description='2D crosscorrelation of the specific substructures '
                    'extracted from the '
                    'regions with differences in their chromatin structures.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'extracted_file',
        type=str,
        help='Output from extract sub-command.')

    parser.add_argument(
        'pairs',
        type=str,
        help='Region pairs that have been identified to '
             'contain structural differences. Expected to be in '
             'bedpe format with chrom1, start1, ... '
             'corresponding to reference and chrom2, start2, ... to query.')

    parser.add_argument(
        'outdir',
        type=str,
        help='Path to output directory.')

    return parser


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
