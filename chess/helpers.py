import io
import gzip
import numpy as np
import os
import intervaltree
from future.utils import string_types
from collections import defaultdict
import logging
logger = logging.getLogger(__name__)


class GenomicRegion(object):
    """
    Class representing a genomic region.
    .. attribute:: chromosome
        Name of the chromosome this region is located on
    .. attribute:: start
        Start position of the region in base pairs
    .. attribute:: end
        End position of the region in base pairs
    .. attribute:: strand
        Strand this region is on (+1, -1)
    .. attribute:: ix
        Index of the region in the context of all genomic
        regions.
    """

    def __init__(self, start, end, chromosome=None, ix=None, strand=None):
        """
        Initialize this object.
        :param start: Start position of the region in base pairs
        :param end: End position of the region in base pairs
        :param chromosome: Name of the chromosome this region is located on
        :param ix: Index of the region in the context of all genomic
                   regions.
        """
        self.start = start
        self.end = end
        self.chromosome = chromosome
        self.ix = ix
        self.strand = strand

    @classmethod
    def from_string(cls, region_string):
        """
        Convert a string into a :class:`~GenomicRegion`.
        This is a very useful convenience function to quickly
        define a :class:`~GenomicRegion` object from a descriptor
        string.
        :param region_string: A string of the form
                              <chromosome>[:<start>-<end>[:<strand>]]
                              (with square brackets indicating optional
                              parts of the string). If any optional
                              part of the string is omitted, intuitive
                              defaults will be chosen.
        :return: :class:`~GenomicRegion`
        """
        chromosome = None
        start = None
        end = None

        # strip whitespace
        no_space_region_string = "".join(region_string.split())
        fields = no_space_region_string.split(':')

        if len(fields) > 3:
            raise ValueError(
                ("Genomic range string must be of the form"
                 " <chromosome>[:<start>-<end>:[<strand>]]"))

        # there is chromosome information
        if len(fields) > 0:
            chromosome = fields[0]

        # there is range information
        if len(fields) > 1 and fields[1] != '':
            start_end_bp = fields[1].split('-')
            if len(start_end_bp) > 0:
                try:
                    start = int(start_end_bp[0])
                except ValueError:
                    raise ValueError("Start of genomic range must be integer")

            if len(start_end_bp) > 1:
                try:
                    end = int(start_end_bp[1])
                except ValueError:
                    raise ValueError("End of genomic range must be integer")

                if not end > start:
                    raise ValueError(
                        "The end coordinate must be bigger than the start.")

        return cls(start=start, end=end, chromosome=chromosome)

    def overlaps(self, region):
        """
        Check if this region overlaps with the specified region.
        :param region: :class:`~GenomicRegion` object or string
        """
        if isinstance(region, str):
            region = GenomicRegion.from_string(region)

        if region.chromosome != self.chromosome:
            return False

        if region.start <= self.end and region.end >= self.start:
            return True
        return False

    def contains(self, region):
        """
        Check if the specified region is completely contained in this region.
        :param region: :class:`~GenomicRegion` object or string
        """
        if isinstance(region, str):
            region = GenomicRegion.from_string(region)

        if region.chromosome != self.chromosome:
            return False

        if region.start >= self.start and region.end <= self.end:
            return True
        return False

    def _equals(self, region):
        if region.chromosome != self.chromosome:
            return False
        if region.start != self.start:
            return False
        if region.end != self.end:
            return False
        return True

    def __eq__(self, other):
        return self._equals(other)

    def __ne__(self, other):
        return not self._equals(other)

    def __str__(self):
        return '{}:{}-{}:{}'.format(
            self.chromosome, self.start, self.end, self.strand)

    def copy(self):
        """
        Make a copy of this GenomicRegion object.
        """
        return GenomicRegion(
            chromosome=self.chromosome, start=self.start, end=self.end)


def is_gzipped(file_name):
    if file_name.endswith('.gz') or file_name.endswith('.gzip'):
        return True
    return False


def _open(file_name, mode='r'):
    if is_gzipped(file_name):
        return io.TextIOWrapper(gzip.open(file_name, mode=mode))
    return open(file_name, mode=mode)


def load_regions(file_name, sep=None, ignore_ix=False):
    """
    Load regions from regions bed file.
    :param file_name: Path to regions bed file.
    :param sep: Delimiter in regions bed, defaults to None (splits by tab)
    :returns: Tuple:
                (
                    List of :class:`~GenomicRegion` objects,
                    ix_converter dict: {region_id: position_in_list},
                    ix2region dict: {region_id: [chromosome, start, end]}
                (
    """
    regions = []
    ix2reg = {}
    ix_converter = None
    with _open(file_name, 'r') as f:
        for i, line in enumerate(f):
            line = line.rstrip()
            fields = line.split(sep)
            if len(fields) > 2:
                chromosome = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                ix = i
                ix2reg[ix] = [chromosome, int(fields[1]), int(fields[2])]
                if not ignore_ix and len(fields) > 3 and fields[3] != '.':  # HicPro
                    if ix_converter is None:
                        ix_converter = dict()
                    if fields[3] in ix_converter:
                        raise ValueError(
                            "name column in region BED must "
                            "only contain unique values! ({})".format(
                                fields[3]))
                    ix_converter[fields[3]] = ix
                regions.append(
                    GenomicRegion(
                        chromosome=chromosome, start=start, end=end, ix=ix))

    if ix_converter is not None:
        ix_converter_rev = {str(v): int(k) for k, v in ix_converter.items()}
        ix2reg = {ix_converter_rev[str(k)]: v for k, v in ix2reg.items()}
    return regions, ix_converter, ix2reg


def region_interval_trees(regions):
    chromosome_intervals = defaultdict(list)
    for region in regions:
        interval = intervaltree.Interval(region.start - 1, region.end, data=region)
        chromosome_intervals[region.chromosome].append(interval)

    region_trees = dict()
    for chromosome, intervals in chromosome_intervals.items():
        region_trees[chromosome] = intervaltree.IntervalTree(intervals)

    return region_trees


def sub_regions(regions, region):
    """
    Get regions from a list the overlap with another region.
    :param regions: List of :class:`~GenomicRegion`
    :param region: :class:`~GenomicRegion` used for overlap calculation
    """
    if isinstance(region, string_types):
        region = GenomicRegion.from_string(region)

    # interval tree?
    try:
        rt = regions[region.chromosome]

        ixs = []
        sr = []
        for interval in rt[region.start - 1: region.end]:
            hit_region = interval.data
            ixs.append(hit_region.ix)
            sr.append(hit_region)
        start_ix = min(ixs)
        end_ix = max(ixs)

    # region list
    except TypeError:
        sr = []
        start_ix = None
        end_ix = None
        for i, r in enumerate(regions):
            if (r.chromosome == region.chromosome
                    and r.start <= region.end and r.end >= region.start):
                if start_ix is None:
                    start_ix = i
                end_ix = i
                sr.append(r.copy())
            else:
                if end_ix is not None:
                    break
    except IndexError:
        sr = []
        start_ix = None
        end_ix = None

    if len(sr) == 0:
        raise ValueError("Region not found in dataset! {}:{}-{}".format(
            region.chromosome, region.start, region.end))

    return sr, start_ix, end_ix


def sub_matrix_regions(hic_matrix, regions, region):
    """
    Get a square sub Hi-C matrix that overlaps a given region.
    :param hic_matrix: A square numpy array
    :param regions: List of :class:`~GenomicRegion`
    :param region: :class:`~GenomicRegion` used for overlap calculation
    """
    sr, start_ix, end_ix = sub_regions(regions, region)

    if start_ix is None:
        return np.empty((0, 0)), sr

    return np.copy(hic_matrix[start_ix:end_ix+1, start_ix:end_ix+1]), sr


def sub_matrix_from_edges_dict(edges, regions, region, default_weight=0.0):
    """
    Get a square sub Hi-C matrix that overlaps a given region.
    :param edges: A dict of the form { source_ix1: { sink_ix1: weight1, sink_ix2: weight2, ...}, ...}
    :param regions: List of :class:`~GenomicRegion` or dict of :class:`~intervaltree.IntervalTree`
    :param region: :class:`~GenomicRegion` used for overlap calculation
    :param default_weight: Default value with which to fill matrix
    """
    sr, start_ix, end_ix = sub_regions(regions, region)

    size = end_ix - start_ix + 1
    m = np.full((size, size), default_weight, dtype='float64')
    for i in range(start_ix, end_ix + 1):
        for j in range(i, end_ix + 1):
            try:
                m[i - start_ix, j - start_ix] = edges[i][j]
                m[j - start_ix, i - start_ix] = edges[i][j]
            except KeyError:
                pass

    return m, sr


def edges_from_sparse_matrix(file_name, ix_converter=None, sep="\t"):
    nan_values = 0
    total_values = 0

    with _open(file_name, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            total_values += 1

            line = line.rstrip()
            fields = line.split(sep)
            if ix_converter is None:
                source, sink, weight = int(fields[0]), int(fields[1]), float(fields[2])
            else:
                source = ix_converter[fields[0]]
                sink = ix_converter[fields[1]]
                weight = float(fields[2])

            if not np.isfinite(weight):
                nan_values += 1
                continue

            if source <= sink:
                yield source, sink, weight
            else:
                yield sink, source, weight

    if nan_values/total_values > 0.05:
        logger.warning("Could not import more than 5% of matrix entries ({}/{}), "
                       "as they were not finite".format(nan_values, total_values))


def edges_dict_from_sparse(edges):
    """Constructs edges_dict from edges parsed
    from sparse format data.

    Args:
        edges (iterator): iterator over (sink, source, weight) tuples

    Returns:
        dict: Nexted dictionary that allows
              queries like: dict[source][sink] = weight
    """
    edges_d = defaultdict(dict)
    for source, sink, weight in edges:
        edges_d[source][sink] = weight
    return edges_d


def edges_dict_from_fanc(hic):
    """Constructs edges_dict from edges
    present in a fanc hic object.

    Args:
        hic (fanc.Hic): Object created by fanc.load()

    Returns:
        dict: Nexted dictionary that allows
              queries like: dict[source][sink] = weight
    """
    edges_dict = defaultdict(dict)
    for e in hic.edges(lazy=True):
        edges_dict[e.source][e.sink] = e.weight
    return edges_dict


def oe_edges_dict_from_fanc(hic):
    """Constructs edges_dict from edges
    present in a fanc hic object. Weights are
    transformed to observed / expected values.

    Args:
        hic (fanc.Hic): Object created by fanc.load()

    Returns:
        dict: Nexted dictionary that allows
              queries like: dict[source][sink] = weight
    """
    edges_dict = defaultdict(dict)
    for e in hic.edges(oe=True, lazy=True):
        edges_dict[e.source][e.sink] = e.weight
    return edges_dict


def load_sparse_matrix(file_name, ix_converter=None, sep="\t"):
    return list(edges_from_sparse_matrix(file_name, ix_converter, sep))


def load_matrix(file_name, size=None, sep=None, ix_converter=None, is_oe=False):
    """Load Numpy matrix from Hi-C matrix file in sparse format.

    :param file_name: Path to Hi-C matrix file in sparse format.
    :param size: Size of the matrix (number of bins / edge length),
                 defaults to None
    :param sep: Delimiter in matrix file, defaults to None (splits by tab)
    :param ix_converter: ix_converter dict corresponding to the matrix file,
                         as produced by :func:`load_regions`, defaults to None
    :returns: Numpy matrix.
    """

    if size is None:
        raise ValueError("Must provide matrix size!")

    if is_oe:
        m = np.ones((size, size))
    else:
        m = np.zeros((size, size))

    for source, sink, weight in edges_from_sparse_matrix(
                                file_name, ix_converter=ix_converter, sep=sep):
        m[source, sink] = weight
        m[sink, source] = weight

    return m


def load_pairs(file_name, sep=None):
    """Load region pairs from bedpe file.

    Read input bedpe file and convert to dict mapping the comparison ID to the
    reference to query regions.
    :param file_name: Path to input bedpe file. Expected columns:
        chromosome_reference, start_reference, end_reference,
        chromosome_query, start_query, end_query,
        comparison_id, dummy column (not used), strand_reference, strand_query
    :param sep: Delimiter in bedpe file, defaults to None (splits by tab)
    :returns: Dict of form {ID: (reference_region, query_region)}, with
              *region being :class:`GenomicRegion` objects.
    """
    pairs = dict()
    for ID, cref, cqry in load_pairs_iter(file_name, sep=sep):
        pairs[ID] = (cref, cqry)
    return pairs


def load_pairs_iter(file_name, sep=None):
    """Load region pairs from bedpe file.

    Read input bedpe file and convert to dict mapping the comparison ID to the
    reference to query regions.
    :param file_name: Path to input bedpe file. Expected columns:
        chromosome_reference, start_reference, end_reference,
        chromosome_query, start_query, end_query,
        comparison_id, dummy column (not used), strand_reference, strand_query
    :param sep: Delimiter in bedpe file, defaults to None (splits by tab)
    :returns: Dict of form {ID: (reference_region, query_region)}, with
              *region being :class:`GenomicRegion` objects.
    """
    pairs = {}
    with _open(file_name, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            sline = line.rstrip()
            fields = sline.split(sep)
            if len(fields) != 10:
                raise ValueError((
                    "10 columns expected but {} found"
                    " in bedpe input"
                    ).format(len(fields)))
            c1, s1, e1, c2, s2, e2, ID, _, str1, str2 = line.split(sep)
            if ID in pairs:
                raise ValueError((
                    "Pair ID {} not unique in bedpe input!"))
            cref = GenomicRegion(
                chromosome=str(c1), start=int(s1) + 1,
                end=int(e1), strand=str(str1))
            cqry = GenomicRegion(
                chromosome=str(c2), start=int(s2) + 1,
                end=int(e2), strand=str(str2))
            yield ID, cref, cqry


def load_pairs_for_chroms(file_name, chromosomes, sep=None):
    """Load region pairs from bedpe file.

    Read input bedpe file and convert to dict mapping the comparison ID to the
    reference to query regions. Return only pairs which are on the specified
    chromosomes.

    :param file_name: Path to input bedpe file. Expected columns:
        chromosome_reference, start_reference, end_reference,
        chromosome_query, start_query, end_query,
        comparison_id, dummy column (not used), strand_reference, strand_query
    :param sep: Delimiter in bedpe file, defaults to None (splits by tab)
    :returns: Dict of form {ID: (reference_region, query_region)}, with
              *region being :class:`GenomicRegion` objects.
    """
    pairs = set()
    res = []
    dropped = 0
    with _open(file_name, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            sline = line.rstrip()
            fields = sline.split(sep)
            if len(fields) != 10:
                raise ValueError((
                    "10 columns expected but {} found"
                    " in bedpe input"
                    ).format(len(fields)))
            c1, s1, e1, c2, s2, e2, ID, _, str1, str2 = line.split(sep)
            if ID in pairs:
                raise ValueError((
                    "Pair ID {} not unique in bedpe input!"))
            pairs.add(ID)
            if any(c not in chromosomes for c in (c1, c2)):
                dropped += 1
            else:
                cref = GenomicRegion(
                    chromosome=str(c1), start=int(s1) + 1,
                    end=int(e1), strand=str(str1))
                cqry = GenomicRegion(
                    chromosome=str(c2), start=int(s2) + 1,
                    end=int(e2), strand=str(str2))
                res.append((ID, cref, cqry))
    return res, dropped


def chunks(l, p):
    n = int(np.ceil(len(l)/p))
    """Yield successive len(l)/p-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def read_chromosome_sizes(file_name):
    chrom_sizes = {}
    with open(os.path.expanduser(file_name), 'r') as chrom_sizes_file:
        for line in chrom_sizes_file:
            line = line.rstrip()
            if line != '':
                chromosome, chromosome_length = line.split("\t")
                chrom_sizes[chromosome] = int(chromosome_length)
    return chrom_sizes


def load_contacts(matrix_file, regions_file=None):
    import fanc
    from chess.oe import observed_expected
    try:
        # try loading via fanc
        reference_loaded = fanc.load(matrix_file)
        edges = oe_edges_dict_from_fanc(reference_loaded)
        regions = reference_loaded.regions
        region_trees = region_interval_trees(
            regions)
    except ValueError as initial_error:
        print(initial_error)
        try:
            assert regions_file is not None, (
                "Regions file needs to be"
                "specified for sparse input.")
            regions, _ix_converter, _ = load_regions(
                regions_file)
            region_trees = region_interval_trees(
                regions)
            _, reference_oe = observed_expected(
                regions_file, matrix_file)
            edges = edges_dict_from_sparse(reference_oe)
        except AssertionError as error:
            raise ValueError(error)
    return edges, region_trees, regions


def load_oe_contacts(matrix_file, regions_file=None):
    import fanc
    try:
        # try loading via fanc
        reference_loaded = fanc.load(matrix_file)
        edges = edges_dict_from_fanc(reference_loaded)
        regions = reference_loaded.regions
        region_trees = region_interval_trees(
            regions)
    except ValueError:
        try:
            assert regions_file is not None
            regions, ix_converter, _ = load_regions(
                regions_file)
            region_trees = region_interval_trees(
                regions)
            edges = edges_dict_from_sparse(
                edges_from_sparse_matrix(
                    matrix_file,
                    ix_converter))
        except AssertionError:
            raise ValueError
    return edges, region_trees, regions
