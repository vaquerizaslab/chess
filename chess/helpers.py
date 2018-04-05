import io
import gzip
import numpy as np
import os
from future.utils import string_types


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


def load_regions(file_name, sep=None):
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
                if len(fields) > 3 and fields[3] != '.':  # HicPro
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


def sub_regions(regions, region):
    """
    Get regions from a list the overlap with another region.
    :param regions: List of :class:`~GenomicRegion`
    :param region: :class:`~GenomicRegion` used for overlap calculation
    """
    if isinstance(region, string_types):
        region = GenomicRegion.from_string(region)

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

    if start_ix is None or end_ix is None:
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


def edges_from_sparse_matrix(file_name, ix_converter=None, sep="\t"):
    with _open(file_name, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            fields = line.split(sep)
            if ix_converter is None:
                source, sink, weight = (
                    int(fields[0]), int(fields[1]), float(fields[2]))
            else:
                source = ix_converter[fields[0]]
                sink = ix_converter[fields[1]]
                weight = float(fields[2])

            if source <= sink:
                yield source, sink, weight
            else:
                yield sink, source, weight


def load_sparse_matrix(file_name, ix_converter=None, sep="\t"):
    return list(edges_from_sparse_matrix(file_name, ix_converter, sep))


def load_matrix(file_name, size=None, sep=None, ix_converter=None):
    """
    Load Hi-C matrix. Must be in sparse matrix format.
    Region information corresponding to the rows in the Hi-C file
    must be given in a separate BED file.
    """

    if size is None:
        raise ValueError("Must provide matrix size!")

    m = np.zeros((size, size))
    for source, sink, weight in edges_from_sparse_matrix(file_name, ix_converter=ix_converter, sep=sep):
        m[source, sink] = weight
        m[sink, source] = weight

    return m


def load_pairs(file_name, sep=None):
    """
    Load reference and query map from bedpe file.
    Expected colums:
        chrm_ref, start_ref, end_ref, chrm_qry, start_qry, end_qry, id,
        score (not used), strand_ref, strand_qry
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
            pairs[ID] = (cref, cqry)
    return pairs


def split_by_chrom(matrix_file, regions, ix2reg, full_matrix_ix_converter=None,
                   sep=None, sampleID=str(), work_dir='./'):
    """
    Split sparse matrix into single chromosome matrix files
    + the corresponding region bed files.
    """

    # sparse matrices
    chrom_handles = {}
    ix_converters = {}
    chrom_ixs = {}
    with _open(matrix_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            s_line = line.rstrip()
            fields = s_line.split(sep)
            ix1, ix2 = (int(v) for v in fields[:2])
            if ix1 > ix2:
                ix1, ix2 = ix2, ix1
            (c1, s1, e1), (c2, s2, e2) = (ix2reg[int(v)] for v in [ix1, ix2])
            # allow only intra-chromosomal contacts
            if c1 != c2:
                continue
            if c1 not in chrom_handles:
                chrom_handles[c1] = open(os.path.join(
                    work_dir, '{}_{}.sparse'.format(sampleID, c1)), 'w')
                chrom_ixs[c1] = set()
            chrom_handles[c1].write(line)
            chrom_ixs[c1].update([ix1, ix2])

    # region.bed files
    chrom_sizes = {}
    chrom_regions = {}
    for region in regions:
        if region.chromosome not in chrom_regions:
            if region.chromosome in chrom_handles:
                chrom_handles[region.chromosome].close()
            chrom_regions[region.chromosome] = [region, region]
            continue
        cs, ce = chrom_regions[region.chromosome]
        if region.ix < cs.ix:
            chrom_regions[region.chromosome][0] = region
        if region.ix > ce.ix:
            chrom_regions[region.chromosome][1] = region

    for chrom, (cs, ce) in chrom_regions.items():
        size = ce.ix - cs.ix + 1
        # map ixs to split matrix size
        if full_matrix_ix_converter is not None:
            rev_full_m_ix = {v: k for k, v in full_matrix_ix_converter.items()}
            cix_s, cix_e = [int(rev_full_m_ix[e]) for e in (cs.ix, ce.ix)]
        else:
            cix_s, cix_e = cs.ix, ce.ix
        ix_converters[chrom] = {
            str(k): v
            for k, v in zip(
                range(cix_s, cix_e + 1),
                range(size))}
        chrom_sizes[chrom] = size

        with open(
            os.path.join(
                work_dir, '{}_{}.regions.bed'.format(
                    sampleID, chrom)), 'w') as f:
            for ix, cix in ix_converters[chrom].items():
                line = '\t'.join(
                    [str(e) for e in ix2reg[int(ix)]] + [str(cix)])
                f.write(line + '\n')

    return chrom_sizes, ix_converters
