from collections import defaultdict
import numpy as np
from future.utils import string_types
from .helpers import load_regions, load_sparse_matrix


def _chromosome_bins(regions):
    chr_bins = {}
    for r in regions:
        if chr_bins.get(r.chromosome) is None:
            chr_bins[r.chromosome] = [r.ix, r.ix + 1]
        else:
            chr_bins[r.chromosome][1] = r.ix + 1
    return chr_bins


def possible_contacts(regions, mappability):
    chromosomes = []
    current_chromosome = None
    for r in regions:
        if r.chromosome != current_chromosome:
            chromosomes.append(r.chromosome)
        current_chromosome = r.chromosome

    cb = _chromosome_bins(regions)
    chromosome_max_distance = defaultdict(int)
    max_distance = 0
    chromosome_subtractions = dict()
    for chromosome in chromosomes:
        start, stop = cb[chromosome]
        max_distance = max(max_distance, stop - start)
        chromosome_max_distance[chromosome] = max(chromosome_max_distance[chromosome], stop - start)
        chromosome_subtractions[chromosome] = np.zeros(stop - start, dtype='int32')

    chromosome_mappable = defaultdict(int)
    chromosome_unmappable = defaultdict(set)
    for i, mappable in enumerate(mappability):
        chromosome = regions[i].chromosome
        if not mappable:  # unmappable
            s = chromosome_subtractions[chromosome]
            o = cb[chromosome][0]
            ix = i - o
            # horizontal
            s[0: len(s) - ix] += 1
            # vertical
            for j in range(1, ix + 1):
                if ix - j not in chromosome_unmappable[chromosome]:
                    s[j] += 1
            chromosome_unmappable[chromosome].add(ix)
        else:
            chromosome_mappable[chromosome] += 1

    inter_total = 0
    intra_total = [0] * max_distance
    chromosome_intra_total = dict()
    for chromosome, d in chromosome_max_distance.items():
        chromosome_intra_total[chromosome] = [0] * d

    for i, chromosome in enumerate(chromosomes):
        start, stop = cb[chromosome]
        count = stop - start

        # intra-chromosomal
        s = chromosome_subtractions[chromosomes[i]]
        for distance in range(0, count):
            intra_total[distance] += count - distance - s[distance]
            chromosome_intra_total[chromosome][distance] += count - distance - s[distance]

        # inter-chromosomal
        for j in range(i + 1, len(chromosomes)):
            count_mappable = chromosome_mappable[chromosomes[i]]
            count2_mappable = chromosome_mappable[chromosomes[j]]
            inter_total += count_mappable * count2_mappable

    return intra_total, chromosome_intra_total, inter_total


def expected_values(regions, sparse_matrix, selected_chromosome=None):
    # get all the bins of the different chromosomes
    chromosome_bins = _chromosome_bins(regions)
    chromosome_dict = defaultdict(list)

    chromosome_max_distance = defaultdict(int)
    max_distance = 0
    for chromosome, (start, stop) in chromosome_bins.items():
        max_distance = max(max_distance, stop - start)
        chromosome_max_distance[chromosome] = max(chromosome_max_distance[chromosome], stop - start)

        for i in range(start, stop):
            chromosome_dict[i] = chromosome

    chromosome_intra_sums = dict()
    chromosome_intra_expected = dict()
    for chromosome, d in chromosome_max_distance.items():
        chromosome_intra_sums[chromosome] = [0.0] * d
        chromosome_intra_expected[chromosome] = [0.0] * d

    # get the sums of edges at any given distance
    marginals = [0.0] * len(regions)
    inter_sums = 0.0
    intra_sums = [0.0] * max_distance
    for source, sink, weight in sparse_matrix:
        if not np.isfinite(weight):
            continue

        source_chromosome = chromosome_dict[source]
        sink_chromosome = chromosome_dict[sink]

        marginals[source] += weight
        marginals[sink] += weight

        if sink_chromosome != source_chromosome:
            inter_sums += weight
        else:
            distance = sink - source
            intra_sums[distance] += weight
            chromosome_intra_sums[source_chromosome][distance] += weight

    mappability = [m > 0 for m in marginals]
    intra_total, chromosome_intra_total, inter_total = possible_contacts(regions, mappability)

    # expected values
    inter_expected = 0 if inter_total == 0 else inter_sums / inter_total

    intra_expected = [0.0] * max_distance
    bin_size = regions[0].end - regions[0].start
    distances = []
    for d in range(max_distance):
        distances.append(bin_size * d)

        # whole genome
        count = intra_total[d]
        if count > 0:
            intra_expected[d] = intra_sums[d] / count

    # chromosomes
    for chromosome in chromosome_intra_expected:
        for d in range(chromosome_max_distance[chromosome]):
            chromosome_count = chromosome_intra_total[chromosome][d]
            if chromosome_count > 0:
                chromosome_intra_expected[chromosome][d] = chromosome_intra_sums[chromosome][d] / chromosome_count

    if selected_chromosome is not None:
        return chromosome_intra_expected[selected_chromosome]

    return intra_expected, chromosome_intra_expected, inter_expected


def _oe_generator(sparse_matrix, ix_to_chromosome, intra_expected, inter_expected, log=False):
    for source, sink, weight in sparse_matrix:
        source_chromosome = ix_to_chromosome[source]
        sink_chromosome = ix_to_chromosome[sink]
        if source_chromosome == sink_chromosome:
            e = intra_expected[source_chromosome][sink - source]
        else:
            e = inter_expected

        try:
            oe = weight/e
        except ZeroDivisionError:
            oe = 1

        if log:
            yield source, sink, np.log2(oe)

        yield source, sink, oe


def observed_expected(regions, sparse_matrix):
    ix_converter = None
    if isinstance(regions, string_types):
        regions, ix_converter, _ = load_regions(regions)

    if isinstance(sparse_matrix, string_types):
        sparse_matrix = load_sparse_matrix(sparse_matrix, ix_converter=ix_converter)

    _, chromosome_intra_expected, inter_expected = expected_values(regions, sparse_matrix)

    ix_to_chromosome = {i: r.chromosome for i, r in enumerate(regions)}

    return regions, list(_oe_generator(sparse_matrix, ix_to_chromosome, chromosome_intra_expected, inter_expected))
