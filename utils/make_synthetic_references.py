#!/usr/bin/env python

import numpy as np
import knight_matrix_balancing as kmb
import random
import sys
import argparse
from copy import deepcopy


def main():
    """
    Create args.n random matrices of dimension args.b x args.b,
    normalize with Knight-Ruiz balancing.

    Store the matrices, boundary positions, loop positions and intensities
    in .npy objects in args.outdir.
    """
    n = args.n
    test_reference_matrices = np.empty((n, args.b, args.b))
    test_reference_boundaries = []
    test_reference_loops = []
    for i in range(n):
        print('Making reference #{0} / {1}'.format(i + 1, n))
        x = make_decay_matrix(bins=args.b)
        b, l = generate_random_features(x)
        x = add_features(x, b, l)
        test_reference_matrices[i] = x
        test_reference_boundaries.append(b)
        test_reference_loops.append(l)
    np.save(
        '{}/{}_test_reference_matrices_{}bins.npy'.format(
            args.outdir, n, args.b),
        test_reference_matrices)
    np.save(
        '{}/{}_test_reference_boundaries_{}bins.npy'.format(
            args.outdir, n, args.b),
        np.array(test_reference_boundaries))
    np.save(
        '{}/{}_test_reference_loops_{}bins.npy'.format(
            args.outdir, n, args.b),
        np.array(test_reference_loops))


def knight_ruiz_norm(m):
    a, b = kmb.correct_matrix(m)
    return a


def kth_diag_indices(n, k):
    """
    Return indices of bins k steps away from the diagonal.
    from:
        http://stackoverflow.com/questions/10925671/numpy-k-th-diagonal-indices
    """

    rows, cols = np.diag_indices(n)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols


def add_features(blank, borders, loops):
    m = deepcopy(blank)
    big_brother_bins = borders[0][-1]
    scaling_factor = np.shape(m)[0] / big_brother_bins
    level_intensities = [75, 100, 150]
    for level, border_list in enumerate(borders):
        for i in range(len(border_list) - 1):
            start = int(border_list[i] * scaling_factor)
            end = int(border_list[i + 1] * scaling_factor)
            m[start:end, start:end] += level_intensities[level]
            loop = loops[level][i]
            if loop is not None:
                rad_rel, intensity = loop
                rad = int(max(((end - start) * rad_rel), 1))
                anchor_a_s = start - rad
                anchor_a_e = start + rad
                anchor_b_s = end - rad
                anchor_b_e = end + rad
                m[anchor_a_s:anchor_a_e, anchor_b_s:anchor_b_e] += intensity

    for i in range(np.shape(m)[0]):
        for j in range(i, np.shape(m)[1]):
            m[j][i] = m[i][j]

    return m


def make_decay_matrix(bins=1000):
    def power_decay(distances):
        res = np.power(distances, -0.85)
        return res
    x = np.full((bins, bins), 0, dtype=np.float_)
    # 1000 bin matrix as reference.
    # n bin matrix will take the first n entries for the decay.
    distances = np.linspace(0.0001, 1, 1000)
    decay = power_decay(distances)
    for i in range(bins):
        x[kth_diag_indices(bins, i)] = decay[i]
    # make x symetric
    for i in range(np.shape(x)[0]):
        for j in range(i, np.shape(x)[1]):
            x[j][i] = x[i][j]
    x = x.clip(min=1)
    return x


def generate_random_features(m, rounds=3, initial_intensity=75,
                             round_intensity_increase=50,
                             first_round_loops=True,
                             borders=None, sizes=None, loops=None):

    from scipy.stats import truncnorm

    def get_truncated_normal(mean, sd, low, upp):
        return truncnorm(
            (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

    lvl_2_TAD_sizes = {
        '2': get_truncated_normal(mean=150, sd=100, low=100, upp=300),
        '1': get_truncated_normal(mean=20, sd=15, low=10, upp=50),
        '0': get_truncated_normal(mean=10, sd=5, low=5, upp=20)
    }

    def add_tads_in_interval(istart, iend):
        new_borders = []
        new_loops = []
        curr_b = istart
        for i in range(iend):
            start = curr_b
            end = start + int(lvl_2_TAD_sizes[str(rounds)].rvs())
            if end >= iend:
                end = iend
            new_borders.append(int(curr_b))
            if first_round_loops and random.choice(loop_chances):
                rad_rel = random.choice(loop_radius_rel)
                intensity = random.choice(loop_intensity_rel) * constant
                new_loops.append((float(rad_rel), float(intensity)))
            else:
                new_loops.append(None)
            curr_b = end
            if end == iend:
                new_borders.append(curr_b)
                new_loops.append(None)
                break
        return new_borders, new_loops

    # loop stuff hardcoded for now, but whateeeevs
    loop_radius_rel = [0.05, 0.07, 0.1]
    loop_chances = [1, 0, 0]
    loop_intensity_rel = [0.9, 1.4, 2.2]
    constant = initial_intensity
    bins = m.shape[0]
    rounds -= 1
    if not borders:
        loops = []
        margin = bins / 100 * 10
        lo, hi = int(0 + margin), int(bins - margin)
        istart_2 = int(random.choice(range(lo, hi)))
        istart_1 = bins - istart_2
        iend = bins
        new_borders_2, new_loops_2 = add_tads_in_interval(istart_2, iend)
        new_borders_1, new_loops_1 = add_tads_in_interval(istart_1, iend)
        new_borders_1 = [abs(bins - b) for b in new_borders_1][::-1]
        borders = [new_borders_1[:-1] + new_borders_2]
        loops = [new_loops_1[:-1] + new_loops_2]
    else:
        last_level = borders[-1]
        new_level = []
        new_loops = []
        for i in range(len(last_level) - 1):
            istart = last_level[i]
            iend = last_level[i+1]
            cnv, cnl = add_tads_in_interval(istart, iend)
            new_level.extend(cnv)
            new_loops.extend(cnl)
        borders.append(new_level)
        loops.append(new_loops)
    if rounds > 0:
        initial_intensity += round_intensity_increase
        res = generate_random_features(
            m, rounds, initial_intensity,
            round_intensity_increase, first_round_loops=True, borders=borders,
            sizes=sizes, loops=loops)
        return res
    else:
        return borders, loops


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


if __name__ == '__main__':
    parser = MyParser(
        description='''
        Produce a set of normalized reference matrices
        for testing CHESS.''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'outdir',
        type=str,
        help='''Write output files to this directory.''')
    parser.add_argument(
        '-n',
        type=int,
        default=100,
        help='''Number of reference matrices to create.''')
    parser.add_argument(
        '-b',
        type=int,
        default=1000,
        help='''Size of reference matrices in bins.''')
    args = parser.parse_args()
    main()
