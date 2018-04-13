#!/usr/bin/env python

import numpy as np
import knight_matrix_balancing as kmb
import sys
import argparse
import pandas as pd
import random
from skimage.transform import resize
from copy import deepcopy
from scipy.stats import zscore
from skimage.measure import compare_ssim
from scipy.stats import truncnorm


def main():
    ref_m = np.load(args.reference_matrices)
    ref_b = np.load(args.reference_boundaries)
    ref_l = np.load(args.reference_loops)
    test_results = run_test(
        ref_m, ref_b, ref_l, size_factor=args.size_factor,
        query_noise=args.query_noise,
        ref_noise=args.ref_noise,
        test_pool=args.test_pool, pool_size=args.pool_size,
        query_size=args.qs)
    base = args.out
    test_results.to_csv(base, sep='\t', index=False)


def run_test(ref_m, ref_b, ref_l, size_factor=1,
             query_noise=0, ref_noise=0, pool_size=1000, test_pool=None,
             query_size=None):

    # 0: norm, 1: OE, 2: log2(OE)
    mtypes = {0: 'norm', 1: 'OE', 2: 'log2(OE)'}
    ssim = []
    z_scores = []
    p_val = []
    win_sizes = []
    matrix_ids = []
    matrix_types = []
    if args.store_background:
        background_scores = []
    n_bins = ref_m[0].shape[0]
    if query_size is None:
        query_size = int(n_bins * size_factor)
    ref_depth = scale_depth(n_bins, args.depth)
    query_depth = scale_depth(query_size, args.depth)
    query_decay = make_decay_matrix(bins=query_size)
    ref_decay = make_decay_matrix(bins=n_bins)

    if test_pool is None:
        test_pool = {'0': np.empty((pool_size, query_size, query_size)),
                     '1': np.empty((pool_size, query_size, query_size)),
                     '2': np.empty((pool_size, query_size, query_size))}
        for i in range(pool_size):
            print('filling test pool: {0} / {1}'.format(i + 1, pool_size))
            r = deepcopy(ref_decay)
            x = deepcopy(query_decay)
            br, lr = generate_random_features(r)
            x = add_features(x, br, lr)
            x = adjust_depth(x, query_depth)
            x = add_noise(x, ref_noise)
            p = add_noise(x, query_noise)
            p_n = knight_ruiz_norm(p)
            test_pool['0'][i] = p_n
            test_pool['1'][i] = make_OE(p_n)
            test_pool['2'][i] = np.log2(make_OE(p_n))
    else:
        test_pool = np.load(test_pool)[()]
    n_references = len(ref_b)

    for i in range(n_references):
        print('current test round: {0} / {1}'.format(i + 1, n_references))
        curr_id = i
        x, b, l = ref_m[i], ref_b[i], ref_l[i]
        x = add_noise(
            adjust_depth(
                x, ref_depth),
            ref_noise)
        sb = add_noise(
            adjust_depth(
                add_features(
                    query_decay, b, l),
                query_depth),
            ref_noise)
        sb = add_noise(sb, query_noise)
        x_n = knight_ruiz_norm(x)
        sb_n = knight_ruiz_norm(sb)
        oe_x = make_OE(x_n)
        oe_sb = make_OE(sb_n)
        oe_x_log2 = np.log2(make_OE(x_n))
        oe_sb_log2 = np.log2(make_OE(sb_n))
        rq_pairs = [(x_n, sb_n), (oe_x, oe_sb), (oe_x_log2, oe_sb_log2)]
        ran_with_7_bins = False

        for rel_window in [0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1]:
            win_size = max([int(n_bins * rel_window), 7])
            if win_size == 7:
                if ran_with_7_bins:
                    continue
                else:
                    ran_with_7_bins = True
            for i, (reference, query) in enumerate(rq_pairs):
                truth_ssim = compare_matrices(reference, query, rel_window)
                control_ssim = [
                    compare_matrices(
                        reference, t, rel_window)
                    for t in test_pool[str(i)]
                ]
                truth_p = np.sum(
                    [1 for v in control_ssim if v >= truth_ssim]) / pool_size
                all_ssim = control_ssim + [truth_ssim]
                zs = zscore(all_ssim)
                z = zs[-1]
                if args.store_background:
                    background_scores[str(i)].extend(control_ssim)
                ssim.append(truth_ssim)
                z_scores.append(z)
                p_val.append(truth_p)
                win_sizes.append(win_size)
                matrix_ids.append(curr_id)
                matrix_types.append(mtypes[i])

    if args.store_background:
        base = args.out
        np.save(
            base + '_BACKGROUND.npy', background_scores)
    results_frame = pd.DataFrame({
        'p': p_val,
        'ssim': ssim,
        'ws': win_sizes,
        'z': z_scores,
        'matrix_id': matrix_ids,
        'matrix_type': matrix_types
    })

    return results_frame


def random_subset(iterator, K):
    """
    Reservior sampling,
    from: https://stackoverflow.com/questions/2612648/reservoir-sampling
    """
    result = []
    N = 0
    for item in iterator:
        N += 1
        if len(result) < K:
            result.append(item)
        else:
            s = int(random.random() * N)
            if s < K:
                result[s] = item

    return result


def make_symetric(m):
    for i in range(np.shape(m)[0]):
        for j in range(i, np.shape(m)[1]):
            m[j][i] = m[i][j]

    return m


def scale_depth(size, depth):
    """
    Return number of reads m should have
    given a depth for a 100^2 (1 mb) matrix.
    """
    size_factor = size / 100
    scaled_depth = depth * size_factor

    return scaled_depth


def adjust_depth(m, depth):
    """
    Adjust number of reads in matrix to specified depth.
    (Random removal of reads.)
    """
    reads = {}
    rows = m.tolist()
    total_reads = 0
    max_pos = m.shape[0]
    # extract read counts
    for i, row in enumerate(rows):
        for j, n_reads in enumerate(row[i:]):
            curr_dot = (i, j + i)
            reads[curr_dot] = n_reads
            total_reads += n_reads
    # this returns non-unique values
    rnd_pos = random_subset(
        (k for k, v in reads.items() for x in range(0, int(v))),
        int(total_reads - depth))
    # remove random reads
    for k in rnd_pos:
        reads[k] -= 1
    # convert reads back to matrix
    mn = np.full((max_pos, max_pos), 0, dtype=np.float_)
    for (i, j), n_reads in reads.items():
        mn[i, j] = n_reads
    mn = make_symetric(mn)
    mn = mn.clip(min=1)

    return mn


def add_noise(m, noise_level, float_adjust=True):
    reads = {}
    rows = m.tolist()
    total_reads = 0
    max_pos = m.shape[0]
    # extract read counts
    for i, row in enumerate(rows):
        for j, n_reads in enumerate(row[i:]):
            curr_dot = (i, j + i)
            reads[curr_dot] = n_reads
            total_reads += n_reads
    # this returns non-unique values
    rnd_pos = random_subset(
        (k for k, v in reads.items() for x in range(0, int(v))),
        int(total_reads * noise_level))
    # remove random reads
    for k in rnd_pos:
        reads[k] -= 1
    # add random reads
    keys = list(reads.keys())
    rnd_pos_add = np.random.randint(
        0, len(keys), int(total_reads * noise_level))
    for i in rnd_pos_add:
        reads[keys[i]] += (1 + random.random() * float_adjust)
    if float_adjust:
        rnd_pos_rm = np.random.randint(
            0, len(keys), int(total_reads * noise_level))
        for i in rnd_pos_rm:
            reads[keys[i]] -= random.random()
    # convert reads back to matrix
    mn = np.full((max_pos, max_pos), 0, dtype=np.float_)
    out_total_reads = 0
    for (i, j), n_reads in reads.items():
        mn[i, j] = n_reads
        out_total_reads += n_reads
    mn = make_symetric(mn)
    mn = mn.clip(min=1)

    return mn


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


def compare_matrices(a, b, rel_window):
    a_size = np.shape(a)[0]
    b_size = np.shape(b)[0]
    if a_size > b_size:
        smaller_m = b
        bigger_m = a
    else:
        smaller_m = a
        bigger_m = b
    if args.scaling_mode == 1:
        smaller_m_resized = resize(
            smaller_m, np.shape(bigger_m),
            order=0, mode='reflect').astype(np.float64)
        a, b = bigger_m, smaller_m_resized
    elif args.scaling_mode == 0:
        bigger_m_resized = resize(
            bigger_m, np.shape(smaller_m),
            order=0, mode='reflect').astype(np.float64)
        a, b = bigger_m_resized, smaller_m
    else:
        raise ValueError(
            ('Invalid value for -scale_mode: {}.'
             'Must be 1 or 0.'))
    winsize = int(np.shape(a)[0] * rel_window)
    if winsize % 2 == 0:
        winsize -= 1
    winsize = max([winsize, 7])
    ssim = compare_ssim(a, b, win_size=winsize)

    return ssim


def reformat_test_results(df):
    ssim = []
    p = []
    matrix_type = []
    win_size = []
    z_scores = []
    for index, row in df.iterrows():
        for t in ['oe_log2', 'oe', 'norm']:
            matrix_type.append(t)
            ssim.append(row['ssim_{0}'.format(t)])
            p.append(row['p_{0}'.format(t)])
            win_size.append(row['ws_{0}'.format(t)])
            z_scores.append(row['z_{0}'.format(t)])
    res = pd.DataFrame({
        'ssim': ssim,
        'p': p,
        'matrix_type': matrix_type,
        'win_size': win_size,
        'z_score': z_scores
    })

    return res


def make_OE(m):
    n = np.shape(m)[0]
    expected = np.full((n, n), 0, dtype=np.float_)
    for i in range(n):
        observed_diag = m[kth_diag_indices(n, i)]
        expected[kth_diag_indices(n, i)] = (
            np.sum(observed_diag) / np.size(observed_diag))
    for i in range(np.shape(expected)[0]):
        for j in range(i, np.shape(expected)[1]):
            expected[j][i] = expected[i][j]

    return m / expected


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
        Test the CHESS on a set of specified synthetic Hi-C matrices.

        Create query matrices and a pool of matrices forming the background
        model according to input parameters.
        Return .csv with raw ssim score, p and z-value, window size,
        matrix type and id (corresponding to the input references indices).''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'reference_matrices',
        type=str,
        help='''.npy file (3D array) with previously generated
        reference matrices.''')
    parser.add_argument(
        'reference_boundaries',
        type=str,
        help='''.npy file with boundary locations for each
        of the reference matrices.''')
    parser.add_argument(
        'reference_loops',
        type=str,
        help='''.npy file with loop locations and intensities for each
        of the reference matrices.''')
    parser.add_argument(
        'out',
        type=str,
        help='''Path to output file.''')
    parser.add_argument(
        '-s', '--size-factor',
        type=float,
        default=1.0,
        help='''Scaling factor for the query matrices
        (s = size(query) / size(reference)).''')
    parser.add_argument(
        '-qs',
        type=int,
        default=None,
        help='''Query size in bins. Overwrites -s, --size_factor.''')
    parser.add_argument(
        '-qn', '--query-noise',
        type=float,
        default=0,
        help='''Noise level for the query matrices, with 0 <= n <= 1.
        Fraction of simulated reads that will be replaced by randomly
        positioned reads. Will be placed on top of reference noise.''')
    parser.add_argument(
        '-rn', '--ref-noise',
        type=float,
        default=0,
        help='''Noise level for the reference matrices, with 0 <= n <= 1.
        Fraction of simulated reads that will be replaced by randomly
        positioned reads.''')
    parser.add_argument(
        '--depth',
        type=int,
        default=1500000,
        help='''Number of reads per 100x100 (1 mb) matrix.''')
    parser.add_argument(
        '-p', '--test-pool',
        type=str,
        default=None,
        help='''Path to a precomputed set of test_matrices. Make sure
        they have the correct scaling and noise level.
        Will produce new test pool if not provided.''')
    parser.add_argument(
        '--store-background',
        action='store_true',
        default=False,
        help='''Save scores calculated for the
        comparisons to the decoy pool.''')
    parser.add_argument(
        '--pool-size',
        type=int,
        default=1000,
        help='''Number of random query matrices that will be generated
        for the test pool.''')
    parser.add_argument(
        '--scaling-mode',
        type=bool,
        default=1,
        help='''Scale smaller matrix up to size of larger matrix (1), or
        scale larger matrix down (0).''')
    args = parser.parse_args()
    main()
