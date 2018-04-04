#!/usr/bin/env python

import os
import logging
import numpy as np
from math import ceil
from scipy.stats import zscore
from skimage.transform import resize
from skimage.measure import compare_ssim
from .helpers import load_matrix, load_regions, sub_matrix_regions, GenomicRegion

logger = logging.getLogger(__name__)


def find_masked_rows(m, masking_value=1):
    s = np.sum(m, 0)
    n_bins = m.shape[0]
    cutoff = n_bins * masking_value
    idxs = np.where(s == cutoff)[0]

    return idxs


def remove_rows(m, target_rows):
    idxs = target_rows
    m_removed = np.delete(m, idxs, 0)
    m_removed = np.delete(m_removed, idxs, 1)

    return m_removed


def compare_structures_genome_scan(reference_ID, query_ID, sampleID2hic, worker_ID,
                                   pairs, min_bins=20, work_dir='./',
                                   keep_unmappable_bins=False, absolute_windowsize=None,
                                   relative_windowsize=1., mappability_cutoff=0.1,
                                   limit_background=False):

    def load_chrom(sample, chrom):
        size, ix = (sampleID2hic[sample][k][chrom]
                    for k in ['sizes', 'ix'])
        m = load_matrix(
            os.path.join(work_dir, '{}_{}.sparse'.format(sample, chrom)),
            ix_converter=ix, size=size)
        r, c_ix_conv, c_ix2reg = load_regions(
            os.path.join(work_dir, '{}_{}.regions.bed'.format(
                sample, chrom)))

        return m, r

    # init
    results = {}  # store the ssim in here
    query_size = None
    no_pass = set()
    refchrm = None
    qrychrm = None

    # transform pairs to tuple representation and sort to avoid
    # uneccessary chromosome loading
    tpairs = [(k, v[0], v[1]) for k, v in pairs.items()]
    tpairs = sorted(tpairs, key=lambda p: p[1].chromosome)
    logger.info("[WORKER #{0}]: Running comparisons.".format(worker_ID))
    # load input matrices and get rounded positions
    len_ids = len(pairs)
    for curr_pos, (ID, refreg, qryreg) in enumerate(tpairs):
        curr_pos = curr_pos + 1
        logger.info("[WORKER #{0}]: Matrix {1} / {2} ".format(
            worker_ID, curr_pos, len_ids))
        if refchrm != refreg.chromosome:
            refchrm_m, refchrm_r = load_chrom(reference_ID, refreg.chromosome)
            refchrm = refreg.chromosome
        if qrychrm != qryreg.chromosome:
            qrychrm_m, qrychrm_r = load_chrom(query_ID, qryreg.chromosome)
            qrychrm = qryreg.chromosome

        reference, ref_rs = sub_matrix_regions(refchrm_m, refchrm_r, refreg)
        query, qry_rs = sub_matrix_regions(qrychrm_m, qrychrm_r, qryreg)

        # binround_qry_start = qry_rs[0].start
        # binround_qry_end = qry_rs[-1].end
        # curr_query_size = binround_qry_end - binround_qry_start
        curr_query_size = len(qry_rs)
        if curr_query_size != query_size and query_size is not None:
            logger.debug("Current region: {}".format(qryreg))
            raise ValueError(
                'Varying query sizes in genome scan!')
        query_size = curr_query_size

        # check size criteria
        if query.shape[0] < min_bins or reference.shape[0] < min_bins:
            no_pass.add(ID)
            continue

        # check content in unmappable bins in input matrices
        idxs = []
        init_filter_pass = True
        if mappability_cutoff > 0:
            ms = [reference, query]
            for m in ms:
                m = np.array(m)
                idx = find_masked_rows(m, masking_value=0)
                idxs.append(idx)
                if len(idx) / np.shape(m)[0] > mappability_cutoff:
                    init_filter_pass = False
        if not init_filter_pass:
            no_pass.add(ID)
            continue

        # remove masked bins
        if not keep_unmappable_bins:
            target_rows = np.union1d(*idxs)
            reference = remove_rows(reference, target_rows)
            query = remove_rows(query, target_rows)

        # compare!
        if absolute_windowsize:
            winsize = absolute_windowsize
        else:
            winsize = int(np.shape(reference)[0] * relative_windowsize)
        if winsize % 2 == 0:
            winsize -= 1
        # prevent from going lower than three.
        # breaks at ws = 1 because of 0 division, and 2 is not odd.
        if winsize < 3:
            winsize = 3
        curr_ssim = compare_ssim(
            reference, query, win_size=winsize)
        results[ID] = curr_ssim
    logger.info(
        ('[WORKER #{0}]:'
         'Done! The following pairs were filtered out: {1}').format(
         worker_ID, no_pass))

    return results


def compare_structures_sliding_window(reference_ID, query_ID, sampleID2hic,
                                      worker_ID, pairs, min_bins=20, work_dir='./',
                                      keep_unmappable_bins=False, absolute_windowsize=None,
                                      relative_windowsize=1., mappability_cutoff=0.1,
                                      limit_background=False):
    print('<-------', work_dir)
    def load_chrom(sample, chrom):
        size, ix = (sampleID2hic[sample][k][chrom]
                    for k in ['sizes', 'ix'])
        m = load_matrix(
            os.path.join(work_dir, '{}_{}.sparse'.format(sample, chrom)),
            ix_converter=ix, size=size)
        r, c_ix_conv, c_ix2reg = load_regions(
            os.path.join(work_dir, '{}_{}.regions.bed'.format(
                sample, chrom)))

        return m, r

    # init
    references = {
        0: {},  # native orientation
        1: {}   # flipped
    }
    reference_masks = {
        0: {},  # native orientation
        1: {}   # flipped
    }
    results = {}                # store the full comparisons in here
    rounded_queries = {}        # store the rounded region pairs in here
    query_sizes = {}            # store the query sizes in here
    size_factors = {}           # store the size factors in here
    no_pass = set()
    chromosomes_to_parse = set()
    resolution = None
    refchrm = None
    qrychrm = None

    # transform pairs to tuple representation and sort to avoid
    # uneccessary chromosome loading
    tpairs = [(k, v[0], v[1]) for k, v in pairs.items()]
    tpairs = sorted(tpairs, key=lambda p: p[1].chromosome)

    logger.info("[WORKER #{0}]: Loading reference matrices into memory".format(
        worker_ID))
    # load input matrices and get rounded positions
    len_ids = len(pairs)
    for curr_pos, (ID, refreg, qryreg) in enumerate(tpairs):
        curr_pos = curr_pos + 1
        logger.info("[WORKER #{0}]: Matrix {1} / {2} ".format(
            worker_ID, curr_pos, len_ids))
        if refchrm != refreg.chromosome:
            refchrm_m, refchrm_r = load_chrom(reference_ID, refreg.chromosome)
            refchrm = refreg.chromosome
        if qrychrm != qryreg.chromosome:
            qrychrm_m, qrychrm_r = load_chrom(query_ID, qryreg.chromosome)
            qrychrm = qryreg.chromosome
        reference, ref_rs = sub_matrix_regions(refchrm_m, refchrm_r, refreg)
        query, qry_rs = sub_matrix_regions(qrychrm_m, qrychrm_r, qryreg)
        binround_qry_start = qry_rs[0].start
        binround_qry_end = qry_rs[-1].end
        syn_size = binround_qry_end - binround_qry_start

        if resolution is None:
            resolution = refchrm_r[0].end - refchrm_r[0].start + 2

        # check size criteria
        if query.shape[0] < min_bins or reference.shape[0] < min_bins:
            no_pass.add(ID)
            continue

        # resize matrices if needed
        size_factor = len(reference) / len(query)
        if size_factor < 1:
            reference = resize(
                reference, np.shape(query),
                order=0, mode='reflect').astype(np.float64)
        elif size_factor > 1:
            query = resize(
                query, np.shape(reference),
                order=0, mode='reflect').astype(np.float64)

        # check content in unmappable bins in input matrices
        idxs = []
        init_filter_pass = True
        if mappability_cutoff > 0:
            ms = [reference, query]
            for m in ms:
                m = np.array(m)
                idx = find_masked_rows(m, masking_value=0)
                idxs.append(idx)
                if len(idx) / np.shape(m)[0] > mappability_cutoff:
                    init_filter_pass = False
        if not init_filter_pass:
            no_pass.add(ID)
            continue

        # get masked values for reference
        if not keep_unmappable_bins:
            reference_mask = idxs[0]
            reference_mask_flipped = np.sort(
                np.shape(reference)[0] - 1 - reference_mask)
            for k, v in enumerate([reference_mask, reference_mask_flipped]):
                reference_masks[k][ID] = v

        # save all the stuff needed for the comparisons
        chromosomes_to_parse.add(qryreg.chromosome)
        references[0][ID] = reference
        references[1][ID] = np.fliplr(np.flipud(reference))
        rounded_queries[ID] = ':'.join(
            (qryreg.chromosome, str(binround_qry_start),
             str(binround_qry_end), qryreg.strand))
        query_sizes[ID] = syn_size
        size_factors[ID] = size_factor
        results[ID] = {}

    logger.info(
        ('[WORKER #{0}]:'
         ' The following pairs were filtered out: {1}.'
         ' Continue with {2} matrices').format(
            worker_ID, no_pass, len(references[0])))

    c_l = {
        k: v * resolution for k, v in sampleID2hic[query_ID]['sizes'].items()}

    # run sliding window comparisons
    if not limit_background:
        chromosomes_to_parse = c_l.keys()
    for chrom in chromosomes_to_parse:
        logger.info(
            ('[WORKER #{0}]:'
             ' Loading new chromosome matrix into memory: {1}').format(
                worker_ID, chrom))
        qrychrm_m, qrychrm_r = load_chrom(query_ID, chrom)
        logger.info("[WORKER #{0}]: Done loading. Running comparisons".format(
            worker_ID))

        for region in qrychrm_r:
            for ID, reference_pos in references[0].items():

                if limit_background:
                    if pairs[ID][1].chromosome != chrom:
                        continue

                reference_neg = references[1][ID]
                query_bins = query_sizes[ID]

                # get query position
                win_chrm = region.chromosome
                win_start = int(region.start)
                win_end = int(win_start + query_bins)
                winreg = GenomicRegion(
                    chromosome=win_chrm, start=win_start, end=win_end)

                # skip if window overlaps end of chromosome
                if win_end > c_l[win_chrm]:
                    continue

                # get query matrix
                query, query_bin_windows = sub_matrix_regions(
                    qrychrm_m, qrychrm_r, winreg)

                # resize
                if size_factors[ID] > 1:
                    query = resize(
                        query, np.shape(reference_pos),
                        order=0, mode='reflect').astype(np.float64)
                if np.shape(reference_pos) != np.shape(query):
                    logger.debug('[WORKER #{0}] ERROR! Block ID: {1}'.format(
                        worker_ID, ID))
                    logger.debug(
                        '[WORKER #{0}] ERROR! Size factor: {1}'.format(
                            worker_ID, size_factors[ID]))
                    logger.debug(
                        '[WORKER #{0}] ERROR! Reference shape: {1}'.format(
                            worker_ID, np.shape(reference_pos)))
                    logger.debug(
                        '[WORKER #{0}] ERROR! Query shape: {1}'.format(
                            worker_ID, np.shape(query)))
                    logger.debug(
                        '[WORKER #{0}] ERROR! Curr window: {1}'.format(
                            worker_ID, str(winreg)))
                    raise RuntimeError(
                        'Reference and query sizes differ after resizing!.')

                # check mappability criteria
                if mappability_cutoff <= 0 and keep_unmappable_bins:
                    pass
                else:
                    idx = find_masked_rows(
                        np.array(query),
                        masking_value=0)
                    if len(idx) / np.shape(query)[0] > mappability_cutoff:
                        continue

                # compare
                for flip, reference in enumerate(
                        [reference_pos, reference_neg]):
                    # create entry in results
                    orient = '-' if flip else '+'
                    curr_win = ':'.join(
                                (win_chrm,
                                 str(win_start),
                                 str(win_end), orient))
                    results[ID][curr_win] = {}

                    # remove unmappable bins
                    if (mappability_cutoff <= 0
                            and keep_unmappable_bins):
                        trunc_reference = reference
                        trunc_query = query
                    else:
                        target_rows = np.union1d(
                            idx, reference_masks[flip][ID])
                        trunc_reference = remove_rows(reference, target_rows)
                        trunc_query = remove_rows(query, target_rows)

                    # compare matrices and store results
                    if absolute_windowsize:
                        winsize = absolute_windowsize
                    else:
                        winsize = int(
                            np.shape(trunc_reference)[0] * relative_windowsize)
                    if winsize % 2 == 0:
                        winsize -= 1
                    # prevent from going lower than three.
                    # breaks at ws = 1 because of 0 division, and 2 is not odd.
                    if winsize < 3:
                        winsize = 3
                    curr_ssim = compare_ssim(
                        trunc_reference, trunc_query, win_size=winsize)
                    results[ID][curr_win] = curr_ssim

    logger.info("[WORKER #{0}]: Done!".format(
        worker_ID))
    return results, rounded_queries


def post_process(raw_results, rounded_queries):
    rows = []
    for ID, results in raw_results.items():
        truth_pos = rounded_queries[ID]
        pair_comp = results[truth_pos]
        background = [v for k, v in results.items() if k != truth_pos]
        all_scores = background + [pair_comp]
        p = (np.sum([1 for s in all_scores if s >= pair_comp])
             / len(all_scores))
        all_zs = zscore(all_scores)
        z = all_zs[-1]

        rows.append({
            'ID': ID,
            'p': p,
            'ssim': pair_comp,
            'z-score': z
            })

    return rows


def post_process_simple(raw_results):
    rows = []
    IDs, scores = zip(
        *[(k, v)
          for k, v in raw_results.items()])
    zscores = zscore(scores)
    for p, ID in enumerate(IDs):
        ssim = scores[p]
        zssim = zscores[p]
        rows.append({
            'ID': ID,
            'ssim': ssim,
            'zssim': zssim
            })

    return rows


def cleanup(chromosome_basenames, converted_oe_files=[], reference_ID='REF',
            query_ID='QRY', work_dir='./'):

    def remove(path):
        try:
            os.remove(path)
        except OSError:
            pass

    for file in converted_oe_files:
        remove(file)
    for sampleID in reference_ID, query_ID:
        remove(os.path.join(work_dir, 'chrom_sizes_' + sampleID))
        remove(os.path.join(work_dir, 'ix_converters_' + sampleID))
        for basename in chromosome_basenames:
            for ext in ['.sparse', '.regions.bed']:
                remove(os.path.join(
                    work_dir, sampleID + '_' + basename + ext))


def distribute_workload(pairs, limit_background, p=1):
    if limit_background:
        logger.warning(
            ('[MAIN]: --limit-background set.'
             ' Will calculate backgrounds '
             'only on syntenic chromosome.'))
        pc_subsets = {}

        for ID, (refreg, qryreg) in pairs.items():
            chrm = qryreg.chromosome
            if chrm not in pc_subsets:
                pc_subsets[chrm] = []
            pc_subsets[chrm].append(ID)

        n_id_pc = {k: len(v) for k, v in pc_subsets.items()}
        n_ids = sum([v for k, v in n_id_pc.items()])
        workers_per_chromosome = {k: ceil(v / n_ids * p)
                                  for k, v in n_id_pc.items()}
        n_workers = sum([v for k, v in workers_per_chromosome.items()])
        pair_subsets = {}
        worker_IDs = iter(range(n_workers))

        for chrm, IDs in pc_subsets.items():
            share = int(n_id_pc[chrm] / workers_per_chromosome[chrm])
            overhang = set(
                range(n_id_pc[chrm] % workers_per_chromosome[chrm]))

            curr_start = 0
            for pos in range(workers_per_chromosome[chrm]):
                worker_ID = next(worker_IDs)
                start = curr_start
                if pos in overhang:
                    stop = start + share + 1
                else:
                    stop = start + share
                curr_start = stop
                if pos == workers_per_chromosome[chrm] - 1:
                    pair_share = IDs[start:]
                else:
                    pair_share = IDs[start:stop]
                pair_subsets[worker_ID] = {
                    ID: pairs[ID] for ID in pair_share
                }
    else:
        pair_subsets = {}
        worker_IDs = list(range(p))
        pair_ids = list(pairs.keys())
        share = int(len(pair_ids) / p)
        overhang = set(range(len(pair_ids) % p))
        curr_start = 0
        for worker_ID in worker_IDs:
            start = curr_start
            if worker_ID in overhang:
                stop = start + share + 1
            else:
                stop = start + share
            curr_start = stop
            if worker_ID == p - 1:
                pair_share = pair_ids[start:]
            else:
                pair_share = pair_ids[start:stop]
            pair_subsets[worker_ID] = {
                 ID: pairs[ID] for ID in pair_share
             }

    return pair_subsets

