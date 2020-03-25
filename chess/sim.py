#!/usr/bin/env python

import logging
import numpy as np

from skimage.transform import resize
from skimage.metrics import structural_similarity
from numpy.lib.stride_tricks import as_strided
from .helpers import sub_matrix_from_edges_dict
import uuid

logger = logging.getLogger(__name__)


def SN(M1, M2, size=7):
    """Compute average of local signal to noise ratios in M1 - M2.

    Slide a window of dimensions size^2 over M1 - M2, compute
    abs(mean(window)) / var(window) in each and take mean of all windows.
    :param M1: Numpy Matrix.
    :param M2: Numpy Matrix.
    :param size: edge length of sliding window (must be odd), defaults to 7
    :returns: float, the average of the local signal to noise ratios.
    """
    M = M1 - M2
    uM = filter_uniform_filter(M, size=size)
    vM = filter_uniform_filter(M, size=size, func=np.nanvar)
    SN = np.nanmean(np.abs(uM) / vM)

    return float(SN)


def filter_uniform_filter(M, size=7, no_data_val=None,
                          func=np.nanmean):
    """Apply func to sliding window on input matrix.

    Adaption from:
    www.programcreek.com/python/example/93943/scipy.ndimage.uniform_filter
    This is equivalent to scipy.ndimage.uniform_filter, but can handle nan's,
    and can use numpy nanmean/median/max/min functions.
    :param M: input data, numpy matrix.
    :param size: edge length of sliding window (must be odd), defaults to 7
    :type size: number, optional
    :param no_data_val: value in matrix that is treated as no data value,
                        defaults to None
    :param func: [description], defaults to np.nanmean
    :returns: func of the matrix A filtered by a uniform kernel of size=size
    """

    assert size % 2 == 1, 'Please supply an odd size'
    assert M.dtype == np.dtype('float64'), 'Input must have dtype(float64)'
    rows, cols = M.shape

    padded_A = np.empty(shape=(rows + size-1,
                               cols + size-1),
                        dtype=M.dtype)
    padded_A[:] = np.nan
    rows_pad, cols_pad = padded_A.shape

    if no_data_val:
        mask = M == no_data_val
        M[mask] = np.nan

    upleft_data_start = int((size - 1) / 2)
    rl, cl = rows_pad - upleft_data_start, cols_pad - upleft_data_start
    padded_A[upleft_data_start: rl, upleft_data_start: cl] = M.copy()

    n, m = M.shape
    strided_data = as_strided(padded_A, (n, m, size, size),
                              padded_A.strides+padded_A.strides)
    strided_data = strided_data.copy().reshape((n, m, size**2))

    return func(strided_data, axis=2)


def find_masked_rows(m, masking_value=0):
    """Find rows in m filled with masking value.

    :param m: Numpy matrix to be searched for masked rows.
    :param masking_value: Value that corresponds to missing contacts,
                          defaults to 0
    :returns: Indices of masked rows.
    """
    s = np.sum(m, 0)
    n_bins = m.shape[0]
    cutoff = n_bins * masking_value
    idxs = np.where(s == cutoff)[0]

    return idxs


def remove_rows(m, target_rows):
    """Remove target bins from m.

    Remove rows and columns with indices indicated in target rows form m.
    :param m: Numpy matrix.
    :param target_rows: List of indices of bins to be removed from m.
    :returns: Truncated numpy matrix.
    """
    idxs = target_rows
    m_removed = np.delete(m, idxs, 0)
    m_removed = np.delete(m_removed, idxs, 1)

    return m_removed


def comparison_worker(input_queue, output_queue,
                      reference_edges, reference_regions,
                      query_edges, query_regions,
                      min_bins=20,
                      keep_unmappable_bins=False,
                      absolute_window_size=None,
                      relative_window_size=1.,
                      mappability_cutoff=0.1,
                      default_value=1):
    worker_id = uuid.uuid4()
    try:
        while True:
            data = input_queue.get(block=True)
            if data is None:
                break
            pair, compute_sn = data

            pair_ix, ssim, sn = compare_pair(pair,
                                             reference_edges, reference_regions,
                                             query_edges, query_regions,
                                             min_bins=min_bins,
                                             keep_unmappable_bins=keep_unmappable_bins,
                                             absolute_window_size=absolute_window_size,
                                             relative_window_size=relative_window_size,
                                             mappability_cutoff=mappability_cutoff,
                                             default_value=default_value,
                                             compute_sn=compute_sn)
            output_queue.put((pair_ix, ssim, sn))
    except Exception as e:
        logger.error("Worker {} encountered a problem: {}".format(worker_id, e))
        output_queue.put(e)

    logger.info("Worker {} exited gracefully".format(worker_id))


def chunk_comparison_worker(input_queue, output_queue,
                            reference_edges, reference_regions,
                            query_edges, query_regions,
                            min_bins=20,
                            keep_unmappable_bins=False,
                            absolute_window_size=None,
                            relative_window_size=1.,
                            mappability_cutoff=0.1,
                            default_value=1):
    worker_id = uuid.uuid4()
    previous_reference = [None, (None, None)]
    previous_query = [None, (None, None)]
    try:
        while True:
            data = input_queue.get(block=True)
            if data is None:
                break
            chunk, compute_sn = data

            results = []
            for pair in chunk:
                pair_ix, reference_region, query_region = pair
                try:
                    if previous_reference[0] is not None and reference_region == previous_reference[0]:
                        reference, ref_rs = previous_reference[1]
                    else:
                        reference, ref_rs = sub_matrix_from_edges_dict(reference_edges,
                                                                       reference_regions,
                                                                       reference_region,
                                                                       default_weight=default_value)
                        previous_reference = [reference_region, (reference, ref_rs)]

                    if previous_query[0] is not None and query_region == previous_query[0]:
                        query, qry_rs = previous_query[1]
                    else:
                        query, qry_rs = sub_matrix_from_edges_dict(query_edges,
                                                                   query_regions,
                                                                   query_region,
                                                                   default_weight=default_value)
                        previous_query = [query_region, (query, qry_rs)]
                except ValueError:
                    results.append([pair_ix, np.nan, np.nan])
                    continue

                pair_ix, ssim, sn = compare_matrix_pair(reference, reference_region,
                                                        query, query_region,
                                                        pair_ix=pair_ix,
                                                        min_bins=min_bins,
                                                        keep_unmappable_bins=keep_unmappable_bins,
                                                        absolute_window_size=absolute_window_size,
                                                        relative_window_size=relative_window_size,
                                                        mappability_cutoff=mappability_cutoff,
                                                        default_value=default_value,
                                                        compute_sn=compute_sn)
                results.append([pair_ix, ssim, sn])
            output_queue.put(results)
    except Exception as e:
        logger.error("Worker {} encountered a problem: {}".format(worker_id, e))
        output_queue.put(e)

    logger.debug("Worker {} exited gracefully".format(worker_id))


def compare_pair(pair, reference_edges, reference_regions,
                 query_edges, query_regions,
                 min_bins=20,
                 keep_unmappable_bins=False,
                 absolute_window_size=None,
                 relative_window_size=1.,
                 mappability_cutoff=0.1,
                 default_value=1,
                 compute_sn=True
                 ):
    """Run comparison of given Hi-C matrices without any background model.

    Compare given submatrices of the full Hi-C matrices in sampleID2hic.
    All submatrices (specified in pairs) have to be of the same size.
    The assignment of reference and query samples is arbitrary in this case,
    but has to be consistent.
    :param edges: edges dictionary
    :param regions: regions intervaltree

    :param worker_ID: String or int used as ID for the process in which this
                      function is run.
    :param pairs: Dictionary indicating the submatrix pairs that
                  will be compared.
                  Must have form:
                  {
                    pair_id: (reference_region, query_region)
                  },
                  where reference_region and query_region
                  are :class:`~GenomicRegion` objects.
    :param min_bins: Minimum size of matrix to be considered for comparison,
                     defaults to 20
    :param keep_unmappable_bins: If True, disables removal of unmappable bins
                                 from matrices prior to comparison,
                                 defaults to False
    :param absolute_windowsize: Absolute window size used in the
                                structural similarity function. Overwrites
                                relative_windowsize, defaults to None
    :param relative_windowsize: Window size relative to the size
                                of the compared matrices,
                                used in the structural similarity function,
                                defaults to 1.
    :param mappability_cutoff: Maximum fraction of unmappable bins allowed
                               in a matrix to be considered for comparison,
                               defaults to 0.1

    :returns: results dictionary of form:
              {
                pair_id: (similarity_score, signal_to_noise_value)
              }
    """
    pair_ix, reference_region, query_region = pair
    try:
        reference, ref_rs = sub_matrix_from_edges_dict(reference_edges, reference_regions,
                                                       reference_region, default_weight=default_value)
        query, qry_rs = sub_matrix_from_edges_dict(query_edges, query_regions,
                                                   query_region, default_weight=default_value)
    except ValueError:
        return pair_ix, np.nan, np.nan

    return compare_matrix_pair(reference, reference_region,
                               query, query_region,
                               pair_ix=pair_ix,
                               min_bins=min_bins,
                               keep_unmappable_bins=keep_unmappable_bins,
                               absolute_window_size=absolute_window_size,
                               relative_window_size=relative_window_size,
                               mappability_cutoff=mappability_cutoff,
                               default_value=default_value,
                               compute_sn=compute_sn)


def compare_matrix_pair(reference, reference_region,
                        query, query_region, pair_ix=None,
                        min_bins=20,
                        keep_unmappable_bins=False,
                        absolute_window_size=None,
                        relative_window_size=1.,
                        mappability_cutoff=0.1,
                        default_value=1,
                        compute_sn=True):
    # check size criteria
    if query.shape[0] < min_bins or reference.shape[0] < min_bins:
        return pair_ix, np.nan, np.nan

    # resize matrices if needed
    size_factor = len(reference) / len(query)
    if size_factor < 1:
        reference = resize(reference, np.shape(query), order=0,
                           mode='reflect').astype(np.float64)
    elif size_factor > 1:
        query = resize(query, np.shape(reference), order=0,
                       mode='reflect').astype(np.float64)

    # check if either matrix needs to be flipped
    if reference_region.strand != query_region.strand:
        query = np.fliplr(np.flipud(query))

    # check content in unmappable bins in input matrices
    idxs = []
    init_filter_pass = True
    if mappability_cutoff > 0:
        ms = [reference, query]
        for m in ms:
            m = np.array(m)
            idx = find_masked_rows(m, masking_value=default_value)
            idxs.append(idx)
            if len(idx) / np.shape(m)[0] > mappability_cutoff:
                init_filter_pass = False
    if not init_filter_pass:
        return pair_ix, np.nan, np.nan

    # remove masked bins
    if not keep_unmappable_bins:
        target_rows = np.union1d(*idxs)
        reference = remove_rows(reference, target_rows)
        query = remove_rows(query, target_rows)

    # compare!
    if absolute_window_size:
        window_size = absolute_window_size
    else:
        window_size = int(np.shape(reference)[0] * relative_window_size)
    if window_size % 2 == 0:
        window_size -= 1

    # prevent from going lower than three.
    # breaks at ws = 1 because of 0 division, and 2 is not odd.
    window_size = max(window_size, 3)

    # actual comparison
    curr_ssim = structural_similarity(reference, query, win_size=window_size)
    if compute_sn:
        curr_sn = SN(reference, query)
    else:
        curr_sn = None

    return pair_ix, curr_ssim, curr_sn


def compare_structures(reference_edges, reference_regions,
                       query_edges, query_regions,
                       pairs, worker_ID,
                       min_bins=20,
                       keep_unmappable_bins=False,
                       absolute_window_size=None,
                       relative_window_size=1.,
                       mappability_cutoff=0.1,
                       default_value=1):
    """Run comparison of given Hi-C matrices without any background model.

    Compare given submatrices of the full Hi-C matrices in sampleID2hic.
    All submatrices (specified in pairs) have to be of the same size.
    The assignment of reference and query samples is arbitrary in this case,
    but has to be consistent.
    :param edges: edges dictionary
    :param regions: regions intervaltree

    :param worker_ID: String or int used as ID for the process in which this
                      function is run.
    :param pairs: Dictionary indicating the submatrix pairs that
                  will be compared.
                  Must have form:
                  {
                    pair_id: (reference_region, query_region)
                  },
                  where reference_region and query_region
                  are :class:`~GenomicRegion` objects.
    :param min_bins: Minimum size of matrix to be considered for comparison,
                     defaults to 20
    :param keep_unmappable_bins: If True, disables removal of unmappable bins
                                 from matrices prior to comparison,
                                 defaults to False
    :param absolute_windowsize: Absolute window size used in the
                                structural similarity function. Overwrites
                                relative_windowsize, defaults to None
    :param relative_windowsize: Window size relative to the size
                                of the compared matrices,
                                used in the structural similarity function,
                                defaults to 1.
    :param mappability_cutoff: Maximum fraction of unmappable bins allowed
                               in a matrix to be considered for comparison,
                               defaults to 0.1
    :returns: results dictionary of form:
              {
                pair_id: (similarity_score, signal_to_noise_value)
              }
    """

    results = {}  # store the ssim in here
    no_pass = set()

    for pair in pairs:
        pair_ix, ssim, sn = compare_pair(pair,
                                         reference_edges, reference_regions,
                                         query_edges, query_regions,
                                         min_bins=min_bins,
                                         keep_unmappable_bins=keep_unmappable_bins,
                                         absolute_window_size=absolute_window_size,
                                         relative_window_size=relative_window_size,
                                         mappability_cutoff=mappability_cutoff,
                                         default_value=default_value)
        if np.isnan(ssim):
            no_pass.add(pair_ix)
            results[pair_ix] = (np.nan, np.nan)
        else:
            results[pair_ix] = (ssim, sn)

    logger.debug('[WORKER #{}]: Done!'.format(worker_ID))

    return results
