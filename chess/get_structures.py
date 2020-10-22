#!/usr/bin/env python

from os import path
import logging
import numpy as np
from scipy import ndimage as ndi
from skimage import restoration
from skimage.morphology import square, closing
import skimage.filters as filters
from skimage.measure import label, regionprops
from numpy import inf
from scipy.ndimage import zoom
import warnings
from future.utils import string_types
from .helpers import (
    load_regions, load_sparse_matrix, sub_matrix_from_edges_dict
    )
warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)


def get_sparse_matrix(regions, sparse_matrix):
    '''Get all the sparse HiC matrix from the bed file'''
    ix_converter = None
    if isinstance(regions, string_types):
        regions, ix_converter, _ = load_regions(regions)

    if isinstance(sparse_matrix, string_types):
        sparse_matrix = load_sparse_matrix(
            sparse_matrix, ix_converter=ix_converter)

    return regions, sparse_matrix


def clipped_zoom(img, zoom_factor, **kwargs):
    h, w = img.shape[:2]
    zoom_tuple = (zoom_factor,) * 2 + (1,) * (img.ndim - 2)
    if zoom_factor < 1:
        zh = int(np.round(h * zoom_factor))
        zw = int(np.round(w * zoom_factor))
        top = (h - zh) // 2
        left = (w - zw) // 2
        out = np.zeros_like(img)
        out[top:top+zh, left:left+zw] = zoom(img, zoom_tuple, **kwargs)
    elif zoom_factor > 1:
        zh = int(np.round(h / zoom_factor))
        zw = int(np.round(w / zoom_factor))
        top = (h - zh) // 2
        left = (w - zw) // 2
        out = zoom(img[top:top+zh, left:left+zw], zoom_tuple, **kwargs)
        trim_top = ((out.shape[0] - h) // 2)
        trim_left = ((out.shape[1] - w) // 2)
        out = out[trim_top:trim_top+h, trim_left:trim_left+w]
    else:
        out = img
    return out


def get_info_feature(labels, submatrix, outfile, position, area, reg):
    for region in regionprops(labels):
        minr, minc, maxr, maxc = region.bbox
        if region.area > area:
            bx = (minc, maxc, maxc, minc, minc)
            by = (minr, minr, maxr, maxr, minr)
            y_min, y_max = np.min(by), np.max(by)
            x_min, x_max = np.min(bx), np.max(bx)
            submat = submatrix[y_min:y_max, x_min:x_max]
            flat_mat = list(submat.flatten())
            flat_mat.insert(0, reg)
            flat_mat.insert(1, position)
            flat_mat.insert(2, x_min)
            flat_mat.insert(3, x_max)
            flat_mat.insert(4, y_min)
            flat_mat.insert(5, y_max)
            outfile.write(','.join(map(str, flat_mat)) + '\n')
            position += 1
    return position


def extract_structures(
    reference_edges,
    reference_regions,
    query_edges,
    query_regions,
    pairs,
    output,
    windowsize,
    sigma_spatial,
    size_medianfilter,
    closing_square
):
    pos_query = 0
    pos_reference = 0

    w1 = open(path.join(output, 'gained_features.tsv'), '+w')
    w2 = open(path.join(output, 'lost_features.tsv'), '+w')

    for pair_ix, reference_region, query_region in pairs:
        reference, ref_rs = sub_matrix_from_edges_dict(
            reference_edges,
            reference_regions,
            reference_region,
            default_weight=0.)
        query, qry_rs = sub_matrix_from_edges_dict(
            query_edges,
            query_regions,
            query_region,
            default_weight=0.)
        area = int((5 * np.shape(query)[0]) / 100)
        size = np.shape(query)[0]
        or_matrix = np.log(np.divide(query, reference))
        where_are_NaNs = np.isnan(or_matrix)
        or_matrix[where_are_NaNs] = 0.
        or_matrix[or_matrix == -inf] = 0.
        or_matrix[or_matrix == inf] = 0.
        std = np.std(or_matrix)
        # thres1, thres2 = mean + 0.5 * std, mean - 0.5 * std
        # positive = np.where(or_matrix > thres1, or_matrix, 0.)
        # negative = np.where(or_matrix < thres2, or_matrix, 0.)
        positive = np.where(or_matrix > (0.5 * std), or_matrix, 0.)
        negative = np.abs(np.where(or_matrix < -(0.5 * std), or_matrix, 0.))

        # denoise
        denoise_positive = restoration.denoise_bilateral(
            positive,
            sigma_color=np.mean(positive),
            win_size=windowsize,
            sigma_spatial=sigma_spatial,
            bins=size,
            multichannel=False)
        denoise_negative = restoration.denoise_bilateral(
            negative,
            sigma_color=np.mean(negative),
            win_size=windowsize,
            sigma_spatial=sigma_spatial,
            bins=size,
            multichannel=False)
        # smooth
        filter_positive = ndi.median_filter(
            denoise_positive, size_medianfilter)
        filter_negative = ndi.median_filter(
            denoise_negative, size_medianfilter)

        # binarise
        if np.all(filter_positive == 0.):
            threshold_pos = positive
        else:
            filter1 = filters.threshold_otsu(filter_positive, nbins=size)
            threshold_pos = filter_positive > filter1

        if np.all(filter_negative == 0.):
            threshold_neg = negative
        else:
            filter2 = filters.threshold_otsu(filter_negative, nbins=size)
            threshold_neg = filter_negative > filter2

        # zoom and rotate 45 degrees
        zm1 = clipped_zoom(threshold_pos, 0.7)
        rot1 = ndi.rotate(zm1, 45, reshape=False)
        zm2 = clipped_zoom(threshold_neg, zoom_factor=0.7)
        rot2 = ndi.rotate(zm2, 45, reshape=False)

        # Close morphology
        img1 = closing(rot1, square(closing_square))
        label_x1 = label(img1)
        img2 = closing(rot2, square(closing_square))
        label_x2 = label(img2)

        zm1_reference = clipped_zoom(reference, 0.7)
        rot_reference = ndi.rotate(zm1_reference, 45, reshape=False)
        zm1_query = clipped_zoom(query, 0.7)
        rot_query = ndi.rotate(zm1_query, 45, reshape=False)

        # get output (file with label and submatrices)
        pos_query = get_info_feature(
            label_x1, rot_query, w1, pos_query, area, pair_ix)
        pos_reference = get_info_feature(
            label_x2, rot_reference, w2, pos_reference, area, pair_ix)
