import logging
import warnings
import numpy as np
logger = logging.getLogger(__name__)


def remove_sparse_rows(m, cutoff=None):
    s = np.sum(m, 0)

    if cutoff is None:
        cutoff = min(s)

    idxs = np.where(s <= cutoff)[0]
    m_removed = np.delete(m, idxs, 0)
    m_removed = np.delete(m_removed, idxs, 1)

    return m_removed, idxs


def restore_sparse_rows(m, idx_sets, rows=None):
    abs_idx = []
    for idxs in reversed(idx_sets):
        for i in sorted(idxs):
            shift = 0
            for j in sorted(abs_idx):
                if j + shift < i:
                    shift += 1
            abs_idx.append(i - shift)
    abs_idx.sort()
    a = np.insert(m, abs_idx, 0, axis=0)
    if len(m.shape) > 1:
        a = np.insert(a, abs_idx, 0, axis=1)
    return a


def correct_matrix(m, max_attempts=50, restore_coverage=False):
    # remove zero-sum rows
    removed_rows = []
    m_nonzero, ixs = remove_sparse_rows(m, cutoff=0)
    removed_rows.append(ixs)

    has_errors = True
    iterations = 0
    x = None
    while has_errors:
        has_errors = False

        try:
            x = get_bias_vector(m_nonzero)
        except ValueError as e:
            logger.debug("Matrix balancing failed (this can happen!), \
                          removing sparsest rows to try again. Error: \
                          %s" % str(e))
            m_nonzero, ixs = remove_sparse_rows(m_nonzero)
            removed_rows.append(ixs)
            has_errors = True

        iterations += 1
        if iterations > max_attempts:
            raise RuntimeError("Exceeded maximum attempts (%d)" % max_attempts)

    if restore_coverage:
        x = x*np.sqrt(np.sum(m_nonzero)/m_nonzero.shape[0])

    logger.debug("Applying bias vector")
    m_nonzero = x*m_nonzero*x[:, np.newaxis]

    logger.debug(removed_rows)
    logger.debug("Restoring {} sets ({} total) sparse rows.".format(
        len(removed_rows), sum(len(x) for x in removed_rows)))
    # restore zero rows
    m_nonzero = restore_sparse_rows(m_nonzero, removed_rows)
    x = restore_sparse_rows(x, removed_rows)

    return m_nonzero, x


def get_bias_vector(A, x0=None, tol=1e-06, delta=0.1, Delta=3, fl=0, high_precision=False, outer_limit=300):
    logger.debug("Starting matrix balancing")

    with warnings.catch_warnings():
        warnings.filterwarnings('error')

        try:
            # basic variables
            # n=size_(A,1)
            if not isinstance(A, np.ndarray):
                try:
                    if high_precision:
                        A = np.array(A, dtype=np.float128)
                    else:
                        A = np.array(A)
                except AttributeError:
                    A = np.array(A)
            n = A.shape[0]
            # e=ones_(n,1)
            try:
                if high_precision:
                    e = np.ones(n, dtype=np.float128)
                else:
                    e = np.ones(n)
            except AttributeError:
                e = np.ones(n)

            if not x0:
                try:
                    if high_precision:
                        x0 = np.ones(n, np.float128)
                    else:
                        x0 = np.ones(n)
                except AttributeError:
                    x0 = np.ones(n)
            else:
                try:
                    if high_precision:
                        x0 = np.array(x0, np.float128)
                    else:
                        x0 = np.array(x0)
                except AttributeError:
                    x0 = np.array(x0)
            res = np.array([])
            g = 0.9
            etamax = 0.1
            eta = etamax
            stop_tol = tol * 0.5
            x = x0.copy()
            rt = tol ** 2
            # v=x.dot((A * x))
            v = x*A.dot(x)
            rk = 1 - v
            # rho_km1=rk.T * rk
            rho_km1 = rk.T.dot(rk)
            rout = rho_km1
            rho_km2 = rho_km1
            rold = rout
            MVP = 0
            i = 0

            n_iterations_outer = 0
            while rout > rt:
                n_iterations_outer += 1

                if n_iterations_outer > outer_limit:
                    raise ValueError("Number of iterations has exceeded the limit (%d)." % outer_limit)

                i += 1
                k = 0
                y = e.copy()
                innertol = max(eta ** 2 * rout, rt)
                n_iterations_inner = 0
                while rho_km1 > innertol:
                    n_iterations_inner += 1

                    k += 1
                    if k == 1:
                        try:
                            Z = rk / v
                        except Warning:
                            raise ValueError("v=0; Remove zero or sparse rows")
                        p = Z.copy()
                        rho_km1 = rk.T.dot(Z)
                    else:
                        beta = rho_km1 / rho_km2
                        p = Z + beta * p
                    # w = x.*(A*(x.*p)) + v.*p;
                    w = x*A.dot(x*p) + v*p
                    alpha = rho_km1 / p.T.dot(w)
                    ap = alpha * p
                    ynew = y + ap
                    if min(ynew) <= delta:
                        if delta == 0:
                            break
                        ind = np.where(ap < 0)[0]
                        # gamma = min((delta  - y(ind))./ap(ind));
                        gamma = min((delta-y[ind])/ap[ind])
                        y = y + gamma * ap
                        break
                    if max(ynew) >= Delta:
                        ind = np.where(ynew > Delta)[0]
                        gamma = min((Delta-y[ind])/ap[ind])
                        y = y + gamma * ap
                        break
                    y = ynew.copy()
                    rk = rk - alpha * w
                    rho_km2 = rho_km1.copy()
                    Z = rk / v
                    rho_km1 = rk.T.dot(Z)

                try:
                    x = x*y
                except Warning:
                    raise ValueError("Value in x or y too small to represent numerically. Try removing sparse rows")

                v = x*A.dot(x)
                rk = 1 - v
                rho_km1 = rk.T.dot(rk)
                rout = rho_km1.copy()
                MVP = MVP + k + 1
                rat = rout / rold
                rold = rout.copy()
                res_norm = np.sqrt(rout)
                eta_o = eta
                eta = g * rat
                if g * eta_o ** 2 > 0.1:
                    eta = max(eta, g * eta_o ** 2)
                eta = max(min(eta, etamax), stop_tol / res_norm)
                if fl == 1:
                    res = np.array([[res], [res_norm]])

            logger.debug("Matrix-vector products = %d\n" % MVP)
            logger.debug("Outer iterations: %d" % n_iterations_outer)
        except Warning as e:
            logger.error(str(e))
            raise ValueError("Generic catch all warnings")

    return x
