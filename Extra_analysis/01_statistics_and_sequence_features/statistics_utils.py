"""Paired rank tests and distribution summaries, without trajectory simulation."""
import itertools

import numpy as np
from scipy import stats


def median_interval(values, confidence=.95):
    """Distribution-free, pointwise order-statistic CI for a population median.

    Coverage uses Binomial(n, .5); ties make the continuous-case bound conservative.
    None denotes an unbounded endpoint when n is too small for a finite interval.
    This interval is for the median, not the signed-rank location parameter.
    """
    x = np.sort(np.asarray(values, dtype=float))
    if not len(x) or not np.isfinite(x).all() or not 0 < confidence < 1:
        raise ValueError('Invalid median interval input')
    tail = (1 - confidence) / 2
    k = int(stats.binom.ppf(tail, len(x), .5))
    while k + 1 <= len(x) // 2 and stats.binom.cdf(k, len(x), .5) <= tail:
        k += 1
    coverage = 1 - 2 * stats.binom.cdf(k - 1, len(x), .5)
    return dict(median=float(np.median(x)), ci_low=float(x[k - 1]) if k else None,
                ci_high=float(x[len(x) - k]) if k else None,
                ci_confidence=confidence, ci_coverage_bound=float(coverage))


def holm_log_adjust(log_pvalues):
    """Holm step-down adjustment in natural-log space, retaining tiny p-values."""
    p = np.asarray(log_pvalues, dtype=float)
    if p.ndim != 1 or not len(p) or not np.isfinite(p).all() or np.any(p > 0):
        raise ValueError('Expected finite nonpositive log p-values')
    order = np.argsort(p, kind='stable')
    adjusted = np.minimum(0., np.maximum.accumulate(p[order] + np.log(np.arange(len(p), 0, -1))))
    result = np.empty_like(p)
    result[order] = adjusted
    return result


def p_fields(log_p, prefix='p'):
    p = float(np.exp(log_p))
    underflow = p == 0.
    return {prefix: p, prefix + '_log10': float(log_p / np.log(10.)),
            prefix + '_underflow': underflow,
            prefix + '_display': ('10^({:.6g})'.format(log_p / np.log(10.))
                                 if underflow else '{:.8g}'.format(p))}


def signed_rank(differences):
    """Two-sided Wilcoxon normal approximation, zero_method=wilcox, no correction.

    Signed ranks and the tie correction match SciPy's asymptotic implementation.
    Differences supplied by the runner are integer path-count sums, avoiding
    spurious floating-point tie breaking on the n-trajectory mean lattice.
    Positive effect means group1 > group2. All-zero differences return p=1.
    """
    d = np.asarray(differences)
    if d.ndim != 1 or not len(d) or not np.isfinite(d).all():
        raise ValueError('Invalid paired differences')
    nonzero = d[d != 0]
    n = len(nonzero)
    if not n:
        return dict(n_nonzero=0, n_zero=len(d), W_min=0., z=0., rank_biserial=0.,
                    **p_fields(0.))
    absolute = np.abs(nonzero)
    ranks = stats.rankdata(absolute, method='average')
    positive = float(ranks[nonzero > 0].sum())
    negative = float(ranks[nonzero < 0].sum())
    _, ties = np.unique(absolute, return_counts=True)
    variance = n * (n + 1.) * (2 * n + 1.) / 24. - np.sum(ties.astype(float)**3 - ties) / 48.
    z = (positive - n * (n + 1.) / 4.) / np.sqrt(variance)
    log_p = min(0., float(np.log(2.) + stats.norm.logsf(abs(z))))
    return dict(n_nonzero=n, n_zero=len(d) - n, W_min=min(positive, negative), z=float(z),
                rank_biserial=(positive - negative) / (positive + negative), **p_fields(log_p))


def to_count_sums(matrix, n_trajectories):
    scaled = np.asarray(matrix, dtype=float) * n_trajectories
    if (not isinstance(n_trajectories, int) or n_trajectories <= 0
            or not np.isfinite(scaled).all() or np.any(scaled < 0)
            or np.any(scaled > np.iinfo(np.int64).max / 2)
            or not np.allclose(scaled, np.rint(scaled), rtol=0, atol=1e-6)):
        raise ValueError('Means are not on the declared complete-trajectory count lattice')
    return np.rint(scaled).astype(np.int64)


def distribution(values):
    x = np.asarray(values, dtype=float)
    q1, median, q3 = np.quantile(x, [.25, .5, .75])
    return dict(n=len(x), mean=float(x.mean()), sd=float(x.std(ddof=1)) if len(x) > 1 else 0.,
                median=float(median), q1=float(q1), q3=float(q3), iqr=float(q3-q1),
                minimum=float(x.min()), maximum=float(x.max()))


def paired_analysis(matrix, names, n_trajectories, progress=None):
    counts = to_count_sums(matrix, n_trajectories)
    n, k = counts.shape
    if n < 2 or k < 3 or len(names) != k or len(set(names)) != k:
        raise ValueError('Friedman analysis requires >=2 blocks and >=3 unique signatures')
    if np.all(counts == counts[:, :1]):
        statistic, pvalue = 0., 1.
    else:
        test = stats.friedmanchisquare(*counts.T)
        statistic, pvalue = float(test.statistic), float(test.pvalue)
    if not np.isfinite(statistic):
        raise ValueError('Undefined Friedman statistic')
    omnibus = dict(test='Friedman, tie-corrected chi-square approximation', n_blocks=n,
                   signatures=k, statistic=statistic, df=k-1, p=pvalue,
                   p_underflow=(pvalue == 0.),
                   p_display=('< smallest representable positive float (underflow)' if pvalue == 0.
                              else '{:.8g}'.format(pvalue)),
                   kendall_W=statistic / (n * (k-1)),
                   asymptotic_size_guidance_met=(n > 10 and k > 6))
    summaries = [dict(Signature=name, **distribution(matrix[:, i])) for i, name in enumerate(names)]
    rows = []
    for i, j in itertools.combinations(range(k), 2):
        d = counts[:, i] - counts[:, j]
        test = signed_rank(d)
        ci = median_interval(d / n_trajectories)
        rows.append(dict(group1=names[i], group2=names[j], n_paired=n,
                         median_group1=float(np.median(matrix[:, i])),
                         median_group2=float(np.median(matrix[:, j])),
                         paired_median_difference=ci.pop('median'), **ci, **test))
        if progress and len(rows) % 500 == 0:
            progress(len(rows), k * (k-1) // 2)
    adjusted = holm_log_adjust([row['p_log10'] * np.log(10.) for row in rows])
    for row, value in zip(rows, adjusted):
        row.update(p_fields(float(value), 'p_holm'))
        row['holm_family_size'] = len(rows)
        row['reject_holm_0_05'] = bool(value < np.log(.05))
    return omnibus, summaries, rows


def independent_rank_effect(x, y):
    """Descriptive rank-biserial effect; no test or biological sampling inference."""
    x, y = np.asarray(x), np.asarray(y)
    rank_sum = stats.rankdata(np.concatenate([x, y]))[:len(x)].sum()
    u = rank_sum - len(x) * (len(x) + 1.) / 2.
    return float(2 * u / (len(x) * len(y)) - 1.)
