"""Independent references and physical invariants of the new climate summaries."""
import numpy as np
import pytest
from scipy.stats import skew

from spea.climate_statistics import statistical_features, STATISTICS, SHAPE, SEASONAL, subset_borders


def test_uniform_climate_conventions():
    x, audit = statistical_features(np.full((1, 12), 10.), np.full((1, 12), 5.))
    assert audit["exact_constant_temperature_rows"] == audit["exact_constant_precipitation_rows"] == 1
    np.testing.assert_allclose(x[["precipitation_normalized_entropy", "precipitation_median_over_mean"]], 1)
    np.testing.assert_allclose(x[["wettest_circular_3month_fraction", "driest_circular_3month_fraction"]], .25)
    zero = [k for k in STATISTICS if k not in ["precipitation_normalized_entropy", "precipitation_median_over_mean", "wettest_circular_3month_fraction", "driest_circular_3month_fraction"]]
    np.testing.assert_allclose(x[zero], 0, atol=1e-7)


def test_quantiles_skew_and_calendar_wrap_against_independent_references():
    p = np.array([[40, 20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 60]], dtype=float)
    t = np.array([[-2, -1, 4, 9, 15, 20, 25, 24, 18, 10, 5, 0]], dtype=float)
    x, _ = statistical_features(p, t)
    assert x.temperature_skew.iloc[0] == pytest.approx(skew(t[0], bias=True), abs=1e-7)
    assert x.precipitation_skew.iloc[0] == pytest.approx(skew(p[0], bias=True), abs=1e-7)
    assert x.temperature_iqr_c.iloc[0] == pytest.approx(np.percentile(t[0], 75) - np.percentile(t[0], 25))
    windows = [sum(p[0, (start + step) % 12] for step in range(3))/sum(p[0]) for start in range(12)]
    assert x.wettest_circular_3month_fraction.iloc[0] == pytest.approx(max(windows)) == 1
    assert x.driest_circular_3month_fraction.iloc[0] == min(windows) == 0


def test_precipitation_scale_and_temperature_translation():
    p = np.arange(1, 25, dtype=float).reshape(2, 12)
    t = np.arange(-10, 14, dtype=float).reshape(2, 12)
    a, _ = statistical_features(p, t)
    b, _ = statistical_features(p * 7, t + 273.15)
    np.testing.assert_allclose(a, b, atol=2e-7)


def test_joint_rotation_preserves_shape_and_coupling_rotates_harmonics():
    p = np.array([[1, 3, 6, 11, 20, 34, 29, 18, 9, 4, 2, 1]], dtype=float)
    t = np.array([[-2, 0, 4, 9, 15, 22, 25, 24, 18, 11, 5, 0]], dtype=float)
    a, _ = statistical_features(p, t)
    b, _ = statistical_features(np.roll(p, 3, axis=1), np.roll(t, 3, axis=1))
    fixed = SHAPE + [SEASONAL[i] for i in [0, 1, 2, 7]]
    np.testing.assert_allclose(a[fixed], b[fixed], atol=2e-7)
    for v in ["precipitation", "temperature"]:
        suffix = "_c" if v == "temperature" else ""
        sin_key, cos_key = v+"_harmonic_sin"+suffix, v+"_harmonic_cos"+suffix
        assert b[sin_key].iloc[0] == pytest.approx(a[cos_key].iloc[0], abs=1e-7)
        assert b[cos_key].iloc[0] == pytest.approx(-a[sin_key].iloc[0], abs=1e-7)


def test_single_wet_month_and_bimodal_cancellation():
    t = np.full((1, 12), 5.)
    p = np.zeros((1, 12)); p[0, 0] = 120
    x, _ = statistical_features(p, t)
    assert x.precipitation_normalized_entropy.iloc[0] == 0
    assert x.precipitation_harmonic_cos.iloc[0] == 1
    p[0, 0] = p[0, 6] = 60
    x, _ = statistical_features(p, t)
    np.testing.assert_allclose(x[["precipitation_harmonic_sin", "precipitation_harmonic_cos"]], 0, atol=1e-7)
    assert x.precipitation_normalized_entropy.iloc[0] < 1


def test_small_real_variance_is_not_replaced_with_zero():
    p = np.arange(1, 13, dtype=float)[None, :]
    t = np.full((1, 12), 1.); t[0, 0] += 1e-10
    x, audit = statistical_features(p, t)
    assert audit["exact_constant_temperature_rows"] == 0
    assert x.temperature_skew.iloc[0] > 2


@pytest.mark.parametrize("p,t", [(np.zeros((1, 12)), np.ones((1, 12))),
                               (-np.ones((1, 12)), np.ones((1, 12))),
                               (np.ones((1, 12)), np.full((1, 12), np.nan)),
                               (np.ones((1, 11)), np.ones((1, 11)))])
def test_invalid_rows_raise_without_dropping(p, t):
    with pytest.raises(ValueError):
        statistical_features(p, t)


def test_crossed_feature_border_remapping(tmp_path):
    source = tmp_path / "union.tsv"
    source.write_text("0\t1.25\n1\t3.5\n2\t0.1\n2\t0.3\n3\t0.5\n")
    dest = tmp_path / "chosen.tsv"
    subset_borders(source, dest, ["a", "b", "c", "d"], ["d", "b"])
    assert dest.read_text() == "0\t0.5\n1\t3.5\n"
    with pytest.raises(FileExistsError):
        subset_borders(source, dest, ["a", "b", "c", "d"], ["d", "b"])
