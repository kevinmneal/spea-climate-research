"""Scientific invariants of the seasonal representations; synthetic data only."""
import numpy as np
import pytest
from spea.seasonal_model import monthly_features, ANNUAL, VARIABILITY, MONTHLY, prefix_borders


def test_annual_and_variability_invariant_to_calendar_permutation():
    p = np.array([[1,2,4,8,16,32,64,32,16,8,4,2]], dtype=float)
    t = np.array([[-2,-1,3,8,14,20,25,24,18,10,4,0]], dtype=float)
    first = monthly_features(p,t)
    shifted = monthly_features(np.roll(p,4,axis=1),np.roll(t,4,axis=1))
    np.testing.assert_array_equal(first[VARIABILITY], shifted[VARIABILITY])
    assert not np.array_equal(first[MONTHLY], shifted[MONTHLY])


def test_physical_shift_and_scale_preserve_monthly_shape():
    p = np.arange(1,25,dtype=float).reshape(2,12)
    t = np.arange(-8,16,dtype=float).reshape(2,12)
    first = monthly_features(p,t)
    changed = monthly_features(p*4,t+5)
    np.testing.assert_allclose(first.iloc[:,2:],changed.iloc[:,2:],atol=1e-7)
    np.testing.assert_allclose(changed[ANNUAL[0]]-first[ANNUAL[0]],5)
    np.testing.assert_allclose(np.expm1(changed[ANNUAL[1]]),4*np.expm1(first[ANNUAL[1]]),rtol=1e-6)


def test_constant_climate_and_fraction_mass():
    features = monthly_features(np.full((2,12),10.),np.full((2,12),5.))
    np.testing.assert_array_equal(features[VARIABILITY[2:]],0)
    np.testing.assert_array_equal(features.iloc[:,4:16],0)
    np.testing.assert_allclose(features.iloc[:,16:].sum(axis=1),1.,atol=2e-7)
    np.testing.assert_allclose(features[ANNUAL[1]],np.log1p(120))


@pytest.mark.parametrize('p,t',[(np.zeros((1,12)),np.ones((1,12))),(-np.ones((1,12)),np.ones((1,12))),
                               (np.ones((1,12)),np.full((1,12),np.nan)),(np.ones((1,11)),np.ones((1,11)))])
def test_invalid_climate_stops_instead_of_dropping_cells(p,t):
    with pytest.raises(ValueError):
        monthly_features(p,t)


def test_shared_prefix_keeps_exact_cutoffs(tmp_path):
    source=tmp_path/'all.tsv';source.write_text('0\t1.25\n1\t3.5\n2\t0.1\n')
    dest=tmp_path/'annual.tsv';prefix_borders(source,dest,2)
    assert dest.read_text() == '0\t1.25\n1\t3.5\n'
    with pytest.raises(FileExistsError):
        prefix_borders(source,dest,2)
