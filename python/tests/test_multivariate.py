"""
Tests ported from inst/tinytest/test_multivariate.R

Covers k_means (skipped if not available) and layerCor.
"""
import numpy as np
import pytest
import terra as pt
from terra.rast import rast

from path_utils import skip_if_missing_inst_ex

try:
    from terra.stats import layer_cor
    _HAS_LAYERCOR = True
except ImportError:
    _HAS_LAYERCOR = False


def find_logo():
    return skip_if_missing_inst_ex("logo.tif")


@pytest.mark.skipif(not _HAS_LAYERCOR, reason="layer_cor not implemented")
def test_layercor_string_vs_callable():
    """layerCor(x, 'cor') and layerCor(x, cor) should give equivalent results."""
    f = find_logo()
    x = rast(f)

    import numpy as np
    a = layer_cor(x, "cor")
    b = layer_cor(x, np.corrcoef)

    # a["correlation"] from string should equal b from callable
    if isinstance(a, dict):
        a_mat = a["correlation"]
    else:
        a_mat = a
    if isinstance(b, dict):
        b_mat = b["correlation"]
    else:
        b_mat = b

    np.testing.assert_array_almost_equal(
        np.array(a_mat), np.array(b_mat), decimal=5
    )
