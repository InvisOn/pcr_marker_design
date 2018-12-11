# Test uMelt functionality

# Really need a function to test if umelt service is
# a) working
# b) working as expected

# mocks to allow routine testing, assuming it is working

import numpy as np
import pytest
from scipy import interpolate
from pcr_marker_design import umelt_service as um


class TestUmeltService:
    @pytest.mark.umelt
    def test_overall_usage(self):
        test_seq = "TATAACCTGACTAACCATGAACCTGGGTAGAATTCCACTCCTCCACCAAATTTTTTAACTTAACCAAG"

        test_seq_hels = np.array([9.73479996e+01, 9.72409973e+01, 9.71320038e+01,
                                  9.70189972e+01, 9.69020004e+01, 9.67819977e+01,
                                  9.66569977e+01, 9.65270004e+01, 9.63919983e+01,
                                  9.62490005e+01, 9.60999985e+01, 9.59400024e+01,
                                  9.57699966e+01, 9.55859985e+01, 9.53860016e+01,
                                  9.51650009e+01, 9.49189987e+01, 9.46409988e+01,
                                  9.43229980e+01, 9.39560013e+01, 9.35270004e+01,
                                  9.30210037e+01, 9.24219971e+01, 9.17080002e+01,
                                  9.08580017e+01, 8.98499985e+01, 8.86640015e+01,
                                  8.72799988e+01, 8.56790009e+01, 8.38259964e+01,
                                  8.16230011e+01, 7.88099976e+01, 7.47559967e+01,
                                  6.82089996e+01, 5.75139999e+01, 4.24399986e+01,
                                  2.64400005e+01, 1.41599998e+01, 6.88999987e+00,
                                  3.20799994e+00, 1.47599995e+00, 6.83000028e-01,
                                  3.19999993e-01, 1.51999995e-01, 7.40000010e-02,
                                  3.70000005e-02, 1.89999994e-02, 9.99999978e-03,
                                  4.99999989e-03, 3.00000003e-03, 2.00000009e-03,
                                  1.00000005e-03, 1.00000005e-03, 0.00000000e+00,
                                  0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                  0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                  0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                  0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                  0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                  0.00000000e+00, 0.00000000e+00],
                                 dtype=np.float32)

        sequence = um.MeltSeq(test_seq)
        umelt = um.UmeltService()
        response = umelt.get_response(sequence)
        helicity = umelt.get_helicity_info(response)

        assert np.equal(helicity.helicity_data, test_seq_hels).all()


def approximately_equal(value_1, value_2, decimal_places):
    text_1 = '{0:.{1}f}'.format(value_1, decimal_places)
    text_2 = '{0:.{1}f}'.format(value_2, decimal_places)

    return text_1 == text_2


def make_spline(spline_x, spline_y):
    """Make the specified spline.

    spline_x and spline_y are array-like
    objects that define the spline on the
    x and y axes.

    Returns the point at which we find the
    steepest descent.
    """

    x = np.array(spline_x)
    y = np.array(spline_y)

    tck = interpolate.splrep(x, y, s=0)
    xnew = np.linspace(x.min(), x.max(), 50)
    dynew_dx = interpolate.splev(xnew, tck, der=1)

    return xnew[dynew_dx.argmin()]


class TestHelicityInfo:
    # TODO: use the lessons learned from here to actually test HelicityInfo
    def test_simple_slope(self):
        """Test that a simple slope is handled properly.
        """

        spline_x = [50, 55, 60, 65, 73, 77, 85, 90, 95, 100]
        spline_y = [100] * 5 + [10, 5, 3, 2, 0]

        melt_point = make_spline(spline_x, spline_y)

        assert approximately_equal(melt_point, 75.510204, 3)

    def test_double_slope_spline(self):
        """Test that the most extreme of two melting points is found.
        """

        spline_x = [50, 55, 60, 62, 65, 70, 72, 75, 80, 90, 95, 100]
        spline_y = [100] * 4 + [50] * 5 + [20, 10, 0]

        melt_point = make_spline(spline_x, spline_y)

        assert approximately_equal(melt_point, 63.265306, 3)
