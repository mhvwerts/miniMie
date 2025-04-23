import unittest
import mie
import numpy as np

class MieTest(unittest.TestCase):
    def test_mie_matzler(self):
        '''Takes the values for the mie function from the 2002 documentation and
        checks against the provided output.'''
        m = 5+0.4j
        x = 1
        result = mie.mie(m, x)
        expected = np.array([5, 0.4, 1, 1.9794, 0.8795, 1.0999, 1.1133, -0.0595,
                            1.2664])

        self.assertTrue(np.isclose(result, expected, atol=5E-4).all())

    def test_mie_large_x(self):
        '''Tests that an x value of 70000 returns non-NaN's'''
        m = 2+0.01j
        x = 70000
        result = mie.mie(m, x)

        self.assertFalse(np.isnan(result.any()))

    def test_mie_ab(self):
        '''Tests mie_ab against MATLAB determined values.'''
        m = 2 + 1j
        x = 5
        result = mie.mie_ab(m, x)
        expected = np.array([[0.6887-0.1051j, 0.4510+0.1921j, 0.3598-0.0839j,
            0.5119-0.0942j, 0.4050-0.0503j, 0.1423-0.0726j, 0.0182-0.0202j,
            0.0016-0.0025j, 0.0001-0.0002j, 0.0000-0.0000j, 0.0000-0.0000j,
            0.0000-0.0000j, 0.0000-0.0000j, 0.0000-0.0000j],
            [0.2994+0.1148j, 0.5542-0.2417j, 0.7419+0.1349j, 0.4040+0.3102j,
            0.1448+0.1685j, 0.0423+0.0430j, 0.0082+0.0054j, 0.0010+0.0004j,
            0.0001+0.0000j, 0.0000+0.0000j, 0.0000-0.0000j, 0.0000-0.0000j,
            0.0000-0.0000j, 0.0000-0.0000j]])

        self.assertTrue(np.isclose(result, expected, atol=5E-4).all())

    def test_mie_abcd(self):
        '''Tests mie_abcd against MATLAB etermined values.'''
        m = 2 + 1j
        x = 5
        result = mie.mie_abcd(m, x)
        expected = np.array([[0.6887-0.1051j, 0.4510+0.1921j, 0.3598-0.0839j,
            0.5119-0.0942j, 0.4050-0.0503j, 0.1423-0.0726j, 0.0182-0.0202j,
            0.0016-0.0025j, 0.0001-0.0002j, 0.0000-0.0000j, 0.0000-0.0000j,
            0.0000-0.0000j, 0.0000-0.0000j, 0.0000-0.0000j],
            [0.2994+0.1148j, 0.5542-0.2417j, 0.7419+0.1349j,  0.4040+0.3102j,
            0.1448+0.1685j, 0.0423+0.0430j, 0.0082+0.0054j, 0.0010+0.0004j,
            0.0001+0.0000j,  0.0000+0.0000j, 0.0000-0.0000j, 0.0000-0.0000j,
            0.0000-0.0000j, 0.0000-0.0000j],
            [0.0030-0.0094j, 0.0005-0.0105j, -0.0040-0.0107j, -0.0102-0.0072j,
            -0.0125+0.0025j, -0.0038+0.0101j, 0.0046+0.0053j, 0.0037-0.0005j,
            0.0009-0.0015j, -0.0002-0.0008j, -0.0003-0.0002j, -0.0001+0.0000j,
            -0.0000+0.0000j, -0.0000+0.0000j],
            [0.0031-0.0095j, 0.0007-0.0110j, -0.0041-0.0121j, -0.0124-0.0094j,
            -0.0183+0.0051j, 0.0004+0.0172j, 0.0088+0.0023j, 0.0029-0.0028j,
            -0.0001-0.0017j, -0.0005-0.0005j, -0.0003-0.0000j, -0.0001+0.0001j,
            -0.0000+0.0001j, 0.0000+0.0000j]])

        self.assertTrue(np.isclose(result, expected, atol=5E-4).all())

    def test_mie_pt_invalid_mu_low(self):
        '''Tests mie_pt with an invalid mu value on the low side. Expect a
        ValueError.'''

        μ = -2
        nmax = 20

        self.assertRaises(ValueError, mie.mie_pt, μ, nmax)

    def test_mie_pt_invalid_mu_hi(self):
        '''Tests mie_pt with an invalid mu value on the high side. Expect a
        ValueError.'''

        μ = 2
        nmax = 20

        self.assertRaises(ValueError, mie.mie_pt, μ, nmax)

    def test_mie_pt_invalid_nmax(self):
        '''Tests mie_pt with an invalid nmax value. Expect a ValueError'''

        μ = 0.5
        nmax = 1

        self.assertRaises(ValueError, mie.mie_pt, μ, nmax)

    def test_mie_pt(self):
        '''Tests mie_pt against MATLAB output.'''

        μ = 0.5
        nmax = 20

        result = mie.mie_pt(μ, nmax)
        expected = np.array([[1.0000, 1.5000, 0.3750, -1.5625, -2.2266, -0.5742,
            1.9756, 2.7729, 0.7237, -2.3171, -3.2291, -0.8481, 2.6147, 3.6286,
            0.9567, -2.8819, -3.9885, -1.0544, 3.1265, 4.3186],
            [0.5000, -1.5000, -5.4375, -5.0000, 3.8086, 13.8633, 11.5083,
            -6.6885, -24.4727, -19.5466, 10.0456, 36.8895, 28.8690, -13.8208,
            -50.8828, -39.3196, 17.9726, 66.2925, 50.7883, -22.4698]])

        self.assertTrue(np.isclose(result, expected, atol=5E-4).all())

    def test_miecoated_ab1(self):
        '''Tests miecoated_ab1 against MATLAB output.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 2
        y = 2.5

        result = mie.miecoated_ab1(m1, m2, x, y)
        expected = np.array([[0.3449-0.2907j, 0.6739-0.1695j, 0.3100-0.1967j,
            0.0170-0.0356j, 0.0007-0.0022j, 0.0000-0.0001j, 0.0000-0.0000j,
            0.0000-0.0000j, 0.0000-0.0000j, 0.0000-0.0000j],
            [0.6967+0.3234j, 0.1775+0.2776j, 0.0283+0.0626j, 0.0033+0.0061j,
            0.0002+0.0003j, 0.0000+0.0000j, 0.0000+0.0000j, 0.0000+0.0000j,
            0.0000+0.0000j, 0.0000+0.0000j]])

        self.assertTrue(np.isclose(result, expected, atol=5E-4).all())

    def test_miecoated_ab1_x_eq_y(self):
        '''Test miecoated_ab1 given the invalid input where x == y.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 2
        y = x

        self.assertRaises(ValueError, mie.miecoated_ab1, m1, m2, x, y)

    def test_miecoated_ab1_x_gt_y(self):
        '''Test miecoated_ab1 given the invalid input where x > y.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 3
        y = 2

        self.assertRaises(ValueError, mie.miecoated_ab1, m1, m2, x, y)

    def test_miecoated_ab2(self):
        '''Tests miecoated_ab2 against MATLAB output.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 2
        y = 2.5

        result = mie.miecoated_ab2(m1, m2, x, y)
        expected = np.array([[0.3449-0.2907j, 0.6739-0.1695j, 0.3100-0.1967j,
            0.0170-0.0356j, 0.0007-0.0022j, 0.0000-0.0001j, 0.0000-0.0000j,
            0.0000-0.0000j, 0.0000-0.0000j, 0.0000-0.0000j],
            [0.6967+0.3234j, 0.1775+0.2776j, 0.0283+0.0626j, 0.0033+0.0061j,
            0.0002+0.0003j, 0.0000+0.0000j, 0.0000+0.0000j, 0.0000+0.0000j,
            0.0000+0.0000j, 0.0000+0.0000j]])

        self.assertTrue(np.isclose(result, expected, atol=5E-4).all())

    def test_miecoated_ab2_x_eq_y(self):
        '''Test miecoated_ab2 given the invalid input where x == y.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 2
        y = x

        self.assertRaises(ValueError, mie.miecoated_ab2, m1, m2, x, y)

    def test_miecoated_ab2_x_gt_y(self):
        '''Test miecoated_ab2 given the invalid input where x > y.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 3
        y = 2

        self.assertRaises(ValueError, mie.miecoated_ab2, m1, m2, x, y)

    def test_miecoated_ab3(self):
        '''Tests miecoated_ab3 against MATLAB output.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 2
        y = 2.5

        result = mie.miecoated_ab3(m1, m2, x, y)
        expected = np.array([[0.344887544998-0.290674498774j,
            0.673933319361-0.169511279685j, 0.310004900520-0.196681430837j,
            0.016972701487-0.035606508715j, 0.000725457871-0.002213356407j,
            0.000027433206-0.000096523054j, 0.000000819366-0.000003117580j,
            0.000000019297-0.000000077175j, 0.000000000364-0.000000001508j,
            0.000000000005-0.000000000023j],
            [0.696725295929+0.323444200126j, 0.177506851256+0.277587920106j,
            0.028310331541+0.062588439288j, 0.003320585625+0.006085965844j,
            0.000222003932+0.000328944968j, 0.000009340496+0.000011871977j,
            0.000000275717+0.000000313029j, 0.000000006098+0.000000006351j,
            0.000000000105+0.000000000102j, 0.000000000001+0.000000000001j]])

        self.assertTrue(np.isclose(result, expected, atol=1E-11).all())

    def test_miecoated_ab3_x_eq_y(self):
        '''Test miecoated_ab3 given the invalid input where x == y.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 2
        y = x

        self.assertRaises(ValueError, mie.miecoated_ab3, m1, m2, x, y)

    def test_miecoated_ab3_x_gt_y(self):
        '''Test miecoated_ab3 given the invalid input where x > y.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 3
        y = 2

        self.assertRaises(ValueError, mie.miecoated_ab3, m1, m2, x, y)

    def test_miecoated_opt1(self):
        '''Test miecoated against MATLAB outputs using opt 1.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 2
        y = 2.5
        opt = 1

        result = mie.miecoated(m1, m2, x, y, opt)
        expected = np.array([3.182016828599, 2.025247626053, 1.156769202546,
            0.637637341412, 0.560757401303, 0.314844137185])

        self.assertTrue(np.isclose(result, expected, atol=1E-11).all())

    def test_miecoated(self):
        '''Test miecoated given the invalid input where x > y.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 3
        y = 2
        opt = 1

        self.assertRaises(ValueError, mie.miecoated, m1, m2, x, y, opt)

    def test_miecoated_opt2(self):
        '''Test miecoated against MATLAB outputs using opt 2.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 2
        y = 2.5
        opt = 2

        result = mie.miecoated(m1, m2, x, y, opt)
        expected = np.array([3.182016828600, 2.025247626053, 1.156769202546,
            0.637637341412, 0.560757401303, 0.314844137185])

        self.assertTrue(np.isclose(result, expected, atol=1E-11).all())

    def test_miecoated_opt3(self):
        '''Test miecoated against MATLAB outputs using opt 3.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 2
        y = 2.5
        opt = 3

        result = mie.miecoated(m1, m2, x, y, opt)
        expected = np.array([3.182016828599, 2.025247626053, 1.156769202546,
            0.637637341412, 0.560757401303, 0.314844137185])

        self.assertTrue(np.isclose(result, expected, atol=1E-11).all())

    def test_miecoated_S12(self):
        '''Test miecoated_S12 against MATLAB outputs.'''

        m1 = 1+2j
        m2 = 1.5+2.5j
        x = 2
        y = 2.5
        u = 0.75

        result = mie.miecoated_S12(m1, m2, x, y, u)
        expected = np.array([3.064118537580-0.772416667613j,
            1.435813972397+1.203900827712j])

        self.assertTrue(np.isclose(result, expected, atol=1E-11).all())

if __name__ == '__main__':
    unittest.main()
