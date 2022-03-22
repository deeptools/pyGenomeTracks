import unittest
import os
from pygenometracks import utilities
import matplotlib.pyplot as plt


ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")


class FakeAxis:
    """Allow Formatter to be called without having a "full" plot set up."""
    def __init__(self, vmin=1, vmax=10):
        self.vmin = vmin
        self.vmax = vmax

    def get_view_interval(self):
        return self.vmin, self.vmax


class TestUilitiesMethods(unittest.TestCase):

    def test_to_string_array(self):
        assert utilities.to_string([bytes("A", 'ascii'),
                                    bytes("B", 'ascii')])[1] == "B"

    def test_count_lines(self):
        with open(os.path.join(ROOT, "bedgraph_withNA.bdg"), 'r') as fh:
            assert utilities.count_lines(fh, asBed=True) == 3

    def test_transform_strange(self):
        input = [0, 1, 2]
        transformed = utilities.transform(input,
                                          "transformation",
                                          0, "file.txt")
        for i, t in zip(input, transformed):
            assert i == t

    def test_invalid_type_coverage_plot(self):
        fig, ax = plt.subplots()
        x_values = [0, 1]
        score_list = [0, 1]
        plot_type = "nothing"
        size = .5
        color = "blue"
        negative_color = color
        alpha = 0.5
        grid = False
        n_children = len(ax.get_children())
        utilities.plot_coverage(ax, x_values, score_list, plot_type, size, color, negative_color, alpha, grid)
        assert len(ax.get_children()) == n_children + 1


class TestFormatter(unittest.TestCase):

    def test_easy_cases_b(self):
        fmt = utilities.MyBasePairFormatter()
        fmt.axis = FakeAxis(0, 20)
        fmt.set_locs([0, 10, 20])
        assert fmt(0) == '0'
        assert fmt(10) == '10 b'
        assert fmt(20) == '20'

    def test_easy_cases_Kb(self):
        fmt = utilities.MyBasePairFormatter()
        fmt.axis = FakeAxis(74500000, 74520000)
        fmt.set_locs([74500000, 74510000, 74520000])
        assert fmt(74500000) == '74,500'
        assert fmt(74510000) == '74,510 Kb'
        assert fmt(74520000) == '74,520'

    def test_easy_cases_Mb(self):
        fmt = utilities.MyBasePairFormatter()
        fmt.axis = FakeAxis(73500000, 74520000)
        fmt.set_locs([73500000, 74000000, 74500000])
        assert fmt(73500000) == '73.5'
        assert fmt(74000000) == '74.0 Mb'
        assert fmt(74500000) == '74.5'

    def test_ill_case(self):
        fmt = utilities.MyBasePairFormatter()
        fmt.axis = FakeAxis(0, 15000)
        fmt.set_locs([0, 10, 20, 15000])
        assert fmt(0) == '0.00'
        assert fmt(10) == '0.01'
        assert fmt(20) == '0.02 Kb'
        assert fmt(15000) == '15.00'

    def test_ill_case_2(self):
        fmt = utilities.MyBasePairFormatter()
        fmt.axis = FakeAxis(0, 0)
        fmt.set_locs([0])
        assert fmt(0) == '0 b'

    def test_ill_case_3(self):
        fmt = utilities.MyBasePairFormatter()
        fmt.axis = FakeAxis(0, 0)
        fmt.set_locs([0, 0])
        assert fmt(0) == '0 b'

    def test_no_axe(self):
        fmt = utilities.MyBasePairFormatter()
        assert fmt(0) == ""
