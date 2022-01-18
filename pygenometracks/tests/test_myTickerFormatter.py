# -*- coding: utf-8 -*-
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pygenometracks.utilities import MyBasePairFormatter
from tempfile import NamedTemporaryFile
import pytest
mpl.use('agg')

outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                             delete=False)


class TestMyBasePairFormat:
    params = [
        (0, 1e8, None),  # no major tick => no lab
        (1, 1e8, ['0 b']),   # a single major tick => 0 b
        (2, 1e8, ['0 Mb', '100']),   # 2 major ticks => 0 Mb and 100
        (3, 2e5, ['0', '100 Kb', '200'])   # 2 major ticks => 0 Mb and 100
    ]

    @pytest.mark.parametrize('nb_majorticks, xmax, expected_text', params)
    def test_starting_at_0(
            self, nb_majorticks, xmax, expected_text):
        fig, ax = plt.subplots()
        ax.set_xlim(0, xmax)
        fig.savefig(outfile.name)
        ax.xaxis.set_major_formatter(MyBasePairFormatter())
        ax.set_xticks(np.linspace(0, xmax, nb_majorticks))
        fig.savefig(outfile.name)
        if nb_majorticks == 0:
            assert len(ax.xaxis.get_majorticklabels()) == 0
        else:
            for text, expected in zip(list(ax.xaxis.get_majorticklabels()), expected_text):
                assert text.get_text() == expected
