import unittest

from igraph import Graph
import richclub
from numpy.testing import assert_array_equal, assert_allclose
from numpy import all


class FirstTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Called once before all tests in this class."""
        ns = [60, ]
        ms = [60, 10, 600]
        directeds = [True, ]  # Not designed to work with undirected graphs
        n_rewiress = [2, 10]
        preserves = ['out', 'in']
        weight_ons = [False, True]

        cls.test_cases = [(n, m, directed, n_rewires, preserve, weight_on)
                          for n in ns
                          for m in ms
                          for directed in directeds
                          for n_rewires in n_rewiress
                          for preserve in preserves
                          for weight_on in weight_ons]
    pass

    def test_directed_spr(self):
        """All methods beginning with 'test' are executed"""

        ns = [60, ]
        ms = [60, 10, 600]
        directeds = [True, ]  # Not designed to work with undirected graphs
        n_rewiress = [2, 10]
        preserves = ['out', 'in']
        weight_ons = [False, True]

        test_cases = [(n, m, directed, n_rewires, preserve, weight_on)
                      for n in ns
                      for m in ms
                      for directed in directeds
                      for n_rewires in n_rewiress
                      for preserve in preserves
                      for weight_on in weight_ons]

        for n, m, directed, n_rewires, preserve, weight_on in test_cases:
            print "%i nodes, %i links, directed: %i, %i rewires, "\
                "preserving %s, weights on: %i"\
                % (n, m, directed, n_rewires, preserve, weight_on)

            g = Graph.Erdos_Renyi(n=n, m=m, directed=directed)

            if weight_on:
                from numpy.random import rand
                g.es["weight"] = rand(m)

            gr = richclub.directed_spr(g, n_rewires=n_rewires,
                                       preserve=preserve)

            self.assertEqual(g.is_weighted(), gr.is_weighted())
            self.assertEqual(len(g.vs), len(gr.vs))
            self.assertEqual(len(g.es), len(gr.es))

            for mode in [1, 2, 3]:
                assert_array_equal(
                    g.strength(mode=mode),
                    gr.strength(mode=mode),
                    err_msg="Degree sequence not equal in mode %i" % mode)

            if weight_on:
                from numpy import sort

                if preserve == 'out':
                    mode = 1
                if preserve == 'in':
                    mode = 2
                if not directed:
                    mode = 3

                assert_allclose(
                    sort(g.strength(mode=mode, weights=g.es["weight"])),
                    sort(gr.strength(mode=mode, weights=gr.es["weight"])),
                    err_msg="Strength sequence not equal")

    def test_rich_nodes(self):
        """Docstrings are printed during executions
        of the tests in the Eclipse IDE"""

        ns = [60, ]
        ms = [60, 10, 600]
        directeds = [True, False]
        highests = [True, False]
        scoress = [None, 'r']

        test_cases = [(n, m, directed, highest, score)
                      for n in ns
                      for m in ms
                      for directed in directeds
                      for highest in highests
                      for score in scoress]

        for n, m, directed, highest, score in test_cases:
            g = Graph.Erdos_Renyi(n=n, m=m, directed=directed)
            print "%i nodes, %i links, directed: %i, highest: %i, "\
                "score: %s"\
                % (n, m, directed, highest, score)

            #Returns all nodes
            self.assertEqual(
                len(richclub.rich_nodes(g, fraction=1, highest=highest)),
                len(g.vs))

            if score == 'r':
                from numpy.random import rand
                from numpy import delete, greater, less
                score = rand(n)

                if highest:
                    op = greater
                else:
                    op = less

                for i in range(n):
                    fraction = float(i + 1) / n
                    rnodes = richclub.rich_nodes(g, fraction=fraction,
                                                 highest=highest)
                    rscores = score[rnodes]
                    poorscores = delete(score, rnodes)

                    self.assertEqual(len(rnodes), i)

                    self.assertTrue(
                        all(
                            [op(a, b) for a in rscores for b in poorscores]
                        )
                    )

    def test_rich_club_coefficient(self):
        """Docstrings are printed during executions
        of the tests in the Eclipse IDE"""
        self.assertTrue(True)

if __name__ == '__main__':
    # execute all TestCases in the module
    unittest.main()
