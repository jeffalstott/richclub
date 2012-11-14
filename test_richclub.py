import unittest

from igraph import Graph
import richclub


class FirstTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Called once before all tests in this class."""
        ns = [6, ]
        ms = [10, ]
        directeds = [True, False]
        n_rewiress = [10, ]
        preserves = ['out', 'in']
        weight_ons = [True, False]

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

        for n, m, directed, n_rewires, preserve, weight_on in self.test_cases:
            print "%i nodes, %i links, directed: %i, %i rewires, "\
                "preserving %s, weights on: %i"\
                % (n, m, directed, n_rewires, preserve, weight_on)

            g = Graph.Erdos_Renyi(n=n, m=m, directed=directed)

            if weight_on:
                from numpy.random import rand
                g.es["weight"] = rand(m)

            gr = richclub.directed_spr(g, n_rewires=n_rewires,
                                       preserve=preserve)

            self.assertTrue(g.is_weighted() == gr.is_weighted())
            self.assertTrue(len(g.vs) == len(gr.vs))
            self.assertTrue(len(g.es) == len(gr.es))

            for mode in [1, 2, 3]:
                self.assertTrue(g.strength(mode=mode)
                                == gr.strength(mode=mode))

            if weight_on:
                from numpy import sort, all
                if preserve == 'out':
                    mode = 1
                if preserve == 'in':
                    mode = 2
                self.assertTrue(
                    all(sort(g.strength(mode=mode, weights=g.es["weight"]))
                        ==
                        sort(gr.strength(mode=mode, weights=gr.es["weight"]))
                        ))

    def test_rich_nodes(self):
        """Docstrings are printed during executions
        of the tests in the Eclipse IDE"""
        self.assertEqual(1, 1)

    def test_rich_club_coefficient(self):
        """Docstrings are printed during executions
        of the tests in the Eclipse IDE"""
        self.assertTrue(True)

if __name__ == '__main__':
    # execute all TestCases in the module
    unittest.main()
