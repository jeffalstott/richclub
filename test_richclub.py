import unittest

from igraph import Graph
import richclub


class FirstTestCase(unittest.TestCase):
    def test_directed_spr(self):
        """All methods beginning with 'test' are executed"""
        n = 6
        m = 10
        directed = True
        g = Graph.Erdos_Renyi(n=n, m=m, directed=directed)
        n_rewires = 10
        weighted = 'out'
        weights = True
        if weights:
            from numpy.random import rand
            g.es["weight"] = rand(m)

        gr = richclub.directed_spr(g, n_rewires=n_rewires, weighted=weighted)

        self.assertTrue(g.is_weighted() == gr.is_weighted())
        self.assertTrue(len(g.vs) == len(gr.vs))
        self.assertTrue(len(g.es) == len(gr.es))

        for mode in [1, 2, 3]:
            self.assertTrue(g.strength(mode=mode) == gr.strength(mode=mode))
            if weights:
                from numpy import sort, all
                self.assertTrue(
                    all(sort(g.strength(mode=mode, weights=g.es["weight"]))
                    == sort(gr.strength(mode=mode, weights=gr.es["weight"]))))

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
