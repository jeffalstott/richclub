import unittest

from igraph import Graph
import richclub
from numpy.testing import assert_array_equal, assert_allclose
from numpy import all, asarray
from numpy.random import rand


class FirstTestCase(unittest.TestCase):

    def test_directed_spr(self):
        print """Testing random control generation for directed weighted and
        unweighted graphs. Note that this code presently is not intended
        to function for undirected graphs, and in fact it does not."""

        ns = [60, ]
        ms = [60, 10, 600]
        directeds = [True, ]  # Not designed to work with undirected graphs
        n_rewiress = [2, 10]
        preserves = ['out', 'in']
        node_weight_averaging = [False, True]
        weight_ons = [False, True]

        test_cases = [(n, m, directed, n_rewires, preserve, nwa, weight_on)
                      for n in ns
                      for m in ms
                      for directed in directeds
                      for n_rewires in n_rewiress
                      for preserve in preserves
                      for nwa in node_weight_averaging
                      for weight_on in weight_ons]

        for n, m, directed, n_rewires, preserve, nwa, weight_on in test_cases:
            print "%i nodes, %i links, directed: %i, %i rewires, "\
                    "preserving %s, node weight averaging: %i, weights on: %i"\
                % (n, m, directed, n_rewires, preserve, nwa, weight_on)

            g = Graph.Erdos_Renyi(n=n, m=m, directed=directed)

            if weight_on:
                from numpy.random import rand
                g.es["weight"] = rand(m)

            gr = richclub.directed_spr(g, n_rewires=n_rewires,
                                       preserve=preserve,
                                       average_weight_by_node=nwa)


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
                if nwa:
                    for v in g.vs:
                        if preserve=='out':
                            weights = gr.es.select(_source_in=[v.index])["weight"]
                        if preserve=='in':
                            weights = gr.es.select(_target_in=[v.index])["weight"]
                        if weights:
                            assert_allclose(weights, weights[0])

    def test_richness_scores(self):
        print """Testing richness score identification on directed and undirected
        graphs."""

        ns = [60, ]
        ms = [60, 10, 600]
        directeds = [True, False]

        test_cases = [(n, m, directed)
                      for n in ns
                      for m in ms
                      for directed in directeds]

        for n, m, directed in test_cases:
            g = Graph.Erdos_Renyi(n=n, m=m, directed=directed)
            print "%i nodes, %i links, directed: %i" % (n, m, directed)

            g.es["weight"] = rand(m)
            ins = richclub.richness_scores(g, richness='in_strength')
            outs = richclub.richness_scores(g, richness='out_strength')
            s = richclub.richness_scores(g, richness='strength')

            if directed:
                assert_allclose(
                    asarray(ins)+asarray(outs),
                    asarray(s),
                    err_msg="Identified strength sequence not the sum of the "
                    "out and in strength sequences")
            else:
                assert_allclose(
                    asarray(ins),
                    asarray(s),
                    err_msg="In and total strength sequences not equal")

                assert_allclose(
                    asarray(outs),
                    asarray(s),
                    err_msg="Out and total strength sequences not equal")

    def test_rich_nodes(self):
        print """Testing rich/poor node selection on directed and undirected
        graphs."""

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
                len(richclub.rich_nodes(g, rank=0, highest=highest)),
                len(g.vs))

            if score == 'r':
                from numpy.random import rand
                from numpy import delete, greater, less
                score = rand(n)

                if highest:
                    op = greater
                else:
                    op = less

                for i in range(0, n - 1):
                    percentile = 100.0 * i / n
                    rnodes = richclub.rich_nodes(g, rank=percentile,
                                                 highest=highest, scores=score)
                    rscores = score[rnodes]
                    poorscores = delete(score, rnodes)

                    self.assertEqual(len(rnodes), n - i)
                    self.assertTrue(
                        all(
                            [op(a, b) for a in rscores for b in poorscores]
                        )
                    )

    def test_rich_club_coefficient(self):
        print """Testing rich_club_coefficient"""

        self.assertTrue(True)
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

        rc_properties = [
            'intensity_P_wm',
            'intensity_L_wm',
            'intensity_P_local',
            'intensity_L_local',
            'intensity_L_global',
            'intensity_P_global']

        for n, m, directed, highest, score in test_cases:
            g = Graph.Erdos_Renyi(n=n, m=m, directed=directed)
#            print "%i nodes, %i links, directed: %i, highest: %i, "\
#                "score: %s"\

            for p in rc_properties:
                print p
                rc = richclub.rich_club_coefficient(g, club_property=p)
            self.assertTrue(True)

if __name__ == '__main__':
    # execute all TestCases in the module
    unittest.main()
