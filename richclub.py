def directed_spr(G, n_rewires=10, preserve='out'):
    from numpy.random import randint

    g = G.copy()
    nes = len(g.es)

    i = 0
    rerolls = 0
    while i < (n_rewires * nes):
        if rerolls > n_rewires * nes:
            print "Halted after %i failed rewirings and"\
                "%i successful rewirings" % (rerolls, i)
            return g

        e1 = randint(nes)
        e2 = randint(nes)
        #In case we select the same edge twice, roll again.
        if e1 == e2:
            rerolls += 1
            continue

        s1 = g.es[e1].source
        t1 = g.es[e1].target
        a1 = g.es[e1].attributes()
        s2 = g.es[e2].source
        t2 = g.es[e2].target
        a2 = g.es[e2].attributes()

        #If s1 and s2 are the same, and we are preserving in strength, roll again.
        #This is because rewiring wouldn't accomplish anything.
        #The reverse is true if t1 and t2 are the same.
        if (s1 == s2 and preserve == 'in') or (t1 == t2 and preserve == 'out'):
            rerolls += 1
            continue

        #If either of the to-be-newly-wired connections already exist, roll again.
        #This prevents multiple edges going in the same direction between two nodes.
        if t2 in g.neighbors(s1, mode=1) or t1 in g.neighbors(s2, mode=1):
            rerolls += 1
            continue

        g.delete_edges([e1, e2])
        if preserve == 'out':  # Rewire the outgoing connections
            g.add_edge(s1, t2, **a1)
            g.add_edge(s2, t1, **a2)
        elif preserve == 'in':  # Rewire the incoming connections
            #Only difference is in the weight assignments
            g.add_edge(s1, t2, **a2)
            g.add_edge(s2, t1, **a1)
        #Note that in the unweighted condition these two methods are equivalent, so either option
        #for 'weighted' produces the same (correct) randomization, while preserving in degree
        # and out degree for each node

        i += 1
    return g


def rich_nodes(graph, rank=90.0, mode='percentile', highest=True, scores=None):
    """Extracts the "rich club" of the given graph, i.e. the subgraph spanned
    between vertices having the top X% of some score.

    Scores are given by the vertex degrees by default.

    @param graph:    the graph to work on
    @param rank: the rank of vertices to extract; must be between 0 and 1.
    @param highest:  whether to extract the subgraph spanned by the highest or
                     lowest scores.
    @param scores:   the scores themselves. C{None} uses the vertex degrees.
    """

    if scores is None:
        scores = graph.degree()

    if not highest:
        rank = 100-rank


    from numpy import where
    if mode=='percentile':
        from scipy.stats import scoreatpercentile
        threshold_score = scoreatpercentile(scores, rank)

    if highest:
        targets = where(scores>=threshold_score)[0]
    else:
        targets = where(scores<=threshold_score)[0]

    return targets

def rich_club_coefficient(graph, rank=None, highest=True, scores_name=None,
                          rewire=10, average=1, control=None, preserve=None):
    if type(rank) == float:
        rank = [rank]

    if rank is None:
        from numpy import arange
        rank = arange(10.0, 90.0, 10.0)

    from numpy import zeros

    if not graph.is_weighted():
        from numpy import ones
        #If given an unweighted graph, "strength"="degree", so all edges
        #have equal weight
        graph = graph.copy()
        graph.es["weight"] = ones(len(graph.es))

        #For unweighted graphs, preserving out strength vs in strength does
        #the same thing, so we'll just assign something to use.
        if not preserve:
            preserve = 'out'

    if scores_name is None or scores_name == 'strength':
        scores = graph.strength(graph.vs, mode=3, weights=graph.es["weight"])
        #Will not assign preserve here, because if the graph is actually weighted
        #the correct normalization is domain specific. ie. Up to the user!
    elif scores_name == 'out_strength':
        scores = graph.strength(graph.vs, mode=1, weights=graph.es["weight"])
        if not preserve:
            preserve = 'out'
    elif scores_name == 'in_strength':
        scores = graph.strength(graph.vs, mode=2, weights=graph.es["weight"])
        if not preserve:
            preserve = 'in'

    rc_coefficient = zeros(len(rank))

    for i in range(len(rank)):

        node_indices = rich_nodes(
            graph, rank=rank[i], highest=highest, scores=scores)

        if scores_name is None or scores_name == 'strength':
            numerator = sum(graph.es.select(_within=node_indices)["weight"])
            denominator = sum(
                graph.es.select(_between=(node_indices, graph.vs))["weight"])
        elif scores_name == 'out_strength':
            numerator = sum(graph.es.select(_within=node_indices)["weight"])
            denominator = sum(
                graph.es.select(_source_in=node_indices)["weight"])
        elif scores_name == 'in_strength':
            numerator = sum(graph.es.select(_within=node_indices)["weight"])
            denominator = sum(
                graph.es.select(_target_in=node_indices)["weight"])

        rc_coefficient[i] = numerator / denominator

    if control is not None:
        from igraph import Graph
        from numpy import ndarray
        if ~(type(control) == list or type(control) == ndarray):
            control = list(control)

        control_rc_coefficient = zeros(len(rank))
        for i in range(len(control)):

            random_graph = Graph.Weighted_Adjacency(
                control[i].toarray().tolist())

            control_rc_coefficient = control_rc_coefficient +\
                rich_club_coefficient(
                    random_graph, rank=rank, highest=highest,
                    scores_name=scores_name, rewire=False)

        control_rc_coefficient = control_rc_coefficient / len(control)

        return rc_coefficient / control_rc_coefficient
    elif rewire:
        control_rc_coefficient = zeros(len(rank))
        for i in range(average):

            random_graph = directed_spr(
                graph, n_rewires=rewire, preserve=preserve)

            control_rc_coefficient = control_rc_coefficient +\
                rich_club_coefficient(
                    random_graph, rank=rank, highest=highest,
                    scores_name=scores_name, rewire=False)

        control_rc_coefficient = control_rc_coefficient / average

        return rc_coefficient / control
    else:
        return rc_coefficient
