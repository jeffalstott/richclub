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


def threshold_score(scores, rank=90.0, mode='percentile', highest=True):
    if not highest:
        rank = 100 - rank

    if mode == 'percentile':
        from scipy.stats import scoreatpercentile
        threshold_score = scoreatpercentile(scores, rank)

    return threshold_score


def rich_nodes(graph, scores=None, highest=True, **kwargs):
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

    ts = threshold_score(scores, highest=highest, **kwargs)

    from numpy import where
    if highest:
        targets = where(scores >= ts)[0]
    else:
        targets = where(scores <= ts)[0]

    return targets


def rich_club_coefficient(graph, richness=None, club_property=None,
                          rank=None, weightmax=1, **kwargs):
    if type(rank) == float:
        rank = [rank]

    if rank is None:
        from numpy import arange
        rank = arange(10.0, 90.0, 10.0)

    if not graph.is_weighted():
        from numpy import ones
        #If given an unweighted graph, "strength"="degree", so all edges
        #have equal weight
        graph = graph.copy()
        graph.es["weight"] = ones(len(graph.es))

    if richness is None or richness == 'strength':
        scores = graph.strength(graph.vs, mode=3, weights=graph.es["weight"])
    elif richness == 'out_strength':
        scores = graph.strength(graph.vs, mode=1, weights=graph.es["weight"])
    elif richness == 'in_strength':
        scores = graph.strength(graph.vs, mode=2, weights=graph.es["weight"])
    else:
        raise ValueError("Unrecognized richness metric.")

    from numpy import zeros
    rc_coefficient = zeros(len(rank))

    for i in range(len(rank)):

        rich_node_indices = rich_nodes(
            graph, rank=rank[i], scores=scores, **kwargs)

        rich_subgraph = graph.subgraph(rich_node_indices)

        if club_property.startswith('intensity'):
            numerator = sum(rich_subgraph.es["weight"])

            if 'local' in club_property:
                target_nodes = rich_node_indices
            elif 'global' in club_property:
                target_nodes = graph.vs
            else:
                raise ValueError("Unrecognized club_property metric.")

            if richness is None or richness == 'strength':
                candidate_edges = graph.es.select(_between=(target_nodes, graph.vs))["weight"]
            elif richness == 'out_strength':
                candidate_edges = graph.es.select(_source_in=target_nodes)["weight"]
            elif richness == 'in_strength':
                candidate_edges = graph.es.select(_target_in=target_nodes)["weight"]
            else:
                raise ValueError("Unrecognized richness metric.")

            if 'total' in club_property:
                denominator = sum(candidate_edges)
            elif 'topN_' in club_property:
                from numpy import sort
                candidate_edges = sort(candidate_edges)[::-1]
                N = len(rich_subgraph.es)
                denominator = sum(candidate_edges[:N])
            elif 'topNp_' in club_property:
                from numpy import sort
                candidate_edges = sort(candidate_edges)[::-1]
                n = len(rich_subgraph.vs)
                Np = n * (n - 1)
                denominator = sum(candidate_edges[:Np])
            elif 'topNpweightmax' in club_property:
                n = len(rich_subgraph.vs)
                Np = n * (n - 1)
                denominator = Np * weightmax
            else:
                raise ValueError("Unrecognized club_property metric.")

            rc_coefficient[i] = numerator / denominator

        elif club_property == 'clustering':
            rc_coefficient[i] = rich_subgraph.transitivity_avglocal_undirected(weights='weight')
        elif club_property == 'n_infomap':
            infomap = rich_subgraph.community_infomap(edge_weights=rich_subgraph.es["weight"])
            rc_coefficient[i] = infomap.cluster_graph().vcount()
        elif club_property == 'q_infomap':
            infomap = rich_subgraph.community_infomap(edge_weights=rich_subgraph.es["weight"])
            rc_coefficient[i] = infomap.q
        elif club_property == 'codelength':
            infomap = rich_subgraph.community_infomap(edge_weights=rich_subgraph.es["weight"])
            rc_coefficient[i] = infomap.codelength
        else:
            raise ValueError("Unrecognized club_property option.")


def normalized_rich_club_coefficient(graph, rewire=10, average=1, control=None,
                                     preserve=None, rank=None, **kwargs):

    if type(rank) == float:
        rank = [rank]

    if rank is None:
        from numpy import arange
        rank = arange(10.0, 90.0, 10.0)
    rc_coefficient = rich_club_coefficient(graph, **kwargs)

    from numpy import zeros

    if not preserve:
        if not graph.is_weighted():
            #For unweighted graphs, preserving out strength vs in strength does
            #the same thing, so we'll just assign something to use.
            preserve = 'out'
        elif kwargs['richness'] == 'out_strength':
            preserve = 'out'
        elif kwargs['richness'] == 'in_strength':
            preserve = 'in'
        elif kwargs['richness'] is None or kwargs['richness'] == 'strength':
            #Will not assign preserve here, because if the graph is actually weighted
            #the correct normalization is domain specific. ie. Up to the user!
            raise ValueError("Must provide explicit control graphs"
                             "for this richness option.")

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
                    random_graph, rank=rank, **kwargs)

        control_rc_coefficient = control_rc_coefficient / len(control)

        return rc_coefficient / control_rc_coefficient
    elif rewire:
        control_rc_coefficient = zeros(len(rank))
        for i in range(average):

            random_graph = directed_spr(
                graph, n_rewires=rewire, preserve=preserve)

            control_rc_coefficient = control_rc_coefficient +\
                rich_club_coefficient(
                    random_graph, rank=rank, **kwargs)

        control_rc_coefficient = control_rc_coefficient / average

        return rc_coefficient / control
    else:
        raise ValueError("Must provide explicit control graphs if"
                         "rewiring option is deactivated.")
