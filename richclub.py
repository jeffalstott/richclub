def directed_spr(G, n_rewires=10, preserve='out', average_weight_by_node=False):
    from numpy.random import randint

    g = G.copy()
    nes = len(g.es)

    if not g.is_directed():
        if 'weight' not in g.es.attributes():
            #Just use igraph's built in rewiring.
            g.rewire(n=n_rewires*nes)
            return g
        else:
            raise ValueError("directed_spr currently does not support rewiring "
            "of graphs that are both undirected and weighted. Have a good think "
            "on the best ways to preserve the degree and strength sequence of an"
            " undirected graph and get back to us with what you come up with.")
            ValueError

    if g.is_weighted() and average_weight_by_node:
        from numpy import mean
        for v in g.vs():
            if preserve=='out':
                edges = g.es.select(_source_in=[v.index])
            elif preserve=='in':
                edges = g.es.select(_target_in=[v.index])
            edges["weight"] = mean(edges["weight"])

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

def plot_rich_club(phis, ranks, ax=None, alpha=.1, **kwargs):

    from numpy import shape
    if len(shape(phis))>1:
        multiple_samples = True
        y = phis.mean(axis=0)
        from numpy import std
        error = std(phis, axis=0)
    else:
        multiple_samples = False
        y = phis
        error = 0

    if not ax:
        import matplotlib.pyplot as plt
        plt.plot(ranks, y, **kwargs)
        ax = plt.gca()
    else:
        ax.plot(ranks, y, **kwargs)

    if multiple_samples:
        ax.fill_between(ranks, y-error, y+error, alpha=alpha, **kwargs)

    xlims = (min(ranks), max(ranks))
    ax.set_xlim(xlims)
    return ax

def threshold_score(scores, rank=90.0, mode='percentile', highest=True,
        **kwargs):
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

def richness_scores(graph, richness=None):
    from types import FunctionType

    if richness is None or richness == 'strength':
        scores = graph.strength(graph.vs, mode=3, weights=graph.es["weight"])
    elif richness == 'out_strength':
        scores = graph.strength(graph.vs, mode=1, weights=graph.es["weight"])
    elif richness == 'in_strength':
        scores = graph.strength(graph.vs, mode=2, weights=graph.es["weight"])
    elif type(richness)==FunctionType:
        scores = richness(graph)
    else:
        try:
            len(richness)
            scores = richness
        except TypeError:
            raise ValueError("Unrecognized richness metric.")
    return scores

def rich_club_coefficient(graph, richness=None,
        club_property='intensity_P_wm',
        rank=None,
        controls=None,
        weightmax='max', candidate_edges_function=None,
        directed_local_drawn_from='out_links',
        **kwargs):

    from types import FunctionType, FloatType
    if type(rank) == FloatType:
        rank = [rank]

    if rank is None:
        from numpy import arange
        rank = arange(10.0, 100.0, 10.0)

    if not graph.is_weighted():
        from numpy import ones
        #If given an unweighted graph, "strength"="degree", so all edges
        #have equal weight
        graph = graph.copy()
        graph.es["weight"] = ones(len(graph.es))

    if weightmax=='max':
        weightmax = max(graph.es["weight"])

    scores = richness_scores(graph, richness=richness)

    from numpy import zeros
    rc_coefficient = zeros(len(rank))

    for i in range(len(rank)):

        rich_node_indices = rich_nodes(
            graph, rank=rank[i], scores=scores, **kwargs)

        rich_subgraph = graph.subgraph(rich_node_indices)

        if club_property==None or club_property=='C':
            rc_coefficient[i] = sum(rich_subgraph.es["weight"])
        elif club_property.startswith('intensity'):
            numerator = sum(rich_subgraph.es["weight"])

            if numerator==0:
                rc_coefficient[i] = 0
                continue

            if 'local' in club_property:
                target_nodes = rich_node_indices
            elif 'global' in club_property:
                target_nodes = graph.vs
            elif 'wm' in club_property:
                target_nodes = rich_node_indices
            else:
                raise ValueError("Unrecognized club_property metric.")

            if richness is None or richness == 'strength' or candidate_edges_function=='strength':
                candidate_edges = graph.es.select(_between=(target_nodes, graph.vs))
            elif richness == 'out_strength' or candidate_edges_function=='out_strength':
                candidate_edges = graph.es.select(_source_in=target_nodes)
            elif richness == 'in_strength' or candidate_edges_function=='in_strength':
                candidate_edges = graph.es.select(_target_in=target_nodes)
            elif candidate_edges_function:
                candidate_edges = candidate_edges_function(graph)
            else:
                raise ValueError("Unrecognized richness metric for identifying "
                "candidate edges to be used in calculating the denominator of "
                "intensity properties. Supply a candidate_edges_function.")

            if 'total' in club_property:
                number_to_count = len(candidate_edges)
            elif 'L' in club_property:
                number_to_count = len(rich_subgraph.es)
            elif 'P' in club_property:
                n = len(rich_subgraph.vs)
                number_to_count = n * (n - 1.0)
                if not graph.is_directed():
                    number_to_count = number_to_count / 2.0
            else:
                raise ValueError("Unrecognized club_property metric.")

            if 'wm' in club_property or 'weightmax' in club_property:
                denominator = number_to_count * weightmax
#            elif number_to_count > len(candidate_edges):
#                print("Fewer links present in the network than are sought"
#                        " for with these settings. Using all available.")
#                        "Try using the 'L'"
#                        " setting instead.")
#                from numpy import nan
#                denominator = nan
            elif 'local' in club_property:
            #The local option includes a requirement that each rich node cannot
            #contribute more links into the club than it can attach to all the
            #other rich nodes.
            #For directed networks, we must identify if the node is contributing
            #out links or in links to the club.
                n_nodes = len(graph.vs)
                n_rich = len(rich_node_indices)
                contribute_total = zeros(n_nodes)
                denominator = 0

                from numpy import argsort
                weights = candidate_edges["weight"]
                link_order = argsort(weights)
                j = 0
                limit = min(number_to_count, len(link_order))
                while j < limit:
                    largest_link = int(link_order[j])
                    e = candidate_edges[largest_link]
                    if directed_local_drawn_from=='out_links':
                        home = e.source
                    elif directed_local_drawn_from=='in_links':
                        home = e.target

                    j += 1
                    if contribute_total[home] >= (n_rich-1):
                        continue
                    else:
                        denominator += e["weight"]
                        contribute_total[home] += 1
                        number_to_count -= 1
            else:
                from numpy import sort
                candidate_edges = sort(candidate_edges["weight"])[::-1]
                denominator = sum(candidate_edges[:number_to_count])

            rc_coefficient[i] = numerator / denominator

        elif club_property == 'n_components':
            rc_coefficient[i] = len(rich_subgraph.components().subgraphs())
        elif club_property == 'diameter_weighted':
            rc_coefficient[i] = rich_subgraph.diameter(weights='weight')
        elif club_property == 'diameter':
            rc_coefficient[i] = rich_subgraph.diameter()
        elif club_property == 'clustering_weighted':
            rc_coefficient[i] = rich_subgraph.transitivity_avglocal_undirected(weights='weight')
        elif club_property == 'clustering':
            rc_coefficient[i] = rich_subgraph.transitivity_avglocal_undirected()
        elif club_property == 'n_infomap':
            infomap = rich_subgraph.community_infomap(edge_weights=rich_subgraph.es["weight"])
            rc_coefficient[i] = infomap.cluster_graph().vcount()
        elif club_property == 'q_infomap':
            infomap = rich_subgraph.community_infomap(edge_weights=rich_subgraph.es["weight"])
            rc_coefficient[i] = infomap.q
        elif club_property == 'codelength':
            infomap = rich_subgraph.community_infomap(edge_weights=rich_subgraph.es["weight"])
            rc_coefficient[i] = infomap.codelength
        elif type(club_property)==FunctionType:
            rc_coefficient[i] = club_property(graph, rich_node_indices)
        else:
            raise ValueError("Unrecognized club_property option.")

    if controls is None:
        return rc_coefficient
    else:
        from igraph import Graph
        from numpy import ndarray
        if ~(type(controls) == list or type(controls) == ndarray):
            controls = list(controls)

        control_rc_coefficient = zeros(len(rank))
        for i in range(len(controls)):

            control_graph = controls[i]

            if type(control_graph) != list:
                if type(control_graph)==ndarray:
                    control_graph = control_graph.tolist()
                else: 
                    from scipy.sparse import csc
                    if type(control_graph) == csc.csc_matrix:
                        control_graph - control_graph.toarray().tolist()
                    else:
                        raise TypeError("Can't parse the control graph type.")

            from numpy import transpose, all
            if all(control_graph==transpose(control_graph)):
                mode = 1
            else:
                mode = 0

            control_graph = Graph.Weighted_Adjacency(control_graph, mode=mode)

            control_rc_coefficient = control_rc_coefficient +\
                rich_club_coefficient(
                    control_graph,
                    club_property=club_property,
                    rank=rank,
                    controls=None,
                    weightmax=weightmax,
                    candidate_edges_function=candidate_edges_function,
                    directed_local_drawn_from=directed_local_drawn_from,
                    **kwargs)

        control_rc_coefficient = control_rc_coefficient / len(controls)

        return rc_coefficient / control_rc_coefficient

def normalized_rich_club_coefficient(graph, rewire=10, average=1, control=None,
                                     preserve=None, rank=None, **kwargs):

    if type(rank) == float:
        rank = [rank]

    if rank is None:
        from numpy import arange
        rank = arange(10.0, 100.0, 10.0)

    rc_coefficient = rich_club_coefficient(graph, rank=rank, **kwargs)

    from numpy import zeros

    if control is None and not preserve:
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
            raise ValueError("Must provide explicit control graphs or rewiring"
                             "preservation for this richness option.")

    if control is not None:
        from igraph import Graph
        from numpy import ndarray
        if ~(type(control) == list or type(control) == ndarray):
            control = list(control)

        control_rc_coefficient = zeros(len(rank))
        for i in range(len(control)):

            control_graph = control[i]

            if type(control_graph) != list:
                if type(control_graph)==ndarray:
                    control_graph = control_graph.tolist()
                else: 
                    from scipy.sparse import csc
                    if type(control_graph) == csc.csc_matrix:
                        control_graph - control_graph.toarray().tolist()
                    else:
                        raise TypeError("Can't parse the control graph type.")

            from numpy import transpose, all
            if all(control_graph==transpose(control_graph)):
                mode = 1
            else:
                mode = 0

            control_graph = Graph.Weighted_Adjacency(control_graph, mode=mode)

            control_rc_coefficient = control_rc_coefficient +\
                rich_club_coefficient(
                    control_graph, rank=rank, **kwargs)

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

        return rc_coefficient / control_rc_coefficient
    else:
        raise ValueError("Must provide explicit control graphs if"
                         "rewiring option is deactivated.")
