from numpy import where

def fast_RC(A, sequence=None):
    from numpy import sum, argsort, arange
    strengths = sum(A, axis=0)
    if sequence is None:
        sequence = argsort(strengths)
    A = A[sequence][:,sequence]
    x = [sum(A[n:,n:]) for n in arange(len(A))]
    return x

class RC(object):
    
    def __init__(self,g,scores=None,mode=2, ranks=None):
	from numpy import argsort
        self.g = g
        if scores is None:
	    from numpy import array
            self.scores = array(g.strength(weights='weight',mode=mode))
        else:
            self.scores = scores
        self.order = argsort(self.scores)
        self.n_nodes = g.vcount()
        
        if ranks is not None:
            self.ranks = ranks
	    from numpy import where
            self.order = [where((self.ranks[r]<=self.scores) * 
                                (self.scores<self.ranks[r+1]))[0] for r in range(len(self.ranks)-1)]

        #    from bisect import bisect_left
        #    self.tiers = [bisect_left(self.ranks, score) for score in self.scores]
        #    order = argsort(self.tiers)
        #    tier_breaks = concatenate(([0],where(diff(self.tiers))[0]+1))
        #    self.order = [order[tier_breaks[i]:tier_breaks[i+1]] for i in arange(len(tier_breaks)-1)]
        
    def _ts(self, rank):
	from numpy import where
        s = where(self.scores<rank)[0]
        self.scores = self.scores[self.scores>rank]
        return s
    
    def _rc(self, n):
	try:
            self.G.delete_vertices(map(str,n))
	except TypeError:
            self.G.delete_vertices(str(n))
        return self.G
    
    def phis(self,**kwargs):
        f = lambda g: sum(g.es['weight'])
        return self.fun(f,**kwargs)
    
    def fun(self,f, ns=None):
        if ns is None:
            ns = self.n_nodes
        self.G = self.g.copy()
        self.G.vs['name'] = map(str,range(self.G.vcount()))
        return [f(self._rc(n)) for n in self.order[:ns]]

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
    if not highest and mode == 'percentile':
        rank = 100 - rank

    if mode == 'percentile':
        from scipy.stats import scoreatpercentile
        threshold_score = scoreatpercentile(scores, rank)
    elif mode == 'raw':
        threshold_score = rank

    return threshold_score


def rich_nodes(graph, scores=None, highest=True, mode='percentile', **kwargs):

    if scores is None:
        scores = graph.degree()

    ts = threshold_score(scores, highest=highest, mode=mode, **kwargs)

    if highest:
        return where(scores >= ts)[0]
    else:
        return where(scores <= ts)[0]

def richness_scores(graph, richness=None):
    from types import FunctionType

    if richness is None or richness == 'strength':
        scores = graph.strength(graph.vs, mode=3, weights=graph.es["weight"])
    elif richness == 'out_strength':
        scores = graph.strength(graph.vs, mode=1, weights=graph.es["weight"])
    elif richness == 'in_strength':
        scores = graph.strength(graph.vs, mode=2, weights=graph.es["weight"])
    elif richness=='degree':
        scores = graph.degree()
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
        club_property='C',
        rank=None,
        controls=None,
        weightmax='max', candidate_edges_function=None,
        directed_local_drawn_from='out_links',
        mode='percentile',
        **kwargs):

    if not graph.is_weighted():
        from numpy import ones
        #If given an unweighted graph, "strength"="degree", so all edges
        #have equal weight
        graph = graph.copy()
        graph.es["weight"] = ones(len(graph.es))

    if weightmax=='max':
        weightmax = max(graph.es["weight"])

    scores = richness_scores(graph, richness=richness)

    from types import FunctionType, FloatType
    if type(rank) == FloatType:
        rank = [rank]

    if rank is None:
        if mode=='percentile':
            from numpy import arange
            rank = arange(10.0, 100.0, 10.0)
        elif mode=='raw':
            from numpy import unique
            rank = unique(scores)

    from numpy import zeros
    rc_coefficient = zeros(len(rank))

    #return  [sum(graph.subgraph(where(
    #    graph.strength(weights='weight',mode=2)>rank[i])[0]
    #    ).es['weight']) for i in range(len(rank))]
    for i in range(len(rank)):

        rich_node_indices = rich_nodes(
            graph, rank=rank[i], scores=scores, mode=mode, **kwargs)

        rich_subgraph = graph.subgraph(rich_node_indices)

        if club_property==None or club_property=='C':
            rc_coefficient[i] = sum(rich_subgraph.es["weight"])
        elif type(club_property)==FunctionType:
            rc_coefficient[i] = club_property(graph, rich_node_indices)
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

            if type(control_graph) != Graph:
#If the graph isn't an igraph Graph, we'll turn it into a list, if it's
#not already one.
                if type(control_graph) != list:
                    if type(control_graph)==ndarray:
                        control_graph = control_graph.tolist()
                    else: 
                        from scipy.sparse import csc
                        if type(control_graph) == csc.csc_matrix:
                            control_graph - control_graph.toarray().tolist()
                        else:
                            raise TypeError("Can't parse the control"
                                    "graph type.")

                from numpy import transpose, all
                if all(control_graph==transpose(control_graph)):
                    directed_mode = 1
                else:
                    directed_mode = 0

                control_graph = Graph.Weighted_Adjacency(control_graph,
                        mode=directed_mode)

            control_rc_coefficient = control_rc_coefficient +\
                rich_club_coefficient(
                    control_graph,
                    club_property=club_property,
                    rank=rank,
                    controls=None,
                    mode=mode,
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

def preserve_strength(g, randomize_topology=False, preserve_mode='estimate_both', permute_strength=True, randomize_method='vl'):
    from numpy import array, ndarray
    from scipy.sparse import csc
    if (type(g)==ndarray or type(g)==csc.csc_matrix) and preserve_mode in ['estimate_both', 3] and not randomize_topology:
        A = g
        from numpy import sum
        out_strength = sum(A, axis=1)
        in_strength = sum(A, axis=0)
        adj = 1*(A>0)
        out_degree = sum(adj, axis=1)
        in_degree = sum(adj, axis=0)
              
        if permute_strength:
            from numpy.random import permutation
            from numpy import unique
            out_in_sequence = array(zip(out_degree, in_degree))
            for k in unique(out_in_sequence):
                k_ind = where((out_in_sequence == k).all(axis=1))[0]
                new_ind = permutation(k_ind)
                out_strength[k_ind] = out_strength[new_ind]
                in_strength[k_ind] = in_strength[new_ind]

        from numpy import mean, outer, logical_not
        mean_k = mean([out_degree,in_degree])
        mean_s = mean([out_strength,in_strength])
        G = (mean_k/mean_s) * outer(out_strength,in_strength) /outer(out_degree,in_degree)
        G[logical_not(adj)] = 0
        return G

    out_degree = g.degree(mode=1)
    in_degree = g.degree(mode=2)
   
    if randomize_topology:
        if g.is_directed():
            #Check if all edges are bidirectional.
            #If so, create a random graph with only bidirectional edges.
            G = g.copy()
            G.to_undirected(mode=False)
            if all(array(G.count_multiple())==2):
                G = g.Degree_Sequence(out_degree, method=randomize_method)
                G.to_directed()
            else:
                G = g.copy()
                G.rewire()
        else:
            G = g.Degree_Sequence(out_degree, method=randomize_method)
    else:
        G = g.copy()

    if preserve_mode in ['estimate_both', 3]:
        
        out_strength = array(g.strength(mode=1, weights='weight'))
        in_strength = array(g.strength(mode=2, weights='weight'))

        if permute_strength:
            from numpy.random import permutation
            from numpy import unique
            out_in_sequence = array(zip(out_degree, in_degree))
            for k in unique(out_in_sequence):
                k_ind = where((out_in_sequence == k).all(axis=1))[0]
                new_ind = permutation(k_ind)
                out_strength[k_ind] = out_strength[new_ind]
                in_strength[k_ind] = in_strength[new_ind]
    
        from numpy import mean
        mean_k = mean([out_degree,in_degree])
        mean_s = mean([out_strength,in_strength])
        
        for e in G.es:
            e["weight"] = ( (mean_k/mean_s) * out_strength[e.source] * 
                in_strength[e.target] /
                (out_degree[e.source]*in_degree[e.target]) )
        return G
    
    elif preserve_mode in ['out', 1]:
        preserve_mode = 1
    elif preserve_mode in ['in', 2]:
        preserve_mode = 2
    
    from numpy import sum
    ind = [g.incident(v,mode=preserve_mode) for v in range(g.vcount())]
    weights = g.es[sum(ind)]['weight']
    from numpy.random import shuffle
    map(shuffle,ind)
    G.es[sum(ind)]['weight'] = weights

    return G

def just_plot_rc(rc, rc_controls, ax=None, x=None, top_limit=None, **kwargs):
    from numpy import array, mean
    from numpy.ma import masked_invalid
    
    if x is None:
        from numpy import arange
        x = arange(len(rc))
    
    if ax is None:
        ax = gca()
            
    from scipy.stats import scoreatpercentile
    upper95 = array([scoreatpercentile(rc_controls[:,i], 97.5) for i in range(len(x))])
    lower95 = array([scoreatpercentile(rc_controls[:,i], 2.5) for i in range(len(x))])
    
    control_mean = mean(rc_controls,axis=0)
    rc_norm = rc/control_mean
    rc_norm = masked_invalid(rc_norm)

    ax.plot(x, rc_norm, linewidth=1, color='k')
    ax.plot(ax.get_xlim(), (1,1), linestyle='--', color='k')
    ax.set_ylim(bottom=0)
    if ax.get_ylim()[1]<2:
        ax.set_ylim(top=2)
    if top_limit:
        ax.set_ylim(top=top_limit)
    ax.set_ylabel(r"$RC_{norm}$")#+"\n (normalized rich club coefficient)")
    ax.fill_between(x, 1,rc_norm, where=(rc>upper95) | (rc<lower95), **kwargs)
    
    return ax

def graph_from_sparse(data, directed=None):
    from igraph import Graph
    sources, targets = data.nonzero()
    
    if directed==None:
        from numpy import all
        directed = not all(data[sources, targets]==data[targets, sources])
    from numpy import array
    g = Graph(zip(sources, targets), directed=directed, edge_attrs={'weight': array(data[sources, targets])[0]})
    if g.is_directed():
        return g
    else:
        return g.simplify(combine_edges="first")
