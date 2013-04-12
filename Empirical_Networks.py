# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import richclub
from igraph import Graph
import pickle

import matplotlib as mpl
mpl.use('Agg')
from matplotlib.pyplot import *
from pylab import *

# <codecell>

def preserve_strength(G):
    from numpy.random import permutation
    g = G.copy()
    degree_sequence = g.strength()
    strength_sequence = array(g.strength(weights='weight'))
    
    
    for k in unique(degree_sequence):
        k_graph = g.vs.select(_degree_eq=k)
        k_ind = k_graph.indices
        strength_sequence[k_ind] = strength_sequence[permutation(k_ind)]
    
    mean_k = mean(degree_sequence)
    mean_s = mean(strength_sequence)
    
    for e in g.es:
        e["weight"] = (mean_k/mean_s) * strength_sequence[e.source]*strength_sequence[e.target] / (degree_sequence[e.source]*degree_sequence[e.target])
    return g

# <codecell>

def preserve_strength_weights(G):
    g = G.copy()
    degree_sequence = g.strength()
    strength_sequence = array(g.strength(weights='weight'))
    weight_sequence = sort(array(g.es['weight']))
    
    assigned_strength_sequence = zeros(g.vcount())
    assigned_weight_sequence = zeros(g.ecount())
    
    def initial_expected_weight(ind):
        e = g.es[int(ind)]
        return 1-1e-9/(strength_sequence[e.source]*strength_sequence[e.target])
 
    initial_expected_weight = vectorize(initial_expected_weight)
    
    expected_weight_sequence = initial_expected_weight(range(g.ecount()))

    #Sort the expected weight sequence, but retain the original identity of each link
    indices = argsort(expected_weight_sequence) 
    expected_weight_sequence = expected_weight_sequence[indices]

    #Get the sources and targets of the links, in index
    sources = array([e.source for e in g.es])
    targets = array([e.target for e in g.es])
    
    for i in range(g.ecount()):
        rank = randint(g.ecount()-i)
        ind = indices[rank]
        weight = weight_sequence[ind]
        
        indices = delete(indices, rank)
        expected_weight_sequence = delete(expected_weight_sequence, rank)
        
        source = sources[ind]
        target = targets[ind]
        source_strength = strength_sequence[source]
        target_strength = strength_sequence[target]
        
        assigned_weight_sequence[ind] = weight
        assigned_strength_sequence[source] += weight
        assigned_strength_sequence[target] += weight
        
        #weight_expectation_changed = where((sources==source) | (sources==target) | (targets==source) |  (targets==target))[0]
        #Could try to do the above line with g.es.select. Might be slow.
        #weight_expectation_changed = (sources==source) | (sources==target) | (targets==source) |  (targets==target)
        attached_to_source = (sources==source) | (targets==source)
        attached_to_target = (sources==target) | (targets==target)
        
        for j in range(g.ecount()):
            if attached_to_source[j]:
                rank = where(indices==j)[0]
                expected_weight_sequence[rank] -= weight/source_strength
            if attached_to_target[j]:
                rank = where(indices==j)[0]
                expected_weight_sequence[rank] -= weight/target_strength
        
        #for i in weight_expectation_changed:
        #for j in range(g.ecount()):
        #    if weight_expectation_changed[j]:
        #        rank = where(indices==ind)[0]
        #        expected_weight_sequence[rank] = expected_weight(int(ind))
        
        #expected_weight_sequence = expected_weight(indices)
        sort_ind = argsort(expected_weight_sequence)
        expected_weight_sequence = expected_weight_sequence[sort_ind]
        indices = indices[sort_ind]
        
    g.es["weight"] = assigned_weight_sequence
    return g
    
    

# <codecell>

def preserve_strength_weights(G, steps_to_recalculate=10):
    g = G.copy()
    strength_sequence = array(g.strength(weights='weight'))
    weight_sequence = sort(array(g.es['weight']))
    
    assigned_strength_sequence = zeros(g.vcount())
    assigned_weight_sequence = zeros(g.ecount())
    
    def initial_expected_weight(ind):
        e = g.es[int(ind)]
        return 1 - 1e-9/(strength_sequence[e.source]*strength_sequence[e.target])
    initial_expected_weight = vectorize(initial_expected_weight)
    expected_weight_sequence = initial_expected_weight(range(g.ecount()))
    
    #Sort the expected weight sequence, but retain the original identity of each link
    indices = argsort(expected_weight_sequence) 
    expected_weight_sequence = expected_weight_sequence[indices]

    #Get the sources and targets of the links, in index order
    sources = array([e.source for e in g.es])
    targets = array([e.target for e in g.es])
    
    for i in range(0, g.ecount(), steps_to_recalculate):
        ranks = argsort(expected_weight_sequence[indices], kind='mergesort')        
        rs = randint(0, g.ecount()-i, steps_to_recalculate)
        
        for s in range(steps_to_recalculate):
            r = rs[s]
            rank = ranks[r]
            ind = indices[r]
            
            source = sources[ind]
            target = targets[ind]
            source_strength = strength_sequence[source]
            target_strength = strength_sequence[target]
                
            assigned_weight = weight_sequence[rank]
            assigned_weight_sequence[ind] = assigned_weight
            
            expected_weight_sequence[(sources==source) | (targets==source)] -= assigned_weight/source_strength
            expected_weight_sequence[(sources==target) | (targets==target)] -= assigned_weight/target_strength
            
        expected_weight_sequence[indices] = expected_weight_sequence[indices][ranks]
        sources[indices] = sources[indices][ranks]
        targets[indices] = targets[indices][ranks]
        indices = delete(indices, rs)
        weight_sequence = delete(weight_sequence, rs)
    g.es["weight"] = assigned_weight_sequence
    return g

# <codecell>

def plot_rc(g, controls, ax, **kwargs):
    
    rc = richclub.rich_club_coefficient(g, **kwargs)
    
    rich_club_coefficient = vectorize(richclub.rich_club_coefficient)
    
    rc_controls = []
    for c in controls:
        rc_controls.append(richclub.rich_club_coefficient(c, **kwargs))
    rc_controls = array(rc_controls)
    control_mean = mean(rc_controls, axis=0)
    
    if 'rank' in kwargs:
        x = kwargs['rank']
    else:
        if 'mode' in kwargs:
            mode = kwargs['mode']
        else:
            mode = 'percentile'
        if mode=='percentile':
            x = arange(10.0,100.0, 10.0)
        else:
            x = unique(g.degree())
            
    from scipy.stats import scoreatpercentile
    upper95 = array([scoreatpercentile(rc_controls[:,i], 95) for i in range(len(x))])
    lower95 = array([scoreatpercentile(rc_controls[:,i], 5) for i in range(len(x))])
    
    ax.plot(x, rc, color='r')
    ax.plot(x, control_mean, color='b')
    ax.fill_between(x, lower95, upper95, where=(lower95>0)&(upper95>0), alpha=.1, color='b')
    ax.set_yscale('log')
    #ax.set_xscale('log')
    
    rc_norm = rc/control_mean
    ax2 = ax.twinx()
    ax2.plot(x, rc_norm, linewidth=4, color='k')
    ax2.plot(ax.get_xlim(), (1,1), linestyle='--', color='k')
    ax2.set_ylim(bottom=0)
    if ax2.get_ylim()[1]<2:
        ax2.set_ylim(top=2)
    ax2.set_ylabel(r"$\phi_{norm}$, normalized rich club coefficient")
    ax2.fill_between(x, 1,rc_norm, where=(rc>upper95) | (rc<lower95), alpha=.05, color='g')
    
    return ax, rc, lower95, upper95

# <codecell>

def rich_club_analysis(g, n_samples=10, weight_controls=None, topology_controls=None, **kwargs):

    g=g.copy()
    
    if weight_controls==None:
        print "Creating weight controls"
        weight_controls = array([preserve_strength(g) for i in range(n_samples)])
    if topology_controls==None:
        print "Creating topology controls"
        topology_controls = []
        for i in range(n_samples):
            topology_controls.append(Graph.Degree_Sequence(g.degree(), method='vl'))
            #c = g.copy()
            #topology_controls.append(c.rewire(c.ecount()*5))
        topology_controls = array(topology_controls)
        
    ######
    f = figure(figsize=(8,11))
    
    ax = f.add_subplot(311)
    g_uw = g.copy()
    g_uw.es["weight"] = ones(g_uw.ecount())
    
    ax, rc_norm, lower95, upper95 = plot_rc(g_uw, topology_controls, ax, mode='raw', **kwargs)
    ax.set_title("Topological Rich Club")
    ax.set_ylabel(r"$C$, total connections in rich club")
    ax.set_xlabel("Degree")
    
    ax = f.add_subplot(312)
    ax, a,b,c = plot_rc(g, weight_controls, ax, **kwargs)
    ax.set_title("Weight Rich Club, Preserve Topology and Strengths")
    ax.set_ylabel(r"$C$, total weight in rich club")
    ax.set_xlabel("Strength Percentile")
    
    ax = f.add_subplot(313)
    for i in range(n_samples):
        ax.scatter(sort(g.strength(weights='weight')), sort(weight_controls[i].strength(weights='weight')))
    ax.set_title("Strength Sequence Preservation in Random Controls")
    ax.set_xlabel("Empirical Strength Sequence")
    ax.set_ylabel("Randomized Strength Sequence")
    
    f.subplots_adjust(hspace=.5)
    
    return f, weight_controls, topology_controls

# <codecell>

#full = Graph.Full(100)
#full.es["weight"] = rand(full.ecount())
#g = preserve_strength_weights(full, steps_to_recalculate=1)
#scatter(sort(full.strength(weights='weight')), sort(g.strength(weights='weight')))
#xlabel('Original')
#ylabel('Randomized')
#gca().set_aspect('equal')
#figure()
#scatter(sort(full.es["weight"]), sort(g.es['weight']))
#gca().set_aspect('equal')

# <codecell>

filename = "Wodehouse.txt"
print filename
g = Graph.Read(filename, format='pajek')
f, weight_controls, topology_controls = rich_club_analysis(g, n_samples=50)
savefig(filename[:-4]+'.pdf')
nets = {'empirical': g, 'weight_controls': weight_controls, 'topology_controls': topology_controls}
pickle.dump(nets, open(filename[:-4]+".p", "wb"))

# <codecell>

filename = "Wodehouse_reweighted.txt"
print filename
g = Graph.Read(filename, format='pajek')
f, weight_controls, topology_controls = rich_club_analysis(g, n_samples=50)
savefig(filename[:-4]+'.pdf')
nets = {'empirical': g, 'weight_controls': weight_controls, 'topology_controls': topology_controls}
pickle.dump(nets, open(filename[:-4]+".p", "wb"))

# <codecell>

filename = "openflights.txt"
print filename
g = Graph.Read(filename, format='ncol')
f, weight_controls, topology_controls = rich_club_analysis(g, n_samples=50)
savefig(filename[:-4]+'.pdf')
nets = {'empirical': g, 'weight_controls': weight_controls, 'topology_controls': topology_controls}
pickle.dump(nets, open(filename[:-4]+".p", "wb"))

US_2010 = Graph.Read("USairport_2010.txt", format='ncol')

# <codecell>

filename = "USairport_2010.txt"
print filename
g = Graph.Read(filename, format='ncol')
f, weight_controls, topology_controls = rich_club_analysis(g, n_samples=50)
savefig(filename[:-4]+'.pdf')
nets = {'empirical': g, 'weight_controls': weight_controls, 'topology_controls': topology_controls}
pickle.dump(nets, open(filename[:-4]+".p", "wb"))

# <codecell>

filename = "cond-mat.gml"
print filename
g = Graph.Read(filename, format='gml')
g.es['weight'] = SCN.es['value']
g.delete_vertices(SCN.vs.select(_degree_eq=0))

f, weight_controls, topology_controls = rich_club_analysis(g, n_samples=50)
savefig(filename[:-4]+'.pdf')
nets = {'empirical': g, 'weight_controls': weight_controls, 'topology_controls': topology_controls}
pickle.dump(nets, open(filename[:-4]+".p", "wb"))

# <codecell>

filename = 'OClinks_w.txt'
print filename
g = Graph.Read(filename, format='ncol')
f, weight_controls, topology_controls = rich_club_analysis(g, n_samples=50)
savefig(filename[:-4]+'.pdf')
nets = {'empirical': g, 'weight_controls': weight_controls, 'topology_controls': topology_controls}
pickle.dump(nets, open(filename[:-4]+".p", "wb"))

# <codecell>

filename = 'taskfmri.mat'
print filename
from scipy.io import loadmat
mat = loadmat(filename)
g = Graph.Weighted_Adjacency(mat['Adj_Matrix'].tolist(), mode=1)
f, weight_controls, topology_controls = rich_club_analysis(g, n_samples=50)
savefig(filename[:-4]+'.pdf')
nets = {'empirical': g, 'weight_controls': weight_controls, 'topology_controls': topology_controls}
pickle.dump(nets, open(filename[:-4]+".p", "wb"))

# <codecell>

mat = loadmat('/data/alstottjd/BMS_Study_20120826.mat')
for i in mat.keys():
    if 'placebo' in i and 'C' in i:
        print i

g = Graph.Weighted_Adjacency(mat['C01_placebo'].tolist(), mode=1)
filename = 'C01_placebo.mat'
print filename
f, weight_controls, topology_controls = rich_club_analysis(g, n_samples=50)
savefig(filename[:-4]+'.pdf')
nets = {'empirical': g, 'weight_controls': weight_controls, 'topology_controls': topology_controls}
pickle.dump(nets, open(filename[:-4]+".p", "wb"))

