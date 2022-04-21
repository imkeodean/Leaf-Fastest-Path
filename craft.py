import math
from itertools import combinations_with_replacement
from collections import defaultdict


SET_SIZE = 8
NUM_SHARDS = 17
BASE_CRIT = 0.346
tiers = ("A", "S", "B", "M", "H")
craft_order = {"A":0, "S":1, "B":2, "M":3, "H":4}
craft_bonus = {"A":256, "S":384, "B":512, "M":640, "H":768}
fusion_req = {"A":0, "S":5, "B":37, "M":233, "H":1417}
ascension_req = {"A":0, "S":6, "B":8, "M":10, "H":14}


# Used to calculate the stat bonuses on a given leaf
def power(tier, ascension, shards):
		return (craft_bonus[tier]**1.6)*(12*(ascension+1)+7)*(1+3*shards)


class CraftSet:

	def __init__(self, leaf_list):
		self.leaf_list = tuple(sorted(leaf_list))
		total_wem = 0
		total_crit = 0
		for leaf in leaf_list:
			total_wem += leaf.wem
			total_crit += leaf.crit
		self.wem = 1+total_wem
		self.epk = math.floor(self.wem*18)
		self.crit = BASE_CRIT+total_crit
		self.DErate = (self.epk/30)*5*(1+self.crit)*(3600/144)

	def __eq__(self, other):
		return self.leaf_list == other.leaf_list

	def __hash__(self):
		return hash(self.leaf_list)

	def __repr__(self):
		return str(list(self.leaf_list))

	def parents(self):
		return self.min_replace[self.leaf_list]

	def time_to_new_leaf(self, leaf, new_craft): # Assumes if not a new craft, then the new leaf only has to be ascended once from an old leaf
		if new_craft:
			cost = (50/(1+self.crit))*fusion_req[leaf.tier]+(100/(1+self.crit))*ascension_req[leaf.tier]*((leaf.ascension)*(leaf.ascension+1)/2)
		else:
			cost = (100/(1+self.crit))*ascension_req[leaf.tier]*leaf.ascension
		return cost/self.DErate

	def time_to_max_hemas(self):
		leaf_set = list(self.leaf_list)
		min_leaf = self.leaf_list[0]
		max_hema = Leaf("H", 10)
		time = 0
		while min_leaf != max_hema:
			for i in range(11):
				if min_leaf.wem < Leaf("H", i).wem:
					new_hema = Leaf("H", i)
					break
			time += CraftSet(leaf_set).time_to_new_leaf(new_hema, True)
			for i in range(new_hema.ascension, 10):
				new_set = CraftSet([Leaf("H", i)]+leaf_set[1:])
				time += new_set.time_to_new_leaf(Leaf("H", i+1), False)
			leaf_set = leaf_set[1:]+[max_hema]
			min_leaf = leaf_set[0]
		return time




class Leaf:

	def __init__(self, tier, ascension):
		self.tier = tier
		self.ascension = ascension
		self.wem = power(tier, ascension, 10)*10**(-7.5)
		self.crit = power(tier, ascension, NUM_SHARDS-10)*10**(-9)

	def __eq__(self, other):
		return self.tier == other.tier and self.ascension == other.ascension

	def __lt__(self, other):
		if self.tier != other.tier:
			return craft_order[self.tier] < craft_order[other.tier]
		else:
			return self.ascension < other.ascension

	def __le__(self, other):
		return self == other or self < other

	def __hash__(self):
		return hash((self.tier, self.ascension))

	def __repr__(self):
		return self.tier + str(self.ascension)

	def parents(self):
		return self.min_replace[(self.tier, self.ascension)]


end = CraftSet([Leaf("H",10)]*SET_SIZE)


"""
Creates list of possible leaf crafts
Only Ancient leaf is A10
Sacreds are not worth crafting, can make better Biotites for cheaper
Biotite starts replacing A10 at B3
Malachite starts replacing A10 at M2
Hematite starts replacing A10 at H1
"""
leaf_pool = [Leaf("A", 10)]
for i in range(3, 11):
	leaf_pool.append(Leaf("B", i))
for i in range(2, 11):
	leaf_pool.append(Leaf("M", i))
#for i in range(1, 11):
#	leaf_pool.append(Leaf("H", i))


"""
Maps each leaf in leaf_pool to the list of leafs that can replace it.
A leaf can be replaced by a leaf of the same rarity one ascension level
higher, or the lowest ascension level of a higher rarity leaf that is
an imporvement.
"""
leaf_min_replace = {("H", 10):()}
for leaf1 in leaf_pool:
	replace_list  = []
	tier_flags = dict(zip(tiers, [False]*5))
	for leaf2 in leaf_pool:
		if leaf1 < leaf2 and leaf1.wem < leaf2.wem and not tier_flags[leaf2.tier]:
			replace_list.append(leaf2)
			tier_flags[leaf2.tier] = True
	leaf_min_replace[(leaf1.tier, leaf1.ascension)] = tuple(replace_list)
Leaf.min_replace = leaf_min_replace


"""
Creates list of possible leaf sets using leafs from leaf_pool
"""
combs = tuple(combinations_with_replacement(leaf_pool, SET_SIZE))
print("Combs done")
set_pool = []
for i in combs:
	set_pool.append(CraftSet(i))
set_pool.append(end)
set_pool = tuple(set_pool)
print("Craft set pool done")


"""
Maps each leaf set in set_pool to the list of sets that can replace it and the
time it takes to reach those sets. Replacing any leaf in a leaf set with something
from that leaf's min_replace list creates a valid replacement set.
"""
set_min_replace = {}
for set1 in set_pool:
	replace_list = []
	processed = []
	for i in range(len(set1.leaf_list)):
		replaced_leaf = set1.leaf_list[i]
		if replaced_leaf not in processed:
			processed.append(replaced_leaf)
			for j in replaced_leaf.parents():
				set2_list = list(set1.leaf_list)
				set2_list[i] = j
				replace_list.append((CraftSet(set2_list), set1.time_to_new_leaf(j, replaced_leaf.tier != j.tier)))
	set_min_replace[set1.leaf_list] = tuple(replace_list)
CraftSet.min_replace = set_min_replace
print("Craft set parents done")


"""
Finds single source shortest paths for
Directed Acyclic Graphs with complexity O(V+E)
"""
class Graph:

    def __init__(self,vertices):
        self.V = vertices
        self.graph = defaultdict(list)
 
    def addEdge(self,u,v,w):
        self.graph[u].append((v,w))

    # A recursive function that sorts the graph and stores it in stack
    def topologicalSortUtil(self,v,visited,stack):
        visited[v] = True
        if v in self.graph.keys():
            for node,weight in self.graph[v]:
                if visited[node] == False:
                    self.topologicalSortUtil(node,visited,stack)
        stack.append(v)

    # Finds shortest path from s to t
    def shortestPath(self, s, t):
        visited = dict(zip(self.V, [False]*len(self.V)))
        stack = []
        self.topologicalSortUtil(s,visited,stack)

        # The path dict will map a vertex v to its shortest distance from s,
        # and the previous node in the shortest path to v from s
        path = {}
        for v in self.V:
        	path[v] = (float("Inf"), v)
        path[s] = (0, s)
 
        # Process vertices in topological order
        while stack:
            i = stack.pop()
            for node,weight in self.graph[i]:
                if path[node][0] > path[i][0] + weight:
                    path[node] = (path[i][0] + weight, i)

        # Generate shortest path from s to t
        path_list = []
        node = t
        while node != s:
        	path_list.append((node, round(path[node][0],2)))
        	node = path[node][1]
        path_list.append((s, 0))
        for i in path_list[::-1]:
        	print(i)
        return path


"""
This sets up and solves the problem
"""
set_graph = Graph(set_pool)
for u in set_pool:
	set_graph.addEdge(u, set_pool[-1], u.time_to_max_hemas())
	for v in u.parents():
		set_graph.addEdge(u, v[0], v[1])
print("Graph done")
path = set_graph.shortestPath(set_pool[0], set_pool[-1])

