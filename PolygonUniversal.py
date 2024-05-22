from pygeos import polygons,delaunay_triangles, get_coordinates
import random
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sys
import time
import sage.all
from sage.graphs.graph import Graph



def longest_list(l):
	res = l[0]
	for elem in l:
		if len(elem)>len(res):
			res = elem
	return res

def random_point_in_triangle(t):
    """
    return random point within a given triangle
    """
    r1 = random.random()
    r2 = random.random()
    if r1+r2>1:
        r1 = 1-r1
        r2 = 1-r2
    a = 1-r1-r2
    x= a*t[0][0] + r1*t[1][0] + r2*t[2][0]
    y= a*t[0][1] + r1*t[1][1] + r2*t[2][1]
    return (x,y)
def random_point_on_edge(e,low = 1/4,high=3/4):
    """
    return random point on an edge
    """
    r1 = random.uniform(low,high)
    x= (1-r1)*e[0][0]+r1*e[1][0]
    y= (1-r1)*e[0][1]+r1*e[1][1]
    return (x,y)

def neighbours(adj,i):
    return [j for j in range(len(adj[i])) if adj[i][j]==1]


class Polygon_Universal:
    def __init__(self,adj,cycle,pos,debug=False):
        self.debug = debug
        self.pos = pos
        self.cycle = cycle
        self.polygon= polygons(pos)
        self.adj = adj
        self.n = len(pos)
        triang = [cycle[pos.index(tuple(i))] for i in get_coordinates(delaunay_triangles(self.polygon)).tolist()]
        self.triangulation = []
        for i in range(0,len(triang),4):
            tmp = [triang[i],triang[i+1],triang[i+2]]
            self.triangulation.append(sorted(tmp))
        #create dual of triangulation
        self.dual = {k:[] for k in range(len(self.triangulation))}
        # for any edge sorted index tuple, give the indices of both triangles containing it
        self.edge_triangle = dict()
        # for any 2 triangles in a sorted tuple, gives their shared edge, sorted by indices
        self.triangle_edge = dict()
        for i in range(len(self.triangulation)):
            for j in range(i,len(self.triangulation)):
                # if combined set ==2 then shared edge
                tmp = set(self.triangulation[i]) & set(self.triangulation[j])
                if len(tmp) == 2:
                    self.dual[i].append(j)
                    self.dual[j].append(i)
                    tmp = tuple(sorted(tmp))
                    self.edge_triangle[tmp] = (i,j)
                    self.triangle_edge[(i,j)] = tmp
                    self.triangle_edge[(j,i)] = tmp

    def intersection(self,*elements):
        """
        compute intersection between 2 geometric elements of the graph triangulation
        input element format:
        (2,x) - triangle index x
        (1,(a,b)) - edge connecting vertices a and b
        (0,a) - point index a
        return intersection in same format
        """
        intersection = elements[0]
        for elem in elements[1:]:
            if elem == intersection:
                continue
            #swap if elem bigger than intersection to reduce number of conditions to handle
            if elem[0]>intersection[0]:
                elem,intersection=intersection,elem
            if intersection[0] == 2: #triangle
                if elem[0]==2: #triangle with triangle -> can share an edge or point
                    shared = set(self.triangulation[elem[1]]) & set(self.triangulation[intersection[1]]) 
                    if elem[1] in self.dual[intersection[1]]:
                        intersection = (1,self.triangle_edge[(intersection[1],elem[1])])
                    elif len(shared) == 1:
                        intersection = (0,shared.pop())
                    else:
                        return None
                elif elem[0]==1: #triangle with edge
                    if len(set(elem[1]) & set(self.triangulation[intersection[1]])) ==2:
                        intersection = elem
                    else:
                        return None
                else:
                    if elem[1] in self.triangulation[intersection[1]]: #point with triangle
                        intersection = elem
                    else:
                        return None
            elif intersection[0] == 1: #edge
                if elem[0]==1: #two edges
                    if elem[1][0] in intersection[1]:
                        intersection = (0,elem[1][0])
                    elif elem[1][1] in intersection[1]:
                        intersection = (0,elem[1][1])
                    else:
                        return None
                else: #point and edge
                    if elem[1] in intersection[1]:
                        intersection = elem
                    else:
                        return None
            else: #points that are not equal
                return None
        return intersection 


    def recursive_traversal(self,e_q,previous_triangle):
        """
        create pocket by assigning all points of the graph to positions inside the given pocket bounded by e_q
        return list of tuples with simplex assigned to each vertex of the graph
        """
        #trivial pocket
        if abs(self.cycle.index(e_q[0])-self.cycle.index(e_q[1]))%(self.n-2) == 1:
            if self.debug:
                print(f"trivial {e_q}")
            res=[]
            for i in range(len(self.adj)):
                if i == e_q[0] or i==e_q[1]:
                    res.append((0,i))
                else: #directly restrict pocket by assigning to t_Q+ instead of e_Q
                    res.append((2,previous_triangle))
            return res
        #non trivial pocket
        else:
            if self.debug:
                print(f"through edge {e_q}")
            #find correct triangle by finding shared edge e_q in dual
            current = self.edge_triangle[e_q][0]
            if current == previous_triangle:
                current = self.edge_triangle[e_q][1]
            #find both subpockets
            p_m = list(set(self.triangulation[current])-set(e_q))[0]
            subpockets = [tuple(sorted([e_q[0],p_m])),tuple(sorted([e_q[1],p_m]))]
            #add trivial pocket edges if they exist
            for i in range(3):
                if abs(self.triangulation[current][(i+1)%3]-self.triangulation[current][i])%(self.n-2) ==1:
                    subpockets.append(tuple(sorted([self.triangulation[current][(i+1)%3],self.triangulation[current][i]])))
            #recursive call on both subpockets of current triangle
            t_l = self.recursive_traversal(subpockets[0],current)
            t_r = self.recursive_traversal(subpockets[1],current)
            if t_l is None or t_r is None:
                return None
            #combine both pockets
            if self.debug:
                print(f"L^+{t_l}")
                print(f"R^+{t_r}")
            res = []
            for i in range(len(self.adj)):
                #check for intersections
                intersection = self.intersection(t_l[i],t_r[i])
                if self.debug:
                    print(f"intersection{t_l[i]}{t_r[i]}:{intersection}")
                if intersection is not None:
                    res += [intersection]
                #take restrictive if other in t_outerq 
                elif t_r[i]==(2,current):
                    res += [t_l[i]]
                elif t_l[i]==(2,current):
                    res += [t_r[i]]
                
                else:
                    return None
            if self.debug:
                print(f"before restriction:{res}")
            #restrict pocket
            coord = []
            for i in range(len(self.adj)):
                tmp = [self.intersection(res[k],(1,e_q)) is not None for k in neighbours(self.adj,i)]
                if self.debug:
                    print(f"neighbours {i}:{neighbours(self.adj,i)}")
                    print(f"intersection{(1,e_q)}:{tmp}")
                if (res[i] == (1,e_q) or res[i] == (2,current)) and all(tmp):
                    coord += [(2,previous_triangle)]
                else:
                    intersection = self.intersection(res[i],(1,e_q))
                    if intersection is not None:
                        coord += [intersection]
                    else:
                        coord+=[res[i]]
            res = coord
            if self.debug:
                print(f"after restriction:{res}")
            return res

    def polygon_universal(self,measure = False):
        """
        Compute the embedding of a graph in polygon following 
        the polygon universal algorithm
        @polygon: the given polygon as a list of tuple coordinates
        @graph: adjacency matrix of the graph, first n points are the polygonal cycle
        """

        #pick random triangle of the triangulation as t_root
        begin = time.time()
        t_root = random.randint(0,len(self.triangulation)-1)
        if self.debug:
            print(f"t_root{self.triangulation[t_root]}")
            print("pocket 1")
        p1 = self.recursive_traversal((self.triangulation[t_root][0],self.triangulation[t_root][1]),t_root)
        if self.debug:
            print("pocket 2")
        p2 = self.recursive_traversal((self.triangulation[t_root][1],self.triangulation[t_root][2]),t_root)
        if self.debug:
            print("pocket 3")
        p3 = self.recursive_traversal((self.triangulation[t_root][0],self.triangulation[t_root][2]),t_root)
        if self.debug:
            print(f"pocket 1:{p1}")
            print(f"pocket 2:{p2}")
            print(f"pocket 3:{p3}")
        #combine 3 pockets
        res=  []
        if p1 is not None and p2 is not None and p3 is not None:
            for i in range(len(self.adj)):
                #compute intersection between 3 pockets
                intersection = self.intersection(p1[i],p2[i],p3[i])
                if intersection is not None:
                    res+=[intersection]
                # 2 pockets in t_root and one not
                elif p1[i] == (2,t_root) and p1[i] == p2[i] and p3[i] !=p1[i]:
                    res+= [p3[i]]
                elif p1[i] == (2,t_root) and p1[i] == p3[i] and p1[i]!=p2[i]:
                    res+=[p2[i]]
                elif p2[i] == (2,t_root) and p2[i] == p3[i] and p1[i]!=p2[i]:
                    res+=[p1[i]]
                else:
                    print("No valid sketch")
                    if measure:
                        return time.time()-begin
                    return None
        else:
            print("No valid sketch")
            if measure:
                return time.time()-begin
            return None
        if debug:
            print(res)
        coord = []
        #assign real positions to points
        for i in range(len(res)):
            if res[i][0] == 0:
                coord += [self.pos[self.cycle.index(res[i][1])]]
            elif res[i][0] == 1: #simplex is an edge -> choose random point on edge
                coord += [random_point_on_edge([self.pos[self.cycle.index(k)] for k in res[i][1]])]
            elif res[i][0] == 2: #simplex is a triangle -> choose random point in triangle
                coord += [random_point_in_triangle([self.pos[self.cycle.index(k)] for k in self.triangulation[res[i][1]]])]
        if measure:
            return time.time()-begin
        return coord





if len(sys.argv) < 3:
    print("PolygonUniversal.py [0-7] debug[0-1]:\n")
else:
    case = int(sys.argv[1])
    debug = bool(int(sys.argv[2]))
    match case:
        case 0:
            adj = np.array([
                [0, 1, 0, 0, 1, 1, 0, 0, 0, 0],
                [1, 0, 1, 0, 0, 0, 0, 0, 1, 0],
                [0, 1, 0, 1, 0, 0, 1, 0, 0, 0],
                [0, 0, 1, 0, 1, 0, 0, 1, 0, 0],
                [1, 0, 0, 1, 0, 0, 0, 0, 0, 1],
                [1, 0, 0, 0, 0, 0, 1, 1, 0, 0],
                [0, 0, 1, 0, 0, 1, 0, 1, 0, 0],
                [0, 0, 0, 1, 0, 1, 1, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0]])
            pos = [(0,3),(2,0),(5,0),(7,3),(3.5,5)]
            cycle = [0,1,2,3,4]
        case 1:
            adj = np.array(
                [[0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], 
            [0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], 
            [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 
            [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], 
            [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0], 
            [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0], 
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0], 
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1], 
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]])
            pos = [(0,0), (4,-3), (9,-4),(13,-4),(14,0),(13,5),(9,9),(6,9),(2,8),(1,5)]
            cycle = [0,1,2,3,4,5,6,7,8,9]
        case 2:
            adj = np.array([[0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
			[1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
			[1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
			[0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
			[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0],
			[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],
			[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1],
			[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0]])
            pos = [(3,0), (1,5), (6,8),(11,5),(9,0)]
            cycle = [0,1,2,3,4]
        case 3:
            adj = np.array(
			[
			[0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0]])

            pos = [(3,0), (1,5), (6,8),(11,5),(9,0)]
            cycle = [0,1,2,3,4]
        case 4:
            adj = np.array([
			[0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0]])

            pos = [(0,0), (4,-3), (9,-4),(13,-4),(14,0),(13,5),(9,9),(6,9),(2,8),(1,5)]
            cycle = [0,1,2,3,4,5,6,7,8,9]
        case 5:
            if len(sys.argv)==4:
                file = open(sys.argv[3])
            else:
                print("PolygonUniversal 5 debug graph_file(graph6 format)")
                quit()
            #choose random graph from file
            line = next(file)
            for num, aline in enumerate(file, 2):
                if random.randrange(num):
                    continue
                line = aline
            g = nx.from_graph6_bytes(line.strip().encode())
            adj = nx.to_numpy_array(g)
            #choose a face of the graph as a peripheral circuit by using sage.faces method on a combinatorial embedding of the graph
            tmp = nx.check_planarity(g)[1] #create a planar embedding
            comb_embedding = {node:list(tmp.neighbors_cw_order(node)) for node in tmp.nodes()}
            tmp = Graph(comb_embedding)
            faces = tmp.faces(comb_embedding) #find faces
            #choose largest face as boundary
            cycle = [e[0] for e in longest_list(faces)]
            #assign position of the points of the peripheral circuit
            tmp_pos = nx.spectral_layout(nx.cycle_graph(len(cycle)))
            pos = [tuple(tmp_pos[i]) for i in range(len(cycle))]
            print(cycle)
            print(pos)
        case 6:
            if len(sys.argv)==5:
                file = sys.argv[3]
                n = int(sys.argv[4])
            else:
                print("PolygonUniversal.py 6 debug[0-1] graph_file(graph6 format) n_graphs:int")
                quit()
            times = []
            with open(file) as f:
                i=0
                while i<n:
                    i+=1
                    g = nx.from_graph6_bytes(f.readline().strip().encode())
                    adj = nx.to_numpy_array(g)
                    #choose a face of the graph as a peripheral circuit by using sage.faces method on a combinatorial embedding of the graph
                    tmp = nx.check_planarity(g)[1] #create a planar embedding
                    comb_embedding = {node:list(tmp.neighbors_cw_order(node)) for node in tmp.nodes()}
                    tmp = Graph(comb_embedding)
                    faces = tmp.faces(comb_embedding) #find faces
                    #choose largest face as boundary
                    cycle = [e[0] for e in longest_list(faces)]
                    #assign position of the points of the peripheral circuit
                    if len(cycle)<4:
                        n+=1
                    else:
                        tmp_pos = nx.spectral_layout(nx.cycle_graph(len(cycle)))
                        pos = [tuple(tmp_pos[i]) for i in range(len(cycle))]
                        #measure execution time
                        embedding = Polygon_Universal(adj,cycle,pos)
                        t = embedding.polygon_universal(measure=True)
                        times.append(t)
            print(np.mean(times))
    if case<6:
        g = nx.from_numpy_array(adj)
        embedding = Polygon_Universal(adj,cycle,pos,debug)
        #color the edges of the graph
        for u,v in g.edges():
            g[u][v]["color"] = "black"
        #color the triangulation edges
        if debug:
            for triangle in embedding.triangulation:
                g.add_edge(triangle[0],triangle[1],color="red")
                g.add_edge(triangle[1],triangle[2],color="red")
                g.add_edge(triangle[2],triangle[0],color="red")
        pos = embedding.polygon_universal()
        if pos is not None:
            fig, ax = plt.subplots()
            ax.set(xlim=(-1,5),ylim=(-1,5))
            for u in g.nodes(): #vertex labels
                plt.text(pos[u][0],pos[u][1]+0.1,s=str(u))
            nx.draw_networkx(g,pos,ax=ax,with_labels=False,node_size=50,edge_color = [g[u][v]["color"] for u,v in g.edges()])
            ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
            plt.show()

