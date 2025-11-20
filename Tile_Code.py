from fontTools.misc.transform import Offset
from manim import *
import networkx as nx
from gcol import edge_coloring
from collections import defaultdict

'''
TODO:
1. Boundary Tiles/Global Grid
    a) Can specify coord on Global Grid and X/Z, will show stabilizer on that tile
    b) Will 'cutoff' Tile/Grid Objects correctly at boundaries
2. Intersection between Grids (Conflict Tree)
    a) Should make a new Grid Object storing conflict
3. Partitioning Grid Objects (bipartite)
4. Coloring Grid Object 
'''
class Grid(Graph):
    def __init__(self,tile,make_stab=True, highlightedges = [], with_back = True, **kwargs):
        self.tile = tile
        vertices = []
        edges = []
        layout = {}
        self.stabilizers = []
        self.highedge = highlightedges
        self.offset = tile.offset
        self.back = with_back
        #use first tile.size^2 entries for vertical, next tile.size^2 for horizontal
        if with_back == True:
            for i in range(tile.size):
                for j in range(tile.size):
                    vertices.append((i,j))
                    vertices.append((i+tile.size,j+tile.size))
                    edges.append(((i,j),(i+tile.size,j+tile.size)))
                    if(j+1<tile.size):
                        edges.append(((i,j),(i+tile.size,j+tile.size+1)))
                    if(i+1<tile.size):
                        edges.append(((i+1,j),(i+tile.size,j+tile.size)) )
                    if(j+1<tile.size and i+1<tile.size):
                        edges.append(((i+1,j),(i+tile.size,j+tile.size+1)))
        else:
            edges = highlightedges

            for hor,ver in edges:
                if vertices.__contains__(hor) == False:
                    vertices.append(hor)
                if vertices.__contains__(ver) == False:
                    vertices.append(ver)

        curpoint = (0,0)
        for point in tile.gridArrayVer:

            ref = tile.gridArrayVer[point]
            layout[curpoint] = [ref[0],ref[1],0]

            horpoint = (curpoint[0]+tile.size,curpoint[1]+tile.size)
            ref = tile.gridArrayHor[curpoint]

            layout[horpoint] = [ref[0],ref[1],0]
            #self.add(Line(self[curpoint],self[horpoint]))

            curpoint = (curpoint[0],curpoint[1]+1)
            if(curpoint[1] == tile.size):
                curpoint = (curpoint[0]+1,0)

        super().__init__(vertices=vertices,edges = edges,layout = layout,**kwargs)
        self.set_stabilizer()

        if make_stab == True:
            self.set_stabilizer_color()

        if with_back == True:
            self.highlight_edges()
        else:
            for edge in self.edges:
                self.edges[edge].set_color(GREEN)

        self.shift(self.offset)



    def get_breadth_Grids(self):
        if(self.back == False):
            layers = list(nx.bfs_layers(self._graph, nx.center(self._graph)))[::-1]
            output = []
            for i in range(1, len(layers)):
                prev_layer = layers[i - 1]
                curr_layer = layers[i]
                H = nx.Graph()
                H.add_nodes_from(prev_layer)
                H.add_nodes_from(curr_layer)
                output_edge_i = []
                for u in curr_layer:
                    for v in self._graph[u]:
                        if v in prev_layer:
                            H.add_edge(v, u)
                for edge in H.edges():
                    if(edge[0][0] >= self.tile.size):
                        output_edge_i.append((edge[1],edge[0]))
                    else:
                        output_edge_i.append(edge)
                output.append(output_edge_i)

            #Deal with last layer case (may need to connect edge between two vertices of last layer)
            H = nx.Graph()
            output_edge_i = []
            for u in layers[len(layers)-1]:
                for v in self._graph[u]:
                    if v in layers[len(layers)-1]:
                        H.add_edge(v, u)
            for u in layers[len(layers)-1]:
                for v in self._graph[u]:
                    if v in layers[len(layers)-1]:
                        H.add_edge(v, u)
            for edge in H.edges():
                if(edge[0][0] >= self.tile.size):
                    output_edge_i.append((edge[1],edge[0]))
                else:
                    output_edge_i.append(edge)
            output.append(output_edge_i)

            output_graph = []

            for breadth in output:
                output_graph.append(Grid(self.tile,make_stab=False,highlightedges = breadth, with_back = False))

            self.tile.breadthGraphs = output_graph

    #Returns new graph with spanning tree lines highlighted (changes the color of all edges in highedges)
    def highlight_edges(self):
        for edge in self.highedge:
            self.edges[edge].set_color(GREEN)
            self.edges[edge].set_z_index(3)

    # computes the approximate steiner tree
    def compute_MST(self):
        steiner_graph = nx.algorithms.approximation.steiner_tree(
            self._graph, # pass the manim graph as a networkx graph (this is the parent graph)
            self.stabilizers, # pass the list of qubits in the support of a stabilizer (list of terminal nodes)
            weight = 'weight', # defaults to weight 1 for all edges
            method = "mehlhorn" # the underlying algorithm
        )
        return_list = []

        #Note that edges are always append with the *lowest* (i.e. hor) qubit first, then higher (i.e. ver) qubit
        #Thus if first coordinate has coord > self.tile.size, the it must be a ver coord and the edge vertices should be swapped in the return list
        for edge in steiner_graph.edges():
            if(edge[0][0] >= self.tile.size):
                return_list.append((edge[1],edge[0]))
            else:
                return_list.append(edge)


        # put the steiner tree into the tile grid
        self.tile.MST = Grid(self.tile,highlightedges = return_list, with_back = False)
        self.tile.MST.get_breadth_Grids()
        return self.tile.MST


    def set_stabilizer(self):
        for tuple in self.tile.stab_ver:
            self.stabilizers.append(tuple)
        for tuple in self.tile.stab_hor:
            new_tuple = (tuple[0]+self.tile.size,tuple[1]+self.tile.size)
            self.stabilizers.append(new_tuple)

    def set_stabilizer_color(self):
        for stab in self.stabilizers:
            self[stab].set_fill(ORANGE, opacity=1)


    """
    Takes a steiner tree `F_a` for a stabilizer of type `a` in {X,Z}
    and a `root` qubit in F_a that will be the ancilla (the qubit that all of the information converges to)
    and computes the layered scheduling for that tree.
    This is almost a direct implementation of Algorithm 3: Route-Pack-PerType
    Returns: `schedule`: a list of nx.Graph 
    NOTE: THIS ONLY WORKS FOR A SINGLE STABILIZER/TREE right now, since it uses only a single root
            We may need a separate algorithm to combine these for the forest of trees. 
            Or perhaps the forest is created, the conflicts are found, and then something is done there?
            May need to take a look at Emerson's code
    """
    def compute_layered_schedule_from_conflict_graph(self, F_a: nx.Graph, root) -> list[nx.Graph]:
        schedule: list[nx.Graph] = []
        # check if root is in F_a
        if root not in F_a:
            raise ValueError(f"Root '{root}' is not present in graph F_a")

        # use BFS (breadth-first search) from `root` to determine the layers (time steps)
        # - Note that layer 0 is the farthest from the root, so we reverse the list from bfs_layers
        layers: list[nx.Graph] = list(nx.bfs_layers(F_a, root)).reverse()
        num_layers = layers.len()
        """
        to do this globally, just need to compute the layer for each tree individually. then take the union of all 
        of the "layers" of a specific level and then apply the below script as-is
        this would change what the inputs to the function are (a list of tuples of trees and roots [(F_a, root)])
        """

        # for each "time step" from the outside in, deal with the conflicts by further partitioning
        for l in num_layers:
            # build the conflict graph C_l

            C_l: nx.Graph = layers[l]
            """ brainstorming about what the conflict graph is and how to make it
            - perhaps put creating the conflict graph in a separate method?
            - I think the conflict graph has as its nodes the edges in the layer
            -   that are not orphans, though the orphans still need to be scheduled, so should all be put in?
            - the conflict graph is structured somewhat differently, or maybe I am overthinking this, since this
            - method takes in the steiner tree, which is already in the dual graph
            - upon further thinking I think I am overthinking this (wow that's a sentence, but I'm leaving it in for now, since I may be wrong)
            - thus I think the conflict graph is just a layer of the steiner tree.
            """

            # TODO later for completeness: add interface reservations for boundaries/overlaps

            # compute delta: the maximum degree of C_l
            degrees = C_l.degree()
            delta = max(degree for node, degree in degrees)

            # determine if C_l is bipartite (allows for computational speedup)
            is_bipartite = nx.is_bipartite(C_l)

            if is_bipartite:
                for c in range(delta):
                    # for c=1 to delta (a matching removes one edge from a node, decreasing its degree by 1,
                    #                   so we need to iterate delta times to account for all edges)
                    M_c = nx.algorithms.bipartite.hopcroft_karp_matching(C_l) # returns a dictionary of corresponding (matched) nodes
                    # note that this is redundant (A:B and B:A are both included in the dictionary, but is likely good for convenience)
                    schedule.append(nx.Graph(M_c))
                    C_l.remove_edges_from(M_c.items()) # remove the matched edges
            else:
                coloring_dict = edge_coloring(C_l) # returns a dictionary of (A,B):color_as_int

                # Swap the keys/values so we now have a dictionary where the key is a color (integer)
                # and the value is a list of edges
                colors_to_edges = defaultdict(list)
                for edge, color in coloring_dict.items():
                    colors_to_edges[color].append(edge)

                # create a new nx.Graph for each color and append each new graph to `schedule`
                for color, edge_list in colors_to_edges.items():
                    G_color = nx.Graph()
                    G_color.add_edges_from(edge_list)
                    schedule.append(G_color)

        return schedule





#Takes in list of tuples denoting horizontal and vertical stabilizers, and type
#True -> X, False -> Z
#level increases layer in which drawn
# A single tile represents a single stabilizer in the tile code,
#   with some number of "vertical" & "Horizontal" qubits (here denoted stabilizers)
#   in its support

class Tile(VMobject):

    #Offset offsets location of tile on global grid
    def __init__(self,size=1, stab_hor = [], stab_ver = [], type = True, runStab = True, level = 0, hor = {}, ver = {}, pointss = {}, gridArrayHor = {}, gridArrayVer = {}, global_setup = False, offset = np.array([0,0,0]), globlesize = 100, symetric_check = False,**kwargs):
        super().__init__(**kwargs)

        self.size = size # linear size of the tile
        self.hor = hor
        self.ver = ver
        self.pointss = pointss
        self.stab_hor = stab_hor # support qubits on "horizontal edges" of the primary lattice
        self.stab_ver = stab_ver # support qubits on "vertical edges" of the primary lattice
        self.type = type # X or Z stabilizer
        self.level = level # which layer is it drawn in
        self.gridArrayHor = gridArrayHor
        self.gridArrayVer = gridArrayVer
        self.grid = None
        self.global_setup = global_setup
        self.offset = offset
        self.globlesize = globlesize

        # Only set once compute_MST() is run on grid object derived from this Tile
        self.MST = None
        #list of grid objects which contain breadth layers graphs for this *single* Tile object
        self.breadthGraphs = None

        if(global_setup == False or symetric_check == True):
            self.createTile()
        else:
            self.setUpTile()

        if(runStab == True):
            self.setStabilizer()

        #self.grid.move_to(self)
    def setPoints(self,dict):
        copied = dict.get_all_points().copy()
        for point in copied:
            point[0] += self.offset[0]
            point[1] += self.offset[1]
        dict.set_points(copied)

    def setUpTile(self):
        for i in range(self.size):
            for j in range(self.size):
                if ((i+self.offset[0] < self.globlesize and j+self.offset[1] < self.globlesize) and (i+self.offset[0] >= 0 and j+self.offset[1] >=0)) and self.global_setup == True:
                    if self.hor.__contains__((i,j)) and self.ver.__contains__((i,j)):
                        self.hor[(i,j)] = self.hor[(i,j)].copy().set_color(WHITE)
                        self.ver[(i,j)] = self.ver[(i,j)].copy().set_color(WHITE)
                        self.add(self.hor[(i,j)])
                        self.add(self.ver[(i,j)])
                        self.setPoints(self.hor[(i,j)])
                        self.setPoints(self.ver[(i,j)])


    # creates a new, empty tile (empty as in containing no support qubits for a stabilizer)
    def createTile(self):
        for i in range(self.size):
            for j in range(self.size):
                if ((i+self.offset[0] < self.globlesize and j+self.offset[1] < self.globlesize) and (i+self.offset[0] >= 0 and j+self.offset[1] >=0))  and self.global_setup == True:
                    curDot = Dot([i,j,0]).set_color(WHITE).scale(1)
                    hor = Line([i,j,0],[i+1,j,0]).set_z_index(1)
                    ver = Line([i,j,0],[i,j+1,0]).set_z_index(1)
                    self.gridArrayHor[(i,j)] = hor.get_center()
                    self.gridArrayVer[(i,j)] = ver.get_center()

                    self.hor[(i,j)] = hor
                    self.ver[(i,j)] = ver
                    self.pointss[(i,j)] = curDot
                    self.add(hor)
                    self.add(ver)
                    self.setPoints(self.hor[(i,j)])
                    self.setPoints(self.ver[(i,j)])



    # set the color of the support qubits of a stabilizer (Tile) and
    # set the level (layer) at which it is drawn
    def setStabilizer(self):
        for coord in self.stab_hor:
            if ((coord[0]+self.offset[0] < self.globlesize and coord[1]+self.offset[1] < self.globlesize) and (coord[0]+self.offset[0] >= 0 and coord[1]+self.offset[1] >=0))  and self.global_setup == True:
                if self.hor.__contains__(coord):
                    if self.type == True: # X type
                        self.hor[coord].set_color(RED).set_z_index(1+self.level)
                    else: # Z type
                        self.hor[coord].set_color(BLUE).set_z_index(1+self.level)
        for coord in self.stab_ver:
            if ((coord[0]+self.offset[0] < self.globlesize and coord[1]+self.offset[1] < self.globlesize) and (coord[0]+self.offset[0] >= 0 and coord[1]+self.offset[1] >=0)) and self.global_setup == True:
                if self.ver.__contains__(coord):
                    if self.type == True: # X type
                        self.ver[coord].set_color(RED).set_z_index(1+self.level)
                    else: # Z type
                        self.ver[coord].set_color(BLUE).set_z_index(1+self.level)

    # set the color of support qubits for two tiles, setting overlapping (potentially conflicting)
    # support qubits to purple
    def setStabilizerWithOverlap(self, tile1, tile2):
        for coord1 in tile1.stab_hor:
            for coord2 in tile2.stab_hor:
                updatedcoord1 = (coord1[0]+tile1.offset[0],coord1[1]+tile1.offset[1])
                updatedcoord2 = (coord2[0]+tile2.offset[0],coord2[1]+tile2.offset[1])
                if(updatedcoord1 == updatedcoord2):
                    if ((coord1[0]+tile1.offset[0] < self.globlesize and coord1[1]+tile1.offset[1] < self.globlesize) and (coord1[0]+tile1.offset[0] >= 0 and coord1[1]+tile1.offset[1] >=0)) and self.global_setup == True:
                        self.hor[coord1].set_color(PURPLE).set_z_index(1+self.level)
        for coord1 in tile1.stab_ver:
            for coord2 in tile2.stab_ver:
                updatedcoord1 = (coord1[0]+tile1.offset[0],coord1[1]+tile1.offset[1])
                updatedcoord2 = (coord2[0]+tile2.offset[0],coord2[1]+tile2.offset[1])
                if(updatedcoord1 == updatedcoord2):
                    if ((coord1[0]+tile1.offset[0] < self.globlesize and coord1[1]+tile1.offset[1] < self.globlesize) and (coord1[0]+tile1.offset[0] >= 0 and coord1[1]+tile1.offset[1] >=0)) and self.global_setup == True:
                        self.ver[coord1].set_color(PURPLE).set_z_index(1+self.level)

    #Returns Corresponding (dual) X/Z tile to current tile
    def makeSymetricTile(self, globle_setup):
        new_stab_hor = []
        new_stab_ver = []
        for point in self.stab_hor:
            new_stab_ver.append(self.findSymPoint(point))
        for point in self.stab_ver:
            new_stab_hor.append(self.findSymPoint(point))
        return Tile(self.size,new_stab_hor,new_stab_ver, not self.type, offset = self.offset, global_setup = globle_setup, globlesize=self.globlesize, symetric_check = True)

    # Create a new tile whose support qubits are the union of the two input tiles (self and otherTile)
    def makeOverlapTile(self,otherTile,level,globle_setup):
        new_stab_hor = self.stab_hor + otherTile.stab_hor
        new_stab_ver = self.stab_ver + otherTile.stab_ver
        overlapTile = Tile(self.size,new_stab_hor,new_stab_ver, not self.type, False,level, offset = self.offset, global_setup = globle_setup, globlesize=self.globlesize )
        overlapTile.setStabilizerWithOverlap(self,otherTile)
        return overlapTile

    # map the coordinates of a qubit in a X(Z) Tile to its dual in a Z(X) tile
    def findSymPoint(self,point):
        return (self.size-1-point[0],self.size-1-point[1])

    def getPoint(self, point):
        return self.pointss[point]
    def getHor(self,point):
        return self.hor[point]
    def getVer(self,point):
        return self.ver[point]
    def getGrid(self):
        return self.grid
    def setGrid(self,grid):
        self.grid = grid

class Global_Tile(VMobject):
    def __init__(self,gridsize=1,tilesize=1, **kwargs):
        super().__init__(**kwargs)

        #tiles has tuple int keys (with key coord cooresponding to bottom left point) and tile object values
        self.tiles = []

        self.size = gridsize
        self.tilesize = tilesize
        self.hor = {}
        self.ver = {}
        self.pointss = {}
        self.gridArrayHor = {}
        self.gridArrayVer = {}
        self.createGrid()
        self.MEGATiles_TRUE = {}
        self.MEGATiles_FALSE = {}
        self.MEGABipartite_TRUE = {}
        self.MEGABipartite_FALSE = {}

    def createGrid(self):
        for i in range(self.size):
            for j in range(self.size):
                curDot = Dot([i,j,0]).set_color(WHITE).scale(1)
                hor = Line([i,j,0],[i+1,j,0], color=GREY).set_z_index(0)
                ver = Line([i,j,0],[i,j+1,0], color=GREY).set_z_index(0)
                self.gridArrayHor[(i,j)] = hor.get_center()
                self.gridArrayVer[(i,j)] = ver.get_center()
                self.hor[(i,j)] = hor
                self.ver[(i,j)] = ver
                self.pointss[(i,j)] = curDot
                self.add(hor)
                self.add(ver)

    #Takes in coord = (i,j) which specifies bottom leftmost point of tile on global grid
    #0<i,j<(gridsize-tilesize)
    #stab_hor and stab_ver are lists of tuples, corresponding to horizontal and vertical sides respectively
    #Type corresponds to type of stabilizer, True -> X, False -> Z
    def createTile(self, coord,stab_hor,stab_ver,type = True, MEGA = False):
        current_hor = {}
        current_ver = {}
        current_pointss = {}
        current_gridArrayHor = {}
        current_gridArrayVer = {}
        offset = np.array([coord[0],coord[1],0])
        if MEGA == False:
            for i in range(self.tilesize):
                for j in range(self.tilesize):
                    if (i+coord[0] < self.size and j+coord[1] < self.size) and (i+coord[0] >= 0 and j+coord[1] >=0):
                        current_hor[(i,j)] = self.hor[(i,j)]
                        current_ver[(i,j)] = self.ver[(i,j)]
                        current_pointss[(i,j)] = self.pointss[(i,j)]
                        current_gridArrayHor[(i,j)] = self.gridArrayHor[(i,j)]
                        current_gridArrayVer[(i,j)] = self.gridArrayVer[(i,j)]
        else:
            for i in range(self.size):
                for j in range(self.size):
                    if (i+coord[0] < self.size and j+coord[1] < self.size) and (i+coord[0] >= 0 and j+coord[1] >=0):
                        current_hor[(i,j)] = self.hor[(i,j)]
                        current_ver[(i,j)] = self.ver[(i,j)]
                        current_pointss[(i,j)] = self.pointss[(i,j)]
                        current_gridArrayHor[(i,j)] = self.gridArrayHor[(i,j)]
                        current_gridArrayVer[(i,j)] = self.gridArrayVer[(i,j)]
        if(MEGA == False):
            self.tiles.append(Tile(self.tilesize,stab_hor,stab_ver,type = type, level = 1, hor = current_hor, ver = current_ver, pointss = current_pointss, gridArrayHor=current_gridArrayHor,gridArrayVer=current_gridArrayVer,global_setup = True, offset = offset, globlesize=self.size))
            return self.tiles[len(self.tiles)-1]
        else:
            return Tile(self.size,stab_hor,stab_ver,type = type, level = 1, hor = current_hor, ver = current_ver, pointss = current_pointss, gridArrayHor=current_gridArrayHor,gridArrayVer=current_gridArrayVer,global_setup = True, offset = offset, globlesize=self.size)

            #MegaGraph is what I've termed each of the two bipartite graphs for given X/Z from given depth level of steiner trees of all tiles present
    #Will include nodes across whole of global grid
    def createMegaGraph(self,depth,type):
        H = nx.Graph()
        for tile in self.tiles:
            if tile.type == type:
                if(tile.breadthGraphs != None and len(tile.breadthGraphs) > depth):
                    newNodes = []
                    for node in tile.breadthGraphs[depth]._graph.nodes():

                        if((int(node[0]+tile.offset[0]),int(node[1]+tile.offset[1])) not in newNodes):
                            newNodes.append((int(node[0]+tile.offset[0]),int(node[1]+tile.offset[1])))
                    newEdges = []
                    for edge in tile.breadthGraphs[depth]._graph.edges():

                        newVertex1 = (int(edge[0][0]+tile.offset[0]),int(edge[0][1]+tile.offset[1]))
                        newVertex2 = (int(edge[1][0]+tile.offset[0] + self.size - tile.size),int(edge[1][1]+tile.offset[1] + self.size - tile.size))
                        newEdge = (newVertex1,newVertex2)
                        if newEdge not in newEdges:
                            newEdges.append(newEdge)

                    H.add_nodes_from(newNodes)
                    H.add_edges_from(newEdges)
        BIGTILE = self.createTile((0,0),[],[],type = type,MEGA = True)
        grid = Grid(BIGTILE,make_stab=False,highlightedges = list(H.edges()), with_back = False)

        self.createBiPartiteMega(depth,grid,type)

        if(type == True):
            self.MEGATiles_TRUE[depth] = grid
            return self.MEGATiles_TRUE[depth]
        else:
            self.MEGATiles_FALSE[depth] = grid
            return self.MEGATiles_FALSE[depth]

    def createBiPartiteMega(self,depth,grid,type):
        nx_MEGA = grid._graph.copy()

        # determine if C_l is bipartite (allows for computational speedup)
        is_bipartite = nx.is_bipartite(nx_MEGA)

        output_edge1 = []
        if is_bipartite:
            M_c = nx.Graph(self.hopcroft_karp_on_disjoint(nx_MEGA)) # returns a dictionary of corresponding (matched) nodes

            output_edge1 = list(M_c.edges())
            nx_MEGA.remove_edges_from(list(M_c.edges())) # remove the matched edges

            output_edge2 = list(nx_MEGA.edges())

            BIGTILE1 = self.createTile((0,0),[],[],type = type,MEGA = True)
            BIGTILE2 = self.createTile((0,0),[],[],type = type,MEGA = True)

            color1 = Grid(BIGTILE1,make_stab=False,highlightedges = output_edge1, with_back = False)
            color2 = Grid(BIGTILE2,make_stab=False,highlightedges = output_edge2, with_back = False)
            if(type == True):
                self.MEGABipartite_TRUE[depth] = (color1,color2)
            else:
                self.MEGABipartite_FALSE[depth] = (color1,color2)

        else:
            print("ERROR")

    def hopcroft_karp_on_disjoint(self, G):
        # final matching edge list
        final_edges = []

        # loop over connected components
        for nodes in nx.connected_components(G):
            # extract subgraph for this component
            H = G.subgraph(nodes).copy()

            # run Hopcroft-Karp on this connected component
            matching = nx.algorithms.bipartite.hopcroft_karp_matching(H)

            # convert matching dict → unique edge list
            edges = [(u, v) for u, v in matching.items() if u < v]
            for edge in edges:
                if(edge[0][0] >= self.size):
                    final_edges.append((edge[1],edge[0]))
                else:
                    final_edges.append(edge)
            final_edges.extend(edges)

        return final_edges






# scene that is currently being used for testing
# the main scene (manim renders scenes)
class Tile_Scene(MovingCameraScene):
    def showTile(self,vertex,a):
        self.play(FadeIn(vertex))
        self.play(FadeIn(a))
        self.play(FadeOut(vertex))
        self.play(FadeOut(a))

    def construct(self):
        self.wait(1)
        self.next_section(skip_animations=True)

        self.camera.frame_width = 28
        self.camera.frame_height = 7.875*2

        global_tile = Global_Tile(10,5);
        self.camera.frame.move_to(global_tile)
        self.add(global_tile)

        stab_hor = [(0,1),(1,1),(2,2),(4,4),(3,4)] # horizontal support qubits
        stab_ver = [(1,0),(1,1),(2,4),(3,3)] # vertical support qubits

        tiles_X = []
        MST_X = []
        tiles_Z = []
        MST_Z = []
        #Make vertices
        for i in range(5):
            for j in range(5):
                if(i%2 == 0 and j%2 == 1):
                    newTile = global_tile.createTile((i,j),stab_hor,stab_ver)
                    newTile.setGrid(Grid(newTile))
                    newMST = newTile.getGrid().compute_MST()
                    tiles_X.append(newTile)
                    MST_X.append(newMST)
                    self.showTile(newTile,newMST)
                if(i%2 == 1 and j%2 == 0):
                    newTile = global_tile.createTile((i,j),stab_hor,stab_ver)
                    newTileZ = newTile.makeSymetricTile(True)
                    newTileZ.setGrid(Grid(newTileZ))
                    newMSTZ = newTileZ.getGrid().compute_MST()
                    tiles_Z.append(newTileZ)
                    MST_Z.append(newMSTZ)
                    self.showTile(newTileZ,newMSTZ)



        self.next_section()



        #plaquette = vertex.makeSymetricTile(True) # dual "plaquette" operator
        #plaquette2 = vertex2.makeSymetricTile(True) # dual "plaquette" operator
        #self.add(plaquette)
        #self.add(vertex.makeOverlapTile(plaquette2,2,True))




        for i in range(10):
            print(i)
            global_tile.createMegaGraph(i,True)
            global_tile.createMegaGraph(i,False)
        #global_tile.createMegaGraph(1,True)
        #self.play(FadeIn(global_tile.MEGATiles[1]))


        sub_colors_X = global_tile.MEGABipartite_TRUE
        sub_colors_Z = global_tile.MEGABipartite_FALSE


        for i in range(len(sub_colors_X)):
            self.play(FadeIn(sub_colors_X[i][0]))
            self.play(FadeIn(sub_colors_X[i][1]))
            self.play(FadeOut(sub_colors_X[i][0]))
            self.play(FadeOut(sub_colors_X[i][1]))

            self.play(FadeIn(sub_colors_Z[i][0]))
            self.play(FadeIn(sub_colors_Z[i][1]))
            self.play(FadeOut(sub_colors_Z[i][0]))
            self.play(FadeOut(sub_colors_Z[i][1]))

        for i in range(len(global_tile.MEGATiles_TRUE)):
            self.play(FadeIn(global_tile.MEGATiles_TRUE[i]))
            self.play(FadeIn(global_tile.MEGATiles_FALSE[i]))

        #self.play(steiner.animate.shift(UP).rotate(PI / 3),run_time=4)

        #self.camera.frame.move_to(steiner)
        #self.add(steiner)
        #plaquette.next_to(vertex, LEFT) # set the two tiles as adjacent to one another

        # add the two tiles to the manim scene
        #self.add(grid)
        #self.add(steiner)
        #self.add(vertex2)
        #self.add(vertex3)
        #self.add(plaquette)
        #self.add(plaquette2)


        # set the initial camera location
        #self.camera.frame.move_to(vertex)
        #self.add(vertex.getGrid())
        #self.wait(0.5)
        #self.combineTwoTileAnimation(vertex,plaquette,3) # manim: move "plaquette" on top of "vertex"

        #self.play(FadeIn(vertex.getGrid()))
        #self.wait(2)


        #self.play(Transform(plaquette,vertex), run_time=5)
        #self.wait()


    def combineTwoTileAnimation(self,tile1,tile2,time):
        overlap = tile1.makeOverlapTile(tile2,1, False)
        overlap.move_to(tile2) # sets the overlap (with purple) to the location of `tile2`. Does not render yet
        group = VGroup(tile1, tile2) # create a group of manim objects (MObjects)
        self.camera.frame.move_to(group)
        self.play(tile1.animate.move_to(tile2),run_time=time) # main sliding animation loop
        #tile1.getGrid().move_to(tile2)
        self.play(FadeIn(overlap),run_time=time) # fade in the purple overlaps


""" Implementation notes:
the networkx steiner approximation is listed as the last algorithm (Mehlhorn) in Emerson's paper

## Graph matching
networkx has hopcroft Karp, which is for graph matching. It requires a bipartite graph
- `networkx.algorithms.bipartite.hopcroft_karp_matching(G)`

For nonbipartite graphs, networkx has the following O(V^3) algorithm 
(although theoretical O(E*sqrt(V)) algorithms exist (eg Micali-Vazirani), they are notoriously hard to implement, 
are not in networkx, and have overhead that would likely not be worth it in our case. If we really
need speed there are python wrappers for the C++ BlossomV, which is still O(V^3) but with many heuristics
that make it quite fast in practice.
- `networkx.algorithms.matching.max_weight_matching(G, maxcardinality=True)`

Note that it is possible that either of these will throw an error (I think)
if there are an odd number of nodes (because an edge connects two nodes)


## Graph coloring
networkx does not implement Misra–Gries or any other graph *edge* coloring algorithm (though it has *vertex* coloring algorithms)
- The package "GCol" is compatible with networkx and implements `gcol.edge_coloring(G)`


"""
"""
Wait, but are there not three types of conflicts? I may have been thinking of the conflict graph incorrectly. 
- for 2 directions of CNOTs and 1 SWAP?
- Actually we can probably limit it to just CNOTs and swaps, or rather this is only an issue for combining multiple tiles I think?
"""
