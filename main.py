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
    def __init__(self,tile,make_stab=True,highlightedges = [],**kwargs):
        self.tile = tile
        vertices = []
        edges = []
        layout = {}
        self.stabilizers = []
        self.highedge = highlightedges
        #use first tile.size^2 entries for vertical, next tile.size^2 for horizontal
        for i in range(tile.size):
            for j in range(tile.size):
                vertices.append((i,j))
                vertices.append((i+tile.size**2,j+tile.size**2))
                edges.append(((i,j),(i+tile.size**2,j+tile.size**2)))
                if(j+1<tile.size):
                    edges.append(((i,j),(i+tile.size**2,j+tile.size**2+1)))
                if(i+1<tile.size):
                    edges.append(((i+1,j),(i+tile.size**2,j+tile.size**2)))
                if((j+1<tile.size) and (i+1<tile.size)):
                    edges.append(((i+1,j),(i+tile.size**2,j+tile.size**2+1)))

        curpoint = (0,0)
        for point in tile.gridArrayVer:
            ref = tile.gridArrayVer[point]
            layout[curpoint] = [ref[0],ref[1],0]

            horpoint = (curpoint[0]+tile.size**2,curpoint[1]+tile.size**2)
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
        self.highlight_edges()

    #Returns new graph with spanning tree lines highlighted (changes the color of all edges in highedges)
    def highlight_edges(self):
        for edge in self.highedge:
            self.edges[edge].set_color(GREEN)

    # computes the approximate steiner tree 
    def compute_MST(self):
        steiner_graph = nx.algorithms.approximation.steiner_tree(
            self._graph, # pass the manim graph as a networkx graph (this is the parent graph)
            self.stabilizers, # pass the list of qubits in the support of a stabilizer (list of terminal nodes)
            weight = 'weight', # defaults to weight 1 for all edges
            method = "mehlhorn" # the underlying algorithm
        )
        return_list = []
        for edge in steiner_graph.edges():
            if(edge[0][0] >= self.tile.size**2):
                return_list.append((edge[1],edge[0]))
            else:
                return_list.append(edge)
        # put the steiner tree into the tile grid
        return Grid(self.tile,highlightedges = return_list)


    def set_stabilizer(self):
        for tuple in self.tile.stab_ver:
            self.stabilizers.append(tuple)
        for tuple in self.tile.stab_hor:
            new_tuple = (tuple[0]+self.tile.size**2,tuple[1]+self.tile.size**2)
            self.stabilizers.append(new_tuple)

    def set_stabilizer_color(self):
        for stab in self.stabilizers:
            self[stab].set_fill(RED, opacity=1)

        
    """
    Takes a steiner tree `F_a` for a stabilizer of type `a` in {X,Z}
    and a `root` qubit in F_a that will be the ancilla (the qubit that all of the information converges to)
    and computes the layered scheduling for that tree.
    This is almost a direct implementation of Algorithm 3: Route-Pack-PerType
    Returns: `schedule`: a list of nx.Graph 
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

        # for each "time step" from the outside in, deal with the conflicts by further partitioning
        for l in num_layers:
            # build the conflict graph C_l
            C_l: nx.Graph = layers[l]

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
    def __init__(self,size=1, stab_hor = [], stab_ver = [], type = True, runStab = True, level = 0,**kwargs):
        super().__init__(**kwargs)

        self.size = size # linear size of the tile
        self.hor = {}
        self.ver = {}
        self.pointss = {}
        self.stab_hor = stab_hor # support qubits on "horizontal edges" of the primary lattice
        self.stab_ver = stab_ver # support qubits on "vertical edges" of the primary lattice
        self.type = type # X or Z stabilizer
        self.level = level # which layer is it drawn in
        self.gridArrayHor = {}
        self.gridArrayVer = {}
        self.grid = None

        self.createTile()
        if(runStab == True):
            self.setStabilizer()

        #self.grid.move_to(self)

    # creates a new, empty tile (empty as in containing no support qubits for a stabilizer)
    def createTile(self):
        for i in range(self.size):
            for j in range(self.size):
                curDot = Dot([i,j,0]).set_color(WHITE).scale(1)
                hor = Line([i,j,0],[i+1,j,0]).set_z_index(0)
                ver = Line([i,j,0],[i,j+1,0]).set_z_index(0)
                self.gridArrayHor[(i,j)] = hor.get_center()
                self.gridArrayVer[(i,j)] = ver.get_center()

                self.hor[(i,j)] = hor
                self.ver[(i,j)] = ver
                self.pointss[(i,j)] = curDot
                self.add(hor)
                self.add(ver)

    # set the color of the support qubits of a stabilizer (Tile) and 
    # set the level (layer) at which it is drawn
    def setStabilizer(self):
        for coord in self.stab_hor:
            if self.type == True: # X type
                self.hor[coord].set_color(RED).set_z_index(1+self.level)
            else: # Z type
                self.hor[coord].set_color(BLUE).set_z_index(1+self.level)
        for coord in self.stab_ver:
            if self.type == True: # X type
                self.ver[coord].set_color(RED).set_z_index(1+self.level)
            else: # Z type
                self.ver[coord].set_color(BLUE).set_z_index(1+self.level)

    # set the color of support qubits for two tiles, setting overlapping (potentially conflicting)
    # support qubits to purple
    def setStabilizerWithOverlap(self, tile1, tile2):
        for coord1 in tile1.stab_hor:
            for coord2 in tile2.stab_hor:
                if(coord1 == coord2):
                    self.hor[coord1].set_color(PURPLE).set_z_index(1+self.level)
        for coord1 in tile1.stab_ver:
            for coord2 in tile2.stab_ver:
                if(coord1 == coord2):
                    self.ver[coord1].set_color(PURPLE).set_z_index(1+self.level)

    #Returns Corresponding (dual) X/Z tile to current tile
    def makeSymetricTile(self):
        new_stab_hor = []
        new_stab_ver = []
        for point in self.stab_hor:
            new_stab_ver.append(self.findSymPoint(point))
        for point in self.stab_ver:
            new_stab_hor.append(self.findSymPoint(point))
        return Tile(self.size,new_stab_hor,new_stab_ver, not self.type)

    # Create a new tile whose support qubits are the union of the two input tiles (self and otherTile)
    def makeOverlapTile(self,otherTile,level):
        new_stab_hor = self.stab_hor + otherTile.stab_hor
        new_stab_ver = self.stab_ver + otherTile.stab_ver
        overlapTile = Tile(self.size,new_stab_hor,new_stab_ver, not self.type, False,level)
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





# scene that is currently being used for testing
# the main scene (manim renders scenes)
class Tile_Scene(MovingCameraScene):
    def construct(self):
        stab_hor = [(0,1),(1,1),(2,2),(4,4),(3,4)] # horizontal support qubits
        stab_ver = [(1,0),(1,1),(2,4),(3,3)] # vertical support qubits
        vertex = Tile(5,stab_hor,stab_ver,True) # "star" operator
        plaquette = vertex.makeSymetricTile() # dual "plaquette" operator

        #vertex.setGrid(Grid(vertex))
        #grid = vertex.getGrid()
        #steiner = grid.compute_MST()

        #self.play(steiner.animate.shift(UP).rotate(PI / 3),run_time=4)

        #self.camera.frame.move_to(steiner)
        plaquette.next_to(vertex, LEFT) # set the two tiles as adjacent to one another

        # add the two tiles to the manim scene
        self.add(vertex)
        self.add(plaquette)


        # set the initial camera location
        self.camera.frame.move_to(vertex)
        #self.add(vertex.getGrid())
        self.wait(0.5)
        self.combineTwoTileAnimation(vertex,plaquette,3) # manim: move "plaquette" on top of "vertex"

        #self.play(FadeIn(vertex.getGrid()))
        self.wait(2)


        #self.play(Transform(plaquette,vertex), run_time=5)
        #self.wait()
    def combineTwoTileAnimation(self,tile1,tile2,time):
        overlap = tile1.makeOverlapTile(tile2,1)
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