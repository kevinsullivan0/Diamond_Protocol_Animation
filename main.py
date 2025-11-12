from manim import *
import networkx as nx

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

    #Returns new graph with spanning tree lines highlighted
    def highlight_edges(self):
        for edge in self.highedge:
            self.edges[edge].set_color(GREEN)

    def compute_MST(self):
        steiner_graph = nx.algorithms.approximation.steiner_tree(
            self._graph,
            self.stabilizers,
            weight = 'weight',
            method = "mehlhorn"
        )
        return_list = []
        for edge in steiner_graph.edges():
            if(edge[0][0] >= self.tile.size**2):
                return_list.append((edge[1],edge[0]))
            else:
                return_list.append(edge)
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





#Takes in list of tuples denoting horizontal and vertical stabilizers, and type
#True -> X, False -> Z
#level increases layer in which drawn
class Tile(VMobject):
    def __init__(self,size=1, stab_hor = [], stab_ver = [], type = True, runStab = True, level = 0,**kwargs):
        super().__init__(**kwargs)

        self.size = size
        self.hor = {}
        self.ver = {}
        self.pointss = {}
        self.stab_hor = stab_hor
        self.stab_ver = stab_ver
        self.type = type
        self.level = level
        self.gridArrayHor = {}
        self.gridArrayVer = {}
        self.grid = None

        self.createTile()
        if(runStab == True):
            self.setStabilizer()

        #self.grid.move_to(self)


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


    def setStabilizer(self):
        for coord in self.stab_hor:
            if self.type == True:
                self.hor[coord].set_color(RED).set_z_index(1+self.level)
            else:
                self.hor[coord].set_color(BLUE).set_z_index(1+self.level)
        for coord in self.stab_ver:
            if self.type == True:
                self.ver[coord].set_color(RED).set_z_index(1+self.level)
            else:
                self.ver[coord].set_color(BLUE).set_z_index(1+self.level)
    def setStabilizerWithOverlap(self, tile1,tile2):
        for coord1 in tile1.stab_hor:
            for coord2 in tile2.stab_hor:
                if(coord1 == coord2):
                    self.hor[coord1].set_color(PURPLE).set_z_index(1+self.level)
        for coord1 in tile1.stab_ver:
            for coord2 in tile2.stab_ver:
                if(coord1 == coord2):
                    self.ver[coord1].set_color(PURPLE).set_z_index(1+self.level)

    #Returns Corresponding X/Z tile to current tile
    def makeSymetricTile(self):
        new_stab_hor = []
        new_stab_ver = []
        for point in self.stab_hor:
            new_stab_ver.append(self.findSymPoint(point))
        for point in self.stab_ver:
            new_stab_hor.append(self.findSymPoint(point))
        return Tile(self.size,new_stab_hor,new_stab_ver, not self.type)

    def makeOverlapTile(self,otherTile,level):
        new_stab_hor = self.stab_hor + otherTile.stab_hor
        new_stab_ver = self.stab_ver + otherTile.stab_ver
        overlapTile = Tile(self.size,new_stab_hor,new_stab_ver, not self.type, False,level)
        overlapTile.setStabilizerWithOverlap(self,otherTile)
        return overlapTile

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






class Tile_Scene(MovingCameraScene):
    def construct(self):
        stab_hor = [(0,1),(1,1),(2,2),(4,4),(3,4)]
        stab_ver = [(1,0),(1,1),(2,4),(3,3)]
        vertex = Tile(5,stab_hor,stab_ver,True)
        plaquette = vertex.makeSymetricTile()

        #vertex.setGrid(Grid(vertex))
        #grid = vertex.getGrid()
        #steiner = grid.compute_MST()

        #self.play(steiner.animate.shift(UP).rotate(PI / 3),run_time=4)

        #self.camera.frame.move_to(steiner)
        plaquette.next_to(vertex, LEFT)

        self.add(vertex)
        self.add(plaquette)



        self.camera.frame.move_to(vertex)
        #self.add(vertex.getGrid())
        self.wait(0.5)
        self.combineTwoTileAnimation(vertex,plaquette,3)

        #self.play(FadeIn(vertex.getGrid()))
        self.wait(2)


        #self.play(Transform(plaquette,vertex), run_time=5)
        #self.wait()
    def combineTwoTileAnimation(self,tile1,tile2,time):
        overlap = tile1.makeOverlapTile(tile2,1)
        overlap.move_to(tile2)
        group = VGroup(tile1, tile2)
        self.camera.frame.move_to(group)
        self.play(tile1.animate.move_to(tile2),run_time=time)
        #tile1.getGrid().move_to(tile2)
        self.play(FadeIn(overlap),run_time=time)
