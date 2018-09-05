from matplotlib import pyplot
from numpy.linalg import det
from numpy import array, concatenate, append
from math import pow
from datetime import datetime
import json

class MeshVertex():
    #Main representation of a point; holds info on x,y position and points to one edge connected to it.
    #If the vertex is on the border of a converx hull, it must point to the interiorboundry edge (turning clockwise)
    def __init__(self, x, y, edge=None):
        self._pos = [x,y]
        self._edge_start = edge

    def pos(self):
        return self._pos

    def set_edge(self, edge):
        self._edge_start = edge

    def edge(self):
        return self._edge_start

class GhostVertex():
    #Representation of a point placed at infinity; used to define ghost faces on the outter perimeter of any mesh/hull
    def __init__(self):
        pass

class MeshEdge():
    #Main representation of a segment in the mesh. The edges are the core component of triangles and thus of the graph.
    #They are defined by an opposite edge, a next edge, a destination vertex and an associated faces(triangle).
    #They have to be defined in a geometricaly coherent process for the overwhole mesh to mesh sense.
    #Defining triangles means that the next,next,next edge should be the same edge, and edges of  triangle should always
    #be defined clockwise.
    #Also, a mesh edge should always oppose a mesh edge
    #Edges can represent obstacles or not. If so the pathfinind algorithm will be forbiden to cross though it.
    def __init__(self, vertex_dest, face=None,  opp_edge=None, next_edge=None, is_obstacle=False):
        self._vertex_dest = vertex_dest
        self._face = face
        self._opp_edge = opp_edge
        self._next_edge= next_edge
        self._is_obstacle = is_obstacle

    def dest(self):
        return self._vertex_dest

    def set_dest(self, dest):
        self._vertex_dest = dest

    def set_face(self, face):
        self._face = face

    def get_face(self):
        return self._face

    def set_opp(self, edge):
        self._opp_edge = edge

    def opp(self):
        return self._opp_edge

    def set_next(self, edge):
        self._next_edge = edge

    def next(self):
        return self._next_edge

    def set_obstacle(self, is_obstacle):
        self._is_obstacle = is_obstacle

    def is_obstacle(self):
        return self._is_obstacle

class GhostEdge():
    #Representation of an externat triangle pointing to a vertex at infinity.
    # A gohst edge should always be opposed to a ghost edge in such a way that it allows turning anti-clockwise arond a mesh/hull
    def __init__(self, vertex, face=None, next_edge=None, opp_edge=None):
        self._opp_edge = opp_edge
        self._vertex = vertex
        self._face = face
        self._next_edge= next_edge

    # def __del__(self):
    #     print('ghost edge deleted')

    def set_face(self, face):
        self._face = face

    def get_face(self):
        return self._face

    def set_opp(self, edge):
        self._opp_edge = edge

    def opp(self):
        return self._opp_edge

    def set_next(self, edge):
        self._next_edge = edge

    def next(self):
        return self._next_edge

    def set_vertex(self, v):
        self._vertex = v

    def vertex(self):
        return self._vertex

class MeshFace():
    #Main representation of a triangle. Counter-intuitively, only holds a reference to one of its edges.
    #with that edge as an entry point in the graph from where to retreive the other edges of the triangle it is enough
    def __init__(self, edge=None):
        self._edge = edge

    def set_edge(self, edge):
        self._edge = edge

    def edge(self):
        return self._edge


class GhostFace():
    #Representation of a triangle with a vertex at infinity. Those tiangles are attached on the exterior of a mesh/hull
    #and allow to circulate around in the anti-clockwise direction
    def __init__(self, edge=None):
        self._edge = edge

    # def __del__(self):
    #     print('ghost face deleted')

    def set_edge(self, edge):
        self._edge = edge

    def edge(self):
        return self._edge

def intersect_check(xa,ya, xb,yb, xk,yk, xl, yl):
    #This function checks for an intersection point beteen two segments
    #It takes the coordinates of the (a,b) points bounding the first segments,
    #then the coordinates of (k,l) of the second one and returns either the x,y coordinates of the intersection point
    #or False if no such point was found
    #It dos not handle parallelism and will return False because in that case
    #there would be either no intersection point or an infinity of them

    #parallel vertical segments
    if xa == xb and xk == xl:
        return False

    #intersection of a vertical segments with a non_vertical segment
    #an vertical segment is a special case because the linear factor a of ax+b is then infinite
    #in that case y values should be checked at the vertical x position
    elif xa != xb and xk == xl:
        if min(xa,xb) > xk or max(xa,xb) < xk:
            return False
        else:
            Aa = (ya-yb)/(xa-xb)
            ba = ya-Aa*xa
            Yx = Aa*xk + ba
            if Yx < max(yk,yl) and Yx > min(yk,yl):
                return [xk, Yx]
            else:
                return False
    elif xk != xl and xa == xb:
        if min(xk,xl) > xa or max(xk,xl) < xa:
            return False
        else:
            Ak = (yk-yl)/(xk-xl)
            bk = yk-Ak*xk
            Yx = Ak*xa + bk
            if Yx < max(ya,yb) and Yx > min(ya,yb):
                return [xa, Yx]
            else:
                return False

    #otherwise, both linear parameters are calculated, then an intersection point is found
    #then returned only if it's defined on both segments

    if max(xa,xb) < min(xk,xl) or max(xk, xl) < min(xa,xb):
        return False
    I1 = [min(xa,xb), max(xa,xb)]
    I2 = [min(xk,xl), max(xk,xl)]

    Aa = (ya-yb)/(xa-xb)
    Ak = (yk-yl)/(xk-xl)
    ba = ya-Aa*xa
    bk = yk-Ak*xk

    if Aa == Ak:
        return False

    Xx = (bk-ba) / (Aa-Ak)
    Yx = Aa*Xx + ba

    if (Xx < max(I1[0],I2[0])) or (Xx > min(I1[1], I2[1])):
        return False
    else:
        return[Xx, Yx]

def clockwise_check(a,b,c):
    #This is a geometric test that checks if the angle between vertices a,b and c in that order is defined clockwise
    #It is used to ensure the coherence of the mesh because any triangle in it must be defined clockwise
    #It calculated the area inside the triangle and returns True if positive (points listed clockwise) and False otherwise
    surface_sum = sum([(b.pos()[0]-a.pos()[0])*(b.pos()[1]+a.pos()[1]),
                       (c.pos()[0]-b.pos()[0])*(c.pos()[1]+b.pos()[1]),
                       (a.pos()[0]-c.pos()[0])*(a.pos()[1]+c.pos()[1])])
    return surface_sum > 0

def colinear_check(a,b,c):
    #This is a geometric test that checks weather the surface area defined by three vertices a,b and c is null
    #If so then the three points are colinear and it is impossible to define anthing else than thwo segments
    surface_sum = sum([(b.pos()[0]-a.pos()[0])*(b.pos()[1]+a.pos()[1]),
                       (c.pos()[0]-b.pos()[0])*(c.pos()[1]+b.pos()[1]),
                       (a.pos()[0]-c.pos()[0])*(a.pos()[1]+c.pos()[1])])
    return surface_sum == 0

def circumcircle_check(v1,v2,v3, x):
    #This function checks if a the evaluated vertex x lies inside the circumcircle defined by the vertices a,b,c
    #This problem is solvable by resolving 4/4 matrix determinant
    M = []
    for v in [v1,v2,v3,x]:
        M.append([v.pos()[0],v.pos()[1],pow(v.pos()[0],2)+pow(v.pos()[1],2),1])
    return det(M) > 0

def cone_check(center_pos, b1_pos, b2_pos, x_pos):
    #This is a geometric test that checks weather or not the evaluated point at given x coordinates fits inside
    #a planar cone defined by its center position and two positions of any point on each outline of the cone
    #this function is often used as a directionnality check with

    if (b1_pos[0]-center_pos[0])*(x_pos[1]-center_pos[1]) - (b1_pos[1]-center_pos[1])*(x_pos[0]-center_pos[0]) > 0:
        return False
    if (b2_pos[0]-center_pos[0])*(x_pos[1]-center_pos[1]) - (b2_pos[1]-center_pos[1])*(x_pos[0]-center_pos[0])  < 0:
        return False
    else:
        return True

def half_plane_check(a_pos, b_pos, x_pos):
    #This is a geometric test that check if a point at the x coordinated lies on the upper half plane
    #defined by the directionnal line a,b
    if (b_pos[0]-a_pos[0])*(x_pos[1]-a_pos[1]) - (b_pos[1]-a_pos[1])*(x_pos[0]-a_pos[0]) >= 0:
        return True
    else:
        return False

def inside_triangle_check(a_pos, b_pos, c_pos, x_pos):
    #This is a geometric thest that checks if a point at the x coordinates lies inside a triangle
    #defined by three points a,b,c given in a clockwise order.
    #It is essensially three consecutive half plane checks
    if (b_pos[0]-a_pos[0])*(x_pos[1]-a_pos[1]) - (b_pos[1]-a_pos[1])*(x_pos[0]-a_pos[0]) > 0:
        return False
    elif (c_pos[0]-a_pos[0])*(x_pos[1]-a_pos[1]) - (c_pos[1]-a_pos[1])*(x_pos[0]-a_pos[0]) < 0:
        return  False
    elif (c_pos[0]-b_pos[0])*(x_pos[1]-b_pos[1]) - (c_pos[1]-b_pos[1])*(x_pos[0]-b_pos[0]) > 0:
        return  False
    else:
        return  True

def dist_point_segment(a_pos, b_pos, x_pos):
    v = [b_pos[0]-a_pos[0], b_pos[1]-a_pos[1]]
    w = [x_pos[0]-a_pos[0], x_pos[1]-a_pos[1]]

    c1 = v[0]*w[0] + v[1]*w[1]
    c2 = v[0]*v[0] + v[1]*v[1]
    if c1 <= 0:
        return pow(x_pos[0]-a_pos[0],2) + pow(x_pos[1]-a_pos[1],2)


    elif c2 <= c1:
        return pow(x_pos[0]-b_pos[0],2) + pow(x_pos[1]-b_pos[1],2)
    else:
        b = c1/c2
        np = [a_pos[0]+b*v[0], a_pos[1]+b*v[1]]
        return pow(x_pos[0]-np[0],2) + pow(x_pos[1]-np[1],2)

def on_segment_check(a_pos, b_pos, x_pos):
    #This is a geometric test that checks if a point at the x_pos coordinates lies on the segment defined by a and b
    if b_pos[0]<a_pos[0]:
        a_pos,b_pos = b_pos, a_pos
    if a_pos[0] == b_pos[0]:
        if x_pos[0] == a_pos[0]:
            if x_pos[1] <= max(a_pos[1], b_pos[1]) and x_pos[1] >= min(a_pos[1], b_pos[1]):
                return True
            else:
                return False
        else:
            return False
    else:
        if x_pos[0] <= max(a_pos[0], b_pos[0]) and x_pos[0] >= min(a_pos[0], b_pos[0]):
            Aa = (a_pos[1]-b_pos[1])/(a_pos[0]-b_pos[0])
            ba = a_pos[1]-Aa*a_pos[0]
            Yx = Aa*x_pos[0] + ba
            if Yx == x_pos[1]:
                return True
            else:
                return False
        else:
            return False

def add_ghost_vertex(edge):
    #Completes an existing mesh edge with a ghost vertex, creating 2 ghost edges defined clockwise and a ghost face
    #should be used with an external edge of the mesh/hull and then the opposites edges of the ghost edges should be
    #linked to other outter ghost edges in a ay that allows to turn anti-clockwise around the mesh/hull

    facex = GhostFace(edge)
    edge.set_face(facex)

    a = edge.opp().dest()
    x = GhostVertex()

    edge_bx = GhostEdge(x, face=facex)
    edge.set_next(edge_bx)

    edge_xa = GhostEdge(a, face=facex, next_edge=edge)
    edge_bx.set_next(edge_xa)

    return {'vertex': x, 'face':facex}

def  create_isolated_edge(a,b):
    #Creates the opposite edges connecting two vertices a and b,
    #then adds gosh vertices on both sides and links their opposite sides
    edge_ab = MeshEdge(b)
    edge_ba = MeshEdge(a, opp_edge=edge_ab)
    edge_ab.set_opp(edge_ba)
    a.set_edge(edge_ab)
    b.set_edge(edge_ba)

    # add ghost vertices, edges, and faces to obtain a quad-faced polygon with two ooposite ghost vertices at infinity
    up_ghost = add_ghost_vertex(edge_ab)

    down_ghost = add_ghost_vertex(edge_ba)

    #link the lateral edges of the ghost faces as opposite edges
    edge_ab.next().set_opp(edge_ba.next().next())
    edge_ba.next().next().set_opp(edge_ab.next())
    edge_ab.next().next().set_opp(edge_ba.next())
    edge_ba.next().set_opp(edge_ab.next().next())

    return {'vertices':[up_ghost['vertex'],down_ghost['vertex']],
            'faces':[up_ghost['face'],down_ghost['face']]}

def create_isolated_face(a,b,c):
    #Creates a triangle graph from three vertices a, b and c by adding edges inbetween in a clockwise manner,
    #then creating opposing edges with ghost faces at infinity that allows to turn around the triangle anti-clockwise

    #add a first edge and its associated opposite edge
    edge_ab = MeshEdge(b)
    edge_ba = MeshEdge(a, opp_edge=edge_ab)
    edge_ab.set_opp(edge_ba)

    #define the current face
    face = MeshFace(edge_ab)
    edge_ab.set_face(face)

    #add second and third edges similarely, linking them by their next and opposite attributes
    edge_bc = MeshEdge(c, face=face)
    edge_cb = MeshEdge(b, opp_edge=edge_bc)
    edge_bc.set_opp(edge_cb)
    edge_ab.set_next(edge_bc)

    edge_ca = MeshEdge(a, face=face)
    edge_ac = MeshEdge(c, opp_edge=edge_ca)
    edge_ca.set_opp(edge_ac)
    edge_bc.set_next(edge_ca)
    edge_ca.set_next(edge_ab)

    a.set_edge(edge_ab)
    b.set_edge(edge_bc)
    c.set_edge(edge_ca)

    #add  3 ghost faces all around
    first_ghost = add_ghost_vertex(edge_ba)
    second_ghost = add_ghost_vertex(edge_cb)
    third_ghost = add_ghost_vertex(edge_ac)

    #then associate the ghosts faces lateral edges by their opposite attribute

    edge_ba.next().next().set_opp(edge_cb.next())
    edge_cb.next().set_opp(edge_ba.next().next())
    edge_cb.next().next().set_opp(edge_ac.next())
    edge_ac.next().set_opp(edge_cb.next().next())
    edge_ac.next().next().set_opp(edge_ba.next())
    edge_ba.next().set_opp(edge_ac.next().next())

    return {'vertices':[first_ghost['vertex'],second_ghost['vertex'],third_ghost['vertex']],
            'faces':[face,first_ghost['face'],second_ghost['face'],third_ghost['face']]}

def edge_flip(edge):
    #Flips an edge by modifying the graph locally, in the edge's face and it's opposing face,
    #by assigning that edge and opposite to point toward the prior external summit vertices
    #This results in the edge 'rotating' clockwise, reordering the edges sequences that define both faces

    #This must be used in the context of an edge between two mesh faces, to change the area defined by both faces
    #often used after a negative circumcircle check of the three vertices of one face against the summit of the opposite one

    if type(edge.get_face()) is MeshFace and type(edge.opp().get_face()) is MeshFace:

        #set the faces's edges to point toward the two external edges that will remain the same
        edge.next().next().get_face().set_edge(edge.next().next())
        edge.opp().next().next().get_face().set_edge(edge.opp().next().next())

        #slip the faces pointers of the other two external edges,
        edge.next().set_face(edge.opp().next().next().get_face())
        edge.opp().next().set_face(edge.next().next().get_face())

        #set the destination of the vertices pointed by the central edge and opposite to exterior edges
        #in case they were pointing to any central edge that will flip
        if edge.opp().dest().edge() == edge:
            edge.opp().dest().set_edge(edge.opp().next())
        if edge.dest().edge() == edge.opp():
            edge.dest().set_edge(edge.next())

        #open one side of both triangles by making the last outter edge clockwise
        # point as next to the first outter edge clockwise of the opposite triangle
        # this is because the vertices around the centrad edge will be the new exterior summmit of the flipped faces
        edge.next().next().set_next(edge.opp().next())
        edge.opp().next().next().set_next(edge.next())

        #first flip the opposite center edge by reassigning destination ertex and next edge
        #it can be used as an entry point in the graph to one of the flipped faces
        edge.opp().set_dest(edge.next().next().next().dest())
        edge.opp().set_next(edge.next().next().next().next())

        #then close the other new face by setting it's second edge next edge to be the central edge
        #at this point all edges are accessible by starting at the opposite central edge, then circulation around
        #all four external edges, and finally the central edge.
        edge.next().next().next().set_next(edge)

        #close the first flipped face completely by  flipping the central edge destination and next edge
        edge.set_dest(edge.opp().next().next().dest())
        edge.set_next(edge.opp().next().next().next())

        #close the second flipped face by setting its second ecternal edge's next edge to the fliped oppposite edge
        edge.opp().next().next().set_next(edge.opp())

class MeshSolid():
    def __init__(self, own_vertices, mesh_vertices, edges=[]):
        self._own_vertices = own_vertices
        self._mesh_vertices = mesh_vertices
        if edges != []:
            self._hull_edges = edges
        else:
            #TODO: find all connecting edges from every vertices
            pass
        self._shared_vertices = []

    def get_hull_edges(self):
        return self._hull_edges

    def set_hull_edges(self, edges):
        self._hull_edges = edges

    def add_hull_edge(self, edge):
        self._hull_edges.append(edge)

    def get_own_vertices(self):
        return self._own_vertices

    def set_own_vertices(self, vertices):
        self._own_vertices = vertices

    def add_own_vertex(self, vertex):
        self._own_vertices.append(vertex)

    def get_mesh_vertices(self):
        return self._mesh_vertices

    def set_mesh_vertices(self, vertices):
        self._mesh_vertices = vertices

    def get_shared_vertices(self):
        shared = []
        for shared_vertex in self._shared_vertices:
            shared.append(shared_vertex['vertex'])
        return shared

    def get_sharing_solids(self, vertex):
        for shared in self._shared_vertices:
            if shared['vertex'] == vertex:
                return shared['solids']

    def set_sharing_solids(self, vertex, solids):
        for share_index in list(range(len(self._shared_vertices))):
            if self._shared_vertices[share_index]['vertex'] == vertex:
                self._shared_vertices[share_index]['solids'] = solids


    def add_sharing_solid(self, vertex, solid):
        shared = False
        for share_index in list(range(len(self._shared_vertices))):
            if self._shared_vertices[share_index]['vertex'] == vertex:
                shared = True
                self._shared_vertices[share_index]['solids'].append(solid)
        if not shared:
            self.add_shared_vertex(vertex, [solid])

    def set_shared_vertices(self, shared):
        self._shared_vertices = shared

    def add_shared_vertex(self, vertex, shared_solids):
        self._shared_vertices.append({'vertex': vertex, 'solids':shared_solids})

class MeshMap():
    def __init__(self, mesh={'faces':[],'verties':[],'grid':[]}, grid_len=None):

        if type(mesh) is list:
            self._vertices = []
            self._faces = []
            self._grid = []
            self._grid_len = grid_len
            self.divide_and_conquer(sorted(mesh, key=lambda k: [k[0], k[1]]))
            self._solids = []
        elif type(mesh) is dict and 'faces' in list(mesh.keys()) and 'vertices' in list(mesh.keys()):
            self._vertices = mesh['faces']
            self._faces = mesh['vertices']
            if 'grid' in list(mesh.keys()):
                self._grid = mesh['grid']
            self._solids = []

        if self._vertices != [] and self._grid_len is not None:
            self.grid_sampling()
            self.fill_grid()

    def get_vertices(self):
        return self._vertices
    def get_faces(self):
        return self._faces
    def get_solids(self):
        return self._solids
    def set_solid(self, solids):
        self._solids = solids

    def divide_and_conquer(self, ordered_points):
        #dichotomic function that takes any number of points ordered by x then y,
        #and applies itself to half those points until there are only less than three points
        #in that case, it creates a single face (3 points) or edge (2 points) with their ghost faces
        #add finlly merges the result of this dichotomic process

        #Case of isolated segment
        if len(ordered_points) == 2:
            vertices = []
            faces = []

            #Create 2 mesh vertices from points coordinates
            a = MeshVertex(ordered_points[0][0],ordered_points[0][1])
            b = MeshVertex(ordered_points[1][0],ordered_points[1][1])

            ghosts = create_isolated_edge(a,b)

            #fill up the resulting lists
            vertices.append(a)
            vertices.append(b)
            vertices.append(ghosts['vertices'][0])
            vertices.append(ghosts['vertices'][1])

            faces.append(ghosts['faces'][0])
            faces.append(ghosts['faces'][1])

            self._vertices = vertices
            self._faces = faces

        #case of triangulation
        elif len(ordered_points) == 3:

            vertices = []
            faces = []

            #Create 3 mesh vertices from points coordinates
            a = MeshVertex(ordered_points[0][0],ordered_points[0][1])
            b = MeshVertex(ordered_points[2][0],ordered_points[2][1])
            c = MeshVertex(ordered_points[1][0],ordered_points[1][1])

            if not clockwise_check(a,b,c):

                #TODO: check coherence of the middle point placment and order of segments

                if colinear_check(a,b,c):

                    if not (b.pos()[0]<= max(a.pos()[0], c.pos()[0]) \
                                    and b.pos()[0] >= min(a.pos()[0], c.pos()[0]) \
                                    and b.pos()[1]<= max(a.pos()[1], c.pos()[1]) \
                                    and b.pos()[1] >= min(a.pos()[1], c.pos()[1])):

                        b,c = c,b

                    first_ghosts = create_isolated_edge(a,b)
                    other_ghosts = create_isolated_edge(b,c)


                    first_ghosts['faces'][0].edge().next().set_opp(other_ghosts['faces'][0].edge().next().next())
                    other_ghosts['faces'][0].edge().next().next().set_opp(first_ghosts['faces'][0].edge().next())

                    first_ghosts['faces'][1].edge().next().next().set_opp(other_ghosts['faces'][1].edge().next())
                    other_ghosts['faces'][1].edge().next().set_opp(first_ghosts['faces'][1].edge().next().next())

                    #fill up the resulting lists
                    vertices.append(a)
                    vertices.append(b)
                    vertices.append(c)


                    vertices.append(first_ghosts['vertices'][0])
                    vertices.append(first_ghosts['vertices'][1])
                    vertices.append(other_ghosts['vertices'][0])
                    vertices.append(other_ghosts['vertices'][1])

                    faces.append(first_ghosts['faces'][0])
                    faces.append(first_ghosts['faces'][1])
                    faces.append(other_ghosts['faces'][0])
                    faces.append(other_ghosts['faces'][1])

                    self._vertices = vertices
                    self._faces = faces
                else:
                    b,c = c,b

                    tri = create_isolated_face(a,b,c)

                    vertices.append(a)
                    vertices.append(b)
                    vertices.append(c)
                    vertices.append(tri['vertices'][0])
                    vertices.append(tri['vertices'][1])
                    vertices.append(tri['vertices'][2])

                    faces.append(tri['faces'][0])
                    faces.append(tri['faces'][1])
                    faces.append(tri['faces'][2])
                    faces.append(tri['faces'][3])

                    self._vertices = vertices
                    self._faces = faces

            else:
                tri = create_isolated_face(a,b,c)

                vertices.append(a)
                vertices.append(b)
                vertices.append(c)
                vertices.append(tri['vertices'][0])
                vertices.append(tri['vertices'][1])
                vertices.append(tri['vertices'][2])

                faces.append(tri['faces'][0])
                faces.append(tri['faces'][1])
                faces.append(tri['faces'][2])
                faces.append(tri['faces'][3])

                self._vertices = vertices
                self._faces = faces

        #case of overpopulation of points; halfen the list then call recursively
        else:
            Lmesh = MeshMap(ordered_points[:len(ordered_points)//2])
            Rmesh = MeshMap(ordered_points[len(ordered_points)//2:])

            merged =merge_meshes(Lmesh,Rmesh)

            self._vertices = merged['vertices']
            self._faces = merged['faces']

    def grid_sampling(self):
        #generates a square grid that holds references to vertices as entry points in the graph
        #it is usefull to quistart the search to find in wich triangle lies a point of coordinates x,y
        #the resolution depends on the square_len and the functions associate each square with any vertex it contains
        #the final grid has one vertex by square, zero if there is none
        sampled_vertices = [[None]]
        for vertex in self._vertices:
            if type(vertex) is MeshVertex:
                i = int(vertex.pos()[0]//self._grid_len)
                j = int(vertex.pos()[1]//self._grid_len)
                if i >= len(sampled_vertices):
                    x_extend = []
                    for k in list(range(len(sampled_vertices[-1]))):
                        x_extend.append(None)
                    for extend_len in list(range(i-len(sampled_vertices)+1)):
                        sampled_vertices.append(x_extend)
                if j >= len(sampled_vertices[-1]):
                    for l in list(range(len(sampled_vertices))):
                        for extend_len in list(range(j- len(sampled_vertices[l])+1)):
                            sampled_vertices[l].append(None)
                if sampled_vertices[i][j] is None:
                    sampled_vertices[i][j] = vertex
        self._grid = sampled_vertices

    def fill_grid(self):
        #an expansion function that fills the grid whenever the value of a square isn't set
        #it will then copy the vertex reference of its closest neighbour until all the grid is filled
        grid_filled = False
        while not grid_filled:
            none_remaining = False
            for i in list(range(len(self._grid))):
                for j in list(range(len(self._grid[i]))):
                    if self._grid[i][j] is None:
                        closest_sample = None
                        if i == 0:
                            ibound = [+1]
                        elif  i == len(self._grid) -1:
                            ibound = [-1]
                        else:
                            ibound = [+1, -1]
                        if j == 0:
                            jbound = [+1]
                        elif  j == len(self._grid[i]) -1:
                            jbound = [-1]
                        else:
                            jbound = [+1, -1]

                        for di in ibound:
                            for dj in jbound:
                                if self._grid[i+di][j+dj] is not None:
                                    if closest_sample is None:
                                        closest_sample = self._grid[i+di][j+dj]
                                    else:
                                        if pow(self._grid[i+di][j+dj].pos()[0]-self._grid_len*(i+1/2), 2) + pow(self._grid[i+di][j+dj].pos()[1]-self._grid_len*(j+1/2), 2) \
                                                < pow(closest_sample.pos()[0]-self._grid_len*(i+1/2),2)+pow(closest_sample.pos()[1]-self._grid_len*(j+1/2),2):
                                            closest_sample = self._grid[i+di][j+dj]
                        if closest_sample is None:
                            none_remaining = True
                        else:
                            self._grid[i][j] = closest_sample
            if not none_remaining:
                grid_filled = True

    def find_face(self, x_pos, start_vertex=None, face=None):
        #given a point's coordinates x,y, find the face in the mesh that point lies in
        #must be given a start vertex as an entry point in the graph
        # (the closer the entry point, thefaster the face will be found)
        if start_vertex is None:
            if self._grid != [] and self._grid_len is not None:
                start_vertex = [int(x_pos[0]//self._grid_len)][int(x_pos[1]//self._grid_len)]
            else:
                start_vertex = self._vertices[0]

        if face is None:
            #if face is not set, pick any face around the starting vertex
            face = start_vertex.edge().get_face()
            return self.find_face(x_pos, start_vertex, face)
        else:
            #TODO: check if the x_pos is a vertex position or on an edge!!
            if type(face) is GhostFace:
                border_edge = face.edge()
                if border_edge.opp().dest() == start_vertex:
                    if not cone_check(start_vertex.pos(), border_edge.dest().pos(),
                                      border_edge.next().next().opp().next().next().opp().dest().pos() , x_pos):
                        return self.find_face(x_pos, start_vertex, border_edge.next().next().opp().next().next().opp().get_face())
                    else:
                        return face
                elif border_edge.dest() == start_vertex:
                    if not cone_check(start_vertex.pos(), border_edge.next().opp().next().opp().dest().pos(),
                                      border_edge.opp().dest().pos() , x_pos):
                        return self.find_face(x_pos, start_vertex, border_edge.opp().get_face())
                    else:
                        return face
                else:
                    print('oopsys')
            else:

                #identify the edge of the current face at the opposite of the starting vertex
                crossing_edge = face.edge()
                if crossing_edge.dest() == start_vertex:
                    crossing_edge = crossing_edge.next().next()
                elif crossing_edge.next().next().dest() == start_vertex:
                    crossing_edge = crossing_edge.next()

                # print('yo')
                # ar = display_faces(self._faces, 'r')
                # tx = []
                # ty = []
                # tx.append(start_vertex.pos()[0])
                # tx.append(x_pos[0])
                # ty.append(start_vertex.pos()[1])
                # ty.append(x_pos[1])
                # ar.plot(tx,ty,'g')
                # ex = []
                # ey = []
                # ex.append(crossing_edge.opp().dest().pos()[0])
                # ex.append(crossing_edge.dest().pos()[0])
                # ey.append(crossing_edge.opp().dest().pos()[1])
                # ey.append(crossing_edge.dest().pos()[1])
                # ar.plot(ex,ey,'k')
                # px,py = [],[]
                # px.append(crossing_edge.next().opp().dest().pos()[0])
                # px.append(crossing_edge.next().dest().pos()[0])
                # py.append(crossing_edge.next().opp().dest().pos()[1])
                # py.append(crossing_edge.next().dest().pos()[1])
                # ar.plot(px,py,'b')
                # # pyplot.show()
                # lx,ly = [],[]
                # lx.append(crossing_edge.next().next().opp().dest().pos()[0])
                # lx.append(crossing_edge.next().next().dest().pos()[0])
                # ly.append(crossing_edge.next().next().opp().dest().pos()[1])
                # ly.append(crossing_edge.next().next().dest().pos()[1])
                # ar.plot(lx,ly,'m')
                # print(crossing_edge.next().next().opp().dest().pos(), crossing_edge.next().next().dest().pos())
                # pyplot.show()
                #verify if the point to localize is in the cone defined by the current face, with the starting vertex as summit
                if cone_check(start_vertex.pos(), crossing_edge.next().next().dest().pos(), crossing_edge.dest().pos(), x_pos):
                    #if the point to localize is in the direction of the cone, it can then be either
                    #inside the triangle or outside and the line start/end will intersect the crossing edge
                    if intersect_check(start_vertex.pos()[0], start_vertex.pos()[1], x_pos[0], x_pos[1],
                                       crossing_edge.opp().dest().pos()[0],
                                       crossing_edge.opp().dest().pos()[1],
                                       crossing_edge.dest().pos()[0],
                                       crossing_edge.dest().pos()[1]):
                        #if the edge is crossed, then the function falls itself
                        # after changing the face to be the one on the other side of the crossed edge,
                        #thus getting the start vertex closer to the end vertex
                        if on_segment_check(crossing_edge.opp().dest().pos(), crossing_edge.dest().pos(), x_pos):
                            return face
                        else:
                            if type(crossing_edge.opp().get_face()) is MeshFace:
                                return self.find_face(x_pos, crossing_edge.opp().next().dest(), None)
                            else:
                                return crossing_edge.opp().get_face()
                    else:
                        #if the point to localize is in the cone and no interection happens,
                        #it means that it should be inside the face, wich will be the reslut to return
                        if inside_triangle_check(face.edge().dest().pos(),
                                                 face.edge().next().dest().pos(),
                                                 face.edge().next().next().dest().pos(), x_pos):
                            return face
                        else:
                            if on_segment_check(crossing_edge.dest().pos(), crossing_edge.next().dest().pos(), x_pos) \
                                    or on_segment_check(crossing_edge.next().dest().pos(), crossing_edge.next().next().opp().dest().pos(), x_pos):
                                return face
                            else:
                                print('oops')
                else:
                    #if the point to localize is not in the cone as defined, it means that the start/end line
                    #will go though another face around the start vertex, so the function
                    #is called again with a face that is the next face by rotation around the starting vertex
                    return self.find_face(x_pos, start_vertex, crossing_edge.next().opp().get_face())

    def vertex_insertion(self, vertex, face, epsilon=1):
        #modifies the graph to insert a vertex in the mesh, given the face it will be located in

        # ar = display_faces(self._faces, 'r')
        # tx = []
        # ty = []
        # tx.append(vertex.pos()[0])
        # tx.append(face.edge().dest().pos()[0])
        # ty.append(vertex.pos()[1])
        # ty.append(face.edge().dest().pos()[1])
        # ar.plot(tx,ty,'g')
        # ex = []
        # ey = []
        # ex.append(face.edge().opp().dest().pos()[0])
        # ex.append(face.edge().dest().pos()[0])
        # ey.append(face.edge().opp().dest().pos()[1])
        # ey.append(face.edge().dest().pos()[1])
        # ar.plot(ex,ey,'k')
        # px,py = [],[]
        # px.append(face.edge().next().opp().dest().pos()[0])
        # px.append(face.edge().next().dest().pos()[0])
        # py.append(face.edge().next().opp().dest().pos()[1])
        # py.append(face.edge().next().dest().pos()[1])
        # ar.plot(px,py,'b')
        # # pyplot.show()
        # lx,ly = [],[]
        # lx.append(face.edge().next().next().opp().dest().pos()[0])
        # lx.append(face.edge().next().next().dest().pos()[0])
        # ly.append(face.edge().next().next().opp().dest().pos()[1])
        # ly.append(face.edge().next().next().dest().pos()[1])
        # ar.plot(lx,ly,'m')
        # pyplot.show()

        superposition = False
        e = face.edge()
        i = 0
        while not superposition and i < 3:
            if pow(e.dest().pos()[0] - vertex.pos()[0],2) + pow(e.dest().pos()[1] - vertex.pos()[1],2) < pow(epsilon,2):
                return e.dest()
            e = e.next()
            i+=1

        if not superposition:
            e = face.edge()
            i = 0
            while not superposition and i < 3 :
                if on_segment_check(e.opp().dest().pos(), e.dest().pos(), vertex.pos()) \
                        or dist_point_segment(e.opp().dest().pos(), e.dest().pos(), vertex.pos()) < pow(epsilon,2):
                    # split the edge in two and the two opposing faces in four (each face is split in two)
                    superposition = True
                    face.set_edge(e.next().next())
                    new_face = MeshFace(e.next())
                    e.next().set_face(new_face)
                    split_new_edge = MeshEdge(e.dest(), new_face, opp_edge=e.opp(), next_edge=e.next())
                    vertex.set_edge(split_new_edge)
                    if e.is_obstacle():
                        split_new_edge.set_obstacle(True)
                        for solid in self.get_solids():
                            if e in solid.get_hull_edges():
                                solid.add_hull_edge(split_new_edge)
                                solid.add_own_vertex(vertex)

                    outward_new_edge = MeshEdge(e.next().dest(), face, next_edge=e.next().next())
                    inward_new_edge = MeshEdge(vertex, new_face, opp_edge=outward_new_edge, next_edge=split_new_edge)
                    outward_new_edge.set_opp(inward_new_edge)
                    e.next().set_next(inward_new_edge)
                    e.set_dest(vertex)
                    e.set_next(outward_new_edge)

                    e.opp().get_face().set_edge(e.opp().next().next())
                    opp_face = MeshFace(e.opp().next())
                    e.opp().next().set_face(opp_face)
                    outward_opp_edge = MeshEdge(e.opp().next().dest(), e.opp().get_face(), next_edge=e.opp().next().next())
                    inward_opp_edge = MeshEdge(vertex, opp_face, opp_edge= outward_opp_edge)
                    outward_opp_edge.set_opp(inward_opp_edge)
                    e.opp().next().set_next(inward_opp_edge)
                    split_opp_edge = MeshEdge(e.opp().dest(), opp_face, opp_edge=e, next_edge=e.opp().next())
                    inward_opp_edge.set_next(split_opp_edge)
                    if e.opp().is_obstacle():
                        split_opp_edge.set_obstacle(True)
                        for solid in self.get_solids():
                            if e.opp() in solid.get_hull_edges():
                                solid.add_hull_edge(split_opp_edge)
                                solid.add_own_vertex(vertex)

                    e.opp().set_dest(vertex)
                    e.opp().set_next(outward_opp_edge)
                    e.opp().set_opp(split_new_edge)
                    e.set_opp(split_opp_edge)

                    self._faces.insert(self._faces.index(face), new_face)
                    self._faces.insert(self._faces.index(face), opp_face)

                    vertex_index = 0
                    index_found = False
                    while vertex_index <= len(self._vertices) and not index_found:
                        if vertex.pos()[0] < self._vertices[vertex_index].pos()[0]:
                            vertex_index += 1
                        else:
                            if vertex.pos()[1] < self._vertices[vertex_index].pos()[1]:
                                vertex_index += 1
                            else:
                                index_found = True
                    self._vertices.insert(vertex_index, vertex)

                    #after point insertion, check if the delaunay criterion is preserved in those three faces;
                    # if not, flip the external edge
                    ext_edge = vertex.edge().next()
                    for fnum in [0,1,2,3]:
                        if type(ext_edge.opp().get_face()) is MeshFace:
                            if not circumcircle_check(ext_edge.dest(), ext_edge.next().dest(), ext_edge.opp().dest(), ext_edge.opp().next().dest()) \
                                    and not ext_edge.is_obstacle():
                                edge_flip(ext_edge)
                                ext_edge = ext_edge.opp().next().next().opp().next()
                            else:
                                ext_edge = ext_edge.next().opp().next()

                e = e.next()
                i+=1

        if not superposition:
            #assuming the point is indeed located in the given face, then that face should be split in three faces
            #so there will be two new faces and three new edges and opposites
            # linking the new inserted vertex with every vertex summit of the face it is inserted in

            #the entry point will be any edge of the face that face is pointing to
            #the other two edges are associted with newly created faces
            first_new_face = MeshFace(face.edge().next().next())
            face.edge().next().next().set_face(first_new_face)
            second_new_face = MeshFace(face.edge().next())
            face.edge().next().set_face(second_new_face)

            #add the first new and opposite edge after the third edge of the face
            first_new_edge = MeshEdge(face.edge().opp().dest(),face, next_edge=face.edge())
            first_opp_edge = MeshEdge(vertex, first_new_face, opp_edge=first_new_edge)
            first_new_edge.set_opp(first_opp_edge)
            vertex.set_edge(first_new_edge)
            face.edge().next().next().set_next(first_opp_edge)

            #add the second new edge and opposite after the second edge of the face, then close the first new face
            second_new_edge = MeshEdge(face.edge().next().dest(), first_new_face,  next_edge=face.edge().next().next())
            first_opp_edge.set_next(second_new_edge)
            second_opp_edge = MeshEdge(vertex, second_new_face, opp_edge=second_new_edge)
            second_new_edge.set_opp(second_opp_edge)
            face.edge().next().set_next(second_opp_edge)

            #add the third new edge and opposite after the first edge of the face, then close the second and third new faces
            third_new_edge = MeshEdge(face.edge().dest(), second_new_face, next_edge=face.edge().next())
            second_opp_edge.set_next(third_new_edge)
            third_opp_edge = MeshEdge(vertex, face, opp_edge=third_new_edge, next_edge= first_new_edge)
            face.edge().set_next(third_opp_edge)
            third_new_edge.set_opp(third_opp_edge)

            #insert the new faces and the new vertex in the mesh where they should be
            vertex_index = 0
            index_found = False
            while vertex_index <= len(self._vertices) and not index_found:
                if vertex.pos()[0] < self._vertices[vertex_index].pos()[0]:
                    vertex_index += 1
                else:
                    if vertex.pos()[1] < self._vertices[vertex_index].pos()[1]:
                        vertex_index += 1
                    else:
                        index_found = True
            self._vertices.insert(vertex_index, vertex)
            self._faces.insert(self._faces.index(face), first_new_face)
            self._faces.insert(self._faces.index(face), second_new_face)

            #after point insertion, check if the delaunay criterion is preserved in those three faces;
            # if not, flip the external edge
            ext_edge = vertex.edge().next()
            for fnum in [0,1,2]:
                if type(ext_edge.opp().get_face()) is MeshFace:
                    if not circumcircle_check(ext_edge.dest(), ext_edge.next().dest(), ext_edge.opp().dest(), ext_edge.opp().next().dest()) \
                            and not ext_edge.is_obstacle():
                        edge_flip(ext_edge)
                        ext_edge = ext_edge.opp().next().next().opp().next()
                    else:
                        ext_edge = ext_edge.next().opp().next()

    def vertex_removal(self, vertex):
        #TODO:add colinearity checks and handle edge merge
        flipped = True
        start_edge = vertex.edge()
        current_edge = start_edge.next().next().opp()

        while not current_edge.next().next().opp().next().next().opp().next().next().opp() == current_edge and flipped:
            flipped = False
            if type(current_edge.get_face()) is MeshFace and type(current_edge.next().next().opp().get_face()) is MeshFace:
                while current_edge != start_edge:
                    if type(current_edge.get_face()) is MeshFace and type(current_edge.next().next().opp().get_face()) is MeshFace:
                        if not inside_triangle_check(current_edge.dest().pos(), current_edge.next().next().opp().dest().pos(), current_edge.next().next().opp().next().next().opp().dest().pos(), vertex.pos())\
                            and not inside_triangle_check(vertex.pos(), current_edge.dest().pos(), current_edge.next().next().opp().next().next().opp().dest().pos(), current_edge.next().next().opp().dest().pos()) \
                            and not current_edge.next().next().opp().is_obstacle() and not current_edge.next().next().is_obstacle():
                            if not on_segment_check(current_edge.dest().pos(), current_edge.next().next().opp().next().next().opp().dest().pos(), current_edge.next().next().opp().dest().pos())\
                                and not on_segment_check(current_edge.dest().pos(), current_edge.next().next().opp().next().next().opp().dest().pos(), vertex.pos()):
                                flipped = True

                                # plt = display_faces([current_edge.get_face(), current_edge.next().next().opp().get_face(), current_edge.next().next().opp().next().next().opp().get_face()])

                                current_edge = current_edge.next().next().opp()
                                edge_flip(current_edge)

                                # plt = display_faces([current_edge.get_face(), current_edge.next().next().opp().get_face(), current_edge.next().next().opp().next().next().opp().get_face()])
                                # plt = display_faces([current_edge.next().opp().get_face(), current_edge.next().opp().next().next().opp().get_face(), current_edge.next().opp().next().next().opp().next().next().opp().get_face()])

                                start_edge = current_edge.next().next()
                                current_edge = current_edge.next().opp()
                                # if disp:
                                #     display_map(self)
                            else:
                                current_edge = current_edge.next().next().opp()
                        else:
                            current_edge = current_edge.next().next().opp()
                    else:
                        current_edge = current_edge.next().next().opp()
            else:
                current_edge = current_edge.next().next().opp()

        if type(current_edge.get_face()) is MeshFace and type(current_edge.next().next().opp().get_face()) is MeshFace and type(current_edge.next().next().opp().next().next().opp().get_face()) is MeshFace:
            if inside_triangle_check(current_edge.dest().pos(), current_edge.next().next().opp().dest().pos(), current_edge.next().next().opp().next().next().opp().dest().pos(), vertex.pos())\
                and current_edge.next().next().opp().next().next().opp().next().next().opp() == current_edge:
                out_edge = current_edge.next()
                first_edge = current_edge
                second_edge = current_edge.next().next().opp()
                third_edge = current_edge.next().next().opp().next().next().opp()
                vertex.set_edge(None)
                current_edge.get_face().set_edge(current_edge.next())
                for edge in [first_edge, second_edge, third_edge]:
                    if edge.dest().edge() == edge.opp():
                        edge.dest().set_edge(edge.next())
                    edge.next().set_face(current_edge.get_face())

                # for edge in [current_edge.next().next().opp().next().next().opp(), current_edge.next().next().opp(), current_edge]:
                for edge in [first_edge, second_edge, third_edge]:
                    edge.opp().next().next().set_next(edge.next())

                # current_edge.next().next().opp().next().next().set_next(current_edge.opp().next().next())
                # current_edge.next().set_next(current_edge.next().next().opp().next())
                # current_edge.opp().next().next().set_next(current_edge.next())

                second_edge.get_face().set_edge(None)
                self._faces.remove(second_edge.get_face())

                third_edge.get_face().set_edge(None)
                self._faces.remove(third_edge.get_face())

                for edge in [first_edge, second_edge, third_edge]:
                    edge.set_dest(None)
                    edge.opp().set_dest(None)
                    edge.set_face(None)
                for edge in [first_edge, second_edge, third_edge]:
                    edge.set_next(None)
                    edge.opp().set_next(None)
                    edge.opp().set_face(None)
                    edge.opp().set_opp(None)
                    # del edge.opp()
                    edge.set_opp(None)
                    # del edge
                self._vertices.remove(vertex)
                first_edge = out_edge
                second_edge = out_edge.next()
                third_edge = out_edge.next().next()

                for edge in [first_edge, second_edge, third_edge]:

                    if type(edge.get_face()) is MeshFace and type(edge.opp().get_face()) is MeshFace:
                        if not circumcircle_check(edge.dest(), edge.opp().dest(), edge.opp().next().dest(), edge.next().dest()) \
                                and not edge.is_obstacle() and not edge.opp().is_obstacle():
                            edge_flip(edge)

                # if disp:
                #     display_map(self)

            elif current_edge.next().next().opp().next().next().opp().next().next().opp().next().next().opp() == current_edge \
                    and type(current_edge.next().next().opp().next().next().opp().next().next().opp().get_face()) is MeshFace:
                if on_segment_check(current_edge.dest().pos(), current_edge.next().next().opp().next().next().opp().dest().pos(), vertex.pos()):
                    first_edge = current_edge.next().next().opp()
                    second_edge = current_edge.next().next().opp().next().next().opp().next().next().opp()
                    out_edge = current_edge.next()

                    for edge in [first_edge, second_edge]:
                        edge.get_face().set_edge(edge.next())
                        edge.opp().next().next().set_face(edge.get_face())
                        if edge.dest().edge() == edge.opp():
                            edge.dest().set_edge(edge.next())
                        edge.opp().next().next().set_next(edge.next())
                        if edge.opp().next().dest().edge() == edge.opp().next().opp():
                            edge.opp().next().dest().set_edge(edge.opp().next().next())
                    for edge in [first_edge, second_edge]:
                        edge.opp().next().get_face().set_edge(None)
                        self._faces.remove(edge.opp().next().get_face())
                        edge.opp().next().set_face(None)
                        edge.next().next().set_dest(edge.opp().next().dest())
                        edge.next().next().set_next(edge.opp().next().next())
                        edge.next().next().set_opp(edge.opp().next().opp())
                        edge.opp().next().set_dest(None)
                        edge.opp().next().set_next(None)
                        edge.opp().next().set_opp(None)

                    for edge in [first_edge, second_edge]:

                        edge.opp().set_dest(None)
                        edge.set_dest(None)
                        edge.opp().set_face(None)
                        edge.set_face(None)
                        edge.opp().set_next(None)
                        edge.set_next(None)
                        edge.opp().set_opp(None)
                        edge.set_opp(None)

                    self._vertices.remove(vertex)

                    if type(out_edge.next().next().get_face()) is MeshFace and type(out_edge.next().next().opp().get_face()) is MeshFace:
                        if not circumcircle_check(out_edge.opp().dest(), out_edge.dest(), out_edge.next().dest(), out_edge.next().next().opp().next().dest()) \
                                and not out_edge.next().next().is_obstacle() and not out_edge.next().next().opp().is_obstacle():
                            edge_flip(out_edge.next().next())

                    # if disp:
                    #     display_map(self)

                elif on_segment_check(current_edge.next().next().opp().next().next().opp().dest().pos(), current_edge.next().next().opp().next().next().opp().next().next().opp().dest().pos(), vertex.pos()):

                    current_edge = current_edge.next().next().opp()

                    first_edge = current_edge.next().next().opp()
                    second_edge = current_edge.next().next().opp().next().next().opp().next().next().opp()
                    out_edge = current_edge.next()

                    for edge in [first_edge, second_edge]:
                        edge.get_face().set_edge(edge.next())
                        edge.opp().next().next().set_face(edge.get_face())
                        if edge.dest().edge() == edge.opp():
                            edge.dest().set_edge(edge.next())
                        edge.opp().next().next().set_next(edge.next())
                        if edge.opp().next().dest().edge() == edge.opp().next().opp():
                            edge.opp().next().dest().set_edge(edge.opp().next().next())
                    for edge in [first_edge, second_edge]:
                        edge.opp().next().get_face().set_edge(None)
                        self._faces.remove(edge.opp().next().get_face())
                        edge.opp().next().set_face(None)
                        edge.next().next().set_dest(edge.opp().next().dest())
                        edge.next().next().set_next(edge.opp().next().next())
                        edge.next().next().set_opp(edge.opp().next().opp())
                        edge.opp().next().set_dest(None)
                        edge.opp().next().set_next(None)
                        edge.opp().next().set_opp(None)

                    for edge in [first_edge, second_edge]:

                        edge.opp().set_dest(None)
                        edge.set_dest(None)
                        edge.opp().set_face(None)
                        edge.set_face(None)
                        edge.opp().set_next(None)
                        edge.set_next(None)
                        edge.opp().set_opp(None)
                        edge.set_opp(None)

                    self._vertices.remove(vertex)

                    if type(out_edge.next().next().get_face()) is MeshFace and type(out_edge.next().next().opp().get_face()) is MeshFace:
                        if not circumcircle_check(out_edge.opp().dest(), out_edge.dest(), out_edge.next().dest(), out_edge.next().next().opp().next().dest()) \
                                and not out_edge.next().next().is_obstacle() and not out_edge.next().next().opp().is_obstacle():
                            edge_flip(out_edge.next().next())

                    # if disp:
                    #     display_map(self)

                elif on_segment_check(current_edge.dest().pos(), current_edge.next().next().opp().next().next().opp().dest().pos(), current_edge.next().next().opp().dest().pos()):
                    current_edge = current_edge.next().next().opp()

                    first_edge = current_edge.next().next().opp()
                    second_edge = current_edge.next().next().opp().next().next().opp().next().next().opp()
                    out_edge = current_edge.next()

                    for edge in [first_edge, second_edge]:
                        edge.get_face().set_edge(edge.next())
                        edge.opp().next().next().set_face(edge.get_face())
                        if edge.dest().edge() == edge.opp():
                            edge.dest().set_edge(edge.next())
                        edge.opp().next().next().set_next(edge.next())
                        if edge.opp().next().dest().edge() == edge.opp().next().opp():
                            edge.opp().next().dest().set_edge(edge.opp().next().next())
                    for edge in [first_edge, second_edge]:
                        edge.opp().next().get_face().set_edge(None)
                        self._faces.remove(edge.opp().next().get_face())
                        edge.opp().next().set_face(None)
                        edge.next().next().set_dest(edge.opp().next().dest())
                        edge.next().next().set_next(edge.opp().next().next())
                        edge.next().next().set_opp(edge.opp().next().opp())
                        edge.opp().next().set_dest(None)
                        edge.opp().next().set_next(None)
                        edge.opp().next().set_opp(None)

                    for edge in [first_edge, second_edge]:

                        edge.opp().set_dest(None)
                        edge.set_dest(None)
                        edge.opp().set_face(None)
                        edge.set_face(None)
                        edge.opp().set_next(None)
                        edge.set_next(None)
                        edge.opp().set_opp(None)
                        edge.set_opp(None)

                    self._vertices.remove(vertex)

                    if type(out_edge.next().next().get_face()) is MeshFace and type(out_edge.next().next().opp().get_face()) is MeshFace:
                        if not circumcircle_check(out_edge.opp().dest(), out_edge.dest(), out_edge.next().dest(), out_edge.next().next().opp().next().dest()) \
                                and not out_edge.next().next().is_obstacle() and not out_edge.next().next().opp().is_obstacle():
                            edge_flip(out_edge.next().next())

                    # if disp:
                    #     display_map(self)

                elif on_segment_check(current_edge.next().next().opp().next().next().opp().dest().pos(), current_edge.next().next().opp().next().next().opp().next().next().opp().dest().pos(), current_edge.next().next().opp().next().next().opp().dest().pos()):

                    current_edge = current_edge.next().next().opp().next().next().opp()
                    first_edge = current_edge.next().next().opp()
                    second_edge = current_edge.next().next().opp().next().next().opp().next().next().opp()
                    out_edge = current_edge.next()

                    for edge in [first_edge, second_edge]:
                        edge.get_face().set_edge(edge.next())
                        edge.opp().next().next().set_face(edge.get_face())
                        if edge.dest().edge() == edge.opp():
                            edge.dest().set_edge(edge.next())
                        edge.opp().next().next().set_next(edge.next())
                        if edge.opp().next().dest().edge() == edge.opp().next().opp():
                            edge.opp().next().dest().set_edge(edge.opp().next().next())
                    for edge in [first_edge, second_edge]:
                        edge.opp().next().get_face().set_edge(None)
                        self._faces.remove(edge.opp().next().get_face())
                        edge.opp().next().set_face(None)
                        edge.next().next().set_dest(edge.opp().next().dest())
                        edge.next().next().set_next(edge.opp().next().next())
                        edge.next().next().set_opp(edge.opp().next().opp())
                        edge.opp().next().set_dest(None)
                        edge.opp().next().set_next(None)
                        edge.opp().next().set_opp(None)

                    for edge in [first_edge, second_edge]:

                        edge.opp().set_dest(None)
                        edge.set_dest(None)
                        edge.opp().set_face(None)
                        edge.set_face(None)
                        edge.opp().set_next(None)
                        edge.set_next(None)
                        edge.opp().set_opp(None)
                        edge.set_opp(None)

                    self._vertices.remove(vertex)

                    if type(out_edge.next().next().get_face()) is MeshFace and type(out_edge.next().next().opp().get_face()) is MeshFace:
                        if not circumcircle_check(out_edge.opp().dest(), out_edge.dest(), out_edge.next().dest(), out_edge.next().next().opp().next().dest()) \
                                and not out_edge.next().next().is_obstacle() and not out_edge.next().next().opp().is_obstacle():
                            edge_flip(out_edge.next().next())

                    # if disp:
                    #     display_map(self)

        # if disp:
        #     display_map(self)
        #     pyplot.show()

    def delete_solid(self, solid):
        #TODO: check vertex edges on shared solid edges, orient them properly on solid removal
        for edge in solid.get_hull_edges():
            edge.set_obstacle(False)
        solid.set_hull_edges(None)
        remain = []
        for own_vertex in solid.get_own_vertices():
            if own_vertex not in solid.get_shared_vertices():
                self.vertex_removal(own_vertex)
        solid.set_own_vertices(None)
        solid.set_mesh_vertices(None)
        for shared_vertex in solid.get_shared_vertices():
            for sharing_solid in solid.get_sharing_solids(shared_vertex):
                sharing_remain = sharing_solid.get_sharing_solids(shared_vertex).remove(solid)
                if sharing_remain is not None:
                    sharing_solid.set_sharing_solids(shared_vertex, sharing_remain)
                else:
                    sharing_solid.get_shared_vertices().remove(shared_vertex)
                    if sharing_solid.get_shared_vertices() is None:
                        sharing_solid.set_shared_vertices([])
        solid.set_shared_vertices(None)
        self._solids.remove(solid)

    def edge_slice(self, start_vertex, end_vertex, face=None, epsilon=1, obstacle=True, slice={'added_vertices':[],'mesh_vertices':[],'edges':[]}):

        #modifies the graph to create a straight line between two vertices
        #it will also add a vertex on each crossed edge, effectively slicing each crossed face in three
        #the obstacle parameter can be set to give the obstacle property to ever created edge
        #the obstacle property is directionnal/oriened so any object should be a hull of edge slices
        # between vertices given in clockwise order

        #the fonction's structure resembles the one of find_face() expect it is not as passive
        #it will go though the graph the same way but will edit it on the way
        if face is None:
            face = start_vertex.edge().get_face()
            return self.edge_slice(start_vertex, end_vertex, face, epsilon, obstacle, slice)

        else:
            if type(face) is GhostFace:
                border_edge = face.edge()
                if border_edge.opp().dest() == start_vertex:
                    if not cone_check(start_vertex.pos(), border_edge.dest().pos(),
                                      border_edge.next().next().opp().next().next().opp().dest().pos() , end_vertex.pos()):
                        return self.edge_slice(start_vertex, end_vertex, border_edge.next().next().opp().next().next().opp().get_face(), epsilon, obstacle, slice)
                elif border_edge.dest() == start_vertex:
                    if not cone_check(start_vertex.pos(), border_edge.next().opp().next().opp().dest().pos(),
                                      border_edge.opp().dest().pos() , end_vertex.pos()):
                        return self.edge_slice(start_vertex, end_vertex, border_edge.opp().get_face(), epsilon, obstacle, slice)
                else:
                    print('oopsys')
            else:
                #if the current face contains the end vertex, one of its edges should be the last edge of the slice
                end_edge = face.edge()
                end_found = False
                i = 0
                while i < 3 and not end_found:
                    if end_edge.dest().pos() == end_vertex.pos():
                        if end_edge.opp().dest() == start_vertex:
                            end_found = True
                        elif end_edge.next().dest() == start_vertex:
                            end_edge = end_edge.next().opp()
                            end_found = True
                    else:
                        end_edge = end_edge.next()
                    i += 1

                if end_found:
                    if not end_edge.opp().is_obstacle():
                        start_vertex.set_edge(end_edge)
                    slice['edges'].append(end_edge)
                    if obstacle:
                        end_edge.set_obstacle(True)

                    # ar = display_faces(self._faces, 'r')
                    # tx = []
                    # ty = []
                    # tx.append(start_vertex.pos()[0])
                    # tx.append(end_vertex.pos()[0])
                    # ty.append(start_vertex.pos()[1])
                    # ty.append(end_vertex.pos()[1])
                    # ar.plot(tx,ty,'g')
                    # ex = []
                    # ey = []
                    # ex.append(end_edge.opp().dest().pos()[0])
                    # ex.append(end_edge.dest().pos()[0])
                    # ey.append(end_edge.opp().dest().pos()[1])
                    # ey.append(end_edge.dest().pos()[1])
                    # ar.plot(ex,ey,'b')
                    # pyplot.show()

                    return slice

                else:
                    #TODO: check if the target line is on an edge
                    #identify the external edge opposed to the current start vertex
                    crossing_edge = face.edge()
                    if crossing_edge.dest() == start_vertex:
                        crossing_edge = crossing_edge.next().next()
                    elif crossing_edge.next().next().dest() == start_vertex:
                        crossing_edge = crossing_edge.next()
                    #check if the diretion of the end vertex is in the cone defined by the current face
                    if cone_check(start_vertex.pos(), crossing_edge.next().next().dest().pos(), crossing_edge.dest().pos(), end_vertex.pos()):

                        intersection = intersect_check(start_vertex.pos()[0], start_vertex.pos()[1],
                                                       end_vertex.pos()[0], end_vertex.pos()[1],
                                                       crossing_edge.opp().dest().pos()[0],
                                                       crossing_edge.opp().dest().pos()[1],
                                                       crossing_edge.dest().pos()[0],
                                                       crossing_edge.dest().pos()[1])

                        #if so, evaluate the position of the intersection point between the crossing edge and the direction line
                        #if that point lies too close to an existing vertex (distance epsilon) then that vertex sould be used instead
                        if intersection != False:

                            if (abs(crossing_edge.dest().pos()[0] - intersection[0]) < epsilon) and (abs(crossing_edge.dest().pos()[1] - intersection[1]) < epsilon):
                                if not crossing_edge.next().is_obstacle():
                                    start_vertex.set_edge(crossing_edge.next().opp())
                                slice['mesh_vertices'].append(crossing_edge.next().opp().dest())
                                slice['edges'].append(crossing_edge.next().opp())
                                if obstacle:
                                    crossing_edge.next().opp().set_obstacle(True)
                                    # crossing_edge.next().set_obstacle(True)
                                return self.edge_slice(crossing_edge.dest(), end_vertex, None, epsilon, obstacle, slice)

                            elif (abs(crossing_edge.opp().dest().pos()[0] - intersection[0]) < epsilon) and (abs(crossing_edge.opp().dest().pos()[1] - intersection[1]) < epsilon):
                                if not crossing_edge.next().next().opp().is_obstacle():
                                    start_vertex.set_edge(crossing_edge.next().next())
                                slice['mesh_vertices'].append(crossing_edge.next().next().dest())
                                slice['edges'].append(crossing_edge.next().next())
                                if obstacle:
                                    crossing_edge.next().next().set_obstacle(True)
                                    # crossing_edge.next().next().opp().set_obstacle(True)
                                return self.edge_slice(crossing_edge.opp().dest(), end_vertex, None, epsilon, obstacle, slice)

                            #if a vertex needs to be inserted in the crossing edge then the crossed face will be split in two
                            #so must be the next face, opposing the crossing edge.
                            #since the next face will be split along the start/end line at the next round, when the next crossing vertex is found
                            #it will now be split simply by edges between the crossing vertex and the next face summit vertex
                            if type(crossing_edge.opp().get_face()) is MeshFace:
                                #TODO:  add a vetex at crossing point, split both triangles on the side of the crossing edge
                                #by adding edges and faces linking opping vertices, then iterate from the crossing point

                                #add the crossing vertex from intersection position
                                crossing_vertex = MeshVertex(intersection[0], intersection[1])
                                # added_vertices.append(crossing_vertex)

                                #add edges to link the starting vertex and the crossing vertexx
                                inward_new_edge = MeshEdge(crossing_edge.next().dest(), face= face, next_edge=crossing_edge.next().next())
                                face.set_edge(inward_new_edge)
                                inward_opp_edge = MeshEdge(crossing_vertex, opp_edge=inward_new_edge)
                                crossing_edge.next().dest().set_edge(inward_opp_edge)
                                #create a new face to split the current one at the newly created edges
                                inward_new_face = MeshFace(inward_opp_edge)
                                inward_opp_edge.set_face(inward_new_face)
                                crossing_edge.next().set_face(inward_new_face)
                                #close the start_vertex side of both faces around the new edges
                                crossing_edge.next().set_next(inward_opp_edge)
                                inward_new_edge.set_opp(inward_opp_edge)
                                #split the crossing edge at the crossing vertex, then assign each edge to its face
                                #(the the face that has be split for the split crossing edge and the newly created one for the new split edge)
                                split_new_edge = MeshEdge(crossing_edge.dest(), inward_new_face, next_edge=crossing_edge.next())
                                #then close both face on that side as well
                                inward_opp_edge.set_next(split_new_edge)
                                crossing_edge.set_next(inward_new_edge)
                                crossing_edge.set_dest(crossing_vertex)

                                #for the face on the opposite edge of the crossing edge to be split,
                                #the split segment is the one between the crossing vertex and the summit vertex of that face
                                #that edge and its opposite are created
                                outward_new_edge = MeshEdge(crossing_edge.opp().next().dest(), crossing_edge.opp().get_face(), next_edge=crossing_edge.opp().next().next())
                                crossing_edge.opp().get_face().set_edge(outward_new_edge)
                                outward_opp_edge = MeshEdge(crossing_vertex, opp_edge=outward_new_edge)
                                outward_new_edge.set_opp(outward_opp_edge)
                                #as the previous face, this one will aslo be split in two, so a new face is created
                                outward_new_face = MeshFace(outward_opp_edge)
                                outward_opp_edge.set_face(outward_new_face)
                                #the second new face is close at its summit end
                                crossing_edge.opp().next().set_next(outward_opp_edge)
                                crossing_edge.opp().next().set_face(outward_new_face)
                                #the the oppoing edge of the crossing edge is also split in two
                                split_opp_edge = MeshEdge(crossing_edge.opp().dest(), outward_new_face, crossing_edge, crossing_edge.opp().next())
                                #the second createdface is closed on its opposed crossing edge side
                                outward_opp_edge.set_next(split_opp_edge)
                                #the two split opposing edges need to be set as opposites of the two split crossig edges
                                crossing_edge.opp().set_opp(split_new_edge)
                                split_new_edge.set_opp(crossing_edge.opp())
                                #the opposing face now split can be closed on the opposed crossing side, completing te process
                                crossing_edge.opp().set_dest(crossing_vertex)
                                crossing_edge.opp().set_next(outward_new_edge)
                                crossing_edge.set_opp(split_opp_edge)

                                #TODO might not stop (check stop condition) about the last jump => no edge crossing, only a link to the opposite vertex

                                crossing_vertex.set_edge(outward_new_edge)
                                start_vertex.set_edge(inward_opp_edge)

                                if obstacle:
                                    inward_opp_edge.set_obstacle(True)

                                #the disposition of the edges once they get split may produce fces that don't
                                #comply to the delaunay criterion, so we check them along their creation
                                #and proceed to edge flips if necessary
                                if type(split_new_edge.next().opp().get_face()) is MeshFace:
                                    if clockwise_check(start_vertex, crossing_vertex, split_new_edge.next().opp().next().dest()):
                                        if not circumcircle_check(start_vertex, crossing_vertex, split_new_edge.dest(), split_new_edge.next().opp().next().dest()) \
                                                and not split_new_edge.next().is_obstacle() and not split_new_edge.next().opp().is_obstacle():

                                            edge_flip(split_new_edge.next())

                                    if (not circumcircle_check(split_new_edge.dest(), split_new_edge.next().dest(), split_new_edge.next().next().dest(), split_new_edge.opp().next().dest()) \
                                                or not circumcircle_check(split_new_edge.opp().dest(),split_new_edge.opp().next().dest(), split_new_edge.opp().next().next().dest(), split_new_edge.next().dest())) \
                                            and not split_new_edge.is_obstacle() and not split_new_edge.opp().is_obstacle():

                                        edge_flip(split_new_edge)

                                #add the new vertices and faces to the mesh
                                vertex_index = 0
                                index_found = False
                                while vertex_index <= len(self._vertices) and not index_found:
                                    if crossing_vertex.pos()[0] < self._vertices[vertex_index].pos()[0]:
                                        vertex_index += 1
                                    else:
                                        if crossing_vertex.pos()[1] < self._vertices[vertex_index].pos()[1]:
                                            vertex_index += 1
                                        else:
                                            index_found = True

                                self._vertices.insert(vertex_index, crossing_vertex)

                                self._faces.insert(0, inward_new_face)
                                self._faces.insert(0, outward_new_face)

                                # TODO: add the crossing vertex to an object representation in the map

                                #the function calls itself with the new crossing edge as starting point
                                #if the next vertex to be identified on the path is not the end vertex,
                                #then there should be a face around the crossing vertex to split in the direction of the end vertex

                                # ar = display_faces(self._faces, 'r')
                                # tx = []
                                # ty = []
                                # tx.append(start_vertex.pos()[0])
                                # tx.append(end_vertex.pos()[0])
                                # ty.append(start_vertex.pos()[1])
                                # ty.append(end_vertex.pos()[1])
                                # # ar.plot(tx,ty,'g')
                                # ex = []
                                # ey = []
                                # ex.append(crossing_edge.opp().dest().pos()[0])
                                # ex.append(crossing_edge.dest().pos()[0])
                                # ey.append(crossing_edge.opp().dest().pos()[1])
                                # ey.append(crossing_edge.dest().pos()[1])
                                # ar.plot(ex,ey,'k')
                                # pyplot.show()
                                slice['added_vertices'].append(crossing_vertex)
                                slice['edges'].append(inward_opp_edge)

                                return self.edge_slice(crossing_vertex, end_vertex, None, epsilon, obstacle, slice)
                            else:
                                return slice
                        else:
                            if on_segment_check(start_vertex.pos(), end_vertex.pos(), crossing_edge.next().opp().dest().pos()):
                                if not crossing_edge.next().is_obstacle():
                                    start_vertex.set_edge(crossing_edge.next().opp())
                                slice['mesh_vertices'].append(crossing_edge.next().opp().dest())
                                slice['edges'].append(crossing_edge.next().opp())
                                if obstacle:
                                    crossing_edge.next().opp().set_obstacle(True)
                                return self.edge_slice(crossing_edge.dest(), end_vertex, None, epsilon, obstacle, slice)
                            elif on_segment_check(start_vertex.pos(), end_vertex.pos(), crossing_edge.next().next().dest().pos()):
                                if not crossing_edge.next().next().opp().is_obstacle():
                                    start_vertex.set_edge(crossing_edge.next().next())
                                slice['mesh_vertices'].append(crossing_edge.next().next().dest())
                                slice['edges'].append(crossing_edge.next().next())
                                if obstacle:
                                    crossing_edge.next().next().set_obstacle(True)
                                return self.edge_slice(crossing_edge.next().next().dest(), end_vertex, None, epsilon, obstacle, slice)
                            else:
                                return slice

                    else:
                        #if the direction of the end vertex is not contained in the cone of this face,
                        #the the function is called with the next face in a rotation around the crossing vertex
                        return self.edge_slice(start_vertex, end_vertex, crossing_edge.next().opp().get_face(), epsilon, obstacle, slice)

    def add_solid(self, hull_points, epsilon = 1):
        #adds the reprensentation of a solid to the mesh
        #it takes a list of points that should define a non-null-area hull with its points given in clockwise order
        #it will then slice the mesh along the segments between each points
        #this will create new edges with the obstacle property on the inner edge of the hull
        #new vertices, faces will be added wherever an edge of the mesh is crossed by those new edges
        #the new vertices are returned to keep track of what represents the solid in the mesh (care for vertex superposition)

        init_hull_vertices = []
        own_vertices = []
        mesh_vertices = []
        hull_edges = []
        shared = []
        for hull_point in hull_points:
            hull_vertex = MeshVertex(hull_point[0], hull_point[1])

            #for each segment to become an edge, each x,y coordinate must be associated
            # with the face it should lie in in the mesh graph
            if self._grid !=[] and self._grid_len is not None:
                insert_face = self.find_face(hull_point, self._grid[int(hull_point[0]//self._grid_len)][int(hull_point[1]//self._grid_len)])
            else:
                insert_face = self.find_face(hull_point, self._vertices[0])

            #once each face is found from the point cordinate, it can be inserted in the mesh, splitting that face in three
            superposed_vertex = self.vertex_insertion(hull_vertex, insert_face, epsilon)

            if superposed_vertex is None:
                init_hull_vertices.append(hull_vertex)
                own_vertices.append(hull_vertex)
                shared_vertex = False
                sharing_solids = []
                for solid in self._solids:
                    if hull_vertex in solid.get_own_vertices():
                        shared_vertex = True
                        sharing_solids.append(solid)
                if shared_vertex:
                    shared.append({'vertex':hull_vertex,'solids':sharing_solids})
            else:
                init_hull_vertices.append(superposed_vertex)
                shared_vertex = False
                sharing_solids = []
                for solid in self._solids:
                    if superposed_vertex in solid.get_own_vertices():
                        shared_vertex = True
                        sharing_solids.append(solid)
                if shared_vertex:
                    own_vertices.append(superposed_vertex)
                    shared.append({'vertex':superposed_vertex,'solids':sharing_solids})

                else:
                    mesh_vertices.append(superposed_vertex)


        for hull_index in list(range(len(init_hull_vertices))):
            #for every pair of newly added vertices, a slice of the mesh between those vertices is done
            #it will create all edges and vertices on a straight line between them and give those edges the obstacle property
            next_index = hull_index+1
            if hull_index == len(init_hull_vertices) -1:
                next_index = 0
            hull = self.edge_slice(init_hull_vertices[hull_index], init_hull_vertices[next_index], epsilon=epsilon, slice={'added_vertices':[],'mesh_vertices':[],'edges':[]})
            hull_edges.extend(hull['edges'])
            own_vertices.extend(hull['added_vertices'])



            for hull_vertex in hull['mesh_vertices']:
                shared_vertex = False
                sharing_solids = []
                for solid in self._solids:

                    if hull_vertex in solid.get_own_vertices():
                        shared_vertex = True
                        sharing_solids.append(solid)
                if shared_vertex:
                    own_vertices.append(hull_vertex)
                    shared.append({'vertex':hull_vertex,'solids':sharing_solids})
                else:
                    mesh_vertices.append(hull_vertex)

                #TODO: sort the added vertices by x then y pos to keep coherence

        new_solid = MeshSolid(own_vertices, mesh_vertices, hull_edges)

        if shared != []:
            new_solid.set_shared_vertices(shared)
            for share_index in list(range(len(shared))):
                for sharing_solid in shared[share_index]['solids']:
                    sharing_solid.add_sharing_solid(shared[share_index]['vertex'], new_solid)

        self._solids.append(new_solid)

    def d_star(self, start_pos, end_pos, iter =0, avoid_vertices=[], epsilon=1):

        start_vertex = MeshVertex(start_pos[0], start_pos[1])
        end_vertex = MeshVertex(end_pos[0], end_pos[1])
        if self._grid != [] and self._grid is not None:
            start_face = self.find_face(start_pos, self._grid[int(start_pos[0]//self._grid_len)][int(start_pos[1]//self._grid_len)])
        else:
            start_face = self.find_face(start_pos, self._vertices[0])

        if self._grid != [] and self._grid_len is not None:
            end_face = self.find_face(end_pos, self._grid[int(start_pos[0]//self._grid_len)][int(start_pos[1]//self._grid_len)])
        else:
            end_face = self.find_face(end_pos, self._vertices[0])

        #TODO: check if start or end vertex is too close to an existing vetex (if so, pick that one instead but without deleting it at the end)

        start_superposed = self.vertex_insertion(start_vertex, start_face)
        if start_superposed is not None:
            start_vertex = start_superposed
        end_superposed = self.vertex_insertion(end_vertex, end_face)
        if end_superposed is not None:
            end_vertex = end_superposed

        upward_explored = [start_vertex]
        upward_explored.extend(avoid_vertices)
        upward_tree = {'vertices':[start_vertex], 'edges':[], 'cost':[], 'split':None}
        upward_connection = {'vertices':[start_vertex], 'edges':[], 'cost':[], 'split':None}
        downward_explored = [end_vertex]
        downward_explored.extend(avoid_vertices)
        downward_tree = {'vertices':[end_vertex], 'edges':[], 'cost':[], 'split': None}
        downward_connection = {'vertices':[end_vertex], 'edges':[], 'cost':[], 'split':None}
        upward_paths = []
        downward_paths = []

        connection = False
        while not connection or iter > 0:
            if upward_tree is not None:
                if downward_tree is not None:
                    upward_pathfinding = tree_crawler(end_vertex, list(upward_explored), dict(upward_tree), dict(upward_connection), dict(downward_connection), dir='up')
                else:
                    upward_pathfinding = tree_crawler(end_vertex, list(upward_explored), dict(upward_tree), dict(upward_connection), {'vertices':[end_vertex], 'edges':[], 'cost':[], 'split': None}, dir='up')
                upward_explored = upward_pathfinding['explored']
                upward_tree = upward_pathfinding['tree']
                upward_connection = upward_pathfinding['full_tree']
                downward_tree = upward_pathfinding['connection']

                if upward_pathfinding['paths'] != []:
                    upward_paths.extend(upward_pathfinding['paths'])
                    if not connection:
                        connection = True
                if connection:
                    iter -= 1

                # if upward_tree is not None:
                #     ar = display_tree(upward_tree, display_map(self))
                #     display_tree(downward_connection, ar, 'b')
                #     pyplot.show()
                #
                # if downward_tree is not None and len(downward_tree['vertices']) > 4:
                #     ar = display_tree(downward_tree, display_map(self))
                #     display_tree(upward_connection, ar, 'b')
                #     pyplot.show()

            else:
                if downward_tree is None:
                    iter=0

            if downward_tree is not None:
                if upward_tree is not None:
                    downward_pathfinding = tree_crawler(start_vertex, list(downward_explored), dict(downward_tree), dict(downward_connection), dict(upward_connection), dir='down')
                else:
                    downward_pathfinding = tree_crawler(start_vertex, list(downward_explored), dict(downward_tree), dict(downward_connection), {'vertices':[start_vertex], 'edges':[], 'cost':[], 'split':None}, dir='down')
                downward_explored = downward_pathfinding['explored']
                downward_tree = downward_pathfinding['tree']
                downward_connection = downward_pathfinding['full_tree']
                upward_tree = downward_pathfinding['connection']


                if downward_pathfinding['paths'] != []:
                    downward_paths.extend(downward_pathfinding['paths'])
                    if not connection:
                        connection = True
                if connection:
                    iter -= 1

                # if upward_tree is not None:
                #     ar = display_tree(upward_tree, display_map(self))
                #     display_tree(downward_connection, ar, 'b')
                #     pyplot.show()
                #
                # if downward_tree is not None and len(downward_tree['vertices']) > 4:
                #     ar = display_tree(downward_tree, display_map(self))
                #     display_tree(upward_connection, ar, 'b')
                #     pyplot.show()
            else:
                if upward_tree is None:
                    iter=0


        # print('up vertices')
        # for path in upward_paths:
        #     display_tree(path, display_map(self))
        #
        # print('up edges')
        # for path in upward_paths:
        #     print('***********')
        #     for nik in path['edges']:
        #         print(nik.dest().pos())
        #     print('**********')
        #     # segz = []
        #     # wegz = []
        #     # for f in self._faces:
        #     #     if type(f) is MeshFace:
        #     #         e = f.edge()
        #     #         for i in [0,1,2,3]:
        #     #             if e.is_obstacle():
        #     #                 wegz.append([e.opp().dest().pos(), e.dest().pos()])
        #     #             else:
        #     #                 segz.append([e.opp().dest().pos(), e.dest().pos()])
        #     #
        #     #             e = e.next()
        #     # fig = pyplot.figure()
        #     # ar = fig.add_subplot(1,1,1)
        #     # sar = array(segz)
        #     # for h in sar:
        #     #     tx = h[:,0]
        #     #     ty = h[:,1]
        #     #     ar.plot(tx,ty, color='r')
        #     # war = array(wegz)
        #     # for w in war:
        #     #     wx = w[:,0]
        #     #     wy = w[:,1]
        #     #     ar.plot(wx,wy, color='k')
        #     # pegz = []
        #     # for vindex in list(range(len(path)-1)):
        #     #     pegz.append([path[vindex].dest().pos(), path[vindex+1].dest().pos()])
        #     # par = array(pegz)
        #     # for g in par:
        #     #     rx = g[:,0]
        #     #     ry = g[:,1]
        #     #     ar.plot(rx,ry, color='g')
        #         # pyplot.show()
        #
        # print('down vertices')
        # for path in downward_paths:
        #     display_tree(path, display_map(self))
        #
        # print('down edges')
        # for path in downward_paths:
        #     print('***********')
        #     for nik in path['edges']:
        #         print(nik.dest().pos())
        #     print('**********')
        # #     segz = []
        # #     wegz = []
        # #     for f in self._faces:
        # #         if type(f) is MeshFace:
        # #             e = f.edge()
        # #             for i in [0,1,2,3]:
        # #                 if e.is_obstacle():
        # #                     wegz.append([e.opp().dest().pos(), e.dest().pos()])
        # #                 else:
        # #                     segz.append([e.opp().dest().pos(), e.dest().pos()])
        # #
        # #                 e = e.next()
        # #     fig = pyplot.figure()
        # #     ar = fig.add_subplot(1,1,1)
        # #     sar = array(segz)
        # #     for h in sar:
        # #         tx = h[:,0]
        # #         ty = h[:,1]
        # #         ar.plot(tx,ty, color='r')
        # #     war = array(wegz)
        # #     for w in war:
        # #         wx = w[:,0]
        # #         wy = w[:,1]
        # #         ar.plot(wx,wy, color='k')
        # #     pegz = []
        # #     for vindex in list(range(len(path)-1)):
        # #         pegz.append([path[vindex].dest().pos(), path[vindex+1].dest().pos()])
        # #     par = array(pegz)
        # #     for g in par:
        # #         rx = g[:,0]
        # #         ry = g[:,1]
        # #         ar.plot(rx,ry, color='g')
        #         # pyplot.show()
        print(len(upward_paths), len(downward_paths))
        return [upward_paths, downward_paths]

    def avoid_verices(self, center_pos, radius):
        vertices_to_avoid = []
        for vertex in self._vertices:
            if type(vertex) is MeshVertex:
                if pow(vertex.pos()[0] - center_pos[0],2) + pow(vertex.pos()[1] - center_pos[1],2) <= pow(radius,2):
                    vertices_to_avoid.append(vertex)
        return vertices_to_avoid

def merge_meshes(meshL, meshR):
    #Main function of the mesh generation that takes two valid delaunay meshes/convex hulls
    #then adds triangles untill both are merged in a single valid delauna mesh/hull

    vertices = []
    for v in meshL.get_vertices():
        if type(v) is MeshVertex:
            vertices.append(v)
    for v in meshR.get_vertices():
        if type(v) is MeshVertex:
            vertices.append(v)
    for v in meshL.get_vertices():
        if type(v) is GhostVertex:
            vertices.append(v)
    for v in meshR.get_vertices():
        if type(v) is GhostVertex:
            vertices.append(v)

    faces = []
    for f in meshL.get_faces():
        if type(f) is MeshFace:
            faces.append(f)
    for f in meshR.get_faces():
        if type(f) is MeshFace:
            faces.append(f)
    for f in meshL.get_faces():
        if type(f) is GhostFace:
            faces.append(f)
    for f in meshR.get_faces():
        if type(f) is GhostFace:
            faces.append(f)

    L_flat = len([f for f in meshL.get_faces() if type(f) is  MeshFace]) == 0
    R_flat = len([f for f in meshR.get_faces() if type(f) is  MeshFace]) == 0

    #Find the upper and lower boundary vertices of each convex hull
    #TODO: check for occlusion between each upper or lower vertices (if there is, select another starting or exit vertex!)
    upper_L_vertex = meshL.get_vertices()[0]
    upper_L_height = meshL.get_vertices()[0].pos()[1]
    for vertex in meshL.get_vertices():
        if type(vertex) is MeshVertex:
            if vertex.pos()[1] > upper_L_height:
                upper_L_vertex = vertex
                upper_L_height = vertex.pos()[1]
            elif vertex.pos()[1] == upper_L_height and vertex.pos()[0] < upper_L_vertex.pos()[0]:
                upper_L_vertex = vertex

    lower_L_vertex = meshL.get_vertices()[0]
    lower_L_height = meshL.get_vertices()[0].pos()[1]
    for vertex in meshL.get_vertices():
        if type(vertex) is MeshVertex:
            if vertex.pos()[1] < lower_L_height:
                lower_L_vertex = vertex
                lower_L_height = vertex.pos()[1]
            elif vertex.pos()[1] == lower_L_height and vertex.pos()[0] < lower_L_vertex.pos()[0]:
                lower_L_vertex = vertex

    upper_R_vertex = meshR.get_vertices()[0]
    upper_R_height = meshR.get_vertices()[0].pos()[1]
    for vertex in meshR.get_vertices():
        if type(vertex) is MeshVertex:
            if vertex.pos()[1] > upper_R_height:
                upper_R_vertex = vertex
                upper_R_height = vertex.pos()[1]
            elif vertex.pos()[1] == upper_R_height and vertex.pos()[0] > upper_R_vertex.pos()[0]:
                upper_R_vertex = vertex

    lower_R_vertex = meshR.get_vertices()[0]
    lower_R_height = meshR.get_vertices()[0].pos()[1]
    for vertex in meshR.get_vertices():
        if type(vertex) is MeshVertex:
            if vertex.pos()[1] < lower_R_height:
                lower_R_vertex = vertex
                lower_R_height = vertex.pos()[1]
            elif vertex.pos()[1] == lower_R_height and vertex.pos()[0] > lower_R_vertex.pos()[0]:
                lower_R_vertex = vertex

    #case of to flat and vertical edges that are colinear, then no face is added; only a linking edge and its ghost faces
    if L_flat and R_flat and lower_L_vertex.pos()[0] == upper_L_vertex.pos()[0] \
            and lower_R_vertex.pos()[0] == upper_R_vertex.pos()[0] \
            and upper_L_vertex.pos()[0] == lower_R_vertex.pos()[0]:

        upper_L_edge = upper_L_vertex.edge()
        lower_R_edge = lower_R_vertex.edge()
        ghosts = create_isolated_edge(upper_L_vertex, lower_R_vertex)

        #after the edge and opposite are created, along with thir ghost faces, they are linked to the next ghosts
        #thatgo around the colinear edges

        upper_L_vertex.edge().next().set_opp(lower_R_edge.next().next())
        lower_R_edge.next().next().set_opp(upper_L_vertex.edge().next())
        lower_R_vertex.edge().next().next().set_opp(lower_R_edge.opp().next())
        lower_R_edge.opp().next().set_opp(lower_R_vertex.edge().next().next())

        upper_L_vertex.edge().next().next().set_opp(upper_L_edge.opp().next())
        upper_L_edge.opp().next().set_opp(upper_L_vertex.edge().next().next())
        lower_R_vertex.edge().next().set_opp(upper_L_edge.next().next())
        upper_L_edge.next().next().set_opp(lower_R_vertex.edge().next())

        lower_R_vertex.set_edge(lower_R_edge)

        #fill up the mesh with the new vertices and faces
        vertices.append(ghosts['vertices'][0])
        vertices.append(ghosts['vertices'][1])

        faces.append(ghosts['faces'][0])
        faces.append(ghosts['faces'][1])

        return {'vertices': vertices, 'faces':faces}
    else:
        #otherwise, at least one side of the merge has an area, or not vertically colinear,
        #meaning that new triangles/faces will have to be created
        start_R = lower_R_vertex

        start_L = lower_L_vertex
        stop_R = upper_R_vertex
        stop_L = upper_L_vertex
        #start by setting the start and stop vertices to the upper and lower bound vertices of each mesh/hull

        #check the case of a flat mesh that needs to be merged with a non flat mesh.
        # If the starting point on the flat mesh is on the same vertical line as a point on the other mesh,
        #then that oint on the other mesh should be the starting point (or ending point depending if left or right)
        if R_flat and not L_flat:
            aligned_L = None
            for vertex in meshL.get_vertices():
                if type(vertex) is MeshVertex:
                    if vertex.pos()[0] == lower_R_vertex.pos()[0]:
                        if aligned_L is None:
                            aligned_L = vertex
                        elif vertex.pos()[1] > aligned_L.pos()[1]:
                            aligned_L = vertex
            if aligned_L is not None:
                start_L = aligned_L

        if L_flat and not R_flat:
            aligned_R = None
            for vertex in meshR.get_vertices():
                if type(vertex) is MeshVertex:
                    if vertex.pos()[0] == upper_L_vertex.pos()[0]:
                        if aligned_R is None:
                            aligned_R = vertex
                        elif vertex.pos()[1] < aligned_R.pos()[1]:
                            aligned_R = vertex
            if aligned_R is not None:
                stop_R = aligned_R

        #then check for occlusion or extension of the field of view
        #an occlusion occurs if the segment between the end or start vertices crosses any mesh
        #some half plane checks are performed.
        #One way is the occlusion, then the start or end vertexmust be displaced
        # towards the center to narrow the field of view
        #The other way is to extnd the field of view if there is no occusion
        # while picking the ext vertex as start or end
        no_occlusion = False
        #check and resolution of the occlusion
        while not no_occlusion:

            no_occlusion = True
            if not R_flat:
                if start_L.pos()[1] >= start_R.pos()[1]:
                    if half_plane_check(start_R.edge().dest().pos(),start_R.pos(), start_L.pos()) \
                            and half_plane_check(start_R.pos(), start_R.edge().opp().next().opp().next().dest().pos(), start_L.pos()):
                        no_occlusion = False
                        start_R = start_R.edge().dest()

                if stop_L.pos()[1] <= stop_R.pos()[1]:
                    if half_plane_check(stop_R.edge().dest().pos(), stop_R.pos(),stop_L.pos()) \
                            and half_plane_check(stop_R.pos(), stop_R.edge().opp().next().opp().next().dest().pos(),stop_L.pos()):
                        no_occlusion = False
                        stop_R = stop_R.edge().opp().next().opp().next().dest()
            if not L_flat:
                if start_L.pos()[1] <= start_R.pos()[1]:
                    if half_plane_check(start_L.edge().dest().pos(), start_L.pos(), start_R.pos()) \
                            and half_plane_check(start_L.pos(), start_L.edge().opp().next().opp().next().dest().pos(),
                                                 start_R.pos()):
                        no_occlusion = False
                        start_L = start_L.edge().opp().next().opp().next().dest()
                if stop_L.pos()[1] >= stop_R.pos()[1]:
                    if half_plane_check(stop_L.pos(), stop_L.edge().opp().next().opp().next().dest().pos(), stop_R.pos()) \
                            and half_plane_check(stop_L.edge().dest().pos(), stop_L.pos(), stop_R.pos()):
                        no_occlusion = False
                        stop_L = stop_L.edge().dest()

        while no_occlusion:
            #check and resolution of the extenion of the field of view, just before any occlusion happens

            no_occlusion = False
            if not R_flat:
                if start_L.pos()[1] < start_R.pos()[1]:
                    if not half_plane_check(start_R.pos(), start_R.edge().opp().next().opp().next().dest().pos(), start_L.pos()):
                        no_occlusion = True
                        start_R = start_R.edge().opp().next().opp().next().dest()
                if stop_L.pos()[1] > stop_R.pos()[1]:
                    if not half_plane_check(stop_R.edge().dest().pos(), stop_R.pos(),stop_L.pos()):
                        no_occlusion = True
                        stop_R = stop_R.edge().dest()
            if not L_flat:
                if start_L.pos()[1] > start_R.pos()[1]:
                    if not half_plane_check(start_L.edge().dest().pos(), start_L.pos(), start_R.pos()):
                        no_occlusion = True
                        start_L = start_L.edge().dest()
                if stop_L.pos()[1] < stop_R.pos()[1]:
                    if not half_plane_check(stop_L.pos(), stop_L.edge().opp().next().opp().next().dest().pos(), stop_R.pos()):
                        no_occlusion = True
                        stop_L = stop_L.edge().opp().next().opp().next().dest()

        #Once the start and end vertices are selected, the merge process cn begin
        #create an edge between both lower points and its opposite
        merge_edge_start = MeshEdge(start_L)
        merge_edge_opp = MeshEdge(start_R, opp_edge=merge_edge_start)

        merge_edge_start.set_opp(merge_edge_opp)

        #add the ghost faces under and over them
        over_ghost_face = add_ghost_vertex(merge_edge_start)
        under_ghost_face = add_ghost_vertex(merge_edge_opp)

        faces.append(over_ghost_face['face'])
        faces.append(under_ghost_face['face'])
        vertices.append(over_ghost_face['vertex'])
        vertices.append(under_ghost_face['vertex'])

        #link the opposites edges of the ghost faces bellow and over the starting edge
        if start_L != stop_L:
            merge_edge_start.next().set_opp(start_L.edge().opp().next().opp())
            start_L.edge().opp().next().opp().set_opp(merge_edge_start.next())
            merge_edge_opp.next().next().set_opp(start_L.edge().opp().next())
            start_L.edge().opp().next().set_opp(merge_edge_opp.next().next())
        else:
            merge_edge_start.next().set_opp(merge_edge_start.opp().next().next())
            merge_edge_start.opp().next().next().set_opp(merge_edge_start.next())

        if start_R != stop_R:
            merge_edge_opp.next().set_opp(start_R.edge().opp().next().opp())
            start_R.edge().opp().next().opp().set_opp(merge_edge_opp.next())
            merge_edge_start.next().next().set_opp(start_R.edge().opp().next())
            start_R.edge().opp().next().set_opp(merge_edge_start.next().next())
        else:
            merge_edge_start.next().next().set_opp(merge_edge_start.opp().next())
            merge_edge_start.opp().next().set_opp(merge_edge_start.next().next())

        current_edge = merge_edge_start

        #after the starting edge is initialized, new edges and faces will be created until the end vertices get linked
        while not (current_edge.dest() == stop_L and current_edge.opp().dest() == stop_R):
            merged = merge_edges(current_edge, faces, vertices)

            faces = merged['faces']
            vertices = merged['vertices']

            #checks the delaunay property and performs edge flips if necessary
            if type(current_edge.opp().get_face()) is MeshFace:
                if not circumcircle_check(current_edge.dest(), current_edge.next().dest(), current_edge.opp().dest(), current_edge.opp().next().dest()):
                    if not colinear_check(current_edge.opp().next().dest(), current_edge.dest(), current_edge.next().dest()) \
                            and not colinear_check(current_edge.opp().next().dest(), current_edge.opp().dest(), current_edge.next().dest()):
                        edge_flip(current_edge)

            if type(current_edge.next().opp().get_face()) is MeshFace:
                if not circumcircle_check(current_edge.dest(), current_edge.next().dest(), current_edge.opp().dest(), current_edge.next().opp().next().dest()):
                    if not colinear_check(current_edge.next().opp().next().dest(), current_edge.dest(), current_edge.opp().dest()) \
                            and not colinear_check(current_edge.next().opp().next().dest(), current_edge.next().dest(), current_edge.opp().dest()):
                        edge_flip(current_edge.next())

            if type(current_edge.next().next().opp().get_face()) is MeshFace:
                if not circumcircle_check(current_edge.dest(), current_edge.next().dest(), current_edge.opp().dest(), current_edge.next().next().opp().next().dest()):
                    if not colinear_check(current_edge.next().next().opp().next().dest(), current_edge.dest(), current_edge.next().dest()) \
                            and not colinear_check(current_edge.next().next().opp().next().dest(), current_edge.dest(), current_edge.opp().dest()):
                        edge_flip(current_edge.next().next())

            current_edge = merged['edge'].opp()

        current_edge.dest().set_edge(current_edge.opp())
        start_L.set_edge(merge_edge_start.opp().next().next().opp().next().next().opp())
        start_R.set_edge(merge_edge_start)

        return {'faces': faces, 'vertices':vertices}

def merge_edges(start_edge, faces, vertices):
    #######
    #One step of the recursive process that generates the mesh
    #Finds and evaluates suitting (clockwise order) candidates amongst the starting edge neighbours
    #elects the vertex that will maximise the minimal angle and creates a new face with it
    #merges the ghost edges of the starting ege and the elected edge where the face is created
    # and associates it to the newly created outter edge
    #checks for collision/crossing with pre-existing edges and proceeds to remove them

    #Identify vertices around the starting edge
    startL = start_edge.dest()
    startR = start_edge.opp().dest()

    #initiate the duo of vertices candidates, starting with the most exterior edge (minimal angle)
    #can be found from the ghost faces of the starting edge
    edgeL = start_edge.next().opp().next().opp()

    if edgeL.opp().dest() == startR:
        #Left hand side one point case; startL is the only left point, so it is also the left candidate
        #then startL == Lcandidate can identify this case for handling
        Lcandidate = startL
        Lfound = True
    else:
        #otherwise, Lcandidate is set to the outter neighbouring vertex
        # and the next edge is set by rotation around the startL vertex
        Lcandidate = edgeL.opp().dest()
        next_edgeL = edgeL.next().opp()
        Lfound = False

    #a counter that holds the number of crossed neighbouring edges that will need to be removed
    Ledges_crossed = 0


    while clockwise_check(startR, startL, Lcandidate) and Lfound==False:
        #the search for a better candidate should stop after the pi angle threshold is crossed
        # or a ghost face signifies the border of the mesh/hull
        if type(next_edgeL) is MeshEdge:
            next_Lcandidate = next_edgeL.opp().dest()
            #check if the candidate vertex complies to the delaunay criterion
            if not circumcircle_check(startR, startL, Lcandidate, next_Lcandidate):
                #if not, change the candidate couple with the next vertex connected to startL by rotation
                edgeL = next_edgeL
                Lcandidate = next_Lcandidate
                next_edgeL = edgeL.next().opp()
                Ledges_crossed += 1
            else:
                #if the delaunay criterion is valid, register a left candidate
                Lfound = True
        else:
            Lfound = True

    #Same sogic applies to the Right hand side
    #only differences comes from the clockwise roation of edges around a mesh/hull and inside faces themselves
    edgeR = start_edge.next().next().opp().next().next().opp()

    Rfound = False

    if edgeR.dest() == startL:
        Rcandidate = startR
        Rfound = True
    else:
        Rcandidate = edgeR.dest()
        next_edgeR = edgeR.next().next().opp()
        Rfound = False

    Redges_crossed = 0

    while clockwise_check(startR, startL, Rcandidate) and Rfound==False:
        if type(next_edgeR) is MeshEdge:
            next_Rcandidate = next_edgeR.dest()
            if not circumcircle_check(startR, startL, Rcandidate, next_Rcandidate):
                edgeR = next_edgeR
                Rcandidate = next_Rcandidate
                next_edgeR = edgeR.next().next().opp()
                Redges_crossed += 1
            else:
                Rfound = True
        else:
            Rfound = True

    #at the point there should be at least one candidate found either on the left or right hand side
    #if both sides have candidates, the better suitted one should be elected
    elected_edge = None
    elected_vertex = None
    if Rfound or Lfound:
        if Rfound and Lfound:
            #the choice of candidate vertex is done with the delaunay criterion using a circumcircle check
            if circumcircle_check(startR, startL, Lcandidate, Rcandidate):
                elected_vertex = Lcandidate
                elected_edge = edgeL
                Rfound = False
            else:
                elected_vertex = Rcandidate
                elected_edge = edgeR
                Lfound = False
        else:
            if Lfound:
                elected_vertex = Lcandidate
                elected_edge = edgeL
            if Rfound:
                elected_vertex = Rcandidate
                elected_edge = edgeR

    if elected_vertex is None:
        pyplot.show()
        print('aie')

    if elected_vertex is not None and elected_edge is not None:


        #finally, connect the elected vertex with edges to form the new face of the mesh,
        # delete (unreference) the edges that were crossed by the new edge,
        #and the new edge's opposite edge should belong to an external ghost edge
        # that is a merge of pre-existing external edges of the start and elected edges because
        #they should be already properly linked with the other external ghost faces around the mesh/hull
        if Lfound:
            if Ledges_crossed > 0:
                #loop that will erase/unreference edges that will be crossed by the new edge
                current_edge = start_edge.next().opp().next().opp()
                a,b = elected_vertex, startR
                for edge_num in list(range(Ledges_crossed)):
                    k,l = current_edge.dest(), current_edge.opp().dest()
                    #the next edge is defined by turning around the assiociated start vertex to find neighbouring vertices
                    next_edge = current_edge.next().opp()
                    #for every face not elected as candidate, the non elected edge should be shecked
                    #for collision with the new edge and if there is, removed
                    if intersect_check(a.pos()[0],a.pos()[1],
                                       b.pos()[0],b.pos()[1],
                                       k.pos()[0],k.pos()[1],
                                       l.pos()[0],l.pos()[1]):

                        #the ghost face of the edge and opposite about to be erased is assigned
                        #to the next neighbouring edge to be evaluated

                        #vertices or faces that could point to the edge to be destroyed are reassigned to point to external edges

                        current_edge.opp().next().next().set_next(current_edge.next().next())
                        current_edge.opp().dest().set_edge(current_edge.next().next().opp())

                        if current_edge.dest().edge() == current_edge.opp():
                            current_edge.dest().set_edge(current_edge.next())

                        if type(current_edge.next().next().opp().get_face()) is MeshFace or current_edge.opp().dest().edge() ==  current_edge:
                            current_edge.opp().dest().set_edge(current_edge.next().next().opp())

                        if type(current_edge.next().opp()) is MeshFace():
                            current_edge.next().dest().set_edge(current_edge.next().opp())
                        else:
                            current_edge.next().dest().set_edge(current_edge.next().next())
                        current_edge.next().next().set_next(current_edge.opp().next())
                        # current_edge.next().next().dest().set_edge(current_edge.next().next().opp())




                        #the mesh face is unreferenced/removed from any list composing the graph
                        current_edge.get_face().set_edge(None)
                        current_edge.next().next().set_face(current_edge.opp().get_face())
                        current_edge.next().set_next(None)
                        current_edge.next().set_face(None)
                        del faces[faces.index(current_edge.get_face())]
                        current_edge.set_face(None)

                        #add ghost vertex and face to the second resulting external edge efter the removal
                        #each new external edge inherists half the ghost edges of the ghost face
                        # that was on the opposite side of the edge to be removed because they hold links
                        #to the rest of the external mesh/hull
                        # both ghost faces will then be completed with opposing sides
                        new_ghost = add_ghost_vertex(current_edge.next())
                        faces.append(new_ghost['face'])
                        vertices.append(new_ghost['vertex'])
                        current_edge.next().next().set_opp(current_edge.opp().next().next())
                        current_edge.next().next().next().set_opp(start_edge.next())

                        #finally, the edge and it's opposite can be removed frm any list in the mesh, removing any reference
                        current_edge.opp().set_next(None)
                        current_edge.opp().set_face(None)
                        current_edge.opp().set_opp(None)
                        current_edge.opp().set_dest(None)


                        current_edge.set_opp(None)
                        current_edge.set_next(None)
                        current_edge.set_dest(None)
                        current_edge.set_face(None)
                        del current_edge


                    current_edge = next_edge

            #create the new edge and face
            new_edge = MeshEdge(startR, next_edge=start_edge)
            new_face = MeshFace(start_edge)

            #add an opposite edge
            new_opp_edge = MeshEdge(elected_vertex, opp_edge=new_edge)
            new_edge.set_opp(new_opp_edge)


            #unset ghost face reference to prepare the merging of the ghost faces
            elected_edge.opp().get_face().set_edge(None)
            del faces[faces.index(elected_edge.opp().get_face())]
            #set the face pointer of the previously exterior elected edge to the new face
            elected_edge.opp().set_face(new_face)
            #associates the external ghost edge to the ghost face of the starting edge, to merge both ghost faces
            elected_edge.opp().next().set_face(start_edge.get_face())
            #delete refernces of the second ghost edge
            elected_edge.opp().next().next().set_next(None)
            elected_edge.opp().next().next().set_opp(None)
            elected_edge.opp().next().next().set_face(None)
            elected_edge.opp().next().next().set_vertex(None)
            #set the next attribute to link both gohst triangles
            elected_edge.opp().next().set_next(start_edge.next().next())

            # if start_edge.next().vertex() in vertices:
            #delete a ghost infinity vertex
            del vertices[vertices.index(elected_edge.opp().next().vertex())]
            elected_edge.opp().next().set_vertex(start_edge.next().vertex())

            #close the new merged ghost triangle
            elected_edge.opp().next().next().set_next(new_opp_edge)
            new_opp_edge.set_next(elected_edge.opp().next())
            new_opp_edge.set_face(start_edge.get_face())
            # new_opp_edge.next().set_vertex(start_edge.next().vertex())
            new_opp_edge.get_face().set_edge(new_opp_edge)

            #delete all references to the first ghost edge of the starting face to complete the merge
            start_edge.next().set_face(None)
            start_edge.next().set_next(None)
            start_edge.next().set_opp(None)
            start_edge.next().set_vertex(None)

            #finally, close the the new face by linking edges clockwise
            start_edge.set_next(elected_edge.opp())
            elected_edge.opp().set_next(new_edge)

            #associate the new face
            start_edge.set_face(new_face)
            new_edge.set_face(new_face)

            #set the vertices so that they point along the edges of the created hull (clockwise)
            if type(elected_edge.get_face()) is MeshFace:
                elected_vertex.set_edge(new_edge)
            else:
                startL.set_edge(elected_edge.opp())

            #add the new face to the mesh
            faces.insert(0,new_face)

            merged = {'faces': faces, 'edge': new_edge, 'vertices':vertices}
            return merged

        if Rfound:

            if Redges_crossed > 0:
                current_edge = start_edge.next().next().opp().next().next().opp()
                a,b = elected_vertex, startL
                for edge_num in list(range(Redges_crossed)):
                    k,l = current_edge.dest(), current_edge.opp().dest()
                    next_edge = current_edge.next().next().opp()
                    if intersect_check(a.pos()[0],a.pos()[1],
                                       b.pos()[0],b.pos()[1],
                                       k.pos()[0],k.pos()[1],
                                       l.pos()[0],l.pos()[1]):

                        current_edge.next().next().set_next(None)
                        current_edge.next().next().set_face(None)

                        if current_edge.dest().edge() == current_edge.opp():
                            current_edge.dest().set_edge(current_edge.next())

                        if type(current_edge.next().opp().get_face()) is MeshFace:
                            current_edge.next().dest().set_edge(current_edge.next().opp())
                        else:
                            current_edge.next().dest().set_edge(current_edge.next().next())

                        if type(current_edge.next().next().opp().get_face()) is MeshFace or current_edge.next().next().dest().edge() == current_edge:

                            # startR.set_edge(current_edge.next().next().opp())
                            current_edge.next().next().dest().set_edge(current_edge.next().next().opp())
                        else:
                            current_edge.next().next().dest().set_edge(current_edge.next().next().opp().next().next().opp().next().next().opp())

                        new_ghost = add_ghost_vertex(current_edge.next().next())
                        faces.append(new_ghost['face'])
                        vertices.append(new_ghost['vertex'])
                        current_edge.next().next().next().next().set_opp(current_edge.opp().next())
                        current_edge.next().next().next().set_opp(start_edge.next().next())
                        current_edge.next().set_next(current_edge.opp().next())
                        current_edge.next().set_face(current_edge.opp().get_face())
                        current_edge.opp().next().next().set_next(current_edge.next())

                        current_edge.get_face().set_edge(None)
                        del faces[faces.index(current_edge.get_face())]

                        current_edge.opp().set_next(None)
                        current_edge.opp().set_opp(None)
                        current_edge.opp().set_face(None)
                        current_edge.opp().set_dest(None)

                        current_edge.set_face(None)
                        current_edge.set_opp(None)
                        current_edge.set_next(None)
                        current_edge.set_dest(None)
                        del current_edge

                    current_edge = next_edge

            #create new edge
            new_edge = MeshEdge(elected_vertex, next_edge=elected_edge.opp())
            new_face = MeshFace(start_edge)

            new_opp_edge = MeshEdge(startL, opp_edge=new_edge)
            new_edge.set_opp(new_opp_edge)

            #delete references to left ghost face
            start_edge.get_face().set_edge(None)
            del faces[faces.index(start_edge.get_face())]

            start_edge.next().set_face(elected_edge.opp().get_face())

            #delete references to left side of linking ghost face

            start_edge.next().next().set_face(None)
            start_edge.next().next().set_opp(None)
            start_edge.next().next().set_next(None)
            start_edge.next().next().set_vertex(None)

            #link both ghost by a next pointer
            start_edge.next().set_next(elected_edge.opp().next().next())
            del vertices[vertices.index(elected_edge.opp().next().vertex())]

            elected_edge.opp().next().set_vertex(None)

            #then delete right hand side of the ghost edge to remove
            elected_edge.opp().next().set_face(None)
            elected_edge.opp().next().set_opp(None)
            elected_edge.opp().next().set_next(None)
            elected_edge.opp().next().set_vertex(None)

            #redefine the last edge of the ghost triangle
            start_edge.next().next().set_next(new_opp_edge)
            new_opp_edge.set_next(start_edge.next())

            new_opp_edge.set_face(elected_edge.opp().get_face())
            new_opp_edge.get_face().set_edge(new_opp_edge)

            start_edge.set_face(new_face)
            start_edge.set_next(new_edge)
            elected_edge.opp().set_next(start_edge)
            elected_edge.opp().set_face(new_face)
            new_edge.set_face(new_face)

            if type(elected_edge.get_face()) is not MeshFace:
                elected_vertex.set_edge(elected_edge.opp())
            faces.insert(0,new_face)

            merged = {'faces': faces, 'edge': new_edge, 'vertices':vertices}
            return merged

def pathfinding_step(end_vertex, explored, vertices_path, edges_path=[], cost_path=[]):
    explored_vertices = explored
    vertices_path = vertices_path
    edges_path = edges_path
    cost_path = cost_path
    first_edge = vertices_path[-1].edge()
    start_edge = first_edge

    while start_edge.dest() in explored_vertices:
        if type(start_edge.get_face()) is MeshFace:
            if start_edge.is_obstacle():
                while not (start_edge.opp().is_obstacle() and not start_edge.is_obstacle()):
                    if type(start_edge.get_face()) is MeshFace:
                        start_edge = start_edge.next().next().opp()
                    else:
                        start_edge = start_edge.next().next().opp().next().next().opp()
                if start_edge == first_edge:
                    if len(vertices_path) > 1:
                        vertices_path.remove(vertices_path[-1])
                        if len(cost_path) > 0:
                            cost_path.remove(cost_path[-1])
                        if len(edges_path) >0:
                            edges_path.remove(edges_path[-1])
                    else:
                        return None

                    start_edge = vertices_path[-1].edge()
                    first_edge = start_edge
            else:
                start_edge = start_edge.next().next().opp()
                if start_edge == first_edge:
                    if len(vertices_path) > 1:
                        vertices_path.remove(vertices_path[-1])
                        if len(cost_path) > 0:
                            cost_path.remove(cost_path[-1])
                        if len(edges_path) >0:
                            edges_path.remove(edges_path[-1])
                    else:
                        return None
                    start_edge = vertices_path[-1].edge()
                    first_edge = start_edge
        else:
            start_edge = start_edge.next().next().opp().next().next().opp()
            if start_edge == first_edge:
                vertices_path.remove(vertices_path[-1])
                if len(cost_path) > 0:
                    cost_path.remove(cost_path[-1])
                if len(edges_path) >0:
                    edges_path.remove(edges_path[-1])
                else:
                    return None

                start_edge = vertices_path[-1].edge()
                first_edge = start_edge

    best_edge = start_edge
    split = False

    edge_vector = [best_edge.dest().pos()[0]-best_edge.opp().dest().pos()[0], best_edge.dest().pos()[1]-best_edge.opp().dest().pos()[1]]
    edge_dist = pow(edge_vector[0],2)+pow(edge_vector[1],2)
    direction = [end_vertex.pos()[0]-best_edge.opp().dest().pos()[0], end_vertex.pos()[1]-best_edge.opp().dest().pos()[1]]
    projection_scalar = (edge_vector[0]*direction[0]) + (edge_vector[1]*direction[1])
    best_ratio = projection_scalar/edge_dist

    if start_edge.dest() == end_vertex:
        vertices_path.append(best_edge.dest())
        explored_vertices.append(best_edge.dest())
        edges_path.append(best_edge)
        cost_path.append(best_ratio)
    else:
        if type(start_edge.get_face()) is MeshFace:
            if start_edge.is_obstacle():
                current_edge = start_edge.next().next().opp()
                while not (current_edge.opp().is_obstacle() and not current_edge.is_obstacle()):
                    if type(current_edge.get_face()) is MeshFace:
                        current_edge = current_edge.next().next().opp()
                    else:
                        current_edge = current_edge.next().next().opp().next().next().opp()
            else:
                current_edge = start_edge.next().next().opp()

        else:
            current_edge = start_edge.next().next().opp().next().next().opp()

        if current_edge.dest() == end_vertex:

            edge_vector = [current_edge.dest().pos()[0]-current_edge.opp().dest().pos()[0], current_edge.dest().pos()[1]-current_edge.opp().dest().pos()[1]]
            edge_dist = pow(edge_vector[0],2)+pow(edge_vector[1],2)
            direction = [end_vertex.pos()[0]-current_edge.opp().dest().pos()[0], end_vertex.pos()[1]-current_edge.opp().dest().pos()[1]]
            projection_scalar = (edge_vector[0]*direction[0]) + (edge_vector[1]*direction[1])
            current_ratio = projection_scalar/edge_dist

            vertices_path.append(current_edge.dest())
            explored_vertices.append(current_edge.dest())
            edges_path.append(current_edge)
            cost_path.append(current_ratio)
        else:
            while current_edge != start_edge and current_edge.dest() != end_vertex:
                if current_edge.dest() not in explored_vertices:
                    edge_vector = [current_edge.dest().pos()[0]-current_edge.opp().dest().pos()[0], current_edge.dest().pos()[1]-current_edge.opp().dest().pos()[1]]
                    edge_dist = pow(edge_vector[0],2)+pow(edge_vector[1],2)
                    direction = [end_vertex.pos()[0]-current_edge.opp().dest().pos()[0], end_vertex.pos()[1]-current_edge.opp().dest().pos()[1]]
                    projection_scalar = (edge_vector[0]*direction[0]) + (edge_vector[1]*direction[1])
                    current_ratio = projection_scalar/edge_dist

                    if current_ratio > best_ratio:
                        best_ratio = current_ratio
                        best_edge = current_edge

                if type(current_edge.get_face()) is MeshFace:
                    if current_edge.is_obstacle():
                        while not (current_edge.opp().is_obstacle() and not current_edge.is_obstacle()):
                            if type(current_edge.get_face()) is MeshFace:
                                current_edge = current_edge.next().next().opp()
                            else:
                                current_edge = current_edge.next().next().opp().next().next().opp()
                    else:
                        current_edge = current_edge.next().next().opp()
                else:
                    current_edge = current_edge.next().next().opp().next().next().opp()

            if best_edge.dest() == end_vertex:
                vertices_path.append(best_edge.dest())
                explored_vertices.append(best_edge.dest())
                edges_path.append(best_edge)
                cost_path.append(best_ratio)

            #TODO: chek this shortcut method's efficiency

            else:
                shortcut = True
                while shortcut:
                    shortcut = False
                    if len(vertices_path) > 2:
                        if type(best_edge.opp().get_face()) is MeshFace:
                            if vertices_path[-2] == best_edge.opp().next().dest() and not best_edge.opp().next().is_obstacle:
                                shortcut = True
                                best_edge = best_edge.opp().next().next()
                                vertices_path.remove(vertices_path[-1])
                                if len(edges_path) > 0:
                                    edges_path.remove(edges_path[-1])
                                if len(cost_path) > 0:
                                    cost_path.remove(cost_path[-1])
                                edge_vector = [current_edge.dest().pos()[0]-current_edge.opp().dest().pos()[0], current_edge.dest().pos()[1]-current_edge.opp().dest().pos()[1]]
                                edge_dist = pow(edge_vector[0],2)+pow(edge_vector[1],2)
                                direction = [end_vertex.pos()[0]-current_edge.opp().dest().pos()[0], end_vertex.pos()[1]-current_edge.opp().dest().pos()[1]]
                                projection_scalar = (edge_vector[0]*direction[0]) + (edge_vector[1]*direction[1])
                                best_ratio = projection_scalar/edge_dist

                            elif type(best_edge.get_face()) is MeshFace:
                                if vertices_path[-2] == best_edge.next().dest() and not best_edge.next().next().is_obstacle():
                                    shortcut = True
                                    best_edge = best_edge.next().opp()
                                    vertices_path.remove(vertices_path[-1])
                                    if len(edges_path) > 0:
                                        edges_path.remove(edges_path[-1])
                                    if len(cost_path) > 0:
                                        cost_path.remove(cost_path[-1])
                                    edge_vector = [current_edge.dest().pos()[0]-current_edge.opp().dest().pos()[0], current_edge.dest().pos()[1]-current_edge.opp().dest().pos()[1]]
                                    edge_dist = pow(edge_vector[0],2)+pow(edge_vector[1],2)
                                    direction = [end_vertex.pos()[0]-current_edge.opp().dest().pos()[0], end_vertex.pos()[1]-current_edge.opp().dest().pos()[1]]
                                    projection_scalar = (edge_vector[0]*direction[0]) + (edge_vector[1]*direction[1])
                                    best_ratio = projection_scalar/edge_dist




                explored_vertices.append(best_edge.dest())
                if best_edge.is_obstacle() or best_edge.opp().is_obstacle():
                    current_edge = best_edge.next().next().opp()
                    while not (current_edge.is_obstacle() and not current_edge.opp().is_obstacle()) and not (current_edge.opp().is_obstacle() and not current_edge.is_obstacle()):
                        current_edge = current_edge.next().next().opp()
                    if current_edge != best_edge and not current_edge in explored_vertices:
                        if len(vertices_path)>0 and len(edges_path)>1:
                            if not edges_path[-1].is_obstacle() and not edges_path[-1].opp().is_obstacle():
                                edge_vector = [current_edge.dest().pos()[0]-current_edge.opp().dest().pos()[0], current_edge.dest().pos()[1]-current_edge.opp().dest().pos()[1]]
                                edge_dist = pow(edge_vector[0],2)+pow(edge_vector[1],2)
                                direction = [end_vertex.pos()[0]-current_edge.opp().dest().pos()[0], end_vertex.pos()[1]-current_edge.opp().dest().pos()[1]]
                                projection_scalar = (edge_vector[0]*direction[0]) + (edge_vector[1]*direction[1])
                                current_ratio = projection_scalar/edge_dist



                                split = True
                                # split_vertices_path=[[vertices_path[-1],best_edge.dest()], [vertices_path[-1],current_edge.dest()]]
                                split_vertices_path=[[best_edge.opp().dest(),best_edge.dest()], [current_edge.opp().dest(),current_edge.dest()]]
                                split_edges_path=[[best_edge], [current_edge]]
                                split_cost_path =[[best_ratio], [current_ratio]]

                              #TODO: check shrotcut on splitting: could collapse on eachside, making the returned vertices list longer et the split
                              #and so the parent vertices list needs to be trimed by as much

                                shortcut = True
                                while shortcut:
                                    shortcut = False
                                    for split_index in list(range(len(split_vertices_path))):
                                        if len(split_vertices_path[split_index]) > 2 :
                                            if type(split_edges_path[split_index][-1].opp().get_face()) is MeshFace:
                                                if split_vertices_path[split_index][-3] == split_edges_path[split_index][-1].opp().next().dest() and not split_edges_path[split_index][-1].opp().next().is_obstacle:
                                                    shortcut = True
                                                    split_edges_path[split_index][-1] = split_edges_path[split_index][-1].opp().next().next()
                                                    split_vertices_path[split_index].remove(split_vertices_path[split_index][-2])
                                                    if len(split_edges_path[split_index]) > 1:
                                                        split_edges_path[split_index].remove(split_edges_path[split_index][-2])
                                                    if len(split_cost_path[split_index]) > 1:
                                                        split_cost_path[split_index].remove(split_cost_path[split_index][-2])
                                                    edge_vector = [split_edges_path[split_index][-1].dest().pos()[0]-split_edges_path[split_index][-1].opp().dest().pos()[0], split_edges_path[split_index][-1].dest().pos()[1]-split_edges_path[split_index][-1].opp().dest().pos()[1]]
                                                    edge_dist = pow(edge_vector[0],2)+pow(edge_vector[1],2)
                                                    direction = [end_vertex.pos()[0]-split_edges_path[split_index][-1].opp().dest().pos()[0], end_vertex.pos()[1]-split_edges_path[split_index][-1].opp().dest().pos()[1]]
                                                    projection_scalar = (edge_vector[0]*direction[0]) + (edge_vector[1]*direction[1])
                                                    split_cost_path[split_index][-1] = projection_scalar/edge_dist

                                                elif type(split_edges_path[split_index][-1].get_face()) is MeshFace:
                                                    if split_vertices_path[split_index][-3] == split_edges_path[split_index][-1].next().dest() and not split_edges_path[split_index][-1].next().next().is_obstacle():
                                                        shortcut = True
                                                        split_edges_path[split_index][-1] = split_edges_path[split_index][-1].next().opp()
                                                        split_vertices_path[split_index].remove(split_vertices_path[split_index][-2])
                                                        if len(split_edges_path[split_index]) > 1:
                                                            split_edges_path[split_index].remove(split_edges_path[split_index][-2])
                                                        if len(split_cost_path[split_index]) > 1:
                                                            split_cost_path[split_index].remove(split_cost_path[split_index][-2])
                                                        edge_vector = [split_edges_path[split_index][-1].dest().pos()[0]-split_edges_path[split_index][-1].opp().dest().pos()[0], split_edges_path[split_index][-1].dest().pos()[1]-split_edges_path[split_index][-1].opp().dest().pos()[1]]
                                                        edge_dist = pow(edge_vector[0],2)+pow(edge_vector[1],2)
                                                        direction = [end_vertex.pos()[0]-split_edges_path[split_index][-1].opp().dest().pos()[0], end_vertex.pos()[1]-split_edges_path[split_index][-1].opp().dest().pos()[1]]
                                                        projection_scalar = (edge_vector[0]*direction[0]) + (edge_vector[1]*direction[1])
                                                        split_cost_path[split_index][-1] = projection_scalar/edge_dist

                                        elif len(vertices_path) - len(split_vertices_path[1-split_index]) > 1:
                                            if type(split_edges_path[split_index][-1].opp().get_face()) is MeshFace:
                                                if vertices_path[-len(split_vertices_path[1-split_index])] == split_edges_path[split_index][-1].opp().next().dest() and not split_edges_path[split_index][-1].opp().next().is_obstacle:
                                                    shortcut = True
                                                    split_vertices_path[split_index][-2] = vertices_path[-len(split_vertices_path[1-split_index])]
                                                    split_edges_path[split_index][-1] = split_edges_path[split_index][-1].opp().next().next()
                                                    split_cost_path[split_index][-1] = cost_path[-len(split_cost_path[1-split_index])]
                                                    split_vertices_path[1-split_index].insert(0,vertices_path[-len(split_vertices_path[1-split_index])])
                                                    if len(edges_path) - len(split_edges_path[1-split_index]) > 1:
                                                        split_edges_path[1-split_index].insert(0,edges_path[-len(split_edges_path[1-split_index])])
                                                    if len(cost_path) - len(split_cost_path[1-split_index]) > 1:
                                                        split_cost_path[1-split_index].insert(0,cost_path[-len(split_cost_path[1-split_index])])
                                                    edge_vector = [split_edges_path[split_index][-1].dest().pos()[0]-split_edges_path[split_index][-1].opp().dest().pos()[0], split_edges_path[split_index][-1].dest().pos()[1]-split_edges_path[split_index][-1].opp().dest().pos()[1]]
                                                    edge_dist = pow(edge_vector[0],2)+pow(edge_vector[1],2)
                                                    direction = [end_vertex.pos()[0]-split_edges_path[split_index][-1].opp().dest().pos()[0], end_vertex.pos()[1]-split_edges_path[split_index][-1].opp().dest().pos()[1]]
                                                    projection_scalar = (edge_vector[0]*direction[0]) + (edge_vector[1]*direction[1])
                                                    split_cost_path[split_index][-1] = projection_scalar/edge_dist

                                                elif type(split_edges_path[split_index][-1].get_face()) is MeshFace:
                                                    if vertices_path[-len(split_vertices_path[1-split_index])] == split_edges_path[split_index][-1].next().dest() and not split_edges_path[split_index][-1].next().next().is_obstacle():
                                                        shortcut = True
                                                        split_vertices_path[split_index][-2] = vertices_path[-len(split_vertices_path[1-split_index])]
                                                        split_edges_path[split_index][-1] = split_edges_path[split_index][-1].next().opp()
                                                        split_cost_path[split_index][-1] = cost_path[-len(split_cost_path[1-split_index])]
                                                        split_vertices_path[1-split_index].insert(0,vertices_path[-len(split_vertices_path[1-split_index])])
                                                        if len(edges_path) - len(split_edges_path[1-split_index]) > 1:
                                                            split_edges_path[1-split_index].insert(0,edges_path[-len(split_edges_path[1-split_index])])
                                                        if len(cost_path) - len(split_cost_path[1-split_index]) > 1:
                                                            split_cost_path[1-split_index].insert(0,cost_path[-len(split_cost_path[1-split_index])])
                                                        edge_vector = [split_edges_path[split_index][-1].dest().pos()[0]-split_edges_path[split_index][-1].opp().dest().pos()[0], split_edges_path[split_index][-1].dest().pos()[1]-split_edges_path[split_index][-1].opp().dest().pos()[1]]
                                                        edge_dist = pow(edge_vector[0],2)+pow(edge_vector[1],2)
                                                        direction = [end_vertex.pos()[0]-split_edges_path[split_index][-1].opp().dest().pos()[0], end_vertex.pos()[1]-split_edges_path[split_index][-1].opp().dest().pos()[1]]
                                                        projection_scalar = (edge_vector[0]*direction[0]) + (edge_vector[1]*direction[1])
                                                        split_cost_path[split_index][-1] = projection_scalar/edge_dist

                                # explored_vertices.append(best_edge.dest())
                                explored_vertices.append(current_edge.dest())

                                vertices_path   =  split_vertices_path
                                edges_path  = split_edges_path
                                cost_path = split_cost_path
                                # vertices_path=[[vertices_path[-1],best_edge.dest()], [vertices_path[-2],current_edge.dest()]]
                                # edges_path=[[best_edge], [current_edge]]
                                # cost_path = [[best_ratio], [current_ratio]]
                                # explored_vertices.append(best_edge.dest())
                                # explored_vertices.append(current_edge.dest())



                if not split:

                    # shortcut = True
                    # while shortcut:
                    #     shortcut = False
                    #     if len(vertices_path) > 2:
                    #         if type(best_edge.opp().get_face()) is MeshFace:
                    #             if vertices_path[-2] == best_edge.opp().next().dest() and not best_edge.opp().next().is_obstacle:
                    #                 shortcut = True
                    #                 best_edge = best_edge.opp().next().next()
                    #                 vertices_path.remove(vertices_path[-1])
                    #                 if len(edges_path) > 0:
                    #                     edges_path.remove(edges_path[-1])
                    #                 if len(cost_path) > 0:
                    #                     cost_path.remove(cost_path[-1])
                    #                 edge_vector = [current_edge.dest().pos()[0]-current_edge.opp().dest().pos()[0], current_edge.dest().pos()[1]-current_edge.opp().dest().pos()[1]]
                    #                 edge_dist = pow(edge_vector[0],2)+pow(edge_vector[1],2)
                    #                 direction = [end_vertex.pos()[0]-current_edge.opp().dest().pos()[0], end_vertex.pos()[1]-current_edge.opp().dest().pos()[1]]
                    #                 projection_scalar = (edge_vector[0]*direction[0]) + (edge_vector[1]*direction[1])
                    #                 best_ratio = projection_scalar/edge_dist
                    #
                    #             elif type(best_edge.get_face()) is MeshFace:
                    #                 if vertices_path[-2] == best_edge.next().dest() and not best_edge.next().next().is_obstacle():
                    #                     shortcut = True
                    #                     best_edge = best_edge.next().opp()
                    #                     vertices_path.remove(vertices_path[-1])
                    #                     if len(edges_path) > 0:
                    #                         edges_path.remove(edges_path[-1])
                    #                     if len(cost_path) > 0:
                    #                         cost_path.remove(cost_path[-1])
                    #                     edge_vector = [current_edge.dest().pos()[0]-current_edge.opp().dest().pos()[0], current_edge.dest().pos()[1]-current_edge.opp().dest().pos()[1]]
                    #                     edge_dist = pow(edge_vector[0],2)+pow(edge_vector[1],2)
                    #                     direction = [end_vertex.pos()[0]-current_edge.opp().dest().pos()[0], end_vertex.pos()[1]-current_edge.opp().dest().pos()[1]]
                    #                     projection_scalar = (edge_vector[0]*direction[0]) + (edge_vector[1]*direction[1])
                    #                     best_ratio = projection_scalar/edge_dist


                    vertices_path.append(best_edge.dest())
                    # explored_vertices.append(best_edge.dest())
                    edges_path.append(best_edge)
                    cost_path.append(best_ratio)

    return {'vertices': vertices_path, 'edges':edges_path, 'explored':explored_vertices, 'split':split, 'cost':cost_path}

def reverse_connection(vertex, tree):
    if vertex in tree['vertices']:
        reverse_vertices_path = []
        reverse_edges_path = []
        reverse_cost_path = []
        for v_index in list(range(tree['vertices'].index(vertex))):
            reverse_vertices_path.insert(0,tree['vertices'][v_index])
            reverse_edges_path.insert(0,tree['edges'][v_index].opp())
            reverse_cost_path.insert(0,tree['cost'][v_index])
        return  {'path':{'vertices':reverse_vertices_path, 'edges':reverse_edges_path, 'cost': reverse_cost_path}, 'tree':None}
    else:
        if tree['split'] is not None:
            if tree['split'][0] is not None:
                split = reverse_connection(vertex, tree['split'][0])
                if split['path'] is not None:
                    split['path']['vertices'] = split['path']['vertices'][:-1]
                    for v_index in list(range(len(tree['vertices']))):
                        split['path']['vertices'].append(tree['vertices'][len(tree['vertices'])-v_index-1])
                    for e_index in list(range(len(tree['edges']))):
                        split['path']['edges'].append(tree['edges'][len(tree['edges'])-e_index-1].opp())
                        split['path']['cost'].append(tree['cost'][len(tree['edges'])-e_index-1])

                    tree['split'][0] = split['tree']

                    # if split['tree'] is not None:
                    #     tree['split'][0] = split['tree']
                    # else:
                    #     tree['split'] = None
                    # if split['tree'] is not None:
                    #     tree['split'][0] = split['tree']
                    # else:
                    #     tree['vertices'].extend(tree['split'][1]['vertices'][1:])
                    #     tree['edges'].extend(tree['split'][1]['edges'])
                    #     tree['cost'].extend(tree['split'][1]['cost'])
                    #     tree['split'] = tree['split'][1]['split']
                    return {'path':split['path'],'tree':tree}
                else:
                    if tree['split'][1] is not None:
                        split = reverse_connection(vertex, tree['split'][1])
                        if split['path'] is not None:
                            split['path']['vertices'] = split['path']['vertices'][:-1]
                            for v_index in list(range(len(tree['vertices']))):
                                split['path']['vertices'].append(tree['vertices'][len(tree['vertices'])-v_index-1])
                            for e_index in list(range(len(tree['edges']))):
                                split['path']['edges'].append(tree['edges'][len(tree['edges'])-e_index-1].opp())
                                split['path']['cost'].append(tree['cost'][len(tree['edges'])-e_index-1])

                            tree['split'][1] = split['tree']


                            # if split['tree'] is not None:
                            #     tree['split'][1] = split['tree']
                            # else:
                            #     tree['split'] = None

                            # if split['tree'] is not None:
                            #     tree['split'][1] = split['tree']
                            # else:
                            #     tree['vertices'].extend(tree['split'][0]['vertices'][:1])
                            #     tree['edges'].extend(tree['split'][0]['edges'])
                            #     tree['cost'].extend(tree['split'][0]['cost'])
                            #     tree['split'] = tree['split'][0]['split']

                            return {'path':split['path'],'tree':tree}
                        else:
                            return {'path':None,'tree':tree}
                    else:
                        return {'path':None,'tree':tree}
            else:
                return {'path':None,'tree':tree}
        else:
            return {'path':None,'tree':tree}

def tree_crawler(target, explored, path_tree, full_tree, connecting_tree, dir='unset'):
    # paths = []
    paths = []
    if path_tree['split'] is None:
        path = pathfinding_step(target, list(explored), list(path_tree['vertices']), list(path_tree['edges']), list(path_tree['cost']))
        if path is not None:
            for explored_vertex in path['explored']:
                if explored_vertex not in explored:
                    explored.append(explored_vertex)
            if path['split'] == False:
                for key in ['vertices', 'edges', 'cost']:
                    path_tree[key] = path[key]
                    full_tree[key] = path_tree[key]

                # path_tree['vertices'] = path['vertices']
                # path_tree['edges'] = path['edges']
                # path_tree['cost'] = path['cost']


                connection = reverse_connection(path['vertices'][-1], connecting_tree)
                if connection['path'] is not None:
                    if connection['tree'] is not None:

                        connecting_tree = connection['tree']
                    else:
                        connecting_tree = {'vertices':[connecting_tree['vertices'][0]], 'edges':[], 'cost': [], 'split':None}

                    connection_path = {'vertices':[], 'edges':[], 'cost':[]}
                    for key in list(connection_path.keys()):
                        connection_path[key].extend(path_tree[key])
                        connection_path[key].extend(connection['path'][key])
                    paths.append(connection_path)
                    path_tree = None

            else:
                path_tree['split'] = [{'vertices':path['vertices'][0], 'edges':path['edges'][0], 'cost':path['cost'][0],'split':None},{'vertices':path['vertices'][1], 'edges':path['edges'][1], 'cost':path['cost'][1],'split':None}]
                full_tree['split'] = path_tree['split']

                if len(path_tree['split'][0]['vertices']) > 2 or len(path_tree['split'][1]['vertices']) > 2:
                    path_tree['vertices'] = path_tree['vertices'][:2-max(len(path_tree['split'][0]['vertices']),len(path_tree['split'][1]['vertices']))]
                    path_tree['edges'] = path_tree['edges'][:1-max(len(path_tree['split'][0]['edges']),len(path_tree['split'][1]['edges']))]

        else:
            path_tree = None

        return {'explored':explored, 'tree':path_tree, 'full_tree':full_tree, 'paths':paths, 'connection':connecting_tree}
    else:
        if path_tree['split'][0] is not None:
            first_split_tree = tree_crawler(target, list(explored), dict(path_tree['split'][0]), dict(full_tree['split'][0]), dict(connecting_tree), dir)
            connecting_tree = first_split_tree['connection']
            for explored_vertex in first_split_tree['explored']:
                if explored_vertex not in explored:
                    explored.append(explored_vertex)
            # explored = sub_tree['explored']
            path_tree['split'][0] = first_split_tree['tree']
            full_tree['split'][0] = first_split_tree['full_tree']

            if first_split_tree['paths'] != []:
                for path_index in list(range(len(first_split_tree['paths']))):
                    connecting_path = {'vertices':[], 'edges':[], 'cost':[]}
                    connecting_path['vertices'].extend(path_tree['vertices'])
                    connecting_path['vertices'].extend(first_split_tree['paths'][path_index]['vertices'][1:])
                    for key in ['edges', 'cost']:
                        connecting_path[key].extend(path_tree[key])
                        connecting_path[key].extend(first_split_tree['paths'][path_index][key])

                    paths.append(connecting_path)

        if path_tree['split'][1] is not None:
            second_split_tree = tree_crawler(target, list(explored), dict(path_tree['split'][1]), dict(full_tree['split'][1]), dict(connecting_tree), dir)
            connecting_tree = second_split_tree['connection']
            for explored_vertex in second_split_tree['explored']:
                if explored_vertex not in explored:
                    explored.append(explored_vertex)
            # explored = sub_tree['explored']
            path_tree['split'][1] = second_split_tree['tree']
            full_tree['split'][1] = second_split_tree['full_tree']

            if second_split_tree['paths'] != []:
                for path_index in list(range(len(second_split_tree['paths']))):
                    connecting_path = {'vertices':[], 'edges':[], 'cost':[]}
                    connecting_path['vertices'].extend(path_tree['vertices'])
                    connecting_path['vertices'].extend(second_split_tree['paths'][path_index]['vertices'][1:])
                    for key in ['edges', 'cost']:
                        connecting_path[key].extend(path_tree[key])
                        connecting_path[key].extend(second_split_tree['paths'][path_index][key])
                    paths.append(connecting_path)

        if path_tree['split'][0] is None and path_tree['split'][1] is None:
            path_tree['split'] = None
        # elif path_tree['split'][0] is None:
        #     print('collapse0')
        #     path_tree['vertices'].extend(path_tree['split'][1]['vertices'][1:])
        #     path_tree['edges'].extend(path_tree['split'][1]['edges'])
        #     path_tree['cost'].extend(path_tree['split'][1]['cost'])
        #     path_tree['split'] = path_tree['split'][1]['split']
        # elif path_tree['split'][1] is None:
        #     print('collapse1')
        #     print(path_tree['vertices'][-1].pos())
        #     print(path_tree['split'][0]['vertices'][-1].pos())
        #     path_tree['vertices'].extend(path_tree['split'][0]['vertices'][1:])
        #     path_tree['edges'].extend(path_tree['split'][0]['edges'])
        #     path_tree['cost'].extend(path_tree['split'][0]['cost'])
        #     path_tree['split'] = path_tree['split'][0]['split']


        return {'explored':explored, 'tree':path_tree, 'full_tree':full_tree, 'paths':paths, 'connection':connecting_tree}

def path_line_of_sight(path):
    pass


def json_load(filename, grid_len=None):
    points_list =  json.loads(open(filename))
    return MeshMap(points_list, grid_len)


if __name__ == '__main__':

    #an utility display function
    def display_faces(faces, color=None):
        fig = pyplot.figure()
        ar = fig.add_subplot(1,1,1, aspect='equal')
        for face in faces:
            if type(face) is MeshFace:
                edge = face.edge()
                if type(edge) is MeshFace:
                    print('nimp!')
                    print(edge.next().dest().pos(),edge.next().opp().dest().pos())
                else:
                    tx = []
                    ty = []
                    for iter in [0,1,2,3]:
                        tx.append(edge.dest().pos()[0])
                        ty.append(edge.dest().pos()[1])

                        edge = edge.next()
            if color is not None:
                ar.plot(tx,ty,color)
            else:
                ar.plot(tx,ty)
        return ar

    def display_map(meshmap):
        segz = []
        wegz = []
        for f in meshmap.get_faces():
            if type(f) is MeshFace:
                e = f.edge()
                for i in [0,1,2,3]:
                    if e.is_obstacle():
                        wegz.append([e.opp().dest().pos(), e.dest().pos()])
                    else:
                        segz.append([e.opp().dest().pos(), e.dest().pos()])
                    e = e.next()
        fig = pyplot.figure()
        ar = fig.add_subplot(1,1,1, aspect='equal')
        sar = array(segz)
        for h in sar:
            tx = h[:,0]
            ty = h[:,1]
            ar.plot(tx,ty, color='r')
        war = array(wegz)
        for w in war:
            wx = w[:,0]
            wy = w[:,1]
            ar.plot(wx,wy, color='k')
        return ar

    def display_tree(tree, fig, color='g'):
        if type(tree) is not dict:
            print(tree)
        pegz = []
        for vindex in list(range(len(tree['vertices'])-1)):
            pegz.append([tree['vertices'][vindex].pos(), tree['vertices'][vindex+1].pos()])
        par = array(pegz)
        for g in par:
            rx = g[:,0]
            ry = g[:,1]
            fig.plot(rx,ry,color)
        if 'split' in list(tree.keys()):
            if tree['split'] is not None:
                if tree['split'][0] is not None:
                    fig = display_tree(tree['split'][0], fig)
                if tree['split'][1] is not None:
                    fig = display_tree(tree['split'][1], fig)
        return fig


    #List of points that will be the base vertices of the triangle mesh
    points_list=[[1.,2.],[1.6,2.3],[2.,1.],[2.2,8.],[2.7,1.8],[3.5,5.],[4.2,8.],[4.6,3.6],[5.,3.],[5.4,2.],[5.6,5.],[6.,4.2], [6.2, 5.8]]

    # map_patch = []
    # for i in list(range(3)):
    #     for j in list(range(3)):
    #         for pnum in list(range(len(points_list))):
    #             map_patch.append([points_list[pnum][0]+i*8,points_list[pnum][1]+j*8])
    # map_points = []
    # for i in list(range(3)):
    #     for j in list(range(3)):
    #         for pnum in list(range(len(map_patch))):
    #             map_points.append([map_patch[pnum][0]+(i)*20,map_patch[pnum][1]+(j)*20])

    map_points = []
    # dist = 200
    # num = 5
    # for i in list(range(3*num+1)):
    #     for j in list(range(2*num+1)):
    #         map_points.append([i*dist,j*dist])
    #         map_points.append([i*dist+100,j*dist+100])

    dist = 200
    num = 3
    for i in list(range(3*num+1)):
        for j in list(range(2*num+1)):
            map_points.append([i*dist,j*dist])
            map_points.append([i*dist+100,j*dist+100])

    print(len(map_points))

    print(datetime.now())
    grid_len = 10

    #To generate the Mesh, the points need to be sorted by x then y
    sorted_map = sorted(map_points, key=lambda k: [k[0], k[1]])

    map=MeshMap(sorted_map,grid_len)

    # display_faces(map.get_faces())
    # pyplot.show()

    #Solids are defined by a list of points given in CLOCKORK ORDER
    # map.add_solid([[-5,0],[3005,0], [3005,-5], [-5,-5]])
    # map.add_solid([[-5,2005],[3005,2005], [3005,2000], [-5,2000]])

    # map.add_solid([[200,200],[200,400],[400,400],[400,200]])
    # map.add_solid([[1400,400],[1400,800],[1800,800],[1800,400]])

    # map.add_solid([[480,50],[480,1040],[560,1040],[560,50]])
    # map.add_solid([[720,500],[720,1900],[920,1900],[920,500]])


    # map.add_solid([[920,500],[920,800],[1300,800],[1300,500]])


    map.add_solid([[480,20],[480,640],[560,640],[560,20]])
    map.add_solid([[720,500],[720,1200],[920,1200],[920,500]])
    map.add_solid([[920,500],[920,800],[1300,800],[1300,500]])

    print(datetime.now())

    # display_map(map)

    # paths = map.d_star([60,80],[1480,1220],15)
    #
    # print(datetime.now())
    #
    # print(len(paths[0]), len(paths[1]))
    # for path in paths[0]:
    #     display_tree(path, display_map(map))
    #
    # for path in paths[1]:
    #     display_tree(path, display_map(map))

    # map.add_solid([[920,500],[920,800],[1300,800],[1300,500]])
    # paths = map.d_star([60,80],[1480,1220],5)

    # display_map(map)
    # pyplot.show()

    # paths = map.d_star([60,80],[1480,1220],15)

    # display_map(map)
    # for path in paths[0]:
    #     display_tree(path, display_map(map))
    # for path in paths[1]:
    #     display_tree(path, display_map(map))
    # pyplot.show()

    # display_map(map)
    # pyplot.show()
    #
    #
    # map.delete_solid(map.get_solids()[1])
    # map.delete_solid(map.get_solids()[0])

    # display_map(map)

    map.add_solid([[200,200],[200,400],[400,400],[400,200]])

    display_map(map)
    pyplot.show()

    paths = map.d_star([60,80],[1480,1220],25)

    for path in paths[0]:
        display_tree(path, display_map(map))
    for path in paths[1]:
        display_tree(path, display_map(map))

    print('avoid')
    avoid = []
    avoid.extend(map.avoid_verices([600, 1000], 200))
    avoid.extend(map.avoid_verices([900, 300], 180))

    paths = map.d_star([60,80],[1480,1220],15, avoid_vertices=avoid)

    # display_map(map)

    for path in paths[0]:
        ar = display_tree(path, display_map(map))
        ar.add_patch(pyplot.Circle((600, 1000), 200))
        ar.add_patch(pyplot.Circle((900, 300), 180))
    for path in paths[1]:
        ar = display_tree(path, display_map(map))
        ar.add_patch(pyplot.Circle((600, 1000), 200))
        ar.add_patch(pyplot.Circle((900, 300), 180))

    pyplot.show()

    print(datetime.now())


