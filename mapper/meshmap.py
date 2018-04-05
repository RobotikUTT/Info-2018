from matplotlib import pyplot
from numpy.linalg import det
from numpy import array, concatenate, append
from math import pow
from datetime import datetime


class MeshPlane():

    def __init__(self, width, length):
        self._width = width
        self._length = length

class MeshObj():
    def __init__(self):
        pass


class MeshVertex():
    def __init__(self, x, y, edge=None):
        self._pos = [x,y]
        self._edge_start = edge

    def get_pos(self):
        return self._pos

    def set_edge(self, edge):
        self._edge_start = edge

    def get_edge(self):
        return self._edge_start

class GhostVertex():
    def __init__(self):
        pass

class GhostEdge():
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

    def get_opp(self):
        return self._opp_edge

    def set_next(self, edge):
        self._next_edge = edge

    def get_next(self):
        return self._next_edge


class MeshEdge():
    def __init__(self, vertex_dest, face=None,  opp_edge=None, next_edge=None):
        self._vertex_dest = vertex_dest
        self._face = face
        self._opp_edge = opp_edge
        self._next_edge= next_edge

    def get_dest(self):
        return self._vertex_dest

    def set_dest(self, dest):
        self._vertex_dest = dest

    def set_face(self, face):
        self._face = face

    def get_face(self):
        return self._face

    def set_opp(self, edge):
        self._opp_edge = edge

    def get_opp(self):
        return self._opp_edge

    def set_next(self, edge):
        self._next_edge = edge

    def get_next(self):
        return self._next_edge

class GhostFace():
    def __init__(self, edge=None):
        self._edge = edge

    # def __del__(self):
    #     print('ghost face deleted')

    def set_edge(self, edge):
        self._edge = edge

    def get_edge(self):
        return self._edge

class MeshFace():
    def __init__(self, edge=None):
        self._edge = edge

    def set_edge(self, edge):
        self._edge = edge

    def get_edge(self):
        return self._edge


def order_points_list(points_list):
    ord_points = [points_list[0]]
    for point_num in list(range(len(points_list)-1)):
        ord_index = 0
        added = False
        while not added and ord_index < len(ord_points):
                if points_list[point_num+1][0] > ord_points[ord_index][0]:

                    ord_index += 1
                elif points_list[point_num+1][0] == ord_points[ord_index][0]:
                    if points_list[point_num+1][1] > ord_points[ord_index][1]:
                        ord_index += 1
                    else:
                        ord_points.insert(ord_index, points_list[point_num+1])
                        added = True
                else:
                    ord_points.insert(ord_index, points_list[point_num+1])
                    added = True
        if not added:
            ord_points.append(points_list[point_num+1])
            added = True
    return ord_points

def clockwise_check(a,b,c):
    surface_sum = sum([(b.get_pos()[0]-a.get_pos()[0])*(b.get_pos()[1]+a.get_pos()[1]),
                       (c.get_pos()[0]-b.get_pos()[0])*(c.get_pos()[1]+b.get_pos()[1]),
                       (a.get_pos()[0]-c.get_pos()[0])*(a.get_pos()[1]+c.get_pos()[1])])
    return surface_sum > 0

def colinear_check(a,b,c):
    surface_sum = sum([(b.get_pos()[0]-a.get_pos()[0])*(b.get_pos()[1]+a.get_pos()[1]),
                       (c.get_pos()[0]-b.get_pos()[0])*(c.get_pos()[1]+b.get_pos()[1]),
                       (a.get_pos()[0]-c.get_pos()[0])*(a.get_pos()[1]+c.get_pos()[1])])
    return surface_sum == 0

def circumcircle_check(v1,v2,v3, x):
    M = []
    for v in [v1,v2,v3,x]:
        M.append([v.get_pos()[0],v.get_pos()[1],pow(v.get_pos()[0],2)+pow(v.get_pos()[1],2),1])
    return det(M) > 0

def add_ghost_vertex(edge):
    #complete an edge with a gohst vertex, creating 2 ghost edges and a ghost face

    facex = GhostFace(edge)
    edge.set_face(facex)

    a = edge.get_opp().get_dest()
    x = GhostVertex()

    edge_bx = GhostEdge(x, face=facex)
    edge.set_next(edge_bx)

    edge_xa = GhostEdge(a, face=facex, next_edge=edge)
    edge_bx.set_next(edge_xa)

    return {'vertex': x, 'face':facex}

def  create_isolated_edge(a,b):
    #Create the opposite edges connecting the veertices
    edge_ab = MeshEdge(b)
    edge_ba = MeshEdge(a, opp_edge=edge_ab)
    edge_ab.set_opp(edge_ba)
    a.set_edge(edge_ab)
    b.set_edge(edge_ba)

    # add ghost vertices, edges, and faces to obtain a quad-faced polygon
    up_ghost = add_ghost_vertex(edge_ab)

    down_ghost = add_ghost_vertex(edge_ba)

    #link the lateral edges of the ghost faces as opposite edges
    edge_ab.get_next().set_opp(edge_ba.get_next().get_next())
    edge_ba.get_next().get_next().set_opp(edge_ab.get_next())
    edge_ab.get_next().get_next().set_opp(edge_ba.get_next())
    edge_ba.get_next().set_opp(edge_ab.get_next().get_next())

    return {'vertices':[up_ghost['vertex'],down_ghost['vertex']],
            'faces':[up_ghost['face'],down_ghost['face']]}

def create_isolated_face(a,b,c):
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

    edge_ba.get_next().get_next().set_opp(edge_cb.get_next())
    edge_cb.get_next().set_opp(edge_ba.get_next().get_next())
    edge_cb.get_next().get_next().set_opp(edge_ac.get_next())
    edge_ac.get_next().set_opp(edge_cb.get_next().get_next())
    edge_ac.get_next().get_next().set_opp(edge_ba.get_next())
    edge_ba.get_next().set_opp(edge_ac.get_next().get_next())

    return {'vertices':[first_ghost['vertex'],second_ghost['vertex'],third_ghost['vertex']],
            'faces':[face,first_ghost['face'],second_ghost['face'],third_ghost['face']]}


def merge_edges(start_edge, faces):
    #######
    #Identify vertices around the starting edge
    startL = start_edge.get_dest()
    startR = start_edge.get_opp().get_dest()
    #initiate the duo of vertices candidates, starting with the most exterior edge (minimal angle)
    #
    edgeL = start_edge.get_next().get_opp().get_next().get_opp()
    Lcandidate = edgeL.get_opp().get_dest()
    next_edgeL = edgeL.get_next().get_opp()

    Lfound = False
    Ledges_crossed = 0

    #Then check if the candidate vertex complies to the delaunay criterion

    while clockwise_check(startR, startL, Lcandidate) and Lfound==False:
        if type(next_edgeL) is MeshEdge:
            next_Lcandidate = next_edgeL.get_opp().get_dest()

            if not circumcircle_check(startR, startL, Lcandidate, next_Lcandidate):
                #if not, change the candidate couple with the next vertices connected to startL and iterate
                edgeL = next_edgeL
                Lcandidate = next_Lcandidate
                next_edgeL = edgeL.get_next().get_opp()
                Ledges_crossed += 1
            else:
                #if the delaunay criterion is valid, register a left candidate
                Lfound = True
        else:
            Lfound = True


    edgeR = start_edge.get_next().get_next().get_opp().get_next().get_next().get_opp()


    Rcandidate = edgeR.get_dest()
    next_edgeR = edgeR.get_next().get_next().get_opp()

    Rfound = False
    Redges_crossed = 0



    while clockwise_check(startR, startL, Rcandidate) and Rfound==False:

        if type(next_edgeR) is MeshEdge:
            next_Rcandidate = next_edgeR.get_dest()
            if not circumcircle_check(startR, startL, Rcandidate, next_Rcandidate):
                edgeR = next_edgeR
                Rcandidate = next_Rcandidate
                next_edgeR = edgeR.get_next().get_next().get_opp()
                Redges_crossed += 1
            else:
                Rfound = True
        else:
            Rfound = True

    elected_edge = None
    elected_vertex = None
    if Rfound or Lfound:
        if Rfound and Lfound:
            #if vertices were found both left and right, one must be elected
            #it is done with the delaunay criterion
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
            else:
                elected_vertex = Rcandidate
                elected_edge = edgeR

    if elected_vertex is None:
        print('aie')
    if elected_vertex is not None and elected_edge is not None:
        #finally, connect the elected vertex with edges to form the new face of the mesh
        #this will replace (unreference) ghost triangles on concerned edges
        if Lfound:

            if Ledges_crossed > 0:
                current_edge = start_edge.get_next().get_opp().get_next().get_opp()
                a,b = elected_vertex, startR
                for edge_num in list(range(Ledges_crossed)):
                    k,l = current_edge.get_dest(), current_edge.get_opp().get_dest()
                    next_edge = current_edge.get_next().get_opp()
                    if intersect_check(a.get_pos()[0],a.get_pos()[1],
                                       b.get_pos()[0],b.get_pos()[1],
                                       k.get_pos()[0],k.get_pos()[1],
                                       l.get_pos()[0],l.get_pos()[1]):
                        # if current_edge.get_next().get_next().get_dest().get_edge() == current_edge:
                        #     current_edge.get_next().get_next().get_dest().set_edge(current_edge.get_next().get_next().get_opp())

                        current_edge.get_opp().get_next().get_next().set_next(current_edge.get_next().get_next())
                        current_edge.get_next().get_next().set_next(current_edge.get_opp().get_next())
                        current_edge.get_next().get_next().get_dest().set_edge(current_edge.get_next().get_next().get_opp())

                        current_edge.get_face().set_edge(None)

                        current_edge.get_next().get_next().set_face(current_edge.get_opp().get_face())
                        current_edge.get_next().set_next(None)
                        current_edge.get_next().set_face(None)
                        del faces[faces.index(current_edge.get_face())]
                        current_edge.set_face(None)
                        new_ghost = add_ghost_vertex(current_edge.get_next())
                        faces.append(new_ghost)
                        current_edge.get_next().get_next().set_opp(current_edge.get_opp().get_next().get_next())
                        current_edge.get_next().get_next().get_next().set_opp(start_edge.get_next())
                        current_edge.get_opp().set_next(None)
                        current_edge.get_opp().set_face(None)
                        current_edge.get_opp().set_opp(None)

                        current_edge.set_opp(None)
                        current_edge.set_next(None)
                        current_edge.set_dest(None)
                        current_edge.set_face(None)
                        del current_edge

                    current_edge = next_edge

            #create the new edge and face
            new_edge = MeshEdge(startR, next_edge=start_edge)
            new_face = MeshFace(start_edge)

            #add an opposite ghost face
            new_opp_edge = MeshEdge(elected_vertex, opp_edge=new_edge)
            new_edge.set_opp(new_opp_edge)


            #unset left ghost face reference to an edge
            elected_edge.get_opp().get_face().set_edge(None)
            #set the face pointer to new face
            elected_edge.get_opp().set_face(new_face)
            #BUG HERE, Face stays in faces list but has no edge (none)
            elected_edge.get_opp().get_next().set_face(start_edge.get_face())
            #delete refeences of the left hand ghost edge to remove
            elected_edge.get_opp().get_next().get_next().set_next(None)
            elected_edge.get_opp().get_next().get_next().set_opp(None)
            elected_edge.get_opp().get_next().get_next().set_face(None)
            #set the next attribute to link both gohst triangles
            elected_edge.get_opp().get_next().set_next(start_edge.get_next().get_next())

            #close the new merged ghost triangle
            elected_edge.get_opp().get_next().get_next().set_next(new_opp_edge)
            new_opp_edge.set_next(elected_edge.get_opp().get_next())
            new_opp_edge.set_face(start_edge.get_face())
            new_opp_edge.get_face().set_edge(new_opp_edge)

            #delete all references ofthe righ hand ghost edge to merge
            start_edge.get_next().set_face(None)
            start_edge.get_next().set_next(None)
            start_edge.get_next().set_opp(None)

            start_edge.set_next(elected_edge.get_opp())
            elected_edge.get_opp().set_next(new_edge)

            start_edge.set_face(new_face)
            new_edge.set_face(new_face)

            #set the vertices so that they point along the edges of the created hull (clockwise)

            if type(elected_edge) is MeshEdge:
                elected_vertex.set_edge(new_edge)
            if type(startL.get_edge()) is GhostEdge:
                startL.set_edge(elected_edge)

            faces.insert(0,new_face)
            merged = {'faces': faces, 'edge': new_edge}
            return merged

        else:

            if Redges_crossed > 0:
                current_edge = start_edge.get_next().get_next().get_opp().get_next().get_next().get_opp()
                a,b = elected_vertex, startL
                for edge_num in list(range(Redges_crossed)):
                    k,l = current_edge.get_dest(), current_edge.get_opp().get_dest()
                    next_edge = current_edge.get_next().get_next().get_opp()
                    if intersect_check(a.get_pos()[0],a.get_pos()[1],
                                       b.get_pos()[0],b.get_pos()[1],
                                       k.get_pos()[0],k.get_pos()[1],
                                       l.get_pos()[0],l.get_pos()[1]):

                        current_edge.get_next().get_next().set_next(None)
                        current_edge.get_next().get_next().set_face(None)
                        new_ghost = add_ghost_vertex(current_edge.get_next().get_next())
                        faces.append(new_ghost['face'])
                        current_edge.get_next().get_next().get_next().get_next().set_opp(current_edge.get_opp().get_next())
                        current_edge.get_next().get_next().get_next().set_opp(start_edge.get_next().get_next())
                        current_edge.get_next().set_next(current_edge.get_opp().get_next())
                        current_edge.get_next().set_face(current_edge.get_opp().get_face())
                        current_edge.get_opp().get_next().get_next().set_next(current_edge.get_next())



                        current_edge.get_face().set_edge(None)
                        del faces[faces.index(current_edge.get_face())]
                        current_edge.set_face(None)
                        current_edge.get_opp().set_next(None)
                        current_edge.get_opp().set_opp(None)
                        current_edge.get_opp().set_face(None)
                        current_edge.get_opp().set_dest(None)
                        current_edge.set_opp(None)
                        current_edge.set_next(None)
                        current_edge.set_dest(None)

                        del current_edge
                    current_edge = next_edge

            #create new edge
            new_edge = MeshEdge(elected_vertex, next_edge=elected_edge.get_opp())
            new_face = MeshFace(start_edge)

            new_opp_edge = MeshEdge(startL, opp_edge=new_edge)
            new_edge.set_opp(new_opp_edge)
            # add_ghost_vertex(new_opp_edge)

            #delete references to left ghost face
            start_edge.get_face().set_edge(None)

            start_edge.get_next().set_face(elected_edge.get_opp().get_face())
            #delete references to left side of linking ghost face
            start_edge.get_next().get_next().set_face(None)
            start_edge.get_next().get_next().set_opp(None)
            start_edge.get_next().get_next().set_next(None)

            #link both ghost by a next pointer
            start_edge.get_next().set_next(elected_edge.get_opp().get_next().get_next())
            #then delete right hand side of the ghodt edge to remove
            elected_edge.get_opp().get_next().set_face(None)
            elected_edge.get_opp().get_next().set_opp(None)
            elected_edge.get_opp().get_next().set_next(None)

            #redefine the last edge of the ghost triangle
            start_edge.get_next().get_next().set_next(new_opp_edge)
            new_opp_edge.set_next(start_edge.get_next())

            new_opp_edge.set_face(elected_edge.get_opp().get_face())
            new_opp_edge.get_face().set_edge(new_opp_edge)

            start_edge.set_face(new_face)
            start_edge.set_next(new_edge)
            elected_edge.get_opp().set_next(start_edge)
            elected_edge.get_opp().set_face(new_face)
            new_edge.set_face(new_face)

            if type(elected_edge) is GhostEdge:
                elected_vertex.set_edge(elected_edge)
            if type(startL.get_edge()) is MeshEdge:
                startL.set_edge(new_edge)

            faces.insert(0,new_face)

            # bugs = 0
            # for fnum in list(range(len(faces))):
            #     if type(faces[fnum]) is MeshFace:
            #         if type(faces[fnum].get_edge()) is not MeshEdge:
            #             bugs +=1
            #             print(fnum)
            #             print('edge ' + str(1))
            #         if type(faces[fnum].get_edge().get_next()) is not MeshEdge:
            #             bugs +=1
            #             print(fnum)
            #             print('edge ' + str(2))
            #             print(faces[fnum].get_edge().get_next())
            #         if type(faces[fnum].get_edge().get_next().get_next()) is not MeshEdge:
            #             bugs +=1
            #             print(fnum)
            #             print('edge ' + str(3))
            #             print(faces[fnum].get_edge().get_next().get_next())
            # print('bugs')
            # print(bugs)

            merged = {'faces': faces, 'edge': new_edge}
            return merged

def merge_mesh(meshL, meshR):

    #TODO use the ghost edges instead of vertices to find uper and lower bound of Lmesh hull and Rmesh hull
    #because vertices cannot point to edges if it's isolated segments

    vertices = []
    for v in meshL['vertices']:
        if type(v) is MeshVertex:
            vertices.append(v)
    for v in meshR['vertices']:
        if type(v) is MeshVertex:
            vertices.append(v)
    for v in meshL['vertices']:
        if type(v) is GhostVertex:
            vertices.append(v)
    for v in meshR['vertices']:
        if type(v) is GhostVertex:
            vertices.append(v)

    faces = []
    for f in meshL['faces']:
        if type(f) is MeshFace:
            faces.append(f)
    for f in meshR['faces']:
        if type(f) is MeshFace:
            faces.append(f)
    for f in meshL['faces']:
        if type(f) is GhostFace:
            faces.append(f)
    for f in meshR['faces']:
        if type(f) is GhostFace:
            faces.append(f)

    L_flat = len([f for f in meshL['faces'] if type(f) is  MeshFace]) == 0
    R_flat = len([f for f in meshR['faces'] if type(f) is  MeshFace]) == 0

    #Find the upper and lower boundary vertices of each convex hull
    upper_L_vertex = meshL['vertices'][0]
    upper_L_height = meshL['vertices'][0].get_pos()[1]
    for vertex in meshL['vertices']:
        if type(vertex) is MeshVertex:
            if vertex.get_pos()[1] > upper_L_height:
                upper_L_vertex = vertex
                upper_L_height = vertex.get_pos()[1]
            elif vertex.get_pos()[1] == upper_L_height and vertex.get_pos()[0] > upper_L_vertex.get_pos()[0]:
                upper_L_vertex = vertex

    lower_L_vertex = meshL['vertices'][0]
    lower_L_height = meshL['vertices'][0].get_pos()[1]
    for vertex in meshL['vertices']:
        if type(vertex) is MeshVertex:
            if vertex.get_pos()[1] < lower_L_height:
                lower_L_vertex = vertex
                lower_L_height = vertex.get_pos()[1]
            elif vertex.get_pos()[1] == lower_L_height and vertex.get_pos()[0] > lower_L_vertex.get_pos()[0]:
                lower_L_vertex = vertex

    upper_R_vertex = meshR['vertices'][0]
    upper_R_height = meshR['vertices'][0].get_pos()[1]
    for vertex in meshR['vertices']:
        if type(vertex) is MeshVertex:
            if vertex.get_pos()[1] > upper_R_height:
                upper_R_vertex = vertex
                upper_R_height = vertex.get_pos()[1]
            elif vertex.get_pos()[1] == upper_R_height and vertex.get_pos()[0] < upper_R_vertex.get_pos()[0]:
                upper_R_vertex = vertex

    lower_R_vertex = meshR['vertices'][0]
    lower_R_height = meshR['vertices'][0].get_pos()[1]
    for vertex in meshR['vertices']:
        if type(vertex) is MeshVertex:
            if vertex.get_pos()[1] < lower_R_height:
                lower_R_vertex = vertex
                lower_R_height = vertex.get_pos()[1]
            elif vertex.get_pos()[1] == lower_R_height and vertex.get_pos()[0] < lower_R_vertex.get_pos()[0]:
                lower_R_vertex = vertex

    if L_flat and R_flat and meshL['vertices'][0].get_pos()[0] == meshR['vertices'][0].get_pos()[0]:
        upper_L_edge = upper_L_vertex.get_edge()
        lower_R_edge = lower_R_vertex.get_edge()
        ghosts = create_isolated_edge(upper_L_vertex, lower_R_vertex)

        #TODO check that: ghost linking on a colinear vertical line
        upper_L_vertex.get_edge().get_next().set_opp(lower_R_edge.get_next().get_next())
        lower_R_edge.get_next().get_next().set_opp(upper_L_vertex.get_edge().get_next())
        lower_R_vertex.get_edge().get_next().get_next().set_opp(lower_R_edge.get_opp().get_next())
        lower_R_edge.get_opp().get_next().set_opp(lower_R_vertex.get_edge().get_next().get_next())

        upper_L_vertex.get_edge().get_next().get_next().set_opp(upper_L_edge.get_opp().get_next())
        upper_L_edge.get_opp().get_next().set_opp(upper_L_vertex.get_edge().get_next().get_next())
        lower_R_vertex.get_edge().get_next().set_opp(upper_L_edge.get_next().get_next())
        upper_L_edge.get_next().get_next().set_opp(lower_R_vertex.get_edge().get_next())

        #fill up the resulting lists
        vertices.append(ghosts['vertices'][0])
        vertices.append(ghosts['vertices'][1])

        faces.append(ghosts['faces'][0])
        faces.append(ghosts['faces'][1])

        return {'vertices': vertices, 'faces':faces}
    else:
        start_R = lower_R_vertex
        start_L = lower_L_vertex
        stop_R = upper_R_vertex
        stop_L = upper_L_vertex

        if R_flat and not L_flat:
            aligned_L = None
            for vertex in meshL['vertices']:
                if type(vertex) is MeshVertex:
                    if vertex.get_pos()[0] == lower_R_vertex.get_pos()[0]:
                        if aligned_L is None:
                            aligned_L = vertex
                        elif vertex.get_pos()[1] > aligned_L.get_pos()[1]:
                            aligned_L = vertex
            if aligned_L is not None:
                start_L = aligned_L

        if L_flat and not R_flat:
            aligned_R = None
            for vertex in meshR['vertices']:
                if type(vertex) is MeshVertex:
                    if vertex.get_pos()[0] == upper_L_vertex.get_pos()[0]:
                        if aligned_R is None:
                            aligned_R = vertex
                        elif vertex.get_pos()[1] < aligned_R.get_pos()[1]:
                            aligned_R = vertex
            if aligned_R is not None:
                stop_R = aligned_R

        #create an edge between both lower points and its opposite
        merge_edge_start = MeshEdge(start_L)
        merge_edge_opp = MeshEdge(start_R,opp_edge=merge_edge_start)

        merge_edge_start.set_opp(merge_edge_opp)

        #add the ghost faces under and over them
        over_ghost_face = add_ghost_vertex(merge_edge_start)
        under_ghost_face = add_ghost_vertex(merge_edge_opp)

        faces.append(over_ghost_face['face'])
        faces.append(under_ghost_face['face'])

        #link the opposites edges of the ghost faces bellow and over
        merge_edge_start.get_next().set_opp(start_L.get_edge().get_opp().get_next().get_opp())
        start_L.get_edge().get_opp().get_next().get_opp().set_opp(merge_edge_start.get_next())
        merge_edge_opp.get_next().get_next().set_opp(start_L.get_edge().get_opp().get_next())
        start_L.get_edge().get_opp().get_next().set_opp(merge_edge_opp.get_next().get_next())
        merge_edge_opp.get_next().set_opp(start_R.get_edge().get_opp().get_next().get_opp())
        start_R.get_edge().get_opp().get_next().get_opp().set_opp(merge_edge_opp.get_next())
        merge_edge_start.get_next().get_next().set_opp(start_R.get_edge().get_opp().get_next())
        start_R.get_edge().get_opp().get_next().set_opp(merge_edge_start.get_next().get_next())

        start_R.set_edge(merge_edge_start)
        current_edge = merge_edge_start

        L_edge = start_L.get_edge()

        while not (current_edge.get_dest() == stop_L and current_edge.get_opp().get_dest() == stop_R):

            merged = merge_edges(current_edge, faces)
            if merged is None:
                print('problem?')
                segz = []
                for fnum in list(range(len(faces))):
                    if type(faces[fnum]) is MeshFace:
                        e = faces[fnum].get_edge()
                        for i in [0,1,2,3]:
                            if type(e) is GhostEdge:
                                print('ge')
                                print(fnum)
                                print(i)
                            segz.append([e.get_dest().get_pos(), e.get_next().get_dest().get_pos()])
                            e = e.get_next()
                fig = pyplot.figure()
                ar = fig.add_subplot(1,1,1)
                sar = array(segz)
                for h in sar:
                    tx = h[:,0]
                    ty = h[:,1]
                    ar.plot(tx,ty)
                pyplot.show()

            faces = merged['faces']

            current_edge = merged['edge'].get_opp()

        start_L.set_edge(L_edge)

        # segz = []
        # for f in faces:
        #     if type(f) is MeshFace:
        #         e = f.get_edge()
        #         for i in [0,1,2,3]:
        #             segz.append([e.get_dest().get_pos(), e.get_next().get_dest().get_pos()])
        #             e = e.get_next()
        #
        # fig = pyplot.figure()
        # ar = fig.add_subplot(1,1,1)
        # sar = array(segz)
        # for h in sar:
        #     tx = h[:,0]
        #     ty = h[:,1]
        #     ar.plot(tx,ty)
        # pyplot.show()

        return {'faces': faces, 'vertices':vertices}

def divide_and_conquer(ordered_points):

    #TODO: Handle the lonesome point
    if len(ordered_points) ==1:
        print('nope')
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

        return {'vertices': vertices, 'faces':faces}

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

                if not (b.get_pos()[0]<= max(a.get_pos()[0], c.get_pos()[0]) \
                                and b.get_pos()[0] >= min(a.get_pos()[0], c.get_pos()[0]) \
                                and b.get_pos()[1]<= max(a.get_pos()[1], c.get_pos()[1]) \
                                and b.get_pos()[1] >= min(a.get_pos()[1], c.get_pos()[1])):

                    b,c = c,b

                first_ghosts = create_isolated_edge(a,b)
                other_ghosts = create_isolated_edge(b,c)


                first_ghosts['faces'][0].get_edge().get_next().set_opp(other_ghosts['faces'][0].get_edge().get_next().get_next())
                other_ghosts['faces'][0].get_edge().get_next().get_next().set_opp(first_ghosts['faces'][0].get_edge().get_next())

                first_ghosts['faces'][1].get_edge().get_next().get_next().set_opp(other_ghosts['faces'][1].get_edge().get_next())
                other_ghosts['faces'][1].get_edge().get_next().set_opp(first_ghosts['faces'][1].get_edge().get_next().get_next())

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

                return {'vertices': vertices, 'faces':faces}

            else:
                b,c = c,b

        tri = create_isolated_face(a,b,c)

        vertices.append(a)
        vertices.append(b)
        vertices.append(c)
        vertices.append(tri['vertices'][0])
        vertices.append(tri['vertices'][0])
        vertices.append(tri['vertices'][0])

        faces.append(tri['faces'][0])
        faces.append(tri['faces'][1])
        faces.append(tri['faces'][2])
        faces.append(tri['faces'][3])

        # segz = []
        # for f in faces:
        #     if type(f) is MeshFace:
        #         e = f.get_edge()
        #         for i in [0,1,2,3]:
        #             segz.append([e.get_dest().get_pos(), e.get_next().get_dest().get_pos()])
        #             e = e.get_next()
        #
        # fig = pyplot.figure()
        # ar = fig.add_subplot(1,1,1)
        # sar = array(segz)
        # for h in sar:
        #     tx = h[:,0]
        #     ty = h[:,1]
        #     ar.plot(tx,ty)
        # pyplot.show()

        return {'vertices': vertices, 'faces':faces}

    #case of overpopulation of points; halfen the list then call recursively
    else:

        Lmesh = divide_and_conquer(ordered_points[:len(ordered_points)//2])
        Rmesh = divide_and_conquer(ordered_points[len(ordered_points)//2:])

        return merge_mesh(Lmesh,Rmesh)

def intersect_check(xa,ya, xb,yb, xk,yk, xl, yl):
    if xa == xb and xk == xl:
        return False
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
            bk = ya-Ak*xk
            Yx = Ak*xa + bk
            if Yx < max(ya,yb) and Yx > min(ya,yb):
                return [xa, Yx]
            else:
                return False

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

    if (Xx < max(I1[0],I2[0])) and (Xx > min(I1[1], I2[1])):
        return False
    else:
        return[Xx, Yx]

if __name__ == '__main__':

    def aff_face(faces):
        fig = pyplot.figure()
        ar = fig.add_subplot(1,1,1)
        for face in faces:
            if type(face) is MeshFace:
                edge = face.get_edge()
                tx = []
                ty = []
                for iter in [0,1,2,3]:
                    tx.append(edge.get_dest().get_pos()[0])
                    ty.append(edge.get_dest().get_pos()[1])

                    edge = edge.get_next()
            ar.plot(tx,ty)
        pyplot.show()


    a = MeshVertex(1,2)
    b = MeshVertex(1.6,2.3)
    c = MeshVertex(1,10)
    print(clockwise_check(a,b,c))

    # tri = create_isolated_face(a,b,c)
    # aff_face([tri['faces'][0]])
    #
    #CAN HANDLE 3 * 2^n VERTICES
    #TODO: add signle point handling to hendle any number of points
    # points_list=[[1,1],[2,3],[3,4],[4,7],[5,6],[6,3]]
    # points_list=[[1.,2.],[2.,1.],[2.2,8.],[5.,3.],[5.4,2.],[5.6,5.]]
    points_list=[[1.,2.],[1.6,2.3],[2.,1.],[2.2,8.],[2.7,1.8],[3.5,5.],[4.2,8],[4.6,3.6],[5.,3.],[5.4,2.],[5.6,5.],[6.,4.2]]
    # points_list=[[4.2,4.8],[4.6,3.6],[5.,3.],[5.4,2.],[5.6,5.],[6.,4.2]]
    # points_list=[[1.,2],[2.,2.],[3.,8.],[4.,4.],[5.,1.],[6.,6.]]
    # points_list=[[1.2,6.],[1.6,4.8],[2.,1.],[2.5,8.],[2.7,1.8]]

    # points_list=[[1.,2.],[1.,4.],[1.,6.],[3.5,1.],[3.5,3.], [3.5,5.]]
    # points_list=[[1.2,2.],[1.2,5.],[1.2,7.]]

    # map=divide_and_conquer(points_list)
    # aff_face(map['faces'])

    translated_map = []
    for i in list(range(2)):
        for j in list(range(2)):
            for pnum in list(range(len(points_list))):
                translated_map.append([points_list[pnum][0]+i*8,points_list[pnum][1]+j*8])
    map = []
    for i in list(range(4)):
        for j in list(range(4)):
            for pnum in list(range(len(translated_map))):
                map.append([translated_map[pnum][0]+(i)*20,translated_map[pnum][1]+(j)*20])

    translated_map = order_points_list(map)
    print(len(translated_map))

    print(datetime.now())
    map=divide_and_conquer(translated_map)
    print(datetime.now())

    #
    aff_face(map['faces'])
