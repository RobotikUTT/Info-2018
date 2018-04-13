from matplotlib import pyplot
from numpy.linalg import det
from numpy import array, concatenate, append
from math import pow
from datetime import datetime
import os


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

    def set_vertex(self, v):
        self._vertex = v

    def get_vertex(self):
        return self._vertex

class MeshEdge():
    def __init__(self, vertex_dest, face=None,  opp_edge=None, next_edge=None, is_obstacle=False):
        self._vertex_dest = vertex_dest
        self._face = face
        self._opp_edge = opp_edge
        self._next_edge= next_edge
        self._is_obstacle = is_obstacle

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

    def set_obstacle(self, is_obstacle):
        self._is_obstacle = is_obstacle

    def is_obstacle(self):
        return self._is_obstacle

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
            bk = yk-Ak*xk
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

    if (Xx < max(I1[0],I2[0])) or (Xx > min(I1[1], I2[1])):
        return False
    else:
        return[Xx, Yx]

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

def insert_point(point, sorted_list):
    insertion_index = len(sorted_list)
    offset = 0
    index_found = False
    while not index_found:
        if insertion_index > 0:
            insertion_index = insertion_index //2
            if point[0] >= sorted_list[insertion_index+offset][0] and point[1] >= sorted_list[insertion_index+offset][1]:
                if sorted_list[insertion_index+offset+1][0] >= point[0] and sorted_list[insertion_index+offset+1][1] >= point[0]:
                    index_found = True
                else:
                    offset += insertion_index
        else:
            index_found = True
    sorted_list.insert(insertion_index+offset+1, point)
    return sorted_list

def concatenate_sorted(list_of_lists):
    sorted_list = list_of_lists[0]
    for points_list in list_of_lists[1:]:
        insertion_index = 0
        for point in points_list:
            added = False
            while not added:
                if insertion_index < len(sorted_list):
                    if point[0] < sorted_list[insertion_index][0]:
                        insertion_index += 1
                    elif point[1] < sorted_list[insertion_index][1]:
                        insertion_index += 1
                    else:
                        sorted_list.insert(insertion_index, point)
                        added = True
                else:
                    sorted_list.append(point)
                    added = True

    return sorted_list

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

def cone_check(center_pos, b1_pos, b2_pos, x_pos):


    # fig = pyplot.figure()
    # ar = fig.add_subplot(1,1,1)
    # tx = []
    # ty = []
    # rx = []
    # ry = []
    # tx.append(center_pos[0])
    # tx.append(b1_pos[0])
    # tx.append(center_pos[0])
    # tx.append(b2_pos[0])
    # rx.append(center_pos[0])
    # rx.append(x_pos[0])
    # ty.append(center_pos[1])
    # ty.append(b1_pos[1])
    # ty.append(center_pos[1])
    # ty.append(b2_pos[1])
    # ry.append(center_pos[1])
    # ry.append(x_pos[1])
    # ar.plot(tx,ty, color='k')
    # ar.plot(rx,ry, color='r')
    # pyplot.show()

    if (b1_pos[0]-center_pos[0])*(x_pos[1]-center_pos[1]) - (b1_pos[1]-center_pos[1])*(x_pos[0]-center_pos[0]) >= 0:
        return False
    if (b2_pos[0]-center_pos[0])*(x_pos[1]-center_pos[1]) - (b2_pos[1]-center_pos[1])*(x_pos[0]-center_pos[0])  <= 0:
        return False
    else:
        return True

def half_plane_check(a_pos, b_pos, x_pos):
    if (b_pos[0]-a_pos[0])*(x_pos[1]-a_pos[1]) - (b_pos[1]-a_pos[1])*(x_pos[0]-a_pos[0]) >= 0:
        return True
    else:
        return False

def inside_triangle_check(a_pos, b_pos, c_pos, x_pos):

    ab_side = (b_pos[0]-a_pos[0])*(x_pos[1]-a_pos[1]) - (b_pos[1]-a_pos[1])*(x_pos[0]-a_pos[0]) >= 0
    if ((c_pos[0]-a_pos[0])*(x_pos[1]-a_pos[1]) - (c_pos[1]-a_pos[1])*(x_pos[0]-a_pos[0])  >= 0) == ab_side:
        return False
    elif ((c_pos[0]-b_pos[0])*(x_pos[1]-b_pos[1]) - (c_pos[1]-b_pos[1])*(x_pos[0]-b_pos[0])  >= 0) != ab_side:
        return False
    else:
        return True

def diamond_angle(x,y):

    if y >= 0:
        if x>=0:
            res = y/(x+y)
        else:
            res = 1-x/(-x+y)
    else:
        if x < 0:
            res = 2-y/(-x-y)
        else:
            res = 3+x/(x-y)

    print('dia')
    print(x,y)
    print(2-res)

    return res-2

def diamond_diff(alpha, beta):
    if alpha >= 0:
        if beta >= 0:
            res = abs(alpha - beta)
        else:
            if beta >= alpha - 2:
                res = abs(alpha - beta)
            else:
                res = abs((2- alpha)-(-2 - beta))
    else:
        if beta < 0:
            res = abs(alpha - beta)
        else:
            if alpha >= beta -2:
                res = abs(beta - alpha)
            else:
                res = abs((2 - beta)-(-2 - alpha))

    print('diff')
    print(alpha, beta, res)

    return res

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

def edge_flip(edge):
    if type(edge.get_face()) is MeshFace and type(edge.get_opp().get_face()) is MeshFace:
        edge.get_next().get_next().get_face().set_edge(edge.get_next().get_next())
        edge.get_opp().get_next().get_next().get_face().set_edge(edge.get_opp().get_next().get_next())
        edge.get_next().set_face(edge.get_opp().get_next().get_next().get_face())
        edge.get_opp().get_next().set_face(edge.get_next().get_next().get_face())
        edge.get_opp().get_dest().set_edge(edge.get_opp().get_next())
        edge.get_dest().set_edge(edge.get_next())

        edge.get_next().get_next().set_next(edge.get_opp().get_next())
        edge.get_opp().get_next().get_next().set_next(edge.get_next())

        edge.get_opp().set_dest(edge.get_next().get_next().get_next().get_dest())
        edge.get_opp().set_next(edge.get_next().get_next().get_next().get_next())


        edge.get_next().get_next().get_next().set_next(edge)
        edge.set_dest(edge.get_opp().get_next().get_next().get_dest())
        edge.set_next(edge.get_opp().get_next().get_next().get_next())
        edge.get_opp().get_next().get_next().set_next(edge.get_opp())


def merge_edges(start_edge, faces, vertices):

    #######
    #Identify vertices around the starting edge
    startL = start_edge.get_dest()
    startR = start_edge.get_opp().get_dest()

    #initiate the duo of vertices candidates, starting with the most exterior edge (minimal angle)
    #
    edgeL = start_edge.get_next().get_opp().get_next().get_opp()

    if edgeL.get_opp().get_dest() == startR:
        Lcandidate = startL
        Lfound = True
    else:
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

    Rfound = False

    if edgeR.get_dest() == startL:
        Rcandidate = startR
        Rfound = True
    else:
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
            if Rfound:
                elected_vertex = Rcandidate
                elected_edge = edgeR

    if elected_vertex is None:
        aff_faces(faces)
        pyplot.show()
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
                        current_edge.get_opp().get_dest().set_edge(current_edge.get_next().get_next().get_opp())

                        if type(current_edge.get_next().get_opp()) is MeshFace():
                            current_edge.get_next().get_dest().set_edge(current_edge.get_next().get_opp())
                        else:
                            current_edge.get_next().get_dest().set_edge(current_edge.get_next().get_next())
                        current_edge.get_next().get_next().set_next(current_edge.get_opp().get_next())
                        # current_edge.get_next().get_next().get_dest().set_edge(current_edge.get_next().get_next().get_opp())

                        current_edge.get_face().set_edge(None)

                        current_edge.get_next().get_next().set_face(current_edge.get_opp().get_face())
                        current_edge.get_next().set_next(None)
                        current_edge.get_next().set_face(None)
                        del faces[faces.index(current_edge.get_face())]
                        current_edge.set_face(None)
                        new_ghost = add_ghost_vertex(current_edge.get_next())
                        faces.append(new_ghost['face'])
                        vertices.append(new_ghost['vertex'])
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
            del faces[faces.index(elected_edge.get_opp().get_face())]
            #set the face pointer to new face
            elected_edge.get_opp().set_face(new_face)
            #BUG HERE, Face stays in faces list but has no edge (none)
            elected_edge.get_opp().get_next().set_face(start_edge.get_face())
            #delete refeences of the left hand ghost edge to remove
            elected_edge.get_opp().get_next().get_next().set_next(None)
            elected_edge.get_opp().get_next().get_next().set_opp(None)
            elected_edge.get_opp().get_next().get_next().set_face(None)
            elected_edge.get_opp().get_next().get_next().set_vertex(None)
            #set the next attribute to link both gohst triangles
            elected_edge.get_opp().get_next().set_next(start_edge.get_next().get_next())

            # if start_edge.get_next().get_vertex() in vertices:
            del vertices[vertices.index(elected_edge.get_opp().get_next().get_vertex())]
            elected_edge.get_opp().get_next().set_vertex(start_edge.get_next().get_vertex())

            #close the new merged ghost triangle
            elected_edge.get_opp().get_next().get_next().set_next(new_opp_edge)
            new_opp_edge.set_next(elected_edge.get_opp().get_next())
            new_opp_edge.set_face(start_edge.get_face())
            # new_opp_edge.get_next().set_vertex(start_edge.get_next().get_vertex())
            new_opp_edge.get_face().set_edge(new_opp_edge)

            #delete all references ofthe righ hand ghost edge to merge
            start_edge.get_next().set_face(None)
            start_edge.get_next().set_next(None)
            start_edge.get_next().set_opp(None)
            start_edge.get_next().set_vertex(None)

            start_edge.set_next(elected_edge.get_opp())
            elected_edge.get_opp().set_next(new_edge)

            start_edge.set_face(new_face)
            new_edge.set_face(new_face)

            #set the vertices so that they point along the edges of the created hull (clockwise)

            if type(elected_edge.get_face()) is MeshFace:
                elected_vertex.set_edge(new_edge)
            else:
                startL.set_edge(elected_edge.get_opp())


            faces.insert(0,new_face)
            merged = {'faces': faces, 'edge': new_edge, 'vertices':vertices}
            return merged

        if Rfound:
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
                        # current_edge.get_dest().set_edge(current_edge.get_next())
                        if type(current_edge.get_next().get_opp().get_face()) is MeshFace:
                            current_edge.get_next().get_dest().set_edge(current_edge.get_next().get_opp())
                        else:
                            current_edge.get_next().get_dest().set_edge(current_edge.get_next().get_next())

                        if type(current_edge.get_next().get_next().get_opp().get_face()) is MeshFace:
                            startR.set_edge(current_edge.get_next().get_next().get_opp())

                        new_ghost = add_ghost_vertex(current_edge.get_next().get_next())
                        faces.append(new_ghost['face'])
                        vertices.append(new_ghost['vertex'])
                        current_edge.get_next().get_next().get_next().get_next().set_opp(current_edge.get_opp().get_next())
                        current_edge.get_next().get_next().get_next().set_opp(start_edge.get_next().get_next())
                        current_edge.get_next().set_next(current_edge.get_opp().get_next())
                        current_edge.get_next().set_face(current_edge.get_opp().get_face())
                        current_edge.get_opp().get_next().get_next().set_next(current_edge.get_next())

                        current_edge.get_face().set_edge(None)
                        del faces[faces.index(current_edge.get_face())]

                        current_edge.get_opp().set_next(None)
                        current_edge.get_opp().set_opp(None)
                        current_edge.get_opp().set_face(None)
                        current_edge.get_opp().set_dest(None)

                        current_edge.set_face(None)
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
            del faces[faces.index(start_edge.get_face())]

            start_edge.get_next().set_face(elected_edge.get_opp().get_face())

            #delete references to left side of linking ghost face

            start_edge.get_next().get_next().set_face(None)
            start_edge.get_next().get_next().set_opp(None)
            start_edge.get_next().get_next().set_next(None)
            start_edge.get_next().get_next().set_vertex(None)

            #link both ghost by a next pointer
            start_edge.get_next().set_next(elected_edge.get_opp().get_next().get_next())
            del vertices[vertices.index(elected_edge.get_opp().get_next().get_vertex())]
            elected_edge.get_opp().get_next().set_vertex(None)

            #then delete right hand side of the ghodt edge to remove
            elected_edge.get_opp().get_next().set_face(None)
            elected_edge.get_opp().get_next().set_opp(None)
            elected_edge.get_opp().get_next().set_next(None)
            elected_edge.get_opp().get_next().set_vertex(None)

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

            if type(elected_edge.get_face()) is not MeshFace:
                # elected_vertex.set_edge(new_opp_edge.get_next().get_next().get_opp().get_next().get_next().get_opp())
                elected_vertex.set_edge(elected_edge.get_opp())
            # else:
            #     elected_vertex.set_edge(elected_edge.get_opp())

            # if type(startL.get_edge().get_face()) is MeshFace:
            #     startL.set_edge(new_edge)
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
            merged = {'faces': faces, 'edge': new_edge, 'vertices':vertices}
            return merged

def merge_mesh(meshL, meshR):

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
    #TODO: check for occlusion between each upper or lower vertices (if there is, select another starting or exit vertex!)
    upper_L_vertex = meshL['vertices'][0]
    upper_L_height = meshL['vertices'][0].get_pos()[1]
    for vertex in meshL['vertices']:
        if type(vertex) is MeshVertex:
            if vertex.get_pos()[1] > upper_L_height:
                upper_L_vertex = vertex
                upper_L_height = vertex.get_pos()[1]
            elif vertex.get_pos()[1] == upper_L_height and vertex.get_pos()[0] < upper_L_vertex.get_pos()[0]:
                upper_L_vertex = vertex

    lower_L_vertex = meshL['vertices'][0]
    lower_L_height = meshL['vertices'][0].get_pos()[1]
    for vertex in meshL['vertices']:
        if type(vertex) is MeshVertex:
            if vertex.get_pos()[1] < lower_L_height:
                lower_L_vertex = vertex
                lower_L_height = vertex.get_pos()[1]
            elif vertex.get_pos()[1] == lower_L_height and vertex.get_pos()[0] < lower_L_vertex.get_pos()[0]:
                lower_L_vertex = vertex

    upper_R_vertex = meshR['vertices'][0]
    upper_R_height = meshR['vertices'][0].get_pos()[1]
    for vertex in meshR['vertices']:
        if type(vertex) is MeshVertex:
            if vertex.get_pos()[1] > upper_R_height:
                upper_R_vertex = vertex
                upper_R_height = vertex.get_pos()[1]
            elif vertex.get_pos()[1] == upper_R_height and vertex.get_pos()[0] > upper_R_vertex.get_pos()[0]:
                upper_R_vertex = vertex

    lower_R_vertex = meshR['vertices'][0]
    lower_R_height = meshR['vertices'][0].get_pos()[1]
    for vertex in meshR['vertices']:
        if type(vertex) is MeshVertex:
            if vertex.get_pos()[1] < lower_R_height:
                lower_R_vertex = vertex
                lower_R_height = vertex.get_pos()[1]
            elif vertex.get_pos()[1] == lower_R_height and vertex.get_pos()[0] > lower_R_vertex.get_pos()[0]:
                lower_R_vertex = vertex

    if L_flat and R_flat and lower_L_vertex.get_pos()[0] == upper_L_vertex.get_pos()[0] \
            and lower_R_vertex.get_pos()[0] == upper_R_vertex.get_pos()[0]\
            and upper_L_vertex.get_pos()[0] == lower_R_vertex.get_pos()[0]:

        upper_L_edge = upper_L_vertex.get_edge()
        lower_R_edge = lower_R_vertex.get_edge()
        ghosts = create_isolated_edge(upper_L_vertex, lower_R_vertex)

        upper_L_vertex.get_edge().get_next().set_opp(lower_R_edge.get_next().get_next())
        lower_R_edge.get_next().get_next().set_opp(upper_L_vertex.get_edge().get_next())
        lower_R_vertex.get_edge().get_next().get_next().set_opp(lower_R_edge.get_opp().get_next())
        lower_R_edge.get_opp().get_next().set_opp(lower_R_vertex.get_edge().get_next().get_next())

        upper_L_vertex.get_edge().get_next().get_next().set_opp(upper_L_edge.get_opp().get_next())
        upper_L_edge.get_opp().get_next().set_opp(upper_L_vertex.get_edge().get_next().get_next())
        lower_R_vertex.get_edge().get_next().set_opp(upper_L_edge.get_next().get_next())
        upper_L_edge.get_next().get_next().set_opp(lower_R_vertex.get_edge().get_next())

        lower_R_vertex.set_edge(lower_R_edge)

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

        no_occlusion = False
        while not no_occlusion:

            # segz = []
            # for f in faces:
            #     if type(f) is MeshFace:
            #         e = f.get_edge()
            #         for i in [0,1,2,3]:
            #             segz.append([e.get_dest().get_pos(), e.get_next().get_dest().get_pos()])
            #             e = e.get_next()
            # fig = pyplot.figure()
            # ar = fig.add_subplot(1,1,1)
            # sar = array(segz)
            # for h in sar:
            #     tx = h[:,0]
            #     ty = h[:,1]
            #     ar.plot(tx,ty, color='r')
            # ar.plot(array([start_L.get_pos()[0],start_R.get_pos()[0]]),array([start_L.get_pos()[1],start_R.get_pos()[1]]), color='k')
            # ar.plot(array([stop_L.get_pos()[0],stop_R.get_pos()[0]]),array([stop_L.get_pos()[1],stop_R.get_pos()[1]]), color='b')
            # pyplot.show()

            no_occlusion = True
            if not R_flat:
                if start_L.get_pos()[1] >= start_R.get_pos()[1]:
                    if half_plane_check(start_R.get_edge().get_dest().get_pos(),start_R.get_pos(), start_L.get_pos()) \
                            and half_plane_check(start_R.get_pos(), start_R.get_edge().get_opp().get_next().get_opp().get_next().get_dest().get_pos(), start_L.get_pos()):
                        no_occlusion = False
                        start_R = start_R.get_edge().get_dest()

                if stop_L.get_pos()[1] <= stop_R.get_pos()[1]:
                    if half_plane_check(stop_R.get_edge().get_dest().get_pos(), stop_R.get_pos(),stop_L.get_pos()) \
                            and half_plane_check(stop_R.get_pos(), stop_R.get_edge().get_opp().get_next().get_opp().get_next().get_dest().get_pos(),stop_L.get_pos()):
                        no_occlusion = False
                        stop_R = stop_R.get_edge().get_opp().get_next().get_opp().get_next().get_dest()
            if not L_flat:
                if start_L.get_pos()[1] <= start_R.get_pos()[1]:
                    if half_plane_check(start_L.get_edge().get_dest().get_pos(), start_L.get_pos(), start_R.get_pos()) \
                            and half_plane_check(start_L.get_pos(), start_L.get_edge().get_opp().get_next().get_opp().get_next().get_dest().get_pos(),
                                                start_R.get_pos()):
                        no_occlusion = False
                        start_L = start_L.get_edge().get_opp().get_next().get_opp().get_next().get_dest()
                if stop_L.get_pos()[1] >= stop_R.get_pos()[1]:
                    if half_plane_check(stop_L.get_pos(), stop_L.get_edge().get_opp().get_next().get_opp().get_next().get_dest().get_pos(), stop_R.get_pos()) \
                            and half_plane_check(stop_L.get_edge().get_dest().get_pos(), stop_L.get_pos(), stop_R.get_pos()):
                        no_occlusion = False
                        stop_L = stop_L.get_edge().get_dest()

        while no_occlusion:

            # segz = []
            # for f in faces:
            #     if type(f) is MeshFace:
            #         e = f.get_edge()
            #         for i in [0,1,2,3]:
            #             segz.append([e.get_dest().get_pos(), e.get_next().get_dest().get_pos()])
            #             e = e.get_next()
            # fig = pyplot.figure()
            # ar = fig.add_subplot(1,1,1)
            # sar = array(segz)
            # for h in sar:
            #     tx = h[:,0]
            #     ty = h[:,1]
            #     ar.plot(tx,ty, color='r')
            # ar.plot(array([start_L.get_pos()[0],start_R.get_pos()[0]]),array([start_L.get_pos()[1],start_R.get_pos()[1]]), color='k')
            # ar.plot(array([stop_L.get_pos()[0],stop_R.get_pos()[0]]),array([stop_L.get_pos()[1],stop_R.get_pos()[1]]), color='b')
            # pyplot.show()

            no_occlusion = False
            if not R_flat:
                if start_L.get_pos()[1] < start_R.get_pos()[1]:
                    if not half_plane_check(start_R.get_pos(), start_R.get_edge().get_opp().get_next().get_opp().get_next().get_dest().get_pos(), start_L.get_pos()):
                            # and not half_plane_check(start_R.get_edge().get_opp().get_next().get_opp().get_next().get_dest().get_pos(), start_R.get_edge().get_opp().get_next().get_opp().get_next().get_next().get_opp().get_next().get_dest().get_pos(), start_L.get_pos()):
                        no_occlusion = True
                        start_R = start_R.get_edge().get_opp().get_next().get_opp().get_next().get_dest()
                if stop_L.get_pos()[1] > stop_R.get_pos()[1]:
                    if not half_plane_check(stop_R.get_edge().get_dest().get_pos(), stop_R.get_pos(),stop_L.get_pos()):
                            # and not half_plane_check(stop_R.get_edge().get_opp().get_next().get_next().get_opp().get_next().get_next().get_opp().get_dest().get_pos(), stop_R.get_edge().get_dest().get_pos(),stop_L.get_pos()):
                        no_occlusion = True
                        stop_R = stop_R.get_edge().get_dest()
            if not L_flat:
                if start_L.get_pos()[1] > start_R.get_pos()[1]:
                    if not half_plane_check(start_L.get_edge().get_dest().get_pos(), start_L.get_pos(), start_R.get_pos()):
                            # and not half_plane_check(start_L.get_edge().get_opp().get_next().get_next().get_opp().get_next().get_next().get_opp().get_dest().get_pos(), start_L.get_edge().get_dest().get_pos(), start_R.get_pos()):
                        no_occlusion = True
                        start_L = start_L.get_edge().get_dest()
                if stop_L.get_pos()[1] < stop_R.get_pos()[1]:
                    if not half_plane_check(stop_L.get_pos(), stop_L.get_edge().get_opp().get_next().get_opp().get_next().get_dest().get_pos(), stop_R.get_pos()):
                            # and not half_plane_check(stop_L.get_edge().get_opp().get_next().get_opp().get_next().get_dest().get_pos(),stop_L.get_edge().get_opp().get_next().get_opp().get_next().get_next().get_opp().get_next().get_dest().get_pos() ,stop_R.get_pos()):
                        no_occlusion = True
                        stop_L = stop_L.get_edge().get_opp().get_next().get_opp().get_next().get_dest()

        # segz = []
        # for f in faces:
        #     if type(f) is MeshFace:
        #         e = f.get_edge()
        #         for i in [0,1,2,3]:
        #             segz.append([e.get_dest().get_pos(), e.get_next().get_dest().get_pos()])
        #             e = e.get_next()
        # fig = pyplot.figure()
        # ar = fig.add_subplot(1,1,1)
        # sar = array(segz)
        # for h in sar:
        #     tx = h[:,0]
        #     ty = h[:,1]
        #     ar.plot(tx,ty, color='r')
        # ar.plot(array([start_L.get_pos()[0],start_R.get_pos()[0]]),array([start_L.get_pos()[1],start_R.get_pos()[1]]), color='k')
        # ar.plot(array([stop_L.get_pos()[0],stop_R.get_pos()[0]]),array([stop_L.get_pos()[1],stop_R.get_pos()[1]]), color='b')
        # pyplot.show()




        #create an edge between both lower points and its opposite
        merge_edge_start = MeshEdge(start_L)
        merge_edge_opp = MeshEdge(start_R,opp_edge=merge_edge_start)

        merge_edge_start.set_opp(merge_edge_opp)

        #add the ghost faces under and over them
        over_ghost_face = add_ghost_vertex(merge_edge_start)
        under_ghost_face = add_ghost_vertex(merge_edge_opp)

        faces.append(over_ghost_face['face'])
        faces.append(under_ghost_face['face'])
        vertices.append(over_ghost_face['vertex'])
        vertices.append(under_ghost_face['vertex'])

        #link the opposites edges of the ghost faces bellow and over
        if start_L != stop_L:
            merge_edge_start.get_next().set_opp(start_L.get_edge().get_opp().get_next().get_opp())
            start_L.get_edge().get_opp().get_next().get_opp().set_opp(merge_edge_start.get_next())
            merge_edge_opp.get_next().get_next().set_opp(start_L.get_edge().get_opp().get_next())
            start_L.get_edge().get_opp().get_next().set_opp(merge_edge_opp.get_next().get_next())
        else:
            merge_edge_start.get_next().set_opp(merge_edge_start.get_opp().get_next().get_next())
            merge_edge_start.get_opp().get_next().get_next().set_opp(merge_edge_start.get_next())

        if start_R != stop_R:
            merge_edge_opp.get_next().set_opp(start_R.get_edge().get_opp().get_next().get_opp())
            start_R.get_edge().get_opp().get_next().get_opp().set_opp(merge_edge_opp.get_next())
            merge_edge_start.get_next().get_next().set_opp(start_R.get_edge().get_opp().get_next())
            start_R.get_edge().get_opp().get_next().set_opp(merge_edge_start.get_next().get_next())
        else:
            merge_edge_start.get_next().get_next().set_opp(merge_edge_start.get_opp().get_next())
            merge_edge_start.get_opp().get_next().set_opp(merge_edge_start.get_next().get_next())

        current_edge = merge_edge_start

        segz = []
        for f in faces:
            if type(f) is MeshFace:
                e = f.get_edge()
                for i in [0,1,2,3]:
                    segz.append([e.get_dest().get_pos(), e.get_next().get_dest().get_pos()])
                    e = e.get_next()

        # fig = pyplot.figure()
        # ar = fig.add_subplot(1,1,1)
        # sar = array(segz)
        # for h in sar:
        #     tx = h[:,0]
        #     ty = h[:,1]
        #     ar.plot(tx,ty, color='r')
        # ar.plot(array([start_L.get_pos()[0],start_R.get_pos()[0]]),array([start_L.get_pos()[1],start_R.get_pos()[1]]), color='k')
        # ar.plot(array([stop_L.get_pos()[0],stop_R.get_pos()[0]]),array([stop_L.get_pos()[1],stop_R.get_pos()[1]]), color='b')
        # pyplot.show()

        while not (current_edge.get_dest() == stop_L and current_edge.get_opp().get_dest() == stop_R):

            merged = merge_edges(current_edge, faces, vertices)
            # if merged is None:
            #     print('problem?')
            #     segz = []
            #     for fnum in list(range(len(faces))):
            #         if type(faces[fnum]) is MeshFace:
            #             e = faces[fnum].get_edge()
            #             for i in [0,1,2,3]:
            #                 if type(e) is GhostEdge:
            #                     print('ge')
            #                     print(fnum)
            #                     print(i)
            #                 segz.append([e.get_dest().get_pos(), e.get_next().get_dest().get_pos()])
            #                 e = e.get_next()
            #     fig = pyplot.figure()
            #     ar = fig.add_subplot(1,1,1)
            #     sar = array(segz)
            #     for h in sar:
            #         tx = h[:,0]
            #         ty = h[:,1]
            #         ar.plot(tx,ty)
            #     pyplot.show()

            faces = merged['faces']
            vertices = merged['vertices']

            if type(current_edge.get_opp().get_face()) is MeshFace:
                if not circumcircle_check(current_edge.get_dest(), current_edge.get_next().get_dest(), current_edge.get_opp().get_dest(), current_edge.get_opp().get_next().get_dest()):
                    if not colinear_check(current_edge.get_opp().get_next().get_dest(), current_edge.get_dest(), current_edge.get_next().get_dest()) \
                            and not colinear_check(current_edge.get_opp().get_next().get_dest(), current_edge.get_opp().get_dest(), current_edge.get_next().get_dest()):
                        edge_flip(current_edge)

            if type(current_edge.get_next().get_opp().get_face()) is MeshFace:
                if not circumcircle_check(current_edge.get_dest(), current_edge.get_next().get_dest(), current_edge.get_opp().get_dest(), current_edge.get_next().get_opp().get_next().get_dest()):
                    if not colinear_check(current_edge.get_next().get_opp().get_next().get_dest(), current_edge.get_dest(), current_edge.get_opp().get_dest()) \
                            and not colinear_check(current_edge.get_next().get_opp().get_next().get_dest(), current_edge.get_next().get_dest(), current_edge.get_opp().get_dest()):
                        edge_flip(current_edge.get_next())

            if type(current_edge.get_next().get_next().get_opp().get_face()) is MeshFace:
                if not circumcircle_check(current_edge.get_dest(), current_edge.get_next().get_dest(), current_edge.get_opp().get_dest(), current_edge.get_next().get_next().get_opp().get_next().get_dest()):
                    if not colinear_check(current_edge.get_next().get_next().get_opp().get_next().get_dest(), current_edge.get_dest(), current_edge.get_next().get_dest()) \
                            and not colinear_check(current_edge.get_next().get_next().get_opp().get_next().get_dest(), current_edge.get_dest(), current_edge.get_opp().get_dest()):
                        edge_flip(current_edge.get_next().get_next())

            current_edge = merged['edge'].get_opp()

        start_L.set_edge(merge_edge_start.get_opp().get_next().get_next().get_opp().get_next().get_next().get_opp())
        start_R.set_edge(merge_edge_start)
        current_edge.get_dest().set_edge(current_edge.get_opp())

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
    if len(ordered_points) == 1:
        print('bra')
        vertices = []
        faces = []
        a = MeshVertex(ordered_points[0][0],ordered_points[0][1])
        vertices.append(a)
        return {'vertices': vertices, 'faces':faces}

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
        vertices.append(tri['vertices'][1])
        vertices.append(tri['vertices'][2])

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

        merged = merge_mesh(Lmesh,Rmesh)
        return merged

def grid_sampling(square_len, sorted_vertices):
    sampled_vertices = [[None]]
    for vertex in sorted_vertices:
        if type(vertex) is MeshVertex:
            i = int(vertex.get_pos()[0]//square_len)
            j = int(vertex.get_pos()[1]//square_len)
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
    return sampled_vertices

def fill_grid(sampled_vertices, square_len):

    grid_filled = False
    while not grid_filled:
        none_remaining = False
        for i in list(range(len(sampled_vertices))):
            for j in list(range(len(sampled_vertices[i]))):
                if sampled_vertices[i][j] is None:
                    closest_sample = None
                    if i == 0:
                        ibound = [+1]
                    elif  i == len(sampled_vertices) -1:
                        ibound = [-1]
                    else:
                        ibound = [+1, -1]
                    if j == 0:
                        jbound = [+1]
                    elif  j == len(sampled_vertices[i]) -1:
                        jbound = [-1]
                    else:
                        jbound = [+1, -1]

                    for di in ibound:
                        for dj in jbound:
                            if sampled_vertices[i+di][j+dj] is not None:
                                if closest_sample is None:
                                    closest_sample = sampled_vertices[i+di][j+dj]
                                else:
                                    if pow(sampled_vertices[i+di][j+dj].get_pos()[0]-square_len*(i+1/2), 2) + pow(sampled_vertices[i+di][j+dj].get_pos()[1]-square_len*(j+1/2), 2) \
                                            < pow(closest_sample.get_pos()[0]-square_len*(i+1/2),2)+pow(closest_sample.get_pos()[1]-square_len*(j+1/2),2):
                                       closest_sample = sampled_vertices[i+di][j+dj]
                    if closest_sample is None:
                        none_remaining = True
                    else:
                        sampled_vertices[i][j] = closest_sample
        if not none_remaining:
            grid_filled = True
    return sampled_vertices

def find_face(x_pos, start_vertex, face=None):
    if face is None:
        face = start_vertex.get_edge().get_face()
        return find_face(x_pos, start_vertex, face)
    else:
        #TODO: check if the x_pos is a vertex position or on an edge!!
        if type(face) is GhostFace:
            border_edge = face.get_edge()
            if border_edge.get_opp().get_dest() == start_vertex:
                if not cone_check(start_vertex.get_pos(), border_edge.get_dest().get_pos(),
                                  border_edge.get_next().get_next().get_opp().get_next().get_next().get_opp().get_dest().get_pos() , x_pos):
                    return find_face(x_pos, start_vertex, border_edge.get_next().get_next().get_opp().get_next().get_next().get_opp().get_face())
                else:
                    return face
            elif border_edge.get_dest() == start_vertex:
                if not cone_check(start_vertex.get_pos(), border_edge.get_next().get_opp().get_next().get_opp().get_dest().get_pos(),
                                  border_edge.get_opp().get_dest().get_pos() , x_pos):
                    return find_face(x_pos, start_vertex, border_edge.get_opp().get_face())
                else:
                    return face
            else:
                print('oopsys')
        else:
            crossing_edge = face.get_edge()
            if crossing_edge.get_dest() == start_vertex:
                crossing_edge = crossing_edge.get_next().get_next()
            elif crossing_edge.get_next().get_next().get_dest() == start_vertex:
                crossing_edge = crossing_edge.get_next()

            if cone_check(start_vertex.get_pos(), crossing_edge.get_next().get_next().get_dest().get_pos(), crossing_edge.get_dest().get_pos(), x_pos):
                # fig = pyplot.figure()
                # ar = fig.add_subplot(1,1,1)
                # tx = []
                # ty = []
                # tx.append(start_vertex.get_pos()[0])
                # tx.append(crossing_edge.get_next().get_next().get_dest().get_pos()[0])
                # tx.append(crossing_edge.get_dest().get_pos()[0])
                # tx.append(start_vertex.get_pos()[0])
                # tx.append(x_pos[0])
                # ty.append(start_vertex.get_pos()[1])
                # ty.append(crossing_edge.get_next().get_next().get_dest().get_pos()[1])
                # ty.append(crossing_edge.get_dest().get_pos()[1])
                # ty.append(start_vertex.get_pos()[1])
                # ty.append(x_pos[1])
                # ar.plot(tx,ty)
                # pyplot.show()

                if intersect_check(start_vertex.get_pos()[0], start_vertex.get_pos()[1], x_pos[0], x_pos[1],
                                   crossing_edge.get_opp().get_dest().get_pos()[0],
                                   crossing_edge.get_opp().get_dest().get_pos()[1],
                                   crossing_edge.get_dest().get_pos()[0],
                                   crossing_edge.get_dest().get_pos()[1]):
                    # fig = pyplot.figure()
                    # ar = fig.add_subplot(1,1,1)
                    # tx = []
                    # ty = []
                    # rx = []
                    # ry = []
                    # sx = []
                    # sy = []
                    # tx.append(x_pos[0])
                    # tx.append(start_vertex.get_pos()[0])
                    # sx.append(start_vertex.get_pos()[0])
                    # sx.append(intersect_check(start_vertex.get_pos()[0], start_vertex.get_pos()[1], x_pos[0], x_pos[1],
                    #                           crossing_edge.get_opp().get_dest().get_pos()[0],
                    #                           crossing_edge.get_opp().get_dest().get_pos()[1],
                    #                           crossing_edge.get_dest().get_pos()[0],
                    #                           crossing_edge.get_dest().get_pos()[1])[0])
                    # rx.append(crossing_edge.get_dest().get_pos()[0])
                    # rx.append(crossing_edge.get_opp().get_dest().get_pos()[0])
                    # ty.append(x_pos[1])
                    # ty.append(start_vertex.get_pos()[1])
                    # sy.append(start_vertex.get_pos()[1])
                    # sy.append(intersect_check(start_vertex.get_pos()[0], start_vertex.get_pos()[1], x_pos[0], x_pos[1],
                    #                           crossing_edge.get_opp().get_dest().get_pos()[0],
                    #                           crossing_edge.get_opp().get_dest().get_pos()[1],
                    #                           crossing_edge.get_dest().get_pos()[0],
                    #                           crossing_edge.get_dest().get_pos()[1])[1])
                    # ry.append(crossing_edge.get_dest().get_pos()[1])
                    # ry.append(crossing_edge.get_opp().get_dest().get_pos()[1])
                    # ar.plot(tx,ty, color='b')
                    # ar.plot(rx,ry, color='r')
                    # ar.plot(sx,sy, color='k')
                    # pyplot.show()

                    if type(crossing_edge.get_opp().get_face()) is MeshFace:
                        return find_face(x_pos, crossing_edge.get_opp().get_next().get_dest(), None)
                    else:
                        return crossing_edge.get_opp().get_face()
                else:
                    if inside_triangle_check(face.get_edge().get_dest().get_pos(),
                                             face.get_edge().get_next().get_dest().get_pos(),
                                             face.get_edge().get_next().get_next().get_dest().get_pos(), x_pos):
                        return face
                    else:
                        print('oops')
            else:
                # if x_pos == [19.2,18.4]:
                #     fig = pyplot.figure()
                #     ar = fig.add_subplot(1,1,1)
                #     tx = []
                #     ty = []
                #     rx = []
                #     ry = []
                #     lx = []
                #     ly = []
                #     tx.append(crossing_edge.get_next().get_dest().get_pos()[0])
                #     tx.append(crossing_edge.get_next().get_next().get_dest().get_pos()[0])
                #     tx.append(crossing_edge.get_next().get_dest().get_pos()[0])
                #     tx.append(crossing_edge.get_dest().get_pos()[0])
                #     tx.append(crossing_edge.get_next().get_next().get_dest().get_pos()[0])
                #     tx.append(crossing_edge.get_dest().get_pos()[0])
                #     rx.append(crossing_edge.get_next().get_dest().get_pos()[0])
                #     rx.append(x_pos[0])
                #     lx.append(face.get_edge().get_dest().get_pos()[0])
                #     lx.append(face.get_edge().get_next().get_dest().get_pos()[0])
                #     lx.append(face.get_edge().get_next().get_dest().get_pos()[0])
                #     lx.append(face.get_edge().get_next().get_next().get_dest().get_pos()[0])
                #     lx.append(face.get_edge().get_next().get_next().get_dest().get_pos()[0])
                #     lx.append(face.get_edge().get_dest().get_pos()[0])
                #     ty.append(crossing_edge.get_next().get_dest().get_pos()[1])
                #     ty.append(crossing_edge.get_next().get_next().get_dest().get_pos()[1])
                #     ty.append(crossing_edge.get_next().get_dest().get_pos()[1])
                #     ty.append(crossing_edge.get_dest().get_pos()[1])
                #     ty.append(crossing_edge.get_next().get_next().get_dest().get_pos()[1])
                #     ty.append(crossing_edge.get_dest().get_pos()[1])
                #     ry.append(crossing_edge.get_next().get_dest().get_pos()[1])
                #     ry.append(x_pos[1])
                #     ly.append(face.get_edge().get_dest().get_pos()[1])
                #     ly.append(face.get_edge().get_next().get_dest().get_pos()[1])
                #     ly.append(face.get_edge().get_next().get_dest().get_pos()[1])
                #     ly.append(face.get_edge().get_next().get_next().get_dest().get_pos()[1])
                #     ly.append(face.get_edge().get_next().get_next().get_dest().get_pos()[1])
                #     ly.append(face.get_edge().get_dest().get_pos()[1])
                #
                #     ar.plot(tx,ty, color='b')
                #     ar.plot(rx,ry, color='g')
                #     ar.plot(lx,ly, color='k')
                    # pyplot.show()
                    # print('wut')

                    # fig = pyplot.figure()
                    # ar = fig.add_subplot(1,1,1)
                    # tx = []
                    # ty = []
                    # rx = []
                    # ry = []
                    # tx.append(crossing_edge.get_next().get_dest().get_pos()[0])
                    # tx.append(crossing_edge.get_next().get_next().get_dest().get_pos()[0])
                    # tx.append(crossing_edge.get_next().get_dest().get_pos()[0])
                    # tx.append(crossing_edge.get_dest().get_pos()[0])
                    # tx.append(crossing_edge.get_next().get_next().get_dest().get_pos()[0])
                    # tx.append(crossing_edge.get_dest().get_pos()[0])
                    # rx.append(crossing_edge.get_next().get_dest().get_pos()[0])
                    # rx.append(x_pos[0])
                    # ty.append(crossing_edge.get_next().get_dest().get_pos()[1])
                    # ty.append(crossing_edge.get_next().get_next().get_dest().get_pos()[1])
                    # ty.append(crossing_edge.get_next().get_dest().get_pos()[1])
                    # ty.append(crossing_edge.get_dest().get_pos()[1])
                    # ty.append(crossing_edge.get_next().get_next().get_dest().get_pos()[1])
                    # ty.append(crossing_edge.get_dest().get_pos()[1])
                    # ry.append(crossing_edge.get_next().get_dest().get_pos()[1])
                    # ry.append(x_pos[1])
                    # ar.plot(tx,ty, color='b')
                    # pyplot.show()

                return find_face(x_pos, start_vertex, crossing_edge.get_next().get_opp().get_face())


def vertex_insertion(vertex, face, faces, vertices):
    #TODO: check if the new vertex falls on an edge that would need to be split

    #TODO: check if new vertex falls close to another vertex by a distance epsilon. if so, no vertex creation, return the closest vertex

    first_new_face = MeshFace(face.get_edge().get_next().get_next())
    face.get_edge().get_next().get_next().set_face(first_new_face)
    second_new_face = MeshFace(face.get_edge().get_next())
    face.get_edge().get_next().set_face(second_new_face)

    first_new_edge = MeshEdge(face.get_edge().get_opp().get_dest(),face, next_edge=face.get_edge())
    first_opp_edge = MeshEdge(vertex, first_new_face, opp_edge=first_new_edge)
    first_new_edge.set_opp(first_opp_edge)
    vertex.set_edge(first_new_edge)
    face.get_edge().get_next().get_next().set_next(first_opp_edge)

    second_new_edge = MeshEdge(face.get_edge().get_next().get_dest(), first_new_face,  next_edge=face.get_edge().get_next().get_next())
    first_opp_edge.set_next(second_new_edge)
    second_opp_edge = MeshEdge(vertex, second_new_face, opp_edge=second_new_edge)
    second_new_edge.set_opp(second_opp_edge)
    face.get_edge().get_next().set_next(second_opp_edge)

    third_new_edge = MeshEdge(face.get_edge().get_dest(), second_new_face, next_edge=face.get_edge().get_next())
    second_opp_edge.set_next(third_new_edge)
    third_opp_edge = MeshEdge(vertex, face, opp_edge=third_new_edge, next_edge= first_new_edge)
    face.get_edge().set_next(third_opp_edge)
    third_new_edge.set_opp(third_opp_edge)

    vertex_index = 0
    index_found = False
    while vertex_index <= len(vertices) and not index_found:
        if vertex.get_pos()[0] < vertices[vertex_index].get_pos()[0]:
            vertex_index += 1
        else:
            if vertex.get_pos()[1] < vertices[vertex_index].get_pos()[1]:
                vertex_index += 1
            else:
                index_found = True
    vertices.insert(vertex_index, vertex)
    faces.insert(faces.index(face), first_new_face)
    faces.insert(faces.index(face), second_new_face)

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

    ext_edge = vertex.get_edge().get_next()
    for fnum in [0,1,2]:
        if type(ext_edge.get_opp().get_face()) is MeshFace:
            if not circumcircle_check(ext_edge.get_dest(), ext_edge.get_next().get_dest(), ext_edge.get_opp().get_dest(), ext_edge.get_opp().get_next().get_dest()) \
                    and not ext_edge.is_obstacle():
                edge_flip(ext_edge)
                ext_edge = ext_edge.get_opp().get_next().get_next().get_opp().get_next()
            else:
                ext_edge = ext_edge.get_next().get_opp().get_next()

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

    return {'faces': faces, 'vertices': vertices}


def edge_slice(start_vertex, end_vertex, faces, vertices, face=None, epsilon=0.001, added_vertices=[], obstacle=True):
    if face is None:
        face = start_vertex.get_edge().get_face()
        return edge_slice(start_vertex, end_vertex, faces, vertices, face, epsilon, added_vertices, obstacle)

    else:
        if type(face) is GhostFace:
            border_edge = face.get_edge()
            if border_edge.get_opp().get_dest() == start_vertex:
                if not cone_check(start_vertex.get_pos(), border_edge.get_dest().get_pos(),
                                  border_edge.get_next().get_next().get_opp().get_next().get_next().get_opp().get_dest().get_pos() , end_vertex.get_pos()):
                    return edge_slice(start_vertex, end_vertex, faces, vertices, border_edge.get_next().get_next().get_opp().get_next().get_next().get_opp().get_face(), epsilon, added_vertices, obstacle)
                else:
                    return {'faces': faces, 'vertices': vertices, 'added_vertices':added_vertices}
            elif border_edge.get_dest() == start_vertex:
                if not cone_check(start_vertex.get_pos(), border_edge.get_next().get_opp().get_next().get_opp().get_dest().get_pos(),
                                  border_edge.get_opp().get_dest().get_pos() , end_vertex.get_pos()):
                    return edge_slice(start_vertex, end_vertex, faces, vertices, border_edge.get_opp().get_face(), epsilon, added_vertices, obstacle)
                else:
                    return {'faces': faces, 'vertices': vertices, 'added_vertices':added_vertices}
            else:
                print('oopsys')
        else:

            end_edge = face.get_edge()
            end_found = False
            i = 0
            while i < 2 and not end_found:
                if end_edge.get_dest().get_pos() == end_vertex.get_pos():
                    if end_edge.get_opp().get_dest() == start_vertex:
                        end_found = True
                    elif end_edge.get_next().get_dest() == start_vertex:
                        end_edge = end_edge.get_next().get_opp()
                        end_found = True
                else:
                    end_edge = end_edge.get_next()
                i += 1

            if end_found:
                print('end')
                start_vertex.set_edge(end_edge)
                print(end_edge.is_obstacle())
                if obstacle:
                    end_edge.set_obstacle(True)
                    # end_edge.get_opp().set_obstacle(True)
                print(end_edge.is_obstacle())
                print(end_edge.get_opp().get_dest().get_pos(), end_edge.get_dest().get_pos())

                return {'faces': faces, 'vertices': vertices, 'added_vertices':added_vertices}

            else:
                #TODO: check if the target line is on an edge
                crossing_edge = face.get_edge()
                if crossing_edge.get_dest() == start_vertex:
                    crossing_edge = crossing_edge.get_next().get_next()
                elif crossing_edge.get_next().get_next().get_dest() == start_vertex:
                    crossing_edge = crossing_edge.get_next()

                if cone_check(start_vertex.get_pos(), crossing_edge.get_next().get_next().get_dest().get_pos(), crossing_edge.get_dest().get_pos(), end_vertex.get_pos()):

                    intersection = intersect_check(start_vertex.get_pos()[0], start_vertex.get_pos()[1],
                                                   end_vertex.get_pos()[0], end_vertex.get_pos()[1],
                                                   crossing_edge.get_opp().get_dest().get_pos()[0],
                                                   crossing_edge.get_opp().get_dest().get_pos()[1],
                                                   crossing_edge.get_dest().get_pos()[0],
                                                   crossing_edge.get_dest().get_pos()[1])

                    if intersection != False:

                        if (abs(crossing_edge.get_dest().get_pos()[0] - intersection[0]) < epsilon) and (abs(crossing_edge.get_dest().get_pos()[1] - intersection[1]) < epsilon):
                            start_vertex.set_edge(crossing_edge.get_next().get_opp())
                            if obstacle:
                                crossing_edge.get_next().get_opp().set_obstacle(True)
                                # crossing_edge.get_next().set_obstacle(True)
                            return edge_slice(crossing_edge.get_dest(), end_vertex, faces, vertices, None, epsilon, added_vertices, obstacle)

                        elif (abs(crossing_edge.get_opp().get_dest().get_pos()[0] - intersection[0]) < epsilon) and (abs(crossing_edge.get_opp().get_dest().get_pos()[1] - intersection[1]) < epsilon):
                            start_vertex.set_edge(crossing_edge.get_next().get_next())
                            if obstacle:
                                crossing_edge.get_next().get_next().set_obstacle(True)
                                # crossing_edge.get_next().get_next().get_opp().set_obstacle(True)
                            return edge_slice(crossing_edge.get_opp().get_dest(), end_vertex, faces, vertices, None, epsilon, added_vertices, obstacle)

                        if type(crossing_edge.get_opp().get_face()) is MeshFace:
                            #TODO:  add a vetex at crossing point, split both triangles on the side of the crossing edge
                            #by adding edges and faces linking opping vertices, then iterate from the crossing point

                            crossing_vertex = MeshVertex(intersection[0], intersection[1])
                            added_vertices.append(crossing_vertex)
                            inward_new_edge = MeshEdge(crossing_edge.get_next().get_dest(), face= face, next_edge=crossing_edge.get_next().get_next())
                            face.set_edge(inward_new_edge)
                            inward_opp_edge = MeshEdge(crossing_vertex, opp_edge=inward_new_edge)
                            crossing_edge.get_next().get_dest().set_edge(inward_opp_edge)
                            inward_new_face = MeshFace(inward_opp_edge)
                            inward_opp_edge.set_face(inward_new_face)
                            crossing_edge.get_next().set_face(inward_new_face)
                            crossing_edge.get_next().set_next(inward_opp_edge)
                            inward_new_edge.set_opp(inward_opp_edge)
                            split_new_edge = MeshEdge(crossing_edge.get_dest(), inward_new_face, next_edge=crossing_edge.get_next())
                            inward_opp_edge.set_next(split_new_edge)
                            crossing_edge.set_next(inward_new_edge)
                            crossing_edge.set_dest(crossing_vertex)

                            outward_new_edge = MeshEdge(crossing_edge.get_opp().get_next().get_dest(), crossing_edge.get_opp().get_face(), next_edge=crossing_edge.get_opp().get_next().get_next())
                            crossing_edge.get_opp().get_face().set_edge(outward_new_edge)
                            outward_opp_edge = MeshEdge(crossing_vertex, opp_edge=outward_new_edge)
                            outward_new_edge.set_opp(outward_opp_edge)
                            outward_new_face = MeshFace(outward_opp_edge)
                            outward_opp_edge.set_face(outward_new_face)
                            crossing_edge.get_opp().get_next().set_next(outward_opp_edge)
                            crossing_edge.get_opp().get_next().set_face(outward_new_face)
                            split_opp_edge = MeshEdge(crossing_edge.get_opp().get_dest(), outward_new_face, crossing_edge, crossing_edge.get_opp().get_next())
                            outward_opp_edge.set_next(split_opp_edge)
                            crossing_edge.get_opp().set_opp(split_new_edge)
                            split_new_edge.set_opp(crossing_edge.get_opp())
                            crossing_edge.get_opp().set_dest(crossing_vertex)
                            crossing_edge.get_opp().set_next(outward_new_edge)
                            crossing_edge.set_opp(split_opp_edge)

                            #TODO might not stop (check stop condition) aout the last jump => no edge crossing, only a link to the pposite vertex

                            crossing_vertex.set_edge(outward_new_edge)

                            if obstacle:
                                # outward_new_edge.set_obstacle(True)
                                # outward_opp_edge.set_obstacle(True)
                                # inward_new_edge.set_obstacle(True)
                                inward_opp_edge.set_obstacle(True)

                                print(outward_new_edge.is_obstacle())

                            if type(split_new_edge.get_next().get_opp().get_face()) is MeshFace:
                                if clockwise_check(start_vertex, crossing_vertex, split_new_edge.get_next().get_opp().get_next().get_dest()):
                                    if not circumcircle_check(start_vertex, crossing_vertex, split_new_edge.get_dest(), split_new_edge.get_next().get_opp().get_next().get_dest()) \
                                            and not split_new_edge.get_next().is_obstacle() and not split_new_edge.get_next().get_opp().is_obstacle():

                                        #NOPE NOPE NOPE

                                        # split_new_edge.get_next().set_dest(outward_new_edge.get_dest())
                                        # split_new_edge.get_next().set_next(outward_new_edge.get_next())
                                        #
                                        # outward_new_edge.get_next().set_face(split_new_edge.get_next().get_opp().get_next().get_face())
                                        # split_new_edge.get_next().set_face(split_new_edge.get_next().get_opp().get_next().get_face())
                                        #
                                        #
                                        # split_new_edge.set_dest(split_new_edge.get_next().get_opp().get_next())
                                        # split_new_edge.set_next(split_new_edge.get_next().get_opp().get_next().get_next())
                                        #
                                        # split_new_edge.get_next().get_opp().get_next().set_next(split_new_edge.get_next())
                                        # split_new_edge.get_next().get_opp().get_next().get_face().set_edge(split_new_edge.get_next().get_opp().get_next())
                                        #
                                        #
                                        # split_new_edge.get_next().get_opp().get_next().get_next().set_next(inward_opp_edge)
                                        # split_new_edge.get_next().get_opp().get_next().get_next().set_face(inward_opp_edge.get_face())
                                        #
                                        # split_new_edge.get_next().get_opp().set_dest(split_new_edge.get_next().get_opp().get_next().get_dest())
                                        # split_new_edge.get_next().get_opp().set_next(split_new_edge.get_opp())
                                        # split_new_edge.get_next().get_opp().set_face(outward_new_edge.get_face())
                                        # outward_new_edge.set_next(split_new_edge.get_opp())
                                        #
                                        # split_new_edge.set_dest(split_new_edge.get_next().get_next().get_next().get_dest())
                                        # # split_new_edge.set_next(split_new_edge.get_next().get_next().get_next())

                                        edge_flip(split_new_edge.get_next())

                                if (not circumcircle_check(split_new_edge.get_dest(), split_new_edge.get_next().get_dest(), split_new_edge.get_next().get_next().get_dest(), split_new_edge.get_opp().get_next().get_dest())\
                                        or not circumcircle_check(split_new_edge.get_opp().get_dest(),split_new_edge.get_opp().get_next().get_dest(), split_new_edge.get_opp().get_next().get_next().get_dest(), split_new_edge.get_next().get_dest()))\
                                        and not split_new_edge.is_obstacle() and not split_new_edge.get_opp().is_obstacle():

                                    edge_flip(split_new_edge)

                            vertex_index = 0
                            index_found = False
                            while vertex_index <= len(vertices) and not index_found:
                                if crossing_vertex.get_pos()[0] < vertices[vertex_index].get_pos()[0]:
                                    vertex_index += 1
                                else:
                                    if crossing_vertex.get_pos()[1] < vertices[vertex_index].get_pos()[1]:
                                        vertex_index += 1
                                    else:
                                        index_found = True
                            vertices.insert(vertex_index, crossing_vertex)
                            faces.insert(0, inward_new_face)
                            faces.insert(0, outward_new_face)

                            # aff_faces([crossing_edge.get_face(), inward_new_face, outward_new_face, split_new_edge.get_opp().get_face()])


                            return edge_slice(crossing_vertex, end_vertex, faces, vertices, None, epsilon, added_vertices, obstacle)
                        else:
                            return {'faces': faces, 'vertices': vertices, 'added_vertices':added_vertices}
                    else:

                        print('oops')
                        return {'faces': faces, 'vertices': vertices, 'added_vertices':added_vertices}
                else:
                    return edge_slice(start_vertex, end_vertex, faces, vertices, crossing_edge.get_next().get_opp().get_face(), epsilon, added_vertices, obstacle)

def add_solid(hull_points, faces, vertices, faces_grid=None, square_len=None):

    init_hull_vertices = []
    for hull_point in hull_points:
        hull_vertex = MeshVertex(hull_point[0], hull_point[1])
        init_hull_vertices.append(hull_vertex)
        if faces_grid is not None and square_len is not None:
            insert_face = find_face(hull_point, faces_grid[int(hull_point[0]//square_len)][int(hull_point[1]//square_len)])
        else:
            insert_face = find_face(hull_point, vertices[0])
        hull = vertex_insertion(hull_vertex, insert_face, faces, vertices)
        faces = hull['faces']
        vertices = hull['vertices']

    hull_vertices = []
    hull_vertices.extend(init_hull_vertices)

    for hull_index in list(range(len(init_hull_vertices))):
        next_index = hull_index+1
        if hull_index == len(init_hull_vertices) -1:
            next_index = 0
        hull_edges = edge_slice(init_hull_vertices[hull_index], init_hull_vertices[next_index], faces, vertices)

        faces = hull_edges['faces']
        vertices = hull_edges['vertices']
        hull_vertices.extend(hull_edges['added_vertices'])

        #TODO: sort the added vertices by x then y pos to keep coherence
    return {'faces': faces, 'vertices': vertices, 'hull_vertices': hull_vertices}


def vertex_removal(vertex, vertices, faces):
    start_edge = vertex.get_edge()
    face = start_edge.get_next().get_face()
    face.set_edge(start_edge.get_next())

    current_edge = start_edge.get_next()
    stop_edge = start_edge.get_opp().get_next().get_next()

    while current_edge != stop_edge:
        current_edge.set_next(current_edge.get_next().get_opp().get_next())
        current_edge.get_next().get_next().get_next().get_opp().set_dest(None)
        current_edge.get_next().get_next().get_next().get_opp().set_opp(None)
        current_edge.get_next().get_next().get_next().get_opp().set_next(None)
        current_edge.get_next().get_next().get_next().get_opp().set_face(None)

        current_edge.get_next().get_next().get_next().set_dest(None)
        current_edge.get_next().get_next().get_next().set_opp(None)
        current_edge.get_next().get_next().get_next().set_next(None)
        current_edge.get_next().get_next().get_next().set_face(None)

        current_edge.get_next().get_face().set_edge(None)
        if current_edge.get_next().get_face() in faces:
            del faces[faces.index(current_edge.get_next().get_face())]

        current_edge.get_next().set_face(face)
        current_edge = current_edge.get_next()

    current_edge.get_next().get_opp().set_dest(None)
    current_edge.get_next().get_opp().set_oppp(None)
    current_edge.get_next().get_opp().set_next(None)
    current_edge.get_next().get_opp().set_face(None)

    current_edge.get_next().set_dest(None)
    current_edge.get_next().set_oppp(None)
    current_edge.get_next().set_next(None)
    current_edge.get_next().set_face(None)

    current_edge.set_next(start_edge)

    if vertex in vertices:
        del vertices[vertices.index(vertex)]

    return {'faces': faces, 'vertices': vertices}

def d_star(start_pos, end_pos, faces, vertices, faces_grid=None, square_len= None):

    start_vertex = MeshVertex(start_pos[0], start_pos[1])
    end_vertex = MeshVertex(end_pos[0], end_pos[1])
    if faces_grid is not None and square_len is not None:
        start_face = find_face(start_pos, faces_grid[int(start_pos[0]//square_len)][int(start_pos[1]//square_len)])
    else:
        start_face = find_face(start_pos, vertices[0])

    if faces_grid is not None and square_len is not None:
        end_face = find_face(end_pos, faces_grid[int(start_pos[0]//square_len)][int(start_pos[1]//square_len)])
    else:
        end_face = find_face(end_pos, vertices[0])

    #TODO: check if start or end vertex is too close to an existing vetex (if so, pick that one instead but without deleting it at the end)

    vertex_insertion(start_vertex, start_face, faces, vertices)
    vertex_insertion(end_vertex, end_face, faces, vertices)

    explored_vertices = [start_vertex]
    start_vertices_path = [start_vertex]
    cost_path = []
    start_edges_path = []
    # end_vertices_path = [end_vertex]

    while start_vertices_path[-1] != end_vertex:
        print('whole loop')
        for sv in start_vertices_path:
            print(sv.get_pos())
        first_edge = start_vertices_path[-1].get_edge()
        start_edge = first_edge

        while start_edge.get_dest() in explored_vertices:
            print('init loop')
            print(start_edge.get_opp().get_dest().get_pos(),start_edge.get_dest().get_pos())
            if type(start_edge.get_face()) is MeshFace:
                if start_edge.is_obstacle():
                    while not start_edge.get_opp().is_obstacle():
                        print('obstacle loop')
                        if type(start_edge.get_face()) is MeshFace:
                            start_edge = start_edge.get_next().get_next().get_opp()
                        else:
                            start_edge = start_edge.get_next().get_next().get_opp().get_next().get_next().get_opp()
                    if start_edge == first_edge:
                        del start_vertices_path[-1]
                        if len(start_edges_path) >0:
                            del start_edges_path[-1]
                        print('start over')
                        print(start_vertices_path[-1].get_edge().get_dest().get_pos())
                        start_edge = start_vertices_path[-1].get_edge()
                        first_edge = start_edge
                else:
                    start_edge = start_edge.get_next().get_next().get_opp()
                    if start_edge == first_edge:
                        del start_vertices_path[-1]
                        if len(start_edges_path) >0:
                            del start_edges_path[-1]
                        print('start over')
                        print(start_vertices_path[-1].get_edge().get_dest().get_pos())
                        start_edge = start_vertices_path[-1].get_edge()
                        first_edge = start_edge
            else:
                start_edge = start_edge.get_next().get_next().get_opp().get_next().get_next().get_opp()
                if start_edge == first_edge:
                    del start_vertices_path[-1]
                    if len(start_edges_path) >0:
                        del start_edges_path[-1]
                    print('start over')
                    print(start_vertices_path[-1].get_edge().get_dest().get_pos())
                    start_edge = start_vertices_path[-1].get_edge()
                    first_edge = start_edge

        best_edge = start_edge
        if start_edge.get_dest() == end_vertex:
            best_angle = 0
            cost_path.append(pow(best_edge.get_dest().get_pos()[0] - start_vertices_path[-1].get_pos()[0],2)+ pow(pow(best_edge.get_dest().get_pos()[1] - start_vertices_path[-1].get_pos()[1],2),2))
            start_vertices_path.append(best_edge.get_dest())
            explored_vertices.append(best_edge.get_dest())
            start_edges_path.append(best_edge)
        else:
            print('yo')
            print(end_vertex.get_pos(), start_vertices_path[-1].get_pos())
            best_angle = diamond_diff(diamond_angle(end_vertex.get_pos()[0]-start_vertices_path[-1].get_pos()[0], end_vertex.get_pos()[1]-start_vertices_path[-1].get_pos()[1]),
                            diamond_angle(best_edge.get_dest().get_pos()[0]-start_vertices_path[-1].get_pos()[0], best_edge.get_dest().get_pos()[1]-start_vertices_path[-1].get_pos()[1]))

            if type(start_edge.get_face()) is MeshFace:
                if start_edge.is_obstacle():
                    current_edge = start_edge.get_next().get_next().get_opp()
                    while not current_edge.get_opp().is_obstacle():
                        if type(current_edge.get_face()) is MeshFace:
                            current_edge = current_edge.get_next().get_next().get_opp()
                        else:
                            current_edge = current_edge.get_next().get_next().get_opp().get_next().get_next().get_opp()
                else:
                    current_edge = start_edge.get_next().get_next().get_opp()

            else:
                current_edge = start_edge.get_next().get_next().get_opp().get_next().get_next().get_opp()

            if current_edge.get_dest() == end_vertex:
                best_angle = 0
                cost_path.append(pow(current_edge.get_dest().get_pos()[0] - start_vertices_path[-1].get_pos()[0],2)+ pow(pow(current_edge.get_dest().get_pos()[1] - start_vertices_path[-1].get_pos()[1],2),2))
                start_vertices_path.append(current_edge.get_dest())
                explored_vertices.append(current_edge.get_dest())
                start_edges_path.append(current_edge)
            else:
                print('current')
                print(start_vertices_path[-1].get_pos())
                print(end_vertex.get_pos())
                print(current_edge.get_dest().get_pos())
                current_angle = diamond_diff(diamond_angle(end_vertex.get_pos()[0]-start_vertices_path[-1].get_pos()[0], end_vertex.get_pos()[1]-start_vertices_path[-1].get_pos()[1]),
                                    diamond_angle(current_edge.get_dest().get_pos()[0]-start_vertices_path[-1].get_pos()[0], current_edge.get_dest().get_pos()[1]-start_vertices_path[-1].get_pos()[1]))

                #TODO: check for obstacles to stop turning around the vertex and ghost faces to ignore and jump though
                while current_edge != start_edge and current_edge.get_dest() != end_vertex:
                    print('search loop')
                    if current_edge.get_dest() not in explored_vertices:
                        if current_angle < best_angle:

                            best_edge = current_edge
                            best_angle = current_angle

                    if type(current_edge.get_face()) is MeshFace:
                        if current_edge.is_obstacle():
                            while not current_edge.get_opp().is_obstacle():
                                if type(current_edge.get_face()) is MeshFace:
                                    current_edge = current_edge.get_next().get_next().get_opp()
                                else:
                                    current_edge = current_edge.get_next().get_next().get_opp().get_next().get_next().get_opp()
                        else:
                            current_edge = current_edge.get_next().get_next().get_opp()
                    else:
                        current_edge = current_edge.get_next().get_next().get_opp().get_next().get_next().get_opp()

                    print('current turn')
                    print(start_vertices_path[-1].get_pos())
                    print(end_vertex.get_pos())
                    print(current_edge.get_dest().get_pos())
                    current_angle = diamond_diff(diamond_angle(end_vertex.get_pos()[0]-start_vertices_path[-1].get_pos()[0], end_vertex.get_pos()[1]-start_vertices_path[-1].get_pos()[1]),
                                        diamond_angle(current_edge.get_dest().get_pos()[0]-start_vertices_path[-1].get_pos()[0], current_edge.get_dest().get_pos()[1]--start_vertices_path[-1].get_pos()[1]))

                if current_edge.get_dest() == end_vertex:
                    best_angle = 0
                    cost_path.append(pow(current_edge.get_dest().get_pos()[0] - start_vertices_path[-1].get_pos()[0],2)+ pow(pow(current_edge.get_dest().get_pos()[1] - start_vertices_path[-1].get_pos()[1],2),2))
                    start_vertices_path.append(current_edge.get_dest())
                    explored_vertices.append(current_edge.get_dest())
                    start_edges_path.append(current_edge)

                print('chosen')
                print(best_edge.get_dest().get_pos())
                cost_path.append(pow(best_edge.get_dest().get_pos()[0] - start_vertices_path[-1].get_pos()[0],2)+ pow(pow(best_edge.get_dest().get_pos()[1] - start_vertices_path[-1].get_pos()[1],2),2))
                start_vertices_path.append(best_edge.get_dest())
                explored_vertices.append(best_edge.get_dest())
                start_edges_path.append(best_edge)

                # segz = []
                # wegz = []
                # for f in faces:
                #     if type(f) is MeshFace:
                #         e = f.get_edge()
                #         for i in [0,1,2,3]:
                #             if e.is_obstacle():
                #                 wegz.append([e.get_opp().get_dest().get_pos(), e.get_dest().get_pos()])
                #             else:
                #                 segz.append([e.get_opp().get_dest().get_pos(), e.get_dest().get_pos()])
                #             e = e.get_next()
                # fig = pyplot.figure()
                # ar = fig.add_subplot(1,1,1)
                # sar = array(segz)
                # for h in sar:
                #     tx = h[:,0]
                #     ty = h[:,1]
                #     ar.plot(tx,ty, color='r')
                # war = array(wegz)
                # for w in war:
                #     wx = w[:,0]
                #     wy = w[:,1]
                #     ar.plot(wx,wy, color='k')
                # pegz = []
                # for vindex in list(range(len(start_vertices_path)-1)):
                #     pegz.append([start_vertices_path[vindex].get_pos(), start_vertices_path[vindex+1].get_pos()])
                # par = array(pegz)
                # for g in par:
                #     rx = g[:,0]
                #     ry = g[:,1]
                #     ar.plot(rx,ry, color='b')
                # pyplot.show()

    segz = []
    wegz = []
    for f in faces:
        if type(f) is MeshFace:
            e = f.get_edge()
            for i in [0,1,2,3]:
                if e.is_obstacle():
                    wegz.append([e.get_opp().get_dest().get_pos(), e.get_dest().get_pos()])
                else:
                    segz.append([e.get_opp().get_dest().get_pos(), e.get_dest().get_pos()])
                e = e.get_next()
    fig = pyplot.figure()
    ar = fig.add_subplot(1,1,1)
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
    pegz = []
    for vindex in list(range(len(start_vertices_path)-1)):
        pegz.append([start_vertices_path[vindex].get_pos(), start_vertices_path[vindex+1].get_pos()])
    par = array(pegz)
    for g in par:
        rx = g[:,0]
        ry = g[:,1]
        ar.plot(rx,ry, color='b')
    # pyplot.show()


    return {'vertices': start_vertices_path, 'edges': start_edges_path, 'costs': cost_path}

def external_crossing(box_pos, square_len, a_pos, b_pos):
    if a_pos[0] != b_pos[0]:
        A = (a_pos[1]-b_pos[1])/(a_pos[0]-b_pos[0])
        b = b_pos[1]-A*b_pos[0]
        borders = []
        if b_pos[0] > a_pos[0]:
            if (A * box_pos[0] + b) >= box_pos[1] and (A * box_pos[0] + b) < box_pos[1]+square_len:
                borders.append('west')
                if (A * (box_pos[0]+square_len) + b) >= box_pos[1] \
                        and (A * (box_pos[0]+square_len) + b) < box_pos[1]+square_len:
                    borders.append('est')
                else:
                    if b_pos[1] >= a_pos[1]:
                        if ((box_pos[1]+square_len - b) / A) >= box_pos[0] and ((box_pos[1]+square_len - b) / A) < box_pos[0]+square_len:
                            borders.append('north')
                    else:
                        if ((box_pos[1] - b) / A) >= box_pos[0] and ((box_pos[1] - b) / A) < box_pos[0]+square_len:
                            borders.append('south')
            else:
                if b_pos[1] >= a_pos[1]:
                    if ((box_pos[1] - b) / A) >= box_pos[0] and ((box_pos[1] - b) / A) < box_pos[0]+square_len:
                        borders.append('south')
                        if ((box_pos[1]+square_len - b) / A) >= box_pos[0] and ((box_pos[1]+square_len - b) / A) < box_pos[0]+square_len:
                            borders.append('north')
                        elif (A * (box_pos[0]+square_len) + b) >= box_pos[1] \
                                and (A * (box_pos[0]+square_len) + b) < box_pos[1]+square_len:
                            borders.append('est')
                    else:
                        return None
                else:
                    if ((box_pos[1]+square_len - b) / A) >= box_pos[0] and ((box_pos[1]+square_len - b) / A) < box_pos[0]+square_len:
                        borders.append('north')
                        if ((box_pos[1] - b) / A) >= box_pos[0] and ((box_pos[1] - b) / A) < box_pos[0]+square_len:
                            borders.append('south')
                        elif (A * (box_pos[0]+square_len) + b) >= box_pos[1] \
                                and (A * (box_pos[0]+square_len) + b) < box_pos[1]+square_len:
                            borders.append('est')
                    else:
                        return None
        else:
            if (A * (box_pos[0]+square_len) + b) >= box_pos[1] \
                    and (A * (box_pos[0]+square_len) + b) < box_pos[1]+square_len:
                borders.append('est')
                if (A * box_pos[0] + b) >= box_pos[1] and (A * box_pos[0] + b) < box_pos[1]+square_len:
                    borders.append('west')
                else:
                    if b_pos[1] >= a_pos[1]:
                        if ((box_pos[1]+square_len - b) / A) >= box_pos[0] and ((box_pos[1]+square_len - b) / A) < box_pos[0]+square_len:
                            borders.append('north')
                    else:
                        if ((box_pos[1] - b) / A) >= box_pos[0] and ((box_pos[1] - b) / A) < box_pos[0]+square_len:
                            borders.append('south')
            else:
                if b_pos[1] >= a_pos[1]:
                    if ((box_pos[1] - b) / A) >= box_pos[0] and ((box_pos[1] - b) / A) < box_pos[0]+square_len:
                        borders.append('south')
                        if (A * box_pos[0] + b) >= box_pos[1] and (A * box_pos[0] + b) < box_pos[1]+square_len:
                            borders.append('west')
                        elif ((box_pos[1]+square_len - b) / A) >= box_pos[0] and ((box_pos[1]+square_len - b) / A) < box_pos[0]+square_len:
                            borders.append('north')
                    else:
                        return None
                else:
                    if ((box_pos[1]+square_len - b) / A) >= box_pos[0] and ((box_pos[1]+square_len - b) / A) < box_pos[0]+square_len:
                        borders.append('north')
                        if (A * box_pos[0] + b) >= box_pos[1] and (A * box_pos[0] + b) < box_pos[1]+square_len:
                            borders.append('west')
                        elif ((box_pos[1] - b) / A) >= box_pos[0] and ((box_pos[1] - b) / A) < box_pos[0]+square_len:
                            borders.append('south')
                    else:
                        return None
    else:
        if b_pos[1] >= a_pos[1]:
            if a_pos[0] >= box_pos[0] and a_pos[0] < box_pos[0]+square_len:
                return ['south','north']
            else:
                return None
        else:
            if a_pos[0] >= box_pos[0] and a_pos[0] < box_pos[0]+square_len:
                return ['north','south']
            else:
                return None

def edge_crossing(box_pos, square_len, a_pos, b_pos):
    if a_pos[0] != b_pos[0]:
        A = (a_pos[1]-b_pos[1])/(a_pos[0]-b_pos[0])
        b = b_pos[1]-A*b_pos[0]
        if b_pos[0] > a_pos[0]:
            if (A * (box_pos[0]+square_len) + b) >= box_pos[1] \
                    and (A * (box_pos[0]+square_len) + b) < box_pos[1]+square_len:
                return 'est'
            else:
                if b_pos[1] >= a_pos[1]:
                    if ((box_pos[1]+square_len - b) / A) >= box_pos[0] and ((box_pos[1]+square_len - b) / A) < box_pos[0]+square_len:
                        return 'north'
                else:
                    if ((box_pos[1] - b) / A) >= box_pos[0] and ((box_pos[1] - b) / A) < box_pos[0]+square_len:
                        return 'south'
        else:
            if (A * box_pos[0] + b) >= box_pos[1] and (A * box_pos[0] + b) < box_pos[1]+square_len:
                return 'west'
            else:
                if b_pos[1] >= a_pos[1]:
                    if ((box_pos[1]+square_len - b) / A) >= box_pos[0] and ((box_pos[1]+square_len - b) / A) < box_pos[0]+square_len:
                        return 'north'
                else:
                    if ((box_pos[1] - b) / A) >= box_pos[0] and ((box_pos[1] - b) / A) < box_pos[0]+square_len:
                        return 'south'
    else:
        if b_pos[1] >= a_pos[1]:
            if b_pos[1] >= box_pos[1]:
                return 'north'
        else:
            if b_pos[1] < box_pos[1]+square_len:
                return 'south'

def round_swipe(start_edge, centers_list, box_pos, square_len, split):

    gone_round = False
    current_edge = start_edge
    crossed_faces = []
    centers = centers_list
    borders_crossed = {'north':[], 'est':[], 'south':[],'west':[]}


    while not (current_edge == start_edge and gone_round):
        if type(current_edge.get_face()) is GhostFace:
            current_edge = current_edge.get_next().get_next().get_opp()
        else:
            centers.append(current_edge.get_opp().get_dest())
            gone_round = True
            if current_edge.get_next().get_next().get_opp().get_dest().get_pos()[0] >= box_pos[0]+square_len \
                    or current_edge.get_next().get_next().get_opp().get_dest().get_pos()[0] < box_pos[0] \
                    or current_edge.get_next().get_next().get_opp().get_dest().get_pos()[1] >= box_pos[1]+square_len \
                    or current_edge.get_next().get_next().get_opp().get_dest().get_pos()[1] < box_pos[1]:

                next_border_crossed = edge_crossing(box_pos, square_len, current_edge.get_next().get_next().get_dest().get_pos(),current_edge.get_next().get_next().get_opp().get_dest().get_pos())

                if next_border_crossed is not None:
                    if current_edge.get_next().get_next().get_face() not in borders_crossed[next_border_crossed]:
                        borders_crossed[next_border_crossed].append(current_edge.get_next().get_next().get_face())

                if current_edge.get_dest().get_pos()[0] >= box_pos[0]+square_len \
                        or current_edge.get_dest().get_pos()[0] < box_pos[0] \
                        or current_edge.get_dest().get_pos()[1] >= box_pos[1]+square_len \
                        or current_edge.get_dest().get_pos()[1] < box_pos[1]:

                    crossed_faces.append(current_edge.get_face())
                    border_crossed =edge_crossing(box_pos, square_len, current_edge.get_opp().get_dest().get_pos(), current_edge.get_dest().get_pos())

                    if border_crossed is not None:
                        if border_crossed != next_border_crossed:
                            if current_edge.get_face() not in borders_crossed[border_crossed]:
                                borders_crossed[border_crossed].append(current_edge.get_face())

                            full_crossing = external_crossing(box_pos, square_len, current_edge.get_dest().get_pos(), current_edge.get_next().get_dest().get_pos())
                            if split is not None:
                                if current_edge.get_next() != split() and current_edge.get_next().get_opp() != split():
                                    if full_crossing is not None:
                                        while full_crossing is not None:
                                            for border in full_crossing:
                                                if current_edge.get_face() not in borders_crossed[border]:
                                                    borders_crossed[border].append(current_edge.get_face())
                                            if type(current_edge.get_next().get_opp().get_face()) is MeshFace:
                                                if current_edge.get_next().get_opp().get_next().get_dest().get_pos()[0] >= box_pos[0] \
                                                        and current_edge.get_next().get_opp().get_next().get_dest().get_pos()[0] < box_pos[0]+square_len \
                                                        and current_edge.get_next().get_opp().get_next().get_dest().get_pos()[1] >= box_pos[1] \
                                                        and current_edge.get_next().get_opp().get_next().get_dest().get_pos()[1] < box_pos[1]+square_len:
                                                    full_crossing = None
                                                    centers.append(current_edge.get_next().get_opp().get_next().get_dest())
                                                    swipe = round_swipe(current_edge.get_next().get_opp().get_next(), centers, box_pos, square_len, current_edge.get_next())
                                                    other_faces = swipe['faces']
                                                    for box_border in list(swipe['borders'].keys()):
                                                        if swipe['borders'][box_border] != []:
                                                            borders_crossed[box_border].extend(swipe['borders'][box_border])
                                                    crossed_faces.append(other_faces)
                                                else:
                                                    full_L_crossing = external_crossing(box_pos, square_len, current_edge.get_next().get_opp().get_dest().get_pos(), current_edge.get_next().get_opp().get_next().get_dest().get_pos())
                                                    full_R_crossing = external_crossing(box_pos, square_len, current_edge.get_next().get_opp().get_next().get_dest().get_pos(), current_edge.get_next().get_opp().get_next().get_next().get_dest().get_pos())
                                                    if full_L_crossing is not None:
                                                        current_edge = current_edge.get_next().get_opp()
                                                        full_crossing = full_L_crossing
                                                    elif full_R_crossing is not None:
                                                        current_edge = current_edge.get_next().get_opp().get_next()
                                                        full_crossing = full_R_crossing
                                                    else:
                                                        if full_crossing == ['north', 'south']:
                                                            borders_crossed['west'].append(current_edge.get_face())
                                                        elif full_crossing == ['est', 'west']:
                                                            borders_crossed['north'].append(current_edge.get_face())
                                                        elif full_crossing == ['south', 'north']:
                                                            borders_crossed['est'].append(current_edge.get_face())
                                                        elif full_crossing == ['west', 'est']:
                                                            borders_crossed['south'].append(current_edge.get_face())
                                                        full_crossing = None
                                            else:
                                                full_crossing = None
                                    else:
                                        if next_border_crossed == 'north' and border_crossed == 'south':
                                            borders_crossed['west'].append(current_edge.get_face())
                                        elif next_border_crossed == 'est' and border_crossed == 'west':
                                            borders_crossed['north'].append(current_edge.get_face())
                                        elif next_border_crossed == 'south' and border_crossed == 'north':
                                            borders_crossed['est'].append(current_edge.get_face())
                                        elif next_border_crossed == 'west' and border_crossed == 'est':
                                            borders_crossed['south'].append(current_edge.get_face())
                                    #!!!  check for full edge crossing !!!
                                    #in that case, if the next vertex is in the square, call a recurtion on it
                                    # otherwise add the opposite face of the external crossing edge and iterate
                                    #until the next edge doesn't cross the box boundaries
                else:
                    if current_edge.get_dest() not in centers:
                        crossed_faces.append(current_edge.get_face())
                        centers.append(current_edge.get_dest())
                        swipe = round_swipe(current_edge.get_opp(), centers, box_pos, square_len)
                        other_faces = swipe['faces']
                        for box_border in list(swipe['borders'].keys()):
                            if swipe['borders'][box_border] != []:
                                borders_crossed[box_border].extend(swipe['borders'][box_border])
                        crossed_faces.append(other_faces)
            else:
                if current_edge.get_next().get_next().get_opp().get_dest() not in centers:
                    if current_edge.get_dest().get_pos()[0] >= box_pos[0]+square_len \
                            or current_edge.get_dest().get_pos()[0] < box_pos[0] \
                            or current_edge.get_dest().get_pos()[1] >= box_pos[1]+square_len \
                            or current_edge.get_dest().get_pos()[1] < box_pos[1]:
                        crossed_faces.append(current_edge.get_face())
                    else:
                        if current_edge.get_dest() not in centers:
                            crossed_faces.append(current_edge.get_face())
                            centers.append(current_edge.get_dest())
                            swipe = round_swipe(current_edge.get_opp(), centers, box_pos, square_len)
                            other_faces = swipe['faces']
                            for box_border in list(swipe['borders'].keys()):
                                if swipe['borders'][box_border] != []:
                                    borders_crossed[box_border].extend(swipe['borders'][box_border])
                            crossed_faces.append(other_faces)

            current_edge = current_edge.get_next().get_next().get_opp()

    return {'faces':crossed_faces, 'borders': borders_crossed}

def box_collision(box_pos, square_len, vertex):

    start_edge = vertex.get_edge()
    if start_edge.get_dest() == vertex:
        start_edge = start_edge.get_opp()

    return round_swipe(start_edge, [vertex], box_pos, square_len)

def propagate_edge_crossing(grid_faces, square_len):
    for i in list(range(len(grid_faces))):
        for j in list(range(len(grid_faces[0]))):
            if grid_faces[i][j] is not None:
                borders = ['north', 'est', 'south', 'west']
                if i == 0:
                    borders.remove('west')
                elif i == len(grid_faces):
                    borders.remove('est')
                if j == 0:
                    borders.remove('south')
                elif j == len(grid_faces[0]):
                    borders.remove('north')
                for border in borders:
                    if grid_faces[i][j][border] != []:
                        if border == 'est':
                            if grid_faces[i+1][j] is None:
                                grid_faces[i+1][j] = {'faces':[], 'borders': {'north':[], 'est':[], 'south':[],'west':[]}}
                                grid_faces[i+1][j]['faces'].extend(grid_faces[i][j]['borders'][border])
                                grid_faces[i+1][j]['borders']['west'].extend(grid_faces[i][j]['borders'][border])
                                single_cell = True
                                if len(grid_faces[i+1][j]['border']['west']) != 1:
                                    single_cell = False
                                expand_edge = []
                                for face in grid_faces[i+1][j]['borders']['west']:
                                    if type(face) is MeshFace:
                                        edge = face.get_edge()

                                        for k in [0,1,2]:
                                            full_crossing = external_crossing([j*square_len, (i+1)*square_len], square_len, edge.get_dest().get_pos(), edge.get_opp().get_dest().get_pos())
                                            if full_crossing is not None:
                                                single_cell = False
                                                if 'west' in full_crossing:
                                                    full_crossing.remove('west')
                                                    if face not in grid_faces[i+1][j]['borders'][full_crossing[0]]:
                                                        grid_faces[i+1][j]['borders'][full_crossing[0]].append(face)
                                                    #Transversal case
                                                    if full_crossing[0] == 'est':
                                                        if type(edge.get_opp().get_face()) is GhostFace:
                                                            if edge.get_dest().get_pos()[0] <  (i+1)*square_len:
                                                                grid_faces[i+1][j]['borders']['south'].append(edge.get_opp().get_face())
                                                            elif edge.get_dest().get_pos()[0] >=  (i+1)*square_len+square_len:
                                                                grid_faces[i+1][j]['borders']['north'].append(edge.get_opp().get_face())
                                                        other_crossing = False
                                                        next_crossing = external_crossing([(i+1)*square_len, j*square_len], square_len, edge.get_next().get_dest().get_pos(), edge.get_next().get_opp().get_dest().get_pos())
                                                        if next_crossing is not None:
                                                            other_crossing = True
                                                        next_crossing = external_crossing([(i+1)*square_len, j*square_len], square_len, edge.get_next().get_next().get_dest().get_pos(), edge.get_next().get_next().get_opp().get_dest().get_pos())
                                                        if next_crossing is not None:
                                                            other_crossing = True
                                                        if not other_crossing:
                                                            if edge.get_dest().get_pos()[0] <  (i+1)*square_len:
                                                                grid_faces[i+1][j]['borders']['north'].append(face)
                                                            elif edge.get_dest().get_pos()[0] >=  (i+1)*square_len+square_len:
                                                                grid_faces[i+1][j]['borders']['south'].append(face)
                                                else:
                                                    if face not in grid_faces[i+1][j]['borders'][full_crossing[0]]:
                                                        grid_faces[i+1][j]['borders'][full_crossing[0]].append(face)
                                                    if face not in grid_faces[i+1][j]['borders'][full_crossing[1]]:
                                                        grid_faces[i+1][j]['borders'][full_crossing[1]].append(face)
                                                    if edge.get_opp not in grid_faces[i+1][j]['faces']:
                                                        expand_edge.append(edge.get_opp())
                                            edge = edge.get_next()
                                if expand_edge == [] and single_cell:
                                    for other_border in ['south', 'est', 'north']:
                                        grid_faces[i+1][j]['borders'][other_border].extend(grid_faces[i+1][j]['borders']['west'])
                                else:
                                    while expand_edge != []:
                                        edge = expand_edge[0]
                                        expand_edge = expand_edge[1:]
                                        crossing = external_crossing([(i+1)*square_len, j*square_len], square_len, edge.get_dest().get_pos(), edge.get_opp().get_dest().get_pos())
                                        if edge.get_face() not in grid_faces[i+1][j]['faces']:
                                            grid_faces[i+1][j]['faces'].append(edge.get_face())
                                        if edge.get_face() not in grid_faces[i+1][j]['borders'][crossing[0]]:
                                            grid_faces[i+1][j]['borders'][crossing[0]].append(edge.get_face())
                                        if edge.get_face() not in grid_faces[i+1][j]['borders'][crossing[1]]:
                                            grid_faces[i+1][j]['borders'][crossing[1]].append(edge.get_face())

                                        other_crossing = False
                                        if type(edge.get_face()) is MeshFace:
                                            next_crossing = external_crossing([(i+1)*square_len, j*square_len], square_len, edge.get_next().get_dest().get_pos(), edge.get_next().get_opp().get_dest().get_pos())
                                            if next_crossing is not None:
                                                expand_edge.append(edge.get_next().get_opp())
                                                other_crossing = True
                                            next_crossing = external_crossing([(i+1)*square_len, j*square_len], square_len, edge.get_next().get_next().get_dest().get_pos(), edge.get_next().get_next().get_opp().get_dest().get_pos())
                                            if next_crossing is not None:
                                                expand_edge.append(edge.get_next().get_next().get_opp())
                                                other_crossing = True
                                        if not other_crossing and (crossing == ['north', 'south'] or crossing == ['south', 'north']):
                                            grid_faces[i+1][j]['borders']['est'].append(edge.get_face())
                        elif border == 'north':
                            if grid_faces[i][j+1] is None:
                                grid_faces[i][j+1] = {'faces':[], 'borders': {'north':[], 'est':[], 'south':[],'west':[]}}
                                grid_faces[i][j+1]['faces'].extend(grid_faces[i][j]['borders'][border])
                                grid_faces[i][j+1]['borders']['south'].extend(grid_faces[i][j]['borders'][border])
                                single_cell = True
                                if len(grid_faces[i][j+1]['border']['south']) != 1:
                                    single_cell = False
                                expand_edge = []
                                for face in grid_faces[i][j+1]['borders']['south']:
                                    if type(face) is MeshFace:
                                        edge = face.get_edge()

                                        for k in [0,1,2]:
                                            full_crossing = external_crossing([i*square_len, (j+1)*square_len], square_len, edge.get_dest().get_pos(), edge.get_opp().get_dest().get_pos())
                                            if full_crossing is not None:
                                                single_cell = False
                                                if 'south' in full_crossing:
                                                    full_crossing.remove('south')
                                                    if face not in grid_faces[i][j+1]['borders'][full_crossing[0]]:
                                                        grid_faces[i][j+1]['borders'][full_crossing[0]].append(face)
                                                    #Transversal case
                                                    if full_crossing[0] == 'north':
                                                        if type(edge.get_opp().get_face()) is GhostFace:
                                                            if edge.get_dest().get_pos()[1] <  (j+1)*square_len:
                                                                grid_faces[i][j+1]['borders']['est'].append(edge.get_opp().get_face())
                                                            elif edge.get_dest().get_pos()[1] >=  (j+1)*square_len+square_len:
                                                                grid_faces[i][j+1]['borders']['west'].append(edge.get_opp().get_face())
                                                        other_crossing = False
                                                        next_crossing = external_crossing([i*square_len, (j+1)*square_len], square_len, edge.get_next().get_dest().get_pos(), edge.get_next().get_opp().get_dest().get_pos())
                                                        if next_crossing is not None:
                                                            other_crossing = True
                                                        next_crossing = external_crossing([i*square_len, (j+1)*square_len], square_len, edge.get_next().get_next().get_dest().get_pos(), edge.get_next().get_next().get_opp().get_dest().get_pos())
                                                        if next_crossing is not None:
                                                            other_crossing = True
                                                        if not other_crossing:
                                                            if edge.get_dest().get_pos()[1] < (j+1)*square_len:
                                                                grid_faces[i][j+1]['borders']['west'].append(face)
                                                            elif edge.get_dest().get_pos()[1] >=  (j+1)*square_len+square_len:
                                                                grid_faces[i][j+1]['borders']['est'].append(face)
                                                else:
                                                    if face not in grid_faces[i][j+1]['borders'][full_crossing[0]]:
                                                        grid_faces[i][j+1]['borders'][full_crossing[0]].append(face)
                                                    if face not in grid_faces[i][j+1]['borders'][full_crossing[1]]:
                                                        grid_faces[i][j+1]['borders'][full_crossing[1]].append(face)
                                                    if edge.get_opp not in grid_faces[i][j+1]['faces']:
                                                        expand_edge.append(edge.get_opp())
                                            edge = edge.get_next()
                                if expand_edge == [] and single_cell:
                                    for other_border in ['west', 'north', 'est']:
                                        grid_faces[i][j+1]['borders'][other_border].extend(grid_faces[i][j+1]['borders']['south'])
                                else:
                                    while expand_edge != []:
                                        edge = expand_edge[0]
                                        expand_edge = expand_edge[1:]
                                        crossing = external_crossing([i*square_len, (j+1)*square_len], square_len, edge.get_dest().get_pos(), edge.get_opp().get_dest().get_pos())
                                        if edge.get_face() not in grid_faces[i][j+1]['faces']:
                                            grid_faces[i][j+1]['faces'].append(edge.get_face())
                                        if edge.get_face() not in grid_faces[i][j+1]['borders'][crossing[0]]:
                                            grid_faces[i][j+1]['borders'][crossing[0]].append(edge.get_face())
                                        if edge.get_face() not in grid_faces[i][j+1]['borders'][crossing[1]]:
                                            grid_faces[i][j+1]['borders'][crossing[1]].append(edge.get_face())

                                        other_crossing = False
                                        if type(edge.get_face()) is MeshFace:
                                            next_crossing = external_crossing([i*square_len, (j+1)*square_len], square_len, edge.get_next().get_dest().get_pos(), edge.get_next().get_opp().get_dest().get_pos())
                                            if next_crossing is not None:
                                                expand_edge.append(edge.get_next().get_opp())
                                                other_crossing = True
                                            next_crossing = external_crossing([i*square_len, (j+1)*square_len], square_len, edge.get_next().get_next().get_dest().get_pos(), edge.get_next().get_next().get_opp().get_dest().get_pos())
                                            if next_crossing is not None:
                                                expand_edge.append(edge.get_next().get_next().get_opp())
                                                other_crossing = True
                                        if not other_crossing and (crossing == ['est', 'west'] or crossing == ['west', 'est']):
                                            grid_faces[i][j+1]['borders']['north'].append(edge.get_face())
                        elif border == 'west':
                            if grid_faces[i-1][j] is None:
                                grid_faces[i-1][j] = {'faces':[], 'borders': {'north':[], 'est':[], 'south':[],'west':[]}}
                                grid_faces[i-1][j]['faces'].extend(grid_faces[i][j]['borders'][border])
                                grid_faces[i-1][j]['borders']['est'].extend(grid_faces[i][j]['borders'][border])
                                single_cell = True
                                if len(grid_faces[i-1][j]['border']['est']) != 1:
                                    single_cell = False
                                expand_edge = []
                                for face in grid_faces[i-1][j]['borders']['est']:
                                    if type(face) is MeshFace:
                                        edge = face.get_edge()

                                        for k in [0,1,2]:
                                            full_crossing = external_crossing([(i-1)*square_len, j*square_len], square_len, edge.get_dest().get_pos(), edge.get_opp().get_dest().get_pos())
                                            if full_crossing is not None:
                                                single_cell = False
                                                if 'est' in full_crossing:
                                                    full_crossing.remove('est')
                                                    if face not in grid_faces[i-1][j]['borders'][full_crossing[0]]:
                                                        grid_faces[i-1][j]['borders'][full_crossing[0]].append(face)
                                                    #Transversal case
                                                    if full_crossing[0] == 'west':
                                                        if type(edge.get_opp().get_face()) is GhostFace:
                                                            if edge.get_dest().get_pos()[0] <  (i-1)*square_len:
                                                                grid_faces[i-1][j]['borders']['south'].append(edge.get_opp().get_face())
                                                            elif edge.get_dest().get_pos()[0] >=  (i-1)*square_len+square_len:
                                                                grid_faces[i-1][j]['borders']['north'].append(edge.get_opp().get_face())
                                                        other_crossing = False
                                                        next_crossing = external_crossing([(i-1)*square_len, j*square_len], square_len, edge.get_next().get_dest().get_pos(), edge.get_next().get_opp().get_dest().get_pos())
                                                        if next_crossing is not None:
                                                            other_crossing = True
                                                        next_crossing = external_crossing([(i-1)*square_len, j*square_len], square_len, edge.get_next().get_next().get_dest().get_pos(), edge.get_next().get_next().get_opp().get_dest().get_pos())
                                                        if next_crossing is not None:
                                                            other_crossing = True
                                                        if not other_crossing:
                                                            if edge.get_dest().get_pos()[0] <  (i-1)*square_len:
                                                                grid_faces[i-1][j]['borders']['north'].append(face)
                                                            elif edge.get_dest().get_pos()[0] >=  (i-1)*square_len+square_len:
                                                                grid_faces[i-1][j]['borders']['south'].append(face)
                                                else:
                                                    if face not in grid_faces[i-1][j]['borders'][full_crossing[0]]:
                                                        grid_faces[i-1][j]['borders'][full_crossing[0]].append(face)
                                                    if face not in grid_faces[i-1][j]['borders'][full_crossing[1]]:
                                                        grid_faces[i-1][j]['borders'][full_crossing[1]].append(face)
                                                    if edge.get_opp not in grid_faces[i-1][j]['faces']:
                                                        expand_edge.append(edge.get_opp())
                                            edge = edge.get_next()
                                if expand_edge == [] and single_cell:
                                    for other_border in ['south', 'west', 'north']:
                                        grid_faces[i-1][j]['borders'][other_border].extend(grid_faces[i-1][j]['borders']['est'])
                                else:
                                    while expand_edge != []:
                                        edge = expand_edge[0]
                                        expand_edge = expand_edge[1:]
                                        crossing = external_crossing([(i-1)*square_len, j*square_len], square_len, edge.get_dest().get_pos(), edge.get_opp().get_dest().get_pos())
                                        if edge.get_face() not in grid_faces[i-1][j]['faces']:
                                            grid_faces[i-1][j]['faces'].append(edge.get_face())
                                        if edge.get_face() not in grid_faces[i-1][j]['borders'][crossing[0]]:
                                            grid_faces[i-1][j]['borders'][crossing[0]].append(edge.get_face())
                                        if edge.get_face() not in grid_faces[i-1][j]['borders'][crossing[1]]:
                                            grid_faces[i-1][j]['borders'][crossing[1]].append(edge.get_face())

                                        other_crossing = False
                                        if type(edge.get_face()) is MeshFace:
                                            next_crossing = external_crossing([(i-1)*square_len, j*square_len], square_len, edge.get_next().get_dest().get_pos(), edge.get_next().get_opp().get_dest().get_pos())
                                            if next_crossing is not None:
                                                expand_edge.append(edge.get_next().get_opp())
                                                other_crossing = True
                                            next_crossing = external_crossing([(i-1)*square_len, j*square_len], square_len, edge.get_next().get_next().get_dest().get_pos(), edge.get_next().get_next().get_opp().get_dest().get_pos())
                                            if next_crossing is not None:
                                                expand_edge.append(edge.get_next().get_next().get_opp())
                                                other_crossing = True
                                        if not other_crossing and (crossing == ['north', 'south'] or crossing == ['south', 'north']):
                                             grid_faces[i-1][j]['borders']['west'].append(edge.get_face())


                        elif border == 'south':
                            if grid_faces[i][j-1] is None:
                                grid_faces[i][j-1] = {'faces':[], 'borders': {'north':[], 'est':[], 'south':[],'west':[]}}
                                grid_faces[i][j-1]['faces'].extend(grid_faces[i][j]['borders'][border])
                                grid_faces[i][j-1]['borders']['north'].extend(grid_faces[i][j]['borders'][border])
                                single_cell = True
                                if len(grid_faces[i][j-1]['border']['north']) != 1:
                                    single_cell = False
                                expand_edge = []
                                for face in grid_faces[i][j-1]['borders']['north']:
                                    if type(face) is MeshFace:
                                        edge = face.get_edge()

                                        for k in [0,1,2]:
                                            full_crossing = external_crossing([i*square_len, (j-1)*square_len], square_len, edge.get_dest().get_pos(), edge.get_opp().get_dest().get_pos())
                                            if full_crossing is not None:
                                                single_cell = False
                                                if 'north' in full_crossing:
                                                    full_crossing.remove('north')
                                                    if face not in grid_faces[i][j-1]['borders'][full_crossing[0]]:
                                                        grid_faces[i][j-1]['borders'][full_crossing[0]].append(face)
                                                    #Transversal case
                                                    if full_crossing[0] == 'south':
                                                        if type(edge.get_opp().get_face()) is GhostFace:
                                                            if edge.get_dest().get_pos()[1] <  (j-1)*square_len:
                                                                grid_faces[i][j-1]['borders']['est'].append(edge.get_opp().get_face())
                                                            elif edge.get_dest().get_pos()[1] >=  (j-1)*square_len+square_len:
                                                                grid_faces[i][j-1]['borders']['west'].append(edge.get_opp().get_face())
                                                        other_crossing = False
                                                        next_crossing = external_crossing([i*square_len, (j-1)*square_len], square_len, edge.get_next().get_dest().get_pos(), edge.get_next().get_opp().get_dest().get_pos())
                                                        if next_crossing is not None:
                                                            other_crossing = True
                                                        next_crossing = external_crossing([i*square_len, (j-1)*square_len], square_len, edge.get_next().get_next().get_dest().get_pos(), edge.get_next().get_next().get_opp().get_dest().get_pos())
                                                        if next_crossing is not None:
                                                            other_crossing = True
                                                        if not other_crossing:
                                                            if edge.get_dest().get_pos()[1] < (j-1)*square_len:
                                                                grid_faces[i][j-1]['borders']['west'].append(face)
                                                            elif edge.get_dest().get_pos()[1] >=  (j-1)*square_len+square_len:
                                                                grid_faces[i][j-1]['borders']['est'].append(face)
                                                else:
                                                    if face not in grid_faces[i][j-1]['borders'][full_crossing[0]]:
                                                        grid_faces[i][j-1]['borders'][full_crossing[0]].append(face)
                                                    if face not in grid_faces[i][j-1]['borders'][full_crossing[1]]:
                                                        grid_faces[i][j-1]['borders'][full_crossing[1]].append(face)
                                                    if edge.get_opp not in grid_faces[i][j-1]['faces']:
                                                        expand_edge.append(edge.get_opp())
                                            edge = edge.get_next()
                                if expand_edge == [] and single_cell:
                                    for other_border in ['west', 'south', 'est']:
                                        grid_faces[i][j-1]['borders'][other_border].extend(grid_faces[i][j-1]['borders']['north'])
                                else:
                                    while expand_edge != []:
                                        edge = expand_edge[0]
                                        expand_edge = expand_edge[1:]
                                        crossing = external_crossing([i*square_len, (j-1)*square_len], square_len, edge.get_dest().get_pos(), edge.get_opp().get_dest().get_pos())
                                        if edge.get_face() not in grid_faces[i][j-1]['faces']:
                                            grid_faces[i][j-1]['faces'].append(edge.get_face())
                                        if edge.get_face() not in grid_faces[i][j-1]['borders'][crossing[0]]:
                                            grid_faces[i][j-1]['borders'][crossing[0]].append(edge.get_face())
                                        if edge.get_face() not in grid_faces[i][j-1]['borders'][crossing[1]]:
                                            grid_faces[i][j-1]['borders'][crossing[1]].append(edge.get_face())

                                        other_crossing = False
                                        if type(edge.get_face()) is MeshFace:
                                            next_crossing = external_crossing([i*square_len, (j-1)*square_len], square_len, edge.get_next().get_dest().get_pos(), edge.get_next().get_opp().get_dest().get_pos())
                                            if next_crossing is not None:
                                                expand_edge.append(edge.get_next().get_opp())
                                                other_crossing = True
                                            next_crossing = external_crossing([i*square_len, (j-1)*square_len], square_len, edge.get_next().get_next().get_dest().get_pos(), edge.get_next().get_next().get_opp().get_dest().get_pos())
                                            if next_crossing is not None:
                                                expand_edge.append(edge.get_next().get_next().get_opp())
                                                other_crossing = True
                                        if not other_crossing and (crossing == ['est', 'west'] or crossing == ['west', 'est']):
                                            grid_faces[i][j-1]['borders']['south'].append(edge.get_face())

    return grid_faces


if __name__ == '__main__':

    def aff_faces(faces):
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
        # pyplot.show()


    # a = MeshVertex(1,2)
    # b = MeshVertex(1.6,2.3)
    # c = MeshVertex(1,10)
    # print(clockwise_check(a,b,c))
    # x = MeshVertex(1.2, 3)
    # print(inside_triangle_check(a.get_pos(), c.get_pos(), b.get_pos(), x.get_pos()))
    #
    # print(cone_check([1,2],[3,0],[3,4],[8,10]))
    # print(cone_check([1,2],[3,0],[3,4],[0,2]))

    # tri = create_isolated_face(a,b,c)
    # aff_face([tri['faces'][0]])
    #
    #CAN HANDLE 3 * 2^n VERTICES
    #TODO: add signle point handling to handle any number of points
    # points_list=[[1,1],[2,3],[3,4],[4,7],[5,6],[6,3]]
    # points_list=[[1.,2.],[2.,1.],[2.2,8.],[5.,3.],[5.4,2.],[5.6,5.]]
    points_list=[[1.,2.],[1.6,2.3],[2.,1.],[2.2,8.],[2.7,1.8],[3.5,5.],[4.2,8.],[4.6,3.6],[5.,3.],[5.4,2.],[5.6,5.],[6.,4.2], [6.2, 5.8]]
    # points_list=[[4.2,4.8],[4.6,3.6],[5.,3.],[5.4,2.],[5.6,5.],[6.,4.2]]
    # points_list=[[1.,2],[2.,2.],[3.,8.],[4.,4.],[5.,1.],[6.,6.]]
    # points_list=[[1.2,6.],[1.6,4.8],[2.,1.],[2.5,8.],[2.7,1.8]]

    # points_list=[[1.,2.],[1.,4.],[1.,6.],[3.5,1.],[3.5,3.], [3.5,5.]]
    # points_list=[[1.2,2.],[1.2,5.],[1.2,7.]]

    # map=divide_and_conquer(points_list)
    # aff_face(map['faces'])

    translated_map = []
    for i in list(range(3)):
        for j in list(range(3)):
            for pnum in list(range(len(points_list))):
                translated_map.append([points_list[pnum][0]+i*8,points_list[pnum][1]+j*8])
    map = []
    for i in list(range(1)):
        for j in list(range(1)):
            for pnum in list(range(len(translated_map))):
                map.append([translated_map[pnum][0]+(i)*20,translated_map[pnum][1]+(j)*20])

    translated_map = sorted(map, key=lambda k: [k[0], k[1]])
    print(len(translated_map))

    print(datetime.now())
    map=divide_and_conquer(translated_map)
    print(datetime.now())

    # mapi=divide_and_conquer(points_list)
    # aff_faces(mapi['faces'])
    # pyplot.show()

    # map2 = divide_and_conquer([[4,4],[5,8],[6,2],[7,6],[8,3]])
    # aff_faces(map2['faces'])
    # pyplot.show()
    # edge_flip(map2['faces'][0].get_edge())
    # aff_faces(map2['faces'])
    # pyplot.show()

    #
    # aff_faces(map['faces'])
    #
    # for v in map['vertices']:
    #     if type(v) is MeshVertex:
    #         if v.get_edge() is None:
    #             print('map')
    #             print(v.get_pos())


    print(diamond_angle(1,0))
    print(diamond_angle(1,1))
    print(diamond_angle(0,1))
    print(diamond_angle(-1,1))
    print(diamond_angle(-1,0))
    print(diamond_angle(-1,-1))
    print(diamond_angle(0,-1))
    print(diamond_angle(1,-1))
    print(diamond_angle(1,-0.05))

    sq_len = 5

    grid = grid_sampling(sq_len, map['vertices'])
    print(datetime.now())
    # print(len(grid),len(grid[0]))

    # for row in grid:
    #     print(row)

    grid = fill_grid(grid, sq_len)
    print(datetime.now())

    # for row in grid:
    #     print(row)

    # p = MeshVertex(18,12)
    # f = find_face(p.get_pos(),grid[1][1])
    # # aff_faces([f])
    # # print(f)
    # print(datetime.now())
    # # print('ready for insertion')
    # inserted = vertex_insertion(p, f, map['faces'], map['vertices'])
    # print(datetime.now())
    #
    # # print('inserted')
    # print(inserted['vertices'][8].get_pos())
    # print(inserted['vertices'][2640].get_pos())
    #
    # # aff_faces(inserted['faces'])
    # sliced = edge_slice(inserted['vertices'][8], inserted['vertices'][2640], map['faces'], inserted['vertices'])
    # print(datetime.now())
    #
    # aff_faces(sliced['faces'])

    # aff_faces(map['faces'])

    # hull_points = [[17.2,18.4],[17.2,21.4],[19.2,21.4],[19.2,18.4]]
    hull_points = [[10.25,5.45],[10.25,18.45],[12.25,18.45],[12.25,5.45]]
    #hull_points = [[10.2,5.4],[10.2,18.4],[13.2,18.4],[13.2,5.4]]
    obstacle_map = add_solid(hull_points, map['faces'], map['vertices'], faces_grid=grid, square_len=sq_len)

    hull_points = [[15.25,12.45],[15.25,20.45],[16.25,20.45],[16.25,12.45]]

    obstacle_map = add_solid(hull_points, obstacle_map['faces'], obstacle_map['vertices'], faces_grid=grid, square_len=sq_len)

    aff_faces(obstacle_map['faces'])

    path = d_star([3.6,5.8],[19.6,16.6],obstacle_map['faces'],obstacle_map['vertices'], grid, sq_len)
    for v in path['vertices']:
        print(v.get_pos())


    # print(diamond_diff(0.4, 0.8))
    # print(diamond_diff(0.8, 1.2))
    # print(diamond_diff(1.2, 1.9))
    # print(diamond_diff(0.2, -0.2))
    # print(diamond_diff(-1.2, -1.9))
    # print(diamond_diff(0.4, -0.4))
    # print(diamond_diff(0.4, -1.4))
    # print(diamond_diff(1.4, -1.4))
    print(diamond_diff(0.8, -0.6))
    print(diamond_diff(0.8, -1.6))

    print('oyak')

    pyplot.show()
