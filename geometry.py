import operator
import math
from LookUpTables import edgeTable, triTable
from DataStructures import BoundingBox, GridCell, Vertex, Triangle


def compute_normal(triangle_vertices):
    v1 = triangle_vertices[1] - triangle_vertices[0]
    v2 = triangle_vertices[2] - triangle_vertices[1]
    normal = Vertex.cross(v1, v2).normalize()
    return normal

def segment_plane_intersection(vertex1, vertex2, normal, point):
        p = vertex2 - vertex1
        denominator = normal * p
        if abs(denominator) < 0.00001:
            print("xabu", vertex1, vertex2, normal)
            return None
        r = normal * (point - vertex1)
        r /= denominator
        if r < 0 or r > 1:
            print("segment problem")
            return None
        return r

def ray_plane_intersection(vertex1, vertex2, normal, point):
    p = vertex2 - vertex1
    denominator = normal * p
    if abs(denominator) < 0.0001:
        return None
    r = normal * (point - vertex1)
    r /= denominator
    # if r < 0 or r > 1:
    #     return None
    intersection_point = vertex1 + p.scalar_mult(r)
    return intersection_point

def triangle_intersection(v1, v2, triangle_vertices, normal = None):
    if normal is None:
        normal = compute_normal(triangle_vertices)
    t = segment_plane_intersection(v1, v2, normal, triangle_vertices[0])
    if t is None:
        return None

    #scalar vector multiplication
    intersection_point = v1 + (v2-v1).scalar_mult(t)
    u = triangle_vertices[1] - triangle_vertices[0]
    v = triangle_vertices[2] - triangle_vertices[1]
    uu = u*u
    uv = u*v
    vv = v*v
    w = intersection_point - triangle_vertices[0]
    wu = w*u
    wv = w*v
    D = uv * uv - uu * vv

    # get and test parametric coords
    s = (uv * wv - vv * wu) / D
    if (s < 0.0 or s > 1.0):         # I is outside T
        return None
    t = (uv * wu - uu * wv) / D
    if (t < 0.0 or (s + t) > 1.0):  # I is outside T
        return None

    return intersection_point                       # I is in T

def fake_triangle_intersection(v1, v2, triangle_vertices, normal = None):
    if normal is None:
        normal = compute_normal(triangle_vertices)
    t = segment_plane_intersection(v1, v2, normal, triangle_vertices[0])
    if t is None:
        return None

    #scalar vector multiplication
    intersection_point = v1 + (v2-v1).scalar_mult(t)

    q0 = intersection_point - triangle_vertices[0]
    q1 = triangle_vertices[1] - triangle_vertices[0]
    q2 = triangle_vertices[2] - triangle_vertices[0]


    return intersection_point                       # I is in T

def moller_triangle_intersection(v1, v2, triangle_vertices, normal = None):
    EPSILON = 0.000001

    e1 = triangle_vertices[1] - triangle_vertices[0]
    e2 = triangle_vertices[2] - triangle_vertices[0]

    D = v2 - v1
    P = Vertex.cross(D, e2)

    det = e1 * P
    if abs(det) < EPSILON:
        return None

    inv_det = 1.0/det
    T = v1 - triangle_vertices[0]
    u = (T*P) * inv_det

    if(u < 0.0 or u > 1.0):
        return None

    Q = Vertex.cross(T, e1)
    v = (D*Q) * inv_det
    if(v < 0.0 or (u + v)  > 1.0):
        return None

    t = (e2*Q) * inv_det

    if(t > EPSILON):
        intersection_point = v1 + (v2-v1).scalar_mult(t)
        return intersection_point

    return None



def compute_side(vertex, triangle_vertices):
    normal = compute_normal(triangle)
    if vertex * normal > 0:
        return 1
    return -1
