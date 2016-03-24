from DataStructures import BoundingBox, Vertex, GridCell, Triangle
import operator
import math
from geometry import *


SENTINEL_VALUE = 72347

def make_grid(box_volume, max_cells):

    min_corner = box_volume.min_corner
    max_corner = box_volume.max_corner
    greatest_dimension = max((map(operator.sub, max_corner, min_corner)))
    print("greatest_dimension: ", greatest_dimension)
    cell_dimension = greatest_dimension/max_cells
    print("cell_dimension: ", cell_dimension)
    print("min_corner:{}, max_corner:{}".format(min_corner, max_corner))

    #default_value = max_corner[2] - min_corner[2]

    z = min_corner[2]
    grid = []
    while(z < max_corner[2]):
        y = min_corner[1]
        square = []
        while(y < max_corner[1]):
            x = min_corner[0]
            row = []
            while(x < max_corner[0]):
                row.append(Vertex(x, y, z, SENTINEL_VALUE))
                x += cell_dimension
            square.append(row)
            y += cell_dimension
        grid.append(square)
        z += cell_dimension

    print("Total cubes: i: {}, j: {}, k: {}".format(len(grid[0][0]), len(grid[0]), len(grid)))
    return grid, cell_dimension


def indices_for_point(vertex, min_corner, cell_dimension):
    i = math.floor((vertex.x - min_corner[0]) / cell_dimension)
    j = math.floor((vertex.y - min_corner[1]) / cell_dimension)
    k = math.floor((vertex.z - min_corner[2]) / cell_dimension)

    return (k, j, i)

def cell_for_indices(grid, indices):
    #TODO: inverter para padronizar
    k = indices[0]
    j = indices[1]
    i = indices[2]
    # vertices = [
    #     grid[k+1][j+1][i],
    #     grid[k+1][j+1][i+1],
    #     grid[k][j+1][i+1],
    #     grid[k][j+1][i],
    #     grid[k+1][j][i],
    #     grid[k+1][j][i+1],
    #     grid[k][j][i+1],
    #     grid[k][j][i]
    # ]
    vertices = [
        grid[k+1][j][i],
        grid[k+1][j][i+1],
        grid[k][j][i+1],
        grid[k][j][i],

        grid[k+1][j+1][i],
        grid[k+1][j+1][i+1],
        grid[k][j+1][i+1],
        grid[k][j+1][i]
    ]
    return GridCell(vertices)

def possible_intersection(vertex, min_corner, cell_dimension):
    i = math.floor((vertex.x - min_corner[0]) / cell_dimension)
    j = math.floor((vertex.y - min_corner[1]) / cell_dimension)
    k = math.floor((vertex.z - min_corner[2]) / cell_dimension)
    return [(i, j), (i, j+1), (i, j-1),
            (i+1, j), (i+1, j+1), (i+1, j-1),
            (i-1, j), (i-1, j+1), (i-1, j-1),

            (i, j+2), (i+1, j+2), (i+1, j-2),
            (i, j-2), (i-1, j+2), (i-1, j-2),

            (i+2, j), (i+2, j+1), (i+2, j-1),
            (i-2, j), (i-2, j+1), (i-2, j-1),

            (i+2, j+2), (i+2, j-2),
            (i-2, j+2), (i-2, j-2)
            ]


def assing_values(grid, vertexBuffer, indicesBuffer, min_corner, cell_dimension):

    inter_points = []
    chosen_pairs = set([])
    for face in indicesBuffer:
        voxels_indices = []
        #print(face)
        triangle = (vertexBuffer[face.a], vertexBuffer[face.b], vertexBuffer[face.c])
        voxels_indices.extend(possible_intersection(triangle[0], min_corner, cell_dimension))
        voxels_indices.extend(possible_intersection(triangle[1], min_corner, cell_dimension))
        voxels_indices.extend(possible_intersection(triangle[2], min_corner, cell_dimension))

        voxels_indices = set(voxels_indices)

        #chosen_pairs = chosen_pairs.union(voxels_indices)

        normal = compute_normal(triangle)
        for indices in voxels_indices:
            try:
                point = grid[0][indices[1]][indices[0]]
            except Exception as e:
                continue

            #already found intersection
            if indices in chosen_pairs:
                continue

            intersection_point = moller_triangle_intersection(point, point + Vertex(0, 0, 7), triangle, normal)
            #intersection_point = ray_plane_intersection(point, point + Vertex(0, 0, 7), triangle[0], normal)
            if intersection_point is not None:
                chosen_pairs.add(indices)
                #d1 = intersection_point.z - point.z
                #d2 = Vertex.distance(intersection_point, point)
                #print(d1, d2)

                inter_points.append("%f %f %f\n"%(intersection_point.x,intersection_point.y,intersection_point.z))
                point.value = (intersection_point.z - point.z)
                #print("intersection:", intersection_point)
                for k in range(1, len(grid)):
                    p = grid[k][indices[1]][indices[0]]
                    p.value = (intersection_point.z - p.z)
                    #break
    inter_file = open("inter_points_moller.ply","w")
    inter_file.write('''ply
format ascii 1.0
element vertex %d
property float x
property float y
property float z
end_header
%s
'''%(len(inter_points),"".join(inter_points)))
    inter_file.close()

    print("all grid points values computed")

    remaining_pairs = [(j, i) for j in range(len(grid[0])) for i in range(len(grid[0][0])) if (j, i) not in chosen_pairs]
    print(remaining_pairs)

    return grid


def assing_all_values(grid, vertexBuffer, indicesBuffer, min_corner, cell_dimension):

    inter_points = []
    for j in range(len(grid[0])):
        for i in range(len(grid[0][0])):
            point = grid[0][j][i]
            for face in indicesBuffer:
                triangle = (vertexBuffer[face.a], vertexBuffer[face.b], vertexBuffer[face.c])
                normal = compute_normal(triangle)

                intersection_point = triangle_intersection(point, point + Vertex(0, 0, 7), triangle, normal)
                if intersection_point is not None:
                    point.value = (intersection_point.z - point.z)
                    #print(j, i, point.value)
                    for k in range(1, len(grid)):
                        p = grid[k][j][i]
                        #p.value = (intersection_point.z - p.z)
                        p.value = point.value - k*cell_dimension
                        #print(p.value)
                    inter_points.append("%f %f %f\n"%(intersection_point.x,intersection_point.y,intersection_point.z))
                    #exit()
                    break
    inter_file = open("inter_points_all.ply","w")
    inter_file.write('''ply
format ascii 1.0
element vertex %d
property float x
property float y
property float z
end_header
%s
'''%(len(inter_points),"".join(inter_points)))
    inter_file.close()
    print("inter points written")


    # for k in range(1, len(grid)):
    #     for j in range(len(grid[0])):
    #         for i in range(len(grid[0][0])):
    #             point = grid[k][j][i]
    #             ref_point = grid[0][j][i]
    #             # floating point arithmetics might be a source of error here
    #             point.value = ref_point.value - k*cell_dimension#- (ref_point.z - point.z)
    print("all grid points values computed")

    return grid
