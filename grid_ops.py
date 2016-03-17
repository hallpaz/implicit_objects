from DataStructures import BoundingBox, Vertex, GridCell, Triangle
import operator
import math
from geometry import *


def make_grid(box_volume, max_cells):

    min_corner = box_volume.min_corner
    max_corner = box_volume.max_corner
    greatest_dimension = max((map(operator.sub, max_corner, min_corner)))
    print("greatest_dimension: ", greatest_dimension)
    cell_dimension = greatest_dimension/max_cells
    print("cell_dimension: ", cell_dimension)

    default_value = max_corner[2] - min_corner[2]

    z = min_corner[2]
    grid = []
    while(z < max_corner[2]):
        y = min_corner[1]
        square = []
        while(y < max_corner[1]):
            x = min_corner[0]
            row = []
            while(x < max_corner[0]):
                row.append(Vertex(x, y, z, default_value))
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
    vertices = [
        grid[k+1][j+1][i],
        grid[k+1][j+1][i+1],
        grid[k][j+1][i+1],
        grid[k][j+1][i],
        grid[k+1][j][i],
        grid[k+1][j][i+1],
        grid[k][j][i+1],
        grid[k][j][i]
    ]
    return GridCell(vertices)

def possible_intersection(vertex, min_corner, cell_dimension):
    i = math.floor((vertex.x - min_corner[0]) / cell_dimension)
    j = math.floor((vertex.y - min_corner[1]) / cell_dimension)
    k = math.floor((vertex.z - min_corner[2]) / cell_dimension)
    return [(i, j, 0), (i, j+1, 0), (i, j-1, 0),
            (i+1, j, 0), (i+1, j+1, 0), (i+1, j-1, 0),
            (i-1, j, 0), (i-1, j+1, 0), (i-1, j-1, 0)
    ]


def assing_values(grid, vertexBuffer, indicesBuffer, min_corner, cell_dimension):

    # inter_points = []
    for indices in indicesBuffer:
        voxels_indices = []

        triangle = (vertexBuffer[indices.a], vertexBuffer[indices.b], vertexBuffer[indices.c])
        voxels_indices.extend(possible_intersection(triangle[0], min_corner, cell_dimension))
        voxels_indices.extend(possible_intersection(triangle[1], min_corner, cell_dimension))
        voxels_indices.extend(possible_intersection(triangle[2], min_corner, cell_dimension))
        # voxels_indices.append(indices_for_point(triangle[0], min_corner, cell_dimension))
        # voxels_indices.append(indices_for_point(triangle[1], min_corner, cell_dimension))
        # voxels_indices.append(indices_for_point(triangle[2], min_corner, cell_dimension))
        voxels_indices = set(voxels_indices)

        normal = compute_normal(triangle)
        for indices in voxels_indices:
            try:
                point = grid[0][indices[1]][indices[0]]
            except Exception as e:
                print(e, indices)
                continue

            intersection_point = triangle_intersection(point, point + Vertex(0, 0, 7), triangle, normal)
            if intersection_point is not None:

                #d1 = intersection_point.z - point.z
                #d2 = Vertex.distance(intersection_point, point)
                #print(d1, d2)
                #exit()
                # inter_points.append("%f %f %f\n"%(intersection_point.x,intersection_point.y,intersection_point.z))
                point.value = intersection_point.z - point.z
                #sign = 1 if (point - intersection_point) * normal > 0 else -1
                # new_value = intersection_point.z - point.z
                # if new_value < point.value:
                #     point.value = new_value
#     inter_file = open("inter_points.ply","w")
#     inter_file.write('''ply
# format ascii 1.0
# element vertex %d
# property float x
# property float y
# property float z
# end_header
# %s
# '''%(len(inter_points),"".join(inter_points)))
#     inter_file.close()
#     print("all intersections computed")
#     exit()
    for k in range(1, len(grid)):
        for j in range(len(grid[0])):
            for i in range(len(grid[0][0])):
                point = grid[k][j][i]
                ref_point = grid[0][j][i]
                # floating point arithmetics might be a source of error here
                point.value = ref_point.value - (ref_point.z - point.z)

    print("all grid points values computed")

    return grid
