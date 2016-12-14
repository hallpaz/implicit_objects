from DataStructures import BoundingBox, Vertex, GridCell, Triangle
import operator
import math
from geometry import *


SENTINEL_VALUE = 72347
BIG = 9999999

def write_point_cloud(vertices, filename):
    points = [str(v) for v in vertices]
    with open(filename, "w") as pcfile:
        pcfile.write(
        '''ply
        format ascii 1.0
        element vertex %d
        property float x
        property float y
        property float z
        end_header
        %s
        '''%(len(points),"".join(points)))

def read_OFF(off_file):
    vertexBuffer = []
    indexBuffer = []
    with open(off_file, "r") as modelfile:
        first = modelfile.readline().strip()
        if first != "OFF":
            raise(Exception("not a valid OFF file ({})".format(first)))

        parameters = modelfile.readline().strip().split()
        min_distance = BIG
        minx = BIG
        miny = BIG
        minz = BIG
        maxx = -BIG
        maxy = -BIG
        maxz = -BIG

        if len(parameters) < 2:
            raise(Exception("OFF file has invalid number of parameters"))

        for i in range(int(parameters[0])):
            coordinates = modelfile.readline().split()
            X, Y, Z = float(coordinates[0]), float(coordinates[1]), float(coordinates[2])
            vertexBuffer.append(Vertex(X, Y, Z))
            if X < minx:
                minx = X
            if X > maxx:
                maxx = X
            if Y < miny:
                miny = Y
            if Y > maxy:
                maxy = Y
            if Z < minz:
                minz = Z
            if Z > maxz:
                maxz = Z

        for i in range(int(parameters[1])):
            indices = modelfile.readline().split()
            indexBuffer.append(Triangle(int(indices[1]), int(indices[2]), int(indices[3])))

        bounding_box = BoundingBox(minx, miny, minz-1, maxx, maxy, maxz+1)

    return vertexBuffer, indexBuffer, bounding_box

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


def assing_values(grid, vertexBuffer, indicesBuffer, min_corner, cell_dimension, filename = None):

    chosen_pairs = set([])
    inter_points = []
    for face in indicesBuffer:
        voxels_indices = []
        triangle = (vertexBuffer[face.a], vertexBuffer[face.b], vertexBuffer[face.c])
        voxels_indices.extend(possible_intersection(triangle[0], min_corner, cell_dimension))
        voxels_indices.extend(possible_intersection(triangle[1], min_corner, cell_dimension))
        voxels_indices.extend(possible_intersection(triangle[2], min_corner, cell_dimension))

        voxels_indices = set(voxels_indices)

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
            if filename is not None:
                inter_points.append(intersection_point)
            if intersection_point is not None:
                chosen_pairs.add(indices)

                point.value = (intersection_point.z - point.z)
                #print("intersection:", intersection_point)
                for k in range(1, len(grid)):
                    p = grid[k][indices[1]][indices[0]]
                    p.value = (intersection_point.z - p.z)
                    #break

    if filename is not None:
        write_point_cloud(inter_points, filename)

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
