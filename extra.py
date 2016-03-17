# def signed_distance(cell, vertices, faces):
#     edges = [
#         (0, 1), (1, 2), (2, 3), (3, 0),
#         (4, 5), (5, 6), (6, 7), (7, 4),
#         (2, 6), (1, 5), (3, 7), (0, 4)
#     ]
#     hasintersection = False
#     for edge in edges:
#         for indices in faces:
#             triangle = (vertices[indices.a], vertices[indices.b], vertices[indices.c])
#             normal = compute_normal(triangle)
#             intersection_point = triangle_intersection(cell.positions[edge[0]], cell.positions[edge[1]], triangle, normal)
#             if intersection_point is not None:
#                 hasintersection = True
#                 #sign = compute_side(cell.positions[edge[0]], triangle)
#                 sign = 1 if cell.positions[edge[0]] * normal > 0 else -1
#                 cell.values[edge[0]] = cell.weights[edge[0]]*cell.values[edge[0]] + sign * Vertex.distance(cell.positions[edge[0]], intersection_point)
#                 cell.weights[edge[0]] += 1
#                 cell.values[edge[0]] /= cell.weights[edge[0]]
#
#                 cell.values[edge[1]] = cell.weights[edge[1]]*cell.values[edge[1]] -sign * Vertex.distance(cell.positions[edge[1]], intersection_point)
#                 cell.weights[edge[1]] += 1
#                 cell.values[edge[1]] /= cell.weights[edge[1]]
#     return hasintersection


def assign_distance_to_grid(grid, vertices, faces):
    flatgrid = []
    # for i in range(len(grid)):
    #     for j in range(len(grid[0])):
    #         for k in range(len(grid[0][0])):
    #             cell = grid[i][j][k]
    #             if signed_distance(cell, vertices, faces):
    #                 flatgrid.append(cell)
    #             print(i, j, k)

    print("all cells with value")
    return flatgrid


def polygonise_cube(cell, isovalue, vertices, triangles):
    vertlist = 12*[None]
    cubeindex = 0
    if (cell.values[0] < isovalue): cubeindex |= 1
    if (cell.values[1] < isovalue): cubeindex |= 2
    if (cell.values[2] < isovalue): cubeindex |= 4
    if (cell.values[3] < isovalue): cubeindex |= 8
    if (cell.values[4] < isovalue): cubeindex |= 16
    if (cell.values[5] < isovalue): cubeindex |= 32
    if (cell.values[6] < isovalue): cubeindex |= 64
    if (cell.values[7] < isovalue): cubeindex |= 128

    # Cube is entirely in/out of the surface
    if (edgeTable[cubeindex] == 0):
       return(0)

    # Find the vertices where the surface intersects the cube
    if (edgeTable[cubeindex] & 1):
       vertlist[0] = vertex_interpolation(isovalue,cell.positions[0],cell.positions[1],cell.values[0],cell.values[1])

    if (edgeTable[cubeindex] & 2):
       vertlist[1] = vertex_interpolation(isovalue,cell.positions[1],cell.positions[2],cell.values[1],cell.values[2])

    if (edgeTable[cubeindex] & 4):
       vertlist[2] = vertex_interpolation(isovalue,cell.positions[2],cell.positions[3],cell.values[2],cell.values[3])

    if (edgeTable[cubeindex] & 8):
       vertlist[3] = vertex_interpolation(isovalue,cell.positions[3],cell.positions[0],cell.values[3],cell.values[0])

    if (edgeTable[cubeindex] & 16):
       vertlist[4] = vertex_interpolation(isovalue,cell.positions[4],cell.positions[5],cell.values[4],cell.values[5])

    if (edgeTable[cubeindex] & 32):
       vertlist[5] = vertex_interpolation(isovalue,cell.positions[5],cell.positions[6],cell.values[5],cell.values[6])

    if (edgeTable[cubeindex] & 64):
       vertlist[6] = vertex_interpolation(isovalue,cell.positions[6],cell.positions[7],cell.values[6],cell.values[7])

    if (edgeTable[cubeindex] & 128):
       vertlist[7] = vertex_interpolation(isovalue,cell.positions[7],cell.positions[4],cell.values[7],cell.values[4])

    if (edgeTable[cubeindex] & 256):
       vertlist[8] = vertex_interpolation(isovalue,cell.positions[0],cell.positions[4],cell.values[0],cell.values[4])

    if (edgeTable[cubeindex] & 512):
       vertlist[9] = vertex_interpolation(isovalue,cell.positions[1],cell.positions[5],cell.values[1],cell.values[5])

    if (edgeTable[cubeindex] & 1024):
       vertlist[10] = vertex_interpolation(isovalue,cell.positions[2],cell.positions[6],cell.values[2],cell.values[6])

    if (edgeTable[cubeindex] & 2048):
       vertlist[11] = vertex_interpolation(isovalue,cell.positions[3],cell.positions[7],cell.values[3],cell.values[7])


    # Create the triangles
    i = 0
    triangles_count = 0
    while triTable[cubeindex][i] != -1:
        vertices.append(vertlist[triTable[cubeindex][i] ])
        vertices.append(vertlist[triTable[cubeindex][i+1]])
        vertices.append(vertlist[triTable[cubeindex][i+2]])
        triangles_count += 1
        i += 3
        triangles.append(Triangle(len(triangles)*3 + 2, len(triangles)*3 + 1, len(triangles)*3))


    return triangles_count
