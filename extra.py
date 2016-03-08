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
