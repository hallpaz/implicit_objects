import operator
import math
import grid_ops
from LookUpTables import edgeTable, triTable
from DataStructures import BoundingBox, GridCell, Vertex, Triangle
from naive_mesh import generateOFF, generatePlaneOFF



def makeGrid(box_volume, max_cells):
    min_corner = box_volume.min_corner
    max_corner = box_volume.max_corner
    greatest_dimension = max((map(operator.sub, max_corner, min_corner)))
    print("greatest_dimension: ", greatest_dimension)
    cell_dimension = greatest_dimension/max_cells
    print("cell_dimension: ", cell_dimension)

    x = min_corner[0]
    volume_of_cubes = []
    while x < max_corner[0]:
        square_of_cubes = []
        y = min_corner[1]
        while y < max_corner[1]:
            row_of_cubes = []
            z = min_corner[2]
            while z < max_corner[2]:
                vertices_positions = []
                for u in range(8):
                    vertices_positions.append(Vertex(x, y, z))
                    vertices_positions.append(Vertex(x + cell_dimension, y, z))
                    vertices_positions.append(Vertex(x + cell_dimension, y + cell_dimension, z))
                    vertices_positions.append(Vertex(x, y + cell_dimension, z))
                    vertices_positions.append(Vertex(x, y, z + cell_dimension))
                    vertices_positions.append(Vertex(x + cell_dimension, y, z + cell_dimension))
                    vertices_positions.append(Vertex(x + cell_dimension, y + cell_dimension, z + cell_dimension))
                    vertices_positions.append(Vertex(x, y + cell_dimension, z + cell_dimension))
                row_of_cubes.append(GridCell(vertices_positions))
                z += cell_dimension
            square_of_cubes.append(row_of_cubes)
            y += cell_dimension
        volume_of_cubes.append(square_of_cubes)
        x += cell_dimension

    return volume_of_cubes, cell_dimension


def assign_values(grid, function):
    for k in range(len(grid)-1):
        for j in range(len(grid[0])-1):
            for i in range(len(grid[0][0])-1):
                cell = grid_ops.cell_for_indices(grid, (k, j, i))
                for u in range(8):
                    value = function(cell.vertices[u].x, cell.vertices[u].y, cell.vertices[u].z)
                    cell.vertices[u].value = value



# special case when s1 <= isovalue <= s2 or s2 <= isovalue <= s1
def vertex_interpolation(isovalue, p1, p2, s1, s2):
    if abs(isovalue - s1) < 0.00001:
        return p1
    if abs(isovalue - s2) < 0.00001:
        return p2
    if abs(s1 - s2) < 0.00001:
        return p1
    alfa = (isovalue - s1) / (s2 - s1)
    r = Vertex(0, 0, 0)

    r.x = (1-alfa)*p1.x + alfa*p2.x
    r.y = (1-alfa)*p1.y + alfa*p2.y
    r.z = (1-alfa)*p1.z + alfa*p2.z
    return r


def compute_normal(triangle_vertices):
    v1 = triangle_vertices[1] - triangle_vertices[0]
    v2 = triangle_vertices[2] - triangle_vertices[1]
    normal = Vertex.cross(v1, v2).normalize()
    return normal

def segment_plane_intersection(vertex1, vertex2, normal, point):
        p = vertex2 - vertex1
        denominator = normal * p
        if denominator < 0.0001:
            return None
        r = normal * (point - vertex1)
        r /= denominator
        if r < 0 or r > 1:
            return None
        return r

def ray_plane_intersection(vertex1, vertex2, normal, point):
    p = vertex2 - vertex1
    denominator = normal * p
    if denominator < 0.0001:
        return None
    r = normal * (point - vertex1)
    r /= denominator
    # if r < 0 or r > 1:
    #     return None
    return r

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


def compute_side(vertex, triangle_vertices):
    normal = compute_normal(triangle)
    if vertex * normal > 0:
        return 1
    return -1

def nearby_voxels(vertex, min_corner, cell_dimension):
    i = math.floor((vertex.x - min_corner[0]) / cell_dimension)
    j = math.floor((vertex.y - min_corner[1]) / cell_dimension)
    k = math.floor((vertex.z - min_corner[2]) / cell_dimension)
    return [(i, j, k), (i, j+1, k), (i, j-1, k), (i, j, k+1), (i, j, k-1), (i, j+1, k+1), (i, j-1, k-1), (i, j+1, k-1), (i, j-1, k+1),
            (i+1, j, k), (i+1, j+1, k), (i+1, j-1, k), (i+1, j, k+1), (i+1, j, k-1), (i+1, j+1, k+1), (i+1, j-1, k-1), (i+1, j+1, k-1), (i+1, j-1, k+1),
            (i-1, j, k), (i-1, j+1, k), (i-1, j-1, k), (i-1, j, k+1), (i-1, j, k-1), (i-1, j+1, k-1), (i-1, j-1, k-1), (i-1, j+1, k-1), (i-1, j-1, k+1)
    ]


def signed_distance(grid, vertices, faces, min_corner, cell_dimension):
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (2, 6), (1, 5), (3, 7), (0, 4)
    ]

    zedges = [(3, 0), (1, 2), (5, 6), (7, 4)]
    yedges = [(3, 7), (2, 6), (1, 5), (0, 4)]
    xedges = [(2, 3), (6, 7), (0, 1), (4, 5)]

    hasintersection = False
    flatgrid = []
    processed_cells = 0
    for indices in faces:
        voxels_indices = []
        triangle = (vertices[indices.a], vertices[indices.b], vertices[indices.c])
        voxels_indices.extend(nearby_voxels(triangle[0], min_corner, cell_dimension))
        voxels_indices.extend(nearby_voxels(triangle[1], min_corner, cell_dimension))
        voxels_indices.extend(nearby_voxels(triangle[2], min_corner, cell_dimension))
        voxels_indices = set(voxels_indices)
        for indices in voxels_indices:
            #print(indices[0], indices[1], indices[2])
            try:
                cell = grid[indices[0]][indices[1]][indices[2]]
            except Exception:
                with open("mesh_log.txt", "a+") as logfile:
                    logfile.write("Out of bounds: {}\n".format(indices))
                continue
            normal = compute_normal(triangle)
            for edge in edges:
                intersection_point = triangle_intersection(cell.positions[edge[0]], cell.positions[edge[1]], triangle, normal)
                if intersection_point is not None:
                    hasintersection = True
                    # Compute Vertex Side
                    sign = 1 if (cell.positions[edge[0]] - intersection_point) * normal > 0 else -1
                    cell.values[edge[0]] = cell.weights[edge[0]]*cell.values[edge[0]] + sign * Vertex.distance(cell.positions[edge[0]], intersection_point)
                    cell.weights[edge[0]] += 1
                    cell.values[edge[0]] /= cell.weights[edge[0]]

                    cell.values[edge[1]] = cell.weights[edge[1]]*cell.values[edge[1]] -sign * Vertex.distance(cell.positions[edge[1]], intersection_point)
                    cell.weights[edge[1]] += 1
                    cell.values[edge[1]] /= cell.weights[edge[1]]
                    flatgrid.append(cell)
                    processed_cells += 1
                    if processed_cells % 200 == 0:
                        print('cells on flat grid: ' ,processed_cells)

    flatgrid = set(flatgrid)
    return flatgrid

def polygonise_cube(cell, isovalue, vertices, triangles):
    vertlist = 12*[None]
    cubeindex = 0
    if (cell.vertices[0].value > isovalue): cubeindex |= 1
    if (cell.vertices[1].value > isovalue): cubeindex |= 2
    if (cell.vertices[2].value > isovalue): cubeindex |= 4
    if (cell.vertices[3].value > isovalue): cubeindex |= 8
    if (cell.vertices[4].value > isovalue): cubeindex |= 16
    if (cell.vertices[5].value > isovalue): cubeindex |= 32
    if (cell.vertices[6].value > isovalue): cubeindex |= 64
    if (cell.vertices[7].value > isovalue): cubeindex |= 128

    # Cube is entirely in/out of the surface
    if (edgeTable[cubeindex] == 0):
       return(0)

    # Find the vertices where the surface intersects the cube
    if (edgeTable[cubeindex] & 1):
       vertlist[0] = vertex_interpolation(isovalue,cell.vertices[0],cell.vertices[1],cell.vertices[0].value,cell.vertices[1].value)

    if (edgeTable[cubeindex] & 2):
       vertlist[1] = vertex_interpolation(isovalue,cell.vertices[1],cell.vertices[2],cell.vertices[1].value,cell.vertices[2].value)

    if (edgeTable[cubeindex] & 4):
       vertlist[2] = vertex_interpolation(isovalue,cell.vertices[2],cell.vertices[3],cell.vertices[2].value,cell.vertices[3].value)

    if (edgeTable[cubeindex] & 8):
       vertlist[3] = vertex_interpolation(isovalue,cell.vertices[3],cell.vertices[0],cell.vertices[3].value,cell.vertices[0].value)

    if (edgeTable[cubeindex] & 16):
       vertlist[4] = vertex_interpolation(isovalue,cell.vertices[4],cell.vertices[5],cell.vertices[4].value,cell.vertices[5].value)

    if (edgeTable[cubeindex] & 32):
       vertlist[5] = vertex_interpolation(isovalue,cell.vertices[5],cell.vertices[6],cell.vertices[5].value,cell.vertices[6].value)

    if (edgeTable[cubeindex] & 64):
       vertlist[6] = vertex_interpolation(isovalue,cell.vertices[6],cell.vertices[7],cell.vertices[6].value,cell.vertices[7].value)

    if (edgeTable[cubeindex] & 128):
       vertlist[7] = vertex_interpolation(isovalue,cell.vertices[7],cell.vertices[4],cell.vertices[7].value,cell.vertices[4].value)

    if (edgeTable[cubeindex] & 256):
       vertlist[8] = vertex_interpolation(isovalue,cell.vertices[0],cell.vertices[4],cell.vertices[0].value,cell.vertices[4].value)

    if (edgeTable[cubeindex] & 512):
       vertlist[9] = vertex_interpolation(isovalue,cell.vertices[1],cell.vertices[5],cell.vertices[1].value,cell.vertices[5].value)

    if (edgeTable[cubeindex] & 1024):
       vertlist[10] = vertex_interpolation(isovalue,cell.vertices[2],cell.vertices[6],cell.vertices[2].value,cell.vertices[6].value)

    if (edgeTable[cubeindex] & 2048):
       vertlist[11] = vertex_interpolation(isovalue,cell.vertices[3],cell.vertices[7],cell.vertices[3].value,cell.vertices[7].value)


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


def sphere(x, y, z):
    return x**2 + y**2 + z**2

def bitorus(x, y, z):
    return ((x**2 + y**2)**2 - x**2 + y**2)**2 + z**2

def coeur(x, y, z):
    return (x**2 + (9/4)*y**2 + z**2 - 1)**3 - x**2 * z**3 - (9/80)*y**2 * z**3



def main():
    vertexBuffer, indexBuffer, box_volume = generateOFF("images/vase_rgb.jpg", "images/vase_depth.png", "models/myvase.off")

    grid, dimension = grid_ops.make_grid(box_volume, 400)

    grid_ops.assing_values(grid, vertexBuffer, indexBuffer, box_volume.min_corner, dimension)
    #grid_ops.assing_all_values(grid, vertexBuffer, indexBuffer, box_volume.min_corner, dimension)

    vertices = []
    triangles = []

    # for k in range(len(grid)-1):
    #     print(grid[k][4][3].value)

    flag = False
    for k in range(len(grid)-1):
        print("K:", k)
        for j in range(len(grid[0])-1):
            for i in range(len(grid[0][0])-1):
                cell = grid_ops.cell_for_indices(grid, (k, j, i))
                for vertex in cell.vertices:
                    if vertex.value == grid_ops.SENTINEL_VALUE:
                        flag = True
                        break
                if(flag):
                    flag = False
                    continue
                polygonise_cube(cell, 0, vertices, triangles)

    points = [str(i) for i in vertices]
    indices = [str(i) for i in triangles]

    meshfile = open("vase200_moller.off","w")
    meshfile.write(
    '''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    meshfile.close()

def plane_test():
    vertexBuffer, indexBuffer, box_volume = generatePlaneOFF("models/plane.off")

    grid, dimension = grid_ops.make_grid(box_volume, 17)

    grid_ops.assing_values(grid, vertexBuffer, indexBuffer, box_volume.min_corner, dimension)
    #grid_ops.assing_all_values(grid, vertexBuffer, indexBuffer, box_volume.min_corner, dimension)

    vertices = []
    triangles = []

    flag = False
    for k in range(len(grid)-1):
        print("K:", k)
        for j in range(len(grid[0])-1):
            for i in range(len(grid[0][0])-1):
                cell = grid_ops.cell_for_indices(grid, (k, j, i))
                for vertex in cell.vertices:
                    if vertex.value == grid_ops.SENTINEL_VALUE:
                        flag = True
                        break
                if(flag):
                    flag = False
                    continue
                polygonise_cube(cell, 0, vertices, triangles)

    points = [str(i) for i in vertices]
    indices = [str(i) for i in triangles]

    meshfile = open("rec_plane_17moller.off","w")
    meshfile.write(
    '''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    meshfile.close()

def poly_function():
    half_size = 3
    box_volume = BoundingBox(-half_size, -half_size, -half_size, half_size, half_size, half_size)

    grid, dimension = grid_ops.make_grid(box_volume, 50)

    assign_values(grid, bitorus)

    vertices = []
    triangles = []
    for k in range(len(grid)-1):
        print("K:", k)
        for j in range(len(grid[0])-1):
            for i in range(len(grid[0][0])-1):
                cell = grid_ops.cell_for_indices(grid, (k, j, i))
                polygonise_cube(cell, 0.01, vertices, triangles)

    points = [str(i) for i in vertices]
    indices = [str(i) for i in triangles]

    meshfile = open("torus100_new.off","w")
    meshfile.write(
    '''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    meshfile.close()
    vertexBuffer, indexBuffer, box_volume = generatePlaneOFF("models/plane.off")

    grid, dimension = grid_ops.make_grid(box_volume, 10)

    grid_ops.assing_values(grid, vertexBuffer, indexBuffer, box_volume.min_corner, dimension)
    #grid_ops.assing_all_values(grid, vertexBuffer, indexBuffer, box_volume.min_corner, dimension)

    vertices = []
    triangles = []

    for k in range(len(grid)-1):
        print("K:", k)
        for j in range(len(grid[0])-1):
            for i in range(len(grid[0][0])-1):
                cell = grid_ops.cell_for_indices(grid, (k, j, i))
                polygonise_cube(cell, 0, vertices, triangles)

    points = [str(i) for i in vertices]
    indices = [str(i) for i in triangles]

    meshfile = open("rec_plane_10.off","w")
    meshfile.write(
    '''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    meshfile.close()

def artificial_valoration():
    vertexBuffer, indexBuffer, box_volume = generatePlaneOFF("models/plane.off")

    grid, dimension = grid_ops.make_grid(box_volume, 100)

    # grid_ops.assing_values(grid, vertexBuffer, indexBuffer, box_volume.min_corner, dimension)
    #grid_ops.assing_all_values(grid, vertexBuffer, indexBuffer, box_volume.min_corner, dimension)

    for j in range(len(grid[0])):
        for i in range(len(grid[0][0])):
            point = grid[0][j][i]
            point.value = 1.0
            for k in range(1, len(grid)):
                p = grid[k][j][i]
                p.value = point.value - k*dimension
                #print(p.value)
            #exit()

    vertices = []
    triangles = []

    for k in range(len(grid)-1):
        print("K:", k)
        for j in range(len(grid[0])-1):
            for i in range(len(grid[0][0])-1):
                cell = grid_ops.cell_for_indices(grid, (k, j, i))
                polygonise_cube(cell, 0, vertices, triangles)

    points = [str(i) for i in vertices]
    indices = [str(i) for i in triangles]

    meshfile = open("artificial_plane_100.off","w")
    meshfile.write(
    '''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    meshfile.close()

def cube_test():
    half_size = 3
    box_volume = BoundingBox(-half_size, -half_size, -half_size, half_size, half_size, half_size)

    # vertexBuffer = [
    #     Vertex(-1.0, -1.0, -2.0),
    #     Vertex(1.0, -1.0, -2.0),
    #     Vertex(1.0, -1.0, 0.0),
    #     Vertex(-1.0, -1.0, 0.0),
    #
    #     Vertex(-1.0, 1.0, -2.0),
    #     Vertex(1.0, 1.0, -2.0),
    #     Vertex(1.0, 1.0, 0.0),
    #     Vertex(-1.0, 1.0, 0.0)
    # ]
    #
    # indexBuffer = [
    #     Triangle(0, 3, 1),
    #     Triangle(3, 2, 1),
    #     Triangle(4, 7, 5),
    #     Triangle(7, 6, 5),
    #     Triangle(6, 2, 5),
    #     Triangle(2, 1, 5),
    #     Triangle(4, 0, 7),
    #     Triangle(0, 3, 7),
    #     Triangle(7, 3, 6),
    #     Triangle(3, 2, 6),
    #     Triangle(4, 0, 5),
    #     Triangle(0, 1, 5)
    # ]


    vertexBuffer = [
        Vertex(2.0, 2.0, -2.0),
        Vertex(-2.0, 2.0, -2.0),
        Vertex(-2.0, -2.0, -2.0),
        Vertex(2.0, -2.0, -2.0),
    ]

    indexBuffer = [
        Triangle(0, 1, 2),
        Triangle(2, 3, 0)
    ]

    grid, dimension = grid_ops.make_grid(box_volume, 100)

    #grid_ops.assing_values(grid, vertexBuffer, indexBuffer, box_volume.min_corner, dimension)
    grid_ops.assing_all_values(grid, vertexBuffer, indexBuffer, box_volume.min_corner, dimension)
    #exit()
    vertices = []
    triangles = []

    # for k in range(len(grid)-1):
    #     print(grid[k][4][3].value)

    # for k in range(len(grid)):
    #     print("K: {}, 16: {}, 17: {}; 83: {}, 84: {}".format(k, grid[k][20][16].value,
    #                                                         grid[k][20][17].value,
    #                                                         grid[k][20][83].value,
    #                                                         grid[k][20][84].value))
    # print("HOP!!!!!!!!")

    flag = False
    for k in range(len(grid)-1):
        print("K:", k)
        for j in range(len(grid[0])-1):
            for i in range(len(grid[0][0])-1):
                cell = grid_ops.cell_for_indices(grid, (k, j, i))
                for vertex in cell.vertices:
                    if vertex.value == grid_ops.SENTINEL_VALUE:
                        flag = True
                        break
                if(flag):
                    flag = False
                    continue
                polygonise_cube(cell, 0, vertices, triangles)

    points = [str(i) for i in vertices]
    indices = [str(i) for i in triangles]

    meshfile = open("plane100_2triangles.off","w")
    meshfile.write(
    '''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    meshfile.close()

def using_off():
    vertexBuffer, indexBuffer, box_volume = grid_ops.read_OFF("models/simplified_vase.off")

    grid, dimension = grid_ops.make_grid(box_volume, 400)

    grid_ops.assing_values(grid, vertexBuffer, indexBuffer, box_volume.min_corner, dimension)
    #grid_ops.assing_all_values(grid, vertexBuffer, indexBuffer, box_volume.min_corner, dimension)

    vertices = []
    triangles = []

    flag = False
    for k in range(len(grid)-1):
        print("K:", k)
        for j in range(len(grid[0])-1):
            for i in range(len(grid[0][0])-1):
                cell = grid_ops.cell_for_indices(grid, (k, j, i))
                for vertex in cell.vertices:
                    if vertex.value == grid_ops.SENTINEL_VALUE:
                        flag = True
                        break
                if(flag):
                    flag = False
                    continue
                polygonise_cube(cell, 0, vertices, triangles)

    points = [str(i) for i in vertices]
    indices = [str(i) for i in triangles]

    meshfile = open("simplified400_moller.off","w")
    meshfile.write(
    '''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    meshfile.close()

if __name__ == '__main__':
    using_off()
    #main()
    #poly_function()
    #cube_test()
    #plane_test()

    #artificial_valoration()
