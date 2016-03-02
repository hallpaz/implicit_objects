import operator
from LookUpTables import edgeTable, triTable
from DataStructures import BoundingBox, GridCell, Vertex, Triangle



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

    return volume_of_cubes


def assign_values(grid, function):
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            for k in range(len(grid[0][0])):
                cell = grid[i][j][k]
                for u in range(8):
                    value = function(cell.positions[u].x, cell.positions[u].y, cell.positions[u].z)
                    cell.values[u] = value

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


def triangle_intersection(v1, v2, triangle_vertices):
    #TODO: code a proper triangle segment intercetion
    has_intersection = False
    intersection_point = Vertex(0, 0, 0)
    if has_intersection:
        return intersection_point

    return None

def signed_distance(cell, triangle):
    

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


def sphere(x, y, z):
    return x**2 + y**2 + z**2

def bitorus(x, y, z):
    return ((x**2 + y**2)**2 - x**2 + y**2)**2 + z**2

if __name__ == '__main__':
    half_size = 2
    box_volume = BoundingBox(-half_size, -half_size, -half_size, half_size, half_size, half_size)

    grid = makeGrid(box_volume, 70)
    #assign_values(grid, sphere)
    assign_values(grid, bitorus)
    vertices = []
    triangles = []
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            for k in range(len(grid[0][0])):
                cell = grid[i][j][k]
                polygonise_cube(cell, 0.01, vertices, triangles)

    points = [str(i) for i in vertices]
    indices = [str(i) for i in triangles]

    meshfile = open("torus.off","w")
    meshfile.write(
    '''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    meshfile.close()
