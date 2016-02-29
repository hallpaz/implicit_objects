class DepthMap:
    def __init__(self, points, width, height):
        self.points = points
        self.width = width
        self.height = height


class Vertex:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "%f %f %f\n"%(self.x,self.y,self.z)


class Triangle:
    #triangle indices
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

    def __str__(self):
        return "3 %d %d %d\n"%(self.a, self.b, self.c)


class Cube:
    def __init__(self, dimension):


class BoundingBox:
    def __init__(self, minx, miny, minz, maxx, maxy, maxz):
        self.min_corner = (minx, miny, minz)
        self.max_corner = (maxx, maxy, maxz)
        self.vertices = [
            Vertex(minx, miny, minz),
            Vertex(maxx, miny, minz),
            Vertex(maxx, maxy, minz),
            Vertex(minx, maxy, minz),
            Vertex(minx, miny, maxz),
            Vertex(maxx, miny, maxz),
            Vertex(maxx, maxy, maxz),
            Vertex(minx, maxy, maxz)
        ]
        self.faces = [
            Triangle(2, 3, 0),
            Triangle(0, 1, 2),
            Triangle(6, 7, 4),
            Triangle(4, 5, 6)

            #Triangle(1, 5, 6),
            #Triangle(6, 2, 1),
            #Triangle(5, 4, 7),
            #Triangle(7, 6, 5),
            #Triangle(4, 1, 3),
            #Triangle(3, 7, 4)
        ]

    def points(self):
        return [str(i) for i in self.vertices]

    def indices(self):
        return [str(i) for i in self.faces]


    def __str__(self):
        return "Min Corner: " + str(self.min_corner) + " Max Corner: " + str(self.max_corner)




def vertex_interpolation(isovalue, p1, p2):
    

XYZ VertexInterp(isolevel,p1,p2,valp1,valp2)
double isolevel;
XYZ p1,p2;
double valp1,valp2;
{
   double mu;
   XYZ p;

   if (ABS(isolevel-valp1) < 0.00001)
      return(p1);
   if (ABS(isolevel-valp2) < 0.00001)
      return(p2);
   if (ABS(valp1-valp2) < 0.00001)
      return(p1);
   mu = (isolevel - valp1) / (valp2 - valp1);
   p.x = p1.x + mu * (p2.x - p1.x);
   p.y = p1.y + mu * (p2.y - p1.y);
   p.z = p1.z + mu * (p2.z - p1.z);

   return(p);
}
