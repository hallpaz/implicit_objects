import argparse
import sys
import os
from PIL import Image
from DataStructures import Vertex, Triangle, DepthMap, BoundingBox

focalLength = 525.0
centerX = 319.5
centerY = 239.5
scalingFactor = 1000.00

BIG = 9999999.99

def generateOFF(rgb_file,depth_file,off_file):
    """
    Generate a colored point cloud in PLY format from a color and a depth image.

    Input:
    rgb_file -- filename of color image
    depth_file -- filename of depth image
    off_file -- filename of off file

    """
    rgb = Image.open(rgb_file)
    depth = Image.open(depth_file)

    print(rgb.mode, "rgb")

    if rgb.size != depth.size:
        raise Exception("Color and depth image do not have the same resolution.")
    if rgb.mode != "RGB":
        raise Exception("Color image is not in RGB format")
    if depth.mode != "I":
        print(depth.mode, "depth")
        raise Exception("Depth image is not in intensity format")

    mean_distance = 0.0
    min_distance = BIG
    minx = BIG
    miny = BIG
    minz = BIG
    maxx = -BIG
    maxy = -BIG
    maxz = -BIG
    max_distance = 0.0
    points = []
    vertices = []
    last_depth = 0.0
    divide_factor = 0
    for v in range(rgb.size[1]):
        for u in range(rgb.size[0]):
            #color = rgb.getpixel((u,v))
            Z = -depth.getpixel((u,v)) / scalingFactor
            if Z==0:
                continue
            X = -(u - centerX) * Z / focalLength
            Y = (v - centerY) * Z / focalLength

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

            distance = abs(last_depth - Z)
            mean_distance += distance
            if distance < min_distance:
                min_distance = distance
            if distance > max_distance:
                max_distance = distance
            last_depth = Z
            divide_factor += 1
            vertex = Vertex(X,Y,Z)
            points.append(str(vertex))
            vertices.append(vertex)

    bounding_box = BoundingBox(minx, miny, minz, maxx, maxy, maxz)
    #print(bounding_box)

    mean_distance = mean_distance / divide_factor
    #print(mean_distance, min_distance, max_distance)
    width = depth.size[0]
    height = depth.size[1]
    indices = []
    faces = []
    distanceThreshold = mean_distance*16
    for v in range(height - 1):
        for u in range(width - 1):
            #superior triangle
            hdistance = abs(depth.getpixel((u, v)) - depth.getpixel((u+1, v))) / scalingFactor
            vdistance = abs(depth.getpixel((u, v)) - depth.getpixel((u, v+1))) / scalingFactor
            if hdistance < distanceThreshold and vdistance < distanceThreshold:
                triangle = Triangle(v*width + u + 1, v*width + u, (v+1)*width + u )
                faces.append(triangle)
                indices.append( str(triangle))
            #inferior triangle
            hdistance = abs(depth.getpixel((u, v+1)) - depth.getpixel((u+1, v+1))) / scalingFactor
            vdistance = abs(depth.getpixel((u+1, v)) - depth.getpixel((u+1, v+1))) /scalingFactor
            if hdistance < distanceThreshold and vdistance < distanceThreshold:
                triangle = Triangle( v*width + u + 1 , (v+1)*width + u, (v+1)*width + u + 1 )
                faces.append(triangle)
                indices.append(str(triangle))
            #faces.append("4 %d %d %d %d\n"%( v*width + u, v*width + u + 1 , (v+1)*width + u + 1 , (v+1)*width + u ))

    meshfile = open(off_file,"w")
    meshfile.write('''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    meshfile.close()

    points = bounding_box.points()
    indices = bounding_box.indices()
    boxfile = open("models/box.off", "w")
    boxfile.write('''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    boxfile.close()
    return vertices, faces, bounding_box


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a registered pair of color and depth images and generates a colored 3D point cloud in the
    PLY format.
    ''')
    parser.add_argument('rgb_file', help='input color image (format: png)')
    parser.add_argument('depth_file', help='input depth image (format: png)')
    parser.add_argument('ply_file', help='output PLY file (format: ply)')
    args = parser.parse_args()

    generateOFF(args.rgb_file,args.depth_file,args.ply_file)
