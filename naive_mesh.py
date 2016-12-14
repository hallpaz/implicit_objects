import argparse
import sys
import os
from PIL import Image
from DataStructures import Vertex, Triangle, DepthMap, BoundingBox
import numpy as np

focalLength = 525.0
centerX = 319.5
centerY = 239.5
scalingFactor = 1000.00

BIG = 9999999.99
# fx = 525.0; fy = 525.0; // default focal length
# cx = 319.5; cy = 239.5; // default optical center
# Eigen::Matrix4d Tk; // k-th transformation matrix from trajectory.log
#
# // translation from depth pixel (u,v,d) to a point (x,y,z)
# z = d / 1000.0;
# x = (u - cx) * z / fx;
# y = (v - cy) * z / fy;
#
# // transform (x,y,z) to (xw,yw,zw)
# Eigen::Vector4d w = Tk * Eigen::Vector4d(x, y, z, 1);
# xw = w(0); yw = w(1); zw = w(2);
#matrix 2343
# matrix = [[-0.5865089568358111, -0.4016671627616591, 0.7033283258268271, -1.78958890904936],
# [3.904256506020038e-17, -0.8683679602305977, -0.49592044285847997, 1.2941331133042304],
# [0.8099427409091144, -0.29086178161448045, 0.5093055865044889, -0.3940599767154947],
# [0.0, 0.0, 0.0, 1.0]]

#matrix 1991
# matrix = [[0.5405581528083065, 0.28655525127420545, -0.7910012461429174, 1.405061042371321],
# [0.0, -0.9402056240170099, -0.3406073759723719, 1.2016302960981653],
# [-0.8413066524356453, 0.18411809398850973, -0.5082358153786161, 1.7523119096340265],
# [0.0, 0.0, 0.0, 1.0]]

#matrix2039
# matrix = [[0.1251551653807354, 0.32596834683277814, -0.9370596680263356, 1.6348786428914903],
# [-6.5021576180652725e-18, -0.9444859910618658, -0.3285516895221894, 1.1836401989655352],
# [-0.992137180322621, 0.041119941038269645, -0.11820730041113557, 0.8277766725924591],
# [0.0, 0.0, 0.0, 1.0]]

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

    width = depth.size[0]
    height = depth.size[1]

    filenumber = int(depth_file[-9:-4])
    matrix = []
    with open("stonewall_trajectory.log") as trajfile:
        number = int(trajfile.readline().split()[2])
        while number != "":
            if number == filenumber:
                matrix_data = []
                for i in range(4):
                    matrix_data.append([float(n) for n in trajfile.readline().split()])
                matrix = np.matrix(matrix_data)
                break
            for i in range(4):
                trajfile.readline()
            number = int(trajfile.readline().split()[2])

    for v in range(rgb.size[1]):
        for u in range(rgb.size[0]):
            #color = rgb.getpixel((u,v))
            #Z = -depth.getpixel((u,v)) / scalingFactor
            Z = depth.getpixel((u,v)) / scalingFactor
            #print(Z)
            if Z==0:
                continue
            #X = -(u - centerX) * Z / focalLength
            X = (u - centerX) * Z / focalLength
            Y = (v - centerY) * Z / focalLength

            vertex = np.array([X, Y, Z, 1])
            vertex = np.array(matrix*np.transpose(np.matrix(vertex)))
            vertex.shape = 4
            #print(vertex)
            X, Y, Z = vertex[0], vertex[1], vertex[2]

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

    bounding_box = BoundingBox(minx, miny, minz-1, maxx, maxy, maxz+1)
    # bounding_box = BoundingBox(minx, miny, minz, maxx, maxy, maxz)
    #print(bounding_box)

    mean_distance = mean_distance / divide_factor
    #print(mean_distance, min_distance, max_distance)

    indices = []
    faces = []
    # distanceThreshold = mean_distance*20
    distanceThreshold = min_distance

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

    meshfile = open(off_file,"w")
    meshfile.write('''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    meshfile.close()

    points = bounding_box.points()
    indices = bounding_box.indices()
    boxfile = open("models/vasebox.off", "w")
    boxfile.write('''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    boxfile.close()
    return vertices, faces, bounding_box


def generatePlaneOFF(off_file):

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
    scale_factor = 1

    width = 640 #depth.size[0]
    height = 480 #depth.size[1]

    for v in range(height):
        for u in range(width):
            #if u > width/6 and u < 5*width/6 and v > height/6 and v < 5*height/6:
            X = scale_factor*3*u/width
            Y = scale_factor*3*v/height
            Z = -3
            vertex = Vertex(X, Y, Z)
            points.append(str(vertex))
            vertices.append(vertex)
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


    bounding_box = BoundingBox(minx, miny, minz-1, maxx, maxy, maxz+1)

    indices = []
    faces = []
    distanceThreshold = mean_distance*8

    for v in range(height - 1):
        for u in range(width - 1):
            #superior triangle
            triangle = Triangle(v*width + u + 1, v*width + u, (v+1)*width + u )
            faces.append(triangle)
            indices.append( str(triangle))
            #inferior triangle
            triangle = Triangle( v*width + u + 1 , (v+1)*width + u, (v+1)*width + u + 1 )
            faces.append(triangle)
            indices.append(str(triangle))

    meshfile = open(off_file,"w")
    meshfile.write('''OFF
    %d %d 0
    %s%s
    '''%(len(points),len(indices), "".join(points), "".join(indices)))
    meshfile.close()

    points = bounding_box.points()
    indices = bounding_box.indices()
    boxfile = open("models/box_plane.off", "w")
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
