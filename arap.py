import sys
import numpy as np
import math
import face
import offfile

# Read file into arrays
class Deformer:
    def __init__(self, filename):
        self.filename = filename

    def readFile(self):
        fr = offfile.OffFile(self.filename)

        first_line = fr.nextLine().split()

        number_of_verticies =   int(first_line[0])
        number_of_faces =       int(first_line[1])
        number_of_edges =       int(first_line[2])

        # Every vertex in the .off
        self.verts = []
        # Every face in the .off
        self.faces = []
        # The ID of the faces related to this vertx ID (i.e. vtf[i] contains faces that contain ID i)
        self.vertsToFaces = []

        for i in range(0, number_of_verticies):
            vert_line = fr.nextLine().split()
            x = float(vert_line[0])
            y = float(vert_line[1])
            z = float(vert_line[2])
            self.verts.append(np.array([x, y, z]))
            self.vertsToFaces.append([])

        for i in range(0, number_of_faces):
            faceLine = fr.nextLine().split()
            v1_id = int(faceLine[1])
            v2_id = int(faceLine[2])
            v3_id = int(faceLine[3])
            self.faces.append(face.Face(v1_id, v2_id, v3_id))
            # Add this face to each vertex face map
            self.vertsToFaces[v1_id].append(i)
            self.vertsToFaces[v2_id].append(i)
            self.vertsToFaces[v3_id].append(i)

        print(str(len(self.verts)) + " verticies")
        print(str(len(self.faces)) + " faces")
        print(str(number_of_edges) + " edges")

    # Returns a set of IDs that are neighbours to this vertexID (not including the input ID)
    def neighboursOf(self, vertID):
        neighbours = []
        for face in self.vertsToFaces[vertID]:
            for vID in face.vertexIDs:
                neighbours.append(vID)
        neighbours = set(neighbours)
        return neighbours.remove(vertID)

    def buildWeightMatrix(self):
        number_of_verticies = len(self.verts)

        self.weightMatrix = [
            [ None for i in range(number_of_verticies) ]
        ] * number_of_verticies

        for vertex_id in number_of_verticies:
            for neighbour_id in neighboursOf(vertex_id):
                assignWeightForPair(vertex_id, neighbour_id)

    def assignWeightForPair(self, i, j):

        self.weightMatrix[i][j]

    def weightForPair(self, i, j):
        faces = []
        for f_id in self.vertsToFaces[i]:
            face = self.faces[f_id]
            if face.containsPointIDs(i, j):
                faces.apped(face)
        #If there's more than two, we have a problem, or mesh is bad
        assert(len(faces) == 2)

        

# MAIN
filename = "data/02-bar-twist/00-bar-original.off"

if(len(sys.argv) > 1):
    filename = sys.argv[1]
d = Deformer(filename)
d.readFile()
d.buildWeightMatrix()
