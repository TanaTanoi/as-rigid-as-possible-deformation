import sys
import numpy as np
import math
import face
import offfile
import othermath as omath

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
        self.vertsPrime = []
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
            self.vertsPrime.append(np.array([x, y, z]))

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

    # Reads the .sel file and keeps track of the selection status of a vertex
    def readSelectionFile(self, filename):
        # The selection status of each vertex, where 0=fixed, 1=deformable-region, 2=handle

        self.vertStatus = open(filename, 'r').read().strip().split("\n")
        # Remove any lines that aren't numbers
        self.vertStatus = [int(line) for line in self.vertStatus if omath.string_is_int(line)]

        # Keep track of the IDs of the selected verts (i.e. verts with handles/status == 2)
        self.selectedVerts = []
        for i in range(0, len(self.vertStatus)):
            if self.vertStatus[i] == 2:
                self.selectedVerts.append(i)
        assert(len(self.vertStatus) == len(self.verts))

    # Reads the .def file and stores the inner matrix
    def readDeformationFile(self, filename):
        defFileLines = open(filename, 'r').read().strip().split("\n")
        # Remove any lines with comments
        defFileLines = [line for line in defFileLines if "#" not in line]

        # Assert its at least a 4 by something matrix just in case
        assert(len(defFileLines) == 4)
        self.deformationMatrix = np.matrix(";".join(defFileLines))
        print("Deformation matrix to apply")
        print(self.deformationMatrix)
        assert(self.deformationMatrix.size == 16)

    # Returns a set of IDs that are neighbours to this vertexID (not including the input ID)
    def neighboursOf(self, vertID):
        neighbours = []
        for faceID in self.vertsToFaces[vertID]:
            face = self.faces[faceID]
            for vID in face.vertexIDs():
                neighbours.append(vID)
        neighbours = set(neighbours)
        neighbours.remove(vertID)
        return neighbours

    def buildWeightMatrix(self):
        number_of_verticies = len(self.verts)

        self.weightMatrix = [
            [ 0 for i in range(number_of_verticies) ] for j in range(number_of_verticies)
        ]
        for vertex_id in range(number_of_verticies):
            for neighbour_id in self.neighboursOf(vertex_id):
                self.assignWeightForPair(vertex_id, neighbour_id)
        print("Matix complete. ", len(self.weightMatrix)**2, " entries")

    def assignWeightForPair(self, i, j):
        if(self.weightMatrix[j][i] == 0):
            # If the opposite weight has not been computed, do so
            weightIJ = self.weightForPair(i, j)
        else:
            weightIJ = self.weightMatrix[j][i]
        self.weightMatrix[i][j] = weightIJ

    def weightForPair(self, i, j):
        local_faces = []
        # For every face associated with vert index I,
        for f_id in self.vertsToFaces[i]:
            face = self.faces[f_id]
            # If the face contains both I and J, add it
            if face.containsPointIDs(i, j):
                local_faces.append(face)

        # Either a normal face or a boundry edge, otherwise bad mesh
        assert(len(local_faces) <= 2)

        vertex_i = self.verts[i]
        vertex_j = self.verts[j]

        # weight equation: 0.5 * (cot(alpha) + cot(beta)

        cot_theta_sum = 0
        for face in local_faces:
            otherVertexID = face.otherPoint(i, j)
            vertex_o = self.verts[otherVertexID]
            theta = omath.angleBetween(vertex_j - vertex_o, vertex_i - vertex_o)
            cot_theta_sum += omath.cot(theta)
        return cot_theta_sum * 0.5

    def apply_deformation(self):
        print("Length of sel verts", len(self.selectedVerts))
        # Apply first deformation
        for vert_id in self.selectedVerts:
            vert = self.verts[vert_id]
            self.vertsPrime[vert_id] = omath.apply_rotation(self.deformationMatrix, vert)

        self.cellRotations = []
        # Calculate and apply rotations for each cell TODO MOVE SOMEWHERE ELSE
        for vert_id in range(0, len(self.verts)):
            if(self.vertStatus[vert_id] != 0):
                rotation = self.calculate_rotation_matrix_for_cell(vert_id)
            else:
                rotation = np.identity(3)
            vert = self.verts[vert_id]

            self.cellRotations.append(rotation)
        for vert_id in range(0, len(self.verts)):
            self.vertsPrime[vert_id] = omath.apply_rotation(self.cellRotations[vert_id], self.verts[vert_id])

    def calculate_rotation_matrix_for_cell(self, vert_id):
        covariance_matrix = self.calculate_covariance_matrix_for_cell(vert_id)

        U, s, V = np.linalg.svd(covariance_matrix, full_matrices=True)
        # U, S, V_transpose
        # V_transpose_transpose * U_transpose
        return V.transpose() * U.transpose()

    def calculate_covariance_matrix_for_cell(self, vert_id):
        #s_i = P_i * D_i * P_i_prime_transpose
        vert_i = self.verts[vert_id]
        vert_i_prime = self.vertsPrime[vert_id]

        neighbour_ids = self.neighboursOf(vert_id)
        D_i = np.diag(np.array([ self.weightMatrix[vert_id][n_id] for n_id in neighbour_ids ]))
        P_i = []
        P_i_prime = []
        for j in neighbour_ids:
            vert_j = self.verts[j]
            P_i.append(vert_i - vert_j)

            vert_j_prime = self.vertsPrime[j]
            P_i_prime.append(vert_i_prime - vert_j_prime)

        P_i = np.matrix(P_i).transpose()
        P_i_prime = np.matrix(P_i_prime)

        return P_i * D_i * P_i_prime

    def output_s_prime_to_file(self):
        f = open('output.off', 'w')
        f.write("OFF\n")
        f.write(str(len(self.verts)) + " " + str(len(self.faces)) + " 0\n")
        for vert in self.vertsPrime:
            for i in vert:
                f.write(str(i) + " ")
            f.write("\n")
        for face in self.faces:
            f.write(face.off_string() + "\n")
        f.close()
        print("Output file to `output.off`")

# MAIN
filename = "data/02-bar-twist/00-bar-original.off"
selection_filename = "data/02-bar-twist/bar.sel"
deformation_file = "data/02-bar-twist/bar.def"

argc = len(sys.argv)

if(argc > 1):
    filename = sys.argv[1]
    selection_filename = ""
    deformation_file = ""
if(argc > 2):
    selection_filename = sys.argv[2]
    deformation_file = ""
if(argc > 3):
    deformation_file = sys.argv[3]

d = Deformer(filename)
d.readFile()
if len(selection_filename) > 0:
    d.readSelectionFile(selection_filename)
    d.readDeformationFile(deformation_file)
d.buildWeightMatrix()
d.apply_deformation()
d.output_s_prime_to_file()