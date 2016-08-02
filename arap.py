import sys
import numpy as np
import math
import face
import offfile
import othermath as omath
np.set_printoptions(precision=2)
# Read file into arrays
class Deformer:
    def __init__(self, filename):
        self.filename = filename

    def read_file(self):
        fr = offfile.OffFile(self.filename)

        first_line = fr.nextLine().split()

        number_of_verticies =   int(first_line[0])
        number_of_faces =       int(first_line[1])
        number_of_edges =       int(first_line[2])

        # Every vertex in the .off
        self.verts = []
        self.verts_prime = []
        # Every face in the .off
        self.faces = []
        # The ID of the faces related to this vertx ID (i.e. vtf[i] contains faces that contain ID i)
        self.verts_to_face = []

        for i in range(0, number_of_verticies):
            vert_line = fr.nextLine().split()
            x = float(vert_line[0])
            y = float(vert_line[1])
            z = float(vert_line[2])
            self.verts.append(np.array([x, y, z]))
            self.verts_prime.append(np.array([x, y, z]))

            self.verts_to_face.append([])

        self.neighbour_matrix = [ [ 0 for x in range(0, number_of_verticies) ] for y in range(0, number_of_verticies) ]
        self.neighbour_matrix = np.matrix(self.neighbour_matrix)

        print("Generating Adjacencies")
        for i in range(0, number_of_faces):
            face_line = fr.nextLine().split()
            v1_id = int(face_line[1])
            v2_id = int(face_line[2])
            v3_id = int(face_line[3])
            self.faces.append(face.Face(v1_id, v2_id, v3_id))
            # Add this face to each vertex face map
            self.assign_values_to_neighbour_matrix(v1_id, v2_id, v3_id)
            self.verts_to_face[v1_id].append(i)
            self.verts_to_face[v2_id].append(i)
            self.verts_to_face[v3_id].append(i)

        print("Generating Edge Matrix")
        self.edge_matrix = np.zeros((number_of_verticies, number_of_verticies))
        # TODO APPLY THIS TO OTHER STUFF ^^
        for row in range(number_of_verticies):
            self.edge_matrix[row][row] = self.neighbour_matrix[row].sum()
        print("Generating Laplacian Matrix")
        self.laplacian_matrix = self.edge_matrix - self.neighbour_matrix

        print(str(len(self.verts)) + " verticies")
        print(str(len(self.faces)) + " faces")
        print(str(number_of_edges) + " edges")

    def assign_values_to_neighbour_matrix(self, v1, v2 ,v3):
        self.neighbour_matrix[v1, v2] = 1
        self.neighbour_matrix[v2, v1] = 1
        self.neighbour_matrix[v1, v3] = 1
        self.neighbour_matrix[v3, v1] = 1
        self.neighbour_matrix[v2, v3] = 1
        self.neighbour_matrix[v3, v2] = 1

    # Reads the .sel file and keeps track of the selection status of a vertex
    def read_selection_file(self, filename):
        # The selection status of each vertex, where 0=fixed, 1=deformable-region, 2=handle

        self.vert_status = open(filename, 'r').read().strip().split("\n")
        # Remove any lines that aren't numbers
        self.vert_status = [int(line) for line in self.vert_status if omath.string_is_int(line)]

        # Keep track of the IDs of the selected verts (i.e. verts with handles/status == 2)
        self.selected_verts = []
        for i in range(0, len(self.vert_status)):
            if self.vert_status[i] == 2:
                self.selected_verts.append(i)
        assert(len(self.vert_status) == len(self.verts))

    # Reads the .def file and stores the inner matrix
    def read_deformation_file(self, filename):
        def_file_lines = open(filename, 'r').read().strip().split("\n")
        # Remove any lines with comments
        def_file_lines = [line for line in def_file_lines if "#" not in line]

        # Assert its at least a 4 by something matrix just in case
        assert(len(def_file_lines) == 4)
        self.deformation_matrix = np.matrix(";".join(def_file_lines))
        print("Deformation matrix to apply")
        print(self.deformation_matrix)
        assert(self.deformation_matrix.size == 16)

    # Returns a set of IDs that are neighbours to this vertexID (not including the input ID)
    def neighbours_of(self, vert_id):
        neighbours = []
        for n_id in range(0, len(self.verts)):
            if(self.neighbour_matrix[vert_id, n_id] == 1):
                neighbours.append(n_id)
        return neighbours

    def build_weight_matrix(self):
        number_of_verticies = len(self.verts)

        self.weight_matrix = np.zeros((number_of_verticies, number_of_verticies))
        for vertex_id in range(number_of_verticies):
            for neighbour_id in self.neighbours_of(vertex_id):
                self.assign_weight_for_pair(vertex_id, neighbour_id)
        print("Matix complete. ", len(self.weight_matrix)**2, " entries")
        print(self.weight_matrix)

    # def build_laplacian_matrix(self):
    #     

    def assign_weight_for_pair(self, i, j):
        if(self.weight_matrix[j, i] == 0):
            # If the opposite weight has not been computed, do so
            weightIJ = self.weight_for_pair(i, j)
        else:
            weightIJ = self.weight_matrix[j, i]
        self.weight_matrix[i, j] = weightIJ

    def weight_for_pair(self, i, j):
        local_faces = []
        # For every face associated with vert index I,
        for f_id in self.verts_to_face[i]:
            face = self.faces[f_id]
            # If the face contains both I and J, add it
            if face.contains_point_ids(i, j):
                local_faces.append(face)

        # Either a normal face or a boundry edge, otherwise bad mesh
        assert(len(local_faces) <= 2)

        vertex_i = self.verts[i]
        vertex_j = self.verts[j]

        # weight equation: 0.5 * (cot(alpha) + cot(beta)

        cot_theta_sum = 0
        for face in local_faces:
            other_vertex_id = face.other_point(i, j)
            vertex_o = self.verts[other_vertex_id]
            theta = omath.angle_between(vertex_j - vertex_o, vertex_i - vertex_o)
            cot_theta_sum += omath.cot(theta)
        return cot_theta_sum * 0.5

    def apply_deformation(self):
        print("Length of sel verts", len(self.selected_verts))
        # Apply first deformation
        for vert_id in self.selected_verts:
            vert = self.verts[vert_id]
            new_vert = omath.apply_rotation(self.deformation_matrix, vert)
            self.verts_prime[vert_id] = new_vert
        self.cell_rotations = []

        for i in range(3):
            # Calculate rotations
            # Calculate and apply rotations for each cell TODO MOVE SOMEWHERE ELSE
            for vert_id in range(0, len(self.verts)):
                if(self.vert_status[vert_id] == 1):
                    rotation = self.calculate_rotation_matrix_for_cell(vert_id)
                else:
                    rotation = np.identity(3)
                vert = self.verts[vert_id]

                self.cell_rotations.append(rotation)

            b_array = []
            for vert_id in range(len(self.verts)):
                b_array.append(self.calculate_b_for(vert_id))
            L = self.laplacian_matrix
            n = len(self.verts)
            for vert_id in range(n):
                if(self.vert_status[vert_id] != 1):
                    L[vert_id] = 0
                    L[:, vert_id] = 0
                    L[vert_id, vert_id] = 1
                    b_array[vert_id] = self.verts_prime[vert_id]

            p_prime = np.linalg.solve(self.laplacian_matrix, np.array(b_array))
            self.verts_prime = p_prime
            print(p_prime)


    def calculate_rotation_matrix_for_cell(self, vert_id):
        covariance_matrix = self.calculate_covariance_matrix_for_cell(vert_id)

        U, s, V_transpose = np.linalg.svd(covariance_matrix, full_matrices=True, compute_uv=True)
        # U, S, V_transpose
        # V_transpose_transpose * U_transpose
        return V_transpose.transpose() * U.transpose()

    def calculate_covariance_matrix_for_cell(self, vert_id):
        #s_i = P_i * D_i * P_i_prime_transpose
        vert_i = self.verts[vert_id]
        vert_i_prime = self.verts_prime[vert_id]

        neighbour_ids = self.neighbours_of(vert_id)
        D_i = np.diag(np.array([ self.weight_matrix[vert_id, n_id] for n_id in neighbour_ids ]))
        P_i = []
        P_i_prime = []
        for j in neighbour_ids:
            vert_j = self.verts[j]
            P_i.append(vert_i - vert_j)

            vert_j_prime = self.verts_prime[j]

            P_i_prime.append(vert_i_prime - vert_j_prime)

        P_i = np.matrix(P_i).transpose()

        P_i_prime = np.matrix(P_i_prime)
        return P_i * D_i * P_i_prime

    def output_s_prime_to_file(self):
        print("Writing to `output.off`")
        f = open('output.off', 'w')
        f.write("OFF\n")
        f.write(str(len(self.verts)) + " " + str(len(self.faces)) + " 0\n")
        for vert in self.verts_prime:
            for i in vert:
                f.write(str(i) + " ")
            f.write("\n")
        for face in self.faces:
            f.write(face.off_string() + "\n")
        f.close()
        print("Output file to `output.off`")

    def calculate_b_for(self, vert_id):
        b = np.zeros(3)
        for n_id in self.neighbours_of(vert_id):
            w_ij = self.weight_matrix[vert_id, n_id] / 2
            r_ij = self.cell_rotations[vert_id] + self.cell_rotations[n_id]
            p_ij = self.verts_prime[vert_id] - self.verts_prime[n_id]
            b += (w_ij * omath.apply_rotation(r_ij, p_ij))
        return b

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
d.read_file()
if len(selection_filename) > 0:
    d.read_selection_file(selection_filename)
    d.read_deformation_file(deformation_file)
d.build_weight_matrix()
# d.build_laplacian_matrix()
d.apply_deformation()
d.output_s_prime_to_file()