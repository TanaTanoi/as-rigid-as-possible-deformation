import time
import sys
import os
import numpy as np
import math
import face
import offfile
import othermath as omath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
np.set_printoptions(precision=2, suppress=True)

SPARSE = False

if SPARSE:
    import scipy.sparse
    import scipy.sparse.linalg

    solve   = scipy.sparse.linalg.spsolve
    matrix  = scipy.sparse.lil_matrix
else:
    solve   = np.linalg.solve
    matrix  = np.zeros

# Read file into arrays
class Deformer:
    max_iterations = 100
    threshold = 0.001
    def __init__(self, filename):
        self.filename = filename
        self.POWER = 8

    def read_file(self):
        fr = offfile.OffFile(self.filename)

        first_line = fr.nextLine().split()

        number_of_verticies =   int(first_line[0])
        number_of_faces =       int(first_line[1])
        number_of_edges =       int(first_line[2])

        self.n = number_of_verticies

        # Every vertex in the .off
        self.verts = []
        self.verts_prime = []
        # Every face in the .off
        self.faces = []
        # The ID of the faces related to this vertx ID (i.e. vtf[i] contains faces that contain ID i)
        self.verts_to_face = []

        for i in range(self.n):
            vert_line = fr.nextLine().split()
            x = float(vert_line[0])
            y = float(vert_line[1])
            z = float(vert_line[2])
            self.verts.append(np.array([x, y, z]))
            self.verts_prime.append(np.array([x, y, z]))

            self.verts_to_face.append([])
        self.verts_prime = np.asmatrix(self.verts_prime)
        self.neighbour_matrix = np.zeros((self.n, self.n))

        print("Generating Adjacencies")
        for i in range(number_of_faces):
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
        self.edge_matrix = np.zeros((self.n, self.n))

        for row in range(self.n):
            self.edge_matrix[row][row] = self.neighbour_matrix[row].sum()
        print("Generating Laplacian Matrix")

        # N size array of 3x3 matricies
        self.cell_rotations = np.zeros((self.n, 3, 3))

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
        self.fixed_verts = []
        for i in range(self.n):
            if self.vert_status[i] == 2:
                self.selected_verts.append(i)
                self.fixed_verts.append((i, omath.apply_rotation(self.deformation_matrix, self.verts[i])))
            elif self.vert_status[i] == 0:
                self.fixed_verts.append((i, self.verts[i]))
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
        for n_id in range(self.n):
            if(self.neighbour_matrix[vert_id, n_id] == 1):
                neighbours.append(n_id)
        return neighbours

    def build_weight_matrix(self):
        print("Generating Weight Matrix")
        self.weight_matrix = matrix((self.n, self.n), dtype=np.float)
        self.weight_sum = matrix((self.n, self.n), dtype=np.float)

        for vertex_id in range(self.n):
            neighbours = self.neighbours_of(vertex_id)
            for neighbour_id in neighbours:
                self.assign_weight_for_pair(vertex_id, neighbour_id)
        print(self.weight_matrix)

    def assign_weight_for_pair(self, i, j):
        if(self.weight_matrix[j, i] == 0):
            # If the opposite weight has not been computed, do so
            weightIJ = self.weight_for_pair(i, j)
        else:
            weightIJ = self.weight_matrix[j, i]
        self.weight_sum[i, i] += weightIJ * 0.5
        self.weight_sum[j, j] += weightIJ * 0.5
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

        # weight equation: 0.5 * (cot(alpha) + cot(beta))

        cot_theta_sum = 0
        for face in local_faces:
            other_vertex_id = face.other_point(i, j)
            vertex_o = self.verts[other_vertex_id]
            theta = omath.angle_between(vertex_i - vertex_o, vertex_j - vertex_o)
            cot_theta_sum += omath.cot(theta)
        return cot_theta_sum * 0.5

    def calculate_laplacian_matrix(self):
        # initial laplacian
        # self.laplacian_matrix = self.edge_matrix - self.neighbour_matrix
        self.laplacian_matrix = self.weight_sum - self.weight_matrix
        fixed_verts_num = len(self.fixed_verts)
        # for each constrained point, add a new row and col
        new_n = self.n + fixed_verts_num
        new_matrix = matrix((new_n, new_n), dtype=np.float)
        # Assign old values to new matrix
        new_matrix[:self.n, :self.n] = self.laplacian_matrix
        # Add 1s in the row and column associated with the fixed point to constain it
        # This will increase L by the size of fixed_verts
        for i in range(fixed_verts_num):
            new_i = self.n + i
            vert_id = self.fixed_verts[i][0]
            new_matrix[new_i, vert_id] = 1
            new_matrix[vert_id, new_i] = 1
        print(self.laplacian_matrix)

        self.laplacian_matrix = new_matrix

    def apply_deformation(self, iterations):
        print("Length of sel verts", len(self.selected_verts))

        if iterations < 0:
            iterations = self.max_iterations

        self.current_energy = 0

        # initialize b and assign constraints
        number_of_fixed_verts = len(self.fixed_verts)

        self.b_array = np.zeros((self.n + number_of_fixed_verts, 3))
         # Constraint b points
        for i in range(number_of_fixed_verts):
            self.b_array[self.n + i] = self.fixed_verts[i][1]

        # Apply following deformation iterations
        for t in range(iterations):
            print("Iteration: ", t)

            self.calculate_cell_rotations()
            self.apply_cell_rotations()
            iteration_energy = self.calculate_energy()
            print("Total Energy: ", self.current_energy)
            # if(self.energy_minimized(iteration_energy)):
            #     print("Energy was minimized at iteration", t, " with an energy of ", iteration_energy)
            #     break
            self.current_energy = iteration_energy

    def energy_minimized(self, iteration_energy):
        return abs(self.current_energy - iteration_energy)  < self.threshold

    def calculate_cell_rotations(self):
        print("Calculating Cell Rotations")
        for vert_id in range(self.n):
            rotation = self.calculate_rotation_matrix_for_cell(vert_id)
            self.cell_rotations[vert_id] = rotation

    def vert_is_deformable(self, vert_id):
        return self.vert_status[vert_id] == 1

    def precompute_p_i(self):
        self.P_i_array = []
        for i in range(self.n):
            vert_i = self.verts[i]
            neighbour_ids = self.neighbours_of(i)
            number_of_neighbours = len(neighbour_ids)

            P_i = np.zeros((3, number_of_neighbours))

            for n_i in range(number_of_neighbours):
                n_id = neighbour_ids[n_i]

                vert_j = self.verts[n_id]
                P_i[:, n_i] = (vert_i - vert_j)
            self.P_i_array.append(P_i)

    def apply_cell_rotations(self):
        print("Applying Cell Rotations")

        # Regular b points
        for i in range(self.n):
            self.b_array[i] = self.calculate_b_for(i)

        print("Printing B")
        print(self.b_array)

        p_prime = solve(self.laplacian_matrix, self.b_array)

        # self.verts = self.verts_prime

        for i in range(self.n):
            self.verts_prime[i] = p_prime[i]

        # print("p prime")
        # print(p_prime)

    def calculate_rotation_matrix_for_cell(self, vert_id):
        covariance_matrix = self.calculate_covariance_matrix_for_cell(vert_id)

        U, s, V_transpose = np.linalg.svd(covariance_matrix)

        # U, s, V_transpose
        # V_transpose_transpose * U_transpose


        rotation = V_transpose.transpose().dot(U.transpose())
        if np.linalg.det(rotation) <= 0:
            U[:0] *= -1
            rotation = V_transpose.transpose().dot(U.transpose())
        return rotation

    def calculate_covariance_matrix_for_cell(self, vert_id):
        # s_i = P_i * D_i * P_i_prime_transpose
        vert_i_prime = self.verts_prime[vert_id]

        neighbour_ids = self.neighbours_of(vert_id)
        number_of_neighbours = len(neighbour_ids)

        D_i = np.zeros((number_of_neighbours, number_of_neighbours))

        P_i =       self.P_i_array[vert_id]
        P_i_prime = np.zeros((3, number_of_neighbours))

        for n_i in range(number_of_neighbours):
            n_id = neighbour_ids[n_i]

            D_i[n_i, n_i] = self.weight_matrix[vert_id, n_id]

            vert_j_prime = self.verts_prime[n_id]
            P_i_prime[:, n_i] = (vert_i_prime - vert_j_prime)

        P_i_prime = P_i_prime.transpose()
        return P_i.dot(D_i).dot(P_i_prime)

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

    def calculate_b_for(self, i):
        b = np.zeros((1, 3))
        neighbours = self.neighbours_of(i)
        for j in neighbours:
            w_ij = self.weight_matrix[i, j] / 2.0
            r_ij = self.cell_rotations[i] + self.cell_rotations[j]
            # print(r_ij)
            p_ij = self.verts[i] - self.verts[j]
            b += (w_ij * r_ij.dot(p_ij))
        return b

    def calculate_energy(self):
        total_energy = 0
        for i in range(self.n):
            total_energy += self.energy_of_cell(i)
        return total_energy

    def power(x):
        x ** 2

    def energy_of_cell(self, i):
        neighbours = self.neighbours_of(i)
        total_energy = 0
        for j in neighbours:
            w_ij = self.weight_matrix[i, j]
            e_ij_prime = self.verts_prime[i] - self.verts_prime[j]
            e_ij = self.verts[i] - self.verts[j]
            r_i = self.cell_rotations[i]
            value = e_ij_prime - r_i.dot(e_ij)
            norm_power = np.power(value, self.POWER)
            # total_energy += w_ij * np.linalg.norm(, ord=self.POWER) ** self.POWER
            total_energy += w_ij * np.sum(norm_power)
        return total_energy

    def hex_color_for_energy(self, energy, max_energy):
        relative_energy = (energy / max_energy) * 255
        relative_energy = max(0, min(int(relative_energy), 255))
        red = hex(relative_energy)[2:]
        blue = hex(255 - relative_energy)[2:]
        if len(red) == 1:
            red = "0" + red
        if len(blue) == 1:
            blue = "0" + blue
        return "#" + red + "00" + blue

    def hex_color_array(self):
        energies = [ self.energy_of_cell(i) for i in range(self.n) ]
        max_value = np.amax(energies)
        return [ self.hex_color_for_energy(energy, max_value) for energy in energies ]

    def show_graph(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        xs = np.squeeze(np.asarray(self.verts_prime[:, 0]))
        ys = np.squeeze(np.asarray(self.verts_prime[:, 1]))
        zs = np.squeeze(np.asarray(self.verts_prime[:, 2]))
        color = self.hex_color_array()
        # Axes3D.scatter(xs, ys, zs=zs, zdir='z', s=1)#, c=None, depthshade=True, *args, **kwargs)
        ax.scatter(xs, ys, zs, c=color)
        plt.show()
# MAIN
t = time.time()
filename            = "data/02-bar-twist/00-bar-original.off"
selection_filename  = "data/02-bar-twist/bar.sel"
deformation_file    = "data/02-bar-twist/bar.def"
iterations = -1
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
if(argc > 4):
    iterations = int(sys.argv[4])

d = Deformer(filename)
d.read_file()
d.build_weight_matrix()
if len(selection_filename) > 0:
    d.read_deformation_file(deformation_file)
    d.read_selection_file(selection_filename)
d.calculate_laplacian_matrix()
d.precompute_p_i()
print("Precomputation time ", time.time() - t)
t = time.time()
d.apply_deformation(iterations)
print("Total iteration time", time.time() - t)
d.output_s_prime_to_file()
d.show_graph()
# os.system("say complete")