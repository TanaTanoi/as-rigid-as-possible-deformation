import numpy as np
import math
def angle_between(vector_a, vector_b):
    costheta = vector_a.dot(vector_b) / (np.linalg.norm(vector_a)*np.linalg.norm(vector_b))
    return math.acos(costheta)

def cot(theta):
    return math.cos(theta) / math.sin(theta)

def string_is_int(string):
    try:
        int(string)
        return True
    except ValueError:
        return False

def inf_norm(matrix):
    return np.amax(np.abs(matrix))

def rotation_matrix_between(vector_a, vector_b):
    theta = angle_between(vector_a, vector_b)

# Apply the 4x4 matrix to the 1x3 vector
def apply_rotation(rotation_matrix, vector):
    is_4_by_4 = rotation_matrix.size == 16
    if(is_4_by_4):
        vector = np.append(vector, 1)
    vector = rotation_matrix.dot(np.matrix(vector).transpose())
    vector = vector.transpose()
    if(is_4_by_4):
        # then remove last element
        vector = np.delete(vector, 3)
    return np.squeeze(np.asarray(vector))