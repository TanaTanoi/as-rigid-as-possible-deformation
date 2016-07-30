import numpy as np
import math
def angleBetween(vector_a, vector_b):
    costheta = vector_a.dot(vector_b) / (np.linalg.norm(vector_a)*np.linalg.norm(vector_b))
    return math.acos(costheta)

def cot(theta):
    return math.cos(theta) / math.sin(theta)