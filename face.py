class Face:
    def __init__(self, v1, v2, v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3

    def containsPointIDs(self, i, j):
        return i in vertexIDs and j in vertexIDs

    def vertexIDs():
        return [self.v1, self.v2, self.v3]