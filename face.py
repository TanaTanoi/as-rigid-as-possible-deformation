class Face:
    def __init__(self, v1, v2, v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3

    def containsPointIDs(self, i, j):
        vIDs = self.vertexIDs()
        return i in vIDs and j in vIDs

    def otherPoint(self, i, j):
        vIDs = set(self.vertexIDs())
        vIDs.remove(i)
        vIDs.remove(j)
        return vIDs.pop()

    def vertexIDs(self):
        return [self.v1, self.v2, self.v3]