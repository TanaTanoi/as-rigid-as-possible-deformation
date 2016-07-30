class Vertex:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class Face:
    def __init__(self, v1, v2, v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3

class OffFileReader:
    def __init__(self, filename):
        self.lines = open(filename, 'r').read().split("\n")
        #Remove the OFF in the first line
        if 'OFF' in self.lines[0]:
            self.lines.pop(0)

    # Grabs the next line, ignoring comments
    def nextLine(self):
        line = self.lines.pop(0)
        while(line[0] == '#'):
            line = self.lines.pop(0)
        return line

# main

filename = "02-bar-twist/00-bar-original.off"

f = OffFileReader(filename)

first_line = f.nextLine().split()

number_of_verticies =   int(first_line[0])
number_of_faces =       int(first_line[1])
number_of_edges =       int(first_line[2])

verts = []
faces = []
for i in range(0, number_of_verticies):
    vert_line = f.nextLine().split()
    x = float(vert_line[0])
    y = float(vert_line[1])
    z = float(vert_line[2])
    verts.append(Vertex(x, y, z))

print(str(len(verts)) + " verticies")

for i in range(0, number_of_faces):
    faceLine = f.nextLine().split()
    v1_id = int(faceLine[1])
    v2_id = int(faceLine[2])
    v3_id = int(faceLine[3])
    faces.append(Face(v1_id, v2_id, v3_id))

print(str(len(faces)) + " faces")
print(str(number_of_edges) + " edges")