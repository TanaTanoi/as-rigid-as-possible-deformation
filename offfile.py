class OffFile:
    def __init__(self, filename):
        self.lines = open(filename, 'r').read().split("\n")
        #Remove the OFF in the first line
        if 'OFF' in self.lines[0]:
            self.lines.pop(0)

    # Grabs the next line, ignoring comments
    def nextLine(self):
        line = self.lines.pop(0)
        while(len(line) == 0 or line[0] == '#'):
            line = self.lines.pop(0)
        return line