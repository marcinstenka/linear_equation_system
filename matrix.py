class Matrix:
    def __init__(self, X, Y, values):
        self.X = X
        self.Y = Y
        self.values = values

    def getDiagonal(self):
        diagonal = [[0 for x in range(self.X)] for y in range(self.Y)]
        for i in range(self.X):
            for j in range(self.Y):
                if i == j:
                    diagonal[i][j] = self.values[i][j]
        return Matrix(self.X, self.Y, diagonal)

    def getLower(self):
        lower = [[0 for x in range(self.X)] for y in range(self.Y)]
        for i in range(self.X):
            for j in range(self.Y):
                if i >= j+1:
                    lower[i][j] = self.values[i][j]
        return Matrix(self.X, self.Y, lower)

    def getUpper(self):
        upper = [[0 for x in range(self.X)] for y in range(self.Y)]
        for i in range(self.X):
            for j in range(self.Y):
                if j >= i+1:
                    upper[i][j] = self.values[i][j]
        return Matrix(self.X, self.Y, upper)

    def signChange(self):
        for i in range(self.Y):
            for j in range(self.X):
                self.values[i][j] = -1 * self.values[i][j]
        return Matrix(self.X, self.Y, self.values)

    def inverse(self):
        for i in range(self.X):
            for j in range(self.Y):
                if self.values[i][j] != 0:
                    self.values[i][j] = 1 / self.values[i][j]

    def mldivide(self, other):
        temp = [[0 for x in range(self.X)] for y in range(self.Y)]
        for x in range(self.X):
            for y in range(self.Y):
                if x == y:
                    temp[x][y] = (1 / self.values[x][y])
        return Matrix(self.X, self.Y, temp) * other

    def getZeros(self, X, Y):
        return [[0 for x in range(X)] for y in range(Y)]

    def __str__(self):
        print('\n'.join([' '.join(['{:4}'.format(item) for item in row]) for row in self.values]))
        return ""

    def __getitem__(self, x):
        return self.values[x]

    def __add__(self, other):
        temp = Matrix.getZeros(self, self.X, self.Y)
        if self.X == other.X:
            for i in range(self.Y):
                for j in range(self.X):
                    temp[i][j] = self.values[i][j] + other.values[i][j]
        return Matrix(self.X, self.Y, temp)

    def __sub__(self, other):
        temp = Matrix.getZeros(self, self.X, self.Y)
        if self.X == other.X:
            for i in range(self.Y):
                for j in range(self.X):
                    temp[i][j] = self.values[i][j] - other.values[i][j]

        return Matrix(self.X, self.Y, temp)

    def __mul__(self, other):
        temp = Matrix.getZeros(self, other.X, other.Y)
        for i in range(other.Y):
            for j in range(other.X):
                value = 0
                for k in range(self.Y):
                    value = value + self.values[i][k] * other.values[k][j]
                temp[i][j] = value
        return Matrix(other.X, other.Y, temp)
