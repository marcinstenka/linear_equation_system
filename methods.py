import math
import matplotlib.pyplot as plt

from matrix import Matrix


def Jacobi(A, b, N, threshold=10 ** -9, maxIterations=300):
    D = A.getDiagonal()
    U = A.getUpper()
    L = A.getLower()
    counter = 0
    tempx = [[1 for x in range(1)] for y in range(N)]
    x = Matrix(1, N, tempx)
    comp1 = L + U
    comp1.signChange()
    end = True

    while counter < maxIterations and end:
        counter = counter + 1
        comp2 = comp1 * x
        comp3 = comp2 + b
        x = D.mldivide(comp3)
        res = A * x - b
        resNormalized = norm(res)

        if resNormalized <= float(threshold):
            end = False

    print("Jacobi iterations: " + str(counter))
    print("Jacobi norm: " + str(resNormalized))
    return x


def Gauss_Seidl(A, b, N, threshold=10 ** -9, maxIterations=300):
    D = A.getDiagonal()
    U = A.getUpper()
    L = A.getLower()
    counter = 0
    tempx = [[1 for n in range(1)] for y in range(N)]
    x = Matrix(1, N, tempx)
    end = True

    while counter < maxIterations and end:
        counter = counter + 1
        for i in range(N):
            x[i][0] = (b[i][0] - sum(U[i][j] * x[j][0] for j in range(i, N)) - sum(L[i][j] * x[j][0] for j in range(i))) / D[i][i]

        res = A * x - b
        resNormalized = norm(res)
        if resNormalized <= float(threshold):
            end = False
    print("Gauss-Sedil iterations: " + str(counter))
    print("Gauss-Sedil norm: "+ str(resNormalized))
    return x


def LU_decomposition(A, N, a1, a2, a3):
    tempL = [[0 for x in range(N)] for y in range(N)]
    for i in range(N):
        for j in range(N):
            if i == j:
                tempL[i][j] = 1
            if i == j + 1 or j == i + 1 or i == j + 2 or j == i + 2:
                tempL[i][j] = 0
    tempU = [[0 for x in range(N)] for y in range(N)]
    for i in range(N):
        for j in range(N):
            if i == j:
                tempU[i][j] = a1
            if i == j + 1 or j == i + 1:
                tempU[i][j] = a2
            if i == j + 2 or j == i + 2:
                tempU[i][j] = a3
    L = Matrix(N, N, tempL)
    U = Matrix(N, N, tempU)

    for i in range(N):
        for j in range(i + 1, N):
            if U[i][i] == 0:
                U[i + 1], U[i] = U[i], U[i + 1]
                L[j][i] = 0
            else:
                L[j][i] = U[j][i] / U[i][i]
            for k in range(i, N):
                U[j][k] = U[j][k] - (L[j][i] * U[i][k])
    return L, U


def factorizationLU(A, b, N, a1, a2, a3):
    L, U = LU_decomposition(A, N, a1, a2, a3)
    tempy = [[0 for x in range(1)] for y in range(N)]
    tempx = [[0 for x in range(1)] for y in range(N)]
    y = Matrix(1, N, tempy)
    x = Matrix(1, N, tempx)
    for i in range(N):
        sum_L = sum(L[i][k] * y[k][0] for k in range(i))
        y[i][0] = (b[i][0] - sum_L)
    for i in reversed(range(N)):
        sum_U = sum(U[i][j] * x[j][0] for j in range(i + 1, N))
        x[i][0] = (y[i][0] - sum_U) / U[i][i]

    res = A * x - b

    resNormalized = norm(res)

    print("Norma residuum " + str(resNormalized))
    return x

def norm(res):
    sum = 0
    for i in range(res.Y):
        sum += res[i][0] * res[i][0]
    result = math.sqrt(sum)
    return result
