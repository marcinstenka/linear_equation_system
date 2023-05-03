import math
import time
import matplotlib.pyplot as plt
from matrix import Matrix
from methods import Jacobi, Gauss_Seidl, factorizationLU


N = 985 # 985
a1 = 5 + 9
a2 = a3 = -1
tempA = [[0 for x in range(N)] for y in range(N)]
for i in range(N):
    for j in range(N):
        if i == j:
            tempA[i][j] = a1
        if i == j + 1 or j == i + 1:
            tempA[i][j] = a2
        if i == j + 2 or j == i + 2:
            tempA[i][j] = a3

tempb = [[0 for x in range(1)] for y in range(N)]
for n in range(N):
    tempb[n][0] = math.sin(n * 9)


# Task A
A = Matrix(N, N, tempA)
b = Matrix(1, N, tempb)

# Task B
start = time.time()
xJacobi = Jacobi(A, b, N)
end = time.time()
timeJacobi = end - start
print("Jacobi method time: " + str(timeJacobi))

start = time.time()
xGauss_Seidl = Gauss_Seidl(A, b, N)
end = time.time()
timeGauss_Seidl = end - start
print("Gauss_Seidl method time: " + str(timeGauss_Seidl))

# Task C
N = 985 # 985
a1 = 3
a2 = a3 = -1
tempA = [[0 for x in range(N)] for y in range(N)]
for i in range(N):
    for j in range(N):
        if i == j:
            tempA[i][j] = a1
        if i == j + 1 or j == i + 1:
            tempA[i][j] = a2
        if i == j + 2 or j == i + 2:
            tempA[i][j] = a3

tempb = [[0 for x in range(1)] for y in range(N)]
for n in range(N):
    tempb[n][0] = math.sin(n * 9)
A = Matrix(N, N, tempA)
b = Matrix(1, N, tempb)

xJacobi = Jacobi(A, b, N)
xGauss_Seidl = Gauss_Seidl(A, b, N)

# Task D
xFact = factorizationLU(A, b, N, a1, a2, a3)

# Task E
sizes = [100, 500, 1000, 2000, 3000]
plots = [[0 for x in range(5)] for y in range(5)]
N = 985 # 985
a1 = 14
a2 = a3 = -1
number = 0
for N in sizes:
    tempA = [[0 for x in range(N)] for y in range(N)]
    for i in range(N):
        for j in range(N):
            if i == j:
                tempA[i][j] = a1
            if i == j + 1 or j == i + 1:
                tempA[i][j] = a2
            if i == j + 2 or j == i + 2:
                tempA[i][j] = a3

    tempb = [[0 for x in range(1)] for y in range(N)]
    for n in range(N):
        tempb[n][0] = math.sin(n * 9)
    A = Matrix(N, N, tempA)
    b = Matrix(1, N, tempb)


    start = time.time()
    x = Jacobi(A, b, N)
    end = time.time()
    JacobiTime = end-start
    plots[0][number] = JacobiTime
    print("Jacobi " + str(JacobiTime))

    start = time.time()
    x = Gauss_Seidl(A, b, N)
    end = time.time()
    Gauss_SeidlTime = end - start
    plots[1][number] = Gauss_SeidlTime
    print("Gauss-Seidl: " + str(Gauss_SeidlTime))

    start = time.time()
    x = factorizationLU(A, b, N, a1, a2, a3)
    end = time.time()
    factorizationLUTime = end - start
    plots[2][number] = factorizationLUTime
    print("factLU time:" + str(factorizationLUTime))
    number = number + 1

plt.plot(sizes, plots[0], 'r', label='Jacobi')
plt.plot(sizes, plots[1], 'g', label='Gauss-Seidel')
plt.plot(sizes, plots[2], 'b', label='Factorization LU')
plt.xticks(sizes)
plt.title("Comparison of time used by algorithms while working with various data size")
plt.legend()
plt.xlabel("Size of matrix [NxN]")
plt.ylabel("Time [s]")
plt.show()