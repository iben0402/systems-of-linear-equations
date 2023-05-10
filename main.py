import math
import time
import matplotlib.pyplot as plt


def createMatrix(sizeX, sizeY):
    M = [[0 for i in range(sizeX)] for j in range(sizeY)]
    return M


def copyMatrix(mat):
    new_matrix = createMatrix(len(mat[0]), len(mat))
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            new_matrix[i][j] = mat[i][j]

    return new_matrix

# Function that multiplies two matrices if possible
# returns 0 if not possible
def multiply(X, Y):
    if len(X[0]) == len(Y):
        cols = len(Y[0])
        rows = len(X)
        C = createMatrix(cols, rows)
        for row in range(rows):
            for col in range(cols):
                for elt in range(len(Y)):
                    C[row][col] += X[row][elt] * Y[elt][col]
        return C
    else:
        return 0


def multiply_by_const(mat, const):
    mult = createMatrix(len(mat[0]), len(mat))
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            mult[i][j] = const * mat[i][j]

    return mult

# Function that substracts two matrix under the condition that A an B are the same shape
# returns 0 if they're not the same shape
def subtract(X, Y):
    if len(X) ==len(Y) and len(X[0])==len(Y[0]):
        C = createMatrix(len(X[0]), len(X))
        for y in range(len(X)):
            for x in range(len(X[0])):
                C[y][x] = X[y][x] - Y[y][x]
        return C
    else:
        return 0


# Function that adds two matrix under the condition that A an B are the same shape
# returns 0 if they're not the same shape
def add(X, Y):
    if len(X) ==len(Y) and len(X[0])==len(Y[0]):
        C = createMatrix(len(X[0]), len(X))
        for i in range(len(X)):
            for j in range(len(X[0])):
                C[i][j] = X[i][j] + Y[i][j]
        return C
    else:
        return 0



# Function that calculates the norm of the residuum vector
def residuum_norm(A, x, b):
    C = multiply(A, x)
    C = subtract(C, b)
    norm = 0
    for i in range(len(C)):
        norm += C[i][0] ** 2

    norm = math.sqrt(norm)
    return norm

def J_getL(A):
    L = createMatrix(len(A), len(A))
    for i in range(len(A)):
        for j in range(i):
            L[i][j] = A[i][j]

    return L

def J_getU(A):
    U = createMatrix(len(A), len(A))
    for i in range(len(A)-1, -1, -1):
        for j in range(i+1, len(A)):
            U[i][j] = A[i][j]

    return U

def J_getD(A):
    D = createMatrix(len(A), len(A))
    for i in range(len(A)):
        D[i][i] = A[i][i]

    return D

def jacobi(A, b):
    counter = 0
    x = createMatrix(1, len(A))
    x_old = copyMatrix(x)
    a = copyMatrix(A)
    B = copyMatrix(b)
    residuum = []

    start = time.time()
    while residuum_norm(A, x_old, b) >= 10**(-9):
        for i in range(len(A)):
            sigma = 0
            for j in range(len(A)):
                if j != i:
                    sigma += a[i][j] * x_old[j][0]
            x[i][0] = (B[i][0] - sigma) / a[i][i]

        for i in range(len(A)):
            x_old[i][0] = x[i][0]
        residuum.append(residuum_norm(A, x, b))
        counter += 1
    end = time.time()

    print("Liczba iteracji dla metody Jacobiego: ", counter)
    print("Czas wykonania algorytmu Jacobiego: ", (end-start))

    # tworzenie wykresu błedu dla kolejnych iteracji
    # plt.plot(range(1, counter+1), residuum)
    # plt.ylabel("Wartość błędu resydualnego")
    # plt.yscale("log")
    # plt.xlabel("numer iteracji")
    # plt.title("Wartość błędu resudualnego dla kolejnych iteracji dla algorytmu Jacobiego")
    # plt.show()

    return (end-start)


def gauss_seidel(A, b):
    counter = 0
    x = createMatrix(1, len(A))
    x_old = copyMatrix(x)
    a = copyMatrix(A)
    B = copyMatrix(b)
    residuum = []


    start = time.time()
    while residuum_norm(A, x_old, b) >= 10**(-9):
        for i in range(len(A)):
            sigma = 0
            for j in range(len(A)):
                if j != i:
                    sigma += a[i][j] * x_old[j][0]
            x[i][0] = (B[i][0] - sigma) / a[i][i]
            x_old[i][0] = x[i][0]
        residuum.append(residuum_norm(A, x, b))
        counter += 1
    end = time.time()

    print("Liczba iteracji dla metody Gaussa-Seidla: ", counter)
    print("Czas wykonania algorytmu Gaussa-Seidla: ", (end-start))

    # tworzenie wykresow bledu resydualnego
    # plt.plot(range(1, counter+1), residuum)
    # plt.ylabel("Wartość błędu resydualnego")
    # plt.yscale("log")
    # plt.xlabel("numer iteracji")
    # plt.title("Wartość błędu resudualnego dla kolejnych iteracji dla algorytmu Gaussa-Seidla")
    # plt.show()

    return (end-start)


def LU_factorization(A, b):
    L = createMatrix(len(A[0]), len(A))
    for i in range(len(A)):
        L[i][i] = 1
    U = copyMatrix(A)

    start = time.time()
    for i in range(len(A)-1):
        for j in range(i+1, len(A)):
            L[j][i] = U[j][i]/U[i][i]
            for k in range(i, len(A)):
                U[j][k] = U[j][k]-L[j][i]*U[i][k]


    # solving Ly = b with forward substitution
    y = createMatrix(1, len(A))
    for i in range(len(y)):
        y[i][0] = b[i][0]
        for j in range(i):
            y[i][0]-=(L[i][j]*y[j][0])
        y[i][0] = y[i][0]/L[i][i]


    # solving Ux = y with backward substitution
    x = createMatrix(1, len(A))
    for i in range(len(A)-1, -1, -1):
        x[i][0] = y[i][0]
        for j in range(i+1, len(A)):
            x[i][0] -= U[i][j] * x[j][0]
        x[i][0] = x[i][0]/U[i][i]
    end = time.time()
    norm = residuum_norm(A, x, b)
    print("Norma z residuum dla faktoryzacji LU wynosi: ", norm)
    return (end-start)


def createDiagrams():
    N = [100, 500, 1000, 2000, 3000]
    times_for_jacobi = []
    times_for_gauss_seidel = []
    times_for_LU = []

    a1 = 5 + 5
    a2 = a3 = -1
    f = 8

    for k in range(len(N)):
        A = createMatrix(N[k], N[k])
        b = createMatrix(1, N[k])
        for i in range(N[k]):
            A[i][i] = a1
            if i + 1 < N[k]:
                A[i + 1][i] = a2
                A[i][i + 1] = a2
                if i + 2 < N[k]:
                    A[i + 2][i] = a3
                    A[i][i + 2] = a3
            alpha = i * (f + 1)
            b[i][0] = math.sin(alpha)

        times_for_jacobi.append(jacobi(A, b))
        times_for_gauss_seidel.append(gauss_seidel(A, b))
        times_for_LU.append(LU_factorization(A, b))

    plt.plot(N, times_for_jacobi)
    plt.ylabel("Czas wykonania algorytmu")
    plt.xlabel("Rozmiar N macierzy")
    plt.title("Wykres czasu wykonania algorytmu Jacobiego w zależności od rozmiaru macierzy")
    plt.show()

    plt.plot(N, times_for_gauss_seidel)
    plt.ylabel("Czas wykonania algorytmu")
    plt.xlabel("Rozmiar N macierzy")
    plt.title("Wykres czasu wykonania algorytmu Gaussa-Seidla w zależności od rozmiaru macierzy")
    plt.show()


    plt.plot(N, times_for_LU)
    plt.ylabel("Czas wykonania algorytmu")
    plt.xlabel("Rozmiar N macierzy")
    plt.title("Wykres czasu wykonania algorytmu faktoryzacji LU w zależności od rozmiaru macierzy")
    plt.show()
    print("\n")


if __name__ == '__main__':

    # inicialization of variables
    N = 993
    # N = 913
    a1 = 5 + 5
    a2 = a3 = -1
    f = 8
    A = createMatrix(N, N)
    b = createMatrix(1, N)

    for i in range(N):
        A[i][i] = a1
        if i + 1 < N:
            A[i + 1][i] = a2
            A[i][i + 1] = a2
            if i + 2 < N:
                A[i + 2][i] = a3
                A[i][i + 2] = a3
        alpha = i*(f+1)
        b[i][0] = math.sin(alpha)


    # ZADANIE B
    jacobi(A, b)
    gauss_seidel(A, b)

    # ZADANIE C
    for i in range(N):
        A[i][i] = 3

    print("Changed values of A matrix")
    jacobi(A, b)
    gauss_seidel(A, b)

    # ZADANIE D
    for i in range(N):
        A[i][i] = 3
    LU_factorization(A,b)
    print("\n\n\n")

    # ZADANIE E
    createDiagrams()
