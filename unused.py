from main import *

def determinant(A):
    mat = createMatrix(len(A), len(A))
    for i in range(len(A)):
        for j in range(len(A)):
            mat[i][j] = A[i][j]
    n = len(mat)
    temp = [0] * n  # temporary array for storing row
    total = 1
    det = 1  # initialize result

    # loop for traversing the diagonal elements
    for i in range(0, n):
        index = i  # initialize the index

        # finding the index which has non zero value
        while (index < n and mat[index][i] == 0):
            index += 1

        if (index == n):  # if there is non zero element
            # the determinant of matrix as zero
            continue

        if (index != i):
            # loop for swapping the diagonal element row and index row
            for j in range(0, n):
                mat[index][j], mat[i][j] = mat[i][j], mat[index][j]

            # determinant sign changes when we shift rows
            # go through determinant properties
            det = det * int(pow(-1, index - i))

        # storing the values of diagonal row elements
        for j in range(0, n):
            temp[j] = mat[i][j]

        # traversing every row below the diagonal element
        for j in range(i + 1, n):
            num1 = temp[i]  # value of diagonal element
            num2 = mat[j][i]  # value of next row element

            # traversing every column of row
            # and multiplying to every row
            for k in range(0, n):
                # multiplying to make the diagonal
                # element and next row element equal

                mat[j][k] = (num1 * mat[j][k]) - (num2 * temp[k])

            total = total * num1  # Det(kA)=kDet(A);

    # multiplying the diagonal elements to get determinant
    for i in range(0, n):
        det = det * mat[i][i]

    return int(det / total)


def transpose(mat):
    N = len(mat)
    transposed = createMatrix(N, N)
    for i in range(N):
        for j in range(N):
            transposed[i][j] = mat[j][i]

    return transposed

def matrixMinor(m,i,j):
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]


# returns inverse matrix of A
def inverse_matrix(mat):
    det = determinant(mat)
    if len(mat) == 2:
        return [[mat[1][1]/det, -1*mat[0][1]/det],
                [-1*mat[1][0]/det, mat[0][0]/det]]

    cofactors = []
    for r in range(len(mat)):
        cofactorRow = []
        for c in range(len(mat)):
            minor = matrixMinor(mat, r, c)
            cofactorRow.append(((-1)**(r+c)) * determinant(minor))
        cofactors.append(cofactorRow)
    cofactors = transpose(cofactors)

    cofactors = multiply_by_const(cofactors, 1/det)
    return cofactors


def inverse_diagonal(mat):
    inversed = createMatrix(len(mat), len(mat))
    for i in range(len(mat)):
        inversed[i][i] = 1/mat[i][i]

    return inversed


# Implementation of Jacobi method, prints number of iterations and time
# tries to get a solution with residuum norm smaller than 10^-9
def jacobi_and_gauss(A, x, b):
    # Jacobi method
    counter = 0

    L = J_getL(A)
    U = J_getU(A)
    D = J_getD(A)

    # calculating N = D^-1
    N = inverse_diagonal(D)

    # calculating O = D^-1 * b
    O = multiply(N, b)

    # calculating M = D^-1 (L + U)
    M = add(L, U)
    M = multiply(multiply_by_const(N, -1), M)
    norm = residuum_norm(A, x, b)

    start = time.time()
    while norm > 10**(-9):
        x = add(multiply(M, x), O)
        counter+=1
        norm = residuum_norm(A, x, b)

    end = time.time()

    print("liczba iteracji dla metody Jacobiego: ", counter)
    print("Czas wykonania algorytmu Jacobiego: ", (end - start))

    # Gauss-Seidel Method

    x2 = createMatrix(1, len(x))

    P = subtract(D, L)
    P = inverse_matrix(P)

    R = multiply(P, b)

    norm = residuum_norm(A, x2, b)
    counter = 0
    start = time.time()
    while norm > 10**(-9):
        S = multiply(U, x2)
        S = multiply(P, S)
        x2 = add(S, R)
        counter += 1
        norm = residuum_norm(A, x2, b)
    end = time.time()
    print("liczba iteracji dla metody Gaussa-Seidla: ", counter)
    print("Czas wykonania algorytmu Gaussa-Seidla: ", (end - start))