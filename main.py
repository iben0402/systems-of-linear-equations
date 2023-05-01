import math
import numpy as np

# Function that multiplies two matrices if possible
# returns 0 if not possible
def multiply(X, Y):
    if X.shape[1] == Y.shape[0]:
        cols = Y.shape[1]
        rows = X.shape[0]
        C = np.zeros((X.shape[0], Y.shape[1]))
        for row in range(rows):
            for col in range(cols):
                for elt in range(len(Y)):
                    C[row, col] += X[row, elt] * Y[elt, col]
        return C
    else:
        return 0


# Function that substracts two matrix under the condition that A an B are the same shape
# returns 0 if they're not the same shape
def subtract(X, Y):
    if X.shape == Y.shape:
        C = np.zeros_like(A)
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                C[i][j] = X[i][j] - Y[i][j]
        return C
    else:
        return 0


# Function that adds two matrix under the condition that A an B are the same shape
# returns 0 if they're not the same shape
def add(X, Y):
    if X.shape == Y.shape:
        C = np.zeros_like(A)
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                C[i][j] = X[i][j] + Y[i][j]
        return C
    else:
        return 0


# Function that calculates the norm of the residuum vector
def residuum_norm(A, x, b):
    C = multiply(A, x)
    C = subtract(C, b)
    norm = 0
    for i in range(C.shape[0]):
        norm += C[i][0] ** 2

    norm = math.sqrt(norm)
    return norm


if __name__ == '__main__':

    # inicialization of variables
    N = 9 * 9 * 3
    a1 = 5 + 5
    a2 = a3 = -1
    f = 8
    A = np.zeros(shape=(N, N))
    b = np.zeros(shape=(N, 1))
    x = np.zeros(shape=(N, 1))

    for i in range(N):
        A[i][i] = a1
        if i + 1 < 243:
            A[i + 1][i] = a2
            A[i][i + 1] = a2
            if i + 2 < 243:
                A[i + 2][i] = a3
                A[i][i + 2] = a3
        b[i][0] = math.sin(i * (f + 1))

    # Jacobi



