import numpy as np  # NumPy: Used for matrix operations and numerical computations
import time  # Used to measure how long calculations take
from scipy.linalg import expm  # SciPy's built-in matrix exponential function for comparison

"""
1. Check Norm: If the norm of matrix A is small enough, skip to step 4.

2. Scale: If necessary, scale down matrix A to improve numerical stability.

3. Compute Components: Determine the matrices N(A) and D(A).

4. Solve: Solve the linear system D(A) * e^(A) = N(A) to find e^(A).

5. Undo Scaling: If scaling was done, recover the original e^(A) by repeatedly squaring the scaled result.

"""

def factorial_custom(n):
    """
    Calculate factorial of a number (n!) using iteration.
    For example: 5! = 5 x 4 x 3 x 2 x 1 = 120
    """
    if n < 0:
        raise ValueError("Factorial is not defined for negative numbers")
    if n == 0 or n == 1:
        return 1
    result = 1
    for i in range(2, n + 1):
        result *= i
    return result

def matrix_power(matrix, p):
    """
    Calculates matrix raised to a power (like A² or A³) efficiently.
    Uses the square-and-multiply method, similar to how you might quickly calculate 2⁸
    by doing 2 → 4 → 16 → 256 instead of multiplying by 2 eight times.
    """
    result = np.eye(matrix.shape[0])
    base = matrix.copy()
    while p > 0:
        if p % 2 == 1:
            result = np.dot(result, base)
        base = np.dot(base, base)
        p //= 2
    return result

    
def pade_exp(A, q, epsilon=1e-6):
    """
    Compute the matrix exponential using Padé approximation with Lagrange remainder.
    - A: The input matrix.
    - q: The degree of the Padé approximation.
    - epsilon: The tolerance for the Lagrange remainder.
    """
    # Step 1: Normalize the matrix
    A_norm = matrix_norm_1(A)
    m = 0
    while A_norm > 0.5:
        A = A / 2
        m += 1
        A_norm = A_norm / 2

    # Step 2: Initialize identity matrix and storage for Nq and Dq
    I = np.eye(A.shape[0])
    Nq = I.copy()
    Dq = I.copy()
    sign = -1

    # Step 3: Iteratively compute the Padé approximation
    for j in range(1, q + 1):
        coef = (factorial_custom(2 * q - j) * factorial_custom(q)) / (
            factorial_custom(2 * q) * factorial_custom(j) * factorial_custom(q - j)
        )
        Aj = matrix_power(A, j)
        Nq += coef * Aj
        Dq += (coef * sign) * Aj
        sign *= -1

        # Compute the Lagrange remainder at this step
        lagrange_remainder = compute_lagrange_remainder(Nq, q)
        if lagrange_remainder < epsilon:
            print(
                f"Lagrange remainder ({lagrange_remainder:.6e}) is within the tolerance ({epsilon:.6e}). "
                f"Stopping calculation early at step {j}."
            )
            break

    # Solve the linear system to get the matrix exponential
    exp_Am = solve_linear_system(Dq, Nq)

    # Undo the scaling by repeated squaring
    for i in range(m):
        exp_Am = np.dot(exp_Am, exp_Am)

    return exp_Am

def compute_lagrange_remainder(A, q):
    """
    Compute the Lagrange remainder for the Padé approximation.
    - A: The input matrix.
    - q: The degree of the Padé approximation.
    """
    # Compute the norm of A
    A_norm = matrix_norm_1(A)
    
    # Compute the maximum coefficient for the remainder
    factorial_2q = factorial_custom(2 * q)
    max_term = factorial_2q / (factorial_custom(q) ** 2)
    
    # Compute the Lagrange remainder
    remainder = max_term * (A_norm ** (2 * q + 1) / factorial_custom(2 * q + 1))
    return remainder

def matrix_norm_1(matrix):
    """
    Finding the maximum impact any column has in the matrix help for scaling.
    matrix.shape[0] gives number of rows 
    matrix.shape[1] gives number of columns 
    Take each column of the matrix -> Take absolute value of each number in that column
    -> Sum up those absolute values -> Do this for all columns -> Find the largest sum
    """
    norm = max(sum(abs(matrix[i, j]) for i in range(matrix.shape[0])) for j in range(matrix.shape[1]))
    return norm

def solve_linear_system(A, B):
    """
    Solves a system of linear equations (AX = B) using Gaussian elimination.
    This is like solving multiple equations simultaneously but with matrices.
    pade - Rpq(A) = [Dpq(A) (A)]^-1 Npq(A) (B)
    """
    n = A.shape[0]
    X = np.zeros_like(B)
        
    # Forward elimination: Make the matrix triangular
    for i in range(n):
        # Find the best row to use as pivot (to avoid numerical problems) - max element in tow
        max_row = max(range(i, n), key=lambda x: abs(A[x, i]))
        #swap rows
        A[[i, max_row]] = A[[max_row, i]]
        B[[i, max_row]] = B[[max_row, i]]
    
        """
            # For each row j below the pivot row i:
         Scale the pivot row i by a factor (A[j, i] / A[i, i]).
         Subtract the scaled pivot row from row j to:
           1. Zero out the element A[j, i].
           2. Update the remaining elements of row j in A and the corresponding element in B.
        """        
        # Eliminate terms below the pivot
        for j in range(i + 1, n):
            factor = A[j, i] / A[i, i]
            A[j, i:] -= factor * A[i, i:]
            B[j] -= factor * B[i]

    # Back substitution: Solve for X from bottom to top
    for i in range(n - 1, -1, -1):
        # AX=B
        # 1x 2x 3x 4x | B1
        X[i] = (B[i] - np.dot(A[i, i + 1:], X[i + 1:])) / A[i, i]
    return X


# Load the matrix from a CSV file
A = np.loadtxt('exp_data.csv', delimiter=',')
q = 13  # How accurate we want our approximation to be (higher = more accurate but slower)
epsilon = 1e-6  # How small a term needs to be before we consider it negligible

# Test our custom implementation and measure its speed
print("Custom Pade Approximation Implementation:")
start_time = time.time()
exp_matrix_pade = pade_exp(A, q, epsilon)
end_time = time.time()
print(f"Time taken: {end_time - start_time:.6f} seconds")
print("Result (truncated to 5x5):\n", exp_matrix_pade[:5, :5])

# Compare with SciPy's built-in implementation
print("\nExisting Implementation (scipy.linalg.expm):")
start_time = time.time()
exp_matrix_existing = expm(A)
end_time = time.time()
print(f"Time taken: {end_time - start_time:.6f} seconds")
print("Result (truncated 5x5):\n", exp_matrix_existing[:5, :5])