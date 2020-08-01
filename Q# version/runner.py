import qsharp
import numpy as np
from typing import *
from RSPRounding import RSPMain

# generates m-dimensional low weight error vector with entries modulo modulus
def generateErrVector (seed : int, m : int, p : float, modulus : int) -> List[int]:
    np.random.seed(seed)
    e = np.zeros([m, 1], dtype=int)
    for i in range(m):
        if (np.random.random() > p):
            rand = np.random.randint(2)
            if (rand == 1):
                rand = modulus - 1
            e[i] = rand
    return e

# computes y = As + e
# A is m x n, s is an n-dimensional vectors (n x 1), e is an m-dimensional vector (m x 1)
def LWEfunEval (A : List[List[int]], s : List[int], e : List[int], modulus : int) -> List[int]:
    y = (np.dot(A, s) + e) % modulus
    return y

# generates a random (A, s) with A having entries mod modulus and s having entries mod 2
# A is m x n, s is an n-dimensional vector (n x 1)
def generateRandomLWEInstance (seed : int, m : int, n : int, modulus : int) -> (List[List[int]], List[int]):
    np.random.seed(seed)
    A = np.random.randint(modulus, size=[m, n])
    s = np.random.randint(2, size=[n, 1])
    return (A, s)

# combines previous three functions to output (A, s, y = As + e) 
def generateRandomLWESample (seed : int, m : int, n : int, modulus : int) -> (List[List[int]], List[int], List[int]):
    (A, s) = generateRandomLWEInstance(seed, m, n, modulus)
    e = generateErrVector(seed, m, 0.1, modulus)
    y = LWEfunEval(A, s, e, modulus)
    return (A, s, y)



# converting between binary and int
def toBinary(num: int, numBits: int) -> List[int]:
    return list(reversed([int((num & (2**i)) > 0) for i in range(numBits)]))

# converts entries in vector to binary and concatenates results
def convertVectToBinary (x : List[int], modulus : int) -> List[int]:
    bits = int(np.log2(modulus))
    lst = []
    n = len(x)
    for i in range(n):
        lst += reversed(toBinary(x[i], bits)) # reversed because of little endian convention
    return lst

# checks if d . (x1 xor x2) = b
def checkEquation (x1 : List[int], x2: List[int], d : List[int], b : int, modulus : int) -> bool:
    x1bin = convertVectToBinary(x1, modulus)
    x2bin = convertVectToBinary(x2, modulus)
    x1xorx2 = [(a ^ b) for (a, b) in zip(x1bin, x2bin)]
    ds = np.dot(d, x1xorx2) % 2
    return (ds == b)

# converts from num in base 10 to base given
def toBase(num : int, numBits : int, base : int) -> List[int]:
    if (num == 0):
        return [0] * numBits
    digits = []
    while (num > 0):
        digits.append(num % base)
        num //= base
    if (len(digits) < numBits):
        digits += [0] * (numBits - len(digits))
    return digits[::-1]

# finds preimage x such that [Ax mod modulus] = y
def findPreimage (y: List[int], A: List[List[int]], modulus: int) -> List[int]:
    n = len(A[0])
    for i in range(0, modulus ** n):
        x = np.array(toBase(i, n, modulus))
        res = (np.dot(A, x) % modulus) >= modulus / 2
        if ((res == y).all()):
            return np.array(list(map(lambda e: [e], x)))
    return []


# parameters
seed = 2000
m = 7
n = 3
modulus = 4
rnd = 2

# (A, s, y) = generateRandomLWESample(seed, m, n, modulus)

# A = np.array([[0, 2, 2], [0, 1, 1], [1, 1, 0], [1, 2, 2], [2, 3, 1], [1, 2, 0]])
# A = np.array([[0,1,2],[3,0,0],[1,1,3],[1,0,2],[1,0,0],[0,3,0]])
A = np.array([[0, 2, 0], [3, 2, 2], [0, 3, 1], [2, 0, 2], [1, 1, 1], [0, 2, 2], [3, 3, 0]])
s = np.array([[1], [0], [1]])
e = np.array([[1], [0], [0], [0], [0], [1], [0]])
y = LWEfunEval(A, s, e, modulus)

successTot = 0
validRuns = 0
numRuns = 1000
bCount = 0

for run in range(numRuns):
    (b, d, roundRes) = RSPMain.simulate(A=A.tolist(), y=(list(map(lambda x: x[0], y.tolist()))), modulus=modulus, rnd=rnd)
    if (sum(d) > 0):
        print(d)
        x = findPreimage(roundRes, A, modulus)
        if (len(x) > 0):
            successTot += (checkEquation(x, (x + s) % modulus, d, b, modulus) or checkEquation(x, (x - s) % modulus, d, b, modulus))
        validRuns += 1
        print(validRuns.__str__() + " " + successTot.__str__())

print(successTot / validRuns)