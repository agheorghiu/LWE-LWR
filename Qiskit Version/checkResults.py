import json
from typing import *
import numpy as np

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

m = 6
n = 3
modulus = 4

#A = np.array([[0, 2, 2], [0, 1, 1], [1, 1, 0], [1, 2, 2], [2, 3, 1], [1, 2, 0]])
# A = np.array([[3,1,0],[0,3,1],[3,2,2],[2,0,0],[1,2,2],[2,3,0]])
##A = np.array([[3,3,0],[1,2,3],[1,0,1],[1,3,1],[3,1,1],[3,0,3]])
##s = np.array([[0], [1], [0]])
##e = np.array([[0], [3], [0], [0], [3], [0]])

# m = 7
A = np.array([[0, 2, 0], [3, 2, 2], [0, 3, 1], [2, 0, 2], [1, 1, 1], [0, 2, 2], [3, 3, 0]])
s = np.array([[1], [0], [1]])
e = np.array([[1], [0], [0], [0], [0], [1], [0]])

counts = {}
with open('ibmq_16_melbourne_m=7_counts') as dataFile:
# with open('Aer_simulator_m=7_counts') as dataFile:
    counts = json.load(dataFile)

results = list(counts.keys())

successTot = 0
validRuns = 0
bCount = 0

for result in results:
    split = result[::-1].split()
    b = int(split[0])
    d = []
    for i in range(1, n + 1):
        d += list(map(int, list(split[i])))
    roundRes = list(map(int, list(split[n + 1])))

    if (sum(d) > 0):
        bCount += b * counts[result]
        x = findPreimage(roundRes, A, modulus)
        if (len(x) > 0):
            if (checkEquation(x, (x + s) % modulus, d, b, modulus) or checkEquation(x, (x - s) % modulus, d, b, modulus)):
                successTot += counts[result]
        validRuns += counts[result]
        # print(validRuns.__str__() + " " + successTot.__str__())

print(successTot / validRuns)
print(bCount / validRuns)
