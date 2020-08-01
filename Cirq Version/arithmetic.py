import cirq
import numpy as np
from typing import *

# depth 3 QFT on 2 qubits
# reg is taken to be a register (list) of qubits
def qft(reg):
    yield cirq.H(reg[0])
    cs = cirq.S.controlled()
    yield cs(reg[1], reg[0])
    yield cirq.H(reg[1])

# inverse QFT on register
def iqft(reg):
    inverse_circuit = cirq.inverse(qft(reg))
    return inverse_circuit

# maps |x>|y> to |cx + y (mod modulus)>
# x, y are lists of qubits and c, modulus are ints
def addRegisters(x, y, c, modulus):
    n = len(x)
    for i in range(n):
        tot = ((2 ** i) * c) % (modulus)
        for j in range(n):
            power = 2 ** ((n - 1) - j)
            if ((tot % power) != 0 or ((tot % power) == 0 and ((tot / power) % 2) != 0)):
                crz = cirq.rz(np.pi * tot / power).controlled()
                yield crz(x[i], y[j])

# maps |x> to |x + c (mod modulus)>
# res is a list of qubits and c, modulus are ints
def addConstantToRegister(c, modulus, res):
    n = len(res)
    for i in range(n):
        power = 2 ** ((n - 1) - i)
        if ((c % power) != 0 or ((c % power) == 0 and ((c / power) % 2) != 0)):
            yield cirq.rz(c * np.pi / power)(res[i])

# controlled version of previous function (controlled on b)
# b is a single qubit, res is a list of qubits, and c, modulus are ints
def controlledAddConstant(b, c, modulus, res):
    n = len(res)
    for i in range(n):
        power = 2 ** ((n - 1) - i)
        if ((c % power) != 0 or (c % power) == 0 and ((c / power) % 2) != 0):
            crz = cirq.rz(c * np.pi / power).controlled()
            yield crz(b, res[i])

# coherently adds <c, x> + b . y mod modulus to result register
# b is a single qubit, res is a list of qubits; y, modulus are ints
# c is a list of ints, and x is a list of lists of qubits
def coherentLWE(b, x, res, y, c, modulus):
    n = len(c)
    for j in range(n):
        if (c[j] != 0):
            yield addRegisters(x[j], res, c[j], modulus)
    yield controlledAddConstant(b, y, modulus, res)

# inverse of previous function
def adjointCoherentLWE(b, x, res, y, c, modulus):
    inverse_circuit = cirq.inverse(coherentLWE(b, x, res, y, c, modulus))
    return inverse_circuit

# ------------------------------------------------------------------------------
# Checking Results Functions

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
