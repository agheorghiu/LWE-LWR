from arithmetic import *
from cirq.contrib.svg import SVGCircuit
# from cirq.contrib.qcircuit import circuit_to_latex_using_qcircuit
from cirq.contrib.svg import circuit_to_svg

# computes y = As + e
# A is m x n, s is an n-dimensional vectors (n x 1), e is an m-dimensional vector (m x 1)
def LWEfunEval (A : List[List[int]], s : List[int], e : List[int], modulus : int) -> List[int]:
    y = (np.dot(A, s) + e) % modulus
    return y

# parameters
seed = 2000
m = 7
n = 3
modulus = 4
bits = int(np.ceil(np.log2(modulus)))

A = np.array([[0, 2, 0], [3, 2, 2], [0, 3, 1], [2, 0, 2], [1, 1, 1], [0, 2, 2], [3, 3, 0]])
A = A.tolist()
s = np.array([[1], [0], [1]])
e = np.array([[1], [0], [0], [0], [0], [1], [0]])
y = LWEfunEval(A, s, e, modulus)
y = (list(map(lambda w: w[0], y.tolist())))

nQubits = bits * n + bits + m

# initialize registers
b = cirq.LineQubit(0)
x = []
for i in range(n):
    x.append([cirq.LineQubit(i) for i in range(1 + bits * i, bits * (i + 1) + 1)])

anc = [cirq.LineQubit(i) for i in range(1 + bits * n, (n + 1) * bits + 1)]

r = [cirq.LineQubit(i) for i in range((n + 1) * bits + 1, bits * n + bits + m)]
r += [anc[0]]

# initialize circuit
circuit = cirq.Circuit()

# prepare superpositions over b and x
# have \sum_{b, x} |b>|x>|0>
circuit.append(cirq.H(b))
for i in range(n):
    for j in range(bits):
        circuit.append(cirq.H(x[i][j]))

# put ancilla in phase
for i in range(bits):
    circuit.append(cirq.H(anc[i]))

# compute rounded result
# use ancilla to store pre-rounded result; clear ancilla each time
# final state should be |b>|x>|0>|round(Ax + b . (As + e'))>
for i in range(m - 1):
    circuit.append(coherentLWE(b, x, anc, y[i], A[i], modulus))
    
    circuit.append(iqft(anc))

    circuit.append(cirq.CNOT(anc[0], r[i]))

    circuit.append(qft(anc))

    circuit.append(adjointCoherentLWE(b, x, anc, y[i], A[i], modulus))

# last iteration MSB of ancilla is rounded result
circuit.append(coherentLWE(b, x, anc, y[m - 1], A[m - 1], modulus))
circuit.append(iqft(anc))

# get an equation in s
circuit.append(cirq.H(b))
for i in range(n):
    for j in range(bits):
        circuit.append(cirq.H(x[i][j]))

circuit.append(cirq.measure(b, key='bMeas'))
for i in range(n):
    circuit.append(cirq.measure(*x[i], key='d{}'.format(i)))
circuit.append(cirq.measure(*r, key='roundRes'))

numRuns = 1000
simulator = cirq.Simulator()
result = simulator.run(circuit, repetitions=numRuns)
meas = result.measurements

# depth
# print(len(cirq.Circuit(circuit.all_operations())))

# circuit svg file
# with open('circuit2.svg', 'w') as file:
#    file.write(circuit_to_svg(circuit))

# ------------------------------------------------------------------------------
# Checking Results

successTot = 0
validRuns = 0
bCount = 0

for i in range(numRuns):
    d = []
    for j in range(n):
        key = 'd{}'.format(j)
        d += list(meas[key][i])
    b = int(meas['bMeas'][i])
    roundRes = list(meas['roundRes'][i])
    if (sum(d) > 0):
        bCount += b
        x = findPreimage(roundRes, A, modulus)
        if (len(x) > 0):
            if (checkEquation(x, (x + s) % modulus, d, b, modulus) or checkEquation(x, (x - s) % modulus, d, b, modulus)):
                successTot += 1
        validRuns += 1

print(successTot / validRuns)
print(bCount / validRuns)



