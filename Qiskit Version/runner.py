from arithmetic import *
from typing import *
import json

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

# og
#A = np.array([[0, 2, 2], [0, 1, 1], [1, 1, 0], [1, 2, 2], [2, 3, 1], [1, 2, 0]])
A = np.array([[0, 2, 0], [3, 2, 2], [0, 3, 1], [2, 0, 2], [1, 1, 1], [0, 2, 2], [3, 3, 0]])
A = A.tolist()
s = np.array([[1], [0], [1]])
e = np.array([[1], [0], [0], [0], [0], [1], [0]])
y = LWEfunEval(A, s, e, modulus)
y = (list(map(lambda w: w[0], y.tolist())))

nQubits = bits * n + bits + m

# initialize registers
b = QuantumRegister(1, 'b')
x = []
for i in range(n):
    x.append(QuantumRegister(bits, 'x{}'.format(i)))
    
anc = QuantumRegister(bits, 'anc') # ancilla
r = QuantumRegister(m - 1, 'r') # rounded, m - 1 now
bMeas = ClassicalRegister(1, 'bM')
xMeas = [] # store result of measuring x
for i in range(n):
    xMeas.append(ClassicalRegister(bits, 'xM{}'.format(i)))
z = ClassicalRegister(m, 'z') # store result of measuring r register

# put registers in circuit
circ = QuantumCircuit(b)
for i in range(n):
    circ.add_register(x[i])
circ.add_register(anc, r, bMeas)
for i in range(n):
    circ.add_register(xMeas[i])
circ.add_register(z)

# prepare superpositions over b and x
# have \sum_{b, x} |b>|x>|0>
circ.h(b[0])
for i in range(n):
    for j in range(bits):
        circ.h(x[i][j])

# put ancilla in phase
for i in range(bits):
    circ.h(anc[i])

# makes list of qubits to apply coherentLWE to
qubits = [b[0]]
for i in range(n):
    for j in range(bits):
        qubits.append(x[i][j])
qubits += anc[0:bits]

# compute rounded result
# use ancilla to store pre-rounded result; clear ancilla each time
# final state should be |b>|x>|0>|round(Ax + b . (As + e'))>
for i in range(m - 1):
    circ.append(coherentLWE(b, x, anc, y[i], A[i], modulus), qubits)

    circ.append(iqft(bits), anc[0:bits])

    circ.cx(anc[0], r[i])

    circ.append(qft(bits), anc[0:bits])

    circ.append(adjointCoherentLWE(b, x, anc, y[i], A[i], modulus), qubits)

# last iteration MSB of ancilla is rounded result
circ.append(coherentLWE(b, x, anc, y[m - 1], A[m - 1], modulus), qubits)

circ.append(iqft(bits), anc[0:bits])

# get an equation in s
circ.h(b[0])
for i in range(n):
    for j in range(bits):
        circ.h(x[i][j])

circ.measure(b, bMeas)
for i in range(n):
    circ.measure(x[i], xMeas[i])

for i in range(m - 1):
    circ.measure(r[i], z[i])
circ.measure(anc[0], z[m - 1])

#circ.draw(output='latex', filename='circ.png')

#decomposed_circ = circ.decompose()
#decomposed_circ.draw(output='latex', filename='decompose_circ.png')

# print(circ.decompose().decompose().depth())
# print(circ.decompose().decompose().decompose().depth())

# backend = Aer.get_backend('qasm_simulator')
provider = IBMQ.load_account()
# backend = provider.get_backend('ibmq_qasm_simulator')
backend = provider.get_backend('ibmq_16_melbourne')
job = execute(circ, backend=backend, shots=8192)
job_monitor(job)
result = job.result()
counts = result.get_counts(circ)

with open('ibmq_16_melbourne_m=7_counts', 'w') as dataFile:
        json.dump(counts, dataFile)

# print(counts)

# plot_histogram(counts)
