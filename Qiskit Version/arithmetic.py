import numpy as np
import matplotlib.pyplot as plt
from qiskit import(
    QuantumCircuit,
    QuantumRegister,
    ClassicalRegister,
    execute,
    Aer,
    IBMQ)
from qiskit.visualization import plot_histogram
from qiskit.tools.monitor import job_monitor

def sGate():
    circ = QuantumCircuit(1, name="S")
    circ.s(0)
    return circ

# see Qiskit documentation, implementation of QFT for n qubits
def qft(n):
    circuit = QuantumCircuit(n, name="QFT")
#    def swap_registers(circuit, n):
#        for qubit in range(n//2):
#            circuit.swap(qubit, n-qubit-1)
#        return circuit
#    def qft_rotations(circuit, n):
#        if n == 0:
#            return circuit
#        n -= 1
#        circuit.h(n)
#        for qubit in range(n):
#            circuit.cu1(np.pi/2**(n-qubit), qubit, n)
#        qft_rotations(circuit, n)
#    
#    qft_rotations(circuit, n)
#    swap_registers(circuit, n)

    # 2 qubit QFT
    circuit.h(0)
    cs = sGate().to_gate().control()
    circuit.append(cs, [1] + [0])
    circuit.h(1)
    return circuit

# inverse QFT for n qubits
def iqft(n):
    inverse_circuit = qft(n).inverse()
    return inverse_circuit

# QFT controlled on 1 qubit, applied to n
def cqft(n):
    return qft(n).to_gate().control(num_ctrl_qubits=1)

# inverse QFT controlled on 1 qubit, applied to n
def ciqft(n):
    return iqft(n).to_gate().control(num_ctrl_qubits=1)

# maps |x>|y> to |cx + y (mod modulus)>
# x, y are QuantumRegisters and c, modulus are ints
def addRegisters(x, y, c, modulus):
    n = x.size
    circuit = QuantumCircuit(x, y, name="CMULT")
    #circuit.append(qft(n), y[0:n])
    for i in range(n):
        tot = ((2 ** i) * c) % (modulus)
        for j in range(n):
            power = 2 ** ((n - 1) - j)
            if ((tot % power) != 0 or ((tot % power) == 0 and ((tot / power) % 2) != 0)):
                circuit.cu1(np.pi * tot / power, x[i], y[j])
    #circuit.append(iqft(n), y[0:n])
    return circuit

# maps |x> to |x + c (mod modulus)>
# res is QuantumRegister and c, modulus are ints
def addConstantToRegister(c, modulus, res):
    n = res.size
    circuit = QuantumCircuit(res, name="ADD")
    #circuit.append(qft(n), res[0:n])
    for i in range(n):
        power = 2 ** ((n - 1) - i)
        if ((c % power) != 0 or ((c % power) == 0 and ((c / power) % 2) != 0)):
            circuit.u1(c * np.pi / power, res[i])
    #circuit.append(iqft(n), res[0:n])
    return circuit

# controlled version of previous function (controlled on b)
# b is a single qubit QuantumRegister, res is QuantumRegister, y and modulus are ints
def controlledAddConstant(b, c, modulus, res):
    n = res.size
    circuit = QuantumCircuit(b, res, name="CADD")
    #circuit.append(cqft(n), [b[0]] + res[0:n])
    for i in range(n):
        power = 2 ** ((n - 1) - i)
        if ((c % power) != 0 or (c % power) == 0 and ((c / power) % 2) != 0):
            circuit.cu1(c * np.pi / power, b[0], res[i])
    #circuit.append(ciqft(n), [b[0]] + res[0:n])
    return circuit

# coherently adds <c, x> + b . y mod modulus to result register
# b, x, res are QuantumRegisters and y, modulus are ints, c is an array of ints
# note x is an array of QuantumRegisters
def coherentLWE(b, x, res, y, c, modulus):
    n = len(c)
    #circuit = QuantumCircuit(b, x, res, name="CoherentLWE")
    circuit = QuantumCircuit(b, name="CoherentLWE")
    for i in range(n):
        circuit.add_register(x[i])
    circuit.add_register(res)
    for j in range(n):
        if (c[j] != 0):
            circuit.append(addRegisters(x[j], res, c[j], modulus),
                           x[j][0:n] + res[0:n])
    circuit.append(controlledAddConstant(b, y, modulus, res), [b[0]] + res[0:n])
    return circuit

# inverse of previous function
def adjointCoherentLWE(b, x, res, y, c, modulus):
    inverse_circuit = coherentLWE(b, x, res, y, c, modulus).inverse()
    return inverse_circuit
