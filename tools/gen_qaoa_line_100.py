# gen_qaoa_line100.py
from qiskit import QuantumCircuit, transpile
import math

def qaoa_line(n=100, p=1, gamma=0.3, beta=0.8):
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for _ in range(p):
        # cost layer: ZZ on edges (i,i+1)
        for i in range(n-1):
            qc.cx(i, i+1)
            qc.rz(2*gamma, i+1)
            qc.cx(i, i+1)
        # mixer layer: Rx on each qubit
        for i in range(n):
            qc.rx(2*beta, i)
    return qc

if __name__ == "__main__":
    circ = qaoa_line(n=100, p=1, gamma=0.3, beta=0.8)
    # transpile to CZ/U3 to match ZAC’s basis
    circ = transpile(circ, basis_gates=["cz", "u3"])
    circ.qasm(filename="benchmark/hpca/qaoa_line_100.qasm")
    print("Wrote benchmark/hpca/qaoa_line_100.qasm")
