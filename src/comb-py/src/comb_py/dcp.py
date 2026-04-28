from pytket.circuit import Circuit
from pauli_hugr.tables import compute_cnot_plus_t_circuit_table

# These two imports below are only needed for the type annotation.
from pauli_hugr.errors import DCPOperator
from pauli_hugr.faults import Fault


circ = Circuit(2).CX(0, 1).T(1).CX(0, 1)
circ.name = "CNOT_plus_T"

faults_and_output_errors: dict[tuple[Fault, ...], DCPOperator] = (
    compute_cnot_plus_t_circuit_table(circ, "test", all_paulis=True)
)

# Interpretation:
# Generate Pauli faults and their corresponding output errors
# 15 Pauli faults for CX gates
# 3 Faults on the output wire of the T gate
# 2 X faults on |0> initialization.
assert len(faults_and_output_errors) == 15 + 3 + 15 + 2

print(faults_and_output_errors)
