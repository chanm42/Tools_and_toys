import numpy as np
##
##
## Tensor ordering is left-to-right in this file, qubit 1 is the leftmost position in the ket
##
##

## Set debug_flag = True for debugging
debug_flag = False

##=====================================================================================##
def pretty_print(M):
##=====================================================================================##
    nrow,ncol = M.shape
    for row in range(nrow):
        row_str = []
        for col in range(ncol):
            x            = M[row,col]
            strx         = str(x)
            row_str.append(strx)
            row_out      = "\t".join(row_str)
        print(row_out)
    return None

def x_gate():
    return np.array([[0, 1], [1, 0]],dtype=complex)

def y_gate():
    return np.array([[0, -1j], [1j, 0]],dtype=complex)

def z_gate():
    return np.array([[1, 0], [0, -1]],dtype=complex)

def hadamard_gate():
    return (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]],dtype=complex)

def phase_gate():
    return np.array([[1, 0], [0, 1j]],dtype=complex)

def t_gate():
    return np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]],dtype=complex)

def cnot_gate():
    return np.array([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 0, 1],
                     [0, 0, 1, 0]],dtype=complex)

def swap_gate():
    return np.array([[1, 0, 0, 0],
                     [0, 0, 1, 0],
                     [0, 1, 0, 0],
                     [0, 0, 0, 1]],dtype=complex)

def controlled_phase_gate():
    return np.array([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1j]],dtype=complex)

def rotation_x(theta):
    return np.array(
        [[np.cos(theta/2, -1j * np.sin(theta/2))],
         [-1j*np.sin(theta/2), np.cos(theta/2)]],dtype=complex
    )

def rotation_y(theta):
    return np.array(
        [[np.cos(theta/2), -np.sin(theta/2)],
         [np.sin(theta/2), np.cos(theta/2)]],dtype=complex
    )

def rotation_z(theta):
    return np.array(
        [[np.exp(-1j*(theta/2)),0.0],
         [0.0, np.exp(-1j*(theta/2))]],dtype=complex
    )

##
##
## Define fermionic operators as Pauli strings
##
##

def creation(j,nqubit):
    a_dag = 1
    I     = np.eye(2,dtype=complex)
    for i in range(nqubit):
        if i < j:
            a_dag = np.kron(a_dag,z_gate())
        elif i == j:
            a_dag = np.kron(a_dag,((0.5)*(x_gate() - (1j * y_gate()))))
        else:
            a_dag = np.kron(a_dag,I)

    return a_dag

def annihilation(j,nqubit):
    a = 1
    I = np.eye(2,dtype=complex)
    for i in range(nqubit):
        if i < j:
            a = np.kron(a,z_gate())
        elif i == j:
            a = np.kron(a,((0.5)*(x_gate() + (1j * y_gate()))))
        else:
            a = np.kron(a,I)

    return a

def number(j,nqubit):
    num = creation(j,nqubit) @ annihilation(j,nqubit)
    return num

def commutator(A,B):
    assert type(A) == type(B), "A and B must both be arrays"
    assert type(A) == np.ndarray
    comm = A@B - B@A
    return comm

def anticommutator(A,B):
    assert type(A) == type(B), "A and B must both be arrays"
    assert type(A) == np.ndarray
    anticom = A@B + B@A
    return anticom

##
##
## Define some useful operations
##
##

def apply_single_qubit_gate(nqubit, qubit, gate):
    assert qubit < nqubit, "Qubit index out of range"
    op = 1
    I  = np.eye(2,dtype=complex)
    for i in range(nqubit):
        if i == qubit:
            op = np.kron(op,gate)
        else:
            op = np.kron(op,I)
    return op

def apply_two_qubit_gate(nqubit, qubit1, qubit2, gate):
    assert qubit1 < nqubit, "Qubit1 index out of range"
    assert qubit2 < nqubit, "Qubit2 index out of range"
    assert qubit1 != qubit2, "Qubit1 and Qubit2 must be different"
    
    if qubit2 not in [qubit1-1, qubit1+1]:
        raise ValueError("Qubit1 and Qubit2 must be adjacent")
    
    op = 1
    q  = 0
    I  = np.eye(2,dtype=complex)
    while q < nqubit:
        if q == qubit1:
            op = np.kron(op,gate)
            q  += 2
        else:
            op = np.kron(op,I)
            q += 1
    return op


##=====================================================================================##
def test1A():
##=====================================================================================##
    x = x_gate()
    y = y_gate()
    z = z_gate()
    h = hadamard_gate()
    p = phase_gate()
    t = t_gate()
    c = cnot_gate()
    s = swap_gate()
    cp = controlled_phase_gate()
    a_dag = creation(2,4)
    a = annihilation(2,4)
    num = number(2,4)
    pretty_print(x)
    print()
    pretty_print(y)
    print()
    pretty_print(z)
    print()
    pretty_print(h)
    print()
    pretty_print(p)     
    print()
    pretty_print(t)
    print()
    pretty_print(c)
    print()
    pretty_print(s)
    print()
    pretty_print(cp)
    print()
    pretty_print(a_dag)
    print()
    pretty_print(a)
    print()
    pretty_print(num)
    print()
    print(">> Checking commutation")
    print("[a,a_dag] = 1 - 2n")
    comm  = commutator(a, a_dag)
    exact = 1 - 2*num
    pretty_print(comm)
    print(">>Checking anticommutation:")
    print("{a,a_dag} = delta_ij")
    anticom = anticommutator(a, a_dag)
    pretty_print(anticom)
    print()
    print(">>Testing off-diagonal entries:")
    for i in range(4):
        for j in range(4):
            comm_ij = commutator(annihilation(i,4), creation(j,4))
            anticom_ij = anticommutator(annihilation(i,4), creation(j,4))
            print(f"i={i}, j={j} => [a_{i}, a_{j}^dag] =")
            pretty_print(comm_ij)
            print(f"i={i}, j={j} => {{a_{i}, a_{j}^dag}} =")
            pretty_print(anticom_ij)
            print()

def main():
    if debug_flag:
        test1A()


main()