import math

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control


def LEA192_Enc(eng):
    # Qubit
    k0 = eng.allocate_qureg(32)  # Key
    k1 = eng.allocate_qureg(32)
    k2 = eng.allocate_qureg(32)
    k3 = eng.allocate_qureg(32)
    k4 = eng.allocate_qureg(32)
    k5 = eng.allocate_qureg(32)

    c0 = eng.allocate_qubit()  # carry qubit
    c1 = eng.allocate_qubit()  # carry qubit
    c2 = eng.allocate_qubit()  # carry qubit
    c3 = eng.allocate_qubit()  # carry qubit
    c4 = eng.allocate_qubit()  # carry qubit
    c5 = eng.allocate_qubit()  # carry qubit

    x0 = eng.allocate_qureg(32)  # plaintext --> ciphertext
    x1 = eng.allocate_qureg(32)
    x2 = eng.allocate_qureg(32)
    x3 = eng.allocate_qureg(32)

    theta0 = eng.allocate_qureg(32)  # theta
    theta1 = eng.allocate_qureg(32)  # theta
    theta2 = eng.allocate_qureg(32)  # theta
    theta3 = eng.allocate_qureg(32)  # theta
    theta4 = eng.allocate_qureg(32)  # theta
    theta5 = eng.allocate_qureg(32)  # theta

    Round_constant_XOR(eng, theta0, 0xc3efe9db, 32)
    Round_constant_XOR(eng, theta1, 0xc3efe9db, 32)
    Round_constant_XOR(eng, theta2, 0xc3efe9db, 32)
    Round_constant_XOR(eng, theta3, 0xc3efe9db, 32)
    Round_constant_XOR(eng, theta4, 0xc3efe9db, 32)
    Round_constant_XOR(eng, theta5, 0xc3efe9db, 32)


    if (resource == 0):
        #key
        Round_constant_XOR(eng, k0, 0x3c2d1e0f, 32)
        Round_constant_XOR(eng, k1, 0x78695a4b, 32)
        Round_constant_XOR(eng, k2, 0xb4a59687, 32)
        Round_constant_XOR(eng, k3, 0xf0e1d2c3, 32)
        Round_constant_XOR(eng, k4, 0xc3d2e1f0, 32)
        Round_constant_XOR(eng, k5, 0x8796a5b4, 32)

        # Plaintext
        Round_constant_XOR(eng, x0, 0x23222120, 32)
        Round_constant_XOR(eng, x1, 0x27262524, 32)
        Round_constant_XOR(eng, x2, 0x2b2a2928, 32)
        Round_constant_XOR(eng, x3, 0x2f2e2d2c, 32)

    i = 0

    k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
    i = i + 1
    x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

    change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0xc3efe9db, 0x44626b02)
    k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
    i = i + 1
    x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

    change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0x44626b02, 0x79e27c8a)
    k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
    i = i + 1
    x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

    change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0x79e27c8a, 0x78df30ec)
    k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
    i = i + 1
    x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

    change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0x78df30ec, 0x715ea49e)
    k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
    i = i + 1
    x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

    change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0x715ea49e, 0xc785da0a)
    k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
    i = i + 1
    x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

    ##
    for j in range(3):
        change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0xc785da0a, 0xc3efe9db)
        k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
        i = i + 1
        x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

        change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0xc3efe9db, 0x44626b02)
        k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
        i = i + 1
        x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

        change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0x44626b02, 0x79e27c8a)
        k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
        i = i + 1
        x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

        change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0x79e27c8a, 0x78df30ec)
        k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
        i = i + 1
        x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

        change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0x78df30ec, 0x715ea49e)
        k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
        i = i + 1
        x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

        change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0x715ea49e, 0xc785da0a)
        k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
        i = i + 1
        x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

    change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0xc785da0a, 0xc3efe9db)
    k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
    i = i + 1
    x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

    change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0xc3efe9db, 0x44626b02)
    k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
    i = i + 1
    x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

    change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0x44626b02, 0x79e27c8a)
    k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
    i = i + 1
    x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)

    change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, 0x79e27c8a, 0x78df30ec)
    k0, k1, k2, k3, k4, k5 = KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i)
    i = i + 1
    x0, x1, x2, x3 = Enc(eng, x0, x1, x2, x3, k0, k1, k2, k3, k4, k5, c0)
    ##

    if (resource == 0):
        print_state(eng, x0, x1, x2, x3)


def Enc(eng, x0, x1, x2, x3, rk0, rk1, rk2, rk3, rk4, rk5, c0):
    CNOT32(eng, rk5, x3)
    CNOT32(eng, rk4, x2)
    improved_adder(eng, x2, x3, c0, 31)
    CNOT32(eng, rk4, x2)  # reverse
    x3 = logical_right(eng, x3, 3)

    CNOT32(eng, rk3, x2)
    CNOT32(eng, rk2, x1)
    improved_adder(eng, x1, x2, c0, 31)
    CNOT32(eng, rk2, x1)  # reverse
    x2 = logical_right(eng, x2, 5)

    CNOT32(eng, rk1, x1)
    CNOT32(eng, rk0, x0)
    improved_adder(eng, x0, x1, c0, 31)
    CNOT32(eng, rk0, x0)  # reverse
    x1 = logical_left(eng, x1, 9)

    return x1, x2, x3, x0


def KeySchedule(eng, k0, k1, k2, k3, k4, k5, theta0, theta1, theta2, theta3, theta4, theta5, c0, c1, c2, c3, c4, c5, i):
    theta0 = logical_left(eng, theta0, i)
    improved_adder(eng, theta0, k0, c0, 31)

    theta1 = logical_left(eng, theta1, i+1)
    improved_adder(eng, theta1, k1, c1, 31)

    theta2 = logical_left(eng, theta2, i+2)
    improved_adder(eng, theta2, k2, c2, 31)

    theta3 = logical_left(eng, theta3, i+3)
    improved_adder(eng, theta3, k3, c3, 31)

    theta4 = logical_left(eng, theta4, i+4)
    improved_adder(eng, theta4, k4, c4, 31)

    theta5 = logical_left(eng, theta5, i+5)
    improved_adder(eng, theta5, k5, c5, 31)

    # reverse

    k0 = logical_left(eng, k0, 1)
    k1 = logical_left(eng, k1, 3)
    k2 = logical_left(eng, k2, 6)
    k3 = logical_left(eng, k3, 11)
    k4 = logical_left(eng, k4, 13)
    k5 = logical_left(eng, k5, 17)

    return k0, k1, k2, k3, k4, k5


def Round_constant_XOR(eng, k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]


def change_theta(eng, theta0, theta1, theta2, theta3, theta4, theta5, a, b):
    a = a ^ b
    Round_constant_XOR(eng, theta0, a, 32)
    Round_constant_XOR(eng, theta1, a, 32)
    Round_constant_XOR(eng, theta2, a, 32)
    Round_constant_XOR(eng, theta3, a, 32)
    Round_constant_XOR(eng, theta4, a, 32)
    Round_constant_XOR(eng, theta5, a, 32)


def logical_right(eng, a, n):
    a_new = []
    for i in range(32):
        a_new.append(a[(i + n) % 32])

    return a_new


def logical_left(eng, a, n):
    a_new = []
    for i in range(32):
        a_new.append(a[(32 - n + i) % 32])

    return a_new


def CNOT32(eng, a, b):
    for i in range(32):
        CNOT | (a[i], b[i])


def improved_adder(eng, a, b, c, n):  # n = n-1

    for i in range(n - 1):
        CNOT | (a[i + 1], b[i + 1])

    CNOT | (a[1], c)
    Toffoli | (a[0], b[0], c)
    CNOT | (a[2], a[1])
    Toffoli | (c, b[1], a[1])
    CNOT | (a[3], a[2])

    for i in range(n - 4):
        Toffoli | (a[i + 1], b[i + 2], a[i + 2])
        CNOT | (a[i + 4], a[i + 3])
    Toffoli | (a[n - 3], b[n - 2], a[n - 2])

    CNOT | (a[n - 1], b[n])
    CNOT | (a[n], b[n])
    Toffoli | (a[n - 2], b[n - 1], b[n])

    for i in range(n - 2):
        X | b[i + 1]

    CNOT | (c, b[1])

    for i in range(n - 2):
        CNOT | (a[i + 1], b[i + 2])

    Toffoli | (a[n - 3], b[n - 2], a[n - 2])

    for i in range(n - 4):
        Toffoli | (a[n - 4 - i], b[n - 3 - i], a[n - 3 - i])
        CNOT | (a[n - 1 - i], a[n - 2 - i])
        X | (b[n - 2 - i])
    Toffoli | (c, b[1], a[1])
    CNOT | (a[3], a[2])
    X | b[2]
    Toffoli | (a[0], b[0], c)
    CNOT | (a[2], a[1])
    X | b[1]
    CNOT | (a[1], c)

    for i in range(n):
        CNOT | (a[i], b[i])

def print_state(eng, x0, x1, x2, x3):

    All(Measure) | x0
    All(Measure) | x1
    All(Measure) | x2
    All(Measure) | x3

    print('Ciphertext : 0x', end='')
    print_hex(eng, x0)
    print_hex(eng, x1)
    print_hex(eng, x2)
    print_hex(eng, x3)
    print('\n')

def print_hex(eng, qubits):

    for i in reversed(range(8)):
        temp = 0
        temp = temp+int(qubits[4*i+3])*8
        temp = temp+int(qubits[4*i+2])*4
        temp = temp+int(qubits[4*i+1])*2
        temp = temp+int(qubits[4*i])

        temp = hex(temp)
        y = temp.replace("0x", "")
        print(y, end='')


global resource

print('Generate Ciphertext...\n')
resource = 0
Simulate = ClassicalSimulator()
eng = MainEngine(Simulate)
LEA192_Enc(eng)

print('Estimate cost...')
resource = 1
Resource = ResourceCounter()
eng = MainEngine(Resource)
LEA192_Enc(eng)
print(Resource)

eng.flush()

