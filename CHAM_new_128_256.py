import math
from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control

def CHAM128_256_ENC(eng, resource):

    k0 = eng.allocate_qureg(32)
    k1 = eng.allocate_qureg(32)
    k2 = eng.allocate_qureg(32)
    k3 = eng.allocate_qureg(32)
    k4 = eng.allocate_qureg(32)
    k5 = eng.allocate_qureg(32)
    k6 = eng.allocate_qureg(32)
    k7 = eng.allocate_qureg(32)

    key_temp  = eng.allocate_qureg(11)
    key_temp1 = eng.allocate_qureg(11)
    key_temp2 = eng.allocate_qureg(11)

    c0 = eng.allocate_qubit()
    c1 = eng.allocate_qubit()
    c2 = eng.allocate_qubit()

    x0 = eng.allocate_qureg(32)
    x1 = eng.allocate_qureg(32)
    x2 = eng.allocate_qureg(32)
    x3 = eng.allocate_qureg(32)

    if(resource == 0):
        Value_XOR(eng, k0, 0x03020100, 32)
        Value_XOR(eng, k1, 0x07060504, 32)
        Value_XOR(eng, k2, 0x0b0a0908, 32)
        Value_XOR(eng, k3, 0x0f0e0d0c, 32)
        Value_XOR(eng, k4, 0xf3f2f1f0, 32)
        Value_XOR(eng, k5, 0xf7f6f5f4, 32)
        Value_XOR(eng, k6, 0xfbfaf9f8, 32)
        Value_XOR(eng, k7, 0xfffefdfc, 32)

        Value_XOR(eng, x0, 0x33221100, 32)
        Value_XOR(eng, x1, 0x77665544, 32)
        Value_XOR(eng, x2, 0xbbaa9988, 32)
        Value_XOR(eng, x3, 0xffeeddcc, 32)

    i = 0

    for j in range(6):
        with Compute(eng):
            KeyGen64_type_0to7(eng, k0, key_temp) #RK0
            KeyGen64_type_0to7(eng, k1, key_temp1)  # RK1
            KeyGen64_type_0to7(eng, k2, key_temp2)  # RK2
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k0, i, c0)
        i = i + 1
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k1, i, c1)
        i = i + 1
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k2, i, c2)
        i = i + 1
        Uncompute(eng)

        with Compute(eng):
            KeyGen64_type_0to7(eng, k3, key_temp) #RK3
            KeyGen64_type_0to7(eng, k4, key_temp1)  # RK0
            KeyGen64_type_0to7(eng, k5, key_temp2)  # RK0
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k3, i, c0)
        i = i + 1
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k4, i, c1)
        i = i + 1
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k5, i, c2)
        i = i + 1
        Uncompute(eng)

        with Compute(eng):
            KeyGen64_type_0to7(eng, k6, key_temp) #RK0
            KeyGen64_type_0to7(eng, k7, key_temp1)  # RK0
            KeyGen64_type_8to15(eng, k1, key_temp2)  # RK8
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k6, i, c0)
        i = i + 1
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k7, i, c1)
        i = i + 1
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k1, i, c2)
        i = i + 1
        Uncompute(eng)

        with Compute(eng):
            KeyGen64_type_8to15(eng, k0, key_temp)  # RK9
            KeyGen64_type_8to15(eng, k3, key_temp1)  # RK10
            KeyGen64_type_8to15(eng, k2, key_temp2)  # RK11
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k0, i, c0)
        i = i + 1
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k3, i, c1)
        i = i + 1
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k2, i, c2)
        i = i + 1
        Uncompute(eng)

        with Compute(eng):
            KeyGen64_type_8to15(eng, k5, key_temp)  # RK8
            KeyGen64_type_8to15(eng, k4, key_temp1)  # RK9
            KeyGen64_type_8to15(eng, k7, key_temp2)  # RK10
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k5, i, c0)
        i = i + 1
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k4, i, c1)
        i = i + 1
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k7, i, c2)
        i = i + 1
        Uncompute(eng)

        with Compute(eng):
            KeyGen64_type_8to15(eng, k6, key_temp)  # RK11
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k6, i, c0)
        i = i + 1
        Uncompute(eng)

    if (resource == 0):
        print_state(eng, x0, x1, x2, x3)


def CNOT32(eng, a, b):

    for i in range(32):
        CNOT | (a[i], b[i])

def Round_odd(eng, x0, x1, x2, x3, k, i, c0):

    Round_constant_XOR(eng, x0, i)
    x1 = logical_left(eng, x1, 1)
    CNOT32(eng, x1, k)
    x1 = logical_right(eng, x1, 1)
    improved_adder(eng, k, x0, c0, 31)
    #reverse
    x1 = logical_left(eng, x1, 1)
    CNOT32(eng, x1, k)
    x1 = logical_right(eng, x1, 1)
    x0 = logical_left(eng, x0, 8)

    return x1, x2, x3, x0

def Round_even(eng, x0, x1, x2, x3, k, i, c0):

    Round_constant_XOR(eng, x0, i)
    x1 = logical_left(eng, x1, 8)
    CNOT32(eng, x1, k)
    x1 = logical_right(eng, x1, 8)
    improved_adder(eng, k, x0, c0, 31)
    #reverse
    x1 = logical_left(eng, x1, 8)
    CNOT32(eng, x1, k)
    x1 = logical_right(eng, x1, 8)
    x0 = logical_left(eng, x0, 1)

    return x1, x2, x3, x0

def KeyGen64_type_0to7(eng, key, key_temp): #72 CNOT

    for i in range(8):
        CNOT | (key[31 - i], key_temp[i])

    for i in range(24):
        CNOT | (key[30 - i], key[31 - i])
        CNOT | (key[23 - i], key[31 - i])

    for i in range(7):
        CNOT | (key[6 - i], key[7 - i])

    for i in range(8):
        CNOT | (key_temp[i], key[7 - i])

    CNOT | (key_temp[0], key[0])


def KeyGen64_type_8to15(eng, key, key_temp): #75CNOT

    for i in range(11):
        CNOT | (key[31 - i], key_temp[i])

    for i in range(21):
        CNOT | (key[30 - i], key[31 - i])
        CNOT | (key[20 - i], key[31 - i])

    for i in range(10):
        CNOT | (key[9 - i], key[10 - i])

    for i in range(11):
        CNOT | (key_temp[i], key[10 - i])

    CNOT | (key_temp[0], key[0])


def Round_constant_XOR(eng, k, rc):
    for i in range(7):
        if (rc >> i & 1):
            X | k[i]

def Value_XOR(eng, k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]

def logical_right(eng, a, n):
    a_new = []
    for i in range(32):
        a_new.append(a[(i+n)%32])

    return a_new

def logical_left(eng, a, n):
    a_new = []
    for i in range(32):
        a_new.append(a[(32-n+i)%32])

    return a_new

def improved_adder(eng, a, b, c, n): #n = n-1

    for i in range(n-1):
        CNOT | (a[i+1], b[i+1])

    CNOT | (a[1], c)
    Toffoli | (a[0], b[0], c)
    CNOT | (a[2], a[1])
    Toffoli | (c, b[1], a[1])
    CNOT | (a[3], a[2])

    for i in range(n-4):
        Toffoli | (a[i+1], b[i+2] ,a[i+2])
        CNOT | (a[i+4], a[i+3])
    Toffoli | (a[n-3], b[n-2], a[n-2])

    CNOT | (a[n-1], b[n])
    CNOT | (a[n], b[n])
    Toffoli | (a[n-2], b[n-1], b[n])

    for i in range(n-2):
        X | b[i+1]

    CNOT | (c, b[1])

    for i in range(n-2):
        CNOT | (a[i+1], b[i+2])

    Toffoli | (a[n-3], b[n-2], a[n-2])

    for i in range(n-4):
        Toffoli | (a[n-4-i], b[n-3-i], a[n-3-i])
        CNOT | (a[n-1-i], a[n-2-i])
        X | (b[n-2-i])
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

print('Generate Ciphertext...\n')
Simulate = ClassicalSimulator()
eng = MainEngine(Simulate)
CHAM128_256_ENC(eng, 0)

print('Estimate cost...')
Resource = ResourceCounter()
eng = MainEngine(Resource)
CHAM128_256_ENC(eng, 1)
print(Resource)
eng.flush()
