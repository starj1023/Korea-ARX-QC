import math
from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control

def CHAM64_128_ENC(eng, resource):

    k0 = eng.allocate_qureg(16)
    k1 = eng.allocate_qureg(16)
    k2 = eng.allocate_qureg(16)
    k3 = eng.allocate_qureg(16)
    k4 = eng.allocate_qureg(16)
    k5 = eng.allocate_qureg(16)
    k6 = eng.allocate_qureg(16)
    k7 = eng.allocate_qureg(16)

    key_temp0  = eng.allocate_qureg(3)
    key_temp1 = eng.allocate_qureg(3)
    key_temp2 = eng.allocate_qureg(3)

    c0 = eng.allocate_qubit()
    c1 = eng.allocate_qubit()
    c2 = eng.allocate_qubit()

    x0 = eng.allocate_qureg(16)
    x1 = eng.allocate_qureg(16)
    x2 = eng.allocate_qureg(16)
    x3 = eng.allocate_qureg(16)

    if(resource == 0):
        Value_XOR(eng, k0, 0x0100, 16)
        Value_XOR(eng, k1, 0x0302, 16)
        Value_XOR(eng, k2, 0x0504, 16)
        Value_XOR(eng, k3, 0x0706, 16)
        Value_XOR(eng, k4, 0x0908, 16)
        Value_XOR(eng, k5, 0x0b0a, 16)
        Value_XOR(eng, k6, 0x0d0c, 16)
        Value_XOR(eng, k7, 0x0f0e, 16)

        Value_XOR(eng, x0, 0x1100, 16)
        Value_XOR(eng, x1, 0x3322, 16)
        Value_XOR(eng, x2, 0x5544, 16)
        Value_XOR(eng, x3, 0x7766, 16)

    i = 0

    for j in range(5):
        with Compute(eng):
            KeyGen64_type_0to7(eng, k0) #RK0
            KeyGen64_type_0to7(eng, k1) #RK1
            KeyGen64_type_0to7(eng, k2) #RK2
            KeyGen64_type_0to7(eng, k3) #RK3
            KeyGen64_type_0to7(eng, k4) #RK4
            KeyGen64_type_0to7(eng, k5) #RK5
            KeyGen64_type_0to7(eng, k6) #RK6
            KeyGen64_type_0to7(eng, k7) #RK7

        k0 = Logical_swap(eng, k0)
        k1 = Logical_swap(eng, k1)
        k2 = Logical_swap(eng, k2)
        k3 = Logical_swap(eng, k3)
        k4 = Logical_swap(eng, k4)
        k5 = Logical_swap(eng, k5)
        k6 = Logical_swap(eng, k6)
        k7 = Logical_swap(eng, k7)

        #ENC (0~7)
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k0, i, c0)
        i = i + 1
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k1, i, c1)
        i = i+1
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k2, i, c2)
        i = i + 1
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k3, i, c0)
        i = i + 1
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k4, i, c1)
        i = i + 1
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k5, i, c2)
        i = i + 1
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k6, i, c0)
        i = i + 1
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k7, i,  c1)
        i = i + 1

        k0 = logical_right(eng, k0, 1)
        k1 = logical_right(eng, k1, 1)
        k2 = logical_right(eng, k2, 1)
        k3 = logical_right(eng, k3, 1)
        k4 = logical_right(eng, k4, 1)
        k5 = logical_right(eng, k5, 1)
        k6 = logical_right(eng, k6, 1)
        k7 = logical_right(eng, k7, 1)

        Uncompute(eng)

        with Compute(eng):
            KeyGen64_type_8to15(eng, k1, key_temp0)  # RK8
            KeyGen64_type_8to15(eng, k0, key_temp1)  # RK9
            KeyGen64_type_8to15(eng, k3, key_temp2)  # RK10
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k1, i, c2)
        i = i + 1
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k0, i, c0)
        i = i + 1
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k3, i, c1)
        i = i + 1
        Uncompute(eng)

        with Compute(eng):
            KeyGen64_type_8to15(eng, k2, key_temp0)  # RK11
            KeyGen64_type_8to15(eng, k5, key_temp1)  # RK12
            KeyGen64_type_8to15(eng, k4, key_temp2)  # RK13
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k2, i, c0)
        i = i + 1
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k5, i, c1)
        i = i + 1
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k4, i, c2)
        i = i + 1
        Uncompute(eng)

        with Compute(eng):
            KeyGen64_type_8to15(eng, k7, key_temp0)  # RK14
            KeyGen64_type_8to15(eng, k6, key_temp1)  # RK15
        x0, x1, x2, x3 = Round_odd(eng, x0, x1, x2, x3, k7, i, c2)
        i = i + 1
        x0, x1, x2, x3 = Round_even(eng, x0, x1, x2, x3, k6, i, c1)
        i = i + 1
        Uncompute(eng)

    if (resource == 0):
        print_state(eng, x0, x1, x2, x3)

def CNOT16(eng, a, b):

    for i in range(16):
        CNOT | (a[i], b[i])

def Round_odd(eng, x0, x1, x2, x3, k, i, c0):
    Round_constant_XOR(eng, x0, i)
    x1 = logical_left(eng, x1, 1)
    CNOT16(eng, x1, k)
    x1 = logical_right(eng, x1, 1)
    improved_adder(eng, k, x0, c0, 15)
    #reverse
    x1 = logical_left(eng, x1, 1)
    CNOT16(eng, x1, k)
    x1 = logical_right(eng, x1, 1)
    x0 = logical_left(eng, x0, 8)

    return x1, x2, x3, x0

def Round_even(eng, x0, x1, x2, x3, k, i, c0): #real odd

    Round_constant_XOR(eng, x0, i)
    x1 = logical_left(eng, x1, 8)
    CNOT16(eng, x1, k)
    x1 = logical_right(eng, x1, 8)
    improved_adder(eng, k, x0, c0, 15)
    #reverse
    x1 = logical_left(eng, x1, 8)
    CNOT16(eng, x1, k)
    x1 = logical_right(eng, x1, 8)
    x0 = logical_left(eng, x0, 1)

    return x1, x2, x3, x0

def Logical_swap(eng, k):

    new_k = []

    new_k.append(k[15])
    for i in range(15):
        new_k.append(k[i])

    return new_k

def KeyGen64_type_0to7(eng, key): #25CNOT
    for i in range(8):
        CNOT | (key[i], key[i+8]) #

    for i in range(7):
        CNOT | (key[i+9], key[i])

    CNOT | (key[7], key[15])
    CNOT | (key[8], key[15])
    CNOT | (key[8], key[7])

    for i in range(7):
        CNOT | (key[i], key[i+8])


def KeyGen64_type_8to15(eng, key, key_temp): #35CNOT

    CNOT | (key[15], key_temp[0])
    CNOT | (key[9], key_temp[0])

    CNOT | (key[5], key_temp[1])
    CNOT | (key[15], key_temp[1])

    CNOT | (key[4], key_temp[2])
    CNOT | (key[10], key_temp[2])

    CNOT | (key[14], key[15])  # 15 14 4
    CNOT | (key[4], key[15])

    CNOT | (key[9], key[4])  # 4 9 3
    CNOT | (key[3], key[4])

    CNOT | (key[14], key[9])  # 9 14 8
    CNOT | (key[8], key[9])

    CNOT | (key[13], key[14])  # 14 13 3
    CNOT | (key[3], key[14])

    CNOT | (key[8], key[3])  # 3 8 2
    CNOT | (key[2], key[3])

    CNOT | (key[13], key[8])  # 8 13 7
    CNOT | (key[7], key[8])

    CNOT | (key[12], key[13])  # 13 12 2
    CNOT | (key[2], key[13])

    CNOT | (key[7], key[2])  # 2 7 1
    CNOT | (key[1], key[2])

    CNOT | (key[12], key[7])  # 7 6 12
    CNOT | (key[6], key[7])

    #
    CNOT | (key[11], key[12])  # 12 11 1
    CNOT | (key[1], key[12])

    CNOT | (key[6], key[1])  # 1 6 0
    CNOT | (key[0], key[1])

    CNOT | (key[11], key[6])  # 6 11 5
    CNOT | (key[5], key[6])

    CNOT | (key[10], key[11])  # 11 10 0
    CNOT | (key[0], key[11])

    CNOT | (key_temp[0], key[10])
    CNOT | (key_temp[1], key[0])
    CNOT | (key_temp[2], key[5])

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
    for i in range(16):
        a_new.append(a[(i+n)%16])

    return a_new

def logical_left(eng, a, n):

    a_new = []
    for i in range(16):
        a_new.append(a[(16-n+i)%16])

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

    for i in reversed(range(4)):
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
CHAM64_128_ENC(eng, 0)

print('Estimate cost...')
Resource = ResourceCounter()
eng = MainEngine(Resource)
CHAM64_128_ENC(eng, 1)
print(Resource)
eng.flush()
