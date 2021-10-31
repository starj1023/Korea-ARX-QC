import math

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, Swap
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control


######################  Encryption ######################

def HIGHT_Enc(eng):
    # key
    key = eng.allocate_qureg(128)
    plaintext = eng.allocate_qureg(64)

    s0 = eng.allocate_qureg(8)
    s1 = eng.allocate_qureg(8)
    s2 = eng.allocate_qureg(8)
    s3 = eng.allocate_qureg(8)

    c0 = eng.allocate_qubit()  # carry qubit
    c1 = eng.allocate_qubit()
    c2 = eng.allocate_qubit()
    c3 = eng.allocate_qubit()
    c4 = eng.allocate_qubit()
    c5 = eng.allocate_qubit()
    c6 = eng.allocate_qubit()
    c7 = eng.allocate_qubit()


    if(resource == 0):
        Round_constant_XOR(eng, key, 0xffeeddccbbaa99887766554433221100, 128)
        Round_constant_XOR(eng, plaintext,  0x0011223344556677, 64)

    k0 = []; k1 = []; k2 = []; k3 = []; k4 = []; k5 = []; k6 = []; k7 = []
    k8 = []; k9 = []; k10 = []; k11 = []; k12 = []; k13 = []; k14 = []; k15 = []
    p0 = []; p1 = [];  p2 = []; p3 = []; p4 = []; p5 = []; p6 = []; p7 = []

    for i in range(8):
        p0.append(plaintext[i])
        p1.append(plaintext[8+i])
        p2.append(plaintext[16+i])
        p3.append(plaintext[24+i])
        p4.append(plaintext[32+i])
        p5.append(plaintext[40+i])
        p6.append(plaintext[48+i])
        p7.append(plaintext[56+i])

    for i in range(8):
        k0.append(key[i])
        k1.append(key[8+i])
        k2.append(key[16+i])
        k3.append(key[24+i])
        k4.append(key[32+i])
        k5.append(key[40+i])
        k6.append(key[48+i])
        k7.append(key[56+i])
        k8.append(key[64+i])
        k9.append(key[72+i])
        k10.append(key[80+i])
        k11.append(key[88+i])
        k12.append(key[96+i])
        k13.append(key[104+i])
        k14.append(key[112+i])
        k15.append(key[120+i])

    X | s0[1]
    X | s0[3]
    X | s0[4]
    X | s0[6]

    X | s1[1]
    X | s1[3]
    X | s1[4]
    X | s1[6]

    X | s2[1]
    X | s2[3]
    X | s2[4]
    X | s2[6]

    X | s3[1]
    X | s3[3]
    X | s3[4]
    X | s3[6]

    # first whitening
    # WK0 = mk12
    # WK1 = mk13
    # WK2 = mk14
    # WK3 = mk15

    improved_adder(eng, k12, p0, c0, 7)
    CNOT8(eng, k13, p2)
    improved_adder(eng, k14, p4, c1, 7)
    CNOT8(eng, k15, p6)

    #1 ~ 31 round
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round_new(eng, s0, s1, s2, s3, k0, k1, k2, k3, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)

    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k4, k5, k6, k7, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k8, k9, k10, k11, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k12, k13, k14, k15, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)

    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k7, k0, k1, k2, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k3, k4, k5, k6, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k15, k8, k9, k10, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k11, k12, k13, k14, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)

    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k6, k7, k0, k1, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k2, k3, k4, k5, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k14, k15, k8, k9, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k10, k11, k12, k13, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)

    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k5, k6, k7, k0, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k1, k2, k3, k4, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k13, k14, k15, k8, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k9, k10, k11, k12, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)

    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k4, k5, k6, k7, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k0, k1, k2, k3, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k12, k13, k14, k15, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k8, k9, k10, k11, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)

    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k3, k4, k5, k6, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k7, k0, k1, k2, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k11, k12, k13, k14, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k15, k8, k9, k10, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)

    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k2, k3, k4, k5, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k6, k7, k0, k1, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k10, k11, k12, k13, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k14, k15, k8, k9, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)

    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k1, k2, k3, k4, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k5, k6, k7, k0, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    p0, p1, p2, p3, p4, p5, p6, p7, s0, s1, s2, s3 = Round(eng, s0, s1, s2, s3, k9, k10, k11, k12, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)
    Last_Round(eng, s0, s1, s2, s3, k13, k14, k15, k8, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3)

    improved_adder(eng, k0, p0, c0, 7)
    CNOT8(eng, k1, p2)
    improved_adder(eng, k2, p4, c1, 7)
    CNOT8(eng, k3, p6)

    if (resource == 0):
        print_state(eng, p0, p1, p2, p3, p4, p5, p6, p7)




def Round_new(eng, s0, s1, s2, s3, rk0, rk1, rk2, rk3, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3):
    with Compute(eng):
        F1(eng, p0)
        improved_adder(eng, s0, rk0, c0, 7)
        CNOT8(eng, rk0, p0)
    improved_adder(eng, p0, p1, c0, 7) #p2

    s1 = generate_s(eng, s1)

    with Compute(eng):
        F0(eng, p2)
        improved_adder(eng, s1, rk1, c1, 7)
        improved_adder(eng, rk1, p2, c1, 7)
    CNOT8(eng, p2, p3) #p4


    s2 = generate_s(eng, s2)
    s2 = generate_s(eng, s2)

    #####################

    with Compute(eng):
        F1(eng, p4)
        improved_adder(eng, s2, rk2, c2, 7)
        CNOT8(eng, rk2, p4)
    improved_adder(eng, p4, p5, c2, 7) #p6


    s3 = generate_s(eng, s3)
    s3 = generate_s(eng, s3)
    s3 = generate_s(eng, s3)

    with Compute(eng):
        F0(eng, p6)
        improved_adder(eng, s3, rk3, c3, 7)
        improved_adder(eng, rk3, p6, c3, 7)
    CNOT8(eng, p6, p7) # p0

    Uncompute(eng)
    Uncompute(eng)
    Uncompute(eng)
    Uncompute(eng)


    new_p0 = []; new_p1 = []; new_p2 = []; new_p3 = []; new_p4 = []; new_p5 = []; new_p6 = []; new_p7 = []

    for i in range(8):
        new_p0.append(p7[i])
        new_p1.append(p0[i])
        new_p2.append(p1[i])
        new_p3.append(p2[i])
        new_p4.append(p3[i])
        new_p5.append(p4[i])
        new_p6.append(p5[i])
        new_p7.append(p6[i])

    return new_p0, new_p1, new_p2, new_p3, new_p4, new_p5, new_p6, new_p7, s0, s1, s2, s3

def Round(eng, s0, s1, s2, s3, rk0, rk1, rk2, rk3, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3):

    s0 = generate_s(eng, s0)
    s0 = generate_s(eng, s0)
    s0 = generate_s(eng, s0)
    s0 = generate_s(eng, s0)
    with Compute(eng):
        F1(eng, p0)
        improved_adder(eng, s0, rk0, c0, 7)
        CNOT8(eng, rk0, p0)
    improved_adder(eng, p0, p1, c0, 7) #p2

    s1 = generate_s(eng, s1)
    s1 = generate_s(eng, s1)
    s1 = generate_s(eng, s1)
    s1 = generate_s(eng, s1)
    with Compute(eng):
        improved_adder(eng, s1, rk1, c1, 7)
        F0(eng, p2)
        improved_adder(eng, rk1, p2, c1, 7)
    CNOT8(eng, p2, p3) #p4

    #####################
    s2 = generate_s(eng, s2)
    s2 = generate_s(eng, s2)
    s2 = generate_s(eng, s2)
    s2 = generate_s(eng, s2)
    with Compute(eng):
        improved_adder(eng, s2, rk2, c2, 7)
        F1(eng, p4)
        CNOT8(eng, rk2, p4)
    improved_adder(eng, p4, p5, c2, 7) #p6

    s3 = generate_s(eng, s3)
    s3 = generate_s(eng, s3)
    s3 = generate_s(eng, s3)
    s3 = generate_s(eng, s3)
    with Compute(eng):
        improved_adder(eng, s3, rk3, c3, 7)
        F0(eng, p6)
        improved_adder(eng, rk3, p6, c3, 7)
    CNOT8(eng, p6, p7) # p0

    Uncompute(eng)
    Uncompute(eng)
    Uncompute(eng)
    Uncompute(eng)

    new_p0 = []; new_p1 = []; new_p2 = []; new_p3 = []; new_p4 = []; new_p5 = []; new_p6 = []; new_p7 = []

    for i in range(8):
        new_p0.append(p7[i])
        new_p1.append(p0[i])
        new_p2.append(p1[i])
        new_p3.append(p2[i])
        new_p4.append(p3[i])
        new_p5.append(p4[i])
        new_p6.append(p5[i])
        new_p7.append(p6[i])

    return new_p0, new_p1, new_p2, new_p3, new_p4, new_p5, new_p6, new_p7, s0, s1, s2, s3

def Last_Round(eng, s0, s1, s2, s3, rk0, rk1, rk2, rk3, p0, p1, p2, p3, p4, p5, p6, p7, c0, c1, c2, c3):

    s0 = generate_s(eng, s0)
    s0 = generate_s(eng, s0)
    s0 = generate_s(eng, s0)
    s0 = generate_s(eng, s0)
    improved_adder(eng, s0, rk0, c0, 7)
    with Compute(eng):
        F1(eng, p0)
        CNOT8(eng, rk0, p0)
    improved_adder(eng, p0, p1, c0, 7) #p2

    s1 = generate_s(eng, s1)
    s1 = generate_s(eng, s1)
    s1 = generate_s(eng, s1)
    s1 = generate_s(eng, s1)
    improved_adder(eng, s1, rk1, c1, 7)
    with Compute(eng):
        F0(eng, p2)
        improved_adder(eng, rk1, p2, c1, 7)
    CNOT8(eng, p2, p3) #p4

    #####################
    s2 = generate_s(eng, s2)
    s2 = generate_s(eng, s2)
    s2 = generate_s(eng, s2)
    s2 = generate_s(eng, s2)

    improved_adder(eng, s2, rk2, c2, 7)
    with Compute(eng):
        F1(eng, p4)
        CNOT8(eng, rk2, p4)
    improved_adder(eng, p4, p5, c2, 7) #p6

    s3 = generate_s(eng, s3)
    s3 = generate_s(eng, s3)
    s3 = generate_s(eng, s3)
    s3 = generate_s(eng, s3)

    improved_adder(eng, s3, rk3, c3, 7)
    with Compute(eng):
        F0(eng, p6)
        improved_adder(eng, rk3, p6, c3, 7)
    CNOT8(eng, p6, p7) # p0

    Uncompute(eng)
    Uncompute(eng)
    Uncompute(eng)
    Uncompute(eng)

def generate_s(eng, s):
    CNOT | (s[3], s[0])
    new_s = []
    for i in range(7):
        new_s.append(s[(1+i)%7])
    new_s.append(s[7])

    return new_s

def F0(eng, x):  # 21cnot

    CNOT | (x[4], x[5])  # 5 = 5+4
    CNOT | (x[3], x[4])  # 4 = 4+3
    CNOT | (x[2], x[3])  # 3 = 3+2
    CNOT | (x[2], x[0])  # 0 = 0+2
    CNOT | (x[4], x[2])  # 2 = 2+4+3
    CNOT | (x[5], x[2])  # 2 = 2+5+3
    CNOT | (x[7], x[5])  # 5 = 5+4+7
    CNOT | (x[6], x[4])  # 4 = 4+3+6
    CNOT | (x[0], x[3])  # 3 = 3+0
    CNOT | (x[1], x[3])  # 3 = 3+0+1
    CNOT | (x[6], x[1])  # 1 = 1+6
    CNOT | (x[7], x[1])  # 1 = 1+6+7
    CNOT | (x[7], x[0])  # 0 = 0+2+7
    CNOT | (x[5], x[6])  # 6 = 6+5+4+7
    CNOT | (x[4], x[6])  # 6 = 3+5+7
    CNOT | (x[3], x[6])  # 6 = 0+1+5+7
    CNOT | (x[1], x[6])  # 6 = 0+5+6
    CNOT | (x[6], x[7])  # 7 = 7+0+5+6
    CNOT | (x[5], x[7])  # 7 = 4+0+6
    CNOT | (x[1], x[7])  # 7 = 1+4+0+7
    CNOT | (x[0], x[7])  # 7 = 1+4+2
    if (resource == 0):
        Swap | (x[0], x[1])
        Swap | (x[3], x[2])
        Swap | (x[3], x[7])
        Swap | (x[7], x[4])
        Swap | (x[7], x[5])
        Swap | (x[7], x[6])

def F1(eng, x):  # 24cnot
    CNOT | (x[3], x[4])
    CNOT | (x[1], x[4])  # 4= 4+3+1
    CNOT | (x[2], x[3])
    CNOT | (x[0], x[3])  # 3= 3+2+0
    CNOT | (x[1], x[2])
    CNOT | (x[7], x[2])  # 2 = 2+1+7
    CNOT | (x[0], x[1])
    CNOT | (x[6], x[1])  # 1 = 1+0+6
    CNOT | (x[7], x[0])  # 0 = 0+7+5
    CNOT | (x[5], x[0])
    CNOT | (x[6], x[7])  # 7 = 7+6+4
    CNOT | (x[4], x[7])
    CNOT | (x[3], x[7])
    CNOT | (x[2], x[7])
    CNOT | (x[0], x[7])
    CNOT | (x[5], x[7])
    CNOT | (x[0], x[6])
    CNOT | (x[7], x[6])  # 4 5 6 0  #6 = 0+5+4
    CNOT | (x[4], x[6])
    CNOT | (x[1], x[6])  # 4 3 6 0  #6 = 5+3+6
    CNOT | (x[7], x[5])
    CNOT | (x[6], x[5])  # 3 5 7 4  #5 = 7+4+3
    CNOT | (x[3], x[5])
    CNOT | (x[0], x[5])  # 3 5 7 2  #5 = 4+2+5
    if(resource == 0):
        Swap | (x[0], x[5])
        Swap | (x[1], x[6])
        Swap | (x[2], x[7])
        Swap | (x[3], x[5])
        Swap | (x[4], x[6])
        Swap | (x[5], x[7])
        Swap | (x[6], x[7])

def Round_constant_XOR(eng, k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]

def CNOT8(eng, a, b):
    for i in range(8):
        CNOT | (a[i], b[i])

def new_adder(eng, a, b, len):
    for i in range(len-1):
        CNOT | (a[i+1], b[i+1])

    for i in range(len-2):
        CNOT | (a[len-2-i], a[len-1-i])

    for i in range(len-1):
        Toffoli | (b[i], a[i], a[i+1])

    for i in range(len-1):
        CNOT | (a[len-1-i], b[len-1-i])
        Toffoli | (b[len-2-i], a[len-2-i], a[len-1-i])

    for i in range(len-2):
        CNOT | (a[i+1], a[i+2])

    for i in range(len):
        CNOT | (a[i], b[i])

def improved_adder(eng, a, b, c, n): #n = n-1

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

def print_state(eng, p0, p1, p2, p3, p4, p5, p6, p7):

    All(Measure) | p0; All(Measure) | p1; All(Measure) | p2; All(Measure) | p3
    All(Measure) | p4; All(Measure) | p5; All(Measure) | p6; All(Measure) | p7
    print('Ciphertext : 0x', end='')
    print_hex(eng, p7)
    print_hex(eng, p6)
    print_hex(eng, p5)
    print_hex(eng, p4)
    print_hex(eng, p3)
    print_hex(eng, p2)
    print_hex(eng, p1)
    print_hex(eng, p0)
    print('\n')

def print_hex(eng, qubits):

    for i in reversed(range(2)):
        temp = 0
        temp = temp + int(qubits[4 * i + 3]) * 8
        temp = temp + int(qubits[4 * i + 2]) * 4
        temp = temp + int(qubits[4 * i + 1]) * 2
        temp = temp + int(qubits[4 * i])

        temp = hex(temp)
        y = temp.replace("0x", "")
        print(y, end='')

global resource

print('Generate Ciphertext...\n')
resource = 0
Simulate = ClassicalSimulator()
eng = MainEngine(Simulate)
HIGHT_Enc(eng)

print('Estimate cost...')
resource = 1
Resource = ResourceCounter()
eng = MainEngine(Resource)
HIGHT_Enc(eng)
print(Resource)
