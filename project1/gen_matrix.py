'''
Author: shawn233
Date: 2021-06-23 20:40:20
LastEditors: shawn233
LastEditTime: 2021-06-26 13:05:40
Description: Generate arbitrary large matrix
'''

import os
import numpy as np


def write_matrix(X:np.ndarray, path:str) -> None:
    with open(path, "w") as f:
        f.write(f"{X.shape[0]} {X.shape[1]}\n")
        for i in range(X.shape[0]):
            f.write(" ".join(X[i].astype(str)))
            f.write("\n")


if __name__ == "__main__":
    os.makedirs("./matrix", exist_ok=True)

    # Exercise 2.1
    if not os.path.exists("./matrix/A_test.txt"):
        A = np.eye(4)
        B = (np.arange(16)/10).reshape(4, 4)
        write_matrix(A, "./matrix/A_test.txt")
        write_matrix(B, "./matrix/B_test.txt")

    for n in [2, 1024]:
        if not os.path.exists(f"./matrix/C_{n}.txt"):
            A = np.random.randn(n, n)
            B = np.random.randn(n, n)

            C = np.matmul(A, B)

            write_matrix(A, f"./matrix/A_{n}.txt")
            write_matrix(B, f"./matrix/B_{n}.txt")
            write_matrix(C, f"./matrix/C_{n}.txt")

    # Exercise 2.2

    if not os.path.exists("./matrix/conv_test.txt"):
        conv_test = (np.arange(64)/10).reshape(8, 8)
        write_matrix(conv_test, "./matrix/conv_test.txt")

    if not os.path.exists("./matrix/conv_test2.txt"):
        conv_test = np.random.randint(0, 10, (8,8))
        write_matrix(conv_test, "./matrix/conv_test2.txt")

    if not os.path.exists("./matrix/conv_1024.txt"):
        conv = np.random.randn(1024, 1024)
        write_matrix(conv, "./matrix/conv_1024.txt")

    if not os.path.exists("./matrix/kernel_test.txt"):
        kernel = np.array([[-1,0,1],[-1,0,1],[-1,0,1]], dtype=float)
        write_matrix(kernel, "./matrix/kernel_test.txt")

    if not os.path.exists("./matrix/kernel.txt"):
        kernel = np.random.randn(4, 4)
        write_matrix(kernel, "./matrix/kernel.txt")

    # Exercise 2.3

    if not os.path.exists("./matrix/pool_test.txt"):
        pool_test = np.random.randint(0, 10, (8,8))
        write_matrix(pool_test, "./matrix/pool_test.txt")

    if not os.path.exists("./matrix/pool_1024.txt"):
        pool = np.random.randn(1024, 1024)
        write_matrix(pool, "./matrix/pool_1024.txt")

    if not os.path.exists("./matrix/pool_kernel_2.txt"):
        pool_kernel = np.ones((2, 2), dtype=float)/4
        write_matrix(pool_kernel, "./matrix/pool_kernel_2.txt")

    if not os.path.exists("./matrix/pool_kernel_4.txt"):
        pool_kernel = np.ones((4, 4), dtype=float)/16
        write_matrix(pool_kernel, "./matrix/pool_kernel_4.txt")
