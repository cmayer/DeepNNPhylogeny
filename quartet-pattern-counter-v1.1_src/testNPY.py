#!/usr/bin/env python3

import numpy as np
import sys

argc = len(sys.argv)
if (argc != 2):
    print("Usage: testNPY.py npy-file\n")
    sys.exit(0);

filename = sys.argv[1]
v = np.load(filename)

print("Data type of object:       ", type(v))
print("Shape of numpy data:       ", v.shape)
print("Data type of numpy data:   ", v.dtype)

print(v)

