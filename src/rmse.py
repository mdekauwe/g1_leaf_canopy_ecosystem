#!/usr/bin/env python
"""
Root Mean Squared Error

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (16.04.2015)"
__email__ = "mdekauwe@gmail.com"


import numpy as np
import sys

def rmse(x, y):
    return np.sqrt(np.mean((x - y)**2))
