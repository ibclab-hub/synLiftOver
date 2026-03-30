#!/usr/bin/env python3.6

import sys
import re
import fileReader as fR
import filterFiles as fF
import mergeAxt as mA
import subprocess
from collections import defaultdict


def allocation (var1):
    var_convey = {}
    for i in var1:
        var_convey[i] = var1[i]
    return var_convey

