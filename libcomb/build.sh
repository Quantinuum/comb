#!/bin/bash
# static linking
gcc demo.c -o demo /usr/local/lib/libm4ri.a -lpng -lc -lm
#dynamic linking
#gcc demo.c -lm4ri
# LD_LIBRARY_PATH=/usr/local/lib at runtime
