# Experimental: GPU kernels

Prototype CUDA kernels (Montgomery multiplication, tiled RNS math) captured
from early experiments. They are **not** built by the Makefile and are not
part of the library API.

Status: dormant. If this direction is picked up, wire a `nvcc`/`hipcc` build
rule here and link a separate `librabin_gpu` target — keep the CPU library
dependency-free.
