# -*- coding: utf-8 -*-
"""
Created on Fri May 22 10:09:29 2020

@author: ravir
"""
import numpy as np
import matplotlib.pyplot as plt
import pytest
# pytest.importorskip('pycuda')
# from lbmpy.session import *
# from lbmpy.relaxationrates import *
import sympy as sp

τ_0, ρ, ω, ω_total, ω_0 = sp.symbols("tau_0 rho omega omega_total omega_0", positive=True, real=True)
ν_0, C_S, S, Π = sp.symbols("nu_0, C_S, |S|, Pi", positive=True, real=True)

Seq = sp.Eq(S, 3 * ω / 2 * Π)
print(Seq)

def relaxation_rate_from_lattice_viscosity(ν_0, C_S,S):
    omega = 2/(6*(ν_0 + C_S ** 2 * S) + 1)
    return omega

ω = relaxation_rate_from_lattice_viscosity(ν_0, C_S,S)

Seq2 = Seq.subs(ω, relaxation_rate_from_lattice_viscosity(ν_0, C_S,S ))
print(Seq2)

solveRes = sp.solve(Seq2, S)
print(solveRes)
assert len(solveRes) == 1
SVal = solveRes[0]

def lattice_viscosity_from_relaxation_rate(τ_0):
    ω_0 = 1/τ_0
    return ω_0

SVal = SVal.subs(ν_0, lattice_viscosity_from_relaxation_rate(τ_0)).expand()

def second_order_moment_tensor(function_values, stencil):
    assert len(function_values) == len(stencil)
    dim = len(stencil[0])
    return sp.Matrix(dim, dim, lambda i, j: sum(c[i] * c[j] * f for f, c in zip(function_values, stencil)))


def frobenius_norm(matrix, factor=1):
    return sp.sqrt(sum(i*i for i in matrix) * factor)

