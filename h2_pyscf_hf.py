import numpy as np
from openfermion import *
from openfermionpyscf import run_pyscf
geometry = [
    ['H', [0,0,0]],
    ['H', [0,0,0.7474]]
]
basis = 'sto-3g'
multiplicity = 1
charge = 0
molecule = MolecularData(geometry, basis, multiplicity, charge)
molecule2 = run_pyscf(molecule, run_scf=True)

one_body = molecule2.one_body_integrals
two_body = molecule2.two_body_integrals

one_body_integrals = []
two_body_integrals = []
for p in range(two_body.shape[0]):
    for q in range(two_body.shape[1]):
        one_body_integrals.append(f'({p} {q})  {one_body[p,q]}')
        for r in range(two_body.shape[2]):
            for s in range(two_body.shape[3]):
                two_body_integrals.append(f'({p} {q} {r} {s})  {two_body[p,q,r,s]}')

one_body_integral_str = "\n".join(one_body_integrals)
two_body_integral_str = "\n".join(two_body_integrals)

print(one_body_integral_str)
print(two_body_integral_str)

