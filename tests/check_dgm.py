#!python3

import os
import sys
import glob
import subprocess
import numpy as np

import dionysus as dion

top_dgm_dim = 2

def dion_to_numpy(dion_dgms):
    dgms = dict()
    for dim in range(top_dgm_dim + 1):
        dgms[dim] = np.array([ [p.birth, p.death] for p in dion_dgms[dim] ])
    return dgms


def numpy_to_dion(dgms):
    dion_dgms = dict()
    for dim in dgms:
        dion_dgm = dion.Diagram()
        for b, d in dgms[dim]:
            if d not in [np.inf, -np.inf]:
                dion_dgm.append(dion.DiagramPoint(b, d))
        dion_dgms[dim] = dion_dgm
    return dion_dgms


def get_dion_dgms(a, negate):
    fil = dion.fill_freudenthal(a, reverse=negate)
    p = dion.homology_persistence(fil)
    dion_dgms = dion.init_diagrams(p, fil)
    dgms = dion_to_numpy(dion_dgms)
    # to remove infinite points
    dion_dgms = numpy_to_dion(dgms)
    return dion_dgms, dgms


def read_cadmus_dgms(cadmus_dgm_fname):
    dgms = { dim : [] for dim in range(top_dgm_dim + 1) }
    for fname in glob.glob(f"{cadmus_dgm_fname}_*"):
        with open(fname) as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                dim, birth, death = line.split(";")
                dim = int(dim)
                birth = float(birth)
                death = float(death)
                if dim in dgms:
                    dgms[dim].append([birth, death])
    for dim in dgms:
        dgms[dim] = np.array(dgms[dim])
    dion_dgms = numpy_to_dion(dgms)
    return dion_dgms, dgms


def test(negate, n_ranks, input_fname, cadmus_dgm_fname):
    a = np.load(input_fname)
    correct_dion_dgms, correct_numpy_dgms = get_dion_dgms(a, negate)

    try:
        os.remove(cadmus_dgm_fname)
    except:
        pass

    cadmus_args = ["-b", str(n_ranks), "-l", "trace"]

    if negate:
        cadmus_args.append("-n")

    cadmus_args.append(input_fname)
    cadmus_args.append(cadmus_dgm_fname)

    cadmus_exec = f"./cadmus_{a.ndim}"
    subprocess.run(["mpirun", "-n", str(n_ranks), cadmus_exec] + cadmus_args)
    cadmus_dion_dgms, cadmus_numpy_dgms = read_cadmus_dgms(cadmus_dgm_fname)

    for dim in range(3):
        print(f"in dim {dim}, Cadmus has {len(cadmus_numpy_dgms[dim])} points, must be {len(correct_numpy_dgms[dim])}")
        print("Bottleneck distance = ", dion.bottleneck_distance(cadmus_dion_dgms[dim], correct_dion_dgms[dim]))
        assert(dion.bottleneck_distance(cadmus_dion_dgms[dim], correct_dion_dgms[dim]) < 1e-4)


if __name__ == "__main__":
    subprocess.run(["make", "-j"], check=True)

    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input_fname dgm_fname")

    input_fname = sys.argv[1]
    cadmus_dgm_fname = sys.argv[2]

    subprocess.run(f"rm -f {cadmus_dgm_fname}*", shell=True)

    if len(sys.argv) > 3:
        n_ranks = int(sys.argv[3])
    else:
        n_ranks = 4

    for negate in [True, False]:
        test(negate, n_ranks, input_fname, cadmus_dgm_fname)
