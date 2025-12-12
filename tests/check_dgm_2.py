#!python3

import subprocess
import os

import numpy as np

import dionysus as dion

top_dgm_dim = 1

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
            dion_dgm.append(dion.DiagramPoint(b, d))
        dion_dgms[dim] = dion_dgm
    return dion_dgms


def get_dion_dgms(a):
    fil = dion.fill_freudenthal(a)
    p = dion.homology_persistence(fil)
    dion_dgms = dion.init_diagrams(p, fil)
    dgms = dion_to_numpy(dion_dgms)
    return dion_dgms, dgms


def read_cadmus_dgms(cadmus_dgm_fname):
    dgms = { dim : [] for dim in range(top_dgm_dim + 1) }
    with open(cadmus_dgm_fname) as f:
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

def run_test(n, seed, clearing, red):
    np.random.seed(seed)
    a = np.random.randint(1, 51, size=(n, n)).astype(np.float32)
    input_fname = f"random_{n}.npy"
    np.save(input_fname, a)
    correct_dion_dgms, correct_numpy_dgms = get_dion_dgms(a)
    cadmus_dgm_fname = "dgm_cadmus.txt"
    try:
        os.remove(cadmus_dgm_fname)
    except:
        pass
    n_blocks = 4
    n_ranks = 1
    compilation_res = subprocess.run(["make", "-j"], check=True)
    if clearing:
        subprocess.run(["mpirun", "-n", str(n_ranks), "./cadmus_2", "-l", "info", "-c", "-b", str(n_blocks), "-r", red, input_fname, cadmus_dgm_fname])
    else:
        subprocess.run(["mpirun", "-n", str(n_ranks), "./cadmus_2", "-l", "info", "-b", str(n_blocks), "-r", red, input_fname, cadmus_dgm_fname])
    cadmus_dion_dgms, cadmus_numpy_dgms = read_cadmus_dgms(cadmus_dgm_fname)

    # for dim in range(top_dgm_dim + 1):
    #     print(f"in dim {dim}, Cadmus has {len(cadmus_numpy_dgms[dim])} points, must be {len(correct_numpy_dgms[dim])}")
    #     # print(dim, correct_numpy_dgms[dim])
    #     sorted_correct_dgm = correct_numpy_dgms[dim][correct_numpy_dgms[dim][:, 0].argsort()]
    #     sorted_cadmus_dgm = cadmus_numpy_dgms[dim][cadmus_numpy_dgms[dim][:, 0].argsort()]
    #     print("Correct:\n", sorted_correct_dgm)
    #     print("Cadmus:\n", sorted_cadmus_dgm)

    for dim in range(top_dgm_dim + 1):
        print(f"in dim {dim}, Cadmus has {len(cadmus_numpy_dgms[dim])} points, must be {len(correct_numpy_dgms[dim])}")
        # # print(dim, correct_numpy_dgms[dim])
        # sorted_correct_dgm = correct_numpy_dgms[dim][correct_numpy_dgms[dim][:, 0].argsort()]
        # sorted_cadmus_dgm = cadmus_numpy_dgms[dim][cadmus_numpy_dgms[dim][:, 0].argsort()]
        # print("Correct:\n", sorted_correct_dgm)
        # print("Cadmus:\n", sorted_cadmus_dgm)
        assert(dion.bottleneck_distance(cadmus_dion_dgms[dim], correct_dion_dgms[dim]) < 1e-6)



def test_2():
    for clearing in [True, False]:
        # for red in ["naive_elz", "elz", "cascade", "cascade_parallel"]:
        for red in ["cascade_parallel"]:
            for n in [4, 8]:
                for seed in [2, 3]:
                    run_test(n, seed, clearing, red)




if __name__ == "__main__":
    test_2()

