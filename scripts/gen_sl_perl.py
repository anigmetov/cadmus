#!/usr/bin/env python3

import glob


def get_text(n_cpus_per_task: int=2,
             dataset_name:str ="isotropic_pressure",
             res: int = 128, # resolution of dataset
             n_blocks=8,
             swap_reduction: bool=True, clearing_opt: bool=True,
             sl_fname=""):

    fname=f"/pscratch/sd/g/greynarn/Data/klacansky_datasets/{dataset_name}_{res}x{res}x{res}_float32.npy"

    options = ""
    opt_part = ""

    n_nodes = max(int(n_blocks * n_cpus_per_task / 256), 1)

    if swap_reduction:
        options += " -s "
        opt_part += "_swap"

    if clearing_opt:
        options += " -c "
        opt_part += "_clear"

    log_fname = f"log_{dataset_name}_{res}{opt_part}_{n_blocks}_%j.txt"

    text = f"""#!/bin/bash -l

#SBATCH --job-name=cadmus-{dataset_name}-{res}
#SBATCH --constraint=cpu
#SBATCH --account=m636
#SBATCH --qos=regular
#SBATCH -L scratch
#SBATCH --time=0:30:00
#SBATCH --nodes={n_nodes}
#SBATCH --ntasks {n_blocks}
#SBATCH --cpus-per-task={n_cpus_per_task}
#SBATCH --error={log_fname}
#SBATCH --output={log_fname}

cd /pscratch/sd/g/greynarn/cadmus/src

export CALI_CONFIG=runtime-report

nblocks={n_blocks}

module load PrgEnv-gnu

srun --cpu-bind=cores ./cadmus_3 -b $nblocks -n {options} {fname} /dev/null
"""
    return text, log_fname


def generate_scripts():
    with open("run_me.sh", "w") as sbatch_f:
        sbatch_f.write("#!/bin/bash\n")
        for n_ranks in [8, 16, 32, 64, 128]:
            for clearing_opt in [True, False]:
                for swap_reduction in [True, False]:
                    res = 128
                    ds_name = "isotropic_pressure"
                    sl_file_name = f"submit_{ds_name}_{res}_b_{n_ranks}_clear_{clearing_opt}_swap_{swap_reduction}.sh"

                    sbatch_f.write(f"sbatch {sl_file_name}\n")

                    if res == 64:
                        n_cpus_per_task = 2
                    elif res == 128:
                        n_cpus_per_task = 4

                    with open(sl_file_name, "w") as f:
                        text, _ = get_text(dataset_name=ds_name, n_cpus_per_task=n_cpus_per_task,
                                        res= res, n_blocks=n_ranks, swap_reduction=swap_reduction, clearing_opt=clearing_opt)
                        f.write(text)


def str_to_bool(s):
    s1 = s.strip().lower()
    if s1 == "false":
        return False
    elif s1 == "true":
        return True
    else:
        raise RuntimeError(f"Unknown bool value: {s}")


def analyze_file(fname):
    fnames = glob.glob(fname.replace("%j", "*"))
    if len(fnames) != 1:
        raise RuntimeError(f"Expected exactly one file with pattern {fname}, found {fnames}")
    with open(fnames[0], "r") as f:
        found = False
        for line in f:
            if found:
                _, n_blocks, n_simplices, clearing_opt, swap_reduction, local_init, local_reduction, sparsify, exchange, final_sort, final_reduction, total = line.split(";")

                n_blocks = int(n_blocks)
                n_simplices = int(n_simplices)
                clearing_opt = str_to_bool(clearing_opt)
                swap_reduction = str_to_bool(swap_reduction)
                local_init = float(local_init)
                local_reduction = float(local_reduction)
                sparsify = float(sparsify)
                exchange = float(exchange)
                final_sort = float(final_sort)
                final_reduction = float(final_reduction)
                total = float(total)
                return n_blocks, n_simplices, local_init, local_reduction, sparsify, exchange, final_sort, final_reduction, total
            found = "clearing_opt; swap_reduction;" in line


def analyze():
    for clearing_opt in [True, False]:
        for swap_reduction in [True, False]:
            print(f"Clearing: {clearing_opt}. Swap reduction: {swap_reduction}")
            for n_ranks in [8, 16, 32, 64, 128]:
                res = 64
                ds_name = "isotropic_pressure"
                _, log_fname = get_text(dataset_name=ds_name,
                                        res= res, n_blocks=n_ranks, swap_reduction=swap_reduction, clearing_opt=clearing_opt)
                n_blocks, n_simplices, local_init, local_reduction, sparsify, exchange, final_sort, final_reduction, total = analyze_file(log_fname)
                print(f"{n_ranks}; {final_reduction}; {total}")


if __name__ == "__main__":
    generate_scripts()
