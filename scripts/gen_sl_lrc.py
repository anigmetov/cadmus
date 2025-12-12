#!python3

import glob


def get_text(partition: str="lr6", n_nodes: int=1, n_tasks_per_node: int=32, 
             dataset_name:str ="isotropic_pressure",
             res: int = 128, # resolution of dataset
             n_blocks=8,
             swap_reduction: bool=True, clearing_opt: bool=True,
             sl_fname=""):

    fname=f"/global/home/users/anigmetov/SCRATCH/Data/selected_klacansky/{dataset_name}_{res}x{res}x{res}_float32.npy"

    options = ""
    opt_part = ""

    if swap_reduction:
        options += " -s "
        opt_part += "_swap"

    if clearing_opt:
        options += " -c "
        opt_part += "_clear"

    log_fname = f"log_{dataset_name}_{res}{opt_part}_{n_blocks}_%j.txt"

    if res == 128:
        maxtime = "00:15:00"
    elif res == 256:
        maxtime = "08:00:00"

    text = f"""#!/bin/bash -l

#SBATCH --job-name=cadmus-{dataset_name}-{res}
#SBATCH --partition={partition}
#SBATCH --account=pc_mlgeometry
#SBATCH --qos=lr_normal
#SBATCH --time={maxtime}
#SBATCH --nodes={n_nodes}
#SBATCH --ntasks-per-node={n_tasks_per_node}
##SBATCH --mail-user=anigmetov@lbl.gov
##SBATCH --mail-type=ALL
#SBATCH --error={log_fname}
#SBATCH --output={log_fname}

source $HOME/code/spack/share/spack/setup-env.sh
spack env activate mpich

cd /global/home/users/anigmetov/SCRATCH/cadmus/src

export CALI_CONFIG=runtime-report,calc.inclusive,mem.highwatermark

nblocks={n_blocks}
nranks={n_blocks}

mpirun -n $nranks ./cadmus_3 -b $nblocks -n {options} {fname} /dev/null
"""
    return text, log_fname


def generate_scripts(res):
    with open("run_me.sh", "w") as sbatch_f:
        sbatch_f.write("#!/bin/bash\n")
        for n_ranks in [8, 16, 32, 64]:
            for clearing_opt in [True]:
                for swap_reduction in [True]:
                    ds_name = "isotropic_pressure"
                    sl_file_name = f"submit_{ds_name}_{res}_b_{n_ranks}_clear_{clearing_opt}_swap_{swap_reduction}.sh"

                    sbatch_f.write(f"sbatch {sl_file_name}\n")

                    if res == 128:

                        if n_ranks <= 32:
                            n_nodes = 1
                        elif n_ranks <= 64:
                            n_nodes = 2
                        elif n_ranks <= 128:
                            n_nodes = 4
                        n_tasks_per_node=32

                    elif res == 256:
                        if n_ranks <= 32:
                            n_nodes = 4
                        elif n_ranks <= 64:
                            n_nodes = 8
                        elif n_ranks <= 128:
                            n_nodes = 16
                        n_tasks_per_node=8

                    with open(sl_file_name, "w") as f:
                        text, _ = get_text(partition="lr6", n_nodes=n_nodes, dataset_name=ds_name,
                                           n_tasks_per_node=n_tasks_per_node,
                                           res=res, n_blocks=n_ranks, swap_reduction=swap_reduction, clearing_opt=clearing_opt)
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
    fnames = os.glob(fname.replace("%j", "*"))
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


def analyze(n_ranks):
    for clearing_opt in [True, False]:
        for swap_reduction in [True, False]:
            print(f"Clearing: {clearing_opt}. Swap reduction: {swap_reduction}")
            for n_ranks in [8, 16, 32, 64, 128]:
                res = 256
                ds_name = "isotropic_pressure"
                _, log_fname = get_text(partition="lr6", dataset_name=ds_name,
                                        res= res, n_blocks=n_ranks, swap_reduction=swap_reduction, clearing_opt=clearing_opt)
                n_blocks, n_simplices, local_init, local_reduction, sparsify, exchange, final_sort, final_reduction, total = analyze_file(log_fname)
                print(f"{n_ranks}; {final_reduction}; {total}")


if __name__ == "__main__":
    generate_scripts(256)
