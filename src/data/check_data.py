import sys
from os import path

import tqdm


def check_files(my_path, n_iter, n_beads, n_swarms):
    for iter in tqdm.tqdm(range(1, n_iter), desc="Iterations checked", total=n_iter):
        if not path.isfile(f"{my_path}/strings/string{iter}.txt"):
            print(f"Seem to be missing strings/string{iter}.txt")
        for bead in range(1, n_beads - 1):
            for swarm in range(0, n_swarms):
                if not path.isfile(
                    f"{my_path}/md/{iter}/{bead}/s{swarm}/traj_comp.xtc"
                ):
                    print(
                        f"Seem to be missing {my_path}/md/{iter}/{bead}/s{swarm}/traj_comp.xtc"
                    )
                if not path.isfile(f"{my_path}/md/{iter}/{bead}/s{swarm}/pullx.xvg"):
                    print(
                        f"Seem to be missing {my_path}/md/{iter}/{bead}/s{swarm}/pullx.xvg"
                    )
            if not path.isfile(f"{my_path}/md/{iter}/{bead}/restrained/pullx.xvg"):
                print(
                    f"Seem to be missing {my_path}/md/{iter}/{bead}/restrained/pullx.xvg"
                )
            if not path.isfile(f"{my_path}/md/{iter}/{bead}/restrained/confout.gro"):
                print(
                    f"Seem to be missing {my_path}/md/{iter}/{bead}/restrained/confout.gro"
                )

    return


def check_input_values(my_path, n_iter, n_beads, n_swarms):
    assert path.isdir(my_path), "The path provided does not exist"
    return


def main(my_path, n_iter, n_beads, n_swarms):
    check_input_values(my_path, n_iter, n_beads, n_swarms)
    check_files(my_path, n_iter, n_beads, n_swarms)
    return


if __name__ == "__main__":
    my_path = sys.argv[1]
    n_iter = int(sys.argv[2])
    n_beads = int(sys.argv[3])
    n_swarms = int(sys.argv[4])
    main(my_path, n_iter, n_beads, n_swarms)
