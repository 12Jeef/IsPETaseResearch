import os

import matplotlib.pyplot as plt
import numpy as np


# LEGACY!!! DO NOT USE!!!


root = os.path.dirname(os.path.abspath(__file__))


def scenario(wd):
    name = os.path.basename(wd)
    inp_path = os.path.join(wd, "Opt.inp")
    scan_bonds = []
    scan_ranges = []
    with open(inp_path, "r") as inp_file:
        for line in inp_file:
            line = line.strip()
            if line.startswith("B"):
                parts = line.split()
                bond = (int(parts[1]), int(parts[2]))
                # Join all parts after the atom indices to reconstruct the bracketed range
                range_str = " ".join(parts[3:])
                # Extract numbers between [ and ]
                range_nums = range_str[
                    range_str.find("[") + 1 : range_str.find("]")
                ].split()
                rng = [float(x) for x in range_nums]
                scan_bonds.append(bond)
                scan_ranges.append(rng)
    drange_x = scan_ranges[0]
    drange_y = scan_ranges[1]
    bond_x = scan_bonds[0]
    bond_y = scan_bonds[1]

    out = os.path.join(wd, "Opt.out")
    with open(out, "r") as file:
        content = file.read()
    lines = content.split("\n")
    Es = [[None for _ in range(len(drange_x))] for _ in range(len(drange_y))]
    total_points = len(drange_x) * len(drange_y)
    point_idx = 0
    for i, line in enumerate(lines):
        if "                                *** OPTIMIZATION RUN DONE ***" in line:
            E = float(lines[i - 3][len("FINAL SINGLE POINT ENERGY") :].strip())
            E_kcal = E * 627.509  # Convert Hartree to kcal/mol
            x = point_idx // len(drange_y)
            y = point_idx % len(drange_y)
            if y < len(drange_y) and x < len(drange_x):
                Es[y][x] = E_kcal
            point_idx += 1
            if point_idx >= total_points:
                break
    # Convert None to np.nan for plotting
    Es = np.array(Es, dtype=float)
    Es[np.array(Es) == None] = np.nan
    # Plot heatmap
    plt.figure(figsize=(8, 6))
    plt.imshow(Es, cmap="viridis", origin="lower", aspect="auto")
    plt.colorbar(label="Energy (kcal/mol)")
    plt.xticks(ticks=range(len(drange_x)), labels=drange_x)
    plt.yticks(ticks=range(len(drange_y)), labels=drange_y)
    plt.xlabel(f"distance (Å) atoms {bond_x[0]}-{bond_x[1]}")
    plt.ylabel(f"distance (Å) atoms {bond_y[0]}-{bond_y[1]}")
    rel_path = os.path.relpath(wd, os.path.join(root, "data"))
    plt.title(f"Heatmap of Energies for {rel_path} scenario")
    plt.tight_layout()
    plt.show()


def run_all_scenarios(data_dir):
    for root, dirs, files in os.walk(data_dir):
        # Only run scenario if Opt.inp and Opt.out exist in the directory
        if "Opt.inp" in files and "Opt.out" in files:
            scenario(root)


data_path = os.path.join(root, "data")
run_all_scenarios(data_path)
