import os
import sys
import shutil
import heapq

import matplotlib.pyplot as plt
import numpy as np


root = os.path.dirname(os.path.abspath(__file__))


do_path = "-p" in sys.argv


def scenario(wd):
    inp_path = os.path.join(wd, "Opt.inp")
    scan_bonds = []
    scan_ranges = []
    with open(inp_path, "r") as inp_file:
        for line in inp_file:
            line = line.strip()
            if line.startswith("B"):
                parts = line.split()
                bond = (int(parts[1]), int(parts[2]))
                range_str = " ".join(parts[3:])
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
    Es: list[list] = [
        [None for _ in range(len(drange_x))] for _ in range(len(drange_y))
    ]
    total_points = len(drange_x) * len(drange_y)
    point_idx = 0
    for i, line in enumerate(lines):
        done = "                                *** OPTIMIZATION RUN DONE ***" in line
        failed = (
            "       The optimization did not converge but reached the maximum number of"
            in line
        )
        if done or failed:
            E = (
                None
                if failed
                else float(lines[i - 3][len("FINAL SINGLE POINT ENERGY") :].strip())
            )
            E_kcal = None if E is None else E * 627.509
            x = point_idx // len(drange_y)
            y = point_idx % len(drange_y)
            if y < len(drange_y) and x < len(drange_x):
                Es[y][x] = E_kcal
            point_idx += 1
            if point_idx >= total_points:
                break
    return {
        "drange_x": drange_x,
        "drange_y": drange_y,
        "bond_x": bond_x,
        "bond_y": bond_y,
        "Es": Es,
    }


def merge_scenarios(s1, s2):
    # Merge dranges
    drange_x = sorted(set(s1["drange_x"]) | set(s2["drange_x"]))
    drange_y = sorted(set(s1["drange_y"]) | set(s2["drange_y"]))
    bond_x = s1["bond_x"]
    bond_y = s1["bond_y"]
    # Create merged matrix filled with None
    Es_merged = [[None for _ in range(len(drange_x))] for _ in range(len(drange_y))]

    # Helper to fill from a scenario
    def fill_matrix(s):
        x_map = {v: i for i, v in enumerate(drange_x)}
        y_map = {v: i for i, v in enumerate(drange_y)}
        for iy, yval in enumerate(s["drange_y"]):
            for ix, xval in enumerate(s["drange_x"]):
                merged_y = y_map[yval]
                merged_x = x_map[xval]
                val = s["Es"][iy][ix]
                if val is not None:
                    Es_merged[merged_y][merged_x] = val

    fill_matrix(s1)
    fill_matrix(s2)
    return {
        "drange_x": drange_x,
        "drange_y": drange_y,
        "bond_x": bond_x,
        "bond_y": bond_y,
        "Es": Es_merged,
    }


def lowest_peak_path(Es, xi, yi, xf, yf, allow_diag=False):
    """
    Finds the path from (xi, yi) to (xf, yf) that minimizes
    1. the maximum energy value along the path (lowest peak)
    2. among those, the shortest path length.
    """
    nrows, ncols = Es.shape

    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    if allow_diag:
        directions += [(-1, -1), (-1, 1), (1, -1), (1, 1)]

    if np.isnan(Es[yi, xi]) or np.isnan(Es[yf, xf]):
        raise ValueError("Start or end position lies on a NaN cell.")

    pq = [(Es[yi, xi], 0.0, xi, yi)]
    best = np.full((nrows, ncols, 2), np.inf)
    best[xi, yi] = (Es[yi, xi], 0.0)

    parent = {(xi, yi): None}

    while pq:
        peak, length, x, y = heapq.heappop(pq)

        if (x, y) == (xf, yf):
            # reconstruct path
            path = []
            curr = (x, y)
            while curr:
                path.append(curr)
                curr = parent[curr]
            path.reverse()
            return {"path": path, "peak": float(peak), "length": length}

        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if not (0 <= nx < nrows and 0 <= ny < ncols):
                continue

            val = Es[ny, nx]
            if np.isnan(val):
                continue

            # compute penalty for diagonal moves
            step_cost = (dx**2 + dy**2) ** 0.5

            npeak = max(peak, val)
            nlen = length + step_cost

            if (npeak < best[nx, ny, 0]) or (
                npeak == best[nx, ny, 0] and nlen < best[nx, ny, 1]
            ):
                best[nx, ny] = (npeak, nlen)
                parent[(nx, ny)] = (x, y)
                heapq.heappush(pq, (npeak, nlen, nx, ny))

    return None  # no path found


def plot_matrix(result, title="Heatmap", file_path=None, show_plot=True):
    drange_x = result["drange_x"]
    drange_y = result["drange_y"]
    bond_x = result["bond_x"]
    bond_y = result["bond_y"]
    Es = np.array(result["Es"], dtype=float)
    Es[np.array(result["Es"]) == None] = np.nan

    # Zero-based energy levels
    min_energy = np.nanmin(Es)  # Find the minimum energy ignoring NaN
    Es -= min_energy  # Subtract the minimum energy to zero the matrix

    plt.figure(figsize=(8, 6))
    plt.imshow(Es, cmap="viridis", origin="lower", aspect="auto")
    plt.colorbar(
        label="Relative Energy (kcal/mol)"
    )  # Update label to reflect zero-based energy
    plt.xticks(ticks=range(len(drange_x)), labels=drange_x)
    plt.yticks(ticks=range(len(drange_y)), labels=drange_y)
    plt.xlabel(f"distance (Å) atoms {bond_x[0]}-{bond_x[1]}")
    plt.ylabel(f"distance (Å) atoms {bond_y[0]}-{bond_y[1]}")

    if do_path:
        # Highlight local minima (including edges)
        minima = []
        for i in range(Es.shape[0]):
            for j in range(Es.shape[1]):
                if not np.isnan(Es[i, j]) and Es[i, j] < 70:
                    neighbors = []
                    for di in [-1, 0, 1]:
                        for dj in [-1, 0, 1]:
                            if di == 0 and dj == 0:
                                continue
                            ni, nj = i + di, j + dj
                            if 0 <= ni < Es.shape[0] and 0 <= nj < Es.shape[1]:
                                neighbors.append(Es[ni, nj])
                    if all(Es[i, j] <= n for n in neighbors if not np.isnan(n)):
                        if all(np.isnan(n) for n in neighbors):
                            continue
                        minima.append((i, j))
        for i, j in minima:
            plt.scatter(j, i, edgecolor="red", facecolor="none", s=100, linewidth=1.5)
            plt.text(
                j,
                i + 0.25,
                str(round(Es[i, j] * 1e4) / 1e4),
                color="red",
                horizontalalignment="center",
                verticalalignment="baseline",
            )
        if len(minima) == 2:
            yi, xi = minima[0]
            yf, xf = minima[1]
            data = lowest_peak_path(Es, xi, yi, xf, yf, allow_diag=True)
            path = data["path"]
            peak = data["peak"]
            plt.plot([x for x, _ in path], [y for _, y in path], color="red")
            for x, y in path:
                if Es[y, x] == peak:
                    plt.scatter(
                        x, y, edgecolor="red", facecolor="none", s=100, linewidth=1.5
                    )
                    plt.text(
                        x,
                        y + 0.25,
                        str(round(Es[y, x] * 1e4) / 1e4),
                        color="red",
                        horizontalalignment="center",
                        verticalalignment="baseline",
                    )

    plt.title(title)
    plt.tight_layout()

    if show_plot:
        plt.show()
    elif file_path:
        plt.savefig(file_path)
    plt.close()


def run_all_scenarios(data_dir, figures_dir, show_plot=True):
    if not show_plot:
        if os.path.exists(figures_dir):
            shutil.rmtree(figures_dir)  # Clear the folder using shutil.rmtree
        os.makedirs(figures_dir)

    manual_paths = [
        "A_Constrain",
        "A_Relax",
        "B_Relax",
        "C_Constrain",
        "C_Relax",
        "N_Constrain",
        "N_Relax",
        "DifferentHbond/2.7",
        "DifferentHbond/2.9",
        "DifferentHbond/3.1",
        "DifferentHbond/3.3",
        "Protein_LongHbond_Constrain",
        "Protein_Relax",
    ]

    # Run scenarios for all manual paths
    for rel_path in manual_paths:
        wd = os.path.join(data_dir, rel_path)
        if os.path.exists(os.path.join(wd, "Opt.inp")) and os.path.exists(
            os.path.join(wd, "Opt.out")
        ):
            result = scenario(wd)
            file_path = (
                os.path.join(figures_dir, f"{rel_path.replace('/', '_')}.png")
                if not show_plot
                else None
            )
            plot_matrix(
                result,
                title=f"Heatmap of Energies for {rel_path}",
                file_path=file_path,
                show_plot=show_plot,
            )

    # Special case: DifferentHbond/3.5 and DifferentHbond/3.7 with Continue subfolders
    for subfolder in ["3.5", "3.7"]:
        base_path = os.path.join(data_dir, "DifferentHbond", subfolder)
        continue_path = os.path.join(base_path, "Continue")
        if (
            os.path.exists(os.path.join(base_path, "Opt.inp"))
            and os.path.exists(os.path.join(base_path, "Opt.out"))
            and os.path.exists(os.path.join(continue_path, "Opt.inp"))
            and os.path.exists(os.path.join(continue_path, "Opt.out"))
        ):
            result1 = scenario(base_path)
            result2 = scenario(continue_path)
            merged = merge_scenarios(result1, result2)
            file_path = (
                os.path.join(figures_dir, f"DifferentHbond_{subfolder}_merged.png")
                if not show_plot
                else None
            )
            plot_matrix(
                merged,
                title=f"Heatmap of Energies for DifferentHbond/{subfolder} (merged)",
                file_path=file_path,
                show_plot=show_plot,
            )


if __name__ == "__main__":
    user_input = input("Do you want to display plots? (y/n): ").strip().lower()
    show_plot = user_input == "y"
    data_dir = os.path.join(root, "data")
    figures_dir = os.path.join(root, "figures")
    run_all_scenarios(data_dir, figures_dir, show_plot=show_plot)
