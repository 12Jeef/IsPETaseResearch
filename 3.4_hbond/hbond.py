import sys
import shutil
import os
import re
import matplotlib.pyplot as plt

root = os.path.dirname(os.path.abspath(__file__))

# double-bonded oxygen and carbon

OXYGEN_N = 13
CARBON_N = 0

# (nitrogen_n, (other_n, ...))

METHYLAMINE_NS = [
    (14, (16, 17, 18, 19, 20, 21)),
    (15, (22, 23, 24, 25, 26, 27)),
]

# get mode

mode_in = "-i" in sys.argv
mode_out = "-o" in sys.argv

# input mode callback

Di = 2.5
Df = 4.5
Dstep = 0.25


def mode_in_fn():
    with open(os.path.join(root, "in", "guess.xyz"), "r") as file:
        content = file.read()
    lines = content.split("\n")
    n = int(lines.pop(0))
    comment = lines.pop(0)

    other_atoms = []
    methylamines_atoms = [[None] for _ in range(len(METHYLAMINE_NS))]

    oxygen = None

    for i in range(n):
        elem, x_s, y_s, z_s = re.split("\\s+", lines[i])
        x = float(x_s)
        y = float(y_s)
        z = float(z_s)

        for j, (nitrogen_n, other_ns) in enumerate(METHYLAMINE_NS):
            if i == nitrogen_n:
                methylamines_atoms[j][0] = (i, elem, x, y, z)
                break
            if i in other_ns:
                methylamines_atoms[j].append((i, elem, x, y, z))
                break
        else:
            other_atoms.append((i, elem, x, y, z))

        if i == OXYGEN_N:
            oxygen = (i, elem, x, y, z)

    if oxygen is None:
        raise Exception("Expected oxygen to be not None")
    for i, methylamine_atoms in enumerate(methylamines_atoms):
        nitrogen = methylamine_atoms[0]
        if nitrogen is None:
            raise Exception(f"Expected methylamine {i+1} to have a nitrogen")

    def modify_atoms(D):
        atoms = [*other_atoms]
        _, __, ox, oy, oz = oxygen
        for methylamine_atoms in methylamines_atoms:
            nitrogen, *others = methylamine_atoms
            ni, nelem, nx, ny, nz = nitrogen
            d = ((ox - nx) ** 2 + (oy - ny) ** 2 + (oz - nz) ** 2) ** 0.5
            sx = (nx - ox) * ((D / d) - 1)
            sy = (ny - oy) * ((D / d) - 1)
            sz = (nz - oz) * ((D / d) - 1)
            atoms.append((ni, nelem, nx + sx, ny + sy, nz + sz))
            atoms.extend(
                [(i, elem, x + sx, y + sy, z + sz) for i, elem, x, y, z in others]
            )
        atoms.sort(key=lambda atom: atom[0])
        return atoms

    try:
        shutil.rmtree(os.path.join(root, "in", "input"))
    except:
        pass
    os.mkdir(os.path.join(root, "in", "input"))

    for i in range(int((Df - Di) / Dstep) + 1):
        D = Di + Dstep * i
        os.mkdir(os.path.join(root, "in", "input", f"{D}"))
        atoms = modify_atoms(D)
        with open(os.path.join(root, "in", "input", f"{D}", "guess.xyz"), "w") as file:
            file.write(f"{len(atoms)}\n{comment}\n")
            for _, elem, x, y, z in atoms:
                file.write(f"{elem} {x} {y} {z}\n")
        shutil.copyfile(
            os.path.join(root, "in", "Opt.inp"), os.path.join(root, "in", "input", f"{D}", "Opt.inp")
        )


def mode_out_fn():
    def parse_row(s):
        return re.split("\\s\\s+", str(s).strip())

    other_net_mulliken_charges = []
    other_c_mulliken_charges = []
    other_o_mulliken_charges = []
    methylamine_net_mulliken_charges = []
    double_bond_strengths = []

    for i in range(int((Df - Di) / Dstep) + 1):
        D = Di + Dstep * i
        with open(os.path.join(root, "out", f"{D}", "Opt.out"), "r") as file:
            content = file.read()
        lines = content.split("\n")
        lines = lines[lines.index("MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS") + 2 :]

        charges = [
            float(line.split(":")[1].strip().split(" ")[0]) for line in lines[:28]
        ]

        other_net_mulliken_charge = sum(charges[:14])
        other_c_mulliken_charge = charges[0]
        other_o_mulliken_charge = charges[13]
        methylamine_net_mulliken_charge = sum(charges[14:])

        double_bond_strength = 0

        lines = lines[lines.index("  Mayer bond orders larger than 0.100000") + 1 :]
        while True:
            line = lines.pop(0)
            if line.strip() == "":
                break
            for part in line.split("B("):
                part = part.strip()
                if len(part) <= 0:
                    continue
                bond, strength = part.split(":")
                atom_a, atom_b = bond.strip()[:-1].strip().split(",")
                atom_ai, *_ = atom_a.strip().split("-")
                atom_bi, *_ = atom_b.strip().split("-")
                atom_ai = int(atom_ai)
                atom_bi = int(atom_bi)
                strength = float(strength.strip())
                if (
                    atom_ai == OXYGEN_N
                    and atom_bi == CARBON_N
                    or atom_ai == CARBON_N
                    and atom_bi == OXYGEN_N
                ):
                    double_bond_strength = strength

        print("=" * 40)
        print(f"Distance            {D}")
        print(f"Other Charge        {other_net_mulliken_charge}")
        print(f"Methylamine Charge  {methylamine_net_mulliken_charge}")
        print(f"Ester Double Bond   {double_bond_strength}")

        other_net_mulliken_charges.append(other_net_mulliken_charge)
        other_c_mulliken_charges.append(other_c_mulliken_charge)
        other_o_mulliken_charges.append(other_o_mulliken_charge)
        methylamine_net_mulliken_charges.append(methylamine_net_mulliken_charge)
        double_bond_strengths.append(double_bond_strength)

    x = [Di + Dstep * i for i in range(int((Df - Di) / Dstep) + 1)]
    y0 = double_bond_strengths
    y1 = other_c_mulliken_charges
    y2 = other_o_mulliken_charges

    print(max(other_c_mulliken_charges) - min(other_c_mulliken_charges))
    print(max(other_o_mulliken_charges) - min(other_o_mulliken_charges))

    # plt.plot(x, y0)
    plt.plot(x, y1)
    plt.plot(x, y2)
    plt.show()


# run callbacks

if mode_in:
    mode_in_fn()

if mode_out:
    mode_out_fn()
