import os
import re

# generates guess.xyz files at each Opt.inp file
# since these input files have their models included in them

def parse(path):
    with open(path, "r") as file:
        content = file.read()
    lines = content.split("\n")
    atoms = []
    while len(lines) > 0:
        line = lines.pop(0)
        if line.startswith("* xyz"):
            while len(lines) > 0:
                line = lines.pop(0)
                end = False
                if "*" in line:
                    line = line.split("*")[0]
                    end = True
                if not line:
                    break
                line = re.split("\\s+", line)
                line[0] = line[0].split("(")[0]
                atoms.append(" ".join(line))
                if end:
                    break
    with open(os.path.join(os.path.dirname(path), "guess.xyz"), "w") as file:
        file.write(f"{len(atoms)}\n\n" + "\n".join(atoms))


def search_dir(wd):
    names = os.listdir(wd)
    for name in names:
        path = os.path.join(wd, name)
        if name == "Opt.inp":
            parse(path)
        if os.path.isdir(path):
            search_dir(path)


search_dir(os.path.join(os.path.dirname(__file__), "data"))
