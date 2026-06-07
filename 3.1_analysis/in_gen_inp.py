import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import util
from util import log, err, wrn, ttl

# getting args

mols_dir = util.get_arg("mols", err=True)

mols_dir = os.path.join(os.getcwd(), mols_dir)

# recursively parse from dir

ttl(f"Generating Opt/Freq input files for ORCA")

lookup, get = util.load_molecules(mols_dir)
mols = get()

for mol in mols:
  log(f"Generating for   {mol.util_name} {' ' * (20 - len(mol.util_name))} @ {mol.path[len(os.getcwd()):]}...")
  comment = f"{mol.name} ({mol.base}.{mol.variant}])"
  util.dump_opt(os.path.dirname(mol.path), [], [], comment=comment, charge=mol.charge)
  util.dump_freq(os.path.dirname(mol.path), comment=comment, charge=mol.charge)
  guess_path = os.path.join(os.path.dirname(mol.path), "guess.xyz")
  if not os.path.exists(guess_path):
    with open(guess_path, "w") as file:
      file.write("1\n\nC 0 0 0")
