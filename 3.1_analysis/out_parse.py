import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import util
from util import log, err, wrn, ttl

# getting args

mols_dir = util.get_arg("mols", err=True)
data_dir = util.get_arg("data", err=True)
out_file = util.get_arg("out", err=True)

mols_dir = os.path.join(os.getcwd(), mols_dir)
data_dir = os.path.join(os.getcwd(), data_dir)
out_file = os.path.join(os.getcwd(), out_file)

# recursively parse from dir

ttl(f"Parsing Opt/Freq output files from ORCA")

lookup, get = util.load_molecules(mols_dir)
mols = get()

energies = {}

for mol in mols:
  data_path = os.path.join(data_dir, os.path.dirname(mol.path[len(mols_dir):]))
  log(f"Parsing   {mol.util_name} {' ' * (20 - len(mol.util_name))} @ {data_path}...")
  E, S, H, G = util.load_energies(data_path)
  if E is None:
    wrn(f"E could not be found!")
  if S is None:
    wrn(f"S could not be found!")
  if H is None:
    wrn(f"H could not be found!")
  if G is None:
    wrn(f"G could not be found!")
  energies[mol.util_name] = E, S, H, G

with open(out_file, "w") as file:
  file.write("Name, E (kcal/mol), S (kcal/mol), H (kcal/mol), G (kcal/mol)\n")
  for name in sorted(energies.keys()):
    E, S, H, G = energies[name]
    file.write(f"{name}, {E}, {S}, {H}, {G}\n")
