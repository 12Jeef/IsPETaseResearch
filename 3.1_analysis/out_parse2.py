import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import util
from util import log, err, wrn, ttl

# getting args

mols_dir = util.get_arg("mols", err=True)
rxns_file = util.get_arg("rxns", err=True)
energies_file = util.get_arg("energies", err=True)
out_dir = util.get_arg("out", err=True)

mols_dir = os.path.join(os.getcwd(), mols_dir)
rxns_file = os.path.join(os.getcwd(), rxns_file)
energies_file = os.path.join(os.getcwd(), energies_file)
out_dir = os.path.join(os.getcwd(), out_dir)

# recursively parse from dir

ttl(f"Parsing Opt/Freq energies")

def floatish(v):
  try:
    return float(v)
  except:
    pass
  return None

lookup, get = util.load_molecules(mols_dir)
rxns = util.load_reactions(rxns_file)

energies = {}

with open(energies_file, "r") as file:
  content = file.read()
lines = [[part.strip() for part in line.split(",")] for line in content.split("\n")[1:] if len(line) > 0]
for name, E, S, H, G in lines:
  E = floatish(E)
  S = floatish(S)
  H = floatish(H)
  G = floatish(G)
  mol = lookup(name)
  if mol is None:
    continue
  energies[mol] = (E, S, H, G)

rxn_energies = {}
contents = []
n_cols_max = 0
for rxn_name, rxn in rxns.items():
  log(f"Building reaction {rxn_name}...")
  rxn_energies[rxn_name] = util.get_reaction_delta_energies(lookup, energies, rxn)

  left, right = rxn
  left_array, right_array = util.get_reaction_delta_energies_array(lookup, energies, rxn)

  id_name = rxn_name.split(":")[0].strip()
  text_name = rxn_name.split(":")[-1].strip()
  with open(os.path.join(out_dir, f"reaction_{id_name}_energies.csv"), "w") as file:
    file.write(f"{text_name},Left{',,'*len(left_array)}→,Right{',,'*len(right_array)},Δ,\n")

    file.write(f"Name,")
    for count, name in left:
      file.write(f"{count},{name},")
    file.write(f",")
    for count, name in right:
      file.write(f"{count},{name},")
    file.write(f",,\n")
    
    for i in range(4):
      file.write(f"{'ESHG'[i]} (kcal/mol),")
      for slice in left_array:
        file.write(f",{slice[i]},")
      file.write(f",")
      for slice in right_array:
        file.write(f",{slice[i]},")
      file.write(f",{rxn_energies[rxn_name][i]},\n")
  
  with open(os.path.join(out_dir, f"reaction_{id_name}_energies.csv"), "r") as file:
    content = file.read()
    for line in content.split("\n"):
      n_cols_max = max(n_cols_max, line.count(","))
    contents.append(content)

with open(os.path.join(out_dir, f"reaction_all_energies.csv"), "w") as file:
  for content in contents:
    file.write(content)
    for i in range(n_cols_max):
      file.write(f"\\{'='*100},")
    file.write("\n")

with open(os.path.join(out_dir, "reaction_energies.csv"), "w") as file:
  file.write("Reaction Name, ΔE (kcal/mol), ΔS (kcal/mol), ΔH (kcal/mol), ΔG (kcal/mol)\n")
  for rxn_name, (E, S, H, G) in rxn_energies.items():
    text_name = rxn_name.split(":")[-1].strip()
    file.write(f"{text_name}, {E}, {S}, {H}, {G}\n")
