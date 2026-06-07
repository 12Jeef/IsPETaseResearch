import sys
import os
import re

# sys utils

args = sys.argv[1:]

def get_args(name):
  values = []
  for i, arg in enumerate(args):
    if arg == f"--{name}":
      values.append(args[i + 1])
  return values

def get_arg(name, default=None, err=False):
  values = get_args(name)
  if len(values) > 0:
    return values[-1]
  if err:
    raise Exception(f"no --{name} argument(s) defined")
  return default

# class refs

class Residue:
  def __init__(self, *, name, index, type):
    self.name = str(name)
    self.index = int(index)
    self.type = str(type)
  
  @property
  def key(self):
    return f"{self.name[:-1]}.{self.index}"
  
  def to_string(self):
    return f"{self.name}.{self.index}-{self.type}"
  
  @classmethod
  def from_string(cls, content):
    content = str(content).split("-")
    type = content.pop()
    name, index = "-".join(content).split(".")
    return cls(name=name, index=index, type=type)
  
  def __str__(self):
    return f"{self.__class__.__name__}<{self.to_string()}>"
  
  def __repr__(self):
    return self.__str__()
  
  def __eq__(self, value):
    return isinstance(value, Residue) and self.to_string() == value.to_string()
  
  def __hash__(self):
    return hash(self.__str__())

class Atom:
  def __init__(self, *, atom, atom_type, x, y, z, residue):
    self.atom = str(atom)
    self.atom_type = str(atom_type)
    self.x = float(x)
    self.y = float(y)
    self.z = float(z)
    assert isinstance(residue, Residue) or residue is None
    self.residue = residue
  
  def to_string(self):
    return f"{self.atom} {self.x} {self.y} {self.z}"
  
  @classmethod
  def from_string(self, content):
    atom, x, y, z = str(content).strip().split()
    return Atom(atom=atom, atom_type="UNKNOWN", x=x, y=y, z=z, residue=None)
  
  def __str__(self):
    return f"{self.__class__.__name__}<{self.to_string()} : {self.atom_type} {self.residue}>"
  
  def __repr__(self):
    return self.__str__()
  
  def __eq__(self, value):
    return isinstance(value, Atom) and self.to_string() == value.to_string()
  
  def __hash__(self):
    return hash(self.__str__())

class Molecule:
  def __init__(self, *, name, base, variant = "", charge = 0, path):
    self.name = str(name)
    self.base = str(base)
    self.variant = str(variant)
    self.charge = int(charge)
    self.path = str(path)
  
  @property
  def util_name(self):
    return f"{self.base}({self.variant})"
  
  @property
  def names(self):
    if self.variant == "":
      return [self.name, self.util_name, self.base]
    return [self.name, self.util_name]
  
  def __str__(self):
    return f"{self.__class__.__name__}<{self.name}({self.base}.{self.variant}):{self.charge}>"
  
  def __repr__(self):
    return self.__str__()
  
  def __hash__(self):
    return hash(self.__str__())

# loading multiple residues

def load_residues(path):
  with open(path, "r") as file:
    content = file.read()
  
  lines = content.strip().split("\n")
  residues = []
  for line in lines:
    try:
      residues.append(Residue.from_string(line.split("#")[0]))
    except:
      pass
  
  return unify_residue_set(residues)

def load_residues_mini(content):
  residue_contents = str(content).split(",")
  residues = []
  for content in residue_contents:
    try:
      residues.append(Residue.from_string(content))
    except:
      residues.append(None)
  
  return unify_residue_set(residues)

def create_residue_map(residues: list[Residue]):
  residue_map = {}
  for residue in residues:
    if residue is None:
      continue
    residue_map[residue.key] = residue
  return residue_map

# loading/dumping xyz files

def load_xyz(path):
  with open(path, "r") as file:
    content = file.read()
  
  lines = content.strip().split("\n")[1:] # remove count, its useless to us
  comment = lines.pop(0)
  try:
    comment, atom_types, residues = comment.split("|")
    atom_types = atom_types.split(",")
    residues = load_residues_mini(residues)
  except:
    atom_types = []
    residues = []
  atoms = []
  for line in lines:
    atom = Atom.from_string(line)
    atom.atom_type = atom_types.pop(0) if len(atom_types) > 0 else "UNKNOWN"
    atom.residue = residues.pop(0) if len(residues) > 0 else None
    atoms.append(atom)
  
  return comment, unify_atom_residue_set(atoms)

def dump_xyz(path, atoms: list[Atom], comment=""):
  with open(path, "w") as file:
    file.write(f"{len(atoms)}\n")
    comment = f"{comment}|{','.join([atom.atom_type for atom in atoms])}|{','.join([atom.residue.to_string() for atom in atoms])}"
    file.write(f"{comment}\n")
    for atom in atoms:
      file.write(f"{atom.to_string()}\n")

# loading pdb files

def load_pdb(path):
  with open(path, "r") as file:
    content = file.read()
  
  lines = content.strip().split("\n")
  atoms = []
  for line in lines:
    line = line.strip()
    if not line.startswith("ATOM"):
      continue

    # type: 13-16
    atom_type = line[13:17].strip()
    # residue: 17-21
    residue = line[17:22].strip()
    # residue index: 22-31
    residue_index = line[22:32].strip()
    # x: 32-39
    x = line[32:40].strip()
    # y: 40-47
    y = line[40:48].strip()
    # z: 48-55
    z = line[48:56].strip()
    # atom: 77
    atom = line[77].strip()

    # never load water or salt
    if residue in ("WAT", "Na+", "Cl-"):
      continue

    atoms.append(Atom(atom=atom, atom_type=atom_type, x=x, y=y, z=z, residue=Residue(name=residue, index=residue_index, type="UNKNOWN")))
  
  return "", unify_atom_residue_set(atoms)

# loading molecule(s)

def load_molecule(path) -> Molecule:
  with open(path, "r") as file:
    content = file.read()
  args = [arg.strip() for arg in content.split(";") if len(arg.strip()) > 0]
  kv = [arg.split("=") for arg in args]
  args = dict(kv)

  name = os.path.basename(os.path.dirname(path))
  base = args["BASE"]
  variant = args.get("VARIANT", "")
  charge = int(args.get("CHARGE", 0))

  return Molecule(name=name, base=base, variant=variant, charge=charge, path=path)

def load_molecules(path):
  if os.path.isfile(path):
    def lookup(name) -> Molecule | None:
      return None
    def get() -> list[Molecule]:
      return []
    name = os.path.basename(path)
    if name != "mol":
      return lookup, get
    molecule = load_molecule(path)
    def lookup(name) -> Molecule | None:
      if name in molecule.names:
        return molecule
      return None
    def get() -> list[Molecule]:
      return [molecule]
    return lookup, get

  lookups = []
  gets = []
  for name in os.listdir(path):
    lookup, get = load_molecules(os.path.join(path, name))
    lookups.append(lookup)
    gets.append(get)
  
  def lookup(name) -> Molecule | None:
    for fn in lookups:
      mol = fn(name)
      if mol is not None:
        return mol
    return None
  def get() -> list[Molecule]:
    molecules = []
    for fn in gets:
      molecules.extend(fn())
    return molecules
  return lookup, get

# loading reactions

rxn_mol = tuple[int, str]
rxn_side = list[rxn_mol]
rxn_sides = tuple[rxn_side, rxn_side]
rxns = dict[str, rxn_sides]

def load_reaction_side(side) -> rxn_side:
  parts = [part.strip() for part in str(side).strip().split(" + ")]
  molecules = []
  for part in parts:
    part = re.split("\\s+", part)
    if len(part) > 2:
      continue
    if len(part) == 2:
      count = int(part[0])
      name = part[1]
    else:
      count = 1
      name = part[0]
    molecules.append((count, name))
  return molecules

def load_reaction_sides(reaction) -> rxn_sides:
  left, right = str(reaction).strip().split(" -> ")
  return load_reaction_side(left), load_reaction_side(right)

def load_reactions(path) -> rxns:
  with open(path, "r") as file:
    content = file.read()
  
  lines = [line.split("#")[0].strip() for line in content.split("\n")]
  lines = [line for line in lines if len(line) > 0]

  reactions = {}
  for reaction_name, reaction in zip(lines[0::2], lines[1::2]):
    reactions[reaction_name] = load_reaction_sides(reaction)
  
  return reactions

# dumping geom optimization files

OPT = """# {{comment}}

! UKS B3LYP ZORA def2-TZVP SARC/J TightSCF
! Normalprint Opt CPCMC(water)

%pal
        nprocs 16
end

%scf
        MaxIter 500
end

%geom
Constraints
{{constraints}}
end
end

* xyzfile {{charge}} {{multiplicity}} guess.xyz
"""

def dump_opt(path, atoms: list[Atom], constrained_atoms: list[Atom], comment="", charge=0, multiplicity=1):
  residues = list(create_residue_map([atom.residue for atom in atoms]).values())
  comment = f"{comment} | {','.join(residue.to_string() for residue in residues)}"
  constraints = "\n".join(["{C " + str(atoms.index(atom)) + " C}" for atom in constrained_atoms])
  charge = str(charge)
  multiplicity = str(multiplicity)
  content = OPT.replace("{{comment}}", comment).replace("{{constraints}}", constraints).replace("{{charge}}", charge).replace("{{multiplicity}}", multiplicity)
  with open(os.path.join(path, "Opt.inp"), "w") as file:
    file.write(content)

# dumping frequency analysis files

FREQ = """# {{comment}}

! B3LYP ZORA ma-ZORA-def2-TZVP SARC/J TightSCF SlowConv  MORead
! LargePrint Freq CPCMC(water)

%pal
        nprocs 16
end

%moinp "Opt.gbw"        # Read scratch from previous Opt results

%scf
        MaxIter 500
end

* xyzfile {{charge}} {{multiplicity}} Opt.xyz
"""

def dump_freq(path, comment="", charge=0, multiplicity=1):
  comment = f"{comment} | performed from geometry optimization"
  charge = str(charge)
  multiplicity = str(multiplicity)
  content = FREQ.replace("{{comment}}", comment).replace("{{charge}}", charge).replace("{{multiplicity}}", multiplicity)
  with open(os.path.join(path, "Freq.inp"), "w") as file:
    file.write(content)

# loading energies

def load_E(path):
  if not os.path.exists(path):
    return None
  if not os.path.isfile(path):
    return None
  # FINAL SINGLE POINT ENERGY       -79.849606724584
  with open(path, "r") as file:
    content = file.read()
  E = re.findall("FINAL SINGLE POINT ENERGY\\s+([-\\.0-9]+)", content)
  if len(E) <= 0:
    return None
  return float(E[-1]) * 627.5

def load_S(path):
  if not os.path.exists(path):
    return None
  if not os.path.isfile(path):
    return None
  # Final entropy term                ...      0.02691711 Eh     16.89 kcal/mol
  with open(path, "r") as file:
    content = file.read()
  S = re.findall("Final entropy term\\s+\\.\\.\\.\\s+([-\\.0-9]+)", content)
  if len(S) <= 0:
    return None
  return float(S[-1]) * 627.5

def load_H(path):
  if not os.path.exists(path):
    return None
  if not os.path.isfile(path):
    return None
  # Total Enthalpy                    ...    -79.77270507 Eh
  with open(path, "r") as file:
    content = file.read()
  H = re.findall("Total Enthalpy\\s+\\.\\.\\.\\s+([-\\.0-9]+)", content)
  if len(H) <= 0:
    return None
  return float(H[-1]) * 627.5

def load_G(path):
  if not os.path.exists(path):
    return None
  if not os.path.isfile(path):
    return None
  # Final Gibbs free energy         ...    -79.79962218 Eh
  with open(path, "r") as file:
    content = file.read()
  G = re.findall("Final Gibbs free energy\\s+\\.\\.\\.\\s+([-\\.0-9]+)", content)
  if len(G) <= 0:
    return None
  return float(G[-1]) * 627.5

def load_energies(path):
  opt_path = os.path.join(path, "Opt.out")
  freq_path = os.path.join(path, "Freq.out")
  E = load_E(opt_path)
  S = load_S(freq_path)
  H = load_H(freq_path)
  G = load_G(freq_path)
  return E, S, H, G

# get reaction energies

def get_reaction_mol_energies(lookup, energies: dict[Molecule, tuple[float, float, float, float]], reaction_mol: rxn_mol):
    count, name = reaction_mol
    mol = lookup(name)
    if mol not in energies:
      return 0, 0, 0, 0
    E, S, H, G = energies[mol]
    return E * count, S * count, H * count, G * count

def get_reaction_side_energies_array(lookup, energies: dict[Molecule, tuple[float, float, float, float]], reaction_side: rxn_side):
  return [get_reaction_mol_energies(lookup, energies, reaction_mol) for reaction_mol in reaction_side]

def get_reaction_side_energies(lookup, energies: dict[Molecule, tuple[float, float, float, float]], reaction_side: rxn_side):
  array = get_reaction_side_energies_array(lookup, energies, reaction_side)
  return tuple([sum(slice[i] for slice in array) for i in range(4)])

def get_reaction_delta_energies_array(lookup, energies: dict[Molecule, tuple[float, float, float, float]], reaction: rxn_sides):
  left, right = reaction
  return get_reaction_side_energies_array(lookup, energies, left), get_reaction_side_energies_array(lookup, energies, right)

def get_reaction_delta_energies(lookup, energies: dict[Molecule, tuple[float, float, float, float]], reaction: rxn_sides):
  left_array, right_array = get_reaction_delta_energies_array(lookup, energies, reaction)
  left = tuple([sum(slice[i] for slice in left_array) for i in range(4)])
  right = tuple([sum(slice[i] for slice in right_array) for i in range(4)])
  return tuple([right[i] - left[i] for i in range(4)])

# ensures that all residues with same key use the same residue object for reference

def unify_residue_set(residues_array: list[Residue]):
  residue_map = create_residue_map(residues_array)
  return [None if residue is None else residue_map[residue.key] for residue in residues_array]

def unify_atom_residue_set(atoms: list[Atom]):
  residues_array = unify_residue_set([atom.residue for atom in atoms])
  for atom, residue in zip(atoms, residues_array):
    atom.residue = residue
  return [*atoms]

# clearing a directory/file until all empty

def rmtree(path):
  if os.path.dirname(path) == path:
    return
  try:
    os.remove(path)
    rmtree(os.path.dirname(path))
  except:
    pass

# making a directory until all full

def mktree(path):
  if os.path.dirname(path) == path:
    return
  mktree(os.path.dirname(path))
  try:
    os.mkdir(path)
  except:
    pass

# vector math

vec = tuple[float, float, float]

def get(atom: Atom) -> vec:
  return (atom.x, atom.y, atom.z)

def set(atom: Atom, vec: vec):
  atom.x, atom.y, atom.z = vec
  return atom

def mag(vec: vec) -> float:
  x, y, z = vec
  return (x*x + y*y + z*z) ** 0.5

def add(a: vec, b: vec) -> vec:
  ax, ay, az = a
  bx, by, bz = b
  return (ax+bx, ay+by, az+bz)

def mul(a: vec, b: float) -> vec:
  x, y, z = a
  return (x*b, y*b, z*b)

def sub(a: vec, b: vec) -> vec:
  return add(a, mul(b, -1))

def div(a: vec, b: float) -> vec:
  return mul(a, 1/b)

def dist(a: vec, b: vec) -> float:
  return mag(sub(a, b))

def dot(a: vec, b: vec) -> float:
  ax, ay, az = a
  bx, by, bz = b
  return ax*bx + ay*by + az*bz

def cross(a: vec, b: vec) -> vec:
  ax, ay, az = a
  bx, by, bz = b
  return (ay*bz - az*by, az*bx - ax*bz, ax*by - ay*bx)

# console output

def log(s):
  print(f"  {s}")

def err(s):
  log(f"  ERROR: {s}")

def wrn(s):
  log(f"  WARNING: {s}")

def ttl(s):
  s = str(s)
  print()
  print(s)
  print("-" * len(s))
