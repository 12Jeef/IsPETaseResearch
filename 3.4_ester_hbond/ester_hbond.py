import os
import matplotlib.pyplot as plt

root = os.path.dirname(os.path.abspath(__file__))

X = int(input("X="))

start, stop, n = [
  (2, 4, 11),
  (2.8, 4.2, 15),
  (3, 5, 21),
  (3, 5, 21),
  (2, 4, 21),
][X-1]

with open(os.path.join(root, "out", "ester_hbond" + ("" if X == 1 else f"_{X}"), "Opt.out"), "r") as file:
  content = file.read()

energies = []

last_was_final = False
while True:
  index = content.find("FINAL")
  if index < 0:
    break
  content = content[index+len("FINAL"):].strip()
  if last_was_final:
    energy = float(content.split("\n")[0][len("SINGLE POINT ENERGY"):].strip()) * 627.5
    energies.append(energy)
  last_was_final = content.startswith("ENERGY")

print(f"Found {len(energies)} scan energies")
for E in energies:
  print(f"{E}")

with open(os.path.join(root, "out", "methylamine", "Opt.out"), "r") as file:
  content = file.read()

methylamine_E = 0
lines = content.split("\n")
for line in lines:
  if line.startswith("FINAL SINGLE POINT ENERGY"):
    methylamine_E = float(line[len("FINAL SINGLE POINT ENERGY"):].strip()) * 627.5

if X > 2:
  with open(os.path.join(root, "out", "ester_ion", "Opt.out"), "r") as file:
    content = file.read()
else:
  with open(os.path.join(root, "../energies.csv"), "r") as file: # TODO: fix when reorg
    content = file.read()

ethyl_benzoate_E = 0
if X > 2:
  lines = content.split("\n")
  for line in lines:
    if line.startswith("FINAL SINGLE POINT ENERGY"):
      ethyl_benzoate_E = float(line[len("FINAL SINGLE POINT ENERGY"):].strip()) * 627.5
else:
  lines = content.split("\n")
  for line in lines:
    if line.startswith("ethyl_benzoate()"):
      ethyl_benzoate_E = float(line.split(",")[1].strip())

print(f"Methylamine E:   \t{methylamine_E}")
print(f"Ethyl Benzoate E:\t{ethyl_benzoate_E}")

zero_E = methylamine_E + ethyl_benzoate_E
print(f"Zero E:          \t{zero_E}")

print(f"Test E:          \t{(-595.672279358055 * 627.5) - zero_E}")

data = [E - zero_E for E in energies]

print(f"Data:")
for i, E in enumerate(data):
  print(f"{start + (stop-start) * i / (n-1)} \t {E}")

plt.plot([start + (stop-start) * i / (n-1) for i in range(n)], data)
plt.xlabel("Distance (Å)")
plt.ylabel("Energy (kcal/mol)")
plt.show()
