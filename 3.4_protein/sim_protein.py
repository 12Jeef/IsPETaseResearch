from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax

init("""
-relax:constrain_relax_to_start_coords
-relax:ramp_constraints false
""")

pose = pose_from_file("af3_model.cif")

scorefxn = create_score_function("ref2015_cart")

relax = FastRelax()
relax.set_scorefxn(scorefxn)
relax.cartesian(True)

relax.apply(pose)

pose.dump_pdb("relaxed.pdb")