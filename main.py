from instance import InstanceRestrictedAssignment
from instances_template import pick_one_instance
from covering import Covering
import numpy as np


#n, m, M = pick_one_instance("1034/1020")
# M = np.array([[3, 6, 12, 24, 17, 33, 66, 132, 264, 528, 160, 640, 576, 320, 288],
#               [3, 6, 12, 24, 17, 33, 66, 132, 264, 528, 160, 640, 576, 320, 288]
#             ])
#
# n = M.shape[1]
# m = M.shape[0]
# instance = InstanceRestrictedAssignment(n, m, generate=False, M = M)
instance = Covering(3,15, [3, 6, 12, 24, 17, 33, 66, 132, 264, 528, 160, 640, 576, 320, 288])
sol, C_max = instance.opt_IP()
print(f"Optimal C_max: {C_max}, solutions {sol}")
C_max_LP = instance.opt_LP(verbose=True, C_max=C_max)
print(f"Optimal C_max LP relaxation: {C_max_LP[0]} , solutions {C_max_LP[0]}")

