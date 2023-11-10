import os
from lrbound.qa_lr import *

r = 2 * 3 + 2
d = 3
T = 3.33
alpha = 1.53

mini1 = 1
mini2 = 1
mini3 = 1
ratio = 0.7020
path1 = ""
path2 = ""

lr_bd = lr_alpha(r, d, T, alpha)
p = 3

"""======================"""

if p < 3:
    step = 10
else:
    step = 1000
for i, path in enumerate(os.listdir(f"balls/{p}/")):
    num = int(path.split("_")[0][6:])
    if i % step == 0:
        print("STEP", i, flush=True)

    g = nx.read_adjlist(f"balls/{p}/{path}")
    energy = edge_energy(g, T, alpha)

    if p == 3:
        if num > 17:
            omega = 3
            num = num - 17
        elif num < 4:
            omega = 1
        else:
            omega = 2
            num = num - 3
    else:
        omega = num
    if omega == 1:
        if energy - lr_bd < mini1:
            mini1 = (energy - lr_bd)[0, 0].real
            path1 = path

    elif omega == 2:
        if energy - lr_bd < mini2:
            mini2 = (energy - lr_bd)[0, 0].real
            path2 = path

    elif omega == 3:
        if energy - lr_bd < ratio:
            print(path, "has energy", (energy - lr_bd)[0, 0].real, flush=True)

            loc_bd = loc_lr_p3(g, T, alpha)
            assert (energy - loc_bd)[
                0, 0
            ].real >= ratio, f"Ball {path} fails to reach the target ratio"

            print("-> with local bound: ", (energy - loc_bd)[0, 0].real, flush=True)
            if (energy - loc_bd)[0, 0].real < mini3:
                mini3 = (energy - loc_bd)[0, 0].real
                path3 = path


print("=" * 90)
print("worst omega 1 :", path1, " with value =", mini1, flush=True)
print("worst omega 2 :", path2, " with value =", mini2, flush=True)
print("worst omega 3 :", path3, " with value =", mini3, flush=True)
print("=" * 90)
