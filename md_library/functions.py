import matplotlib as mpl
import numpy as np

def ColorFade(c1, c2, mix, temp, T_c):
    c1 = np.array(mpl.colors.to_rgb(c1))
    c2 = np.array(mpl.colors.to_rgb(c2))
    if T_c and temp == int(T_c):
        return "black"
    else:
        return mpl.colors.to_hex((1-mix) * c1 + mix * c2)

def GetLabel(temp, T_c):
    if T_c is not None and temp == int(T_c):
        return r"$T_\mathrm{C} = $" + str(temp)
    else:
        return str(temp)
