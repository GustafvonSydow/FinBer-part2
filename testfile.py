import numpy as np
lst = np.linspace(0, 0.5, 10000)

the_s = 0.27



prev_val = 0
for val in lst:
    if val < the_s:
        prev_val = val
    elif val == the_s:
        index = lst.index(val)
    elif np.abs(val - the_s) <= np.abs(prev_val - the_s):
        index = 

print(ind)