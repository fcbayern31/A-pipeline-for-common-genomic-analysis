import sys
import straw
import numpy as np
import pandas as pd 

IN=sys.argv[1]

mat = straw.straw("KR", IN, "as", "as", "BP", rslu=500000)
x = np.hstack([mat[0], mat[1]])
y = np.hstack([mat[1], mat[0]])

lb, xb, yb  = np.log(np.hstack([mat[2], mat[2]])) , np.hstack([mat[0], mat[1]]) , np.hstack([mat[1], mat[0]])

bm = pd.DataFrame({"x": xb, "y": yb, "counts": lb},columns=["x", 'y', 'counts'])
fm = bm.pivot("x", "y", "counts").fillna(0)
fm.to_csv("draw.matrix", sep='\t')
