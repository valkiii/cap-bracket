"""Chrome tube cap with S-bent mounting flange - parametric draft (all mm, ESTIMATED from photos).
Cylinder axis = Z, open end at z=0, flange on the +Y side. Profile (YZ): sloped top on the far
half, flat ledge on the near half, then up a short outer wall, jog inward, up a tall plate."""
import numpy as np, trimesh
from trimesh.transformations import rotation_matrix as rot

OD, ID, H = 18.5, 15.5, 14.0      # outer dia, bore dia, height at the high (far) edge
LEDGE_DROP = 2.7                  # flat ledge level is this far below the high edge
SLOPE_END = -1.3                  # y where the slope meets the flat ledge
CEIL = 2.0                        # closed-top thickness under the ledge level
TOP_HOLE = 3.4                    # vertical hole through the top (at the centre)
BOSS_OD, BOSS_LEN = 5.6, 2.0      # boss around that hole, inside the bore
SIDE_HOLE, SIDE_Z = 3.4, 7.0      # cross hole through the +Y wall, height above the open end
LEDGE_T = 1.3                     # ledge sheet thickness (rectangular, OD wide, flush with the wall edge)
WALL_T, WALL_H = 1.5, 5.2         # outer step wall thickness, height above the ledge top
JOG_T, JOG_IN = 1.5, 6.5          # horizontal jog thickness, how far it runs inward from the edge
PLATE_T, PLATE_TOP = 2.5, 16.0    # tall plate thickness, its top above the ledge top
PLATE_IN = 4.3                    # plate outer face distance inward from the edge (plate = 4.3..6.8 in)
FL_HOLE, FL_HOLE_DROP = 3.8, 4.0  # plate hole dia, distance below the plate top
R = OD / 2; S = 96
zl = H - LEDGE_DROP               # ledge top level

def cyl(d, h, z0=0):
    m = trimesh.creation.cylinder(radius=d/2, height=h, sections=S); m.apply_translation([0,0,z0+h/2]); return m
def box(x0,x1,y0,y1,z0,z1):
    m = trimesh.creation.box(extents=[x1-x0,y1-y0,z1-z0]); m.apply_translation([(x0+x1)/2,(y0+y1)/2,(z0+z1)/2]); return m
def u(*m): return trimesh.boolean.union(list(m), engine="manifold")
def d(a,*m): return trimesh.boolean.difference([a]+list(m), engine="manifold")
def i(a,b): return trimesh.boolean.intersection([a,b], engine="manifold")

# body: cylinder trimmed by the sloped-then-flat top profile (extruded along X)
def hull(pts):  # pts are (y,z); extruded across x
    P=[(x,y,z) for (y,z) in pts for x in (-R-1,R+1)]
    return trimesh.convex.convex_hull(np.array(P))
prof = u(hull([(-R-1,-1),(R+1,-1),(R+1,zl),(SLOPE_END,zl)]+[(-R-1,-1)]),
         hull([(-R-1,-1),(SLOPE_END,-1),(SLOPE_END,zl),(-R,H),(-R-1,H)]))
body = i(cyl(OD, H+1, z0=-0.5), prof)

# S-bend flange (all OD wide in X)
ledge = box(-R,R, SLOPE_END, R, zl-LEDGE_T, zl)
wall  = box(-R,R, R-0.3-WALL_T, R-0.3, zl, zl+WALL_H)
jog   = box(-R,R, R-JOG_IN-0.3, R-0.3, zl+WALL_H-JOG_T, zl+WALL_H)
plate = box(-R,R, R-PLATE_IN-PLATE_T, R-PLATE_IN, zl+WALL_H-JOG_T, zl+PLATE_TOP)
body = u(body, ledge, wall, jog, plate)

# bore, boss, holes
body = d(body, cyl(ID, zl - CEIL))
body = u(body, cyl(BOSS_OD, BOSS_LEN, z0=zl-CEIL-BOSS_LEN+0.01))
body = d(body, cyl(TOP_HOLE, 60, z0=-1))
def ycyl(dia, length, y0, z):
    m = trimesh.creation.cylinder(radius=dia/2, height=length, sections=S)
    m.apply_transform(rot(np.pi/2,[1,0,0])); m.apply_translation([0,y0+length/2,z]); return m
body = d(body, ycyl(SIDE_HOLE, 10, R-5, SIDE_Z))
body = d(body, ycyl(FL_HOLE, 10, R-PLATE_IN-PLATE_T-2, zl+PLATE_TOP-FL_HOLE_DROP))

body.export("cap_bracket.stl")
print("watertight:", body.is_watertight, "bbox mm:", np.round(body.bounds[1]-body.bounds[0],1), "vol mm3: %.0f" % body.volume)

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig = plt.figure(figsize=(15,5))
for k,(el,az,t) in enumerate([(20,-60,"3D"),(0,0,"front (flange side)"),(0,90,"side profile (flange right)")]):
    ax = fig.add_subplot(1,3,k+1,projection="3d")
    ax.plot_trisurf(*body.vertices.T[:2], body.vertices[:,2], triangles=body.faces, color="#b8c4cc", edgecolor="none", shade=True)
    ax.view_init(el,az); ax.set_title(t); ax.set_box_aspect(np.ptp(body.vertices,axis=0)); ax.set_axis_off()
plt.savefig("preview.png", dpi=90, bbox_inches="tight")
