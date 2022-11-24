from pybinding.repository import graphene
import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
def rectangle(width, height):
    x0 = width / 2
    y0 = height / 2
    return pb.Polygon([[x0, y0], [x0, -y0], [-x0, -y0], [-x0, y0]])

shape = rectangle(1.6, 1.2)
shape.plot()


model = pb.Model(
    graphene.monolayer(),
     rectangle(width=1.6, height=1.2)
)
fig, ax = plt.subplots()
ax = model.plot()

print(np.array(model.lattice))
#model.lattice.plot_vectors(position=[0.6, -0.25])