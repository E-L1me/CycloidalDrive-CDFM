import numpy as np
import matplotlib.pyplot as plt
import cycloidgeometry

angle = 0

#creating objects
cycloid = cycloidgeometry.Cycloid()

outputpinholes = cycloidgeometry.OutputPinHoles()

rollerpins = cycloidgeometry.RollerPins()

outputpins = cycloidgeometry.OutputPins()


#plotting
plt.figure(figsize=(10,10))
for i in [cycloid, outputpinholes, rollerpins]:
    i.set_pos(angle)
    i.plot()

outputpins.plot

plt.savefig("plot")