import numpy as np
import matplotlib.pyplot as plt
import cfdmclasses

angle = 0

#creating objects
cycloid = cfdmclasses.Cycloid()

outputpinholes = cfdmclasses.OutputPinHoles()

rollerpins = cfdmclasses.RollerPins()

outputpins = cfdmclasses.OutputPins()


#plotting
plt.figure(figsize=(10,10))
for i in [cycloid, outputpinholes, rollerpins, outputpins]:
    i.set_pos(angle)
    i.plot()

plt.savefig("plot")