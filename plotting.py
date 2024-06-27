import numpy as np
import matplotlib.pyplot as plt
import math
import cfdmclasses


#creating objects
mycycloid = cfdmclasses.Cycloid()
mycycloid.set_tlist()
mycycloid.get_points()

myoutputpinholes = cfdmclasses.OutputPinHole()
myoutputpinholes.set_tlist()
myoutputpinholes.get_points()

myrollerpins = cfdmclasses.RollerPin()
myrollerpins.set_tlist()
myrollerpins.get_points()

myoutputpins = cfdmclasses.OutputPin()
myoutputpins.set_tlist()
myoutputpins.get_points()

myrollerpins = cfdmclasses.RollerPin()
myrollerpins.set_tlist()
myrollerpins.get_points()

#plotting
plt.figure(figsize=(10,10))
plt.plot(mycycloid.points[:,0], mycycloid.points[:,1])
for i in range(cfdmclasses.zw):
    plt.plot(myoutputpinholes.points[i,:,0],myoutputpinholes.points[i,:,1], color="Orange")
    plt.plot(myoutputpins.points[i,:,0],myoutputpins.points[i,:,1], color="Red")
for i in range(cfdmclasses.zp):
    plt.plot(myrollerpins.points[i,:,0],myrollerpins.points[i,:,1], color="Green")

plt.savefig("plot")