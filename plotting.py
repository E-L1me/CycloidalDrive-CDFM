import numpy as np
import matplotlib.pyplot as plt
import cfdmclasses


#creating objects
mycycloid = cfdmclasses.Cycloid()
mycycloid.set_tlist()
mycycloid.get_points()
mycycloid.translate_points(15)

myoutputpinholes = cfdmclasses.OutputPinHole()
myoutputpinholes.set_tlist()
myoutputpinholes.get_points()
myoutputpinholes.translate_points(15)


myrollerpins = cfdmclasses.RollerPin()
myrollerpins.set_tlist()
myrollerpins.get_points()
myrollerpins.translate_points(15)

myoutputpins = cfdmclasses.OutputPin()
myoutputpins.set_tlist()
myoutputpins.get_points()


#plotting
plt.figure(figsize=(10,10))
plt.plot(mycycloid.position[:,0], mycycloid.position[:,1], color="Blue")
for i in range(cfdmclasses.zw):
    plt.plot(myoutputpinholes.position[i,:,0],myoutputpinholes.position[i,:,1], color="Orange")
    plt.plot(myoutputpins.points[i,:,0],myoutputpins.points[i,:,1], color="Red")
for i in range(cfdmclasses.zp):
    plt.plot(myrollerpins.position[i,:,0],myrollerpins.position[i,:,1], color="Green")

plt.plot(mycycloid.xaxis[:,0], mycycloid.xaxis[:,1], color="Blue")
plt.plot(mycycloid.yaxis[:,0], mycycloid.yaxis[:,1], color="Blue")

plt.savefig("plot")