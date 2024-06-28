import numpy as np
import matplotlib.pyplot as plt
import cfdmclasses

angle = 50

#creating objects
mycycloid = cfdmclasses.Cycloid()
mycycloid.set_tlist()
mycycloid.get_points()
mycycloid.translate_points(angle)

myoutputpinholes = cfdmclasses.OutputPinHole()
myoutputpinholes.set_tlist()
myoutputpinholes.get_points()
myoutputpinholes.translate_points(angle)


myrollerpins = cfdmclasses.RollerPin()
myrollerpins.set_tlist()
myrollerpins.get_points()
myrollerpins.translate_points(angle)

myoutputpins = cfdmclasses.OutputPin()
myoutputpins.set_tlist()
myoutputpins.get_points()


#plotting
plt.figure(figsize=(10,10))
plt.plot(mycycloid.position[:,0], mycycloid.position[:,1], color="Blue")
for i in range(cfdmclasses.zw):
    plt.plot(myoutputpinholes.position[i,:,0],myoutputpinholes.position[i,:,1], color="Orange")
    plt.plot(myoutputpinholes.tcenters[i,0], myoutputpinholes.tcenters[i,1], color="Orange", marker="o")    
    plt.plot(myoutputpins.points[i,:,0],myoutputpins.points[i,:,1], color="Red")
    plt.plot(myoutputpins.centers[i,0], myoutputpins.centers[i,1], color="Red", marker="o")
for i in range(cfdmclasses.zp):
    plt.plot(myrollerpins.position[i,:,0],myrollerpins.position[i,:,1], color="Green")
    plt.plot(myrollerpins.tcenters[i,0],myrollerpins.tcenters[i,1], color="Green", marker="o")

plt.plot(mycycloid.txaxis[:,0], mycycloid.txaxis[:,1], color="Blue")
plt.plot(mycycloid.tyaxis[:,0], mycycloid.tyaxis[:,1], color="Blue")
plt.plot(myoutputpinholes.txaxis[:,0], myoutputpinholes.txaxis[:,1], color="Orange")
plt.plot(myoutputpinholes.tyaxis[:,0], myoutputpinholes.tyaxis[:,1], color="Orange")
plt.plot(myrollerpins.txaxis[:,0], myrollerpins.txaxis[:,1], color="Green")
plt.plot(myrollerpins.tyaxis[:,0], myrollerpins.tyaxis[:,1], color="Green")
plt.plot(myoutputpins.txaxis[:,0], myoutputpins.txaxis[:,1], color="Red")
plt.plot(myoutputpins.tyaxis[:,0], myoutputpins.tyaxis[:,1], color="Red")


plt.savefig("plot")