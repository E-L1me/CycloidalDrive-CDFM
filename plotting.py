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

outputpins.plot()

plt.grid()

plt.savefig("plot")

#a few ideas on why the graphing could be causing problems: 
# - linespace arguments are wrong?
# - other arguments are wrong
# - one of the intermediary functions crush something to a smaller size
# - not all points are passed through?
# - pyplot can't handle all the points?
# - something large is divinding everything
# - the limiting/shrinking is happening in both the graphs of the circles and the graphs of the axes, so somthing fundemental is wrong?
# - so linespace problem?