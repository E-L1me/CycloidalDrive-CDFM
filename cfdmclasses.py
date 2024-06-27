import numpy as np
import math

e = .8 #eccentricity distance
rp = 1.5 #roller pin radius
Rp = 30 #roller pin distribution radius
zc = 29 #number of cycloidal teeth
zp = zc + 1 #number of roller pins
rw = 3.5 #output pin radius
Rw = 21 #output pin distribution radius
rwh = rw + e #output pin hole radius
Rwh = Rw #output pin hole distribution radius
zw = 8 #number of output pins
# omegain = #input drive speed
# omegaout = #output drive speed

"""
if configuration is "pin gear fixed, output disk outputs":
    Iew = omegain/omegaout = -zp/(zp - zc) = -zp #theoretical transmission ratio

if configuration is "out put disk is fixed, pin gear outputs":
    Iep = omegaout/omegain = zp/(zp - zc) = zp #theoretical transmission ratio
"""

re = 8 #radius of input shaft
rec = 15 #baring radius of eccentric shaft in contact with cycloidal gear


de = 0 #eccentric shaft deviation

drwh = 0 #deviation in radius of pin-hole
dRwh = 0 #deviation in radius of distribution of pin-hole

drp = 0 #deviation in radius of roller pin
dRp = 0 #deviation in radius of distribution of roller pin

#translation matricies
def Rc2w(x):
    return np.array([[np.cos(x), -np.sin(x)], [np.sin(x), np.cos(x)]])

def Pc2w(x):
    return np.array([[-(e + de) * np.sin(x)],[(e + de)*np.cos(x)]])

def Rp2w(x):
    return np.array([[np.cos(x), -np.sin(x)], [np.sin(x), np.cos(x)]])

def get_axes(r, steps=100):
    xs = []
    for t in np.linspace(-r, r,steps):
        xs.append([t, 0])
    ys = []
    for t in np.linspace(-r, r,steps):
        ys.append([0, t])
    return xs, ys

class Cycloid:
    def __init__(self):
        self.points = np.array([]) #points of the graph
        self.normdir = np.array([]) #normal direction for each point
        self.tlist = np.array([]) #independent variable for the parametric equations
        self.roc = np.array([]) #radius of curvature
        self.position = np.array([]) #position of the cycloid
        xs, ys = get_axes(Rp) #axes of the cycloid
        self.xaxis = np.array(xs) #x axis of the cycloid
        self.yaxis = np.array(ys) #y axis

    def set_tlist(self, steps=5000):
        self.tlist = np.linspace(0,2*np.pi,steps)
        return self.tlist

    def get_points(self):
        points = []
        for t in self.tlist:
            K1 = (e * zp)/(Rp + dRp)
            x = (Rp + dRp)* np.sin(t) - e * np.sin(zp*t) + ((K1*np.sin(zp*t)-np.sin(t))*(rp + drp))/(math.sqrt(1+K1**2 - 2*K1*np.cos(zc*t)))
            y = (Rp + dRp)* np.cos(t) - e * np.cos(zp*t) + ((K1*np.cos(zp*t)-np.cos(t))*(rp + drp))/(math.sqrt(1+K1**2 - 2*K1*np.cos(zc*t)))
            points.append([x,y,t])
        self.points = np.array(points)
        return self.points

    def get_normdir(self):
        for i in self.points:
            x = (i[0,0] - zc * e * np.sin(zp*i[1]))/(math.sqrt((i[0,0] - zc * e * np.sin(zp * i[1]))**2 + (i[0,0] - zc * e * np.cos(zp * i[1]))**2))
            y = (i[0,1] - zc * e * np.cos(zp*i[1]))/(math.sqrt((i[0,1] - zc * e * np.sin(zp * i[1]))**2 + (i[0,1] - zc * e * np.cos(zp * i[1]))**2))
            np.append(self.normdir, (x,y, i[2]))
        return self.normdir

    def get_roc(self):
        for t in self.tlist:
            K1 = (e * zp)/(Rp + dRp)
            pc = (1 + K1**2 -2*K1*np.cos(zc * t))**1.5/(K1*(1 + zp)*np.cos(zc*t - zp * K1**2 - 1))
            np.append(self.roc, (pc, t))
        return self.roc

    def translate_points(self, degrees):
        position = []
        for i in self.points:
            vector = np.array([[i[0]], [i[1]]])
            new = np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)
            position.append([new[0,0], new[1,0], i[2]])
        self.position = np.array(position)
        return self.position
        """  
        xs = []
        ys = []
        for i in self.xaxis:
            vector = np.array([[i[0]], [i[1]]])
            new = np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)
            xs.append([new[0,0], new[1,0], i[2]])
        for i in self.yaxis:
            vector = np.array([[i[0]], [i[1]]])
            new = np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)
            ys.append([new[0,0], new[1,0], i[2]])
        self.xaxis = np.array(xs)
        self.yaxis = np.array(ys)
        """


    def translate_normdir(self, degrees):
        position = []
        for i in self.points:
            vector = np.array([[i[0]], [i[1]]])
            new = np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)
            position.append([new[0,0], new[1,0], i[2]])
        self.normdir = np.array(position)
        return self.normdir

class OutputPinHole:
    def __init__(self):
        self.npins = zw
        self.points = np.array([]) #points of the graph
        self.normdir = np.array([]) #normal direction for each point
        self.tlist = np.array([]) #independent variable for the parametric equations
        self.position = np.array([]) #position of the output pin holes
        self.pos_normdir = np.array([]) #norm dir for position of the output pin holes

    def set_tlist(self, steps = 5000):
        self.tlist = np.linspace(0,2*np.pi,steps)
        return self.tlist

    def get_points(self):
        points = []
        for z in range(self.npins):
            pin = []
            for t in self.tlist:
                i = z + 1
                x =  -(rwh + drwh) * np.sin(t) - (Rwh + dRwh) * np.sin(((2*i+1)*np.pi)/(zw))
                y = (rwh + drwh) * np.cos(t) + (Rwh + dRwh) * np.cos(((2*i+1)*np.pi)/(zw))
                pin.append([x,y,t])
            points.append(pin)
        self.points = np.array(points)
        return self.points

    def get_normdir(self):
        for i in self.points:
            pin = np.array([])
            for t in i:
                x = - np.sin(t[1])
                y = np.cos(t[1])
                np.append(pin, (x,y,t))
            np.append(self.normdir, pin)
        return self.normdir

    def translate_points(self, degrees):
        position = []
        for z in self.points:
            pin = []
            for i in z:
                vector = np.array([[i[0]], [i[1]]])
                new = np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)
                pin.append([new[0,0], new[1,0], i[2]])
            position.append(pin)
        self.position = np.array(position)
        return self.position

    def translate_normdir(self, degrees):
        position = []
        for z in self.points:
            pin = []
            for i in z:
                vector = np.array([[i[0]], [i[1]]])
                new = np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)
                pin.append([new[0,0], new[1,0], i[2]])
            position.append(pin)
        self.pos_normdir = np.array(position)
        return self.pos_normdir


class RollerPin:
    def __init__(self):
        self.nrollers = zp
        self.points = np.array([]) #points of the graph
        self.normdir = np.array([]) #normal direction for each point
        self.tlist = np.array([]) #independent variable for the parametric equations
        self.position = np.array([]) #position of the roller pins
        self.pos_normdir = np.array([]) #norm dir for position of the roller pins

    def set_tlist(self, steps = 5000):
        self.tlist = np.linspace(0,2*np.pi,steps)
        return self.tlist

    def get_points(self):
        points = []
        for z in range(self.nrollers):
            roller = []
            i = z + 1
            for t in self.tlist:
                x = -(rp + drp) * np.sin(t) - (Rp + dRp) * np.sin(2 * np.pi * (i - 1)/zp)
                y = (rp + drp) * np.cos(t) + (Rp + dRp) * np.cos(2 * np.pi * (i - 1)/zp)
                roller.append([x,y,t])
            points.append(roller)
        self.points = np.array(points)
        return self.points

    def get_normdir(self):
        for i in self.points:
            roller = np.array([])
            for t in i:
                x = - np.sin(t[1])
                y = np.cos(t[1])
                np.append(roller, (x,y,t))
            np.append(self.normdir, roller)
        return self.normdir

    def translate_points(self, degrees):
        position = []
        for z in self.points:
            pin = []
            for i in z:
                vector = np.array([[i[0]], [i[1]]])
                new = np.matmul(Rp2w(degrees), vector)
                pin.append([new[0,0], new[1,0], i[2]])
            position.append(pin)
        self.position = np.array(position)
        return self.position

    def translate_normdir(self, degrees):
        position = []
        for z in self.normdir:
            pin = []
            for i in z:
                vector = np.array([[i[0]], [i[1]]])
                new = np.matmul(Rp2w(degrees), vector)
                pin.append([new[0,0], new[1,0], i[2]])
            position.append(pin)
        self.pos_normdir = np.array(position)
        return self.pos_normdir

class OutputPin:
    def __init__(self):
        self.noutputs = zw
        self.points = np.array([]) #points of the graph
        self.normdir = np.array([]) #normal direction for each point
        self.tlist = np.array([]) #independent variable for the parametric equations

    def set_tlist(self, steps = 5000):
        self.tlist = np.linspace(0,2*np.pi,steps)
        return self.tlist

    def get_points(self):
        points = []
        for z in range(self.noutputs):
            pin = []
            i = z + 1
            for t in self.tlist:
                x = -rw * np.sin(t) - Rw *np.sin(((2*i+1)*np.pi)/(zw))
                y = rw * np.cos(t) + Rw * np.cos(((2*i+1)*np.pi)/(zw))
                pin.append([x,y,t])
            points.append(pin)
        self.points = np.array(points)
        return self.points

    def get_normdir(self):
        for i in self.points:
            pin = np.array([])
            for t in i:
                x = - np.sin(t[1])
                y = np.cos(t[1])
                np.append(pin, (x,y,t))
            np.append(self.normdir, pin)
        return self.normdir


    def translate_points(self):
        pass