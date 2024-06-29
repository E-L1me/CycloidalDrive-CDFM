import numpy as np
import math
from dataclasses import dataclass

e = 0.8  # eccentricity distance
rp = 1.5  # roller pin radius
Rp = 30  # roller pin distribution radius
zc = 29  # number of cycloidal teeth
zp = zc + 1  # number of roller pins
rw = 3.5  # output pin radius
Rw = 21  # output pin distribution radius
rwh = rw + e  # output pin hole radius
Rwh = Rw  # output pin hole distribution radius
zw = 8  # number of output pins
# omegain = #input drive speed
# omegaout = #output drive speed


# if configuration is "pin gear fixed, output disk outputs":
#     Iew = omegain/omegaout = -zp/(zp - zc) = -zp #theoretical transmission ratio

# if configuration is "out put disk is fixed, pin gear outputs":
#     Iep = omegaout/omegain = zp/(zp - zc) = zp #theoretical transmission ratio


re = 8  # radius of input shaft
rec = 15  # baring radius of eccentric shaft in contact with cycloidal gear


de = 0  # eccentric shaft deviation

drwh = 0  # deviation in radius of pin-hole
dRwh = 0  # deviation in radius of distribution of pin-hole

drp = 0  # deviation in radius of roller pin
dRp = 0  # deviation in radius of distribution of roller pin


# translation matricies
def Rc2w(x):
    return np.array([[np.cos(x), -np.sin(x)], [np.sin(x), np.cos(x)]])


def Pc2w(x):
    return np.array([[-(e + de) * np.sin(x)], [(e + de) * np.cos(x)]])


def Rp2w(x):
    return np.array([[np.cos(x), -np.sin(x)], [np.sin(x), np.cos(x)]])


def get_axes(r, steps=100):
    xs = []
    for t in np.linspace(-r, r, steps):
        xs.append([t, 0])
    ys = []
    for t in np.linspace(-r, r, steps):
        ys.append([0, t])
    return xs, ys


def start_cycloid(t: float):
    K1 = (e * zp) / (Rp + dRp)
    x = (
        (Rp + dRp) * np.sin(t)
        - e * np.sin(zp * t)
        + ((K1 * np.sin(zp * t) - np.sin(t)) * (rp + drp))
        / (math.sqrt(1 + K1**2 - 2 * K1 * np.cos(zc * t)))
    )
    y = (
        (Rp + dRp) * np.cos(t)
        - e * np.cos(zp * t)
        + ((K1 * np.cos(zp * t) - np.cos(t)) * (rp + drp))
        / (math.sqrt(1 + K1**2 - 2 * K1 * np.cos(zc * t)))
    )
    x_dir = (x - zc * e * np.sin(zp * t)) / (
        math.sqrt(
            (x - zc * e * np.sin(zp * t)) ** 2 + (x - zc * e * np.cos(zp * t)) ** 2
        )
    )
    y_dir = (y - zc * e * np.cos(zp * t)) / (
        math.sqrt(
            (y - zc * e * np.sin(zp * t)) ** 2 + (y - zc * e * np.cos(zp * t)) ** 2
        )
    )
    K1 = (e * zp) / (Rp + dRp)
    pc = (1 + K1**2 - 2 * K1 * np.cos(zc * t)) ** 1.5 / (
        K1 * (1 + zp) * np.cos(zc * t - zp * K1**2 - 1)
    )
    return [t, x, y, x_dir, y_dir, pc]


def c_to_o(p: list, degrees: float):
    vector = np.array([[p[0]], [p[1]]])
    return np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)


def cycloid_to_output_rf(i: list, degrees: float):
    new = c_to_o(i[1:3], degrees)
    newnormvec = c_to_o(i[3:5], degrees)
    return [i[0], new[0, 0], new[1, 0], newnormvec[0, 0], newnormvec[0, 1]]


class Cycloid:
    num_points: int = 5000

    def __init__(self):
        self.tlist = np.linspace(
            0, 2 * np.pi, self.num_points
        )  # independent variable for the parametric equations
        self.basis = np.array(
            list(map(start_cycloid, self.tlist))
        )  # setting fundemental points for the cycloid with the sturucture: [points, normal vector, radius of curvature]
        xs, ys = get_axes(Rp)  # axes of the cycloid
        self.xaxis = np.array(xs)  # x axis of the cycloid
        self.yaxis = np.array(ys)  # y axis
        self.txaxis = self.xaxis  # position of x axis
        self.tyaxis = self.yaxis  # position of y axis
        self.pos = np.empty((-1, 5))  # position of the cycloid

    def set_pos(self, degrees):
        self.pos = map(cycloid_to_output_rf, self.basis, degrees)
        self.txaxis = map(c_to_o, self.xaxis, degrees)
        self.tyaxis = map(c_to_o, self.yaxis, degrees)
        return [self.pos, self.taxis, self.tyaxis]


class OutputPinHole:
    def __init__(self):
        self.npins = zw
        self.points = np.array([])  # points of the graph
        self.normdir = np.array([])  # normal direction for each point
        self.tlist = np.array([])  # independent variable for the parametric equations
        self.position = np.array([])  # position of the output pin holes
        self.pos_normdir = np.array([])  # norm dir for position of the output pin holes
        xs, ys = get_axes(Rp)  # axes of the cycloid
        self.xaxis = np.array(xs)  # x axis of the cycloid
        self.yaxis = np.array(ys)  # y axis
        self.txaxis = self.xaxis  # position of x axis
        self.tyaxis = self.yaxis  # position of y axis
        self.centers = np.array([])
        self.tcenters = np.array([])

    def set_tlist(self, steps=5000):
        self.tlist = np.linspace(0, 2 * np.pi, steps)
        return self.tlist

    def get_points(self):
        points = []
        centers = []
        for z in range(self.npins):
            pin = []
            i = z + 1
            centers.append(
                [
                    -(Rwh + dRwh) * np.sin(((2 * i + 1) * np.pi) / (zw)),
                    +(Rwh + dRwh) * np.cos(((2 * i + 1) * np.pi) / (zw)),
                ]
            )
            for t in self.tlist:
                x = -(rwh + drwh) * np.sin(t) - (Rwh + dRwh) * np.sin(
                    ((2 * i + 1) * np.pi) / (zw)
                )
                y = (rwh + drwh) * np.cos(t) + (Rwh + dRwh) * np.cos(
                    ((2 * i + 1) * np.pi) / (zw)
                )
                pin.append([x, y, t])
            points.append(pin)
        self.centers = np.array(centers)
        self.points = np.array(points)
        return self.points

    def get_normdir(self):
        for i in self.points:
            pin = np.array([])
            for t in i:
                x = -np.sin(t[1])
                y = np.cos(t[1])
                np.append(pin, (x, y, t))
            np.append(self.normdir, pin)
        return self.normdir

    def translate_points(self, degrees):
        position = []
        xs = []
        ys = []
        for z in self.points:
            pin = []
            for i in z:
                vector = np.array([[i[0]], [i[1]]])
                new = np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)
                pin.append([new[0, 0], new[1, 0], i[2]])
            position.append(pin)
        self.position = np.array(position)
        for i in self.xaxis:
            vector = np.array([[i[0]], [i[1]]])
            new = np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)
            xs.append([new[0, 0], new[1, 0]])
        self.txaxis = np.array(xs)
        for i in self.yaxis:
            vector = np.array([[i[0]], [i[1]]])
            new = np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)
            ys.append([new[0, 0], new[1, 0]])
        self.tyaxis = np.array(ys)
        centers = []
        for i in self.centers:
            vector = np.array([[i[0]], [i[1]]])
            new = np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)
            centers.append([new[0, 0], new[1, 0]])
        self.tcenters = np.array(centers)
        return self.position

    def translate_normdir(self, degrees):
        position = []
        for z in self.points:
            pin = []
            for i in z:
                vector = np.array([[i[0]], [i[1]]])
                new = np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)
                pin.append([new[0, 0], new[1, 0], i[2]])
            position.append(pin)
        self.pos_normdir = np.array(position)
        return self.pos_normdir


class RollerPin:
    def __init__(self):
        self.nrollers = zp
        self.points = np.array([])  # points of the graph
        self.normdir = np.array([])  # normal direction for each point
        self.tlist = np.array([])  # independent variable for the parametric equations
        self.position = np.array([])  # position of the roller pins
        self.pos_normdir = np.array([])  # norm dir for position of the roller pins
        xs, ys = get_axes(Rp)  # axes of the cycloid
        self.xaxis = np.array(xs)  # x axis of the cycloid
        self.yaxis = np.array(ys)  # y axis
        self.txaxis = self.xaxis  # position of x axis
        self.tyaxis = self.yaxis  # position of y axis
        self.centers = np.array([])
        self.tcenters = np.array([])

    def set_tlist(self, steps=5000):
        self.tlist = np.linspace(0, 2 * np.pi, steps)
        return self.tlist

    def get_points(self):
        points = []
        centers = []
        for z in range(self.nrollers):
            roller = []
            i = z + 1
            centers.append(
                [
                    -(Rp + dRp) * np.sin(2 * np.pi * (i - 1) / zp),
                    (Rp + dRp) * np.cos(2 * np.pi * (i - 1) / zp),
                ]
            )
            for t in self.tlist:
                x = -(rp + drp) * np.sin(t) - (Rp + dRp) * np.sin(
                    2 * np.pi * (i - 1) / zp
                )
                y = (rp + drp) * np.cos(t) + (Rp + dRp) * np.cos(
                    2 * np.pi * (i - 1) / zp
                )
                roller.append([x, y, t])
            points.append(roller)
        self.centers = np.array(centers)
        self.points = np.array(points)
        return self.points

    def get_normdir(self):
        for i in self.points:
            roller = np.array([])
            for t in i:
                x = -np.sin(t[1])
                y = np.cos(t[1])
                np.append(roller, (x, y, t))
            np.append(self.normdir, roller)
        return self.normdir

    def translate_points(self, degrees):
        position = []
        xs = []
        ys = []
        centers = []
        for z in self.points:
            pin = []
            for i in z:
                vector = np.array([[i[0]], [i[1]]])
                new = np.matmul(Rp2w(degrees), vector)
                pin.append([new[0, 0], new[1, 0], i[2]])
            position.append(pin)
        self.position = np.array(position)
        for i in self.xaxis:
            vector = np.array([[i[0]], [i[1]]])
            new = np.matmul(Rp2w(degrees), vector)
            xs.append([new[0, 0], new[1, 0]])
        self.txaxis = np.array(xs)
        for i in self.yaxis:
            vector = np.array([[i[0]], [i[1]]])
            new = np.matmul(Rp2w(degrees), vector)
            ys.append([new[0, 0], new[1, 0]])
        self.tyaxis = np.array(ys)
        for i in self.centers:
            vector = np.array([[i[0]], [i[1]]])
            new = np.matmul(Rp2w(degrees), vector)
            centers.append([new[0, 0], new[1, 0]])
        self.tcenters = np.array(centers)
        return self.position

    def translate_normdir(self, degrees):
        position = []
        for z in self.normdir:
            pin = []
            for i in z:
                vector = np.array([[i[0]], [i[1]]])
                new = np.matmul(Rp2w(degrees), vector)
                pin.append([new[0, 0], new[1, 0], i[2]])
            position.append(pin)
        self.pos_normdir = np.array(position)
        return self.pos_normdir


class OutputPin:
    def __init__(self):
        self.noutputs = zw
        self.points = np.array([])  # points of the graph
        self.normdir = np.array([])  # normal direction for each point
        self.tlist = np.array([])  # independent variable for the parametric equations
        xs, ys = get_axes(Rp)  # axes of the cycloid
        self.xaxis = np.array(xs)  # x axis of the cycloid
        self.yaxis = np.array(ys)  # y axis
        self.txaxis = self.xaxis  # position of x axis
        self.tyaxis = self.yaxis  # position of y axis
        self.centers = np.array([])

    def set_tlist(self, steps=5000):
        self.tlist = np.linspace(0, 2 * np.pi, steps)
        return self.tlist

    def get_points(self):
        points = []
        centers = []
        for z in range(self.noutputs):
            pin = []
            i = z + 1
            centers.append(
                [
                    -Rw * np.sin(((2 * i + 1) * np.pi) / (zw)),
                    Rw * np.cos(((2 * i + 1) * np.pi) / (zw)),
                ]
            )
            for t in self.tlist:
                x = -rw * np.sin(t) - Rw * np.sin(((2 * i + 1) * np.pi) / (zw))
                y = rw * np.cos(t) + Rw * np.cos(((2 * i + 1) * np.pi) / (zw))
                pin.append([x, y, t])
            points.append(pin)
        self.centers = np.array(centers)
        self.points = np.array(points)
        return self.points

    def get_normdir(self):
        for i in self.points:
            pin = np.array([])
            for t in i:
                x = -np.sin(t[1])
                y = np.cos(t[1])
                np.append(pin, (x, y, t))
            np.append(self.normdir, pin)
        return self.normdir

    def translate_points(self):
        pass


@dataclass
class Material:
    v: float  # Poisson's Ratio
    e: float  # Young's Modulus
    p: float  # Contact Curvature Radius


# Hertzian contact theory
def a(
    Fc, m_1: Material, m_2: Material, concavity
):  # half the contact deformation width
    p_star = 0
    if concavity == "Convex and concave":
        p_star = abs((m_1.p * m_2.p) / (m_1.p - m_2.p))
    if concavity == "Convex and convex":
        p_star = abs((m_1.p * m_2.p) / (m_1.p + m_2.p))
    return np.sqrt(
        (4 * Fc * p_star * (m_1.e * (1 - m_1.v**2) + m_2.e(1 - m_2.v**2)))
        / (np.pi * b * m_1.e * m_2.e)
    )


def delta(Fc, m_1: Material, m_2: Material, b, concavity):  # deformation
    a = a(Fc, m_1.p, m_2.p, m_1.e, m_2.e, m_1.v, m_2.v, b, concavity)
    return (
        (
            ((1 - m_1.v**2) / m_1.e) * (np.log((4 * abs(m_1.p)) / a) - 1 / 2)
            + ((1 - m_2.v**2) / m_2.e) * (np.log((4 * abs(m_2.p)) / a) - 1 / 2)
        )
        * (2 * Fc)
        / (np.pi * b)
    )


def sigma(Fc, m_1: Material, m_2: Material, b, concavity):  # contact stress
    a = a(Fc, m_1.p, m_2.p, m_1.e, m_2.e, m_1.v, m_2.v, b, concavity)
    return 2 * Fc / (np.pi * a * b)
