import numpy as np
import math
from dataclasses import dataclass
import matplotlib.pyplot as plt
import json

f = open('constants.json')

e, rp, Rp, zc, rw, Rw, zw, re, rec, de, drwh, dRwh, drp, dRp = json.load(f)


e = int(e)  # eccentricity distance
rp = int(rp)  # roller pin radius
Rp = int(Rp)  # roller pin distribution radius
zc = int(zc)  # number of cycloidal teeth
zp = zc + 1  # number of roller pins
rw = int(rw)  # output pin radius
Rw = int(Rw)  # output pin distribution radius
rwh = rw + e  # output pin hole radius
Rwh = Rw # output pin hole distribution radius
zw = int(zw)  # number of output pins
# omegain = #input drive speed
# omegaout = #output drive speed


# if configuration is "pin gear fixed, output disk outputs":
#     Iew = omegain/omegaout = -zp/(zp - zc) = -zp #theoretical transmission ratio

# if configuration is "out put disk is fixed, pin gear outputs":
#     Iep = omegaout/omegain = zp/(zp - zc) = zp #theoretical transmission ratio


re = int(re)  # radius of input shaft
rec = int(rec)  # baring radius of eccentric shaft in contact with cycloidal gear


de = int(de)  # eccentric shaft deviation

drwh = int(drwh)  # deviation in radius of pin-hole
dRwh = int(dRwh)  # deviation in radius of distribution of pin-hole

drp = int(drp)  # deviation in radius of roller pin
dRp = int(dRp)  # deviation in radius of distribution of roller pin


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


def c_to_o(p: list, degrees: float):
    vector = np.array([[p[0]], [p[1]]])
    ret = np.matmul(Rc2w(degrees), vector) + Pc2w(degrees)
    return ret.flatten()


def r_to_o(p: list, degrees: float):
    vector = np.array([[p[0]], [p[1]]])
    ret = np.matmul(Rp2w(degrees), vector)
    return ret.flatten()


def cycloid_to_output_rf(i: list, degrees: float):
    new = c_to_o(i[1:3], degrees)
    newnormvec = c_to_o(i[3:5], degrees)
    return [i[0], new[0], new[1], newnormvec[0], newnormvec[1]]


def rollerpins_to_output_rf(i: list, degrees: float):
    new = r_to_o(i[1:3], degrees)
    newnormvec = r_to_o(i[3:5], degrees)
    return [i[0], new[0], new[1], newnormvec[0], newnormvec[1]]

def start_cycloid(t: float):
    K1 = (e * zp) / (Rp + dRp)
    x = (
        (Rp + dRp) * np.sin(t)
        - e * np.sin(zp * t)
        + ((K1 * np.sin(zp * t) - np.sin(t)) * (rp + drp)) / (math.sqrt(1 + K1**2 - 2 * K1 * np.cos(zc * t)))
    )
    y = (
        (Rp + dRp) * np.cos(t)
        - e * np.cos(zp * t)
        + ((K1 * np.cos(zp * t) - np.cos(t)) * (rp + drp)) / (math.sqrt(1 + K1**2 - 2 * K1 * np.cos(zc * t)))
    )
    x_dir = (x - zc * e * np.sin(zp * t)) / (
        math.sqrt((x - zc * e * np.sin(zp * t)) ** 2 + (x - zc * e * np.cos(zp * t)) ** 2)
    )
    y_dir = (y - zc * e * np.cos(zp * t)) / (
        math.sqrt((y - zc * e * np.sin(zp * t)) ** 2 + (y - zc * e * np.cos(zp * t)) ** 2)
    )
    K1 = (e * zp) / (Rp + dRp)
    pc = (1 + K1**2 - 2 * K1 * np.cos(zc * t)) ** 1.5 / (K1 * (1 + zp) * np.cos(zc * t - zp * K1**2 - 1))
    return [t, x, y, x_dir, y_dir, pc]


def start_outputpinhole(t: float, i: int):
    x = -(rwh + drwh) * np.sin(t) - (Rwh + dRwh) * np.sin(((2 * i + 1) * np.pi) / (zw))
    y = (rwh + drwh) * np.cos(t) + (Rwh + dRwh) * np.cos(((2 * i + 1) * np.pi) / (zw))
    x_dir = -np.sin(t)
    y_dir = np.cos(t)
    return [t, x, y, x_dir, y_dir]


def start_rollerpin(t: float, i: int):
    x = -(rp + drp) * np.sin(t) - (Rp + dRp) * np.sin(2 * np.pi * (i - 1) / zp)
    y = (rp + drp) * np.cos(t) + (Rp + dRp) * np.cos(2 * np.pi * (i - 1) / zp)
    x_dir = -np.sin(t)
    y_dir = np.cos(t)
    return [t, x, y, x_dir, y_dir]


def start_outputpin(t: float, i: int):
    x = -rw * np.sin(t) - Rw * np.sin(((2 * i + 1) * np.pi) / (zw))
    y = rw * np.cos(t) + Rw * np.cos(((2 * i + 1) * np.pi) / (zw))
    x_dir = -np.sin(t)
    y_dir = np.cos(t)
    return [t, x, y, x_dir, y_dir]


class Cycloid:
    num_points: int = 5000
    color: str = "Blue"

    def __init__(self):
        self.tlist: np.ndarray = np.linspace(0, 2 * np.pi, self.num_points)  # independent variable for the parametric equations
        self.basis: np.ndarray = np.array(
            list(map(start_cycloid, self.tlist))
        )  # setting fundemental points for the cycloid with the sturucture: [points, normal vector, radius of curvature]
        xs, ys = get_axes(Rp)  # axes of the cycloid
        self.xaxis: np.ndarray = np.array(xs)  # x axis of the cycloid
        self.yaxis: np.ndarray = np.array(ys)  # y axis
        self.txaxis: np.ndarray = self.xaxis  # position of x axis
        self.tyaxis: np.ndarray = self.yaxis  # position of y axis
        self.pos: np.ndarray = np.empty(2)  # position of the cycloid

    def set_pos(self, degrees):
        d = [degrees for j in range(self.basis.shape[0])]
        self.pos = np.array(list(map(cycloid_to_output_rf, self.basis, d)))
        self.txaxis = np.array(list(map(c_to_o, self.xaxis, d)))
        self.tyaxis = np.array(list(map(c_to_o, self.yaxis, d)))
        return [self.pos, self.txaxis, self.tyaxis]

    def plot(self):
        plt.plot(self.pos[:, 1], self.pos[:, 2], color=self.color)
        plt.plot(self.txaxis[:, 0], self.txaxis[:, 1], color=self.color)
        plt.plot(self.tyaxis[:, 0], self.tyaxis[:, 1], color=self.color)


class OutputPinHoles:
    num_points: int = 5000
    color: str = "Blue"

    def __init__(self):
        self.npins: int = zw
        self.tlist: np.ndarray = np.linspace(
            0, 2 * np.pi, self.num_points
        )  # independent variable for the parametric equations
        self.basis = np.array(
            [
                list(
                    map(
                        start_outputpinhole,
                        self.tlist,
                        [i for j in range(self.num_points)],
                    )
                )
                for i in range(1, self.npins + 1)
            ]
        )  # setting fundemental points for the cycloid with the sturucture: [points, normal vector, radius of curvature]
        xs, ys = get_axes(Rp)  # axes of the cycloid
        self.xaxis = np.array(xs)  # x axis of the cycloid
        self.yaxis = np.array(ys)  # y axis
        self.txaxis = self.xaxis  # position of x axis
        self.tyaxis = self.yaxis  # position of y axis
        self.pos = np.empty(2)  # position of the cycloid
        self.centers = np.array(
            [
                [
                    -(Rwh + dRwh) * np.sin(((2 * i + 1) * np.pi) / (zw)),
                    (Rwh + dRwh) * np.cos(((2 * i + 1) * np.pi) / (zw)),
                ]
                for i in range(1, self.npins + 1)
            ]
        )
        self.tcenters = self.centers

    def set_pos(self, degrees):
        d = [degrees for j in range(self.num_points)]
        self.pos = np.array([np.array(list(map(cycloid_to_output_rf, hole, d))) for hole in self.basis])
        self.txaxis = np.array(list(map(c_to_o, self.xaxis, d)))
        self.tyaxis = np.array(list(map(c_to_o, self.yaxis, d)))
        self.tcenters = np.array(list(map(c_to_o, self.centers, d)))
        return [self.pos, self.txaxis, self.tyaxis]

    def plot(self):
        for i in range(self.npins):
            plt.plot(self.pos[i, :, 1], self.pos[i, :, 2], color=self.color)
        plt.plot(self.txaxis[:, 0], self.txaxis[:, 1], color=self.color)
        plt.plot(self.tyaxis[:, 0], self.tyaxis[:, 1], color=self.color)
        plt.plot(self.tcenters[:, 0], self.tcenters[:, 1], "o", color=self.color)


class RollerPins:
    num_points: int = 5000
    color: str = "Green"

    def __init__(self):
        self.npins: int = Rp
        self.tlist: np.ndarray = np.linspace(
            0, 2 * np.pi, self.num_points
        )  # independent variable for the parametric equations
        self.basis = np.array(
            [
                    list(
                        map(
                            start_rollerpin,
                            self.tlist,
                            [i for j in range(self.tlist.shape[0])],
                        )
                    )
                for i in range(1, self.npins + 1)
            ]
        )  # setting fundemental points for the cycloid with the sturucture: [points, normal vector, radius of curvature]
        xs, ys = get_axes(Rp)  # axes of the cycloid
        self.xaxis = np.array(xs)  # x axis of the cycloid
        self.yaxis = np.array(ys)  # y axis
        self.txaxis = self.xaxis  # position of x axis
        self.tyaxis = self.yaxis  # position of y axis
        self.pos = np.empty(2)  # position of the cycloid
        self.centers = np.array(
            [
                [
                    -(Rp + dRp) * np.sin(2 * np.pi * (i - 1) / zp),
                    (Rp + dRp) * np.cos(2 * np.pi * (i - 1) / zp),
                ]
                for i in range(1, self.npins + 1)
            ]
        )
        self.tcenters = self.centers

    def set_pos(self, degrees):
        d = [degrees for j in range(self.num_points)]
        self.pos = np.array([np.array(list(map(rollerpins_to_output_rf, pin, d))) for pin in self.basis])
        self.txaxis = np.array(list(map(r_to_o, self.xaxis, d)))
        self.tyaxis = np.array(list(map(r_to_o, self.yaxis, d)))
        self.tcenters = np.array(list(map(r_to_o, self.centers, d)))
        return [self.pos, self.txaxis, self.tyaxis]

    def plot(self):
        for i in range(self.npins):
            plt.plot(self.pos[i, :, 1], self.pos[i, :, 2], color=self.color)
        plt.plot(self.txaxis[:, 0], self.txaxis[:, 1], color=self.color)
        plt.plot(self.tyaxis[:, 0], self.tyaxis[:, 1], color=self.color)
        plt.plot(self.tcenters[:, 0], self.tcenters[:, 1], "o", color=self.color)


class OutputPins:
    num_points: int = 5000
    color: str = "Black"

    def __init__(self):
        self.npins: int = zw
        self.tlist: np.ndarray = np.linspace(
            0, 2 * np.pi, self.num_points
        )  # independent variable for the parametric equations
        self.basis = np.array(
                    [list(
                        map(
                            start_outputpin,
                            self.tlist,
                            [i for j in range(self.num_points)],
                        )
                    )
                    for i in range(1, self.npins + 1)]
        )  # setting fundemental points for the cycloid with the sturucture: [points, normal vector, radius of curvature]
        xs, ys = get_axes(Rp)  # axes of the cycloid
        self.xaxis = np.array(xs)  # x axis of the cycloid
        self.yaxis = np.array(ys)  # y axis
        self.txaxis = self.xaxis
        self.tyaxis = self.yaxis
        self.centers = np.array(
            [
                [
                    -Rw * np.sin(((2 * i + 1) * np.pi) / (zw)),
                    Rw * np.cos(((2 * i + 1) * np.pi) / (zw)),
                ]
                for i in range(1, self.npins + 1)
            ]
        )

    def plot(self):
        for i in range(self.npins):
            plt.plot(self.basis[i-1,:, 1], self.basis[i-1, :, 2], color=self.color)
        plt.plot(self.xaxis[:, 0], self.xaxis[:, 1], color=self.color)
        plt.plot(self.yaxis[:, 0], self.xaxis[:, 1], color=self.color)
        plt.plot(self.centers[:, 0], self.centers[:, 1], "o", color=self.color)


@dataclass
class Material:
    v: float  # Poisson's Ratio
    e: float  # Young's Modulus
    p: float  # Contact Curvature Radius


# Hertzian contact theory
def a(Fc, m_1: Material, m_2: Material, b, concavity):  # half the contact deformation width
    p_star = 0
    if concavity == "Convex and concave":
        p_star = abs((m_1.p * m_2.p) / (m_1.p - m_2.p))
    if concavity == "Convex and convex":
        p_star = abs((m_1.p * m_2.p) / (m_1.p + m_2.p))
    return np.sqrt((4 * Fc * p_star * (m_1.e * (1 - m_1.v**2) + m_2.e(1 - m_2.v**2))) / (np.pi * b * m_1.e * m_2.e))


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
