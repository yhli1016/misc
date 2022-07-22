#! /usr/bin/env python

import math
import numpy as np
import matplotlib.pyplot as plt


# Shared constants
H_EYE = 1.65
H_ARROW = 1.50
L_ARROW = 1.0
DISTANCE = 10


def solve_hit(height_aim, hg):
    # Shortcuts for lengths
    AB = H_EYE - H_ARROW
    BG = L_ARROW

    # Solve right-angled triangle BGH
    BH = math.sqrt(BG**2 + hg**2)
    HBG = math.atan(hg / BG)

    # Solve triangle ABH
    BAH = math.pi / 2 + math.atan((height_aim - H_EYE) / DISTANCE)
    AHB = math.asin(math.sin(BAH) * AB / BH)
    ABH = math.pi - BAH - AHB

    # Evaluate height of hit
    EBI = math.pi / 2 - ABH - HBG
    height_hit = H_ARROW + DISTANCE * math.tan(EBI)
    return height_hit


def main():
    height_aim = np.linspace(1.5, 2.0, 6)
    hg = np.linspace(0.10, 0.20, 11)
    for h1 in height_aim:
        result = []
        for h2 in hg:
            result.append(solve_hit(height_aim=h1, hg=h2))
        plt.plot(hg, result)
    plt.show()

    height_aim = np.linspace(1.5, 2.0, 6)
    hg = np.linspace(0.10, 0.20, 11)
    for h2 in hg:
        result = []
        for h1 in height_aim:
            result.append(solve_hit(height_aim=h1, hg=h2))
        plt.plot(height_aim, result)
    plt.show()


if __name__ == "__main__":
    main()
