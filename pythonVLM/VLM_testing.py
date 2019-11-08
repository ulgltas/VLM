#!/usr/bin/env python
# -*- coding: utf8 -*-
def test(file, iterations, references, tolerances):
    global totalforce
    execfile(file, globals())

    force_x = totalforce[iterations-2]
    force_y = totalforce[iterations*2-2]
    force_z = totalforce[iterations*3-2]
    induced_drag = totalforce[-2]
    forces = [force_x, force_y, force_z, induced_drag]
    print(u"X-axis force: {} (expected {}+-{})".format(force_x, references[0], tolerances[0]))
    print(u"Y-axis force: {} (expected {}+-{})".format(force_y, references[1], tolerances[1]))
    print(u"Z-axis force: {} (expected {}+-{})".format(force_z, references[2], tolerances[2]))
    print(u"Induced drag: {} (expected {}+-{})".format(induced_drag, references[3], tolerances[3]))
    
    if any([abs(i-j)>k for i, j, k in zip(forces, references, tolerances)]):
        raise Exception("Incorrect results")
    else:
        print("Test passed")