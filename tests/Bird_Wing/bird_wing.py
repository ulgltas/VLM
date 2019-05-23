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

        
def properties(yaw, aileron, flap, elevator, rudder):
    import wing
    import pythonVLM.VLM_inputs as VLMinputs

    wing_file = "models/Wing-4Sct.arp"
    wing_chordwise = 8
    wing_spanwise = 13
    w = VLMinputs.VLMDesc(wing_chordwise, wing_spanwise, wing_file)
    w.set_aileron(aileron)
    w.set_flap(flap)

    htail_file = "models/HTail.arp"
    htail_chordwise = 4
    htail_spanwise = 8
    htail = VLMinputs.VLMDesc(htail_chordwise, htail_spanwise, htail_file)
    htail.set_elevator(elevator)

    vtail_file = "models/VTail.arp"
    vtail_chordwise = 4
    vtail_spanwise = 8
    vtail = VLMinputs.VLMDesc(vtail_chordwise, vtail_spanwise, vtail_file)
    vtail.set_rudder(rudder)

    properties = VLMinputs.VLMProperties(w, htail, vtail)

    properties.u = 30
    properties.rho = 1.225
    properties.yaw = yaw
    properties.AoA = 5.0
    properties.timesteps = 25
    properties.denominator = 2.0

    properties.write_infile()

def run():
    import PyVLM

    PyVLM.run("infile.arp", "outfile.m")
