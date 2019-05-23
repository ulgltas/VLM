#!/usr/bin/env python
# -*- coding: utf8 -*-

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
