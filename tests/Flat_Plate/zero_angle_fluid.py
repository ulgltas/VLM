import wing
import pythonVLM.VLM_inputs as inputs
import pythonVLM.VLM_testing as tests

airfoils = ["NACA 0001.dat", "NACA 0001.dat"]
span = [4.0]
taper = [1.0, 1.0]
twist = [0.0, 0.0]
sweep = [0.0, 0.0]
dihedral = [0.0, 0.0]
root_chord = 1.0

offset = [0.0, 0.0]

spanwise_panels = 5

w = inputs.VLMWing(airfoils, span, taper, sweep, dihedral, twist, root_chord, offset)

v_airfoils = ["NACA 0009.dat", "NACA 0009.dat", "NACA 0009.dat"]

v_span = [0.5, 0.5]
v_taper = [0.5, 0.5, 0.5]
v_sweep = [45.0, 20.0, 20.0]
v_dihedral = [0.0, 0.0, 0.0]
v_twist = [0.0, 0.0, 0.0]
v_root_chord = 0.5
v_offset = [5.0, 6.0]

vtail = inputs.VLMVTail(v_airfoils, v_span, v_taper, v_sweep, v_dihedral, v_twist, v_root_chord, v_offset)

h_span = [1.0]
h_taper = [0.5, 0.5]
h_sweep = [45.0, 45.0]
h_dihedral = [0.0, 0.0]
h_twist = [0.0, 0.0]
h_root_chord = 0.5
h_offset = [5.0, 6.0]


htail = inputs.VLMHTail(airfoils, h_span, h_taper, h_sweep, h_dihedral, h_twist, h_root_chord, h_offset)

w.write_geofile("Wing.arp")
htail.write_geofile("HTail.arp")
vtail.write_geofile("VTail.arp")

properties = inputs.VLMProperties(w, htail, vtail)
properties.u = 30
properties.rho = 1.225
properties.AoA = 0.0
properties.timesteps = 20
properties.wing.chordwise_panels = 5
properties.wing.spanwise_panels = 8
properties.wing.geometry_file = "Wing.arp"

properties.htail.chordwise_panels = 5
properties.htail.spanwise_panels = 8
properties.htail.geometry_file = "HTail.arp"

properties.vtail.chordwise_panels = 5
properties.vtail.spanwise_panels = 8
properties.vtail.geometry_file = "VTail.arp"


properties.write_infile()

import CVLM

CVLM.run("infile.arp","outfile.m")

tests.test("outfile.py", properties.timesteps, [0.0, 0.0, 0.0, 0.0], [0.00001, 0.00001, 0.00001, 0.00001])