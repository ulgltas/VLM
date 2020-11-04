#!/usr/bin/env python3
# -*- coding: utf8 -*-
# Copyright 2020 Université de Liège
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Test case for a sample light general aviation aircraft

import pythonVLM.VLM_inputs as inputs
import pythonVLM.VLM_driver as driver
import pythonVLM.VLM_testing as tests

# --- Wing ---
root_chord_w = 2.0
wing_file = "Wing.arp"
offset = [0.0, 0.0]

span = [5.0] # half-wing spans
taper = [0.9] # taper ratios
twist = [0.0, -1.0] # twist angle of each airfoil
sweep_LE = [5.0] # sweeps at leading edge
dihedral = [5.0] # dihedral angles
airfoils = ["NACA 2412.dat", "NACA 2412.dat"] # airfoils

wing = inputs.VLMWing(airfoils, span, taper, sweep_LE, dihedral, twist, root_chord_w, offset)

wing.chordwise_panels = 6 # number of panels along x
wing.spanwise_panels = 24 # number of panels along y
wing.geometry_file = wing_file

# Aileron
aileron_span = 2.0 # length
aileron_root_location = 3.0 # y-location
aileron_root_chord = 0.5 # x_location

# Add aileron and deflect it
wing.aileron.add_control(wing, aileron_span, aileron_root_location, aileron_root_chord)
wing.aileron.angle = 5.0 # degrees

# Flap (similar to aileron)
flap_span = 2.0
flap_root_location = 0.0
flap_root_chord = 0.5

# Add flap and deflect it
wing.flap.add_control(wing, flap_span, flap_root_location, flap_root_chord)
wing.flap.angle = 10.0

# Winglet (only one section)
winglet_span = 0.5 # length
winglet_chord = [1.8, 1.2] # chord at root and tip
winglet_sweep_LE = 30.0 # sweep at leading edge
winglet_dihedral = 60.0 # dihedral angle

# Add winglet
wing.winglet.add_winglet(winglet_span, winglet_chord, winglet_sweep_LE, winglet_dihedral, "NACA 0012")

# --- Horizontal tail (similar to wing) ---
root_chord_h = 1.0
htail_file = "HTail.arp"
offset_h = [5.0, 0.5]

span_h = [2.0]
taper_h = [0.9]
twist_h = [0.0, 0.0]
sweep_LE_h = [5.0]
dihedral_h = [0.0]
airfoils_h = ["NACA 0012.dat", "NACA 0012.dat"]

htail = inputs.VLMHTail(airfoils_h, span_h, taper_h, sweep_LE_h, dihedral_h, twist_h, root_chord_h, offset_h)

htail.chordwise_panels = 4
htail.spanwise_panels = 12
htail.geometry_file = htail_file

# Elevator (similar to aileron)
elev_span = 3.8 # length for FULL SPAN elevator
elev_root_location = 0.0
elev_root_chord = 0.3

# Add elevator and deflect it
htail.elevator.add_control(htail, elev_span, elev_root_location, elev_root_chord)
htail.elevator.angle = 10.0

# --- Vertical tail (similar to wing) ---
root_chord_v = 1.0
vtail_file = "VTail.arp"
offset_v = [5.0, 0.5]

span_v = [1.5]
taper_v = [0.7]
twist_v = [0.0, 0.0]
sweep_LE_v = [20.0]
dihedral_v = [0.0]
airfoils_v = ["NACA 0012.dat", "NACA 0012.dat"]

vtail = inputs.VLMVTail(airfoils_v, span_v, taper_v, sweep_LE_v, dihedral_v, twist_v, root_chord_v, offset_v)

vtail.chordwise_panels = 4
vtail.spanwise_panels = 12
vtail.geometry_file = vtail_file

# Rudder (similar to elevator)
rud_span = 1.4
rud_root_location = 0.0
rud_root_chord = 0.3

# Add rudder and deflect it
vtail.rudder.add_control(vtail, rud_span, rud_root_location, rud_root_chord)
vtail.rudder.angle = 10.0

# --- Flow properties ---
properties = inputs.VLMProperties(wing, htail, vtail)
properties.u = 50
properties.rho = 1.225
properties.AoA = 2.0
properties.yaw = 0.0
properties.timesteps = 20
properties.denominator = 2.0
properties.freewake = 0

# --- Write geometry and input files ---
wing.write_geofile(wing_file)
htail.write_geofile(htail_file)
vtail.write_geofile(vtail_file)
properties.write_infile()

# --- Execute ---
VLM = driver.VLMDriver(properties.infile)
VLM.run()
VLM.save()

# --- Test loads ---

references = [616.497859, -803.381703, 13038.832337, 655.011561]
tolerances = [0.01, 0.0001, 0.1, 0.01]

tests.test("outfile.py", properties.timesteps, references, tolerances)
