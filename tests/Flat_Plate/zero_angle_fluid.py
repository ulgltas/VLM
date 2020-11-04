#!/usr/bin/env python3
# -*- coding: utf8 -*-
# Copyright 2019 Université de Liège
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

import wing
import pythonVLM.VLM_inputs as inputs
import pythonVLM.VLM_testing as tests

airfoils = ["NACA 0001.dat", "NACA 0001.dat"] # Symmetric airfoil
span = [4.0]
taper = [1.0]
twist = [0.0, 0.0]
sweep = [0.0]
dihedral = [0.0]
root_chord = 1.0

offset = [0.0, 0.0]

w = inputs.VLMWing(airfoils, span, taper, sweep, dihedral, twist, root_chord, offset)

# Empty VTail
v_airfoils = []

v_span = []
v_taper = []
v_sweep = []
v_dihedral = []
v_twist = []
v_root_chord = 0.5
v_offset = [0.0, 0.0]

vtail = inputs.VLMVTail(v_airfoils, v_span, v_taper, v_sweep, v_dihedral, v_twist, v_root_chord, v_offset)

# Empty HTail
h_airfoils = []
h_span = []
h_taper = []
h_sweep = []
h_dihedral = []
h_twist = []
h_root_chord = 0.5
h_offset = [0., 0.]


htail = inputs.VLMHTail(h_airfoils, h_span, h_taper, h_sweep, h_dihedral, h_twist, h_root_chord, h_offset)

w.write_geofile("Wing.arp")

properties = inputs.VLMProperties(w, htail, vtail)
properties.u = 30
properties.rho = 1.225
properties.AoA = 0.0
properties.timesteps = 20
properties.wing.chordwise_panels = 5
properties.wing.spanwise_panels = 8
properties.wing.geometry_file = "Wing.arp"


properties.write_infile()

import CVLM

CVLM.run("infile.arp","outfile.m")

tests.test("outfile.py", properties.timesteps, [0.0, 0.0, 0.0, 0.0], [0.00001, 0.00001, 0.00001, 0.00001])
