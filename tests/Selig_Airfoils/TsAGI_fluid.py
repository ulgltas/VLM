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

import pythonVLM.VLM_inputs as inputs
import pythonVLM.VLM_driver as driver
import pythonVLM.VLM_testing as tests

## Wing
root_chord_w = 1.0
wing_file = "Wing.arp"
offset = [0.0, 0.0]

span = [2.0, 2.0]
taper = [0.5, 0.5]
twist = [0.0, 2.5, 5.0]
sweep_LE = [0.0, -25.0]
dihedral = [0.0, 5.0]
airfoils = ["TsAGI_R3a.dat", "TsAGI_R3a.dat", "TsAGI_R3a.dat"]

wing = inputs.VLMWing(airfoils, span, taper, sweep_LE, dihedral, twist, root_chord_w, offset)

wing.chordwise_panels = 6
wing.spanwise_panels = 24
wing.geometry_file = wing_file


## Empty horizontal tail
root_chord_h = 0.5
htail_file = "HTail.arp"
offset_h = [5.0, 1.0]

span_h = []
taper_h = []
twist_h = []
sweep_LE_h = []
dihedral_h = []
airfoils_h = []

htail = inputs.VLMHTail(airfoils_h, span_h, taper_h, sweep_LE_h, dihedral_h, twist_h, root_chord_h, offset_h)

htail.chordwise_panels = 0
htail.spanwise_panels = 0
htail.geometry_file = htail_file

## Empty vertical tail
root_chord_v = 0.5
vtail_file = "VTail.arp"
offset_v = [5.0, 1.0]

span_v = []
taper_v = []
twist_v = []
sweep_LE_v = []
dihedral_v = []
airfoils_v = []

vtail = inputs.VLMVTail(airfoils_v, span_v, taper_v, sweep_LE_v, dihedral_v, twist_v, root_chord_v, offset_v)

vtail.chordwise_panels = 0
vtail.spanwise_panels = 0
vtail.geometry_file = vtail_file

## Write geometry files
wing.write_geofile(wing_file)

## Flow properties
properties = inputs.VLMProperties(wing, htail, vtail)
properties.u = 50.0
properties.rho = 1.225
properties.AoA = 5.0
properties.yaw = 0.0
properties.timesteps = 20
properties.denominator = 2.0
properties.freewake = 0

## Input file
properties.write_infile()

VLM = driver.VLMDriver(properties.infile)
VLM.run()
VLM.save()

tests.test("outfile.py", properties.timesteps, [29.047813, 0.000000, 5286.971696, 213.431801], [0.00001, 0.00001, 0.00001, 0.00001])
