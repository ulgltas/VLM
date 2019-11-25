#!/usr/bin/env python
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
