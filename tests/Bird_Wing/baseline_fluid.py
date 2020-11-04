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

import bird_wing
import pythonVLM.VLM_testing as VLMtest

yaw = 5.0
aileron = 0.0
flap = 0.0
elevator = 0.0
rudder = 0.0

bird_wing.properties(yaw, aileron, flap, elevator, rudder)
bird_wing.run()

references = [-97.57, -563.3, 2807.5, 58.85]
tolerances = [0.01, 0.1, 0.1, 0.01]

VLMtest.test("outfile.py", 25, references, tolerances)