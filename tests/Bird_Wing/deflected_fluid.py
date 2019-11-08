#!/usr/bin/env python
# -*- coding: utf8 -*-

import bird_wing
import pythonVLM.VLM_testing as VLMtest

yaw = 0.0
aileron = 0.0
flap = 45.0
elevator = -15.0
rudder = 0.0

bird_wing.properties(yaw, aileron, flap, elevator, rudder)
bird_wing.run()

references = [847.91, 0.0, 7409.0, 970.44]
tolerances = [0.01, 0.0001, 0.1, 0.01]

VLMtest.test("outfile.py", 25, references, tolerances)