#!/usr/bin/env python
# -*- coding: utf8 -*-

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