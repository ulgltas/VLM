#!/usr/bin/env python
# -*- coding: utf8 -*-

import bird_wing

yaw = 0.0
aileron = 0.0
flap = 0.0
elevator = 0.0
rudder = 0.0

bird_wing.properties(yaw, aileron, flap, elevator, rudder)
bird_wing.run()

references = [-98.26, 0.0, 2826.6, 56.38]
tolerances = [0.01, 0.0001, 0.1, 0.01]

bird_wing.test("outfile.py", 25, references, tolerances)