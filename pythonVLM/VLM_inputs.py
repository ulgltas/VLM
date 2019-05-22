#!/usr/bin/env python
# -*- coding: utf8 -*-

import wing
import PyVLM


class VLMProperties:
    def __init__(self, wing, htail, vtail, infile="infile.arp"):
        self.u = 0.0
        self.rho = 0.0
        self.AoA = 0.0
        self.yaw = 0.0
        self.infile = infile
        self.timesteps = 0
        self.denominator = 1.0
        self.freewake = 1
        self.wing = wing
        self.htail = htail
        self.vtail = vtail

    def write_infile(self):
        f = open(self.infile, "w")
        f.write("INP1\t{}\t{}\t{}\t{}\n".format(
            self.u, self.rho, self.AoA, self.yaw))
        f.write("INP2\t{}\t{}\t{}\n".format(self.wing.chordwise_panels,
                                            self.htail.chordwise_panels, self.vtail.chordwise_panels))
        f.write("INP3\t{}\t{}\t{}\n".format(self.wing.spanwise_panels,
                                            self.htail.spanwise_panels, self.vtail.spanwise_panels))
        f.write("INP4\t{}\t{}\t{}\n".format(
            self.timesteps, self.denominator, self.freewake))
        f.write("INP5\t{}\t{}\t{}\t{}\n".format(self.wing.aileron.angle,
                                                self.wing.flap.angle, self.htail.elevator.angle, self.vtail.rudder.angle))
        f.write("INP6\t{}\n".format(self.wing.geometry_file))
        f.write("INP7\t{}\n".format(self.htail.geometry_file))
        f.write("INP8\t{}\n".format(self.vtail.geometry_file))
        f.close()

class VLMDesc:
    def __init__(self, chordwise_panels, spanwise_panels, geometry_file):
        self.chordwise_panels = chordwise_panels
        self.spanwise_panels = spanwise_panels
        self.geometry_file = geometry_file
    def set_aileron(self, angle=0.0):
        self.aileron = VLMControl(1, angle)
    def set_flap(self, angle=0.0):
        self.flap = VLMControl(1, angle)
    def set_elevator(self, angle=0.0):
        self.elevator = VLMControl(1, angle)
    def set_rudder(self, angle=0.0):
        self.rudder = VLMControl(1, angle)


class VLMSurface(wing.Wing):
    def initData(self, filenames, span, twist, sweep, dihedral, offset):
        import numpy as np
        import os
        import csv
        self.airfoils = []
        self.sweep_le = sweep
        self.dihedral = dihedral
        self.offset = offset
        self.twist = twist
        self.x_offset = [self.offset[0]]
        self.z_offset = [self.offset[1]]
        j = -1
        for filename in filenames:
            m = isNACA4or5(filename)
            j+=1
            if m:
                n = 101
                xpline = np.linspace(0, 1, n)
                airfoil = "NACA " + m.group("number")
                self.airfoils.append(airfoil)
                ycamber = PyVLM.nacafourfivedigit(xpline, n, airfoil)
                fname = os.path.join("models", airfoil+".dat")
                if not os.path.isfile(fname):
                    f = open(fname, "w")
                    f.write(airfoil + "\n")
                    for i in range(n-1, -1, -1):
                        f.write("\t{}\t{}\n".format(xpline[i], ycamber[i]))
                    for i in range(0, n):
                        f.write("\t{}\t{}\n".format(xpline[i], ycamber[i]))
                    f.close()
            else:
                fnamelist = filename.split(".")
                airfoil = fnamelist[0]
                fname = os.path.join("models", airfoil+".dat")
                fname_arf = os.path.join("models", airfoil+".arf")
                self.airfoils.append(os.path.join("models", airfoil))
                if not os.path.isfile(fname) and os.path.isfile(fname_arf):
                    arf_reader = csv.DictReader(fname_arf, delimiter=" ")
                    f = open(fname, "w")
                    f.write(airfoil + "\n")
                    for row in arf_reader:
                        if len(row)>2:
                            f.write("\t{}\t{}\n".format(float(row[1])/100, float(row[2])/100))
                    for row in arf_reader:
                        if len(row)>2:
                            f.write("\t{}\t{}\n".format(float(row[3])/100, float(row[4])/100))
                    f.close()
                elif os.path.isfile(fname) and not os.path.isfile(fname_arf):
                    data = self.read(fname)
                    n = len(data)
                    l_data = (n+1)/2
                    f = open(fname_arf, "w")
                    f.write("GMD401 {}\n".format(l_data))
                    f.write("GMD402 {}\n".format(l_data))
                    for i in range(l_data):
                        f.write("GM15 {} {} {} {}\n".format(data[l_data-i-1,0]*100, data[l_data-i-1,1]*100, data[l_data+i-1,0]*100, data[l_data+i-1,1]*100))
                    f.write("GMD6 \n")
                    f.close()
                elif not os.path.isfile(fname) and not os.path.isfile(fname_arf):
                    raise Exception("Neither arf nor dat file found for airfoil {}".format(airfoil))
        filenames = map("models/{}".format, filenames)
        wing.Wing.initData(self, filenames, span, twist, sweep, dihedral, offset)
        self.y_offset = self.spanPos
        for j in range(1,self.n):
            a = self.specPts(j)
            self.x_offset.append(self.pts[j][a[3],0])
            self.z_offset.append(self.pts[j][a[3],2])
        self.calculate_mac()
    def calculate_mac(self):
        mac_integrand = 0.0 # Integrand of the mean aerodynamic chord

        for i in range(1,self.n):
            dy = self.spanPos[i]-self.spanPos[i-1]
            dc = self.chord[i-1]-self.chord[i]
            mac_integrand += dy*(self.chord[i-1]**2.0-dc*self.chord[i-1]+dc**2.0/3)
        self.mac = 1/self.S*mac_integrand
        # Find MAC location
        self.mac_x = 0.0
        self.mac_y = 0.0
        self.mac_z = 0.0
        for i in range(1,self.n):
            taper = self.chord[i]/self.chord[i-1]
            if taper<1.0 and self.chord[i]<self.mac and self.chord[i-1]>self.mac:
                rel_pos = (self.mac-self.chord[i])/(self.chord[i-1]-self.chord[i])
                self.mac_x = self.x_offset[i]*rel_pos+self.x_offset[i-1]*(1-rel_pos)
                self.mac_y = self.y_offset[i]*rel_pos+self.y_offset[i-1]*(1-rel_pos)
                self.mac_z = self.z_offset[i]*rel_pos+self.z_offset[i-1]*(1-rel_pos)
            elif taper==1.0 and self.chord[i]==self.mac:
                self.mac_x = self.x_offset[i-1]
                self.mac_y = self.y_offset[i-1]
                self.mac_z = self.z_offset[i-1]
            elif taper>1.0:
                rel_pos = (self.mac-self.chord[i-1])/(self.chord[i]-self.chord[i-1])
                self.mac_x = self.x_offset[i]*(1-rel_pos)+self.x_offset[i-1]*rel_pos
                self.mac_y = self.y_offset[i]*(1-rel_pos)+self.y_offset[i-1]*rel_pos
                self.mac_z = self.z_offset[i]*(1-rel_pos)+self.z_offset[i-1]*rel_pos
    def get_chord(self, span):
        i = 0
        while i<self.n-1 and self.spanPos[i+1]>span:
            i+=1
        rel_pos = (span-self.spanPos[i-1])/(self.spanPos[i]-self.spanPos[i-1])
        return self.chord[i]*rel_pos+self.chord[i-1]*(1-rel_pos)



class VLMWing(VLMSurface):
    kind = "WNG"

    def initData(self, filenames, span, twist, sweep, dihedral, offset):
        VLMSurface.initData(self, filenames, span, twist, sweep, dihedral, offset)
        self.aileron = VLMControl(0, 0.0)
        self.flap = VLMControl(0, 0.0)
        self.winglet = VLMWinglet(0)


class VLMHTail(VLMSurface):
    kind = "HTL"

    def initData(self, filenames, span, twist, sweep, dihedral, offset):
        VLMSurface.initData(self, filenames, span, twist, sweep, dihedral, offset)
        self.elevator = VLMControl(0, 0.0)


class VLMVTail(VLMSurface):
    kind = "VTL"

    def initData(self, filenames, span, twist, sweep, dihedral, offset):
        VLMSurface.initData(self, filenames, span, twist, sweep, dihedral, offset)
        self.rudder = VLMControl(0, 0.0)



class VLMControl:
    def __init__(self, exists, angle):
        self.angle = angle
        self.exists = exists
    def add_control(self, wing, span, root, chord):
        self.exists = 1
        self.span = span
        self.root = root
        self.rel_chord = chord/wing.get_chord(root)*100.0
        self.rel_span = self.span/wing.b*100.0

class VLMWinglet:
    kind = "WGL"
    def __init__(self, exists):
        self.exists = exists
    def add_winglet(self, span, chords, sweep, dihedral, airfoil):
        self.exists = 1
        self.span = span
        self.root_chord = chords[0]
        self.tip_chord = chords[-1]
        self.area = self.span*(self.root_chord+self.tip_chord)/2
        self.sweep_le = sweep
        self.dihedral = dihedral
        self.airfoil = airfoil


def write_geofile(filename, wing):
    import numpy as np
    offset_p = np.argmin(wing.pts[-1][:, 0])
    end_x = wing.pts[-1][offset_p, 0]-wing.offset[0]
    end_z = wing.pts[-1][offset_p, 2]-wing.offset[1]

    sweep = np.degrees(np.arctan2(end_x, wing.b))
    dihedral = np.degrees(np.arctan2(end_z, wing.b))
    taper = wing.chord[-1]/wing.chord[0]

    f = open(filename, "w")
    f.write("WNG101:\t{}\t-01\n".format(wing.n-1))
    f.write("WNG401:\t000\t{}\n".format(wing.n-1))

    for i in range(0, wing.n-1):
        dy = wing.spanPos[i+1]-wing.spanPos[i]
        dx = dy*np.tan(wing.sweep_le[i])-wing.chord[i]/4+wing.chord[i+1]/4
        sweep_c4 = np.degrees(np.arctan2(dx, dy))
        dx = dx-wing.chord[i]/4+wing.chord[i+1]/4
        sweep_c2 = np.degrees(np.arctan2(dx, dy))
        f.write("WN{}501:\t{}\t{}\n".format(i+1, wing.airfoils[i], wing.airfoils[i+1]))
        f.write("WN{}502:\t{}\t{}\t{}\n".format(
            i+1, dy, wing.chord[i], wing.chord[i+1]))
        f.write("WN{}503:\t{}\t{}\t{}\t{}\t{}\n".format(
            i+1, 0.0, 0.0, np.rad2deg(wing.sweep_le[i]), np.rad2deg(wing.dihedral[i]), np.rad2deg(wing.twist[i])))
        f.write("WN{}504: \n".format(i+1))
        f.write("WN{}505: \n".format(i+1))
        f.write("WN{}601: \n".format(i+1))
        f.write("WN{}602: \n".format(i+1))
        f.write("WN{}603:\t{}\t{}\n".format(i+1, sweep_c4, sweep_c2))
    f.write("WNG601:\t{}\t{}\t{}\t{}\t{}\n".format(
        2*wing.b, wing.offset[0], wing.offset[1], wing.chord[0], wing.chord[-1]))
    f.write("WNG602:\t{}\n".format(2*wing.S))
    f.write("WNG603: \n")
    f.write("WNG604:\t{}\t{}\t{}\t{}\n".format(sweep, 0.0, dihedral, 0.0))
    f.write("WNG605:\t{}\t{}\t{}\n".format(wing.AR, taper, wing.AR))
    f.write("WNG606: \n")
    f.write("WNG607:\t{}\t{}\t{}\t{}\n".format(wing.mac, wing.mac_x, wing.mac_y, wing.mac_z))

    f.write("AIL201:\t{}\n".format(wing.aileron.exists))
    if wing.aileron.exists==1:
        f.write("AIL601:\t{}\t{}\n".format(wing.aileron.span, wing.aileron.root))
        f.write("AIL603:\t{}\t{}\t{}\n".format(100.0-wing.aileron.rel_chord, wing.aileron.rel_chord, wing.aileron.rel_span))
    f.write("LED201:\t0\n")
    f.write("TED201:\t{}\n".format(wing.flap.exists))
    if wing.flap.exists==1:
        f.write("TED601:\t{}\t{}\n".format(wing.flap.span, wing.flap.root))
        f.write("TED603:\t{}\t{}\t{}\t{}\n".format(100.0-wing.flap.rel_chord, wing.flap.rel_chord, wing.flap.rel_span, 1.0))
    f.write("ABK201:\t0\n")
    f.write("WGL201:\t{}\n".format(wing.winglet.exists))
    if wing.winglet.exists==1:
        f.write("WGL102:\t{}\n".format(wing.winglet.airfoil))
        f.write("WGL601:\t{}\t{}\t{}\n".format(wing.winglet.span, wing.winglet.root_chord, wing.winglet.tip_chord))
        f.write("WGL602:\t{}\t{}\t{}\n".format(wing.winglet.area, wing.winglet.sweep_le, wing.winglet.dihedral))
    else:
            f.write("WGL601:\t0\n")
    f.write("FLA201:\t0\n")
    f.close()


def isNACA4or5(airfoil):
    import re
    return re.match("NACA ?(?P<number>\d{4,5})\D*", airfoil)
