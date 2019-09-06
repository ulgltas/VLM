import CVLM
import numpy as np

class VLMDriver(object):
    def __init__(self, infile):
        self.iteration = 0
        self.data = CVLM.VLMData()
        CVLM.setup(infile, self.data)
        self.m = self.data.wing.nshed+2 # Assumes only wing
        self.n = self.data.wing.nvert/self.m
        CVLM.geometry_setup(self.data)
        CVLM.cycleliftsurf(self.data)
        CVLM.memory_setup(self.data)
        self.x = np.zeros(self.data.wing.nvert)
        self.y = np.zeros(self.data.wing.nvert)
        self.z = np.zeros(self.data.wing.nvert)
        self.xv = np.zeros(self.data.wing.nvert)
        self.yv = np.zeros(self.data.wing.nvert)
        self.zv = np.zeros(self.data.wing.nvert)
        for ii in range(self.data.wing.nvert):
            self.x[ii] = self.getX(ii)
            self.y[ii] = self.getY(ii)
            self.z[ii] = self.getZ(ii)
            self.xv[ii] = self.getXv(ii)
            self.yv[ii] = self.getYv(ii)
            self.zv[ii] = self.getZv(ii)
        self.S = 0. # Surface of the wing, imperfect
        for i in range(self.data.wing.nface):
            self.S += CVLM.nsurf_getitem(self.data.wing.nsurf, i)
    def run(self):
        # Run steady simulation
        for i in range(self.data.ntimes-1):
            self.iteration = i
            CVLM.iteration(self.data, self.iteration)
    def __getCoord(self, index, delta):
        # Get the x/y/z (delta=0/1/2, respectively) of vertex index
        if index>self.data.wing.nvert+self.data.flap.nvert:
            c = CVLM.vertices_getitem(self.data.aileron.vertices,
                index-self.data.wing.nvert-self.data.flap.nvert+delta*self.data.aileron.nvert)
        elif index>self.data.wing.nvert:
            c = CVLM.vertices_getitem(self.data.flap.vertices,
                index-self.data.wing.nvert+delta*self.data.flap.nvert)
        else:
            c = CVLM.vertices_getitem(self.data.wing.vertices,
                index+delta*self.data.wing.nvert)
        return c
    def __setCoord(self, index, c, delta):
        # Set the x/y/z (delta=0/1/2, respectively) of vertex index
        if index>self.data.wing.nvert+self.data.flap.nvert:
            CVLM.vertices_setitem(self.data.aileron.vertices,
                index-self.data.wing.nvert-self.data.flap.nvert+delta*self.data.aileron.nvert, c)
        elif index>self.data.wing.nvert:
            CVLM.vertices_setitem(self.data.flap.vertices,
                index-self.data.wing.nvert+delta*self.data.flap.nvert, c)
        else:
            CVLM.vertices_setitem(self.data.wing.vertices,
                index+delta*self.data.wing.nvert, c)
    def __getVortexCoord(self, index, delta):
        return CVLM.vortex_getitem(self.data.wing.vortex, index+delta*self.data.wing.nvert)
    def __setVortexCoord(self, index, c, delta):
            CVLM.vortex_setitem(self.data.wing.vortex, index+delta*self.data.wing.nvert, c)
    def getX(self, index):
        return self.__getCoord(index, 0)
    def getY(self, index):
        return self.__getCoord(index, 1)
    def getZ(self, index):
        return self.__getCoord(index, 2)
    def getXv(self, index):
        return self.__getVortexCoord(index, 0)
    def getYv(self, index):
        return self.__getVortexCoord(index, 1)
    def getZv(self, index):
        return self.__getVortexCoord(index, 2)
    def setX(self, index, x):
        self.__setCoord(index, x, 0)
    def setY(self, index, y):
        self.__setCoord(index, y, 1)
    def setZ(self, index, z):
        self.__setCoord(index, z, 2)
    def dX(self, index, dx):
        x = self.getX(index)
        self.__setCoord(index, self.x[index]+dx, 0)
    def dY(self, index, dy):
        y = self.getY(index)
        self.__setCoord(index, self.y[index]+dy, 1)
    def dZ(self, index, dz):
        z = self.getZ(index)
        self.__setCoord(index, self.z[index]+dz, 2)
    def dXv(self, index, dx):
        x = self.__getVortexCoord(index, 0)
        self.__setVortexCoord(index, self.xv[index]+dx, 0)
    def dYv(self, index, dy):
        y = self.__getVortexCoord(index, 1)
        self.__setVortexCoord(index, self.yv[index]+dy, 1)
    def dZv(self, index, dz):
        z = self.__getVortexCoord(index, 2)
        self.__setVortexCoord(index, self.zv[index]+dz, 2)
    def getVertices(self, panel):
        v = [-1,-1,-1,-1]
        if panel < 10000:
            for j in range(4):
                v[j] = CVLM.faces_getitem(self.data.wing.faces, panel+j*self.data.wing.nface)
        return v
    def __getLoads(self, delta):
        return CVLM.aeroforce_getitem(self.data.wing.aeroforce, delta)
    def getQ(self):
        Q = 0.5*(CVLM.aeroforce_getitem(self.data.UVW, 0)**2+CVLM.aeroforce_getitem(self.data.UVW, 1)**2+
                CVLM.aeroforce_getitem(self.data.UVW, 2)**2)*self.data.rho
        return Q
    def getCl(self):
        x_force = self.__getLoads(0)
        z_force = self.__getLoads(2)
        lift = z_force*np.cos(self.data.aoa)-x_force*np.sin(self.data.aoa)
        return lift/(self.getQ()*self.S)
    def getCd(self):
        x_force = self.__getLoads(0)
        z_force = self.__getLoads(2)
        induced_drag = self.__getLoads(3)
        drag = z_force*np.sin(self.data.aoa)+x_force*np.cos(self.data.aoa)+induced_drag
        return drag/(self.getQ()*self.S)
    def getdeltaP(self, panel):
        if panel < 10000:
            dP = CVLM.Deltap_getitem(self.data.wing.Deltap, panel)
            normal = [CVLM.normal_getitem(self.data.wing.normal, panel), CVLM.normal_getitem(self.data.wing.normal, panel+self.data.wing.nface), CVLM.normal_getitem(self.data.wing.normal, panel+2*self.data.wing.nface)]
            dPnormal = [comp * dP for comp in normal]
        return dPnormal
    def getSurface(self, panel):
        if panel < 10000:
            return CVLM.nsurf_getitem(self.data.wing.nsurf, panel)
    def getForce(self, panel, weight):
        if panel < 10000:
            dP = CVLM.Deltap_getitem(self.data.wing.Deltap, panel)
            normal = [CVLM.normal_getitem(self.data.wing.normal, panel), CVLM.normal_getitem(self.data.wing.normal, panel+self.data.wing.nface), CVLM.normal_getitem(self.data.wing.normal, panel+2*self.data.wing.nface)]
            S = self.getSurface(panel)
            forceNormal = [comp * dP * S * weight for comp in normal]
        return forceNormal
    def update(self):
        CVLM.reset_wake(self.data)
        CVLM.geometry_setup(self.data)
        CVLM.colvec(self.data.wing)
        CVLM.cycleliftsurf(self.data)
    def save(self):
        CVLM.exportTextOutput("outfile.m",self.iteration-1,self.data)