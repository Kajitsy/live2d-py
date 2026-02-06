import math

from ..type import Array
from ..util import UtMath
from .physics_constants import *
from .physics_point import PhysicsPoint
from .physics_src import PhysicsSrc
from .physics_target import PhysicsTarget


class PhysicsHair:
    SRC_TO_X = SRC_TO_X
    SRC_TO_Y = SRC_TO_Y
    SRC_TO_G_ANGLE = SRC_TO_G_ANGLE

    TARGET_FROM_ANGLE = TARGET_FROM_ANGLE
    TARGET_FROM_ANGLE_V = TARGET_FROM_ANGLE_V

    def __init__(self):
        self.p1 = PhysicsPoint()
        self.p2 = PhysicsPoint()
        self.Fo_ = 0
        self.Db_ = 0
        self.L2_ = 0
        self.M2_ = 0
        self.ks_ = 0
        self._9b = 0
        self.iP_ = 0
        self.iT_ = 0
        self.lL_ = Array()
        self.qP_ = Array()
        self.setup(0.3, 0.5, 0.1)

    def setup(self, length=None, stiffness=None, mass=None):
        self.ks_ = self.calcAngle()
        self.p2.setupLast()
        if mass is not None:
            self.Fo_ = length
            self.L2_ = stiffness
            self.p1.mass = mass
            self.p2.mass = mass
            self.p2.y = length
            self.setup()

    def getPhysicsPoint1(self):
        return self.p1

    def getPhysicsPoint2(self):
        return self.p2

    def getAngle(self):
        return self.Db_

    def setAngle(self, angle):
        self.Db_ = angle

    def getLastAngle(self):
        return self.M2_

    def getLastAngleVelocity(self):
        return self._9b

    def calcAngle(self):
        return -180 * (math.atan2(self.p1.x - self.p2.x, -(self.p1.y - self.p2.y))) / math.pi

    def addSrcParam(self, paramId, type, scale, weight):
        srcParam = PhysicsSrc(paramId, type, scale, weight)
        self.lL_.append(srcParam)

    def addTargetParam(self, paramId, type, scale, weight):
        targetParam = PhysicsTarget(paramId, type, scale, weight)
        self.qP_.append(targetParam)

    def update(self, model, timeMSec):
        if self.iP_ == 0:
            self.iP_ = self.iT_ = timeMSec
            self.Fo_ = (math.sqrt((self.p1.x - self.p2.x) * (self.p1.x - self.p2.x) + (self.p1.y - self.p2.y) * (
                    self.p1.y - self.p2.y)))
            return

        deltaTime = (timeMSec - self.iT_) / 1000
        if deltaTime != 0:
            for i in range(len(self.lL_) - 1, 0 - 1, -1):
                srcParam = self.lL_[i]
                srcParam.update(model, self)

            self.updatePhysics(model, deltaTime)
            self.M2_ = self.calcAngle()
            self._9b = (self.M2_ - self.ks_) / deltaTime
            self.ks_ = self.M2_

        for i in range(len(self.qP_) - 1, 0 - 1, -1):
            targetParam = self.qP_[i]
            targetParam.update(model, self)

        self.iT_ = timeMSec

    def updatePhysics(self, model, deltaTime):
        if deltaTime < 0.033:
            deltaTime = 0.033

        invTime = 1 / deltaTime
        self.p1.vx = (self.p1.x - self.p1.lastX) * invTime
        self.p1.vy = (self.p1.y - self.p1.lastY) * invTime
        self.p1.ax = (self.p1.vx - self.p1.lastVX) * invTime
        self.p1.ay = (self.p1.vy - self.p1.lastVY) * invTime
        self.p1.fx = self.p1.ax * self.p1.mass
        self.p1.fy = self.p1.ay * self.p1.mass
        self.p1.setupLast()
        angle = -(math.atan2((self.p1.y - self.p2.y), self.p1.x - self.p2.x))
        cosAngle = math.cos(angle)
        sinAngle = math.sin(angle)
        gravity = 9.8 * self.p2.mass
        angleRad = (self.Db_ * UtMath.DEG_TO_RAD)
        forceGravity = (gravity * math.cos(angle - angleRad))
        forceX = (forceGravity * sinAngle)
        forceY = (forceGravity * cosAngle)
        forceDragX = (-self.p1.fx * sinAngle * sinAngle)
        forceDragY = (-self.p1.fy * sinAngle * cosAngle)
        forceStiffnessX = (-self.p2.vx * self.L2_)
        forceStiffnessY = (-self.p2.vy * self.L2_)
        self.p2.fx = (forceX + forceDragX + forceStiffnessX)
        self.p2.fy = (forceY + forceDragY + forceStiffnessY)
        self.p2.ax = self.p2.fx / self.p2.mass
        self.p2.ay = self.p2.fy / self.p2.mass
        self.p2.vx += self.p2.ax * deltaTime
        self.p2.vy += self.p2.ay * deltaTime
        self.p2.x += self.p2.vx * deltaTime
        self.p2.y += self.p2.vy * deltaTime
        distance = (math.sqrt(
            (self.p1.x - self.p2.x) * (self.p1.x - self.p2.x) + (self.p1.y - self.p2.y) * (self.p1.y - self.p2.y)))
        self.p2.x = self.p1.x + self.Fo_ * (self.p2.x - self.p1.x) / distance
        self.p2.y = self.p1.y + self.Fo_ * (self.p2.y - self.p1.y) / distance
        self.p2.vx = (self.p2.x - self.p2.lastX) * invTime
        self.p2.vy = (self.p2.y - self.p2.lastY) * invTime
        self.p2.setupLast()
