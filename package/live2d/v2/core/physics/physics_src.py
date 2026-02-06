from .iphysics_param import IPhysicsParam
from .physics_constants import SRC_TO_Y, SRC_TO_G_ANGLE, SRC_TO_X


class PhysicsSrc(IPhysicsParam):

    def __init__(self, paramId, aK, scale, weight):
        super().__init__(aK, scale, weight)
        self.tL_ = None
        self.tL_ = paramId

    def update(self, model, physicsContext):
        scaledValue = self.scale * model.getParamFloat(self.paramId)
        physicsPoint1 = physicsContext.getPhysicsPoint1()

        if self.tL_ == SRC_TO_X:
            physicsPoint1.x = physicsPoint1.x + (scaledValue - physicsPoint1.x) * self.weight
        elif self.tL_ == SRC_TO_Y:
            physicsPoint1.y = physicsPoint1.y + (scaledValue - physicsPoint1.y) * self.weight
        elif self.tL_ == SRC_TO_G_ANGLE:
            angle = physicsContext.getAngle()
            angle = angle + (scaledValue - angle) * self.weight
            physicsContext.setAngle(angle)
