from enum import Enum
class Instrument(object):
    def __init__(self,instrument_name,axes):
        self.i_name = instrument_name
        self.axes = axes
    def instrument_info(self):
        return f"Instrument Name: {self.i_name}, Axes: {self.axes}"

axes = {
    "Axis0" : "Rot,0,1,0,1"
}
MAPS = Instrument("MAPS",axes)
MERLIN = Instrument("MERLIN", axes)
LET = Instrument("LET", axes)

class reflection_condition(Enum):
    PRIMITIVE = "Primitive"
    BODY_CENTERED= "Body Centered"
    C_FACE_CENTERED = "C-face centered"
    B_FACE_CENTERED = "B-face centered"
    A_FACE_CENTERED = "A-face centered"
    ALL_FACE_CENTERED = "All-face centered"
    RHOMBO_CENTERED_REV = "Rhombohedrally centred reverse"
    HEX_CENTERED_REV = "Hexagonally centred, reverse"



class UB_METHOD(Enum):
    ONE_PEAK = 1
    TWO_PEAKS = 2
    THREE_PEAKS = 3
    FIVEPLUS_PEAKS = 4

class Peak_picking_methods(Enum):
    MANUAL = 1
    ROI = 2
