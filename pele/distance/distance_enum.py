from enum import Enum, unique  # Package enum34


@unique
class Distance(Enum):
    CARTESIAN = 'cartesian'
    PERIODIC = 'periodic'
    LEES_EDWARDS = 'lees-edwards'
