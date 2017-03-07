from enum import Enum, unique  # Package enum34


@unique
class Distance(Enum):
    CARTESIAN = 1
    PERIODIC = 2
    LEES_EDWARDS = 3
