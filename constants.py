import numpy as np


# Параметры орбиты
APOCENTRE_HEIGHT = 350                         # км
PERICENTRE_HEIGHT = 240                        # км
INCLINATION = 10 * np.pi / 180                 # радиан
ASC_NODE_LONGITUDE = 5 * np.pi / 180           # радиан
PERICENTRE_ARGUMENT = 0 * np.pi / 180          # радиан
MEDIUM_ANOMALY = 240 * np.pi / 180              # радиан

# Параметры Земли
EARTH_MASS = 5.9742 * 1e24                     # кг
EARTH_RADIUS = 6378                            # км
GRAVITATIONAL_PARAMETER = 398600               # км^3 / c^2
GRAVITATIONAL_CONSTANT = 6.6743 * 1e-11        # м^3 / (кг * с^2)
ANGULAR_EARTH_VELOCITY = 7.29 * 1e-5           # рад / c
GEA_SEMI_MAJOR_AXIS = 6378.136                 # км
e_EARTH_SQUARED = 0.0067385254                 # без-ая величина

# Параметры атмосферы
NIGHT_ATMOSPHERE_DENSITY_120 = 1.58868 * 1e-8  # кг/м^3
DRAG_FORCE_COEFFICIENT = 3.5                   # без-ая величина

# Параметры космического аппарата
SPACECRAFT_MASS = 1650                         # кг
SPACECRAFT_AREA = 23 * 1e-6                    # км^2
