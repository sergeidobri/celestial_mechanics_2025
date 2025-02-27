import numpy as np
from constants import e_EARTH_SQUARED


def find_eccentric_anomaly(medium_anomaly: float,  # радианы
                           eccentricity: float,
                           precision=1e-3) -> float:
    previous_eccentric_anomaly = 0
    new_eccentric_anomaly = (medium_anomaly + 
                            eccentricity * np.sin(previous_eccentric_anomaly))
    while(abs(new_eccentric_anomaly - previous_eccentric_anomaly) > precision):
        previous_eccentric_anomaly = new_eccentric_anomaly
        new_eccentric_anomaly = (medium_anomaly + 
                            eccentricity * np.sin(previous_eccentric_anomaly))

    return new_eccentric_anomaly  # радианы

def find_true_anomaly(eccentric_anomaly: float,
                      eccentricity: float,
                      medium_anomaly: float) -> float:
    result = 2 * np.atan(
        ((1 + eccentricity) / (1 - eccentricity)) ** 0.5 * np.tan(eccentric_anomaly / 2)
    )
    if medium_anomaly > np.pi:
        result += 2 * np.pi
    return result


def find_geodetic_coordinates(x: float,
                              y: float,
                              z: float, 
                              a: float, 
                              e_earth_sq: float) -> tuple:
    D = (x**2 + y**2) ** 0.5  # км

    if D == 0:
        B = (np.pi/2) * (z / abs(z))  # рад
        L = 0
        H = z * np.sin(B) - a * ((1 - e_EARTH_SQUARED * (np.sin(B) ** 2)) ** 0.5)
        return L, B, H
    elif D > 0:
        L_a = np.arcsin(y/D)
        L = None
        if y < 0 and x > 0:
            L = 2 * np.pi - L_a
        elif y < 0 and x < 0:
            L = np.pi + L_a
        elif y > 0 and x < 0:
            L = np.pi - L_a
        elif y > 0 and x > 0:
            L = L_a
        
        if z == 0:
            B = 0
            H = D - a
            return L, B, H
        
        r = (x**2 + y**2 + z**2) ** 0.5
        c = np.arcsin(z/r)
        p = (e_earth_sq * a) / (2*r)
        s_1 = 1e3
        s_2 = 0
        precision = 1e-4
        b = None
        while abs(s_1 - s_2) >= precision:
            s_1 = s_2
            b = c + s_1
            s_2 = np.arcsin((p * np.sin(2*b))/(1-e_earth_sq*np.sin(b)**2)**0.5)
        B = b
        H = D * np.cos(B) + z * np.sin(B) - a * (1 - e_earth_sq * np.sin(B) ** 2) ** 0.5
        return L, B, H        
    