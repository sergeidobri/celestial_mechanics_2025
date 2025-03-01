import numpy as np
import csv
import matplotlib.pyplot as plt
from utils import find_eccentric_anomaly, find_geodetic_coordinates, find_true_anomaly
from constants import *


# 1. Рассчет координат в АГЭСК (x_a, y_a, z_a):
apocentre_radius = APOCENTRE_HEIGHT + EARTH_RADIUS      # км
pericentre_radius = PERICENTRE_HEIGHT + EARTH_RADIUS    # км

semi_major_axis = (apocentre_radius + pericentre_radius) / 2                    # км 
eccentricity = (apocentre_radius - pericentre_radius) / (2 * semi_major_axis)   # безразмерная величина

eccentric_anomaly = find_eccentric_anomaly(MEDIUM_ANOMALY, eccentricity)            # радиан
true_anomaly = find_true_anomaly(eccentric_anomaly, eccentricity, MEDIUM_ANOMALY)   # радиан

latitude_argument = true_anomaly + PERICENTRE_ARGUMENT                              # радиан
radius_vector = semi_major_axis * (1 - eccentricity * np.cos(eccentric_anomaly))    # км

x_a = radius_vector * (  # км
    np.cos(latitude_argument) * np.cos(ASC_NODE_LONGITUDE) - 
    np.sin(latitude_argument) * np.sin(ASC_NODE_LONGITUDE) * np.cos(INCLINATION)
)
y_a = radius_vector * (  # км
    np.cos(latitude_argument) * np.sin(ASC_NODE_LONGITUDE) + 
    np.sin(latitude_argument) * np.cos(ASC_NODE_LONGITUDE) * np.cos(INCLINATION)
)
z_a = radius_vector * np.sin(latitude_argument) * np.sin(INCLINATION)  # км

AGECS_coords = np.matrix([  # координаты в АГЭСК
    [x_a], 
    [y_a], 
    [z_a]
])

# 2. Рассчет трансверсальной и радиальной скоростей (v_t, v_r):
orbital_parameter = semi_major_axis * (1 - eccentricity ** 2)  # км

v_r = (GRAVITATIONAL_PARAMETER / orbital_parameter) ** 0.5 * eccentricity * np.sin(true_anomaly)  # км/c
v_t = (GRAVITATIONAL_PARAMETER / orbital_parameter) ** 0.5 * (1 + eccentricity * np.cos(true_anomaly))  # км/c
v = (v_r**2 + v_t**2) ** 0.5  # км/c

# 3. Рассчет положения точки в ГСК (x, y, z) и рассчитать геодезические координаты (L, B, H):
angular_spacecraft_velocity = (GRAVITATIONAL_PARAMETER / (semi_major_axis ** 3)) ** 0.5  # радиан/с
time_moment = MEDIUM_ANOMALY / angular_spacecraft_velocity  # c

sidereal_time = ANGULAR_EARTH_VELOCITY * time_moment  # рад

rotation_matrix = np.matrix([
    [np.cos(sidereal_time), np.sin(sidereal_time), 0],
    [-np.sin(sidereal_time), np.cos(sidereal_time), 0],
    [0, 0, 1]
])

GCS_coords = rotation_matrix @ AGECS_coords  # координаты в ГСК
x, y, z = GCS_coords.item(0), GCS_coords.item(1), GCS_coords.item(2)

L, B, H = geodetic_coords = find_geodetic_coordinates(  # геодезические координаты
                                    x, y, z, GEA_SEMI_MAJOR_AXIS, e_EARTH_SQUARED)

# 4. Рассчёт плотности атмосферы. Высота космического аппарата - 268.2 км, что соответствует первому высотному диапазону
densities = []

with open("data/coefficients_of_atmosphere_density_model.csv") as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',')
    for row in reader:
        force  = float(row['force'])
        a0     = float(row['a0'])
        a1     = float(row['a1'])
        a2     = float(row['a2'])
        a3     = float(row['a3'])
        a4     = float(row['a4'])
        a5     = float(row['a5'])
        a6     = float(row['a6'])
        
        degree = a0 + a1*H + a2*H**2 + a3*H**3 + a4*H**4 + a5*H**5 + a6*H**6
        atmosphere_density = NIGHT_ATMOSPHERE_DENSITY_120 * np.e ** degree * 1e9  # кг / км^3
        densities.append((force, atmosphere_density))
    
force_array = []
S_array = []
T_array = []
W_array = []
w_array = []

for force, density in densities:
    ballistic_coefficient = DRAG_FORCE_COEFFICIENT * SPACECRAFT_AREA / (2 * SPACECRAFT_MASS)
    S = -ballistic_coefficient * density * v * v_r
    T = -ballistic_coefficient * density * v * v_t
    W = 0
    force_array.append(force)
    S_array.append(S)
    T_array.append(T)
    W_array.append(W)
    w_array.append((S**2 + T**2 + W**2) ** 0.5)

g = GRAVITATIONAL_CONSTANT * EARTH_MASS / ((EARTH_RADIUS + H) * 1e3)**2

# 5. Соответствующий форматированный вывод результатов
print(f"\t---Орбита---")
print(f"Большая полуось: {semi_major_axis:.2f} км")
print(f"Эксцентриситет: {eccentricity:.6f}")
print(f"Орбитальный параметр: {orbital_parameter:.2f} км")

print(f"\n ---Положение космического аппарата---")
print(f"Момент времени после прохождения аппаратом перицентра: {time_moment:.2f} с")
print(f"Эксцентрическая аномалия: {eccentric_anomaly:.4f} радиан")
print(f"Истинная аномалия: {true_anomaly:.4f} радиан")
print(f"Координаты в АГЭСК [x_a, y_a, z_a]: ({x_a:.2f} км, {y_a:.2f} км, {z_a:.2f} км)")
print(f"Координаты в ГСК [x, y, z]: ({x:.2f} км, {y:.2f} км, {z:.2f} км)")
print(f"Геодезические координаты [L, B, H]: ({L:.2f} рад, {B:.2f} рад, {H:.2f} км)")

print(f"\n ---Скорость космического аппарата---")
print(f"Радиальная состовляющая скорости: {v_r:.4f} км/c")
print(f"Трансверсальная состовляющая скорости: {v_t:.4f} км/c")
print(f"Модуль скорости: {v:.4f} км/c")
print(f"Гравитационное ускорение космического аппарата: {g:.4f} м/c^2")

print(f"\n\t---Атмосфера---")
print(f"При средней аномалии:", "\033[33m{}\033[0m".format(round(MEDIUM_ANOMALY * 180 / np.pi)) + " градусов")
print(f"{'Ур-нь солн. акт-ности, 10^-22 Вт/м^2*Гц':>40s} | {'Плотность атмосферы, кг/км^3'}")
print("-"*41 + "|" + "-"*30)
for force, density in densities:
    print(f"{force:>40} | {density:.3f}")
print("-"*71)

print(f"\n\t---Возмущения---")
print(f"При средней аномалии:", "\033[33m{}\033[0m".format(round(MEDIUM_ANOMALY * 180 / np.pi)) + " градусов")
print(f"{'Ур-нь солн. акт-ности, 10^-22 Вт/м^2*Гц':>40s} | {'S, км/c^2':>10s} | {'T, км/c^2':>10s} | {'W, км/c^2':>10s} | {'Полное возмущение, км/c^2':>10s}")
print("-"*41 + "|" + "-"*12 + "|" +  "-"*12 + "|" + "-"*12 + "|" + "-"*30)
for force, S, T, W in zip(force_array, S_array, T_array, W_array):
    print(f"{force:>40} | {S:>10.3e} | {T:>10.3e} | {W:>10.3e} | {(S**2 + T**2 + W**2)**0.5:>10.3e}")
print("-"*110)

# 6. Изображение графиков
fig, (ax_top, ax_bot) = plt.subplots(nrows=2, ncols=1, figsize=(8, 10))

ax_top.grid(True)
ax_bot.grid(True)

ax_top.set_xlabel("Уровень солнечной активности F₀, 10⁻²² Вт/(м²·Гц)")
ax_top.set_ylabel("Компоненты возмущения, км/c²")
ax_top.plot(force_array, S_array, color='yellow', label="S(F₀)")
ax_top.plot(force_array, T_array, color='blue', label="T(F₀)")
ax_top.plot(force_array, W_array, color='green', label="W(F₀)")
ax_top.legend()

ax_bot.set_xlabel("Уровень солнечной активности F₀, 10⁻²² Вт/(м²·Гц)")
ax_bot.set_ylabel("Полное возмущение, км/c²")
ax_bot.plot(force_array, w_array, color='red', label="w(F₀)")
ax_bot.legend()

plt.tight_layout()
plt.show()
