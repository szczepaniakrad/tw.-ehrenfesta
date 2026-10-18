import numpy as np
import matplotlib.pyplot as plt


K = int(input("Podaj liczbę kroków czasowych (K): "))
M = int(input("Podaj liczbę punktów przestrzennych (M): "))
L = float(input("Podaj długość przedziału (L): "))
f = float(input("Podaj parametr potencjału f: "))
sigma = float(input("Podaj wartość sigma: "))
x0 = float(input("Podaj początkową pozycję x0: "))


def V(m, L, M, f):
    return (m * L / M) * f #definicja potencjału

T = np.sqrt((x0)/f) #czas dla teorii klasycznej
dt = T / K 

r = -(1j * T / K) / (2 * (L / M) ** 2) #definicja parametru r

def psi_m(m, L, M, x0, sigma): #definicja funkcji psi wraz z warunkiem brzegowym
    x = m * L / M
    if x <= 0 or x >= L:
        return 0.0
    return np.exp(-((x - x0) ** 2) / (2 * sigma ** 2)) / ((np.pi ** 0.25) * (sigma ** 0.5))

def a_m(m, L, M, T, K, r, f): #definicja a_m
    return (1 - 2 * r + 1j * V(m, L, M, f) * (T / (2 * K)))


def b_m(m, L, M, T, K, r, f): #definicja b_m
    return (1 + 2 * r + 1j * V(m, L, M, f) * (T / (2 * K)))

G = np.zeros((M, M), dtype=complex) # macierze G i H rownania schroedingera
H = np.zeros((M, M), dtype=complex)

for m in range(M): #latwo zauwazyc ze powstale macierze sa trojdiagonalne
    G[m, m] = a_m(m+1, L, M, T, K, r, f)
    H[m, m] = b_m(m+1, L, M, T, K, r, f)
    if m > 0:
        G[m, m-1] = r
        H[m, m-1] = -r
    if m < M - 1:
        G[m, m+1] = r
        H[m, m+1] = -r

U = np.linalg.solve(G, H) # GU=H -> U=G^-1H


psi_0 = np.array([psi_m(m+1, L, M, x0, sigma) for m in range(M)], dtype=complex) #wektor poczatkowy funkcji falowej


psi_all = np.zeros((K+1, M), dtype=complex) #tablica ze wszystkimi wartosciami psi^k_m
psi_all[0, :] = psi_0

for k in range(K):
    psi_all[k+1, :] = U @ psi_all[k, :] #obliczamy nastepne wektory poprzez macierz U


x_m = np.array([m * L / M for m in range(M)], dtype=float) #tablica punktow x_m

y_all = np.zeros((K+1, M), dtype=float) #tablica y_m 

for k in range(K+1):
    for m in range(M):
        y_all[k, m] = np.abs(psi_all[k, m])**2 * x_m[m] #wartosci wszystkich y^k_m


integral_k = np.zeros(K+1, dtype=float) #tablica na calki

for k in range(K+1):
    y = y_all[k, :]

    y_odd = y[1:M-1:2]   # y_2m-1
    y_even = y[2:M-2:2]  # y_2m

    integral = L / (3 * M) * (y[0] + 4 * np.sum(y_odd) + 2 * np.sum(y_even) + y[M-1]) #wzor simpsona
    integral_k[k] = integral #zapisywanie wyniku

    
t_values = np.array([k * dt for k in range(K+1)])

plt.figure(figsize=(10, 6))
plt.plot(t_values, integral_k, marker='o', linestyle='-', color='b')
plt.xlabel('Czas', fontsize=12)
plt.ylabel('Wartość oczekiwana położenia <x>', fontsize=12)
plt.title('Wykres dla teorii kwantowej', fontsize=14)
plt.grid(True)
plt.show()


x_classical = x0 - f * t_values ** 2 #polozenie w klasycznej teorii

plt.figure(figsize=(10, 6))
plt.plot(t_values, x_classical, marker='o', linestyle='-', color='r')
plt.xlabel('Czas', fontsize=12)
plt.ylabel('Położenie x', fontsize=12)
plt.title('Wykres dla teorii klasycznej', fontsize=14)
plt.grid(True)
plt.legend()
plt.show()


plt.figure(figsize=(10, 6))
plt.plot(t_values, x_classical - integral_k, color='purple', linestyle='-', marker='o')
plt.xlabel('Czas', fontsize=12)
plt.ylabel('Różnica położenia x - <x>', fontsize=12)
plt.title('Różnica położeń teorii klasycznej i kwantowej', fontsize=14)
plt.ylim(-5, 5)  
plt.legend()
plt.grid(True)
plt.show()

