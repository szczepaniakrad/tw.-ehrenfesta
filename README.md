# Numeryczna weryfikacja twierdzenia Ehrenfesta

Projekt pokazuje numerycznie, że wartość oczekiwana położenia cząstki kwantowej, wyznaczona z równania Schrödingera, podąża za trajektorią klasyczną. Porównujemy je dla cząstki w jednorodnym polu o stałej sile (potencjał liniowy) w jednym wymiarze.

## Idea

Twierdzenie Ehrenfesta mówi, że wartości oczekiwane położenia i pędu spełniają równania zbliżone do klasycznych równań ruchu:

```math
\frac{d\langle x\rangle}{dt} = \frac{\langle p\rangle}{m}, \qquad \frac{d\langle p\rangle}{dt} = -\left\langle \frac{dV}{dx}\right\rangle
```

Dla potencjału liniowego $V(x) = f\,x$ siła $-dV/dx = -f$ jest stała, więc $\langle dV/dx\rangle = dV/d\langle x\rangle$ i równania są identyczne z klasycznymi. Wartość oczekiwana $\langle x\rangle(t)$ powinna więc pokrywać się z klasyczną trajektorią $x(t)$.

## Model

**Jednostki:** $\hbar = 1$, $2m = 1$ (hamiltonian ma postać $\hat H = -\partial_x^2 + V(x)$).

**Równanie Schrödingera:**

```math
i\,\frac{\partial \psi}{\partial t} = -\frac{\partial^2 \psi}{\partial x^2} + f\,x\,\psi
```

**Potencjał:** $V(x) = f\,x$.

**Warunek początkowy:** pakiet gaussowski o zerowym średnim pędzie

```math
\psi(x,0) = \frac{1}{\pi^{1/4}\sqrt{\sigma}}\exp\left(-\frac{(x-x_0)^2}{2\sigma^2}\right)
```

**Warunki brzegowe:** $\psi(0,t) = \psi(L,t) = 0$ (nieskończone ściany na końcach przedziału $[0, L]$).

**Rozwiązanie klasyczne** (cząstka o masie $m = 1/2$, siła $-f$, start w spoczynku z $x_0$):

```math
x_{kl}(t) = x_0 - f\,t^2
```

Czas dojścia do $x = 0$ wynosi $T = \sqrt{x_0/f}$. Symulacja trwa dokładnie tyle.

## Metoda numeryczna

### 1. Siatka

Przedział $[0, L]$ dzielimy na $M$ punktów, $x_m = mL/M$. Czas $[0, T]$ dzielimy na $K$ kroków, $\Delta t = T/K$.

### 2. Schemat Cranka–Nicolson

Równanie Schrödingera dyskretyzujemy tak, że w każdym kroku czasowym spełniony jest układ

```math
G\,\psi^{k+1} = H\,\psi^{k}
```

gdzie $G$ i $H$ są macierzami trójdiagonalnymi z elementami:

- na przekątnej: $a_m = 1 - 2r + \frac{i}{2}V_m\Delta t$ (dla $G$) oraz $b_m = 1 + 2r + \frac{i}{2}V_m\Delta t$ (dla $H$),
- poza przekątną: $r$ (dla $G$) oraz $-r$ (dla $H$),
- parametr $r$ dany wzorem:

```math
r = -\frac{i\,\Delta t}{2(\Delta x)^2}
```

### 3. Macierz ewolucji

$U = G^{-1}H$, a kolejne kroki liczymy jako $\psi^{k+1} = U\psi^{k}$.

### 4. Wartość oczekiwana położenia

```math
\langle x\rangle_k = \int_0^L |\psi^k(x)|^2\,x\,dx
```

Całkę obliczamy metodą Simpsona dla każdego kroku czasowego.

### 5. Porównanie

Porównujemy $\langle x\rangle(t)$ z $x_{kl}(t)$ i rysujemy wykres różnicy $x_{kl}(t) - \langle x\rangle(t)$.

## Parametry wejściowe

Program poprosi o podanie parametrów:

| Parametr | Znaczenie |
|----------|-----------|
| `K` | liczba kroków czasowych |
| `M` | liczba punktów przestrzennych |
| `L` | długość przedziału $[0, L]$ |
| `f` | parametr potencjału $V(x) = f\,x$ |
| `sigma` | szerokość początkowego pakietu gaussowskiego |
| `x0` | początkowe położenie środka pakietu |

Przykładowy zestaw parametrów (do dostosowania): `K = 1000`, `M = 500`, `L = 20`, `f = 1`, `sigma = 1`, `x0 = 10`.

**Wskazówki doboru parametrów:**

- $\sigma \ll x_0$ oraz $x_0 - 4\sigma > 0$ i $x_0 + 4\sigma < L$, aby pakiet na starcie nie dotykał ścian,
- $\Delta x = L/M$ powinno być wyraźnie mniejsze niż $\sigma$,
- większe $K$ i $M$ dają dokładniejszy wynik, ale dłuższy czas obliczeń (macierz $U$ jest gęsta, a układ $GU = H$ rozwiązujemy bezpośrednio).

## Wyniki

Program generuje trzy wykresy:

1. **Teoria kwantowa**: $\langle x\rangle(t)$ otrzymane z rozwiązania równania Schrödingera.
2. **Teoria klasyczna**: $x_{kl}(t) = x_0 - f\,t^2$.
3. **Różnica** $x_{kl}(t) - \langle x\rangle(t)$.

Jeśli twierdzenie Ehrenfesta jest spełnione, wykresy 1 i 2 się pokrywają, a różnica jest bliska zeru i zbliża się do niego wraz ze wzrostem $K$ i $M$.

## Ograniczenia

- Na końcach przedziału narzucone są ściany nieskończone. Gdy pakiet zbliża się do $x = 0$ (końcówka symulacji), odbicie od ściany i rozpływanie się pakietu mogą powodować odchylenia od trajektorii klasycznej. Nie jest to błąd twierdzenia, tylko skutek wpływu ściany.
- Dla potencjałów nieliniowych $\langle V'(x)\rangle \neq V'(\langle x\rangle)$, więc zgodność z mechaniką klasyczną jest tylko przybliżona.
