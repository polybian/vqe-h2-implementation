# Simplified explanation

*El algoritmo VQE y la molécula de hidrógeno, explicados sin
matemáticas que asusten. Si puedes explicárselo a alguien sin
formación técnica, puedes explicarlo en una entrevista.*

---

## ⚛️ Capítulo 1: La pregunta de los mil millones

Los químicos quieren saber una cosa de cada molécula: **¿cuánta energía
tiene en su estado más tranquilo?** (su "estado base"). Con ese número
puedes predecir si una reacción ocurre, si un material conduce
electricidad, si una batería funciona.

El problema: calcular esa energía exactamente es brutalmente difícil
para moléculas grandes, el costo crece exponencialmente. Ahí entra la
computación cuántica.

> 🧒 *"¿Y por qué el hidrógeno, la molécula más simple?"*
> Porque es la única donde sabemos la respuesta exacta de antemano.
> Es como practicar con la partitura enfrente: si tu algoritmo no
> reproduce lo conocido, no puedes confiar en él para lo desconocido.

---

## 🎻 Capítulo 2: El afinador (qué es VQE)

Imagina que afinas una guitarra a oído:

1. Tocas una cuerda 🎸
2. Escuchas si suena alta o baja 👂
3. Giras la clavija 🔧
4. Repites hasta que suena perfecta ✨

VQE hace exactamente eso, pero con moléculas:

1. **La computadora cuántica** prepara un "estado de prueba" de la
   molécula (toca la cuerda).
2. **Mide su energía** (escucha).
3. **Una computadora normal** ajusta unas perillas matemáticas θ
   (gira la clavija).
4. Repite hasta que la energía no baja más.

Hay una garantía hermosa: la física cuántica demuestra que la energía
que mides **nunca** puede estar por debajo de la verdadera. Es como
afinar sabiendo que nunca te pasarás de la nota correcta: solo te
acercas desde arriba. A eso le dicen el **principio variacional**.

> 🧒 *"¿Y las perillas θ qué son?"*
> Son los ángulos de unas compuertas cuánticas (rotaciones Ry). El
> circuito de prueba de este proyecto tiene 4 perillas y una capa que
> entrelaza los qubits. Girarlas cambia el estado de prueba hasta
> encontrar el más tranquilo.

---

## 🧱 Capítulo 3: Traducir química a qubits (el mapeo)

Aquí está el truco más bonito del proyecto. Los electrones obedecen
reglas raras (son *fermiones*: no pueden compartir estado). Los qubits
no. Para simular unos con otros hay que **traducir**.

El proyecto implementa el **mapeo de paridad** desde cero: cada qubit
guarda si el número de electrones hasta cierto punto es par o impar.
Y como la molécula conserva el número total de electrones y el de
"ciertos spines", dos qubits quedan fijos por simetría y se eliminan
sin aproximación, de 4 qubits a 2. A eso le dicen *tapering*.

El resultado es el Hamiltoniano de 2 qubits del proyecto:

```
H = −0.332404·II + 0.397937·ZI − 0.397937·IZ − 0.011280·ZZ + 0.180931·XX
```

Cada término es una "pregunta" que le haces al estado cuántico, y la
energía es la suma ponderada de las respuestas.

---

## 📉 Capítulo 4: Estirar la molécula (la curva de disociación)

Si juntas demasiado los dos átomos de hidrógeno, se repelen (energía
alta). Si los separas muchísimo, cada uno se va por su lado (energía
alta también, pero menos). En el medio hay un punto dulce donde la
molécula es más feliz: **~0.74 Å**. Esa es la longitud de enlace real
del H₂.

El script 03 calcula la energía exacta en 13 distancias y dibuja la
curva. El mínimo cae donde la química dice que debe caer. Cuando tu
simulación reproduce la realidad sin que la fuerces, sabes que la
máquina funciona.

---

## ▶️ Cómo correrlo

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
bash run_all.sh
```

Verás el Hamiltoniano con su energía exacta, el VQE convergiendo
iteración por iteración hasta el mismo número, y la curva de
disociación con su mínimo en el lugar correcto. 🎉
