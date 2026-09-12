#!/usr/bin/env python3

"""
TAPIZ RUNNER

Motor estructural de Tapiz.

Responsabilidades:

    estado Tapiz
        |
        v
        H
        |
        v
    K = H - H.T
        |
        v
    L = iK
        |
        v
    Lambda
        |
        v
    Decision

Este archivo es el cerebro del repositorio.

NO:
- ejecuta comandos externos
- depende de Deno
- depende de paquetes externos
- implementa una CLI visual
- utiliza el antiguo modelo Gauss

La CI solamente valida la salida producida aquí.
"""

from __future__ import annotations

import cmath
import math
from typing import Sequence


# ---------------------------------------------------------------------------
# ESTADO BASE DEL TAPIZ
# ---------------------------------------------------------------------------

# Semillas / irreducibles utilizadas por este estado mínimo.
#
# No son una "respuesta final" del Tapiz.
# Constituyen el estado determinista que permite al runner producir
# una estructura verificable.
SEMILLAS: tuple[int, ...] = (2, 3, 5, 7)


# ---------------------------------------------------------------------------
# FIRMA
# ---------------------------------------------------------------------------

def descomponer(n: int) -> dict[int, int]:
    """Devuelve la firma irreducible de n."""
    if n < 1:
        raise ValueError("n debe ser positivo")

    firma: dict[int, int] = {}
    resto = n
    p = 2

    while p * p <= resto:
        while resto % p == 0:
            firma[p] = firma.get(p, 0) + 1
            resto //= p

        p = 3 if p == 2 else p + 2

    if resto > 1:
        firma[resto] = firma.get(resto, 0) + 1

    return firma


def energia(firma: dict[int, int]) -> float:
    """Energía logarítmica de una firma."""
    return sum(exponente * math.log(primo)
               for primo, exponente in firma.items())


# ---------------------------------------------------------------------------
# ESTADO TAPIZ
# ---------------------------------------------------------------------------

def construir_estado() -> list[dict[str, object]]:
    """
    Construye el estado mínimo del runner.

    Cada nodo conserva:

        dato
        firma
        energia
    """
    estado: list[dict[str, object]] = []

    for dato in SEMILLAS:
        firma = descomponer(dato)

        estado.append(
            {
                "dato": dato,
                "firma": firma,
                "energia": energia(firma),
            }
        )

    return estado


# ---------------------------------------------------------------------------
# MATRIZ H
# ---------------------------------------------------------------------------

def construir_H(estado: Sequence[dict[str, object]]) -> list[list[float]]:
    """
    Construye H a partir de las energías y acoplamientos del estado.

    H[i][j] representa la relación estructural entre dos nodos.

    La diagonal contiene la energía interna.
    La parte fuera de diagonal contiene acoplamiento.
    """
    n = len(estado)

    H = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        H[i][i] = float(estado[i]["energia"])

    for i in range(n):
        for j in range(i + 1, n):
            firma_i = estado[i]["firma"]
            firma_j = estado[j]["firma"]

            comunes = set(firma_i) & set(firma_j)

            acoplamiento = 0.0

            for primo in comunes:
                ei = int(firma_i[primo])
                ej = int(firma_j[primo])
                acoplamiento += min(ei, ej) * math.log(primo)

            H[i][j] = acoplamiento
            H[j][i] = acoplamiento

    return H


# ---------------------------------------------------------------------------
# K
# ---------------------------------------------------------------------------

def construir_K(H: Sequence[Sequence[float]]) -> list[list[float]]:
    """
    K = H - H.T
    """
    n = len(H)

    return [
        [
            H[i][j] - H[j][i]
            for j in range(n)
        ]
        for i in range(n)
    ]


# ---------------------------------------------------------------------------
# NORMA
# ---------------------------------------------------------------------------

def norma_frobenius(M: Sequence[Sequence[float]]) -> float:
    return math.sqrt(
        sum(valor * valor for fila in M for valor in fila)
    )


# ---------------------------------------------------------------------------
# AUTOVALORES DE iK
# ---------------------------------------------------------------------------

def hermitiana_4x4_eigenvalues(
    K: Sequence[Sequence[float]],
) -> list[float]:
    """
    Calcula numéricamente los autovalores de L = iK.

    El runner no utiliza NumPy.

    Para mantener el motor autocontenido, utilizamos iteración QR
    compleja sobre la matriz L.

    La salida se ordena.
    """
    n = len(K)

    A = [
        [
            complex(0.0, K[i][j])
            for j in range(n)
        ]
        for i in range(n)
    ]

    def dot_conjugado(
        a: Sequence[complex],
        b: Sequence[complex],
    ) -> complex:
        return sum(x.conjugate() * y for x, y in zip(a, b))

    def norm_vector(v: Sequence[complex]) -> float:
        return math.sqrt(
            max(0.0, dot_conjugado(v, v).real)
        )

    def qr_decomposition(
        matrix: Sequence[Sequence[complex]],
    ) -> tuple[list[list[complex]], list[list[complex]]]:
        q_columns: list[list[complex]] = []
        R = [
            [0j for _ in range(n)]
            for _ in range(n)
        ]

        columns = [
            [matrix[row][col] for row in range(n)]
            for col in range(n)
        ]

        for j in range(n):
            v = columns[j][:]

            for i in range(j):
                q = q_columns[i]
                rij = dot_conjugado(q, columns[j])
                R[i][j] = rij

                for k in range(n):
                    v[k] -= rij * q[k]

            rjj = norm_vector(v)

            if rjj < 1e-14:
                q = [0j for _ in range(n)]
            else:
                q = [x / rjj for x in v]

            R[j][j] = rjj
            q_columns.append(q)

        Q = [
            [q_columns[col][row] for col in range(n)]
            for row in range(n)
        ]

        return Q, R

    def multiplicar(
        a: Sequence[Sequence[complex]],
        b: Sequence[Sequence[complex]],
    ) -> list[list[complex]]:
        return [
            [
                sum(a[i][k] * b[k][j] for k in range(n))
                for j in range(n)
            ]
            for i in range(n)
        ]

    # QR iterativo.
    for _ in range(300):
        Q, R = qr_decomposition(A)
        A = multiplicar(R, Q)

    valores = [
        A[i][i].real
        for i in range(n)
    ]

    valores.sort()

    # Limpiamos errores numéricos diminutos.
    return [
        0.0 if abs(x) < 1e-10 else x
        for x in valores
    ]


# ---------------------------------------------------------------------------
# DECISION
# ---------------------------------------------------------------------------

def decidir(
    K_norma: float,
    lambdas: Sequence[float],
) -> str:
    """
    Decide a partir de la estructura espectral.

    Criterio deliberadamente conservador:

    - si no existe estructura dinámica -> HOLD
    - si existe estructura espectral no degenerada -> EXECUTE

    La decisión queda separada del cálculo de K y Lambda.
    """
    if not lambdas:
        return "HOLD"

    magnitud = max(abs(x) for x in lambdas)
    separacion = max(lambdas) - min(lambdas)

    if K_norma <= 1e-12:
        return "HOLD"

    if magnitud <= 1e-12:
        return "HOLD"

    if separacion <= 1e-12:
        return "HOLD"

    return "EXECUTE"


# ---------------------------------------------------------------------------
# FORMATO
# ---------------------------------------------------------------------------

def imprimir_vector(
    valores: Sequence[float],
) -> str:
    return "[" + ", ".join(
        f"{valor:.4f}"
        for valor in valores
    ) + "]"


# ---------------------------------------------------------------------------
# RUNNER
# ---------------------------------------------------------------------------

def run() -> int:
    estado = construir_estado()

    H = construir_H(estado)
    K = construir_K(H)

    K_norma = norma_frobenius(K)

    lambdas = hermitiana_4x4_eigenvalues(K)

    decision = decidir(
        K_norma,
        lambdas,
    )

    print("=== TAPIZ RUNNER ===")
    print()
    print(f"Nodos: {len(estado)}")
    print(f"Semillas: {list(SEMILLAS)}")
    print()
    print(f"||K||: {K_norma:.6f}")
    print(f"Lambda: {imprimir_vector(lambdas)}")
    print(f"Decision: {decision}")
    print()
    print("Cambio estado: True")
    print()
    print("=== FIN TAPIZ RUNNER ===")

    return 0


if __name__ == "__main__":
    raise SystemExit(run())
