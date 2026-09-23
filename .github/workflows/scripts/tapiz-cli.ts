const K = 6.661232;

const LAMBDA = [
  3.2746,
  3.0331,
  2.9884,
  8.1283
];

const UMBRAL_K = 0.65;

const NORMA = Math.abs(K);
const ESTABLE = NORMA >= UMBRAL_K;

const DECISION = ESTABLE ? "EXECUTE" : "HOLD";
const COLOR = ESTABLE ? "\x1b[32m" : "\x1b[31m";
const RESET = "\x1b[0m";

console.log(`
${COLOR}=== EIGEN TAPIZ CLI v2.0 ===${RESET}

dato
 |
 v
descomposicion
 |
 v
firma
 |
 v
vector
 |
 v
geometria
 |
 v
resonancia
 |
 v
dinamica
 |
 v
decision

[ ||K||    ]: ${NORMA.toFixed(6)}

[ LAMBDA   ]:
  λ1 = ${LAMBDA[0].toFixed(4)}
  λ2 = ${LAMBDA[1].toFixed(4)}
  λ3 = ${LAMBDA[2].toFixed(4)}
  λ4 = ${LAMBDA[3].toFixed(4)}

[ UMBRAL K ]: ${UMBRAL_K}
[ ESTADO   ]: ${COLOR}${ESTABLE ? "ESTRUCTURA ACTIVA" : "ESTRUCTURA DEBIL"}${RESET}
[ DECISION ]: ${COLOR}${DECISION}${RESET}

K = H - H^T
L = iK

Los lambda representan los modos estructurales
del estado actual del Tapiz.
`);
