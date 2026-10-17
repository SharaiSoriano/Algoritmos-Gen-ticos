import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import random
import copy
from typing import List, Dict

from CromosomaBin import CromosomaBinario
from interfaz import AptitudMochila, SeleccionTorneo, ReproduccionUnCruce, MutacionFlip

class GA(object):
    def __init__(self, num_inv: int, num_genes: int, hyperparams: Dict):
        self.poblacion: List[CromosomaBinario] = [CromosomaBinario(num_genes) for i in range(num_inv)]
        self.num_inv: int = num_inv
        self.hyperparams: Dict = hyperparams
        self.evaluar = hyperparams["Aptitud"]
        self.seleccion = hyperparams["Seleccion"]
        self.reproducion = hyperparams["Reproduccion"]
        self.mutacion = hyperparams["Mutacion"]
        self.num_remplazo: int = int(hyperparams["por_remplazo"] * num_inv)

        #inicializar
        for i in range(num_inv):
            self.poblacion[i].llenar_random_data()

    def buscar_solucion(self, tolerancia: float, var: Dict):
        avg_tol: float = 99999999.99
        sum_score_prev: float = 0.0
        generacion: int = 1

        fail: bool = False
        fail_count: int = 0
        fail_score: float = 0.0

        while not fail:
            self.evaluar.evaluar_poblacion(self.poblacion, var)
            sum_scores: float = sum(cromosoma.score for cromosoma in self.poblacion)
            sum_scores /= self.num_inv

            if sum_scores != 0:
                avg_tol = abs(sum_score_prev - sum_scores) / sum_scores

            if (generacion + 1) % 100 == 0:
                print(f"[{generacion:05d}] avg_err : {avg_tol:.5f} avg_score: {sum_scores:.5f}")

            if fail_score == avg_tol:
                fail_count += 1

            if fail_count >= 10:
                fail = True
                break

            if avg_tol < tolerancia:
                break

            fail_score = avg_tol
            sum_score_prev = sum_scores

            self.poblacion = sorted(self.poblacion, key=lambda c: c.score, reverse=True)
            mating_pool = self.seleccion.seleccion_poblacion(self.poblacion)
            generacion_nueva = copy.deepcopy(self.poblacion[: self.num_remplazo])

            self.reproducion.reproduccion_poblacion(generacion_nueva, mating_pool, self.num_inv)
            self.mutacion.mutacion_poblacion(generacion_nueva)

            generacion += 1
            self.poblacion = copy.deepcopy(generacion_nueva)


def sum_pesos_test(pesos, cromosoma):
    s = 0
    for i in range(cromosoma.n_genes):
        if cromosoma.get_bits(i + 1) == 1:
            s += pesos[i]
    return s


if __name__ == "__main__":
    variables = {}
    variables["pesos"] = [random.randint(5, 25) for _ in range(32)]
    variables["valores"] = [random.randint(50, 300) for _ in range(32)]
    variables["capacidad"] = 150  # kg

    num_genes = 32
    cant_ind = 100      

    gaparams = {}
    gaparams["Aptitud"] = AptitudMochila()  
    gaparams["Seleccion"] = SeleccionTorneo(3)
    gaparams["Reproduccion"] = ReproduccionUnCruce(num_genes)
    gaparams["Mutacion"] = MutacionFlip(0.05)
    gaparams["por_remplazo"] = 0.5

    algoritmoGA = GA(cant_ind, num_genes, gaparams)
    algoritmoGA.buscar_solucion(1e-4, variables)

    solucion = sorted(algoritmoGA.poblacion, key=lambda c: c.score, reverse=True)
    
    print(f"{'id':<7} {'Cromosoma':<34} {'Score':<5}")
    for i, cromosoma in enumerate(solucion[:15]):
        print(f"[{i + 1:>2} ] {cromosoma.bits:032b} : {cromosoma.score:.0f} -> {sum_pesos_test(variables['pesos'], cromosoma):>2} kg")