import numpy as np
import random
import time
from math import sqrt
from tqdm import tqdm
import csv
import os
from dataclasses import dataclass
from typing import Dict, Tuple, List, Optional, Set
import matplotlib.pyplot as plt

# --- Configuração ---
@dataclass
class ConfigEVRP:
    pop_size: int = 100
    mutation_rate: float = 0.30
    crossover_rate: float = 0.90
    elitism_rate: float = 0.10
    tournament_size: int = 5
    max_stagnation: int = 250
    runs: int = 20
    local_search_rate: float = 0.20

# --- Classes InstanciaEVRP e CromossomoEVRP (sem alterações) ---
class InstanciaEVRP:
    def __init__(self, filename: str):
        self.filename = filename
        self.name = os.path.basename(filename)
        self.coords: Dict[int, Tuple[float, float]] = {}
        self.demands: Dict[int, int] = {}
        self.competition_optimal = {'E-n23-k3.evrp': 571.94, 'E-n51-k5.evrp': 529.90}
        self.optimal_value: float = self.competition_optimal.get(self.name, 0.0)
        self.vehicles: int = 0
        self.customers: int = 0
        self.num_stations: int = 0
        self.station_nodes: Set[int] = set()
        self.depot: int = 1
        self.dist_matrix: np.ndarray = np.array([])
        
        self._load_instance()
        self._build_distance_matrix()
    
    def _load_instance(self):
        try:
            with open(self.filename, 'r') as f:
                lines = [line.strip() for line in f if line.strip()]

            section = None
            for line in lines:
                if line.startswith("VEHICLES"): self.vehicles = int(line.split()[-1])
                elif line.startswith("DIMENSION"): self.customers = int(line.split()[-1]) - 1
                elif line.startswith("STATIONS:"): self.num_stations = int(line.split()[-1])
                elif line in ["NODE_COORD_SECTION", "DEMAND_SECTION", "STATIONS_COORD_SECTION", "DEPOT_SECTION"]:
                    section = line
                    continue
                elif line == "EOF": break
                elif section == "NODE_COORD_SECTION":
                    parts = line.split()
                    if len(parts) >= 3:
                        idx, x, y = int(parts[0]), float(parts[1]), float(parts[2])
                        self.coords[idx] = (x, y)
                elif section == "STATIONS_COORD_SECTION":
                    self.station_nodes.add(int(line.split()[0]))
                elif section == "DEPOT_SECTION":
                     if line != "-1": self.depot = int(line)

            print(f"\nInstância carregada: {self.name}")
            print(f"Clientes: {self.customers}, Veículos: {self.vehicles}, Estações: {self.num_stations}")
            print(f"Valor Ótimo da Competição para comparação: {self.optimal_value}")

        except (IOError, ValueError) as e:
            raise ValueError(f"Erro ao carregar ou processar o arquivo da instância '{self.filename}': {e}")
    
    def _build_distance_matrix(self):
        nodes = sorted(self.coords.keys())
        max_node_id = max(nodes)
        self.dist_matrix = np.full((max_node_id + 1, max_node_id + 1), fill_value=np.inf)
        coords_array = np.array([self.coords[i] for i in nodes])
        diff = coords_array[:, np.newaxis, :] - coords_array[np.newaxis, :, :]
        self.dist_matrix[np.ix_(nodes, nodes)] = np.sqrt(np.sum(diff**2, axis=-1))

    def calculate_route_distance(self, route: List[int]) -> float:
        if len(route) < 2: return 0.0
        return self.dist_matrix[route[:-1], route[1:]].sum()

class CromossomoEVRP:
    def __init__(self, solver: 'AlgoritmoGeneticoEVRP', genes: Optional[List[int]] = None):
        self.solver = solver
        self.genes = genes if genes is not None else self._generate_random_genes()
        self.tours: List[List[int]] = []
        self.fitness: float = float('inf')
        self.is_valid: bool = False
        
        self._decode_and_validate()
        self._evaluate()
    
    def _generate_random_genes(self) -> List[int]:
        customers = list(range(2, 2 + self.solver.instance.customers))
        random.shuffle(customers)
        if self.solver.instance.vehicles > 1:
            num_splits = self.solver.instance.vehicles - 1
            if len(customers) > num_splits:
                split_points = sorted(random.sample(range(1, len(customers)), num_splits))
                offset = 0
                for pos in split_points:
                    customers.insert(pos + offset, self.solver.instance.depot)
                    offset += 1
        return customers
    
    def _decode_and_validate(self):
        depot = self.solver.instance.depot
        required_customers = set(range(2, 2 + self.solver.instance.customers))
        
        decoded_tours = []
        current_tour = [depot]
        for gene in self.genes:
            if gene == depot:
                if len(current_tour) > 1:
                    current_tour.append(depot)
                    decoded_tours.append(current_tour)
                current_tour = [depot]
            else:
                current_tour.append(gene)
        
        if len(current_tour) > 1:
            current_tour.append(depot)
            decoded_tours.append(current_tour)
            
        self.tours = decoded_tours
        
        if len(self.tours) > self.solver.instance.vehicles:
            self.is_valid = False
            return
            
        visited_customers = {node for tour in self.tours for node in tour if node != depot}
                    
        if visited_customers != required_customers:
            self.is_valid = False
            return

        self.is_valid = True

    def _evaluate(self):
        if not self.is_valid:
            self.fitness = float('inf')
            return
        self.fitness = sum(self.solver.instance.calculate_route_distance(tour) for tour in self.tours)
    
    def crossover(self, other: 'CromossomoEVRP') -> 'CromossomoEVRP':
        if random.random() >= self.solver.config.crossover_rate:
            return self
        
        p1_genes = [g for g in self.genes if g != self.solver.instance.depot]
        p2_genes = [g for g in other.genes if g != self.solver.instance.depot]

        size = len(p1_genes)
        child_genes = [None] * size
        start, end = sorted(random.sample(range(size), 2))
        
        child_genes[start:end] = p1_genes[start:end]
        segment = set(p1_genes[start:end])
        
        p2_genes_to_fill = [gene for gene in p2_genes if gene not in segment]
        
        fill_idx = 0
        for i in range(size):
            if child_genes[i] is None:
                child_genes[i] = p2_genes_to_fill[fill_idx]
                fill_idx += 1
        
        if self.solver.instance.vehicles > 1:
            num_splits = self.solver.instance.vehicles - 1
            if len(child_genes) > num_splits:
                split_points = sorted(random.sample(range(1, len(child_genes)), num_splits))
                for i, pos in enumerate(split_points):
                    child_genes.insert(pos + i, self.solver.instance.depot)

        return CromossomoEVRP(self.solver, child_genes)

    def mutate(self) -> 'CromossomoEVRP':
        if random.random() >= self.solver.config.mutation_rate:
            return self
        
        mutated_genes = self.genes.copy()
        if len(mutated_genes) < 2: return self

        mutation_type = random.random()
        
        if mutation_type < 0.4: # Swap
            idx1, idx2 = random.sample(range(len(mutated_genes)), 2)
            mutated_genes[idx1], mutated_genes[idx2] = mutated_genes[idx2], mutated_genes[idx1]
        elif mutation_type < 0.8: # Inversion
            start, end = sorted(random.sample(range(len(mutated_genes)), 2))
            if start < end:
                mutated_genes[start:end] = mutated_genes[start:end][::-1]
        else: # Scramble
            start, end = sorted(random.sample(range(len(mutated_genes)), 2))
            if start < end:
                segment = mutated_genes[start:end]
                random.shuffle(segment)
                mutated_genes[start:end] = segment

        return CromossomoEVRP(self.solver, mutated_genes)

    def repair(self) -> 'CromossomoEVRP':
        if self.is_valid: return self
        
        depot = self.solver.instance.depot
        required_customers = set(range(2, 2 + self.solver.instance.customers))
        
        genes = self.genes
        present_customers = []
        seen = set()
        for gene in genes:
            if gene != depot and gene in required_customers and gene not in seen:
                present_customers.append(gene)
                seen.add(gene)
        
        missing = list(required_customers - seen)
        random.shuffle(missing)
        repaired_genes = present_customers + missing
        
        if self.solver.instance.vehicles > 1 and len(repaired_genes) > 1:
            num_splits = self.solver.instance.vehicles - 1
            if len(repaired_genes) > num_splits:
                split_points = sorted(random.sample(range(1, len(repaired_genes)), num_splits))
                for i, pos in enumerate(split_points):
                    repaired_genes.insert(pos + i, depot)
                
        return CromossomoEVRP(self.solver, repaired_genes)
        
    def apply_local_search(self) -> 'CromossomoEVRP':
        if not self.is_valid:
            return self

        new_tours = [self.solver.local_search_2opt(tour) for tour in self.tours]
        
        new_genes = []
        for i, tour in enumerate(new_tours):
            new_genes.extend(tour[1:-1])
            if i < len(new_tours) - 1:
                new_genes.append(self.solver.instance.depot)
        
        return CromossomoEVRP(self.solver, new_genes)

# --- Classe do Algoritmo Memético ---
class AlgoritmoGeneticoEVRP:
    def __init__(self, filename: str, config: Optional[ConfigEVRP] = None):
        self.instance = InstanciaEVRP(filename)
        self.config = config if config else ConfigEVRP()
        self.evaluation_count = 0
        os.makedirs('results', exist_ok=True)
        os.makedirs('plots', exist_ok=True)
    
    def _initialize_population(self) -> List[CromossomoEVRP]:
        population = []
        for _ in range(self.config.pop_size * 20):
            if len(population) >= self.config.pop_size:
                break
            c = CromossomoEVRP(self)
            if c.is_valid:
                population.append(c)
        
        if len(population) < self.config.pop_size:
             raise RuntimeError(f"Não foi possível inicializar a população com o tamanho desejado. Apenas {len(population)} indivíduos válidos foram criados.")
        return population
        
    def _restart_population(self, best_known_solution: CromossomoEVRP) -> List[CromossomoEVRP]:
        print("\n--- REINICIALIZAÇÃO ADAPTATIVA ATIVADA ---")
        new_population = self._initialize_population()
        new_population.sort(key=lambda ind: ind.fitness, reverse=True)
        new_population[0] = best_known_solution
        return new_population

    def local_search_2opt(self, tour: List[int]) -> List[int]:
        if len(tour) <= 4:
            return tour

        best_tour = tour
        improved = True
        while improved:
            improved = False
            for i in range(1, len(best_tour) - 2):
                for j in range(i + 1, len(best_tour)):
                    if j - i == 1: continue
                    
                    current_dist = self.instance.dist_matrix[best_tour[i-1], best_tour[i]] + \
                                   self.instance.dist_matrix[best_tour[j-1], best_tour[j]]
                    
                    new_dist = self.instance.dist_matrix[best_tour[i-1], best_tour[j-1]] + \
                               self.instance.dist_matrix[best_tour[i], best_tour[j]]
                               
                    if new_dist < current_dist:
                        new_tour = best_tour[:i] + best_tour[i:j][::-1] + best_tour[j:]
                        best_tour = new_tour
                        improved = True
                        break
                if improved:
                    break
        return best_tour

    def _select_parent(self, population: List[CromossomoEVRP]) -> CromossomoEVRP:
        tournament_size = min(self.config.tournament_size, len(population))
        tournament = random.sample(population, tournament_size)
        return min(tournament, key=lambda ind: ind.fitness)

    def run(self) -> Dict:
        all_runs_results = []
        n = self.instance.customers + 1 + self.instance.num_stations
        max_evaluations = 25000 * n
        
        for run_num in range(1, self.config.runs + 1):
            print(f"\n--- Execução {run_num}/{self.config.runs} para a instância {self.instance.name} ---")
            print(f"Orçamento de avaliações: {max_evaluations:,}")
            start_time = time.time()
            
            try:
                population = self._initialize_population()
                self.evaluation_count = len(population)
                
                best_overall_solution = min(population, key=lambda ind: ind.fitness)
                stagnation_counter, run_history = 0, []
                
                progress_bar = tqdm(total=max_evaluations, desc="Avaliações", unit="eval", initial=self.evaluation_count, position=0, leave=True)
                
                while self.evaluation_count < max_evaluations:
                    population.sort(key=lambda ind: ind.fitness)
                    
                    elite_size = max(1, int(self.config.elitism_rate * self.config.pop_size))
                    next_generation = population[:elite_size]
                    
                    num_children_to_create = self.config.pop_size - elite_size
                    children = []
                    for _ in range(num_children_to_create):
                        if self.evaluation_count >= max_evaluations: break
                        p1, p2 = self._select_parent(population), self._select_parent(population)
                        
                        child = p1.crossover(p2).mutate()
                        if not child.is_valid: child = child.repair()

                        children.append(child)
                        self.evaluation_count += 1
                        progress_bar.update(1)

                    children.sort(key=lambda c: c.fitness)
                    num_ls = int(len(children) * self.config.local_search_rate)
                    for i in range(num_ls):
                        children[i] = children[i].apply_local_search()

                    next_generation.extend(children)
                    population = next_generation
                    
                    current_best_in_gen = min(population, key=lambda ind: ind.fitness)
                    if current_best_in_gen.fitness < best_overall_solution.fitness:
                        best_overall_solution = current_best_in_gen
                        stagnation_counter = 0
                    else:
                        stagnation_counter += 1
                    
                    run_history.append((self.evaluation_count, best_overall_solution.fitness))
                    
                    if stagnation_counter >= self.config.max_stagnation:
                        population = self._restart_population(best_overall_solution)
                        self.evaluation_count += len(population)
                        progress_bar.update(len(population))
                        stagnation_counter = 0
                
                progress_bar.close()
                exec_time = time.time() - start_time
                gap = (best_overall_solution.fitness - self.instance.optimal_value) / self.instance.optimal_value * 100 if self.instance.optimal_value > 0 else 0
                
                run_result = {'run': run_num, 'fitness': best_overall_solution.fitness, 'gap': gap, 'time': exec_time, 'route': best_overall_solution.tours, 'history': run_history}
                all_runs_results.append(run_result)
                
                self._plot_convergence(run_result)
                self._plot_solution_routes(run_result)
                self._log_best_route_details(run_result)

            except (RuntimeError, ValueError) as e:
                print(f"\nErro crítico na execução {run_num}: {e}")
                continue
        
        # >>>>> ALTERAÇÃO AQUI: Chamando as duas funções de salvamento de arquivos <<<<<
        self._save_summary_results(all_runs_results)
        self._save_detailed_results(all_runs_results) # Salva o novo CSV com as rotas
        
        return self._analyze_final_results(all_runs_results)

    def _log_best_route_details(self, result: Dict, is_final_best: bool = False):
        if is_final_best:
            title = f" MELHORES DETALHES DA ROTA GERAL (Execução {result['run']}) "
            print("\n" + "=" * 25 + title + "=" * 25)
        else:
            print(f"\nMelhor Rota da Execução {result['run']}:")

        print(f"Fitness: {result['fitness']:.4f}")
        if self.instance.optimal_value > 0:
            print(f"Optimal: {self.instance.optimal_value:.2f} (Gap: {result['gap']:.2f}%)")
        
        vehicles_used = len(result['route'])
        vehicles_total = self.instance.vehicles
        print(f"Vehicles used: {vehicles_used}/{vehicles_total}\n")
        print("Route:")
        
        for i, tour in enumerate(result['route'], 1):
            tour_str = ' -> '.join(map(str, tour))
            dist = self.instance.calculate_route_distance(tour)
            print(f"Vehicle {i}: {tour_str} (Distance: {dist:.2f})")
        
        if is_final_best:
            print("=" * (50 + len(title)))

    def _plot_solution_routes(self, result: Dict):
        fig, ax = plt.subplots(figsize=(14, 10))
        customer_nodes = set(range(2, 2 + self.instance.customers))
        cust_coords = np.array([self.instance.coords[c] for c in customer_nodes])
        ax.scatter(cust_coords[:, 0], cust_coords[:, 1], c='skyblue', label='Clientes', s=50, zorder=3)
        if self.instance.station_nodes:
            stat_coords = np.array([self.instance.coords[s] for s in self.instance.station_nodes])
            ax.scatter(stat_coords[:, 0], stat_coords[:, 1], c='lightgreen', marker='s', label='Estações', s=60, zorder=2)
        depot_coord = self.instance.coords[self.instance.depot]
        ax.scatter(depot_coord[0], depot_coord[1], c='red', marker='*', label='Depósito', s=200, zorder=5)

        vehicle_tours = result['route']
        colors = plt.cm.jet(np.linspace(0, 1, len(vehicle_tours)))

        for i, tour in enumerate(vehicle_tours):
            tour_coords = np.array([self.instance.coords[node] for node in tour])
            ax.plot(tour_coords[:, 0], tour_coords[:, 1], color=colors[i], 
                    label=f'Veículo {i+1} (Dist: {self.instance.calculate_route_distance(tour):.2f})',
                    zorder=4, marker='o', markersize=4)
        
        ax.set_title(f"Melhor Rota da Execução {result['run']} - {self.instance.name}\nFitness Total: {result['fitness']:.2f}", fontsize=16)
        ax.set_xlabel("Coordenada X"), ax.set_ylabel("Coordenada Y")
        ax.legend(loc='best'), ax.grid(True, linestyle='--', linewidth=0.5)
        plt.tight_layout()
        filename = f"plots/route_{self.instance.name}_run{result['run']}.png"
        plt.savefig(filename), plt.close(fig)
        print(f"Gráfico da melhor rota salvo em: {filename}")

    def _plot_convergence(self, result: Dict):
        history = result['history']
        evals, fitnesses = zip(*history)
        plt.figure(figsize=(12, 7))
        plt.plot(evals, fitnesses, label='Melhor Fitness por Avaliação')
        if self.instance.optimal_value > 0:
            plt.axhline(y=self.instance.optimal_value, color='r', linestyle='--', label=f'Ótimo Conhecido ({self.instance.optimal_value:.2f})')
        plt.xlabel('Avaliações da Função Objetivo')
        plt.ylabel('Distância Total (Fitness)')
        plt.title(f"Convergência do AM - {self.instance.name} (Execução {result['run']})")
        plt.legend(), plt.grid(True, linestyle='--', linewidth=0.5)
        plt.tight_layout()
        plt.savefig(f"plots/convergence_{self.instance.name}_run{result['run']}.png")
        plt.close()
    
    # Função do CSV de resumo (sem alterações)
    def _save_summary_results(self, results: List[Dict]):
        if not results: return
        filename = f"results/summary_{self.instance.name}.csv"
        with open(filename, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['run', 'fitness', 'gap', 'time'])
            writer.writeheader()
            for res in results:
                writer.writerow({'run': res['run'], 'fitness': f"{res['fitness']:.4f}", 'gap': f"{res['gap']:.2f}", 'time': f"{res['time']:.2f}"})
        print(f"\nResumo estatístico salvo em: {filename}")

    # >>>>> NOVA FUNÇÃO PARA SALVAR O CSV DETALHADO COM AS ROTAS <<<<<
    def _save_detailed_results(self, results: List[Dict]):
        """Salva um CSV com os detalhes completos de cada melhor rota por execução."""
        if not results: return
        filename = f"results/detailed_results_{self.instance.name}.csv"
        fieldnames = ['run', 'fitness', 'gap', 'vehicles_used', 'routes_str']
        
        with open(filename, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for res in results:
                # Formata a lista de tours em uma única string, separada por "|"
                route_str = " | ".join([str(tour) for tour in res['route']])
                row = {
                    'run': res['run'],
                    'fitness': f"{res['fitness']:.4f}",
                    'gap': f"{res['gap']:.2f}",
                    'vehicles_used': len(res['route']),
                    'routes_str': route_str
                }
                writer.writerow(row)
        print(f"Resultados detalhados com as rotas salvos em: {filename}")

    def _analyze_final_results(self, results: List[Dict]) -> Dict:
        if not results: return {}
        
        fitnesses, gaps, times = [r['fitness'] for r in results], [r['gap'] for r in results], [r['time'] for r in results]
        best_run = min(results, key=lambda r: r['fitness'])
        
        self._log_best_route_details(best_run, is_final_best=True)
        
        print("\n" + "="*25 + " Análise Final Consolidada (Algoritmo Memético) " + "="*25)
        print(f"Instância: {self.instance.name}")
        print(f"Execuções: {len(results)}")
        print(f"Melhor fitness (Mínimo): {min(fitnesses):.4f} (Execução {best_run['run']})")
        print(f"Pior fitness (Máximo): {max(fitnesses):.4f}")
        print(f"Média de fitness: {np.mean(fitnesses):.4f}")
        print(f"Desvio Padrão (stdev): {np.std(fitnesses):.4f}")
        print(f"Média de Gap: {np.mean(gaps):.2f}%")
        print(f"Média de tempo/execução: {np.mean(times):.2f}s")
        print("="*85)
        
        return {'best_overall_fitness': best_run['fitness'], 'best_route': best_run['route']}

def main():
    config = ConfigEVRP(runs=20)
    instance_files = ["E-n23-k3.evrp", "E-n51-k5.evrp"] 
    
    for instance_file in instance_files:
        if not os.path.exists(instance_file):
            print(f"AVISO: Arquivo da instância '{instance_file}' não encontrado. Pulando...")
            continue
        try:
            solver = AlgoritmoGeneticoEVRP(filename=instance_file, config=config)
            solver.run()
        except (ValueError, RuntimeError) as e:
            print(f"\nERRO FATAL ao processar {instance_file}: {e}")

if __name__ == "__main__":
    main()