# 🚚 Algoritmo Memético para o Problema de Roteamento de Veículos (VRP)

## 🔍 Resumo

Este projeto apresenta a implementação de um **Algoritmo Memético (AM)** de alta performance em Python para resolver o **Problema de Roteamento de Veículos (VRP)**. Um Algoritmo Memético é uma abordagem híbrida que combina a exploração global de um **Algoritmo Genético** com o poder de refinamento de uma **Busca Local**, resultando em soluções de maior qualidade e consistência.

Esta implementação evoluiu de um AG simples para uma solução robusta, incorporando estratégias de **intensificação**, **diversificação adaptativa** e **análise automatizada** de resultados.

---

## 🚀 Principais Características

- ✅ **Abordagem Híbrida**: Combinação de Algoritmo Genético com Busca Local (2-Opt).
- ✅ **Intensificação**: Uso da heurística **2-Opt** para refinar rotas localmente.
- ✅ **Diversificação**: Estratégia de **Reinicialização Adaptativa** contra estagnação.
- ✅ **Representação por Permutação**: Natural e eficiente para o VRP.
- ✅ **Operadores Avançados**:
  - Seleção por torneio
  - Crossover de Ordem (OX1)
  - Mutações: Swap, Inversion, Scramble
- ✅ **Análise Automatizada**: Geração de gráficos, CSVs e relatórios com estatísticas e visualizações das soluções.

---

## 📁 Estrutura do Projeto

```
.
├── E-n23-k3.evrp                # Instância 1 do problema
├── E-n51-k5.evrp                # Instância 2 do problema
├── memetic_solver.py            # Script principal do algoritmo memético
├── plots/                       # Gráficos gerados
│   ├── convergence_.png
│   ├── route_.png
│   └── (outros gráficos)
├── results/                     # Relatórios de execução
│   ├── summary_.csv             # Estatísticas resumidas
│   └── detailed_results_.csv    # Rotas completas por execução
└── README.md                    # Este arquivo
```

---

## ⚙️ Pré-requisitos

- ✅ **Python 3.8 ou superior**  
  🔗 [https://www.python.org/downloads/](https://www.python.org/downloads/)

> ⚠️ No Windows, marque a opção **"Add Python to PATH"** durante a instalação.

---

## 🧪 Como Executar (Passo a Passo)

### 1. Clone ou baixe o repositório

```bash
git clone https://github.com/cristiano-almeida/memetic_solver
cd memetic_solver
```

Ou baixe manualmente via **Code > Download ZIP**

### 2. Crie e ative um ambiente virtual

```bash
# Criar o ambiente
python -m venv venv

# Ativar no Windows
venv\Scripts\activate

# Ativar no Linux/Mac
source venv/bin/activate
```

### 3. Instale as dependências

Crie um arquivo `requirements.txt` com o conteúdo:

```
numpy
matplotlib
tqdm
```

E instale com:

```bash
pip install -r requirements.txt
```

### 4. Execute o algoritmo

```bash
python memetic_solver.py
```

Os resultados serão salvos automaticamente nas pastas `results/` e `plots/`.

---

## 📊 Resultados e Análise

A execução gera:

- 📄 **Saída no console** com o detalhamento da melhor rota por execução
- 📈 `summary_*.csv`: Estatísticas de desempenho (média, desvio padrão, melhor e pior solução)
- 🗺️ `detailed_results_*.csv`: Rotas completas por execução
- 📊 Gráficos `.png` em `plots/`:
  - Evolução da convergência do algoritmo
  - Visualização das melhores rotas geradas

O algoritmo foi testado para extrair o máximo do orçamento computacional, com ganhos notáveis em robustez e qualidade da solução comparado ao AG simples.

---

## 📚 Referências e Conceitos Aplicados

- **Algoritmos Meméticos**: Híbridos propostos por P. Moscato (1989), combinando algoritmos genéticos com refinamento local.
- **Busca Local 2-Opt**: Heurística eficiente para o TSP, aplicada aqui à etapa de intensificação.
- **Reinicialização Adaptativa**: Estratégia para escapar de ótimos locais e manter a diversidade da população.

---

🔧 Projeto desenvolvido para fins de estudo e experimentação em **otimização de rotas com algoritmos bio-inspirados**.

