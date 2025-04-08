#Importação de bibliotecas necessárias  
import subprocess
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline

def read_data(filename):
    """
    Carrega os dados de um arquivo CSV e retorna um DataFrame com as colunas "E", "P" e "test_value".
    Também extrai o valor do parâmetro de teste (armazenado na coluna "test_value") e o retorna separadamente
    """
    data = pd.read_csv(filename, delim_whitespace=True, header=None, names=["E", "P", "test_value"])
    test_value = data["test_value"].iloc[0]
    return data, test_value

def process_files_by_prefix(directory, prefix):
    """
    Processa todos os arquivos em um diretório que começam com um prefixo específico, carrega os dados de cada arquivo
    usando a função "read_data", e organiza os resultados em duas listas: uma com os DataFrames e outra com os valores
    de teste extraídos de cada arquivo.
    """
    files = sorted(
        [file for file in os.listdir(directory) if file.startswith(prefix)],
        key=lambda x: int(x.split("_file")[1].split(".")[0])
    )
    data_frames = []
    free_pars = []
    for file in files:
        file_path = os.path.join(directory, file)
        data, test_value = read_data(file_path)
        data_frames.append(data)
        free_pars.append(test_value)
    return data_frames, free_pars
    
def chi_squared(directory, prefixes):
    """
    Calcula o valor de chi-quadrado (χ²) para cada conjunto de dados observado e teórico, 
    e retorna um DataFrame contendo os parâmetros de teste e os valores de χ² para cada canal.
    """

    # Dicionários para armazenar os resultados
    chi2_channels = {prefix: [] for prefix in prefixes}
    all_data_frames = []  # Lista para armazenar todos os DataFrames

    # Para cada prefixo
    for prefix in prefixes:
        free_pars = []
        # Arquivo teórico correspondente
        # true_file = f"prob/data/a/{prefix.split('_file')[0]}_STD_bsm.dat"
        true_file = f"prob/data/a/{prefix.split('_file')[0]}_STD_bsm.dat"
        data_true = pd.read_csv(true_file, delim_whitespace=True, header=None, names=["E", "P"])

        # Processar arquivos do canal atual
        data_frames, free_pars_channel = process_files_by_prefix(directory, prefix)

        # Acumular os valores de free_pars para garantir que sejam únicos
        for df in data_frames:
            free_pars.append(df["test_value"].iloc[0])  # Usando o primeiro valor da coluna "test_value"

        """
        # Debug: Verificar os valores de free_pars
        print(f"Valores de 'free_pars' para o prefixo '{prefix}':")
        print(free_pars)  # Imprime todos os valores de free_pars acumulados
        """

        # Calcular χ² para cada DataFrame
        for df in data_frames:
            observed = df["P"]  # Valores observados
            expected = data_true["P"]  # Valores esperados

            # Calcular χ²
            chi_square = np.sum(((observed - expected) ** 2) / expected)
            # chi_square = np.sum(((observed - expected) ** 2) / expected * 0.1)

            
            # Adicionar o valor de χ² à lista do canal correspondente
            chi2_channels[prefix].append(chi_square)
        
        all_data_frames.append((data_frames, free_pars_channel))  # Armazena os DataFrames e os parâmetros livres

        """
        # Debug: Verificar os valores de χ² calculados para o prefixo
        print(f"Valores de χ² para o prefixo '{prefix}':")
        print(chi2_channels[prefix])  # Imprime os valores de χ² para o prefixo atual
        """

    """
    # Verificar se o número de elementos nas listas é o mesmo
    if len(free_pars) != len(chi2_channels[prefixes[0]]):
        print(f"Erro: o número de valores de free_pars ({len(free_pars)}) não é igual ao número de valores de χ² para o prefixo '{prefixes[0]}' ({len(chi2_channels[prefixes[0]])})")
        return None  # Ou você pode lançar um erro, dependendo de como quiser lidar com isso.
    """
    
    # Criar o DataFrame unificado
    data = {"free_pars": free_pars}
    for prefix in prefixes:
        data[prefix] = chi2_channels[prefix]
        
    # Adicionar a coluna chi_total com a soma de todos os chi-quadrados
    df = pd.DataFrame(data)
    df["chi_total"] = df[prefixes].sum(axis=1)
    
    # Encontrar o menor valor de chi² total
    chi_minimo = df["chi_total"].min()

    # Adicionar a coluna delta_chi
    df["delta_chi"] = df["chi_total"] - chi_minimo

    # Criar e retornar o DataFrame
    return df, all_data_frames
    
def plot_probability(data_frames, free_pars_channel, ylabel, test_par):
    """
    Cria um gráfico de probabilidade (P por E) para todos os DataFrames fornecidos, com legendas associadas aos valores de "test_value".
    """
    fig, ax = plt.subplots(figsize=(9, 6))

    # Plota os dados para cada parâmetro livre
    ax.plot(data_true["E"], data_true["P"], color="black", label="True (SI)", linewidth=2)
    for data, test_value in zip(data_frames, free_pars_channel):
        ax.plot(data["E"], data["P"], linestyle='--', label=f"{test_par} = {test_value:.2e}", linewidth=2)

    # Configurações do gráfico
    ax.set_xlabel("Energy [GeV]", fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.legend(fontsize=16)
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=16)

    # Exibir gráfico
    plt.show()

def plot_chi_squared_by_free_par(chi2_df, prefixes, test_par):
    """
    Cria gráficos de χ² por parâmetro livre para cada prefixo fornecido.
    """
    sigma_levels = {"1σ": 1.00, "95% C.L.": 3.84, "3σ": 9.00, "5σ": 25.00}

    for prefix in prefixes:
        fig, ax = plt.subplots(figsize=(8, 6))

        # Valores de parâmetro livre e χ² para o prefixo atual
        free_pars = chi2_df["free_pars"]
        chi2_values = chi2_df[prefix]

        # Plotar χ² por parâmetro livre
        ax.plot(free_pars, chi2_values, color='red', linestyle='-', linewidth=2)

        # Adicionar linhas de nível de confiança
        for label, value in sigma_levels.items():
            ax.axhline(y=value, color='gray', linestyle='--', linewidth=1)
            ax.text(x=free_pars.iloc[-1], y=value, s=label, fontsize=12, color='gray', ha='left', va='bottom')

        # Configurações do gráfico
        ax.set_xlabel(f"{test_par}", fontsize=14)
        ax.set_ylabel(r"$\chi^2$", fontsize=14)
        ax.set_title(f"$\chi^2$ por {test_par} - {prefix}", fontsize=16)
        ax.grid(True)
        ax.set_xscale("log")
        ax.tick_params(axis='both', which='major', labelsize=10)
        plt.ylim(0,12   )

        # Exibir gráfico
        plt.show()
        
# def plot_total_chi_squared(chi2_df, test_par, tech_note_best_fit, delta=None):
def plot_total_chi_squared(chi2_df, test_par, delta=None):

    """
    Cria um gráfico do χ² total por parâmetro livre, considerando a soma dos canais.
    """
    
    # Se delta for especificado (ou seja, não for None), usar a coluna "delta_chi"
    if delta is True:
        total_chi2 = chi2_df["delta_chi"].values
        ylabel = r"$\Delta \chi$"
        title = f"$\Delta \chi$ por {test_par} (Todos os Canais)"
    else:
        total_chi2 = chi2_df["chi_total"].values
        ylabel = r"Total $\chi^2$"
        title = f"Total $\chi^2$ por {test_par} (Todos os Canais)"
        
    free_pars = chi2_df["free_pars"].values
    sigma_levels = {"1σ": 1.00, "95% C.L.": 3.84, "3σ": 9.00, "5σ": 25.00}

    interp_func = interp1d(total_chi2, free_pars, kind="linear", fill_value="extrapolate")
    free_par_95 = interp_func(confindence_level)
    
    spl = UnivariateSpline(free_pars, total_chi2, s=0.02)
    x_smooth = np.linspace(free_pars.min(), free_pars.max(), 500)
    y_smooth = spl(x_smooth)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(free_pars, total_chi2, color='blue', linestyle='-', linewidth=2)
    for label, value in sigma_levels.items():
        ax.axhline(y=value, color='gray', linestyle='--', linewidth=1)
        ax.text(x=free_pars[-1], y=value, s=label, fontsize=12, color='gray', ha='left', va='bottom')
    plt.scatter([free_par_95], [confindence_level], color="black", label=f"95% C.L. {test_par} ={free_par_95:.2e}")
    # plt.scatter([tech_note_best_fit], [confindence_level], color="green", label=f" Tech note value {test_par} = {tech_note_best_fit}")
    ax.set_xlabel(f"{test_par}", fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_title(title, fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.ylim(0, 12)
    ax.legend(fontsize=12)
    ax.set_xscale("log")
    plt.show()
    
def compile_and_run():
    """
    Compila e executa o programa C++ responsável por gerar os arquivos de probabilidade
    """
    try:
        result = subprocess.run(["bash", "prob_bsm_loop_make.sh"], capture_output=True, text=True)
        if result.returncode != 0:
            print("Erro durante a compilação ou execução do C++:")
            print(result.stderr)
            return False
        print("Compilação e execução bem-sucedidas.")
        return True
    except FileNotFoundError:
        print("O arquivo prob_bsm_make.sh não foi encontrado.")
        return False


"""
            Script 1 - Utiliza todas as funções criadas.                
"""


"""
Seção Editável - Personalize esta parte com:
1. O caminho onde seus arquivos de probabilidade estão salvos.
2. O nome do parâmetro que está sendo estudado (ex: "$a_{e \mu}$", "$a_{e \tau}$").
3. O valor do parâmetro que vc queira comparar com o seu resultado.
4. O nome do eixo Y para os gráficos (ex: Probabilidade para diferentes canais).
5. Os prefixos usados para nomear os arquivos de dados (como "prob_numu_file", "prob_nue_file", etc.).
6. O nível de confiança (por exemplo, 95%) para determinar o "melhor valor ajustado" (best fit).
"""

# 1, 2, 3.
test_par, directory = r"$\mathring{a}^{(5)}$", "./prob/data/a"
# test_par, directory = r"$\mathring{c}^{(8)}$", "./prob/data/c"

# test_par, directory, tech_note_best_fit = r"$a_{ee}$", "./prob/a/ee", 1.5e-23

# 4. Dicionário para mapear prefixos aos rótulos do eixo Y
prefix_to_ylabel = {
    "prob_numu_file": r"P($\nu_{\mu} \rightarrow \nu_{\mu}$)", 
    "prob_nue_file": r"P($\nu_{\mu} \rightarrow \nu_{e}$)",
    "prob_numubar_file": r"P($\bar{\nu}_{\mu} \rightarrow \bar{\nu}_{\mu}$)",
    "prob_nuebar_file": r"P($\bar{\nu}_{\mu} \rightarrow \bar{\nu}_{e}$)"
}
        
# 5. Prefixos dos arquivos para os diferentes canais
prefixes = ["prob_nue_file", "prob_nuebar_file", "prob_numu_file", "prob_numubar_file"]
# prefixes = ["prob_nue_file", "prob_nuebar_file"]

# 6. melhor valor ajustado
confindence_level = 3.84 

"""
Análise de dados
"""
# Executar o código C++ para gerar os dados
compile_and_run()

chi2_df, all_data_frames = chi_squared(directory, prefixes)
# chi2_df = chi2_df.sort_values(by="free_pars", ascending=True)

#Verifique o DataFrame resultante
if chi2_df is not None:
    print(chi2_df)

# Plotar os gráficos para cada prefixo
for prefix, (data_frames, free_pars_channel) in zip(prefixes, all_data_frames):
    ylabel = prefix_to_ylabel.get(prefix, "Probability")
    data_true = pd.read_csv(f"prob/data/c/{prefix.split('_file')[0]}_STD_bsm.dat", delim_whitespace=True, header=None, names=["E", "P"])
    plot_probability(data_frames, free_pars_channel, ylabel, test_par)
       
# Plotar os gráficos de χ² por parâmetro livre
# plot_chi_squared_by_free_par(chi2_df, prefixes, test_par)

#Plotar os gráficos de χ² total (somando todos os canais)
# plot_total_chi_squared(chi2_df, test_par, tech_note_best_fit, delta=True)