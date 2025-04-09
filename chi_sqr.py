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
    Carrega os dados e retorna APENAS o DataFrame com E e P, 
    e o test_value separadamente (sem coluna redundante)
    """
    # Ler todas as colunas do arquivo
    raw_data = pd.read_csv(filename, delim_whitespace=True, header=None)
    
    # Extrair test_value (3ª coluna)
    test_value = raw_data.iloc[0, 2]  # Índice [linha 0, coluna 2]
    
    # Criar DataFrame apenas com E e P
    data = pd.DataFrame({
        "E": raw_data.iloc[:, 0],
        "P": raw_data.iloc[:, 1]
    })
    
    return data, test_value

def process_files_by_prefix(directory, prefix):
    """
    Retorna um dicionário {test_value: DataFrame} sem coluna redundante
    """
    files = sorted(
        [file for file in os.listdir(directory) if file.startswith(prefix)],
        key=lambda x: int(x.split("_file")[1].split(".")[0])
    )
    
    data_dict = {}
    for file in files:
        file_path = os.path.join(directory, file)
        data, test_value = read_data(file_path)
        data_dict[test_value] = data  # DataFrame já não tem test_value
    
    return data_dict
    
def get_all_data_frames(directory, prefixes):
    """
    Retorna lista de tuplas (data_dict, data_true) para cada prefixo
    """
    all_data = []
    
    for prefix in prefixes:
        true_file = f"prob/data/a/{prefix.split('_file')[0]}_STD_bsm.dat"
        data_true = pd.read_csv(true_file, delim_whitespace=True, header=None, names=["E", "P"])
        
        data_dict = process_files_by_prefix(directory, prefix)
        all_data.append((data_dict, data_true))
    
    return all_data


def plot_probability(data_dict, data_true, ylabel, test_par):
    """
    Cria um gráfico de probabilidade (P por E) para todos os DataFrames fornecidos.
    """
    fig, ax = plt.subplots(figsize=(9, 6))
    
    ax.plot(data_true["E"], data_true["P"], color="black", label="True (SI)", linewidth=2)
    
    for test_value, df in data_dict.items():
        ax.plot(df["E"], df["P"], linestyle='--', 
                label=f"{test_par} = {test_value:.2e}", linewidth=2)
    
    ax.set_xlabel("Energy [GeV]", fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.legend(fontsize=16)
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=16)
    plt.show()

def calculate_chi_squared(all_data, prefixes):
    """
    Calcula o χ² usando dados pré-processados de get_all_data_frames.
    Retorna o mesmo DataFrame mas sem reprocessar arquivos.
    """
    chi_data = {prefix: {} for prefix in prefixes}
    free_pars = []

    # Processar cada conjunto de dados
    for (data_dict, data_true), prefix in zip(all_data, prefixes):
        if not free_pars:  # Coleta os parâmetros uma vez
            free_pars = sorted(data_dict.keys())  # Garante ordem consistente
            
        for test_value, df in data_dict.items():
            observed = df["P"]
            expected = data_true["P"]
            chi_square = np.sum(((observed - expected) ** 2) / expected)
            chi_data[prefix][test_value] = chi_square

    # Criar DataFrame
    df = pd.DataFrame({'free_pars': free_pars})
    
    # Adicionar colunas para cada canal
    for prefix in prefixes:
        df[prefix] = df['free_pars'].map(chi_data[prefix])
    
    # Calcular totais
    df['chi_total'] = df[prefixes].sum(axis=1)
    df['delta_chi'] = df['chi_total'] - df['chi_total'].min()
    
    return df


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

    # interp_func = interp1d(total_chi2, free_pars, kind="linear", fill_value="extrapolate")
    # free_par_95 = interp_func(confindence_level)
    
    # spl = UnivariateSpline(free_pars, total_chi2, s=0.02)
    # x_smooth = np.linspace(free_pars.min(), free_pars.max(), 500)
    # y_smooth = spl(x_smooth)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(free_pars, total_chi2, color='blue', linestyle='-', linewidth=2)
    for label, value in sigma_levels.items():
        ax.axhline(y=value, color='gray', linestyle='--', linewidth=1)
        ax.text(x=free_pars[-1], y=value, s=label, fontsize=12, color='gray', ha='left', va='bottom')
    # plt.scatter([free_par_95], [confindence_level], color="black", label=f"95% C.L. {test_par} ={free_par_95:.2e}")
    ax.set_xlabel(f"{test_par}", fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_title(title, fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.ylim(0, 12)
    ax.legend(fontsize=12)
    ax.set_xscale("log")
    plt.show()
 
def save_chi_squared(df, filename="prob\data\chi_squared_results.csv"):
    """Salva o DataFrame em um arquivo CSV"""
    df.to_csv(filename, index=False)
    print(f"Resultados salvos em {filename}")

def load_chi_squared(filename="prob\data\chi_squared_results.csv"):
    """Carrega o DataFrame de um arquivo CSV"""
    try:
        df = pd.read_csv(filename)
        print(f"Resultados carregados de {filename}")
        return df
    except FileNotFoundError:
        print("Arquivo não encontrado. Execute o cálculo primeiro.")
        return None
 
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
2. O nome do parâmetro que está sendo estudado (ex: r"$\mathring{a}^{(5)}$", "$a_{e \tau}$").
3. O nome do eixo Y para os gráficos (ex: Probabilidade para diferentes canais).
4. Os prefixos usados para nomear os arquivos de dados (como "prob_numu_file", "prob_nue_file", etc.).
5. O nível de confiança (por exemplo, 95%) para determinar o "melhor valor ajustado" (best fit).
"""

# 1, 2, 3.
test_par, directory = r"$\mathring{a}^{(5)}$", "./prob/data/a"
# test_par, directory = r"$\mathring{c}^{(8)}$", "./prob/data/c"

# 3. Dicionário para mapear prefixos aos rótulos do eixo Y
prefix_to_ylabel = {
    "prob_numu_file": r"P($\nu_{\mu} \rightarrow \nu_{\mu}$)", 
    "prob_nue_file": r"P($\nu_{\mu} \rightarrow \nu_{e}$)",
    "prob_numubar_file": r"P($\bar{\nu}_{\mu} \rightarrow \bar{\nu}_{\mu}$)",
    "prob_nuebar_file": r"P($\bar{\nu}_{\mu} \rightarrow \bar{\nu}_{e}$)"
}

# 4. Prefixos dos arquivos para os diferentes canais
prefixes = ["prob_nue_file", "prob_nuebar_file", "prob_numu_file", "prob_numubar_file"]

# 5. melhor valor ajustado
confindence_level = 3.84 

"""
Análise de dados
"""
# Opção para recarregar resultados existentes
RECARREGAR_DADOS = False  # Mude para True para usar dados salvos

if not RECARREGAR_DADOS:
    # Executar o código C++ para gerar os dados
    compile_and_run()
    all_data = get_all_data_frames(directory, prefixes)
    chi2_df = calculate_chi_squared(all_data, prefixes)
    save_chi_squared(chi2_df)  # Salva automaticamente
else:
    chi2_df = load_chi_squared()

#Verifique o DataFrame resultante
if chi2_df is not None:
    print(chi2_df)

# Plotar os gráficos de probabilidade:
for (data_dict, data_true), prefix in zip(all_data, prefixes):
    ylabel = prefix_to_ylabel.get(prefix, "Probability")
    plot_probability(data_dict, data_true, ylabel, test_par)
       
# Plotar os gráficos de χ² por parâmetro livre
plot_chi_squared_by_free_par(chi2_df, prefixes, test_par)

# Plotar os gráficos de χ² total (somando todos os canais) -  Algumas funcinalidades da função que eu deixei comentado: 
# free_par_95: Serve para calcular o ponto exato onde teriamos o valor de a para o nível de  cl desejado.
# x_smooth, y_smooth: Faz uma curva interpolada com os pontos, se necessário extrapolar a curva ou suaviza-la
plot_total_chi_squared(chi2_df, test_par, delta=True)
