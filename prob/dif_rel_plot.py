import numpy as np
import matplotlib.pyplot as plt

# Função para carregar dados do arquivo
def load_data(filename):
    data = np.loadtxt(filename)
    x = data[:, 0]  # Primeira coluna: E [GeV]
    y = data[:, 1]  # Segunda coluna: P(nu_mu -> nu_e)
    return x, y

# Função para plotar gráficos com dois eixos
def plot_graph(x_std, y_std, x_bsm, y_bsm, title, ylabel, ylim_max):
    # Calcular a diferença relativa
    y_diff = (y_bsm - y_std) / y_std

    # Criar a figura e os eixos
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Configurar o gráfico principal
    ax1.grid(True)
    ax1.plot(x_std, y_std, label='Oscilação padrão', color='blue', linewidth=2)
    ax1.plot(x_bsm, y_bsm, label='Oscilação padrão no arquivo bsm', color='black', linestyle='--', linewidth=2)
    ax1.set_xlabel("E [GeV]", fontsize=16)
    ax1.set_ylabel(ylabel, fontsize=16)
    ax1.legend(loc="upper left")
    plt.ylim(0, ylim_max)
    # plt.xlim(0.5, 5)



    # Eixo secundário para a diferença relativa
    ax2 = ax1.twinx()  # Adiciona um segundo eixo Y
    ax2.plot(x_std, y_diff, label='Diferença Relativa', color='red', linestyle=':', linewidth=2)
    ax2.set_ylabel("Diferença Relativa", fontsize=16)
    ax2.legend(loc="upper right")

    # Configurar título do gráfico
    plt.title(title, fontsize=18)
    plt.show()

# Carregar e plotar cada gráfico

# Gráfico 1: prob_nue_STD.dat e prob_nue_STD_bsm.dat
x_std, y_std = load_data("prob_nue_STD.dat")
x_bsm, y_bsm = load_data("prob_nue_STD_bsm.dat")
plot_graph(x_std, y_std, x_bsm, y_bsm, "Grafico 1: prob_nue_STD vs prob_nue_STD_bsm", r"P($\nu_{\mu} \rightarrow \nu_{e}$)", 0.25)

# Gráfico 2: prob_nuebar_STD.dat e prob_numu_STD.dat
x_std, y_std = load_data("prob_nuebar_STD.dat")
x_bsm, y_bsm = load_data("prob_nuebar_STD_bsm.dat")
plot_graph(x_std, y_std, x_bsm, y_bsm, "Grafico 2: prob_nuebar_STD vs prob_nuebar_STD_bsm", r"P($\overline{\nu}_{\mu} \rightarrow \overline{\nu}_{e}$)", 0.1)

# Gráfico 3: prob_numu_STD.dat e prob_numu_STD_bsm.dat
x_std, y_std = load_data("prob_numu_STD.dat")
x_bsm, y_bsm = load_data("prob_numu_STD_bsm.dat")
plot_graph(x_std, y_std, x_bsm, y_bsm, "Grafico 3: prob_numu_STD vs prob_numu_STD_bsm", r"P($\nu_{\mu} \rightarrow \nu_{\mu}$)", 1.2)

# Gráfico 4: prob_numubar_STD.dat e prob_numu_STD_bsm.dat
x_std, y_std = load_data("prob_numubar_STD.dat")
x_bsm, y_bsm = load_data("prob_numubar_STD_bsm.dat")
plot_graph(x_std, y_std, x_bsm, y_bsm, "Grafico 4: prob_numubar_STD vs prob_numubar_STD_bsm", r"P($\overline{\nu}_{\mu} \rightarrow \overline{\nu}_{\mu}$)", 1.2)
