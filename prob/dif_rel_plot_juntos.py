import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Função para carregar dados do arquivo
def load_data(filename):
    data = pd.read_csv(filename, delim_whitespace=True, header=None, names=["E", "P"])
    return data["E"].values, data["P"].values

# Função para plotar gráficos com dois eixos
def plot_graph(ax1, x_std, y_std, x_bsm, y_bsm, ylabel, ylim_max, show_legend=False):
    # Calcular a diferença relativa
    y_diff = np.sqrt((y_bsm - y_std) ** 2) / y_std

    # Configurar o gráfico principal
    ax1.grid(True)
    line1, = ax1.plot(x_std, y_std, label='Oscilação padrão', color='blue', linewidth=2)
    line2, = ax1.plot(x_bsm, y_bsm, label='BSM, todos os par. zero (OP)', color='black', linestyle='--', linewidth=2)
    ax1.tick_params(axis='both', which='major', labelsize=12)  # Aumenta a fonte dos números do eixo x e y
    ax1.set_xlabel("E [GeV]", fontsize=14)
    ax1.set_ylabel(ylabel, fontsize=14)
    ax1.set_ylim(0, ylim_max)
    ax1.set_xlim(0.5, 5)

    # Eixo secundário para a diferença relativa
    ax2 = ax1.twinx()  # Adiciona um segundo eixo Y
    line3, = ax2.plot(x_std, y_diff, label='Diferença Relativa', color='red', linestyle=':', linewidth=2)
    ax2.set_ylabel("Diferença Relativa", fontsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=12)  # Aumenta a fonte dos números do eixo secundário também

    
    # Exibir a legenda apenas no primeiro gráfico
    if show_legend:
        lines = [line1, line2, line3]
        labels = [line.get_label() for line in lines]
        ax1.legend(lines, labels, loc="upper right", fontsize=10)

# Criar a figura e os eixos para 4 gráficos lado a lado
fig, axs = plt.subplots(2, 2, figsize=(14, 10))
fig.tight_layout(pad=10.0)  # Ajusta o espaçamento entre os gráficos

# Carregar e plotar cada gráfico

# Gráfico 1: prob_nue_STD.dat e prob_nue_STD_bsm.dat (com legenda)
x_std, y_std = load_data("prob_nue_STD.dat")
x_bsm, y_bsm = load_data("prob_nue_STD_bsm.dat")
plot_graph(axs[0, 0], x_std, y_std, x_bsm, y_bsm, r"P($\nu_{\mu} \rightarrow \nu_{e}$)", 0.15, show_legend=True)

# Gráfico 2: prob_nuebar_STD.dat e prob_nuebar_STD_bsm.dat (sem legenda)
x_std, y_std = load_data("prob_nuebar_STD.dat")
x_bsm, y_bsm = load_data("prob_nuebar_STD_bsm.dat")
plot_graph(axs[0, 1], x_std, y_std, x_bsm, y_bsm, r"P($\overline{\nu}_{\mu} \rightarrow \overline{\nu}_{e}$)", 0.03)

# Gráfico 3: prob_numu_STD.dat e prob_numu_STD_bsm.dat (sem legenda)
x_std, y_std = load_data("prob_numu_STD.dat")
x_bsm, y_bsm = load_data("prob_numu_STD_bsm.dat")
plot_graph(axs[1, 0], x_std, y_std, x_bsm, y_bsm, r"P($\nu_{\mu} \rightarrow \nu_{\mu}$)", 1.2)

# Gráfico 4: prob_numubar_STD.dat e prob_numubar_STD_bsm.dat (sem legenda)
x_std, y_std = load_data("prob_numubar_STD.dat")
x_bsm, y_bsm = load_data("prob_numubar_STD_bsm.dat")
plot_graph(axs[1, 1], x_std, y_std, x_bsm, y_bsm, r"P($\overline{\nu}_{\mu} \rightarrow \overline{\nu}_{\mu}$)", 1.2)

# Salvar e exibir o gráfico completo
plt.savefig('dif_rel.png')
plt.show()
