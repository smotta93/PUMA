import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Carregar os dados do primeiro arquivo
file_path1 = './prob/a/mutau/prob_nue_file2.dat'
data1 = pd.read_csv(file_path1, delim_whitespace=True, header=None, names=["Energy", "Probability", "test_value"])

# Carregar os dados do segundo arquivo
file_path2 = './prob/a/prob_nue_STD_bsm.dat'
data2 = pd.read_csv(file_path2, delim_whitespace=True, header=None, names=["Energy", "Probability", "test_value"])

# Extrair energia e probabilidade do primeiro arquivo
x1 = data1['Energy']
y1 = data1['Probability']

# Extrair energia e probabilidade do segundo arquivo
x2 = data2['Energy']
y2 = data2['Probability']

# Definir a janela para a média móvel
window_size = 20  # Tamanho da janela para a média móvel (ajuste conforme necessário)
window = np.ones(window_size) / window_size  # Filtro de média simples

# Aplicar a média móvel ao primeiro conjunto de dados
y1_smoothed = np.convolve(y1, window, mode='same')

# Plotando
plt.figure(figsize=(10, 6))

# Dados e curva suavizada do primeiro arquivo
plt.scatter(x1, y1, color='blue', label=r'Dados Originais considerando $a_5 = 1 \times 10^{-18}$', alpha=0.7, s=10)
plt.plot(x1, y1_smoothed, 'g--', label=r'Curva Suavizada para $a_5$', linewidth=2)

# Dados do segundo arquivo (sem suavização)
plt.plot(x2, y2, color='black', label='Oscilação Padrão', linewidth=2)

# Configurações do gráfico
plt.xlabel('Energy (GeV)')
plt.ylabel(r"P($\nu_{\mu} \rightarrow \nu_{e}$)")
plt.title('Comparação entre Dados Originais e Curva Suavizada')
plt.grid(True)
plt.ylim(0, 0.2)
plt.xlim(0.1,0.8)
plt.legend()
plt.show()
