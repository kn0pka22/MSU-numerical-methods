import pandas as pd
import matplotlib.pyplot as plt

def read_data(file_address):
    df = pd.read_csv(file_address, delim_whitespace=True, header=None, names=['Iteration', 'Residual', 'Theoretical_Residual'])
    return df

file_name = "output1.txt"
data = read_data(file_name)


plt.figure(figsize=(6, 6))
plt.plot(data['Iteration'], data['Residual'], color='g', label='Process Error', marker='+', markersize=2, alpha=0.7, linewidth=0.8, markeredgecolor='red')
plt.plot(data['Iteration'], data['Theoretical_Residual'], color='b', label='Theoretical Error', marker='+', markersize=2, alpha=0.7, linewidth=0.8, markeredgecolor='purple')

plt.xlabel('Iteration')
plt.ylabel('Error')
plt.title('Convergence Rate')
plt.grid(alpha=0.2, linestyle='--')
plt.legend()
plt.tight_layout()
plt.show()

