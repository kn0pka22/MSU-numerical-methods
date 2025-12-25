import matplotlib.pyplot as plt
import numpy as np

def read_data(filename):
    data = np.loadtxt(filename)
    N_values = data[:, 0]
    errors = data[:, 1:]
    return N_values, errors

def plot_convergence(N_values, errors):
    plt.figure(figsize=(12, 6))
    
    # График 1: Логарифмический масштаб для определения порядка сходимости
    plt.subplot(1, 2, 1)
    plt.loglog(N_values, errors, 'o-', linewidth=2)
    plt.xlabel('Number of points (log scale)')
    plt.ylabel('Error (log scale)')
    plt.title('Log-Log Plot for Convergence Order')
    plt.grid(True, which="both", ls="-")
    
    # Добавление теоретических линий сходимости для сравнения
    for i, color in enumerate(plt.rcParams['axes.prop_cycle'].by_key()['color'][:4]):
        order = 1 if i < 2 else 2  # Порядок для схем 1-2 и 3-4
        x_ref = N_values[-1] / 2
        y_ref = errors[-1, i] * 2
        plt.loglog([x_ref, x_ref*2], [y_ref, y_ref/(2**order)], '--', color=color, 
                  label=f'Theoretical order {order}')
    
    # График 2: Нормированные ошибки
    plt.subplot(1, 2, 2)
    for i in range(4):
        normalized_errors = errors[:, i] / errors[0, i]
        plt.plot(N_values, normalized_errors, 'o-', linewidth=2, 
                label=f'Scheme {i+1}')
    
    plt.xlabel('Number of points')
    plt.ylabel('Normalized error')
    plt.title('Normalized Error Comparison')
    plt.xscale('log')
    plt.grid(True, which="both", ls="-")
    
    plt.legend()
    plt.tight_layout()
    plt.savefig('improved_convergence_plot.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    N_values, errors = read_data('out.txt')
    plot_convergence(N_values, errors)