import matplotlib.pyplot as plt
import numpy as np

def read_data(filename):
    data = np.loadtxt(filename)
    log_N = data[:, 0]
    log_errors = data[:, 1:]
    return log_N, log_errors

def plot_results(log_N, log_errors):
    plt.figure(figsize=(10, 10))
    plt.xlabel('log(number of points)')
    plt.ylabel('log(error)')
    plt.title('Convergence of Different Schemes')
    
    schemes = ['Scheme 1', 'Scheme 2', 'Scheme 3', 'Scheme 4']
    markers = ['o-', 's-', '^-', 'd-']
    
    for i in range(4):
        plt.plot(log_N, log_errors[:, i], markers[i], label=schemes[i])
    
    plt.legend()
    plt.grid(True)
    # plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

    # plt.minorticks_on()
    # plt.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.5)

    # plt.xlim(0, 1)
    # plt.ylim(0, 1)

    plt.savefig('convergence_plot.png')
    plt.show()

if __name__ == "__main__":
    log_N, log_errors = read_data('out.txt')
    plot_results(log_N, log_errors)