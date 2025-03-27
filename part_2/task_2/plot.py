import numpy as np
import matplotlib.pyplot as plt

def main():
    try:
        data = np.loadtxt('data.txt')
        log_N = data[:, 0]
        log_y1 = data[:, 1]
        log_y2 = data[:, 2]
        log_y3 = data[:, 3]
        log_y4 = data[:, 4]
        
        plt.figure(figsize=(10, 8))
        plt.plot(log_N, log_y1, 'o-', label='Scheme 1', linewidth=2, markersize=8)
        plt.plot(log_N, log_y2, 's-', label='Scheme 2', linewidth=2, markersize=8)
        plt.plot(log_N, log_y3, 'd-', label='Scheme 3', linewidth=2, markersize=8)
        plt.plot(log_N, log_y4, '^-', label='Scheme 4', linewidth=2, markersize=8)
        
        plt.xlabel('log(N)', fontsize=12)
        plt.ylabel('log(error)', fontsize=12)
        plt.title('Numerical Schemes Comparison', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend(fontsize=12, bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        plt.savefig('schemes_comparison.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()