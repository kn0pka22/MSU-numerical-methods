import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def plot_solutions():
    try:
        x = np.loadtxt('x_grid.txt')
        t = np.loadtxt('t_grid.txt')
        u_num = np.loadtxt('solution_matrix.txt')
        
        X, T = np.meshgrid(x, t)
        u_exact = np.exp(-np.pi**2 * T) * X * (X - 1)
        
        fig = plt.figure(figsize=(18, 6))
        
        # Numerical solution
        ax1 = fig.add_subplot(131, projection='3d')
        ax1.plot_surface(X, T, u_num, cmap='viridis')
        ax1.set_title('Numerical Solution')
        
        # Exact solution
        ax2 = fig.add_subplot(132, projection='3d')
        ax2.plot_surface(X, T, u_exact, cmap='plasma')
        ax2.set_title('Exact Solution')
        
        # Error
        ax3 = fig.add_subplot(133, projection='3d')
        error = u_num - u_exact
        ax3.plot_surface(X, T, error, cmap='coolwarm')
        ax3.set_title('Error')
        
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print(f"Error: {e}")

plot_solutions()

# !make && ./proga convergence implicit 2


# implicit_data = pd.read_csv("convergence_results.txt", sep='\s+')

# implicit_data['h_ratio'] = implicit_data['h'].shift(1) / implicit_data['h'] 

# implicit_data['error_ratio'] = implicit_data['error'].shift(1) / implicit_data['error'] 
# implicit_data['log(error_ratio)'] = (np.log(implicit_data['error'].shift(1) / implicit_data['error'])) / \
#                                     (np.log(implicit_data['M'] / implicit_data['M'].shift(1))) 



# !make && ./proga convergence implicit 1


# implicit_data = pd.read_csv("convergence_results.txt", sep='\s+')

# implicit_data['h_ratio'] = implicit_data['h'].shift(1) / implicit_data['h'] 

# implicit_data['error_ratio'] = implicit_data['error'].shift(1) / implicit_data['error'] 
# implicit_data['log(error_ratio)'] = (np.log(implicit_data['error'].shift(1) / implicit_data['error'])) / \
#                                     (np.log(implicit_data['M'] / implicit_data['M'].shift(1))) 


# ! make && ./proga convergence explicit 3


# explicit_data = pd.read_csv("convergence_results.txt", sep='\s+')

# explicit_data['h_ratio'] = explicit_data['h'].shift(1) / explicit_data['h'] 

# explicit_data['error_ratio'] = explicit_data['error'].shift(1) / explicit_data['error'] 
# explicit_data['log(error_ratio)'] = (np.log(explicit_data['error'].shift(1) / explicit_data['error'])) / \
#                                     (np.log(explicit_data['M'] / explicit_data['M'].shift(1))) 
#                                     # ((explicit_data['tau'].shift(1) + explicit_data['h'].shift(1) * explicit_data['h'].shift(1))/ (explicit_data['tau'] + explicit_data['h'] * explicit_data['h']) )



plt.figure(figsize=(12, 5))

# Сходимость по x
plt.subplot(121)
plt.plot(implicit_data['log_h'], implicit_data['log_error'], 'bo-', label='Неявная схема')
plt.plot(explicit_data['log_h'], explicit_data['log_error'], 'ro-', label='Явная схема')
plt.plot(implicit_data['log_h'], 2*implicit_data['log_h'] - 8, 'k--', label='Теоретическая O(h²)')
plt.xlabel('log(h)')
plt.ylabel('log(error)')
plt.title('Пространственная сходимость')
plt.legend()
plt.grid(True)

# Сходимость по t
plt.subplot(122)
plt.plot(np.log(implicit_data['tau']), np.log(implicit_data['error']), 'b^-', label='Неявная схема')
plt.plot(np.log(explicit_data['tau']), np.log(explicit_data['error']), 'r^-', label='Явная схема')
plt.plot(np.log(implicit_data['tau']), np.log(implicit_data['tau']) - 5, 'k--', label='Теоретическая O(τ)')
plt.plot(np.log(explicit_data['tau']), 2*np.log(explicit_data['tau']) - 9, 'k:', label='Теоретическая O(τ²)')
plt.xlabel('log(τ)')
plt.ylabel('log(error)')
plt.title('Временная сходимость')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()