import matplotlib.pyplot as plt

def read_error_data(filename):
    """
    Чтение данных из файла.
    Возвращает два списка: значения N и соответствующие ошибки.
    """
    N_values = []
    err_values = []
    with open(filename, 'r') as file:
        for line in file:
            N, err = map(float, line.split())  
            N_values.append(N)
            err_values.append(err)
    return N_values, err_values

def plot_error(N_values, err_values):
    """
    Построение графика ошибки.
    """
    plt.figure(figsize=(10, 6))  
    plt.plot(N_values, err_values, marker='o', linestyle='-', color='b', label='Error')  
    plt.xlabel('N (Number of grid nodes)')  
    plt.ylabel('Error (|True - Approx|)')  
    plt.title('Error vs Number of Grid Nodes')  
    plt.grid(True) 
    plt.legend()  
    plt.show()  

if __name__ == "__main__":
    N_values, err_values = read_error_data("/Users/kn0pka22/Desktop/msu/MSU-numerical-methods/part_2/error_data.txt")
    
    plot_error(N_values, err_values)