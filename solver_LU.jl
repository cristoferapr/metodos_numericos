using LinearAlgebra

# Función para la factorización LU
function factorizacion_lu(matriz)
    # Obtener el tamaño de la matriz (n x n)
    n = size(matriz, 1)

    # Inicializar una matriz identidad L de tipo Float64
    L = Matrix{Float64}(I, n, n)

    # Crear una copia de la matriz original y almacenarla en U
    U = copy(matriz)

    # Variable para realizar un seguimiento de los pasos
    paso = 1

    # Iterar a través de las filas
    for k in 1:n-1
        # Iterar a través de las filas subsiguientes
        for i in k+1:n
            # Calcular el factor de multiplicación
            factor = U[i, k] / U[k, k]

            # Almacenar el factor en la matriz L
            L[i, k] = factor

            # Actualizar la matriz U
            for j in k:n
                U[i, j] -= factor * U[k, j]
            end

            # Imprimir los pasos y las matrices L y U
            println("Paso $paso: Fila $i - $factor Fila $k")
            paso += 1
            println("Matriz L:")
            for row in 1:n
                for col in 1:n
                    print(L[row, col], " ")
                end
                println()
            end
            println("Matriz U:")
            for row in 1:n
                for col in 1:n
                    print(U[row, col], " ")
                end
                println()
            end
            println()
        end
    end

    # Devolver las matrices L y U
    return L, U
end

# Función para resolver un sistema de ecuaciones lineales utilizando las matrices LU y el vector de solucion B
function solver_sistemaEcu(L, U, B)
    # Se calcula a, b, c en funcion del vector de solucion B y la matriz L
    a = B[1]
    b = (B[2] - L[2, 1] * a) / L[2, 2]
    c = (B[3] - L[3, 1] * a - L[3, 2] * b) / L[3, 3]

    # Imprimir los resultados intermedios
    println("a = $a, b = $b, c = $c")

    # Calcular x,y,z en función de la matriz U y los datos a, b y c obtenidos anteriormente
    z = c / U[3, 3]
    y = (b - z * U[2, 3]) / U[2, 2]
    x = (a - y * U[1, 2] - z * U[1, 3]) / U[1, 1]

    # Imprimir la solución final
    println("\nSolución:")
    println("x = $x")
    println("y = $y")
    println("z = $z")
end

# Definir el sistema de ecuaciones como una matriz A 3x3 y un vector B
A = [3.0 2.0 1.0;
     2.0 8.0 7.0;
     1.0 7.0 14.0]
B = [10, 11, 15]

# Realizar la factorización LU y resolver el sistema de ecuaciones
L, U = factorizacion_lu(A)
solver_sistemaEcu(L, U, B)
