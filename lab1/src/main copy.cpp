#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <stack>
#include <iomanip>
#include <unordered_map>
#include </usr/include/eigen3/Eigen/Dense>

template<typename T>
class OrderIndexMap {
public:
    OrderIndexMap(size_t expected_size) {
        map.reserve(expected_size); // заранее зарезервировать размер
    }

    // Возвращает индекс элемента: либо уже существующий, либо новый
    int get_or_add(const T& value) {
        auto it = map.find(value);
        if (it != map.end()) {
            return it->second;
        } else {
            int index = current_index++;
            map[value] = index;
            return index;
        }
    }

    void clear() {
        map.clear();       // Удаляем все элементы из хэш-таблицы
        current_index = 0; // Сбрасываем счётчик индексов
    }

private:
    std::unordered_map<T, int> map;
    int current_index = 0;
};

// Граничные условия:
// type 0 - начальные условия, T(t=0) = val
// type 1 - Dirichlet,         T      = val
// type 2 - Neumann,           dT/dn  = val
// type 3 - Robin,             dT/dn  = val * T
class BoundaryCondition {
    public:
        int type;     
        double val;
        int n_i, n_j;
        BoundaryCondition(int _type, double _val, int _n_i, int _n_j) : type(_type), val(_val), n_i(_n_i), n_j(_n_j) {}
    };
    
// Класс узла сетки
class Node {
    public:
        BoundaryCondition bc; 
        int i, j;       // Индексы в КЭ сетке
        double x, y;    // Координаты узла
        double t;       // Температура узла
        Node(int _i, int _j, double _x, double _y, int _bc_type, double _bc_val, int _t = 0., int _n_i = 0, int _n_j = 0) : i(_i), j(_j), x(_x), y(_y), bc(_bc_type, _bc_val, _n_i, _n_j), t(_t) {}
    };

class PlateMeshGrid {
    public:
    std::vector<Node> nodes;
    int nRows, nCols;   // nRows - число строк (индекс i), nCols - число столбцов (индекс j)
    double dx, dy;
       
    // Чтение сетки из CSV файла
    int read_mesh_from_csv(const std::string& filename) {
        std::ifstream fin(filename);
        std::vector<std::string> lines;
        std::string line;  
        
        // прочтение метаданных о сетке
        while (std::getline(fin, line)) {
            const auto first = line.find_first_not_of(" \t");
            if (first == std::string::npos || line[first] == '#') 
                continue;
            std::istringstream ss(line);
            char delim;
            ss >> this->nCols >> delim
               >> this->nRows >> delim
               >> this->dx >> delim
               >> this->dy >> delim;
            break; //метаданные о сетке прочитаны
        }

        while (std::getline(fin, line)) {
            const auto first = line.find_first_not_of(" \t");
            // если строка пустая или первый непробельный символ — '#', пропускаем её
            if (first == std::string::npos || line[first] == '#') 
                continue;
            std::istringstream ss(line);
            int i, j, bc_type;
            double x, y, bc_val;
            int n_i, n_j;
            char delim;
            ss >> i >> delim
               >> j >> delim
               >> x >> delim
               >> y >> delim
               >> bc_type >> delim
               >> bc_val;
            if (bc_type == 2 || bc_type == 3) {
               ss >> delim >> n_i >> delim
                  >> n_j;
               nodes.emplace_back(i, j, x, y, bc_type, bc_val, 0., n_i, n_j);
            }
            else{
               nodes.emplace_back(i, j, x, y, bc_type, bc_val, bc_val);
            }

        }
        fin.close();

        return 0;
    }
 
    // Функция отрисовки (записи) распределения температуры
    int print_mesh(const std::string& filename) {
        std::ofstream fout(filename);
        fout << std::setprecision(6);
        
        for (const auto& node : nodes) {
            // Формат: X  Y  Temperature  dX/2  dY/2
            fout << node.x << " " << node.y << " " << node.t 
                 << " " << dx/2 << " " << dy/2 << "\n";
        }
        fout.close();
        return 0;
    }
};

class PlateTemperatureSolver {
    public:
        int solve_mesh(PlateMeshGrid& mesh, double time_final) {
            const double alpha = 1.0;                    // Коэффициент теплопроводности
            const int N = mesh.nodes.size();             // Число уравнений в СЛАУ будет соответсвтовать числу узлов, в которых не граничное условие 1 рода
            const double dt = 1;                         // Шаг по времени
            int timesteps = (int) time_final / dt;       // Инициализация начального момента времени
            //int timesteps = 1;
            const int nodes_count = mesh.nodes.size();

            // A B = C, B = A^-1 C
            std::vector<std::vector<double>> A(N, std::vector<double>(N));
            std::vector<double> C(N);
            OrderIndexMap<int> B(N); 
            for (auto& node : mesh.nodes) {
                B.get_or_add(node.j * mesh.nCols + node.i);
            }

            for (int k = 0; k < timesteps; k++) {
                for (auto& row : A) std::fill(row.begin(), row.end(), 0.0);
                std::fill(C.begin(), C.end(), 0.0);

                for (int n = 0; n < nodes_count; n++) { // составляем СЛАУ
                    Node node = mesh.nodes[n];
                    int b_poz_i_j, b_poz_im_j, b_poz_ip_j, b_poz_i_jm, b_poz_i_jp, b_poz, b_poz_n;
                    switch (node.bc.type) { // см описание класса
                        case 0:
                            b_poz_i_j = B.get_or_add(node.j * mesh.nCols + node.i);        // T_{i,j}
                            b_poz_im_j = B.get_or_add(node.j * mesh.nCols + node.i - 1);   // T_{i-1,j}
                            b_poz_ip_j = B.get_or_add(node.j * mesh.nCols + node.i + 1);   // T_{i+1,j}
                            b_poz_i_jm = B.get_or_add((node.j - 1) * mesh.nCols + node.i); // T_{i,j-1}
                            b_poz_i_jp = B.get_or_add((node.j + 1) * mesh.nCols + node.i); // T_{i,j+1}

                            A[n][b_poz_i_j] = (mesh.dx*mesh.dx*mesh.dy*mesh.dy) \
                                              + 2*alpha*mesh.dy*mesh.dy*dt \
                                              + 2*alpha*mesh.dx*mesh.dx*dt;

                            A[n][b_poz_im_j] = -alpha*mesh.dy*mesh.dy*dt;
                            A[n][b_poz_ip_j] = -alpha*mesh.dy*mesh.dy*dt;

                            A[n][b_poz_i_jm] = -alpha*mesh.dx*mesh.dx*dt;
                            A[n][b_poz_i_jp] = -alpha*mesh.dx*mesh.dx*dt;                           

                            C[n] = mesh.dx*mesh.dx*mesh.dy*mesh.dy * mesh.nodes[b_poz_i_j].t; 
                            break;
                        case 1:
                            b_poz = B.get_or_add(node.j * mesh.nCols + node.i);
                            A[n][b_poz] = 1.;
                            C[n] = node.bc.val;
                            break;
                        case 2:
                            b_poz = B.get_or_add(node.j * mesh.nCols + node.i);
                            b_poz_n = B.get_or_add(node.bc.n_j * mesh.nCols + node.bc.n_i);
                            A[n][b_poz] = -1.;
                            A[n][b_poz_n] = 1.;
                            C[n] = node.bc.val * sqrt((mesh.dx*mesh.dx*(node.bc.n_i-node.i)*(node.bc.n_i-node.i) \
                                        + mesh.dy*mesh.dy*(node.bc.n_j-node.j)*(node.bc.n_j-node.j)));
                            break;
                        case 3:
                            b_poz = B.get_or_add(node.j * mesh.nCols + node.i);
                            b_poz_n = B.get_or_add(node.bc.n_j * mesh.nCols + node.bc.n_i);
                            A[n][b_poz] = -(node.bc.val * sqrt((mesh.dx*mesh.dx*(node.bc.n_i-node.i)*(node.bc.n_i-node.i) \
                                                                + mesh.dy*mesh.dy*(node.bc.n_j-node.j)*(node.bc.n_j-node.j))) + 1);
                            A[n][b_poz_n] = 1.;
                            C[n] = 0.;
                            break;
                    }
                }

                Eigen::MatrixXd eigenA(N, N);
                Eigen::VectorXd eigenC(N);
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        eigenA(i, j) = A[i][j];
                    }
                    eigenC(i) = C[i];
                }

                Eigen::VectorXd solution = eigenA.partialPivLu().solve(eigenC);

                for (auto& node : mesh.nodes) {
                    if (node.bc.type != 1) {
                        node.t = solution(B.get_or_add(node.j * mesh.nCols + node.i));
                    }
                }
            }
            return 0;
        }
};

// Определяет параметры моделирования, задаваемые пользователем
struct ModelParameters {
    double duration;         // Продолжительность моделирования
    bool service_mode;       // Флаг служебного режима без файлового вывода (для benchmark)
    std::string output_file; // Название файла для вывода результатов моделирования
    std::string mesh_file;   // Название файла-источника сетки для моделирования
};

class ParameterParser {
    public:
    ParameterParser(int c, char** v) : _argument_count(c), _argument_strings(v) {}
  
    // Разбирает аргументы, переданные при создании, для определения параметров моделирования.
    // выбрасывает std::invalid_argument при некорректных значениях параметров.
    ModelParameters getParameters() const {
      // Значения по умолчанию.
      double duration = 25.;
      std::string output_file = "data.txt";
      bool service_mode = false;
      std::string mesh_file;

      // Разбор аргументов при помощи getopt.
      int option = 0;
      while ((option = getopt(_argument_count, _argument_strings, "i:o:t:s")) != -1) {
        switch (option) {
        case 't':
          duration = atof(optarg);
          if (duration <= 0.) {
            throw std::invalid_argument("ModelParameters: duration must be positive");
          }
          break;
        case 'o':
          output_file = optarg;
          break;
        case 'i':
          mesh_file = optarg;
        case 's':
          service_mode = true;
          break;
        default:
          break;
        }
      }

      if (mesh_file.empty()) {
        std::cout << "Parameter -i (mesh_file) is required" << std::endl;
        throw std::invalid_argument("Parameter -i (mesh_file) is required");
      }

      ModelParameters parameters;
      parameters.output_file = output_file;
      parameters.mesh_file = mesh_file;
      parameters.duration  = duration;
      parameters.service_mode = service_mode;

      return parameters;
    }
  
    private:
      const int _argument_count; // Число аргументов запуска программы.
      char** _argument_strings;  // Массив строковых значений аргументов запуска программы.
};

int main(int argc, char** argv) {
    try {
        ParameterParser parser(argc, argv);
        ModelParameters params = parser.getParameters();
        PlateMeshGrid mesh;
        mesh.read_mesh_from_csv(params.mesh_file);

        PlateTemperatureSolver solver;
        solver.solve_mesh(mesh, 25);
        mesh.print_mesh(params.output_file);
    } catch (...) {
        return 1;
    }
    std::cout << "Simulation completed successfully.\n";
    return 0;
}




/*
int gauss_elimination(Matrix* mat) {
    gauss_forward_elimination(mat);
    int res = gauss_check_solutions(mat);
    if (res != 0)
        gauss_backward_elimination(mat);
    return res;
}

void gauss_forward_elimination(Matrix* mat) {
    for (int k = 0; k < mat->rows; k++) {
        /* поиск строки с максимальным абсолютным значением в текущем столбце 
           (от строки k вниз), чтобы избежать деления на слишком маленькие числа 
           и улучшить численную стабильность алгоритма. */
/*           int max_row = k;
           for (int i = k + 1; i < mat->rows; i++) {
               if (fabs(mat->matrix[i][k]) > fabs(mat->matrix[max_row][k])) {
                   max_row = i;
               }
           }
   
           // Меняем строки, если необходимо
           if (max_row != k) {
               swap_rows(mat, k, max_row);
           }
           // Нормализация основной строки
           double pivot = mat->matrix[k][k];
           for (int j = k; j < mat->cols; j++) {
               mat->matrix[k][j] /= pivot;
           }
           
           // Обнуление элементов под основной строкой
           for (int i = k + 1; i < mat->rows; i++) {
               double factor = mat->matrix[i][k];
               mat->matrix[i][k] = 0.;
               for (int j = k+1; j < mat->cols; j++) {
                   mat->matrix[i][j] -= factor * mat->matrix[k][j];
               }
           }
       }
   }
   
   void gauss_backward_elimination(Matrix* mat) {
       for (int k = mat->rows - 1; k >= 0; k--) {
           for (int i = k - 1; i >= 0; i--) {
               double factor = mat->matrix[i][k];
               mat->matrix[i][k] = 0.;
               for (int j = k + 1; j < mat->cols; j++) {
                   mat->matrix[i][j] -= factor * mat->matrix[k][j];
               }
           }
       }
   }
   
   int gauss_check_solutions(Matrix* mat) {
       int rank = 0;
    
       for (int i = 0; i < mat->rows; i++) {
           bool row_nonzero = false;
           for (int j = 0; j < mat->cols - mat->is_augmented; j++) {
               if (fabs(mat->matrix[i][j]) > 1e-9) {  // Пороговое значение для сравнения с нулем
                   row_nonzero = true;
                   break;
               }
           }
           // 0 0 0 0 | 4
           if (mat->is_augmented && !row_nonzero && fabs(mat->matrix[i][mat->cols]) > 1e-9)
               return 0;
           if (row_nonzero) rank++;
       }
   
       if (rank == mat->rows)
           return rank;
       else
           return -rank;
   }
*/

// -------------------------------

/*
Matrix create_matrix(int rows, int cols, bool is_augmented, double (*reference)[cols]) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.is_augmented = is_augmented;
    mat.matrix = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        mat.matrix[i] = (double*)malloc(cols * sizeof(double));
    }
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            mat.matrix[i][j] = reference[i][j];
        }
    }
    return mat;
}

Matrix copy_matrix_with_cut_cols(Matrix* mat, int* sorted_delete_cols_indexes, int del_cols_count) {
    Matrix new_mat;
    new_mat.rows = mat->rows;
    new_mat.cols = mat->cols - del_cols_count;
    new_mat.is_augmented = mat->is_augmented;
    new_mat.matrix = (double**)malloc(mat->rows * sizeof(double*));
    for (int i = 0; i < mat->rows; i++) {
        new_mat.matrix[i] = (double*)malloc(new_mat.cols * sizeof(double));
    }
    for (int i = 0; i < mat->rows; i++) {
        int k = 0;
        for (int j = 0; j < mat->cols; j++) {
            if (k < del_cols_count && j == sorted_delete_cols_indexes[k]) {
                k++;
            }
            else {
                new_mat.matrix[i][j-k] = mat->matrix[i][j];
            }
        }
    }
    return new_mat;
}

void matrix_free(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++) {
        free(mat->matrix[i]);
    }
    free(mat->matrix);
}

void matrix_print(Matrix matrix_to_print) {
    for (int i = 0; i < matrix_to_print.rows; i++) {
        for (int j = 0; j < matrix_to_print.cols - matrix_to_print.is_augmented; j++) {
            printf("%5.2f ", matrix_to_print.matrix[i][j]);
        }
        if (matrix_to_print.is_augmented) {
            printf(" | %5.2f", matrix_to_print.matrix[i][matrix_to_print.cols-1]);
        }       
        putchar('\n');
    }
    putchar('\n');
}

void swap_rows(Matrix* mat, int row1, int row2) {
    for (int i = 0; i < mat->cols; i++) {
        double temp = mat->matrix[row1][i];
        mat->matrix[row1][i] = mat->matrix[row2][i];
        mat->matrix[row2][i] = temp;
    }
}
*/



/*
// 0 - СЛАУ несовместна
// -n - бесконечное количество решений (n - установленный ранг матрицы)
// n - СЛАУ имеет единственное решение (n - установленный ранг матрицы)
int gauss_elimination(Matrix* mat);

// Прямой ход для приведения к верхнетреугольному виду
void gauss_forward_elimination(Matrix* mat);

// Обратный ход для приведения к диагональному виду
void gauss_backward_elimination(Matrix* mat);

// Функция для проверки количества решений СЛАУ
// Принимает на вход матрицу верхнетреугольного вида
// 0 - СЛАУ несовместна (!!определяется только если переданная матрица была расширенной!!)
// -m - бесконечное количество решений (m - установленный ранг матрицы)
// n - СЛАУ имеет единственное решение (n - установленный ранг матрицы)
int gauss_check_solutions(Matrix* mat);
*/

/*
// структура матрицы 
typedef struct {
    int rows;             // количество строк матрицы
    int cols;             // количество колонок матрицы
    bool is_augmented;    // является ли матрица расширенной
    double** matrix;            
} Matrix;

// Функция для создания и инициализации матрицы
// !! требует освобождения памяти !!
Matrix create_matrix(int rows, int cols, bool is_augmented, double (*reference)[cols]);

// функция копирования матрицы с вырезанием указанных столбцов
// массив удаляемых столбцов должен быть отсортирован 
// !! требует освобождения памяти !!
Matrix copy_matrix_with_cut_cols(Matrix* mat, int* sorted_delete_cols_indexes, int delete_cols_count);

// Функция для освобождения динамической памяти, выделенной под матрицу
void matrix_free(Matrix* mat);

// функция форматированной печати матрицы
void matrix_print(Matrix matrix_to_print);

// Функция для обмена строк матрицы
void swap_rows(Matrix* mat, int row1, int row2);
*/