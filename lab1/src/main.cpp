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

void print_augmented_matrix(const std::vector<std::vector<double>>& A, 
    int max_show = 22, 
    int precision = 3) {
const int rows = A.size();
if(rows == 0) {
std::cout << "Empty matrix!\n";
return;
}
const int cols = A[0].size();
const int shown_cols = std::min(cols, max_show + 1); // +1 для RHS

// Заголовок
std::cout << "\n\033[1;36m" << std::setw(12) << "Equation │";
for(int j = 0; j < shown_cols - 1; ++j) {
std::cout << std::setw(10) << "a" + std::to_string(j);
}
std::cout << std::setw(12) << "│ RHS";
std::cout << "\n\033[0m";

// Вывод данных
for(int i = 0; i < std::min(rows, max_show); ++i) {
// Номер строки
std::cout << std::setw(10) << i << " │";

// Коэффициенты
for(int j = 0; j < shown_cols; ++j) {
if(j == cols - 1) std::cout << " │"; // Разделитель для RHS

if(fabs(A[i][j]) < 1e-12) {
std::cout << std::setw(10) << "·";
} else {
std::cout << std::fixed << std::setprecision(precision) 
  << std::setw(10) << A[i][j];
}
}

// Многоточие для больших матриц
if(cols > shown_cols) std::cout << " ...";
std::cout << "\n";
}

// Футер
if(rows > max_show || cols > shown_cols) {
std::cout << "\n\033[3mShowing first " << max_show 
<< " rows/columns. Total size: " 
<< rows << "x" << cols << "\033[0m\n";
}
}


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
               //std::cout << i << "," << j << "," << x << "," << y << "," << bc_type << "," << bc_val << "," << std::endl;
            }
            else{
               nodes.emplace_back(i, j, x, y, bc_type, bc_val, bc_val);
            }
            // std::cout << i << "," << j << std::endl;
        }
        fin.close();
/*
        for (auto& node : nodes) {
            int i = node.i;
            int j = node.j;
            std::cout << i << "," << j << std::endl;
        }
  */ 
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
            const double alpha = 1;                    // Коэффициент теплопроводности
            const int N = mesh.nodes.size();             // Число уравнений в СЛАУ будет соответсвтовать числу узлов, в которых не граничное условие 1 рода
            const double dt = 0.1;                         // Шаг по времени
            int timesteps = (int) time_final / dt;       // Инициализация начального момента времени
            //int timesteps = 1;
            const int nodes_count = mesh.nodes.size();

            // A B = C, B = A^-1 C
            std::vector<std::vector<double>> A(N, std::vector<double>(N+1));
            OrderIndexMap<int> B(N); 
            for (auto& node : mesh.nodes) {
                int tmp = B.get_or_add(node.j * mesh.nCols + node.i);
            }

            for (int k = 0; k < timesteps; k++) {
                for (auto& row : A) std::fill(row.begin(), row.end(), 0.0);

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

                            A[n][N] = mesh.dx*mesh.dx*mesh.dy*mesh.dy * mesh.nodes[b_poz_i_j].t; 

                            A[n][b_poz_i_j] = (mesh.dx*mesh.dx*mesh.dy*mesh.dy) \
                                              + 2*alpha*mesh.dy*mesh.dy*dt \
                                              + 2*alpha*mesh.dx*mesh.dx*dt;

                            A[n][b_poz_im_j] = -alpha*mesh.dy*mesh.dy*dt;
                            A[n][b_poz_ip_j] = -alpha*mesh.dy*mesh.dy*dt;

                            A[n][b_poz_i_jm] = -alpha*mesh.dx*mesh.dx*dt;
                            A[n][b_poz_i_jp] = -alpha*mesh.dx*mesh.dx*dt;                          
                            break;
                        case 1:
                            b_poz = B.get_or_add(node.j * mesh.nCols + node.i);
                            A[n][b_poz] = 1.;
                            A[n][N] = node.bc.val;
                            break;
                        case 2:
                            b_poz = B.get_or_add(node.j * mesh.nCols + node.i);
                            b_poz_n = B.get_or_add(node.bc.n_j * mesh.nCols + node.bc.n_i);
                            A[n][b_poz] = -1.;
                            A[n][b_poz_n] = 1.;
                            A[n][N] = node.bc.val * sqrt((mesh.dx*mesh.dx*(node.bc.n_i-node.i)*(node.bc.n_i-node.i) \
                                        + mesh.dy*mesh.dy*(node.bc.n_j-node.j)*(node.bc.n_j-node.j)));
                            break;
                        case 3:
                            b_poz = B.get_or_add(node.j * mesh.nCols + node.i);
                            b_poz_n = B.get_or_add(node.bc.n_j * mesh.nCols + node.bc.n_i);
                            A[n][b_poz] = -(node.bc.val * sqrt((mesh.dx*mesh.dx*(node.bc.n_i-node.i)*(node.bc.n_i-node.i) \
                                                                + mesh.dy*mesh.dy*(node.bc.n_j-node.j)*(node.bc.n_j-node.j))) + 1);
                            A[n][b_poz_n] = 1.;
                            A[n][N] = 0.;
                            break;
                    }
                }

                // решение расширенной СЛАУ
                gauss_elimination(A);

                for (auto& node : mesh.nodes) {
                    if (node.bc.type != 1) {
                        node.t = A[B.get_or_add(node.j * mesh.nCols + node.i)][N];
                    }
                }

                /*альтернативный решатель СЛАУ (быстрее)*/
                /*
                Eigen::MatrixXd eigenA(N, N);
                Eigen::VectorXd eigenC(N);
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        eigenA(i, j) = A[i][j];
                    }
                    eigenC(i) = A[i][N];
                }

                Eigen::VectorXd solution = eigenA.partialPivLu().solve(eigenC);

                for (auto& node : mesh.nodes) {
                    if (node.bc.type != 1) {
                        node.t = solution(B.get_or_add(node.j * mesh.nCols + node.i));
                    }
                }
                */

            }
            return 0;
        }

    private:
        // 0 - СЛАУ несовместна
        // -n - бесконечное количество решений (n - установленный ранг матрицы)
        // n - СЛАУ имеет единственное решение (n - установленный ранг матрицы)
        int gauss_elimination(std::vector<std::vector<double>>& mat) {
            gauss_forward_elimination(mat);
            int res = gauss_check_solutions(mat);
            if (res != 0)
                gauss_backward_elimination(mat);
            return res;
        }

        // Прямой ход для приведения к верхнетреугольному виду
        void gauss_forward_elimination(std::vector<std::vector<double>>& mat) {
            for (int k = 0; k < mat.size(); k++) {
                /* поиск строки с максимальным абсолютным значением в текущем столбце 
                   (от строки k вниз), чтобы избежать деления на слишком маленькие числа 
                   и улучшить численную стабильность алгоритма. */
                   int max_row = k;
                   for (int i = k + 1; i < mat.size(); i++) {
                       if (fabs(mat[i][k]) > fabs(mat[max_row][k])) {
                           max_row = i;
                       }
                   }
           
                   // Меняем строки, если необходимо
                   if (max_row != k) {
                       mat[k].swap(mat[max_row]);
                   }

                   // Нормализация основной строки
                   double pivot = mat[k][k];                 
                   for (int j = k; j < mat[k].size(); j++) {
                       mat[k][j] /= pivot;
                   }
                   
                   // Обнуление элементов под основной строкой
                   for (int i = k + 1; i < mat.size(); i++) {
                       double factor = mat[i][k];
                       mat[i][k] = 0.;
                       for (int j = k+1; j < mat[k].size(); j++) {
                           mat[i][j] -= factor * mat[k][j];
                       }
                   }
               }
           }
           
           // Обратный ход для приведения к диагональному виду
           void gauss_backward_elimination(std::vector<std::vector<double>>& mat) {
               for (int k = mat.size() - 1; k >= 0; k--) {
                   for (int i = k - 1; i >= 0; i--) {
                       double factor = mat[i][k];
                       mat[i][k] = 0.;
                       for (int j = k + 1; j < mat[i].size(); j++) {
                           mat[i][j] -= factor * mat[k][j];
                       }
                   }
               }
           }

            // Функция для проверки количества решений СЛАУ
            // Принимает на вход матрицу верхнетреугольного вида
            // 0 - СЛАУ несовместна (!!определяется только если переданная матрица была расширенной!!) (mat[i].size() - 1 поэтому всегда считается расширенной всегда)
            // -m - бесконечное количество решений (m - установленный ранг матрицы)
            // n - СЛАУ имеет единственное решение (n - установленный ранг матрицы)
        int gauss_check_solutions(std::vector<std::vector<double>>& mat) {
            int rank = 0;
         
            for (int i = 0; i < mat.size(); i++) {
                bool row_nonzero = false;
                for (int j = 0; j < mat[i].size() - 1; j++) {
                    if (fabs(mat[i][j]) > 1e-12) {  // Пороговое значение для сравнения с нулем
                        row_nonzero = true;
                        break;
                    }
                }
                // 0 0 0 0 | 4
                if (!row_nonzero && fabs(mat[i].back()) > 1e-12)
                    return 0;

                if (row_nonzero) rank++;
            }
        
            if (rank == mat.size())
                return rank;
            else
                return -rank;     
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