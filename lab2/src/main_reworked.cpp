#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <iomanip>
#include <limits>
#include <string>

// Определяет класс аппроксимируемой функции и её аналитическое решение
// A * d^2u/dx^2 + B * du/dx + C = 0
// Изменения параметров требуют пересчёта аналитического решения!!!
class ApproxFunc {
  public:
    double a = 3.;  // du2/d2x
    double b = 2.;  // du/dx
    double c = -7.; // свободный член

    // Использование иных ГУ, кроме 1 и 2 требует доработки FEM::apply_border_conditions
    int bcl = 1; // Слева задано граничнео условие 1 рода (Дирихле)
    double bcl_value = -10.; // u(-2) = -10
    int bcr = 2; // Справа задано ГУ 2 рода (Неймана)
    double bcr_value = 1.; // du/dx (x=5) = u

    double xl = -2.; // координата левой границы
    double xr = 5.;  // координата правой границы

    double analytical_solution(double x) {
      return -1040242697543338869./73384984777978471.+216172782113783808./73384984777978471.*exp(-x*2./3.) + x*7./2.;
      //return -(21*exp(14./3.)*(x-4.)+66*exp(-2./3. * (x-5.))-35.*x + 30.) / (10.-6*exp(14./3.));
    }
};

class GaussSolver {
  public:
    // 0 - СЛАУ несовместна
    // -n - бесконечное количество решений (n - установленный ранг матрицы)
    // n - СЛАУ имеет единственное решение (n - установленный ранг матрицы)
    static int gauss_elimination(std::vector<std::vector<double>>& mat) {
        gauss_forward_elimination(mat);
        int res = gauss_check_solutions(mat);
        if (res != 0)
            gauss_backward_elimination(mat);
        return res;
    }

    // Прямой ход для приведения к верхнетреугольному виду
    static void gauss_forward_elimination(std::vector<std::vector<double>>& mat) {
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
    static void gauss_backward_elimination(std::vector<std::vector<double>>& mat) {
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
    static int gauss_check_solutions(std::vector<std::vector<double>>& mat) {
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

// Для порядка 2 и выше используется операция редукции, которая требует 
// TODO обратного хода Гаусса для получения значений внутри
class FEM {
  public:
    FEM(ApproxFunc f_, long fe_count_, int fe_type_)
      : f(f_), fe_count(fe_count_), fe_type(fe_type_),
      rigidity_matrix(static_cast<std::size_t>(fe_count_*fe_type_ + 1),
      std::vector<double>(static_cast<std::size_t>(fe_count_*fe_type_ + 2), 0.0))
      //L((f.xr - f.xl) / fe_count)
    {}

    long fe_count; // Количество конечных элементов
    int fe_type;   // Степень функции формы (1 - линейная, 3 - кубическая)
    std::vector<std::vector<double>> rigidity_matrix;
    //double L; // L = (f.xr - f.xl) / fe_count
    ApproxFunc f;

    int Solve() {
      ensembling();
      apply_border_conditions();
      //printLatexAugmented(std::out);
      // решения будут располагаться в дополненной части матрицы
      int solving_result_status = GaussSolver::gauss_elimination(rigidity_matrix);
      //printLatexAugmented(std::out);
      return solving_result_status;
    }

    // Печать для gnuplot
    void write_gnuplot(const std::string &prefix) {
      int N = static_cast<int>(rigidity_matrix.size());       
      int aug_col = static_cast<int>(rigidity_matrix[0].size()) - 1; 
      double dx = (f.xr - f.xl) / double(N - 1);
      std::ofstream fem_ofs(prefix + "_fem.dat");
      std::ofstream an_ofs(prefix + "_analytic.dat");
      fem_ofs << std::setprecision(12);
      an_ofs  << std::setprecision(12);
      for (int i = 0; i < N; ++i) {
          double x = f.xl + i * dx;
          double u_num = rigidity_matrix[i][aug_col];       
          double u_exact = f.analytical_solution(x);
          fem_ofs << x << " " << u_num << "\n";
          an_ofs  << x << " " << u_exact << "\n";
      }
      fem_ofs.close();
      an_ofs.close();
    }

    // Печать в .csv
    void write_results_csv(const std::string &filename) {
      double max_ae = 0.;
      int N = static_cast<int>(rigidity_matrix.size());
      if (N == 0) return;
      int aug_col = static_cast<int>(rigidity_matrix[0].size()) - 1;
      double dx = (f.xr - f.xl) / double(N - 1);

      std::ofstream ofs(filename + std::string("_table.csv"));
      if (!ofs.is_open()) return;

      // Заголовок (в кавычках, как вы просили)
      ofs << "\"x\",\"точное решение\",\"МКЭ\",\"абсолютная погрешность\"\n";
      ofs << std::setprecision(14);

      for (int i = 0; i < N; ++i) {
          double x = f.xl + i * dx;
          double u_exact = f.analytical_solution(x);
          double u_num   = rigidity_matrix[i][aug_col];
          double abs_err = std::isnan(u_num) ? std::numeric_limits<double>::quiet_NaN()
                                            : std::fabs(u_exact - u_num);
          if (abs_err > max_ae) {
            max_ae = abs_err;
          }

          auto write_quoted = [&](double v) {
              if (std::isnan(v)) {
                  ofs << "\"nan\"";
              } else {
                  // формируем строку с фиксированной точностью, затем заключаем в кавычки
                  std::ostringstream ss;
                  ss << std::setprecision(14) << v;
                  ofs << "" << ss.str() << "";
              }
          };

          write_quoted(x);        ofs << ';';
          write_quoted(u_exact);  ofs << ';';
          write_quoted(u_num);    ofs << ';';
          write_quoted(abs_err);  ofs << '\n';
      }
      ofs.close();
      //std::cout << std::fixed << std::setprecision(14) << max_ae << std::endl;

  }

    // Печать в LaTeX: дополненная матрица [A | b]
    void printLatexAugmented(std::ostream &out, int precision = 6) const {
        const size_t rows = rigidity_matrix.size();
        if (rows == 0) {
            out << "% Пустая матрица\n";
            return;
        }
        const size_t cols = rigidity_matrix[0].size(); // включая последний столбец (вектор b)

        // Формируем LaTeX-спецификацию колонок: "ccc|c"
        std::string colspec;
        colspec.reserve(cols + 2);
        for (size_t j = 0; j < cols - 1; ++j) colspec += 'c';
        colspec += "|c";

        out << "\\[\n";
        out << "\\begin{array}{" << colspec << "}\n";

        out << std::fixed << std::setprecision(precision);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                out << rigidity_matrix[i][j];
                if (j + 1 < cols) out << " & ";
            }
            if (i + 1 < rows) out << " \\\\ \n";
            else out << "\n";
        }

        out << "\\end{array}\n";
        out << "\\]\n";
    }


  private:
    std::vector<std::vector<double>> linear_rigidity_matrix (double L) {
      std::vector<std::vector<double>> lrm = {
            { - f.a / L - f.b / 2., f.a / L + f.b / 2. },
            { f.a / L - f.b / 2., - f.a / L + f.b / 2. },
      };
      //std::cout << lrm[0][0] << " " << lrm[0][1] << std::endl <<  lrm[1][0] << " " << lrm[1][1] << std::endl;
      return lrm;
    }
    std::vector<double> linear_loads_vector (double L) {
      std::vector<double> llv = {
        - f.c * L / 2,
        - f.c * L / 2,
      };
      return llv;
    }

    std::vector<std::vector<double>> cubic_rigidity_matrix (double L) {
      std::vector<std::vector<double>> crm = {
        {     - f.b/2 - (37*f.a)/(10*L), (57*f.b)/80 + (189*f.a)/(40*L), - (3*f.b)/10 - (27*f.a)/(20*L),   (7*f.b)/80 + (13*f.a)/(40*L)},
        {(189*f.a)/(40*L) - (57*f.b)/80,              -(54*f.a)/(5*L), (81*f.b)/80 + (297*f.a)/(40*L), - (3*f.b)/10 - (27*f.a)/(20*L)},
        {  (3*f.b)/10 - (27*f.a)/(20*L), (297*f.a)/(40*L) - (81*f.b)/80,              -(54*f.a)/(5*L), (57*f.b)/80 + (189*f.a)/(40*L)},
        {  (13*f.a)/(40*L) - (7*f.b)/80,   (3*f.b)/10 - (27*f.a)/(20*L), (189*f.a)/(40*L) - (57*f.b)/80,        f.b/2 - (37*f.a)/(10*L)}
      };
      return crm;
    }
    std::vector<double> cubic_loads_vector (double L) {
      std::vector<double> clv = {
        - f.c * L / 8.,
        - 3. * f.c * L / 8.,
        - 3. * f.c * L / 8.,
        - f.c * L / 8.,
      };
      return clv;
    }

    // Ансамблирование
    void ensembling () {
      std::vector<std::vector<double>> rm;
      std::vector<double> lv;
      if (fe_type == 1) {
        rm = linear_rigidity_matrix((f.xr - f.xl) / fe_count);
        lv = linear_loads_vector((f.xr - f.xl) / fe_count);
        for (long i = 0; i < fe_count; i++) {
          rigidity_matrix[i][i] += rm[0][0]; rigidity_matrix[i][i+1] += rm[0][1];
          rigidity_matrix[i+1][i] += rm[1][0]; rigidity_matrix[i+1][i+1] += rm[1][1];
          rigidity_matrix[i][fe_count+1] += lv[0];
          rigidity_matrix[i+1][fe_count+1] += lv[1];
        }
      }
      else {
        rm = cubic_rigidity_matrix((f.xr - f.xl) / fe_count);
        lv = cubic_loads_vector((f.xr - f.xl) / fe_count);
        for (long i = 0; i < fe_count*fe_type; i+=fe_type) {
          rigidity_matrix[i][i] += rm[0][0]; rigidity_matrix[i][i+1] += rm[0][1]; rigidity_matrix[i][i+2] += rm[0][2]; rigidity_matrix[i][i+3] += rm[0][3];
          rigidity_matrix[i+1][i] += rm[1][0]; rigidity_matrix[i+1][i+1] += rm[1][1]; rigidity_matrix[i+1][i+2] += rm[1][2]; rigidity_matrix[i+1][i+3] += rm[1][3];
          rigidity_matrix[i+2][i] += rm[2][0]; rigidity_matrix[i+2][i+1] += rm[2][1]; rigidity_matrix[i+2][i+2] += rm[2][2]; rigidity_matrix[i+2][i+3] += rm[2][3];
          rigidity_matrix[i+3][i] += rm[3][0]; rigidity_matrix[i+3][i+1] += rm[3][1]; rigidity_matrix[i+3][i+2] += rm[3][2]; rigidity_matrix[i+3][i+3] += rm[3][3];
          rigidity_matrix[i][fe_count*fe_type+1] += lv[0];
          rigidity_matrix[i+1][fe_count*fe_type+1] += lv[1];
          rigidity_matrix[i+2][fe_count*fe_type+1] += lv[2];
          rigidity_matrix[i+3][fe_count*fe_type+1] += lv[3];
        }
      } 
      return;
    }

    // Реализованы ГУ 1 и 2 рода для левой и правой границы рассматриваемого отрезка
    void apply_border_conditions() {
      if (f.bcl == 1) {
        rigidity_matrix[0][0] = 1.;
        rigidity_matrix[0][1] = 0.;
        if (fe_type == 3) {
          rigidity_matrix[0][2] = 0.;
          rigidity_matrix[0][3] = 0.;
        }
        rigidity_matrix[0][fe_count*fe_type+1]  = f.bcl_value;
      } 
      //TODO протестировать
      else{;}

      if (f.bcr == 2) {
        rigidity_matrix[fe_count*fe_type][fe_count*fe_type] += f.a * f.bcr_value;
      } 
      //TODO протестировать
      else{;} 
      return;
    }
};

// Определяет параметры моделирования, задаваемые пользователем
struct ModelParameters {
    long fe_count;          // Количество конечных элементов
    int fe_type;            // Степень функции формы (1 - линейная, 3 - кубическая)
    std::string output_file; // Название файла для вывода результатов моделирования
};

// Парсит параметры, задаваемые пользователем и записывает в объект класса ModelParameters
class ParameterParser {
    public:
    ParameterParser(int c, char** v) : _argument_count(c), _argument_strings(v) {}
  
    // Разбирает аргументы, переданные при создании, для определения параметров моделирования.
    // выбрасывает std::invalid_argument при некорректных значениях параметров.
    ModelParameters getParameters() const {
      // Значения по умолчанию.
      long fe_count = 3;
      int fe_type = 1;
      std::string output_file = "fe_result";

      // Разбор аргументов при помощи getopt.
      int option = 0;
      while ((option = getopt(_argument_count, _argument_strings, "c:t:o:")) != -1) {
        switch (option) {
        case 't':
          fe_type = atoi(optarg);
          if (fe_type != 1 && fe_type != 3) {
            throw std::invalid_argument("ModelParameters: -t (finite element type) must be 1 (linear) or 3 (cubic).");
          }
          break;
        case 'c':
            fe_count = atoi(optarg);
            if (fe_count < 1) {
                throw std::invalid_argument("ModelParameters: -c (the number of finite elements) must be at least 1");
              }
            break;
        case 'o':
          output_file = optarg;
          break;
        default:
          break;
        }
      }

      ModelParameters parameters;
      parameters.output_file = output_file;
      parameters.fe_count = fe_count;
      parameters.fe_type  = fe_type;

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
        ApproxFunc f;
        FEM fem_model(f, params.fe_count, params.fe_type);
        fem_model.Solve();
        fem_model.write_gnuplot(params.output_file);
        fem_model.write_results_csv(params.output_file);
    } catch (...) {
        return 1;
    }
    std::cout << "Расчёт проведен успешно\n";
    return 0;
}


/*
{
-(3*f.c*L*(414*f.a + 19*f.b*L))/(8*(126*f.a + 19*f.b*L)),
(1080*duLi*f.a^3 + 360*duLi*f.a^2*f.b*L - 2160*f.c*f.a*f.a*L + 38*duLi*f.a*f.b^2*L^2 - 180*f.c*f.a*f.b*L^2 - 31*f.c*f.b^2*L^3)/(2*(540*f.a^2 + 180*f.a*f.b*L + 19*f.b^2*L^2)),
-(1080*duLh*f.a^3 + 360*duLh*f.a^2*f.b*L - 1548*f.c*f.a*f.a*L + 38*duLh*f.a*f.b^2*L^2 + 110*f.c*f.a*f.b*L^2 + 7*f.c*f.b^2*L^3)/(2*(540*f.a^2 + 180*f.a*f.b*L + 19*f.b^2*L^2))
}
*/

/*
{(57*f.b)/80 + (189*f.a)/(40*L),                                      - (3*f.b)/10 - (27*f.a)/(20*L),                                                                               - f.b/2 - (37*f.a)/(10*L),                                                                           (7*f.b)/80 + (13*f.a)/(40*L)},
{                         0, (81*(540*f.a^2 + 180*f.a*f.b*L + 19*f.b^2*L^2))/(80*L*(126*f.a + 19*f.b*L)),                                 -(3*(12540*f.a^2 + 3840*f.a*f.b*L + 361*f.b^2*L^2))/(80*L*(126*f.a + 19*f.b*L)),                                    -(3*(510*f.a^2 + 255*f.a*f.b*L + 38*f.b^2*L^2))/(20*L*(126*f.a + 19*f.b*L))},
{                         0,                                                               0,      -(27*(120*f.a^3 + 60*f.a^2*f.b*L + 12*f.a*f.b^2*L^2 + f.b^3*L^3))/(2*L*(540*f.a^2 + 180*f.a*f.b*L + 19*f.b^2*L^2)),      (27*(120*f.a^3 + 60*f.a^2*f.b*L + 12*f.a*f.b^2*L^2 + f.b^3*L^3))/(2*L*(540*f.a^2 + 180*f.a*f.b*L + 19*f.b^2*L^2))},
{                         0,                                                               0, (55380*f.a^3 + 13020*f.a^2*f.b*L + 417*f.a*f.b^2*L^2 - 110*f.b^3*L^3)/(20*L*(540*f.a^2 + 180*f.a*f.b*L + 19*f.b^2*L^2)), -(16620*f.a^3 + 2820*f.a^2*f.b*L - 181*f.a*f.b^2*L^2 - 55*f.b^3*L^3)/(10*L*(540*f.a^2 + 180*f.a*f.b*L + 19*f.b^2*L^2))}
*/

/*

*/