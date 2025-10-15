#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <string>
#include <thread>
#include <chrono>

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
        else {
            std::cout << "SOS! СЛАУ НЕСОВМЕСТНА" << std::endl;
        }
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
            if (!row_nonzero && fabs(mat[i].back()) > 1e-24)
                return 0;

            if (row_nonzero) rank++;
        }
    
        if (rank == mat.size())
            return rank;
        else
            return -rank;     
    }       
};

class Schema9 {
    public:
        // Параметры схемы
        const double I = 1.;
        const double L = 0.001;
        const double C = 1e-8;
        const double rb = 20.;
        const double ru = 1000000.;
        const double Cb = 2e-12; 
        const double It = 1e-12;
        const double R = 1000.;
        const double E = 100.;
        const double IL_0 = 0.;
        const double UC_0 = 0.;
        const double UCb_0 = 0.;
        const double mft = 0.026;
        
        // Параметры времени расчёта
        const double t_start = 0.;
        const double t_end = 0.001;

        // Гиперпараметры решателя
        const double delta_t_min = 1e-13;
        const double delta_t_start = 1e-9;
        const double eps = 1e-2;
        const double eps1 = 1e-6;
        const double eps2 = 1;
        const int max_step = 7;

                // ---------- заранее выделяемые (статические для экземпляра) матрицы/векторы ----------
        // Матрица 6x7 (расширенная матрица) — хранится как vector<vector<double>> для совместимости с GaussSolver
        std::vector<std::vector<double>> augm;

        std::vector<double> _prev_state;      // size 3: IL_nm1, UC_nm1, UCb_nm1
        std::vector<double> _init_approx;     // size 6: IL, phi1, phi2, phi3, phi4, IE
        std::vector<double> _X;               // size 6
        std::vector<double> _X_prev;          // size 6
        std::vector<double> _X_prev_prev;     // size 6
        std::vector<double> _Xdtdt;           // size 6
        std::vector<double> _dX;              // size 6

        // создание матрицы
        // принимает: double deltaT, 
                   // <double IL, double phi1, double phi2, double phi3, double phi4, double IE>, 
                   // <double IL_nm1, double UC_nm1, UCb_nm1>
        // инициализирует Y|I матрицу (см. отчёт)
        void calculate_augmented_matrix(double deltaT, std::vector<double>& init_approx, std::vector<double>& prev_state_variables, std::vector<std::vector<double>>& augm) {
            double IL_nm1 = prev_state_variables[0];
            double UC_nm1 = prev_state_variables[1];
            double UCb_nm1 = prev_state_variables[2];
            double IL = init_approx[0], IE = init_approx[5];
            double phi1 = init_approx[1]; double phi2 = init_approx[2]; double phi3 = init_approx[3]; double phi4 = init_approx[4]; 
            double alpha = It*(std::exp((phi2-phi3)/mft))/mft; 
            double Id = It*(std::exp((phi2-phi3)/mft)-1.);
            augm[0][0] = -L/deltaT; 
            augm[0][1] = 1.; 
            augm[0][2] = 0.; 
            augm[0][3] = 0.; 
            augm[0][4] = 0.; 
            augm[0][5] = 0.;
            augm[0][6] = -(phi1 - L/deltaT*(IL-IL_nm1));

            augm[1][0] = 1.; 
            augm[1][1] = C/deltaT + 1./rb; 
            augm[1][2] = -1./rb; 
            augm[1][3] = 0.; 
            augm[1][4] = 0.; 
            augm[1][5] = 0.;
            augm[1][6] = -(-I + IL + C/deltaT*(phi1 - UC_nm1) + (phi1-phi2)/rb);

            augm[2][0] = 0.; 
            augm[2][1] = -1./rb;
            augm[2][2] = Cb/deltaT + 1./rb + 1./ru + alpha;
            augm[2][3] = -Cb/deltaT - 1./ru - alpha;
            augm[2][4] = 0.; 
            augm[2][5] = 0.;
            augm[2][6] = -( -(phi1-phi2)/rb + Cb/deltaT*(phi2-phi3-UCb_nm1) + (phi2-phi3)/ru + Id );

            augm[3][0] = 0.; 
            augm[3][1] = 0.;
            augm[3][2] = -Cb/deltaT - 1./ru - alpha;
            augm[3][3] = Cb/deltaT + 1./ru + alpha;
            augm[3][4] = 0.; 
            augm[3][5] = 1.;
            augm[3][6] = - ( - Cb/deltaT*(phi2-phi3-UCb_nm1) - Id - (phi2-phi3)/ru + IE );

            augm[4][0] = 0.; 
            augm[4][1] = 0.; 
            augm[4][2] = 0.; 
            augm[4][3] = 0.;
            augm[4][4] = 1./R; 
            augm[4][5] = -1.;
            augm[4][6] = - ( -IE + phi4/R );

            augm[5][0] = 0.; 
            augm[5][1] = 0.; 
            augm[5][2] = 0.;
            augm[5][3] = -1.; 
            augm[5][4] = 1.; 
            augm[5][5] = 0.;
            augm[5][6] = - ( E - (phi3 - phi4) );


            /*
            std::cout << "augm = \n";
            for (int i = 0; i < augm.size(); i++) {
                std::cout << "<";
                for (int j = 0; j < augm[0].size(); j++) {
                    std::cout << augm[i][j] << ",";
                }
                std::cout << ">\n";
            }
            std::cout << std::fflush;
            */
        
            return;
        }

        // Принимает double deltaT, 
                    // <double IL, double phi1, double phi2, double phi3, double phi4, double IE>, 
                    // <double IL_nm1, double UC_nm1, UCb_nm1>
        // Возвращает вектор deltas <IL, phi1, phi2, phi3, phi4, IE>
        bool newton_iteration (double DeltaT, std::vector<double>& init_approx, std::vector<double>& prev_state_variables, std::vector<std::vector<double>>& augm, std::vector<double>& dX) {
            calculate_augmented_matrix(DeltaT, init_approx, prev_state_variables, augm);
            int rang = GaussSolver::gauss_elimination(augm);
            if (rang != augm.size()) {
                return false;
            }
            for (size_t i = 0; i < 6; ++i) 
                dX[i] = augm[i][6];
            return true;
        }

        double calculate_norm(std::vector<double>& deltas) {
            if (deltas.empty()) {
                return 0.0;
            }
            
            double max_abs = std::abs(deltas[0]);
            for (size_t i = 1; i < deltas.size(); ++i) {
                max_abs = std::max(max_abs, std::abs(deltas[i]));
            }
            return max_abs;
        }

    bool solve_schema() {
        double t = t_start;
        double t_prev = t_start;
        double dt = delta_t_start;
        double dt_prev = delta_t_start; // Нельзя 0, иначе при расчёте q будет деление на 0
        double norm_dX;

        // подготовка лог-файла для вывода данных, удобных для GnuPlot
        std::ofstream flog("data.dat");
        if (!flog.is_open()) {
            std::cerr << "Не удалось открыть файл data.dat для записи\n";
            return false;
        }
        //flog << "# idx t phi1 phi2 phi3 phi4\n";
        flog << "# idx t phi1 phi3 phi4\n";

        int count = 0;

        // Инициализация состояний
        std::vector<std::vector<double>> augm(6, std::vector<double>(7, 0.0));
        std::vector<double> dX(6, 0.0);
        std::vector<double> X(6, 0.0);
        std::vector<double> X_prev(6, 0.0); // <IL, phi1, phi2, phi3, phi4, IE> (-1) - приемняется при расчёте второй производной 
        std::vector<double> X_prev_prev(6, 0.0); // <IL, phi1, phi2, phi3, phi4, IE> (-2)
        std::vector<double> init_approx(6, 0.0); // IL, phi1, phi2, phi3, phi4, IE
        std::vector<double> prev_state = {IL_0, UC_0, UCb_0}; // IL, UC, UCb
        std::vector<double> Xdtdt(6, 0.0); // используется при расчёте q

        while (t < t_end) {
            std::cout << t << std::endl;
            
            //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            bool step_accepted = false;
            bool restart_with_smaller_dt = false;
            // Пытаемся выполнить текущий временной шаг, пока он не будет принят
            while (!step_accepted) {
                // TODO ВАЖНО! при превышлении максимального числа итераций методом Ньютона мы изменяем шаг и следовательно начальное приближение пересчитываем!!!
                for (size_t i = 0; i < init_approx.size(); ++i) {
                    init_approx[i] = 2.0 * X_prev[i] - X_prev_prev[i];
                }
                int n_iter = 0;

                do {
                    n_iter++;
                    if (!newton_iteration(dt, init_approx, prev_state, augm, dX)) {
                        std::cerr << "Ошибка во время решения СЛАУ\n";
                        flog.close();
                        return false;
                    }
                    //std::cout << "dX=<" << dX[0] << "," << dX[1] << "," << dX[2] << "," << dX[3] << "," << dX[4] << "," << dX[5] << ">" << std::endl << std::fflush;
                    for (size_t i = 0; i < init_approx.size(); i++)
                        init_approx[i] += dX[i];
                    norm_dX = calculate_norm(dX);
                    //std::cout << "norm_dX=" << norm_dX << std::endl << std::fflush;

                    // Превышение числа итераций — уменьшаем dt и перезапускаем шаг
                    if (n_iter > max_step) {                
                        dt /= 2.0;
                        if (dt < delta_t_min) {
                            std::cerr << "Ошибка: сходимость не достигается (dt < dt_min)\n";
                            std::cout << "t=" << t << std::endl << std::fflush;
                            flog.close();
                            return false;
                        }
                        restart_with_smaller_dt = true;
                        break; // выход из цикла do...while
                    }
                } while (norm_dX > eps);

                if (restart_with_smaller_dt) {
                    // попытаемся снова с уменьшенным dt (внешний while не инкрементирует t)
                    continue;
                }

                for (size_t i = 0; i < X.size(); i++) {
                    Xdtdt[i] = ((init_approx[i]-X_prev[i])/dt-(X_prev[i]-X_prev_prev[i])/dt_prev)/dt;
                }

                double q = 0.5*dt*calculate_norm(Xdtdt)*dt;
                std::cout << "dt=" << dt << std::endl;
                std::cout << "q=" << q << std::endl;

                // === Адаптивное изменение шага ===
                if (q < eps1) { 
                    //std::cout << "eps1" << std::endl << std::fflush;
                    t_prev = t;
                    t += dt;
                    dt_prev = dt;
                    dt *= 2.0;
                } else {
                    if (q < eps2) {
                        //std::cout << "sigma2" << std::endl << std::fflush;
                        t_prev = t;
                        t += dt;
                        dt_prev = dt;
                        dt = dt;
                    }
                    else {
                        // TODO?????? КОСТЫЛЬ
                            t_prev = t;
                            t += dt;
                            dt_prev = dt;
                            dt = dt;
                        /*
                        dt /= 2;
                        if (dt < delta_t_min) {
                            std::cerr << "Ошибка: сходимость не достигается (dt < dt_min)\n";
                            std::cout << "t=" << t << std::endl << std::fflush;
                            flog.close();
                            return false;
                        }
                        continue;*/ 
                    }
                }

                // === Успешное завершение шага ===
                for (size_t i = 0; i < X.size(); i++) {
                    X_prev_prev[i] = X_prev[i];
                    X_prev[i] = X[i];
                    X[i] = init_approx[i];
                }
            
                // Запоминаем значения для следующего шага
                prev_state[0] = X[0]; // IL
                prev_state[1] = X[1]; // UC = phi1
                prev_state[2] = X[2] - X[3]; // UCb = phi2 - phi3

                flog << std::fixed << std::setprecision(12)
                << count++ << ' '
                << t      << ' '
                << X[1]   << ' '
                << X[3]   << ' '
                << X[4]   << '\n';

                step_accepted = true; // выходим из внутреннего while и переходим к следующему временному шагу
            } // конец попыток одного временного шага
        } // конец цикла по времени
        flog.close();
        return true;
    }
};

int main(int argc, char** argv) {
    Schema9 dz;
    if(dz.solve_schema())
        std::cout << "Расчёт проведен успешно\n";
    return 0;
}