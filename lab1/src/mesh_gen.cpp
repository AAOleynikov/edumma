#include <iostream> 
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <string>

// Структура узла КЭ сетки
struct Node {
    int i, j;       // индексы узла: i – строка (по y), j – столбец (по x)
    double x, y;    // координаты узла в пространстве
    // 0 - начальные условия, T(t=0) = val
    // 1 - Dirichlet,         T      = val
    // 2 - Neumann,           dT/dn  = val
    // 3 - Robin,             dT/dn  = val * T
    int bc_type;  
    double bc_val;
    int n_i, n_j;   // определяет "узел нормали" при наличии такового 
};

// Функция возвращает true, если точка (x,y) лежит строго внутри треугольника с вершинами (1,1), (1,3) и (3,1)
bool insideTriangle(double x, double y) {
    const double Ax = 1.0, Ay = 1.0;
    const double Bx = 1.0, By = 3.0;
    const double Cx = 3.0, Cy = 1.0;
    auto area = [](double ax, double ay, double bx, double by, double cx, double cy) -> double {
        return std::fabs((ax*(by - cy) + bx*(cy - ay) + cx*(ay - by)) / 2.0);
    };
    double A  = area(Ax, Ay, Bx, By, Cx, Cy);
    double A1 = area(x, y, Bx, By, Cx, Cy);
    double A2 = area(Ax, Ay, x, y, Cx, Cy);
    double A3 = area(Ax, Ay, Bx, By, x, y);
    const double tol = 1e-8;
    if (std::fabs((A1 + A2 + A3) - A) > tol)
        return false;
    if (A1 < tol || A2 < tol || A3 < tol)
        return false;
    return true;
}

// Функция возвращает true, если точка (x,y) лежит на границе треугольника.
int onTriangleEdge(double x, double y) {
    const double tol = 1e-8;
    if (std::fabs(x - 1.0) < tol && y >= 1.0 - tol && y <= 3.0 + tol)
        return 1;
    if (std::fabs(y - 1.0) < tol && x >= 1.0 - tol && x <= 3.0 + tol)
        return 2;
    if (std::fabs(x + y - 4.0) < tol && x >= 1.0 - tol && x <= 3.0 + tol && y >= 1.0 - tol && y <= 3.0 + tol)
        return 3;
    return 0;
}

// Класс генератора сетки для пластины.
class MeshGenerator {
public:
    MeshGenerator(double Lx, double Ly, int Nx, int Ny)
        : Lx(Lx), Ly(Ly), Nx(Nx), Ny(Ny)
    {
        dx = Lx / Nx;
        dy = Ly / Ny;
    }
    
    // Генерирует узлы сетки с автоматическим назначением граничных условий
    void generateMesh() {
        int counter = 0;
        // Перебор строк (i от 0 до Ny) и столбцов (j от 0 до Nx)
        for (int j = 0; j <= Ny; j++) {
            for (int i = 0; i <= Nx; i++) {
                double x = i * dx;
                double y = j * dy;

                // Если точка находится на границе треугольника, назначаем условие Robin: dT/dn = 1 * T
                if (int norm = onTriangleEdge(x, y)) {
                    Node node;
                    node.i = i;
                    node.j = j;
                    node.x = x;
                    node.y = y;
                    node.bc_type = 3;
                    node.bc_val = 1;
                    switch (norm){
                        case 1:
                            node.n_i = i - 1;
                            node.n_j = j;
                            break;
                        case 2:
                            node.n_i = i;
                            node.n_j = j - 1;     
                            break;
                        case 3:
                            node.n_i = i+1;
                            node.n_j = j+1;
                            break;
                        default:
                            break;
                    }

                    nodes.push_back(node);
                } else
                // Если вне треугольника и не на границах
                if (!insideTriangle(x, y) && (i!=0 && i!=Ny) && (j!=0 && j!=Nx)) {
                    Node node;
                    node.i = i;
                    node.j = j;
                    node.x = x;
                    node.y = y;
                    // По умолчанию назначаем начальное условие: T=30
                    node.bc_type = 0;
                    node.bc_val = 30.0;
                    nodes.push_back(node);
                }
            
                // граничные условия второго рода
                else if (j == 0 && i!=0 && i!=Nx) {
                    Node node;
                    node.i = i;
                    node.j = j;
                    node.x = x;
                    node.y = y;
                    node.bc_type = 2;
                    node.bc_val = 0;
                    node.n_i = i;
                    node.n_j = j+1;
                    nodes.push_back(node);                    
                }
                // граничные условия первого рода
                else if (i == 0 && j!=0 && j!=Ny) {
                    Node node;
                    node.i = i;
                    node.j = j;
                    node.x = x;
                    node.y = y;
                    node.bc_type = 1;
                    node.bc_val = 70;
                    nodes.push_back(node);                    
                }
                else if (i == Nx && j!=0 && j!=Ny) {
                    Node node;
                    node.i = i;
                    node.j = j;
                    node.x = x;
                    node.y = y;
                    node.bc_type = 1;
                    node.bc_val = 70;
                    nodes.push_back(node);                    
                }
                else if (j == Ny && i!=0 && i!=Nx) {
                    Node node;
                    node.i = i;
                    node.j = j;
                    node.x = x;
                    node.y = y;
                    node.bc_type = 1;
                    node.bc_val = 70;
                    nodes.push_back(node);                    
                }
            }
        }
    }

    
    // Экспортирует сетку в CSV файл "mesh_nodes.csv"
    // Первая строка: метаданные max_i, max_j, dx, dy
    // Перед каждой строкой данных добавляется комментарий
    void exportToFile() {
        std::ofstream fout("mesh_nodes.csv");
        fout << std::setprecision(8);
        // Метаданные: максимальные индексы и шаги
        fout << "# метаданные: nRows, nCols, dx, dy" << std::endl;
        fout << Nx + 1 << ", " << Ny + 1 << ", " << dx << ", " << dy << std::endl;

        fout << "# узел: i, j, x, y, bc_type, bc_val, n_i, n_j" << std::endl;
        // Данные узлов: i, j, x, y, bc_type, bc_val
        for (size_t k = 0; k < nodes.size(); ++k) {
            fout << nodes[k].i << ", "
                 << nodes[k].j << ", "
                 << nodes[k].x << ", "
                 << nodes[k].y << ", "
                 << nodes[k].bc_type << ", "
                 << nodes[k].bc_val;
            if (nodes[k].bc_type == 2 || nodes[k].bc_type == 3) {
                fout << ", "<< nodes[k].n_i << ", "
                     << nodes[k].n_j;
            }
            fout << std::endl;
        }
        fout.close();
    }
       
private:
    double Lx, Ly; // размеры области
    int Nx, Ny;    // число разбиений по x и y
    double dx, dy; // шаги по осям
    std::vector<Node> nodes;                   // вектор узлов сетки
};

int main() {
    double Lx = 4.0, Ly = 4.0; // длина и ширина пластины
    int Nx = 40, Ny = 40; // число узлов КЭ по ширине и высоте

    MeshGenerator meshGen(Lx, Ly, Nx, Ny);
    meshGen.generateMesh();
    meshGen.exportToFile();
    
    std::cout << "Генерация сетки завершена. Результаты записаны в mesh_nodes.csv" << std::endl;
    return 0;
}
