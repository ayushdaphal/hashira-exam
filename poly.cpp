#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <algorithm>

// Simple JSON parser for our specific format
class JSONParser {
private:
    std::string removeSpaces(const std::string& str) {
        std::string result;
        for (char c : str) {
            if (c != ' ' && c != '\n' && c != '\t') {
                result += c;
            }
        }
        return result;
    }

public:
    struct Point {
        int x;
        int base;
        std::string value;
        long long decimal_value;
    };
    
    struct JSONData {
        int n, k;
        std::vector<Point> points;
    };

    JSONData parseJSON(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        std::string content, line;
        while (std::getline(file, line)) {
            content += line;
        }
        file.close();

        content = removeSpaces(content);
        JSONData data;
        
        // Parse n and k
        size_t nPos = content.find("\"n\":");
        size_t kPos = content.find("\"k\":");
        
        if (nPos != std::string::npos) {
            size_t start = content.find(':', nPos) + 1;
            size_t end = content.find(',', start);
            data.n = std::stoi(content.substr(start, end - start));
        }
        
        if (kPos != std::string::npos) {
            size_t start = content.find(':', kPos) + 1;
            size_t end = content.find_first_of(",}", start);
            data.k = std::stoi(content.substr(start, end - start));
        }

        // Parse points
        for (int i = 1; i <= data.n; i++) {
            std::string key = "\"" + std::to_string(i) + "\":";
            size_t keyPos = content.find(key);
            
            if (keyPos != std::string::npos) {
                Point point;
                point.x = i;
                
                // Parse base
                size_t basePos = content.find("\"base\":", keyPos);
                size_t baseStart = content.find('\"', content.find(':', basePos)) + 1;
                size_t baseEnd = content.find('\"', baseStart);
                point.base = std::stoi(content.substr(baseStart, baseEnd - baseStart));
                
                // Parse value
                size_t valuePos = content.find("\"value\":", keyPos);
                size_t valueStart = content.find('\"', content.find(':', valuePos)) + 1;
                size_t valueEnd = content.find('\"', valueStart);
                point.value = content.substr(valueStart, valueEnd - valueStart);
                
                data.points.push_back(point);
            }
        }
        
        return data;
    }
};

class PolynomialSolver {
private:
    // Convert number from given base to decimal
    long long convertBaseToDecimal(const std::string& value, int base) {
        long long decimal_value = 0;
        long long power = 1;
        
        for (int i = value.length() - 1; i >= 0; i--) {
            int digit_value;
            char digit = value[i];
            
            if (digit >= '0' && digit <= '9') {
                digit_value = digit - '0';
            } else if (digit >= 'a' && digit <= 'z') {
                digit_value = digit - 'a' + 10;
            } else if (digit >= 'A' && digit <= 'Z') {
                digit_value = digit - 'A' + 10;
            } else {
                throw std::invalid_argument("Invalid character in number: " + std::string(1, digit));
            }
            
            if (digit_value >= base) {
                throw std::invalid_argument("Digit value exceeds base");
            }
            
            decimal_value += digit_value * power;
            power *= base;
        }
        
        return decimal_value;
    }

    // Solve 2x2 system using simple substitution/elimination
    std::pair<double, double> solve2x2(double a1, double b1, double c1,
                                       double a2, double b2, double c2) {
        // System: a1*x + b1*y = c1
        //         a2*x + b2*y = c2
        
        // Use elimination: multiply first equation by a2/a1, subtract from second
        if (std::abs(a1) < 1e-10) {
            throw std::runtime_error("Cannot solve: a1 coefficient is zero");
        }
        
        double multiplier = a2 / a1;
        
        // New second equation: (a2 - a1*multiplier)*x + (b2 - b1*multiplier)*y = (c2 - c1*multiplier)
        double new_b = b2 - b1 * multiplier;
        double new_c = c2 - c1 * multiplier;
        
        if (std::abs(new_b) < 1e-10) {
            throw std::runtime_error("System has no unique solution");
        }
        
        double a0 = new_c / new_b;  // This is our constant term!
        
        // Back substitute to find a1
        double a1_coeff = (c1 - b1 * a0) / a1;
        
        return {a1_coeff, a0};
    }

    // Solve 3x3 system using substitution method
    std::vector<double> solve3x3(const std::vector<std::vector<double>>& matrix,
                                 const std::vector<double>& constants) {
        
        // Simple Gaussian elimination for 3x3
        std::vector<std::vector<double>> aug = matrix;
        for (int i = 0; i < 3; i++) {
            aug[i].push_back(constants[i]);
        }
        
        // Forward elimination
        for (int i = 0; i < 3; i++) {
            // Find pivot
            int maxRow = i;
            for (int k = i + 1; k < 3; k++) {
                if (std::abs(aug[k][i]) > std::abs(aug[maxRow][i])) {
                    maxRow = k;
                }
            }
            std::swap(aug[i], aug[maxRow]);
            
            // Make all rows below this one 0 in current column
            for (int k = i + 1; k < 3; k++) {
                if (std::abs(aug[i][i]) < 1e-10) continue;
                
                double factor = aug[k][i] / aug[i][i];
                for (int j = i; j < 4; j++) {
                    aug[k][j] -= factor * aug[i][j];
                }
            }
        }
        
        // Back substitution
        std::vector<double> solution(3);
        for (int i = 2; i >= 0; i--) {
            solution[i] = aug[i][3];
            for (int j = i + 1; j < 3; j++) {
                solution[i] -= aug[i][j] * solution[j];
            }
            solution[i] /= aug[i][i];
        }
        
        return solution;
    }

public:
    long long substitutionMethod(const std::vector<std::pair<int, long long>>& points, int degree) {
        int n = points.size();
        
        if (degree == 1) {
            // Linear: f(x) = a0 + a1*x
            int x1 = points[0].first, y1 = points[0].second;
            int x2 = points[1].first, y2 = points[1].second;
            
            // Solve: a1*x1 + a0*1 = y1
            //        a1*x2 + a0*1 = y2
            auto result = solve2x2(x1, 1, y1, x2, 1, y2);
            return static_cast<long long>(std::round(result.second));
        }
        else if (degree == 2) {
            // Quadratic: f(x) = a0 + a1*x + a2*x²
            std::vector<std::vector<double>> matrix(3, std::vector<double>(3));
            std::vector<double> constants(3);
            
            for (int i = 0; i < 3; i++) {
                int x = points[i].first;
                long long y = points[i].second;
                
                matrix[i][0] = 1;        // coefficient of a0
                matrix[i][1] = x;        // coefficient of a1
                matrix[i][2] = x * x;    // coefficient of a2
                constants[i] = y;
            }
            
            auto solution = solve3x3(matrix, constants);
            return static_cast<long long>(std::round(solution[0]));
        }
        else {
            // Higher degree - use general substitution
            std::vector<std::vector<double>> A(n, std::vector<double>(n));
            std::vector<double> b(n);
            
            for (int i = 0; i < n; i++) {
                int x = points[i].first;
                long long y = points[i].second;
                b[i] = y;
                
                for (int j = 0; j < n; j++) {
                    A[i][j] = std::pow(x, j);
                }
            }
            
            // Solve using simple Gaussian elimination
            for (int i = 0; i < n; i++) {
                // Find pivot
                int maxRow = i;
                for (int k = i + 1; k < n; k++) {
                    if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
                        maxRow = k;
                    }
                }
                
                // Swap rows
                if (maxRow != i) {
                    std::swap(A[i], A[maxRow]);
                    std::swap(b[i], b[maxRow]);
                }
                
                // Elimination
                for (int k = i + 1; k < n; k++) {
                    if (std::abs(A[i][i]) < 1e-10) continue;
                    
                    double factor = A[k][i] / A[i][i];
                    for (int j = i; j < n; j++) {
                        A[k][j] -= factor * A[i][j];
                    }
                    b[k] -= factor * b[i];
                }
            }
            
            // Back substitution
            std::vector<double> solution(n);
            for (int i = n - 1; i >= 0; i--) {
                solution[i] = b[i];
                for (int j = i + 1; j < n; j++) {
                    solution[i] -= A[i][j] * solution[j];
                }
                solution[i] /= A[i][i];
            }
            
            return static_cast<long long>(std::round(solution[0])); // a0 is the constant term
        }
    }

    long long solvePolynomialConstant(const JSONParser::JSONData& data) {
        std::cout << "Roots provided: " << data.n << ", Polynomial degree: " << (data.k - 1) << "\n";
        
        // Convert all points to decimal and display
        std::vector<std::pair<int, long long>> points;
        
        for (const auto& point : data.points) {
            long long decimal_val = convertBaseToDecimal(point.value, point.base);
            points.push_back({point.x, decimal_val});
            
            std::cout << "Root " << point.x << ": base " << point.base 
                      << " value '" << point.value << "' = " << decimal_val << "\n";
        }
        
        // Sort points by x-coordinate
        std::sort(points.begin(), points.end());
        
        // Use first k points for the system
        std::vector<std::pair<int, long long>> selected_points(points.begin(), points.begin() + data.k);
        
        long long constant_c = substitutionMethod(selected_points, data.k - 1);
        
        std::cout << "Constant coefficient c = " << constant_c << "\n";
        
        return constant_c;
    }
};

int main() {
    try {
        JSONParser parser;
        PolynomialSolver solver;
        
        std::cout << "POLYNOMIAL CONSTANT FINDER\n";
        std::cout << std::string(40, '=') << "\n";
        
        // Test case 1
        std::cout << "\nTEST CASE 1:\n";
        try {
            auto data1 = parser.parseJSON("test_case_1.json");
            long long constant1 = solver.solvePolynomialConstant(data1);
            
            std::cout << "\nTEST CASE 2:\n";
            auto data2 = parser.parseJSON("test_case_2.json");
            long long constant2 = solver.solvePolynomialConstant(data2);
            
            std::cout << "\n" << std::string(40, '=') << "\n";
            std::cout << "FINAL RESULTS:\n";
            std::cout << "Test case 1: c = " << constant1 << "\n";
            std::cout << "Test case 2: c = " << constant2 << "\n";
            
        } catch (const std::exception& e) {
            std::cout << "Error: " << e.what() << "\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}