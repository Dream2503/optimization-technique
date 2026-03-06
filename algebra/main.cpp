#include "algebra.hpp"

using namespace algebra;

int main() {
    Matrix<Fraction> matrix({{1, 2}, {2, 1}});
    const Matrix<Variable> variables({{Variable("X"), Variable("Y")}, {Variable("A"), Variable("B")}});

    // matrix *= matrix;
    // std::cout << matrix << std::endl;
    //
    // const auto res = matrix * variables;
    // std::cout << res << std::endl;

    matrix /= 10;
    std::cout << matrix << std::endl;

    const auto res1 = variables / Fraction(10);
    std::cout << res1 << std::endl;

    std::cout << matrix.inverse();

    const Graph graph;
    Variable v = (Variable("x") ^ 1) * (Variable("y") ^ -1);
    std::cout << v << std::endl;
    Variable x("x");
    std::cout << (x + 1) * (x + 2) << std::endl;
    graph.plot({
        x,
        x + 2,
        2 * x - 5,
    });

    Variable y("y");
    Inequation inequation = 5 * x + 2 * y - 7 <= 5;
    Inequation res = inequation.solve_for(y);
    std::cout << res;
    return 0;
}
