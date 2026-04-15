#include "optimization.hpp"

using namespace algebra;
using namespace tensor;
using namespace optimization;

void test(LPP&& lpp, const std::string& method = "simplex", const Variable& var = {}, const Matrix<Fraction>& coefficients = {}) {
    if (method == "simplex" || method == "two phase") {
        auto x = lpp.tabular_optimize(LPP::Method::SIMPLEX, method == "two phase" ? LPP::ArtificialMethod::TWO_PHASE : LPP::ArtificialMethod::BIG_M);
    } else if (method == "dual") {
        auto x = lpp.tabular_optimize(LPP::Method::DUAL_SIMPLEX);
    } else if (method.starts_with("Var")) {
        lpp = lpp.standardize();
        ComputationalTable table(lpp);
        table.optimize();

        if (method.ends_with("add") || method.ends_with("remove")) {
            if (method.ends_with("add")) {
                table.add_variable(var, coefficients);
            } else {
                table.remove_variable(var);
            }
            table.get_solutions(LPP::Method::SIMPLEX);
        } else {
            method.ends_with('C') ? table.cost_variation() : table.RHS_variation();
        }
    } else if (method == "graphical") {
        auto x = lpp.optimize_graphical(var.variables.front().name);
    } else {
        auto x = lpp.dual(method);
    }
}

void test(std::vector<Equation>&& equations) { basic_feasible_solutions(equations); }

void test(ComputationalTable&& table, const std::string& method = "simplex", const Variable& var = {}, const Inequation& constraint = {},
          const Matrix<Fraction>& coefficients = {}) {
    optimization::GLOBAL_FORMATTING << table << std::endl;

    if (method.starts_with("Var")) {
        if (method.ends_with("add") || method.ends_with("remove")) {
            if (method.ends_with("add")) {
                table.add_variable(var, coefficients);
            } else {
                table.remove_variable(var);
            }
            table.get_solutions(LPP::Method::SIMPLEX);
        } else {
            method.ends_with('C') ? table.cost_variation() : table.RHS_variation();
        }
    } else if (method.starts_with("Constraint")) {
        if (method.ends_with("add")) {
            table.add_constraint(constraint);
        }
        table.get_solutions(LPP::Method::SIMPLEX);
    }
}

void test(IPP&& ipp, const std::string& path) { auto x = ipp.optimize_branch_bound(path); }

void test(NLPP&& nlpp) {
    auto res = nlpp.optimize();

    for (const auto& [variable, fraction] : res) {
        optimization::GLOBAL_FORMATTING << Equation(variable, fraction) << "  ";
    }
    optimization::GLOBAL_FORMATTING << std::endl;
}

void test(QPP&& qpp) {
    auto res = qpp.optimize();

    for (const auto& [variable, fraction] : res) {
        optimization::GLOBAL_FORMATTING << Equation(variable, fraction) << "  ";
    }
    optimization::GLOBAL_FORMATTING << std::endl;
}

int main() {
    const Variable x("x"), y("y"), z("z"), x1("x1"), x2("x2"), x3("x3"), x4("x4"), x5("x5"), s1("s1"), s2("s2"), s3("s3");
    // optimization::GLOBAL_FORMATTING.toggle_file("output.txt");
    // optimization::GLOBAL_FORMATTING.toggle_latex("latex.tex");
    test(LPP(Optimization::MAXIMIZE, 2 * x + 7 * y,
             {
                 3 * x + 5 * y <= 15,
                 7 * x + 3 * y <= 21,
             },
             {x >= 0, y >= 0}),
         "graphical", Variable("outputs/graph1.png"));
    test(LPP(Optimization::MAXIMIZE, 3 * x + 5 * y,
             {
                 4 * x + 3 * y <= 12,
                 5 * x + 4 * y >= 20,
             },
             {x >= 0, y >= 0}),
         "graphical", Variable("outputs/graph2.png"));
    test(LPP(Optimization::MAXIMIZE, x + 2 * y,
             {
                 3 * x + 2 * y <= 6,
                 2 * x + 5 * y >= 10,
             },
             {x >= 0, y >= 0}),
         "graphical", Variable("outputs/graph3.png"));
    test(LPP(Optimization::MINIMIZE, 3 * x - 10 * y,
             {
                 5 * x + 2 * y >= 10,
                 4 * x + 3 * y <= 12,
             },
             {x >= 0, y >= 0}),
         "graphical", Variable("outputs/graph4.png"));
    test(LPP(Optimization::MAXIMIZE, 5 * x + 4 * y,
             {
                 2 * x + 5 * y >= 10,
                 3 * x + 4 * y >= 12,
             },
             {x >= 0, y >= 0}),
         "graphical", Variable("outputs/graph5.png"));
    test(LPP(Optimization::MINIMIZE, 5 * x + 4 * y,
             {
                 2 * x + 5 * y >= 10,
                 3 * x + 4 * y >= 12,
             },
             {x >= 0, y >= 0}),
         "graphical", Variable("outputs/graph6.png"));
    test(LPP(Optimization::MAXIMIZE, 10 * x + 4 * y,
             {
                 5 * x + 2 * y <= 100,
                 3 * x + 2 * y <= 90,
                 x + 2 * y <= 50,
             },
             {x >= 0, y >= 0}),
         "graphical", Variable("outputs/graph7.png"));
    test(LPP(Optimization::MAXIMIZE, 3 * x + 2 * y,
             {
                 x + y <= 4,
                 x - y <= 2,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, x1 + x2 + 3 * x3,
             {
                 3 * x1 + 2 * x2 + x3 <= 3,
                 2 * x1 + x2 + 2 * x3 <= 2,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(LPP(Optimization::MINIMIZE, x1 - 3 * x2 + 2 * x3,
             {
                 3 * x1 - x2 + 2 * x3 <= 7,
                 -2 * x1 + 4 * x2 <= 12,
                 -4 * x1 + 3 * x2 + 8 * x3 <= 10,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(LPP(Optimization::MAXIMIZE, 2 * x + y,
             {
                 4 * x + 3 * y <= 12,
                 4 * x + y <= 8,
                 4 * x - y <= 8,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, 2 * x + y,
             {
                 x - y <= 10,
                 2 * x - y <= 40,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, 3 * x + 2 * y,
             {
                 x - y <= 1,
                 3 * x - 2 * y <= 6,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, x1 + 2 * x2 + x3,
             {
                 2 * x1 + x2 - x3 >= -2,
                 -2 * x1 + x2 - 5 * x3 <= 6,
                 4 * x1 + x2 + x3 <= 6,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    // Big M Method
    test(LPP(Optimization::MAXIMIZE, -4 * x1 - x2,
             {
                 3 * x1 + x2 == 3,
                 4 * x1 + 3 * x2 >= 6,
                 x1 + 2 * x2 <= 3,
             },
             {x1 >= 0, x2 >= 0}));
    test(LPP(Optimization::MAXIMIZE, -x - y,
             {
                 3 * x + 2 * y >= 30,
                 -2 * x + 3 * y <= -30,
                 x + y <= 5,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, -4 * x - y,
             {
                 3 * x + y == 3,
                 4 * x + 3 * y >= 6,
                 x + 2 * y <= 3,
             },
             {x >= 0, y >= 0}),
         "two phase");
    test(LPP(Optimization::MAXIMIZE, 3 * x1 + 2 * x2 + x3,
             {
                 -3 * x1 + 2 * x2 + 2 * x3 == 8,
                 -3 * x1 + 4 * x2 + x3 == 7,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    // BFS
    test({
        x1 + 2 * x2 + x3 == 4,
        2 * x1 + x2 + 5 * x3 == 5,
    });
    test({
        2 * x1 + x2 - x3 == 2,
        3 * x1 + 2 * x2 + x3 == 3,
    });
    test(LPP(Optimization::MAXIMIZE, 2 * x1 + 3 * x2 + 10 * x3,
             {
                 x1 + 2 * x3 == 0,
                 x2 + x3 == 1,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(LPP(Optimization::MAXIMIZE, 2 * x1 + 3 * x2 + 10 * x3,
             {
                 x1 - 2 * x3 == 0,
                 x2 + x3 == 1,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(LPP(Optimization::MAXIMIZE, 2 * x + y,
             {
                 4 * x + 3 * y <= 12,
                 4 * x + y <= 8,
                 4 * x - y <= 8,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, 5 * x + 3 * y,
             {
                 x + y <= 2,
                 5 * x + 2 * y <= 10,
                 3 * x + 8 * y <= 12,
             },
             {x >= 0, y >= 0}));
    // Surprise Test
    test(LPP(Optimization::MAXIMIZE, 2 * x + 3 * y,
             {
                 x + 2 * y >= 2,
                 3 * x + y >= 3,
                 4 * x + 3 * y <= 6,
             },
             {x >= 0, y >= 0}));
    test({
        x + 2 * y - z == 3,
        3 * x - y + 2 * z == 4,
        2 * x + 3 * y - 5 * z == 7,
    });
    test(LPP(Optimization::MAXIMIZE, x1 + 2 * x2 + x3,
             {
                 2 * x1 + x2 - x3 >= -2,
                 -2 * x1 + x2 - 5 * x3 <= 6,
                 4 * x1 + x2 + x3 <= 6,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(LPP(Optimization::MAXIMIZE, 2 * x + 3 * y,
             {
                 x + 2 * y <= 4,
                 x + y == 3,
             },
             {x >= 0, y >= 0}));
    // Dual
    test(LPP(Optimization::MAXIMIZE, 3 * x + 2 * y,
             {
                 4 * x - 3 * y <= 10,
                 x + 2 * y <= 5,
                 3 * x - 5 * y >= 15,
             },
             {x >= 0, y >= 0}),
         "w");
    test(LPP(Optimization::MINIMIZE, 5 * x1 + 4 * x2 - 3 * x3,
             {
                 x1 + x2 + x3 >= 5,
                 2 * x1 + 3 * x2 - 5 * x3 <= 4,
                 x1 + 2 * x2 - 3 * x3 <= 6,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}),
         "w");
    test(LPP(Optimization::MAXIMIZE, x + 2 * y,
             {
                 3 * x + 2 * y <= 4,
                 2 * x - 5 * y == 10,
             },
             {x >= 0, y >= 0}),
         "w");
    test(LPP(Optimization::MINIMIZE, 3 * x1 + 4 * x2 + 7 * x3,
             {
                 x1 + 2 * x2 - 5 * x3 >= 4,
                 2 * x1 + 5 * x2 - x3 == 7,
                 2 * x1 + 3 * x2 - x3 <= -8,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}),
         "w");
    test(LPP(Optimization::MAXIMIZE, 3 * x + 2 * y,
             {
                 x + y <= 5,
                 2 * x + 3 * y >= 4,
                 x - y <= 2,
             },
             {LPP::unrestrict(x), y >= 0}),
         "w");
    test(LPP(Optimization::MAXIMIZE, x + 2 * y,
             {
                 2 * x + 3 * y <= 4,
                 3 * x + 4 * y == 5,
             },
             {x >= 0, y >= 0}),
         "w");
    // Dual Simplex
    test(LPP(Optimization::MAXIMIZE, -5 * x - 6 * y,
             {
                 x + y >= 2,
                 4 * x + y >= 4,
             },
             {x >= 0, y >= 0}),
         "dual");
    test(LPP(Optimization::MINIMIZE, 10 * x1 + 6 * x2 + 2 * x3,
             {
                 -x1 + x2 + x3 >= 1,
                 3 * x1 + x2 - x3 >= 2,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}),
         "dual");
    // Alternate Optimal Solution
    test(LPP(Optimization::MAXIMIZE, 2 * x + 4 * y,
             {
                 x + 2 * y <= 5,
                 x + y <= 4,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, -20 * x - 30 * y,
             {
                 2 * x + 3 * y >= 120,
                 x + y >= 40,
                 2 * x + 3 * y / 2 >= 90,
             },
             {x >= 0, y >= 0}));
    // Variation in C
    test(LPP(Optimization::MAXIMIZE, 3 * x + 5 * y,
             {
                 x + y <= 1,
                 2 * x + 3 * y <= 1,
             },
             {x >= 0, y >= 0}),
         "Var C");
    test(LPP(Optimization::MAXIMIZE, 15 * x + 45 * y,
             {
                 y <= 50,
                 x + 16 * y <= 240,
                 5 * x + 2 * y <= 162,
             },
             {x >= 0, y >= 0}),
         "Var C");
    test(LPP(Optimization::MAXIMIZE, x1 + x2 + 3 * x3,
             {
                 3 * x1 + 2 * x2 + x3 <= 3,
                 2 * x1 + x2 + 2 * x3 <= 2,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}),
         "Var C");
    // Variation in B
    test(LPP(Optimization::MAXIMIZE, -x1 + 2 * x2 - x3,
             {
                 3 * x1 + x2 - x3 <= 10,
                 -x1 + 4 * x2 + x3 >= 6,
                 x2 + x3 <= 4,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}),
         "Var B");
    test(LPP(Optimization::MAXIMIZE, 2 * x + y,
             {
                 3 * x + 5 * y <= 15,
                 6 * x + 2 * y <= 24,
             },
             {x >= 0, y >= 0}),
         "Var B");
    test(LPP(Optimization::MAXIMIZE, 15 * x + 10 * y,
             {
                 4 * x + 6 * y <= 360,
                 3 * x <= 180,
                 5 * y <= 200,
             },
             {x >= 0, y >= 0}),
         "Var B");
    test(ComputationalTable(
             {
                 {x, Variable(2)},
                 {y, Variable(3)},
                 {z, Variable(1)},
                 {s1, Variable()},
                 {s2, Variable()},
             },
             {x, y},
             {
                 {LPP::B, {1, 2}},
                 {x, {1, 0}},
                 {y, {0, 1}},
                 {z, {-1, 2}},
                 {s1, {4, -1}},
                 {s2, {-1, 1}},
             },
             Solution::OPTIMIZED, LPP{Optimization::MAXIMIZE, {}, {x <= 3, x <= 7}, {}}),
         "Var B");
    // Addition of new variable
    test(LPP(Optimization::MAXIMIZE, 3 * x + 5 * y,
             {
                 x <= 4,
                 3 * x + 2 * y <= 18,
             },
             {x >= 0, y >= 0}),
         "Var add", 7 * x1, {1, 2});
    // Deletion of a variable
    test(LPP(Optimization::MAXIMIZE, 3 * x + 5 * y,
             {
                 x <= 4,
                 3 * x + 2 * y <= 18,
             },
             {x >= 0, y >= 0}),
         "Var remove", x);
    test(LPP(Optimization::MAXIMIZE, 3 * x + 5 * y,
             {
                 x <= 4,
                 3 * x + 2 * y <= 18,
             },
             {x >= 0, y >= 0}),
         "Var remove", y);
    test(ComputationalTable(
             {
                 {x1, Variable(2)},
                 {x2, Variable(4)},
                 {x3, Variable(1)},
                 {x4, Variable(3)},
                 {x5, Variable(2)},
                 {s1, Variable()},
                 {s2, Variable()},
                 {s3, Variable()},
             },
             {x1, x2, x3},
             {
                 {LPP::B, {3, 1, 7}},
                 {x1, {1, 0, 0}},
                 {x2, {0, 1, 0}},
                 {x3, {0, 0, 1}},
                 {x4, {-1, 2, -1}},
                 {x5, {0, 1, -2}},
                 {s1, {Fraction(1, 2), -1, 5}},
                 {s2, {Fraction(1, 5), 0, Fraction(-3, 10)}},
                 {s3, {-1, Fraction(1, 2), 2}},
             },
             Solution::OPTIMIZED),
         "Var remove", x2);
    test(ComputationalTable(
             {
                 {x1, Variable(2)},
                 {x2, Variable(4)},
                 {x3, Variable(1)},
                 {x4, Variable(3)},
                 {x5, Variable(2)},
                 {s1, Variable()},
                 {s2, Variable()},
                 {s3, Variable()},
             },
             {x1, x2, x3},
             {
                 {LPP::B, {3, 1, 7}},
                 {x1, {1, 0, 0}},
                 {x2, {0, 1, 0}},
                 {x3, {0, 0, 1}},
                 {x4, {-1, 2, -1}},
                 {x5, {0, 1, -2}},
                 {s1, {Fraction(1, 2), -1, 5}},
                 {s2, {Fraction(-1, 5), 0, Fraction(-2, 5)}},
                 {s3, {-1, Fraction(1, 2), 2}},
             },
             Solution::OPTIMIZED, {}),
         "Constraint add", {}, 2 * x1 + 3 * x2 - x3 + 2 * x4 + 4 * x5 <= 5);
    test(ComputationalTable(
             {
                 {x1, Variable(2)},
                 {x2, Variable(4)},
                 {x3, Variable(1)},
                 {x4, Variable(3)},
                 {x5, Variable(2)},
                 {s1, Variable()},
                 {s2, Variable()},
                 {s3, Variable()},
             },
             {x1, x2, x3},
             {
                 {LPP::B, {3, 1, 7}},
                 {x1, {1, 0, 0}},
                 {x2, {0, 1, 0}},
                 {x3, {0, 0, 1}},
                 {x4, {-1, 2, -1}},
                 {x5, {0, 1, -2}},
                 {s1, {Fraction(1, 2), -1, 5}},
                 {s2, {Fraction(-1, 5), 0, Fraction(-2, 5)}},
                 {s3, {-1, Fraction(1, 2), 2}},
             },
             Solution::OPTIMIZED, {}),
         "Constraint add", {}, 3 * x1 + x2 + 2 * x3 + x4 + 9 * x5 <= 19);
    // Mid Term Examination
    test(LPP(Optimization::MAXIMIZE, 3 * x - 5 * y,
             {
                 4 * x + 3 * y >= 5,
                 2 * x - 5 * y <= 3,
             },
             {x >= 0, y >= 0}),
         "w");
    test(LPP(Optimization::MAXIMIZE, 2 * x + y,
             {
                 x - y <= 10,
                 2 * x - y <= 40,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, 5 * x - y,
             {
                 3 * x + 5 * y <= 15,
                 4 * x + 3 * y <= 12,
             },
             {x >= 0, y >= 0}),
         "graphical", Variable("outputs/graph8.png"));
    test(LPP(Optimization::MAXIMIZE, 5 * x - y,
             {
                 3 * x + 5 * y <= 15,
                 4 * x + 3 * y <= 12,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, x1 + x2,
             {
                 x1 + x2 <= 8,
                 2 * x1 + x2 <= 10,
             },
             {x1 >= 0, x2 >= 0}),
         "Var B");
    test(LPP(Optimization::MAXIMIZE, -4 * x - y,
             {
                 3 * x + y == 3,
                 4 * x + 3 * y >= 6,
                 x + 2 * y <= 3,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MINIMIZE, x1 - 3 * x2 + 2 * x3,
             {
                 3 * x1 - x2 + 2 * x3 <= 7,
                 -2 * x1 + 4 * x2 <= 12,
                 -4 * x1 + 3 * x2 + 8 * x3 <= 10,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(LPP(Optimization::MINIMIZE, 5 * x + 6 * y,
             {
                 x + y >= 2,
                 4 * x + y >= 4,
             },
             {x >= 0, y >= 0}),
         "dual");
    test(LPP(Optimization::MAXIMIZE, 3 * x1 + 5 * x2,
             {
                 x1 + x2 <= 1,
                 2 * x1 + 3 * x2 <= 1,
             },
             {x >= 0, y >= 0}),
         "Var C");
    // IPP
    test(IPP(Optimization::MAXIMIZE, x + 4 * y,
             {
                 2 * x + 4 * y <= 7,
                 5 * x + 3 * y <= 15,
             },
             {x >= 0, y >= 0}),
         "outputs/ipp1");
    test(IPP(Optimization::MAXIMIZE, 7 * x + 9 * y,
             {
                 -x + 3 * y <= 6,
                 7 * x + y <= 35,
                 y <= 7,
             },
             {x >= 0, y >= 0}),
         "outputs/ipp2");
    test(IPP(Optimization::MAXIMIZE, 2 * x + 3 * y,
             {
                 6 * x + 5 * y <= 12,
                 4 * x + 3 * y <= 14,
             },
             {x >= 0, y >= 0}),
         "outputs/ipp3");
    test(IPP(Optimization::MAXIMIZE, 3 * x + 4 * y,
             {
                 3 * x - y <= 12,
                 3 * x + 11 * y <= 66,
             },
             {x >= 0, y >= 0}),
         "outputs/ipp4");
    test(IPP(Optimization::MINIMIZE, 2 * x + 3 * y,
             {
                 2 * x + 3 * y <= 7,
                 x <= 2,
                 y <= 2,
             },
             {x >= 0, y >= 0}),
         "outputs/ipp5");
    test(NLPP(Optimization::MAXIMIZE, 6 * (x1 ^ 2) + 5 * (x2 ^ 2),
              {
                  x1 + 5 * x2 == 3,
              },
              {x1 >= 0, x2 >= 0}));
    test(NLPP(Optimization::MAXIMIZE, 2 * (x1 ^ 2) + 2 * (x2 ^ 2) - 24 * x1 - 8 * x2 + 2 * (x3 ^ 2) - 12 * x3 + 200,
              {
                  x1 + x2 + x3 == 11,
              },
              {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(NLPP(Optimization::MINIMIZE, (x1 ^ 2) + (x2 ^ 2) + (x3 ^ 2),
              {
                  x1 + x2 + 3 * x3 == 2,
                  5 * x1 + 2 * x2 + x3 == 5,
              },
              {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(NLPP(Optimization::MAXIMIZE, 6 * x1 + 8 * x2 - (x1 ^ 2) - (x2 ^ 2),
              {
                  4 * x1 + 3 * x2 == 16,
                  3 * x1 + 5 * x2 == 15,
              },
              {x1 >= 0, x2 >= 0}));
    test(NLPP(Optimization::MINIMIZE, 2 * (x1 ^ 2) + 2 * (x2 ^ 2) - 24 * x1 - 8 * x2 + 2 * (x3 ^ 2) - 12 * x3 + 200,
              {
                  x1 + x2 + x3 == 11,
              },
              {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(NLPP(Optimization::MINIMIZE, 2 * (x1 ^ 2) + 12 * x1 * x2 - 7 * (x2 ^ 2),
              {
                  2 * x1 + 5 * x2 <= 98,
              },
              {x1 >= 0, x2 >= 0}));
    test(NLPP(Optimization::MAXIMIZE, 8 * x1 + 10 * x2 - (x1 ^ 2) - (x2 ^ 2),
              {
                  3 * x1 + 2 * x2 <= 6,
              },
              {x1 >= 0, x2 >= 0}));
    test(NLPP(Optimization::MAXIMIZE, -(x1 ^ 2) - (x2 ^ 2) - (x3 ^ 2) + 4 * x1 + 6 * x2,
              {
                  x1 + x2 <= 2,
                  2 * x1 + 3 * x2 <= 12,
              },
              {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(NLPP(Optimization::MAXIMIZE, 2 * x1 + 3 * x2 - 2 * (x1 ^ 2),
              {
                  x1 + 4 * x2 <= 4,
                  x1 + x2 <= 2,
              },
              {x1 >= 0, x2 >= 0}));
    test(QPP(Optimization::MAXIMIZE, 2 * x1 + 3 * x2 - 2 * (x1 ^ 2),
             {
                 x1 + 4 * x2 <= 4,
                 x1 + x2 <= 2,
             },
             {x1 >= 0, x2 >= 0}));
    test(QPP(Optimization::MAXIMIZE, +2 * x1 + x2 - (x1 ^ 2),
             {
                 2 * x1 + 3 * x2 <= 6,
                 2 * x1 + x2 <= 4,
             },
             {x1 >= 0, x2 >= 0}));
    return 0;
}
