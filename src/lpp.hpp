#pragma once

class optimization::LPP {
    static constexpr auto serial_class = detail::SerialClass::LPP;

public:
    enum class Method : bool { SIMPLEX, DUAL_SIMPLEX };
    enum class ArtificialMethod : bool { BIG_M, TWO_PHASE };
    Optimization type;
    algebra::SimplePolynomial objective;
    std::vector<algebra::Inequation> constraints, restrictions;
    inline static const algebra::Variable B{"@"}, M{"M"}, Z{"Z"}; // '@' for ascii arrangement

    static algebra::Inequation unrestrict(const algebra::Variable& variable) { return algebra::Inequation(variable, {}, algebra::inf); }

    LPP() = default;

    LPP(const Optimization type, const algebra::SimplePolynomial& objective, const std::vector<algebra::Inequation>& constraints,
        const std::vector<algebra::Inequation>& restrictions) :
        type(type), objective(objective), constraints(constraints), restrictions(restrictions) {} // need to add variables integrity check

    LPP standardize(const Method method = Method::SIMPLEX) const {
        LPP lpp = *this;
        int i = 1;
        GLOBAL_FORMATTING << lpp;

        if (lpp.type == Optimization::MINIMIZE) {
            lpp.objective *= -1;
            lpp.type = Optimization::MAXIMIZE;
        }
        for (algebra::Inequation& constraint : lpp.constraints) {
            if (method == Method::DUAL_SIMPLEX && constraint.opr == algebra::RelationalOperator::GE ||
                method == Method::SIMPLEX && static_cast<algebra::Fraction>(constraint.rhs) < 0) {
                constraint = constraint.invert();
            }
            if (constraint.opr != algebra::RelationalOperator::EQ) {
                algebra::Variable variable("s" + std::to_string(i++));
                constraint =
                    algebra::Equation(constraint.lhs + (constraint.opr == algebra::RelationalOperator::LE ? variable : -variable), constraint.rhs);
            }
        }
        GLOBAL_FORMATTING << "Standard Form:" << std::endl << lpp;
        return lpp;
    }

    LPP canonicalize() const {
        LPP lpp = *this;
        GLOBAL_FORMATTING << *this;

        for (const algebra::Inequation& constraint : constraints) {
            if (constraint.opr == algebra::RelationalOperator::EQ) {
                lpp.constraints.push_back(algebra::Inequation(constraint.lhs, algebra::RelationalOperator::LE, constraint.rhs));
                lpp.constraints.push_back(algebra::Inequation(constraint.lhs, algebra::RelationalOperator::GE, constraint.rhs));
            }
        }
        std::erase_if(lpp.constraints,
                      [](const algebra::Inequation& inequation) -> bool { return inequation.opr == algebra::RelationalOperator::EQ; });

        for (algebra::Inequation& constraint : lpp.constraints) {
            if (type == Optimization::MAXIMIZE && constraint.opr == algebra::RelationalOperator::GE ||
                type == Optimization::MINIMIZE && constraint.opr == algebra::RelationalOperator::LE) {
                constraint = constraint.invert();
            }
        }
        GLOBAL_FORMATTING << "Canonical Form:" << std::endl << lpp;
        return lpp;
    }

    std::variant<std::map<algebra::Variable, algebra::Fraction>, Solution> optimize_graphical(const std::string& path) const {
        assert(objective.terms.size() <= 2);
        const int size = constraints.size();
        algebra::Point res;
        algebra::Fraction limit, optimal = type == Optimization::MAXIMIZE ? -algebra::inf : algebra::inf, second_optimal = optimal;
        algebra::Graph graph;
        std::vector<algebra::SimplePolynomial> polynomials;
        std::vector<algebra::Point> points{{0, 0}};
        const std::vector<std::vector<int>> combinations = algebra::detail::generate_combinations(size, 2);
        polynomials.reserve(size);
        GLOBAL_FORMATTING << *this;

        for (const algebra::Inequation& constraint : constraints) {
            if (static_cast<algebra::Fraction>(constraint.rhs) != 0) {
                polynomials.push_back(constraint.lhs / static_cast<algebra::Fraction>(constraint.rhs));
            }
            for (const algebra::Variable& variable : polynomials.back().terms) {
                if (variable.variables == algebra::Variable("x").variables) {
                    points.emplace_back(variable.coefficient.reciprocate(), 0);
                } else if (variable.variables == algebra::Variable("y").variables) {
                    points.emplace_back(0, variable.coefficient.reciprocate());
                }
            }
        }
        for (const std::vector<int>& combination : combinations) {
            std::map<algebra::Variable, algebra::Fraction> solution =
                tensor::solve_linear_system({algebra::Equation(constraints[combination[0]]), algebra::Equation(constraints[combination[1]])});
            points.emplace_back(solution[algebra::Variable("x")], solution[algebra::Variable("y")]);
        }
        std::ranges::sort(points);
        points.erase(std::ranges::unique(points).begin(), points.end());
        erase_if(points, [](const algebra::Point& point) -> bool { return point.x < 0 || point.y < 0; });

        for (const algebra::Point& point : points) {
            std::map<algebra::Variable, algebra::Fraction> substituent = {{algebra::Variable("x"), point.x}, {algebra::Variable("y"), point.y}};

            if (point.x >= 0 && point.y >= 0) {
                limit = std::max({limit, point.x, point.y});
            }
            if (point.x >= 0 && point.y >= 0 &&
                std::ranges::all_of(std::array{constraints, restrictions} | std::views::join,
                                    [&substituent](const algebra::Inequation& constraint) -> bool {
                                        return static_cast<bool>(constraint.substitute(substituent, false));
                                    })) {
                const auto value = static_cast<algebra::Fraction>(objective.substitute(substituent, false));

                if (type == Optimization::MAXIMIZE && optimal < value || type == Optimization::MINIMIZE && optimal > value) {
                    second_optimal = optimal;
                    optimal = value;
                    res = point;
                } else if (type == Optimization::MAXIMIZE && second_optimal < value || type == Optimization::MINIMIZE && second_optimal > value) {
                    second_optimal = value;
                }
            }
        }
        if (graph.plot(constraints, points, limit, path) && type == Optimization::MAXIMIZE) {
            GLOBAL_FORMATTING << "Unbounded Solution" << std::endl;
            return Solution::UNBOUNDED;
        }
        if (type == Optimization::MAXIMIZE && optimal == -algebra::inf || type == Optimization::MINIMIZE && optimal == algebra::inf) {
            GLOBAL_FORMATTING << "Infeasible Solution" << std::endl;
            return Solution::INFEASIBLE;
        }
        if (second_optimal == optimal) {
            GLOBAL_FORMATTING << "Infinitely Many Solutions" << std::endl;
            return Solution::ALTERNATE;
        }
        GLOBAL_FORMATTING << (Z == optimal) << (algebra::Variable("x") == res.x) << (algebra::Variable("y") == res.y) << std::endl;
        return std::map{std::pair{Z, optimal}, {algebra::Variable("x"), res.x}, {algebra::Variable("y"), res.y}};
    }

    std::string to_latex() const {
        const int size = restrictions.size();
        std::string res("\\text{");
        res.append(type == Optimization::MAXIMIZE ? "Maximize" : "Minimize")
            .append("} \\quad & ")
            .append(objective.to_latex())
            .append("\\\\\n\\text{subject to} \\quad\n");

        for (const algebra::Inequation& constraint : constraints) {
            res.append("& ").append(constraint.to_latex()).append("\\\\\n");
        }
        res.append("& ");

        for (int i = 0; i < size; i++) {
            if (static_cast<algebra::Fraction>(restrictions[i].rhs) == algebra::inf) {
                res.append(restrictions[i].lhs.to_latex()).append("\\text{ is unrestricted}");
            } else {
                res.append(restrictions[i].to_latex());
            }
            if (i < size - 1) {
                res.append(",\\; ");
            }
        }
        return res;
    }

    void serialize(std::ofstream& out) const {
        out.write(reinterpret_cast<const char*>(&serial_class), sizeof(serial_class));
        out.write(reinterpret_cast<const char*>(&type), sizeof(type));
        objective.serialize(out);
        size_t size = constraints.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));

        for (const algebra::Inequation& inequation : constraints) {
            inequation.serialize(out);
        }
        size = restrictions.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));

        for (const algebra::Inequation& restriction : restrictions) {
            restriction.serialize(out);
        }
    }

    static LPP deserialize(std::ifstream& in) {
        detail::SerialClass type;
        in.read(reinterpret_cast<char*>(&type), sizeof(type));
        assert(type == serial_class);

        LPP res;
        size_t size;
        in.read(reinterpret_cast<char*>(&res.type), sizeof(res.type));
        res.objective = algebra::SimplePolynomial::deserialize(in);
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        res.constraints.reserve(size);

        for (int i = 0; i < size; i++) {
            res.constraints.push_back(algebra::Inequation::deserialize(in));
        }
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        res.restrictions.reserve(size);

        for (int i = 0; i < size; i++) {
            res.restrictions.push_back(algebra::Inequation::deserialize(in));
        }
        return res;
    }

    std::variant<std::vector<std::map<algebra::Variable, algebra::Fraction>>, Solution>
    tabular_optimize(Method = Method::SIMPLEX, ArtificialMethod = ArtificialMethod::BIG_M,
                     const std::map<algebra::Variable, algebra::Variable>& = {}) const;

    LPP dual(const std::string& = "w") const;
};

namespace std {
    inline string to_string(const optimization::LPP& lpp) {
        const int size = lpp.restrictions.size();
        string res(lpp.type == optimization::Optimization::MAXIMIZE ? "Maximize" : "Minimize");
        res.append("    ").append(to_string(lpp.objective)).append("\nsubject to  ").append(to_string(lpp.constraints.front())).push_back('\n');

        for (const algebra::Inequation& constraint : lpp.constraints | std::views::drop(1)) {
            res.append("            ").append(to_string(constraint)).push_back('\n');
        }
        res.append("            ");

        for (int i = 0; i < size; i++) {
            if (static_cast<algebra::Fraction>(lpp.restrictions[i].rhs) == algebra::inf) {
                res.append(to_string(lpp.restrictions[i].lhs)).append(" is unrestricted");
            } else {
                res.append(to_string(lpp.restrictions[i]));
            }
            if (i < size - 1) {
                res.append(", ");
            }
        }
        return res.append("\n");
    }
} // namespace std

inline std::ostream& optimization::operator<<(std::ostream& out, const LPP& lpp) { return out << std::to_string(lpp); }

inline std::vector<std::map<algebra::Variable, algebra::Fraction>>
optimization::basic_feasible_solutions(const std::vector<algebra::Equation>& equations) {
    std::map<algebra::Variable, std::vector<algebra::Fraction>> variables;

    for (const algebra::Equation& equation : equations) {
        GLOBAL_FORMATTING << equation << std::endl;
    }
    for (const algebra::Equation& equation : equations) {
        for (const algebra::Variable& variable : equation.lhs.terms) {
            variables[variable.basis()].push_back(variable.coefficient);
        }
        variables[LPP::B].push_back(static_cast<algebra::Fraction>(equation.rhs));
    }
    const int col = std::ranges::max(variables | std::views::drop(1) | std::views::values, {}, &std::vector<algebra::Fraction>::size).size();
    const int row = variables.size() - 1, n = std::max(row, col), k = std::min(row, col);
    std::vector<std::map<algebra::Variable, algebra::Fraction>> result;

    for (const std::vector<int>& combination : algebra::detail::generate_combinations(n, k)) {
        tensor::Matrix<algebra::Fraction> B(k, k), C(k, 1);
        std::vector<algebra::Variable> X;
        X.reserve(k);

        for (int i = 0; i < k; i++) {
            const auto itr = std::next(variables.begin(), combination[i] + 1); // B

            for (int j = 0; j < k; j++) {
                B[j, i] = itr->second[j];
            }
            X.push_back(itr->first);
            C[i, 0] = variables[LPP::B][i];
        }
        tensor::Matrix<algebra::Fraction> res = B.inverse() * C;
        std::map<algebra::Variable, algebra::Fraction> element;

        for (int i = 0; i < k; i++) {
            element[X[i]] = res[i, 0];
        }
        result.push_back(element);
    }
    GLOBAL_FORMATTING << "Basic Feasible Solutions:" << std::endl;

    for (const std::map<algebra::Variable, algebra::Fraction>& res : result) {
        for (const auto& [variable, fraction] : res) {
            GLOBAL_FORMATTING << algebra::Equation(variable, fraction) << "  ";
        }
        GLOBAL_FORMATTING << std::endl;
    }
    return result;
}
