#pragma once

class optimization::NLPP : public LPP {
public:
    enum class Method : uint8_t { LAGRANGE, KKT, WOLFS };
    using LPP::LPP;

    std::tuple<algebra::SimplePolynomial, std::vector<algebra::Variable>, std::vector<algebra::Equation>, int>
    lagrange_multiplier(const Method method) const {
        const int size = constraints.size();
        algebra::SimplePolynomial multiplier = objective;
        std::vector<algebra::Variable> variables;
        std::vector<algebra::Equation> equations;
        std::set<algebra::Variable> seen;
        variables.reserve(size);

        for (int i = 0; i < size; i++) {
            variables.emplace_back("L" + std::to_string(i + 1));
            multiplier -= variables.back() * (constraints[i].lhs - constraints[i].rhs);
        }
        GLOBAL_FORMATTING << (algebra::Variable("L") == multiplier) << std::endl;

        for (const algebra::Variable& variable :
             std::array{objective.terms, variables,
                        method == Method::WOLFS
                            ? constraints | std::views::transform([](const algebra::Inequation& constraint) -> std::vector<algebra::Variable> {
                                  return constraint.lhs.terms;
                              }) |
                                std::views::join | std::ranges::to<std::vector>()
                            : std::vector<algebra::Variable>()} |
                 std::views::join | std::views::filter([&seen](const algebra::Variable& var) -> bool {
                     return !var.is_fraction() && !seen.contains(var.variables[0].name);
                 })) {
            equations.push_back(multiplier.differentiate(variable.variables[0].name) == 0);
            const auto itr = std::ranges::find_if(equations.back().lhs.terms, [](const algebra::Variable& var) -> bool { return var.is_fraction(); });

            if (itr != equations.back().lhs.terms.end()) {
                equations.back() -= *itr;
            }
            seen.emplace(variable.variables[0].name);
        }
        for (const algebra::Variable& variable : variables) {
            seen.erase(variable);
        }
        variables.insert(variables.begin(), seen.begin(), seen.end());
        return {multiplier, variables, equations, seen.size()};
    }

    NLPP standardize(const Method method) const {
        NLPP nlpp = *this;
        int i = 1;
        GLOBAL_FORMATTING << nlpp << std::endl;
        nlpp.constraints.clear();

        for (const algebra::Inequation& constraint :
             std::array{constraints, method == Method::WOLFS ? restrictions : std::vector<algebra::Inequation>()} | std::views::join) {
            if (static_cast<algebra::Fraction>(constraint.rhs) < 0) {
                nlpp.constraints.push_back(constraint.invert());
            } else {
                nlpp.constraints.push_back(constraint);
            }
            if (nlpp.constraints.back().opr != algebra::RelationalOperator::EQ) {
                algebra::Variable variable("s" + std::to_string(i++));
                nlpp.constraints.back() = type == Optimization::MAXIMIZE && constraint.opr == algebra::RelationalOperator::GE ||
                        type == Optimization::MINIMIZE && constraint.opr == algebra::RelationalOperator::LE
                    ? nlpp.constraints.back().invert()
                    : nlpp.constraints.back();
                nlpp.constraints.back().lhs += variable ^ 2;
                nlpp.constraints.back().opr = algebra::RelationalOperator::EQ;
            }
        }
        GLOBAL_FORMATTING << "Standard Form:" << std::endl << nlpp << std::endl;
        return nlpp;
    }

    std::map<algebra::Variable, algebra::Fraction> optimize() const {
        if (std::ranges::all_of(constraints,
                                [](const algebra::Inequation& constraint) -> bool { return constraint.opr == algebra::RelationalOperator::EQ; })) {
            return lagrange_multiplier_method();
        }
        return kkt_conditions();
    }

    std::map<algebra::Variable, algebra::Fraction> lagrange_multiplier_method() const {
        GLOBAL_FORMATTING << *this << std::endl;
        NLPP temp_nlpp = *this;

        for (algebra::Inequation& constraint : temp_nlpp.constraints) {
            constraint -= constraint.rhs;
        }
        const auto [multiplier, variables, equations, size] = temp_nlpp.lagrange_multiplier(Method::LAGRANGE);
        std::map<algebra::Variable, algebra::Fraction> res = tensor::solve_linear_system(equations);
        tensor::Matrix<algebra::Fraction> sufficient(size, size);
        std::vector<uint8_t> signs;
        signs.reserve(size - 2);

        if (size > 2) {
            GLOBAL_FORMATTING << "Sufficiency Test:" << std::endl;
        }

        for (int i = 1; i < size; i++) {
            sufficient[i, 0] = sufficient[0, i] = static_cast<algebra::Fraction>(constraints[0].lhs.differentiate(variables[i - 1], false));
        }
        for (int i = 1; i < size; i++) {
            for (int j = 1; j <= i; j++) {
                sufficient[i, j] = sufficient[j, i] =
                    static_cast<algebra::Fraction>(multiplier.differentiate(variables[i - 1], false).differentiate(variables[j - 1], false));
            }
        }
        for (int i = 3; i <= size; i++) {
            sufficient.row = sufficient.column = i;
            auto x = sufficient.determinant();
            signs.push_back(x > 0);
        }
        if (std::ranges::all_of(signs, [](const bool sign) -> bool { return !sign; }) && type == Optimization::MINIMIZE ||
            [&signs] {
                const int sign_size = signs.size();

                if (sign_size <= 1) {
                    return true;
                }
                const bool even = signs[0], odd = signs[1];

                for (int i = 2; i < sign_size; i++) {
                    if (i % 2 == 0 && even != signs[i] || i % 2 == 1 && odd != signs[i]) {
                        return false;
                    }
                }
                return true;
            }() &&
                type == Optimization::MAXIMIZE) {
            std::erase_if(res, [&variables, size](const std::pair<algebra::Variable, algebra::Fraction>& element) -> bool {
                return std::ranges::contains(variables.begin() + size, variables.end(), element.first);
            });
            res[Z] = static_cast<algebra::Fraction>(objective.substitute(res));
            return res;
        }
        return {};
    }

    std::map<algebra::Variable, algebra::Fraction> kkt_conditions() const {
        const NLPP standard = standardize(Method::KKT);
        auto [multiplier, variables, equations, size] = standard.lagrange_multiplier(Method::KKT);
        const int cases = constraints.size();
        std::map<algebra::Variable, algebra::Fraction> substituent;
        std::vector<algebra::Inequation> lagrange_constraints;
        equations.resize(size);
        lagrange_constraints.reserve(variables.size() - size);

        for (int i = size; i < variables.size(); i++) {
            lagrange_constraints.push_back(variables[i] >= 0);
        }
        for (int i = 0; i < 1 << cases; i++) {
            std::vector<algebra::Equation> temp_equations;
            temp_equations.reserve(size);
            substituent.clear();
            GLOBAL_FORMATTING << std::endl << "Case " << i + 1 << ": ";

            for (int j = 0; j < cases; j++) {
                if (i >> j & 1) {
                    temp_equations.push_back(constraints[j]);
                    GLOBAL_FORMATTING << algebra::Inequation(variables[size + j], algebra::RelationalOperator::NE, 0) << ' ';
                } else {
                    substituent.emplace(variables[size + j], 0);
                    GLOBAL_FORMATTING << algebra::Equation(variables[size + j], 0) << ' ';
                }
            }
            GLOBAL_FORMATTING << std::endl;

            for (const algebra::Equation& equation : equations) {
                temp_equations.push_back(equation.substitute(substituent));
            }
            GLOBAL_FORMATTING << std::endl;
            std::map<algebra::Variable, algebra::Fraction> res = tensor::solve_linear_system(temp_equations);

            if (res.empty()) {
                continue;
            }
            GLOBAL_FORMATTING << std::endl;
            res.insert(substituent.begin(), substituent.end());

            for (const algebra::Variable& variable : variables) {
                res.emplace(variable, 0);
            }
            if (std::ranges::all_of(
                    std::array{constraints, restrictions, lagrange_constraints} | std::views::join,
                    [&res](const algebra::Inequation& constraint) -> bool { return static_cast<bool>(constraint.substitute(res)); })) {
                for (int j = size; j < variables.size(); j++) {
                    res.erase(variables[j]);
                }
                res[Z] = static_cast<algebra::Fraction>(objective.substitute(res));
                return res;
            }
        }
        return {};
    }
};
