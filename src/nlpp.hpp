#pragma once

class optimization::NLPP : public LPP {
public:
    NLPP(const Optimization type, const algebra::SimplePolynomial& objective, const std::vector<algebra::Inequation>& constraints,
         const std::vector<algebra::Inequation>& restrictions) : LPP(type, objective, constraints, restrictions) {
        for (algebra::Inequation& constraint : this->constraints) {
            constraint -= constraint.rhs;
        }
    }

    std::map<algebra::Variable, algebra::Fraction> optimize() const {
        GLOBAL_FORMATTING << *this;
        int size = constraints.size();
        algebra::SimplePolynomial lagrange_multiplier = objective;
        std::vector<algebra::Variable> lagrange_variables;
        std::vector<algebra::Equation> equations;
        std::set<algebra::Variable> seen;
        lagrange_variables.reserve(size);

        for (int i = 0; i < size; i++) {
            lagrange_variables.emplace_back("L" + std::to_string(i + 1));
            lagrange_multiplier -= lagrange_variables.back() * constraints[i].lhs;
        }
        GLOBAL_FORMATTING << lagrange_multiplier;

        for (const algebra::Variable& variable :
             std::array{objective.terms, lagrange_variables} | std::views::join | std::views::filter([&seen](const algebra::Variable& var) -> bool {
                 return !var.is_fraction() && !seen.contains(var.variables[0].name);
             })) {
            equations.push_back(lagrange_multiplier.differentiate(variable.variables[0].name) == 0);
            const auto itr = std::ranges::find_if(equations.back().lhs.terms, [](const algebra::Variable& var) -> bool { return var.is_fraction(); });

            if (itr != equations.back().lhs.terms.end()) {
                equations.back() -= *itr;
            }
            seen.emplace(variable.variables[0].name);
        }
        std::map<algebra::Variable, algebra::Fraction> res = tensor::solve_linear_system(equations);

        auto range = seen | std::views::filter([&lagrange_variables](const algebra::Variable& variable) -> bool {
                         return !std::ranges::contains(lagrange_variables, variable);
                     });
        const std::vector variables(range.begin(), range.end());
        size = variables.size() + 1;
        tensor::Matrix<algebra::SimplePolynomial> sufficient(size, size);

        for (int i = 1; i < size; i++) {
            sufficient[i, 0] = sufficient[0, i] = constraints[0].lhs.differentiate(variables[i - 1]);
        }
        for (int i = 1; i < size; i++) {
            for (int j = 1; j <= i; j++) {
                sufficient[i, j] = sufficient[j, i] = lagrange_multiplier.differentiate(variables[i - 1]).differentiate(variables[j - 1]);
            }
        }
        std::erase_if(res, [&lagrange_variables](const std::pair<algebra::Variable, algebra::Fraction>& element) -> bool {
            return std::ranges::contains(lagrange_variables, element.first);
        });
        res[Z] = static_cast<algebra::Fraction>(objective.substitute(res));
        GLOBAL_FORMATTING << sufficient << std::endl << algebra::Equation(Z, res[Z]) << std::endl << std::endl;
        return res;
    }
};
