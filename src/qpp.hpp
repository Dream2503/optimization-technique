#pragma once

class optimization::QPP : public NLPP {
public:
    using NLPP::Method;
    using NLPP::NLPP;

    std::map<algebra::Variable, algebra::Fraction> optimize() const {
        const NLPP standard = standardize(Method::WOLFS);
        auto [multiplier, variables, equations, size] = standard.lagrange_multiplier(Method::WOLFS);
        std::map<algebra::Variable, algebra::Variable> complementary_slack, reversed;

        for (int i = 0; i < constraints.size(); i++) {
            complementary_slack.emplace(algebra::Variable("L" + std::to_string(i + 1)), algebra::Variable("s" + std::to_string(i + 1)));
        }
        for (int i = 0; i < restrictions.size(); i++) {
            complementary_slack.emplace(algebra::Variable("L" + std::to_string(constraints.size() + i + 1)),
                                        variables[constraints.size() + restrictions.size() + i]);
        }
        for (const auto& [key, value] : complementary_slack) {
            reversed.emplace(value, key);
        }
        complementary_slack.merge(reversed);

        const LPP lpp(
            type, objective,
            equations | std::views::take(size - restrictions.size()) | std::views::transform([](algebra::Equation equation) -> algebra::Inequation {
                for (algebra::Variable& variable : equation.lhs.terms) {
                    if (variable.variables[0].name[0] == 's') {
                        variable.variables[0].exponent = 1;
                    }
                }
                return equation;
            }) | std::ranges::to<std::vector>(),
            std::views::join(std::array{std::ranges::subrange(variables.begin(), variables.begin() + constraints.size()),
                                        std::ranges::subrange(variables.begin() + constraints.size() + restrictions.size(), variables.end())}) |
                std::views::transform([](const algebra::Variable& variable) -> algebra::Inequation { return variable >= 0; }) |
                std::ranges::to<std::vector>());
        lpp.tabular_optimize(LPP::Method::SIMPLEX, ArtificialMethod::TWO_PHASE, complementary_slack);
        return {};
    }
};
