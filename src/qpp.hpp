#pragma once

class optimization::QPP : public NLPP {
public:
    using NLPP::Method;
    using NLPP::NLPP;

    std::map<algebra::Variable, algebra::Fraction> optimize() const {
        const NLPP standard = standardize(Method::WOLFS);
        auto [multiplier, variables, equations, size] = standard.lagrange_multiplier(Method::WOLFS);
        const LPP lpp(type, {},
                      equations | std::views::take(size - restrictions.size()) |
                          std::views::transform([](const algebra::Equation& equation) -> algebra::Inequation { return equation; }) |
                          std::ranges::to<std::vector>(),
                      variables | std::views::drop(size - constraints.size()) |
                          std::views::transform([](const algebra::Variable& variable) -> algebra::Inequation { return variable >= 0; }) |
                          std::ranges::to<std::vector>());
        GLOBAL_FORMATTING << lpp;
        return {};
    }
};
