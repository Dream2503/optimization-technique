#pragma once

class optimization::IPP : public LPP {
public:
    using LPP::LPP;

    std::variant<std::map<algebra::Variable, algebra::Fraction>, Solution> optimize_branch_bound(const std::string& path) const {
        int i = 1;
        algebra::Fraction optimal = type == Optimization::MAXIMIZE ? -algebra::inf : algebra::inf;
        std::queue<IPP> queue;
        std::map<algebra::Variable, algebra::Fraction> res;
        queue.push(*this);
        std::filesystem::create_directories(path);

        while (!queue.empty()) {
            IPP current = queue.front();
            std::variant<std::map<algebra::Variable, algebra::Fraction>, Solution> result =
                current.optimize_graphical(path + "/graph" + std::to_string(i++) + ".png");
            queue.pop();

            if (const std::map<algebra::Variable, algebra::Fraction>* ans = std::get_if<std::map<algebra::Variable, algebra::Fraction>>(&result)) {
                bool is_fractional = false;
                auto range = *ans | std::views::filter([](const std::pair<algebra::Variable, algebra::Fraction>& element) -> bool {
                    return element.second.denominator != 1 && element.first.variables != Z.variables;
                });
                auto itr =
                    std::ranges::max_element(range, {}, [](const std::pair<algebra::Variable, algebra::Fraction>& element) -> algebra::Fraction {
                        return element.second.numerator % element.second.denominator;
                    });
                if (itr != range.end()) {
                    const auto [variable, fraction] = *itr;
                    current.constraints.push_back(variable <= static_cast<int64_t>(fraction));
                    queue.push(current);
                    current.constraints.back() = variable >= static_cast<int64_t>(fraction) + 1;
                    queue.push(current);
                    is_fractional = true;
                }
                if (!is_fractional) {
                    const algebra::Fraction& value = ans->at(Z);

                    if (type == Optimization::MAXIMIZE && optimal < value || type == Optimization::MINIMIZE && optimal > value) {
                        optimal = value;
                        res = *ans;
                    }
                }
            }
        }
        for (const auto& [variable, fraction] : res) {
            GLOBAL_FORMATTING << algebra::Equation(variable, fraction) << std::endl;
        }
        return res;
    }
};
