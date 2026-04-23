#pragma once
#include "lpp.hpp"

class optimization::ComputationalTable {
    static constexpr auto serial_class = detail::SerialClass::COMPUTATIONAL_TABLE;

    ComputationalTable() = default;

    static algebra::Fraction extract_coefficient_M(const algebra::SimplePolynomial& polynomial) {
        const auto itr = std::ranges::find(polynomial.terms, LPP::M, &algebra::Variable::basis);

        if (itr != polynomial.terms.end()) {
            return itr->coefficient;
        }
        return static_cast<algebra::Fraction>(polynomial);
    }

    void compute_zj_cj() {
        const int size = coefficient_matrix[LPP::B].size();
        zj_cj.clear();

        for (const algebra::Variable& variable : cost | std::views::keys) {
            algebra::SimplePolynomial polynomial;

            for (int i = 0; i < size; i++) {
                polynomial += cost[basis_vector[i]] * coefficient_matrix[variable][i];
            }
            zj_cj.push_back(polynomial - cost[variable]);
        }
    }

    std::map<algebra::Variable, std::vector<algebra::Fraction>> compute_next_table(const std::pair<algebra::Variable, int>& pivot, const int size) {
        std::map<algebra::Variable, std::vector<algebra::Fraction>> new_coefficient_matrix = coefficient_matrix;
        auto unit_matrix = tensor::Matrix<algebra::Fraction>::make_identity(size);
        basis_vector.erase(basis_vector.begin() + pivot.second);
        basis_vector.insert(basis_vector.begin() + pivot.second, pivot.first);

        for (int i = 0; i < size; i++) {
            new_coefficient_matrix[basis_vector[i]] = unit_matrix[i];
        }
        for (const algebra::Variable& variable : coefficient_matrix | std::views::keys |
                 std::views::filter([this](const algebra::Variable& element) -> bool { return !std::ranges::contains(basis_vector, element); })) {
            for (int i = 0; i < size; i++) {
                if (i == pivot.second) {
                    new_coefficient_matrix[variable][i] /= coefficient_matrix[pivot.first][pivot.second];
                } else {
                    new_coefficient_matrix[variable][i] = (coefficient_matrix[variable][i] * coefficient_matrix[pivot.first][pivot.second] -
                                                           coefficient_matrix[variable][pivot.second] * coefficient_matrix[pivot.first][i]) /
                        coefficient_matrix[pivot.first][pivot.second];
                }
            }
        }
        return new_coefficient_matrix;
    }

public:
    LPP lpp;
    Solution solution;
    std::vector<algebra::Variable> basis_vector;
    std::map<algebra::Variable, algebra::Variable> cost;
    std::map<algebra::Variable, std::vector<algebra::Fraction>> coefficient_matrix;
    std::vector<algebra::SimplePolynomial> zj_cj;
    std::vector<algebra::Fraction> mr;
    static constexpr int TAB_SIZE = 13;

    ComputationalTable(const LPP& lpp) : ComputationalTable(std::move(lpp), lpp.type) {}

    ComputationalTable(const LPP&& lpp, const Optimization type) : lpp(lpp), solution(Solution::UNOPTIMIZED) {
        const int size = lpp.constraints.size();
        const auto unit_matrix = tensor::Matrix<algebra::Fraction>::make_identity(size);
        this->lpp.type = type;

        for (const algebra::Variable& variable : lpp.objective.terms) {
            cost.emplace(algebra::Variable(variable.variables[0].name), variable.coefficient);
        }
        for (const algebra::Inequation& constraint : lpp.constraints) {
            for (const algebra::Variable& variable : constraint.lhs.terms) {
                cost.emplace(variable.basis(), 0);
            }
            coefficient_matrix[LPP::B].push_back(static_cast<algebra::Fraction>(constraint.rhs));
        }
        for (const algebra::Variable& variable : cost | std::views::keys) {
            coefficient_matrix.emplace(variable, std::vector<algebra::Fraction>());
        }
        for (const algebra::Inequation& constraint : lpp.constraints) {
            for (std::vector<algebra::Fraction>& fractions : coefficient_matrix | std::views::drop(1) | std::views::values) { // B
                fractions.push_back(0);
            }
            for (const algebra::Variable& variable : constraint.lhs.terms) {
                coefficient_matrix[variable.basis()].back() = variable.coefficient;
            }
        }
        for (int i = 0, j = 1; i < size; i++) {
            const auto itr =
                std::ranges::find_if(coefficient_matrix | std::views::drop(1), // B
                                     [&unit_matrix, i](const std::pair<algebra::Variable, std::vector<algebra::Fraction>>& element) -> bool {
                                         return element.second == unit_matrix[i];
                                     });

            if (itr != coefficient_matrix.end()) {
                basis_vector.push_back(itr->first);
            } else {
                const algebra::Variable variable("A" + std::to_string(j++));
                coefficient_matrix[variable] = unit_matrix[i];
                cost.emplace(variable, -LPP::M);
                basis_vector.push_back(variable);
            }
        }
    }

    ComputationalTable(const std::map<algebra::Variable, algebra::Variable>& cost, const std::vector<algebra::Variable>& basis_vector,
                       const std::map<algebra::Variable, std::vector<algebra::Fraction>>& coefficient_matrix, const Solution solution,
                       const LPP& lpp = LPP()) :
        lpp(lpp), solution(solution), basis_vector(basis_vector), cost(cost), coefficient_matrix(coefficient_matrix) {
        compute_zj_cj();
    }

    std::variant<std::vector<std::map<algebra::Variable, algebra::Fraction>>, Solution>
    get_solutions(const LPP::Method& method = LPP::Method::SIMPLEX, const LPP::ArtificialMethod artificial_method = LPP::ArtificialMethod::BIG_M,
                  const std::map<algebra::Variable, algebra::Variable>& complementary_slack = {}) {
        const auto get_solution = [this] -> std::map<algebra::Variable, algebra::Fraction> {
            const int size = basis_vector.size();
            auto range = cost |
                std::views::transform(
                             [](const std::pair<algebra::Variable, algebra::Variable>& element) -> std::pair<algebra::Variable, algebra::Fraction> {
                                 return {element.first, static_cast<algebra::Fraction>(element.second)};
                             });
            std::map<algebra::Variable, algebra::Fraction> res, substituent(range.begin(), range.end());

            for (int i = 0; i < size; i++) {
            }
            for (const algebra::Variable& variable :
                 cost | std::views::keys | std::views::filter([](const algebra::Variable& var) -> bool { return var.variables[0].name[0] != 's'; })) {
                const int idx = std::ranges::find(basis_vector, variable) - basis_vector.begin();
                substituent[variable] = res[variable] = idx < size ? coefficient_matrix[LPP::B][idx] : 0;
            }
            GLOBAL_FORMATTING << *this << std::flush;
            res[LPP::Z] = static_cast<algebra::Fraction>(lpp.objective.substitute(substituent));
            res[LPP::Z] *= lpp.type == Optimization::MINIMIZE ? -1 : 1;
            return res;
        };
        bool loop = true;
        std::vector<std::map<algebra::Variable, algebra::Fraction>> res;

        while (loop) {
            solution = optimize(LPP::Method::SIMPLEX, artificial_method, complementary_slack);

            switch (solution) {
            case Solution::OPTIMIZED:
                res.push_back(get_solution());
                loop = false;
                break;

            case Solution::INFEASIBLE:
                GLOBAL_FORMATTING << "Infeasible Solution" << std::endl;
                return Solution::INFEASIBLE;

            case Solution::UNBOUNDED:
                GLOBAL_FORMATTING << "Unbounded Solution" << std::endl;
                return Solution::UNBOUNDED;

            case Solution::ALTERNATE:
                {
                    std::map<algebra::Variable, algebra::Fraction> sol = get_solution();

                    if (!std::ranges::contains(res, sol)) {
                        res.push_back(sol);
                    } else {
                        loop = false;
                    }
                    break;
                }

            case Solution::UNOPTIMIZED:
                std::unreachable();
            }
        }
        for (const std::map<algebra::Variable, algebra::Fraction>& sol : res) {
            for (const auto& [variable, fraction] : sol) {
                GLOBAL_FORMATTING << algebra::Equation(variable, fraction) << std::endl;
            }
            GLOBAL_FORMATTING << std::endl;
        }
        return res;
    }

    Solution optimize(const LPP::Method method = LPP::Method::SIMPLEX, const LPP::ArtificialMethod artificial_method = LPP::ArtificialMethod::BIG_M,
                      const std::map<algebra::Variable, algebra::Variable>& complementary_slack = {}) {
        switch (method) {
        case LPP::Method::SIMPLEX:
            return optimize_simplex(artificial_method, complementary_slack);

        case LPP::Method::DUAL_SIMPLEX:
            return optimize_dual_simplex();
        }
        std::unreachable();
    }

    Solution optimize_simplex(const LPP::ArtificialMethod artificial_method = LPP::ArtificialMethod::BIG_M,
                              const std::map<algebra::Variable, algebra::Variable>& complementary_slack = {}) {
        if (solution != Solution::UNOPTIMIZED && solution != Solution::ALTERNATE) {
            return solution;
        }
        bool phase1 = true;
        const int size = coefficient_matrix[LPP::B].size();
        std::map<algebra::Variable, algebra::Variable> temp_cost = cost;

        if (artificial_method == LPP::ArtificialMethod::TWO_PHASE &&
            std::ranges::contains(basis_vector, 'A', [](const algebra::Variable& variable) -> char { return variable.variables[0].name[0]; })) {
            for (auto& [key, value] : cost) {
                value = algebra::Variable(key.variables[0].name[0] == 'A' ? -1 : 0);
            }
            std::erase_if(temp_cost, [](const std::pair<algebra::Variable, algebra::Variable>& element) -> bool {
                return element.first.variables[0].name[0] == 'A';
            });
        } else {
            phase1 = false;
        }
        while (true) {
            if (solution != Solution::ALTERNATE) {
                compute_zj_cj();

                if (std::ranges::all_of(zj_cj,
                                        [](const algebra::SimplePolynomial& polynomial) -> bool { return extract_coefficient_M(polynomial) >= 0; })) {
                    if (phase1 && artificial_method == LPP::ArtificialMethod::TWO_PHASE) {
                        cost = temp_cost;
                        compute_zj_cj();
                        phase1 = false;
                    } else {
                        solution = std::ranges::contains(basis_vector, 'A',
                                                         [](const algebra::Variable& variable) -> char { return variable.variables[0].name[0]; })
                            ? Solution::INFEASIBLE
                            : Solution::OPTIMIZED;

                        if (artificial_method != LPP::ArtificialMethod::TWO_PHASE) {
                            const int zj_cj_size = zj_cj.size();
                            auto itr = std::next(coefficient_matrix.begin()); // B

                            for (int i = 0; i < zj_cj_size; i++) {
                                if (zj_cj[i].is_fraction() && static_cast<algebra::Fraction>(zj_cj[i]) == 0 &&
                                    !std::ranges::contains(basis_vector, itr->first)) {
                                    solution = Solution::ALTERNATE;
                                    break;
                                }
                                ++itr;
                            }
                        }
                        return solution;
                    }
                }
            }
            auto ev = std::next(coefficient_matrix.begin()); // B
            mr.clear();

            if (solution == Solution::ALTERNATE) {
                const int zj_cj_size = zj_cj.size();

                for (int i = 0; i < zj_cj_size; i++) {
                    if (static_cast<algebra::Fraction>(zj_cj[i]) == 0 && !std::ranges::contains(basis_vector, ev->first)) {
                        solution = Solution::UNOPTIMIZED;
                        break;
                    }
                    ++ev;
                }
            } else {
                if (complementary_slack.empty()) {
                    ev = std::next(coefficient_matrix.begin(), std::ranges::min_element(zj_cj, {}, extract_coefficient_M) - zj_cj.begin() + 1); // B
                } else {
                    auto x = std::views::zip(cost | std::views::keys, zj_cj | std::views::transform(extract_coefficient_M)) |
                        std::views::transform([](auto&& element) -> std::pair<algebra::Variable, algebra::Fraction> {
                                 return std::pair{std::get<0>(element), std::get<1>(element)};
                             }) |
                        std::ranges::to<std::vector>();
                    std::ranges::sort(x, {}, &std::pair<algebra::Variable, algebra::Fraction>::second);
                    ev = coefficient_matrix.find(*std::ranges::begin(
                        x | std::views::keys | std::views::drop_while([this, &complementary_slack](const algebra::Variable& variable) -> bool {
                            return std::ranges::contains(basis_vector, variable) ||
                                std::ranges::contains(basis_vector, complementary_slack.at(variable));
                        })));
                    phase1 = false;
                }
            }
            for (int i = 0; i < size; i++) {
                mr.push_back(ev->second[i] <= 0 ? algebra::inf : coefficient_matrix[LPP::B][i] / ev->second[i]);
            }
            int lv = std::ranges::min_element(mr) - mr.begin();
            bool is_unbounded = true;
            GLOBAL_FORMATTING << *this << std::flush;

            if (!std::ranges::contains(basis_vector, 'A', [](const algebra::Variable& variable) -> char { return variable.variables[0].name[0]; })) {
                for (int k = 0; k < size && mr[lv] != algebra::inf; k++) {
                    std::vector<int> candidates;

                    for (int i = 0; i < size; i++) {
                        if (mr[i] == mr[lv]) {
                            candidates.push_back(i);
                        }
                    }
                    if (candidates.size() > 1) {
                        for (const int candidate : candidates) {
                            mr[candidate] = coefficient_matrix[basis_vector[k]][candidate] / coefficient_matrix[ev->first][k];
                        }
                        lv = std::ranges::min_element(mr) - mr.begin();
                    } else {
                        is_unbounded = false;
                        break;
                    }
                }
            } else {
                is_unbounded = false;
            }
            if (is_unbounded || mr[lv] == algebra::inf) {
                return solution = Solution::UNBOUNDED;
            }
            if (basis_vector[lv].variables[0].name[0] == 'A' || cost[basis_vector[lv]].variables == LPP::M.variables) {
                coefficient_matrix.erase(basis_vector[lv]);
                cost.erase(basis_vector[lv]);
                auto itr = std::ranges::find(lpp.objective.terms, basis_vector[lv], &algebra::Variable::basis);

                if (itr != lpp.objective.terms.end()) {
                    lpp.objective.terms.erase(itr);
                }
            }
            coefficient_matrix = compute_next_table({ev->first, lv}, size);
        }
    }

    Solution optimize_dual_simplex() {
        const int size = coefficient_matrix[LPP::B].size();
        solution = Solution::UNBOUNDED;

        while (true) {
            compute_zj_cj();

            if (std::ranges::all_of(
                    zj_cj, [](const algebra::SimplePolynomial& polynomial) -> bool { return static_cast<algebra::Fraction>(polynomial) >= 0; }) &&
                std::ranges::all_of(coefficient_matrix[LPP::B], [](const algebra::Fraction& fraction) -> bool { return fraction >= 0; })) {
                return solution = Solution::OPTIMIZED;
            }
            const int lv = std::ranges::min_element(coefficient_matrix[LPP::B]) - coefficient_matrix[LPP::B].begin();
            const auto begin = std::next(coefficient_matrix.begin()); // B
            const auto ev = std::ranges::max_element(
                coefficient_matrix | std::views::drop(1),
                [this, lv, &begin](const std::pair<algebra::Variable, std::vector<algebra::Fraction>>& lhs,
                                   const std::pair<algebra::Variable, std::vector<algebra::Fraction>>& rhs) -> bool {
                    const algebra::Fraction &lhs_value = lhs.second[lv], rhs_value = rhs.second[lv];

                    if (lhs_value >= 0) {
                        return true;
                    }
                    if (rhs_value >= 0) {
                        return false;
                    }
                    const int lhs_idx = std::distance(begin, coefficient_matrix.find(lhs.first));
                    const int rhs_idx = std::distance(begin, coefficient_matrix.find(rhs.first));
                    return static_cast<algebra::Fraction>(zj_cj[lhs_idx]) / lhs_value < static_cast<algebra::Fraction>(zj_cj[rhs_idx]) / rhs_value;
                });
            GLOBAL_FORMATTING << *this;
            coefficient_matrix = compute_next_table({ev->first, lv}, size);
        }
    }

    std::vector<algebra::Interval> cost_variation() {
        std::vector<algebra::Interval> res;

        if (solution == Solution::UNOPTIMIZED) {
            optimize_simplex();
        }
        for (const auto& [variable, value] : cost | std::views::filter([](const std::pair<algebra::Variable, algebra::Variable>& element) -> bool {
                                                 return element.first.variables[0].name[0] != 's';
                                             })) {
            const auto itr = std::ranges::find(basis_vector, variable);
            algebra::Variable var("C" + variable.variables[0].name);

            if (itr != basis_vector.end()) {
                const int idx = itr - basis_vector.begin();
                int i = 0;
                algebra::Fraction min = algebra::inf, max = -algebra::inf;

                for (const auto& [name, fractions] : coefficient_matrix | std::views::drop(1)) { // B
                    if (!std::ranges::contains(basis_vector, name)) {
                        if (fractions[idx] > 0) {
                            max = std::max(max, -static_cast<algebra::Fraction>(zj_cj[i] / fractions[idx]));
                        } else if (fractions[idx] < 0) {
                            min = std::min(min, -static_cast<algebra::Fraction>(zj_cj[i] / fractions[idx]));
                        }
                    }
                    i++;
                }
                res.push_back(max + value < var < min + value);
            } else {
                res.push_back(-algebra::inf < var <
                              zj_cj[std::distance(coefficient_matrix.begin(), coefficient_matrix.find(variable)) - 1] + value); // B
            }
        }
        for (const algebra::Interval& interval : res) {
            GLOBAL_FORMATTING << interval << std::endl;
        }
        return res;
    }

    std::vector<algebra::Interval> RHS_variation() {
        const int size = basis_vector.size();
        int i = 0;
        std::vector<algebra::Interval> res;

        for (const auto& [name, fractions] :
             coefficient_matrix | std::views::filter([](const std::pair<algebra::Variable, std::vector<algebra::Fraction>>& element) -> bool {
                 return element.first.variables[0].name[0] == 's';
             })) {
            algebra::Variable var("B" + std::to_string(i + 1));
            algebra::Fraction min = algebra::inf, max = -algebra::inf;

            for (int j = 0; j < size; j++) {
                if (fractions[j] > 0) {
                    max = std::max(max, -coefficient_matrix[LPP::B][j] / fractions[j]);
                } else if (fractions[j] < 0) {
                    min = std::min(min, -coefficient_matrix[LPP::B][j] / fractions[j]);
                }
            }
            res.push_back(max + lpp.constraints[i].rhs < var < min + lpp.constraints[i].rhs);
            i++;
        }
        for (const algebra::Interval& interval : res) {
            GLOBAL_FORMATTING << interval << std::endl;
        }
        return res;
    }

    void add_variable(const algebra::Variable& variable, const tensor::Matrix<algebra::Fraction>& coefficients) {
        int i = 0;
        const int size = basis_vector.size();
        tensor::Matrix<algebra::Fraction> res(size, size);
        lpp.objective += variable;
        GLOBAL_FORMATTING << *this;

        for (const auto& [name, fractions] :
             coefficient_matrix | std::views::filter([](const std::pair<algebra::Variable, std::vector<algebra::Fraction>>& element) -> bool {
                 return element.first.variables[0].name[0] == 's';
             })) {
            for (int j = 0; j < size; j++) {
                res[i, j] = fractions[j];
            }
            i++;
        }
        res = coefficients * res;
        coefficient_matrix[variable.basis()] = res[0];
        cost[variable.basis()] = variable.coefficient;
        solution = Solution::UNOPTIMIZED;
        GLOBAL_FORMATTING << *this;
    }

    void remove_variable(const algebra::Variable& variable) {
        GLOBAL_FORMATTING << *this;

        if (std::ranges::contains(basis_vector, variable)) {
            cost[variable.basis()] = -LPP::M;
            solution = Solution::UNOPTIMIZED;
        } else {
            const int idx = std::ranges::distance(cost.begin(), cost.find(variable));
            lpp.objective -= cost[variable] * variable;
            zj_cj.erase(zj_cj.begin() + idx);
            cost.erase(variable);
            coefficient_matrix.erase(variable);
        }
        GLOBAL_FORMATTING << *this;
    }

    void add_constraint(const algebra::Inequation& inequation) {
        std::map<algebra::Variable, algebra::Fraction> substituent =
            std::get<std::vector<std::map<algebra::Variable, algebra::Fraction>>>(get_solutions())[0];
        GLOBAL_FORMATTING << *this;

        if (static_cast<bool>(inequation.substitute(substituent))) {
            return;
        }
        auto range = coefficient_matrix | std::views::keys |
            std::views::filter([](const algebra::Variable& variable) -> bool { return variable.variables[0].name[0] == 's'; }) |
            std::views::transform([](const algebra::Variable& variable) -> int { return std::stoi(variable.variables[0].name.substr(1)); });
        const int size = basis_vector.size() + 1;
        const algebra::Variable slack("s" + std::to_string(*std::ranges::max_element(range) + 1));
        algebra::SimplePolynomial transformation;
        coefficient_matrix[LPP::B].push_back(static_cast<algebra::Fraction>(inequation.rhs));

        for (const algebra::Variable& variable : inequation.lhs.terms) {
            coefficient_matrix[variable.basis()].push_back(variable.coefficient);

            if (std::ranges::contains(basis_vector, variable.basis())) {
                transformation += variable;
            }
        }
        for (algebra::Variable& variable : transformation.terms) {
            const algebra::Variable var("R" + std::to_string(std::ranges::find(basis_vector, variable.basis()) - basis_vector.begin() + 1));
            substituent[var] = 0;
            variable.variables = std::move(var.variables);
        }
        for (std::vector<algebra::Fraction>& fractions : coefficient_matrix | std::views::values) {
            if (fractions.size() < size) {
                fractions.push_back(0);
            }
        }
        basis_vector.push_back(slack);
        cost.emplace(slack, 0);
        coefficient_matrix.emplace(slack, std::vector<algebra::Fraction>(size, 0));
        coefficient_matrix[slack].back() = 1;

        for (std::vector<algebra::Fraction>& fractions : coefficient_matrix | std::views::values) {
            for (int i = 0; i < size - 1; i++) {
                substituent[basis_vector[i]] = fractions[i];
            }
            fractions.back() -= static_cast<algebra::Fraction>(transformation.substitute(substituent, false));
        }
        solution = Solution::UNOPTIMIZED;
        GLOBAL_FORMATTING << *this;
    }

    std::string to_latex() const {
        const int size = basis_vector.size(), columns = 3 + coefficient_matrix.size();
        std::string res("\\begin{array}{");

        for (int i = 0; i < columns; i++) {
            res.append("c|");
        }
        res.pop_back();
        res.append("}\n &  &  & ");

        for (const algebra::Variable& variable : coefficient_matrix | std::views::drop(1) | std::views::keys) { // B
            res.append(cost.at(variable).to_latex()).append(" & ");
        }
        res.append("\\\\\n\\hline\n");

        for (const char* str : {"BV", "C", "B"}) {
            res.append("\\text{").append(str).append("} & ");
        }
        for (const algebra::Variable& variable : coefficient_matrix | std::views::drop(1) | std::views::keys) { // B
            res.append(variable.to_latex()).append(" & ");
        }
        res.append("\\text{MR} \\\\\n\\hline\n");

        for (int i = 0; i < size; i++) {
            res.append(basis_vector[i].to_latex()).append(" & ").append(cost.at(basis_vector[i]).to_latex()).append(" & ");

            for (const std::vector<algebra::Fraction>& fractions : coefficient_matrix | std::views::values) {
                res.append(fractions[i].to_latex()).append(" & ");
            }
            if (!mr.empty()) {
                res.append(mr[i] == algebra::inf ? "-" : mr[i].to_latex()).append(" & ");
            }
            res.pop_back();
            res.pop_back();
            res.append("\\\\\n");
        }
        res.append("\\hline\n &  & Z_j-C_j & ");

        for (const algebra::SimplePolynomial& polynomial : zj_cj) {
            res.append(polynomial.to_latex()).append(" & ");
        }
        return res.append("\\\\\n\\end{array}\n");
    }

    std::string to_html() const {
        const int size = basis_vector.size(), columns = 3 + coefficient_matrix.size();
        std::string res("<table style='border-collapse: collapse; margin-left: auto; margin-right: auto;'><style>td:first-child { border-left: none; "
                        "} td { padding: 8px 15px;text-align: center;border-left: 1px solid black; }</style><tr><td></td><td></td><td></td>");

        for (const algebra::Variable& variable : coefficient_matrix | std::views::drop(1) | std::views::keys) { // B
            res.append("<td>").append(cost.at(variable).to_html()).append("</td>");
        }
        res.append("<td></td></tr><tr style='border-top: 1px solid black; border-bottom: 1px solid black;'><td>BV</td><td>C</td><td>B</td>");

        for (const algebra::Variable& variable : coefficient_matrix | std::views::drop(1) | std::views::keys) { // B
            res.append("<td><math>").append(variable.to_html()).append("</math></td>");
        }
        res.append("<td>MR</td></tr>");

        for (int i = 0; i < size; i++) {
            res.append(i == size - 1 ? "<tr style='border-bottom: 1px solid black;'><td><math>" : "<tr><td><math>")
                .append(basis_vector[i].to_html())
                .append("</math></td><td><math>")
                .append(cost.at(basis_vector[i]).to_html())
                .append("</math></td>");

            for (const std::vector<algebra::Fraction>& fractions : coefficient_matrix | std::views::values) {
                res.append("<td><math>").append(fractions[i].to_html()).append("</math></td>");
            }
            if (!mr.empty()) {
                res.append("<td><math>").append(mr[i] == algebra::inf ? "<mo>&ndash;</mo>" : mr[i].to_html()).append("</math></td>");
            }
        }
        res.append("<tr><td></td><td></td><td><math><msub><mi>Z</mi><mi>j</mi></msub><mo>-</mo><msub><mi>C</mi><mi>j</mi></msub></math></td>");

        for (const algebra::SimplePolynomial& polynomial : zj_cj) {
            res.append("<td><math>").append(polynomial.to_html()).append("</math></td>");
        }
        return res.append("<td></td></table>");
    }

    void serialize(std::ofstream& out) const {
        out.write(reinterpret_cast<const char*>(&serial_class), sizeof(serial_class));
        lpp.serialize(out);
        out.write(reinterpret_cast<const char*>(&solution), sizeof(solution));
        size_t size = basis_vector.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));

        for (const algebra::Variable& variable : basis_vector) {
            variable.serialize(out);
        }
        size = cost.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));

        for (const auto& [variable, value] : cost) {
            variable.serialize(out);
            value.serialize(out);
        }
        size = coefficient_matrix.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));

        for (const auto& [variable, fractions] : coefficient_matrix) {
            variable.serialize(out);

            size = fractions.size();
            out.write(reinterpret_cast<const char*>(&size), sizeof(size));

            for (const algebra::Fraction& fraction : fractions) {
                fraction.serialize(out);
            }
        }
        size = zj_cj.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));

        for (const algebra::SimplePolynomial& polynomial : zj_cj) {
            polynomial.serialize(out);
        }
        size = mr.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));

        for (const algebra::Fraction& fraction : mr) {
            fraction.serialize(out);
        }
    }

    static ComputationalTable deserialize(std::ifstream& in) {
        detail::SerialClass type;
        in.read(reinterpret_cast<char*>(&type), sizeof(type));
        assert(type == serial_class);

        ComputationalTable res;
        size_t size, len;
        res.lpp = LPP::deserialize(in);
        in.read(reinterpret_cast<char*>(&res.solution), sizeof(res.solution));
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        res.basis_vector.reserve(size);

        for (int i = 0; i < size; i++) {
            res.basis_vector.push_back(algebra::Variable::deserialize(in));
        }
        in.read(reinterpret_cast<char*>(&size), sizeof(size));

        for (int i = 0; i < size; i++) {
            algebra::Variable variable = algebra::Variable::deserialize(in), value = algebra::Variable::deserialize(in);
            res.cost.emplace(variable, value);
        }
        in.read(reinterpret_cast<char*>(&size), sizeof(size));

        for (int i = 0; i < size; i++) {
            algebra::Variable variable = algebra::Variable::deserialize(in);
            std::vector<algebra::Fraction> fractions;
            in.read(reinterpret_cast<char*>(&len), sizeof(len));
            fractions.reserve(len);

            for (int j = 0; j < len; j++) {
                fractions.push_back(algebra::Fraction::deserialize(in));
            }
            res.coefficient_matrix.emplace(variable, fractions);
        }
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        res.zj_cj.reserve(size);

        for (int i = 0; i < size; i++) {
            res.zj_cj.push_back(algebra::SimplePolynomial::deserialize(in));
        }
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        res.mr.reserve(size);

        for (int i = 0; i < size; i++) {
            res.mr.push_back(algebra::Fraction::deserialize(in));
        }
        return res;
    }
};

namespace std {
    inline string to_string(const optimization::ComputationalTable& computational_table) {
        const int size = computational_table.basis_vector.size(), columns = 3 + computational_table.coefficient_matrix.size();
        auto get_partition = [columns] -> std::string {
            std::string str("+");

            for (int i = 0; i < columns; i++) {
                str.append(optimization::ComputationalTable::TAB_SIZE, '-').push_back('+');
            }
            return str.append("\n");
        };
        auto get_format = [](const auto& value) -> std::string {
            const std::string str = std::to_string(value);
            const int remaining = optimization::ComputationalTable::TAB_SIZE - str.size();
            return std::string(remaining / 2, ' ').append(str).append(remaining - remaining / 2, ' ');
        };
        std::string res(" ");

        for (int i = 0; i < 3; i++) {
            res.append(optimization::ComputationalTable::TAB_SIZE + 1, ' ');
        }
        for (const algebra::Variable& variable : computational_table.coefficient_matrix | std::views::drop(1) | std::views::keys) { // B
            res.append(get_format(computational_table.cost.at(variable))).push_back(' ');
        }
        res.append("\n").append(get_partition()).push_back('|');

        for (const std::string string : {"BV", "C", "B"}) {
            res.append(get_format(string)).push_back('|');
        }
        for (const algebra::Variable& variable : computational_table.coefficient_matrix | std::views::drop(1) | std::views::keys) { // B
            res.append(get_format(variable)).push_back('|');
        }
        res.append(get_format(std::string("MR")).append("|\n")).append(get_partition());

        for (int i = 0; i < size; i++) {
            res.append("|")
                .append(get_format(computational_table.basis_vector[i]))
                .append("|")
                .append(get_format(computational_table.cost.at(computational_table.basis_vector[i])))
                .push_back('|');

            for (const std::vector<algebra::Fraction>& fractions : computational_table.coefficient_matrix | std::views::values) {
                res.append(get_format(fractions[i])).push_back('|');
            }
            if (!computational_table.mr.empty()) {
                res.append(get_format(computational_table.mr[i] == algebra::inf ? "-ve" : std::to_string(computational_table.mr[i]))).push_back('|');
            } else {
                res.append(optimization::ComputationalTable::TAB_SIZE, ' ').push_back('|');
            }
            res.push_back('\n');
        }
        res.append(get_partition()).push_back('|');

        for (int i = 0; i < 2; i++) {
            res.append(optimization::ComputationalTable::TAB_SIZE, ' ').push_back('|');
        }
        res.append(get_format(std::string("Zj-Cj"))).push_back('|');

        for (const algebra::SimplePolynomial& polynomial : computational_table.zj_cj) {
            res.append(get_format(polynomial)).push_back('|');
        }
        return res.append(optimization::ComputationalTable::TAB_SIZE, ' ').append("|\n").append(get_partition()).append("\n");
    }
} // namespace std

inline std::ostream& optimization::operator<<(std::ostream& out, const ComputationalTable& computational_table) {
    return out << std::to_string(computational_table);
}

inline optimization::LPP optimization::LPP::dual(const std::string& basis) const {
    LPP canonical = *this;
    GLOBAL_FORMATTING << canonical;

    for (algebra::Inequation& constraint : canonical.constraints) {
        if (type == Optimization::MAXIMIZE && constraint.opr == algebra::RelationalOperator::GE ||
            type == Optimization::MINIMIZE && constraint.opr == algebra::RelationalOperator::LE) {
            constraint = constraint.invert();
        }
    }
    const int objective_size = canonical.objective.terms.size(), constraints_size = canonical.constraints.size();
    LPP res;
    ComputationalTable computational_table(canonical);
    std::erase_if(computational_table.coefficient_matrix, [](const std::pair<algebra::Variable, std::vector<algebra::Fraction>>& element) -> bool {
        return element.first.variables[0].name[0] == 'A';
    });
    std::erase_if(computational_table.cost, [](const std::pair<algebra::Variable, algebra::Variable>& element) -> bool {
        return !element.second.variables.empty() && element.second.variables[0].name == "M";
    });
    res.type = canonical.type == Optimization::MAXIMIZE ? Optimization::MINIMIZE : Optimization::MAXIMIZE;
    res.constraints.resize(objective_size);
    res.restrictions.reserve(constraints_size);

    for (int i = 0; i < constraints_size; i++) {
        auto itr = std::next(computational_table.coefficient_matrix.begin()); // B
        algebra::Variable variable(basis + std::to_string(i + 1));
        res.objective += computational_table.coefficient_matrix[B][i] * variable;

        for (int j = 0; j < objective_size; j++) {
            res.constraints[j].lhs += itr->second[i] * variable;
            ++itr;
        }
        res.restrictions.push_back(canonical.constraints[i].opr == algebra::RelationalOperator::EQ ? unrestrict(variable) : variable >= 0);
    }
    auto itr = computational_table.cost.begin();

    for (int i = 0; i < objective_size; i++) {
        res.constraints[i].opr =
            canonical.restrictions[i].lhs.is_fraction() && static_cast<algebra::Fraction>(canonical.restrictions[i].lhs) == algebra::inf ||
                canonical.restrictions[i].rhs.is_fraction() && static_cast<algebra::Fraction>(canonical.restrictions[i].rhs) == algebra::inf
            ? algebra::RelationalOperator::EQ
            : res.type == Optimization::MAXIMIZE ? algebra::RelationalOperator::LE
                                                 : algebra::RelationalOperator::GE;
        res.constraints[i].rhs = itr->second;
        ++itr;
    }
    GLOBAL_FORMATTING << "Dual:" << std::endl << res;
    return res;
}

inline std::variant<std::vector<std::map<algebra::Variable, algebra::Fraction>>, optimization::Solution>
optimization::LPP::tabular_optimize(const Method method, const ArtificialMethod artificial_method,
                                    const std::map<algebra::Variable, algebra::Variable>& complementary_slack) const {
    switch (method) {
    case Method::SIMPLEX:
        return ComputationalTable(standardize(method), type).get_solutions(method, artificial_method, complementary_slack);

    case Method::DUAL_SIMPLEX:
        return ComputationalTable(canonicalize().standardize(method), type).get_solutions(method);
    }
    std::unreachable();
}
