#pragma once

class algebra::Variable {
    struct Var {
        std::string name;
        Fraction exponent;

        constexpr std::strong_ordering operator<=>(const Var& value) const {
            return std::tuple(name, -exponent) <=> std::tuple(value.name, -value.exponent);
        }

        constexpr bool operator==(const Var&) const = default;
    };

public:
    Fraction coefficient;
    std::vector<Var> variables;

    Variable() = default;

    Variable(const std::string& name) : coefficient(1), variables({{name, 1}}) {}

    Variable(const Fraction& coefficient) : coefficient(coefficient) {}

    Variable operator-() const {
        Variable res = *this;
        res.coefficient = -res.coefficient;
        return res;
    }

    Variable& operator*=(const Fraction& value) { return *this *= Variable(value); }

    Variable operator*(const Fraction& value) const { return *this * Variable(value); }

    Variable& operator*=(const Variable& value) {
        coefficient *= value.coefficient;

        for (const Var& var : value.variables) {
            const auto itr = std::ranges::lower_bound(variables, var);

            if (itr != variables.end() && itr->name == var.name) {
                if ((itr->exponent += var.exponent) == 0) {
                    variables.erase(itr);
                }
            } else {
                variables.insert(itr, var);
            }
        }
        return *this;
    }

    Variable operator*(const Variable& value) const {
        Variable res = *this;
        res *= value;
        return res;
    }

    Variable& operator/=(const Fraction& value) { return *this /= Variable(value); }

    Variable operator/(const Fraction& value) const { return *this / Variable(value); }

    Variable& operator/=(const Variable& value) {
        coefficient /= value.coefficient;

        for (const Var& var : value.variables) {
            const auto itr = std::ranges::lower_bound(variables, var);

            if (itr != variables.end() && itr->name == var.name) {
                if ((itr->exponent -= var.exponent) == 0) {
                    variables.erase(itr);
                }
            } else {
                variables.emplace(itr, var.name, -var.exponent);
            }
        }
        return *this;
    }

    Variable operator/(const Variable& value) const {
        Variable res = *this;
        res /= value;
        return res;
    }

    Variable& operator^=(const Fraction& value) {
        coefficient ^= value;

        for (auto& [_, exponent] : variables) {
            exponent *= value;
        }
        return *this;
    }

    Variable operator^(const Fraction& value) const {
        Variable variable = *this;
        variable ^= value;
        return variable;
    }

    constexpr std::strong_ordering operator<=>(const Variable& value) const {
        const bool is_const = variables.empty(), value_const = value.variables.empty();

        if (is_const != value_const) {
            return is_const ? std::strong_ordering::greater : std::strong_ordering::less;
        }
        return std::tie(variables, coefficient) <=> std::tie(value.variables, value.coefficient);
    }

    constexpr bool operator==(const Variable&) const = default;

    Variable substitute(const std::vector<std::pair<std::string, Fraction>>& values) const {
        Variable res = *this;

        for (const auto& [name, value] : values) {
            const auto itr = std::ranges::lower_bound(res.variables, name, {}, &Var::name);

            if (itr != res.variables.end() && itr->name == name) {
                res.coefficient *= value ^ itr->exponent;
                res.variables.erase(itr);
            }
        }
        return res;
    }

    Variable basis() const {
        Variable res = *this;
        res.coefficient = 1;

        for (auto& [_, exponent] : res.variables) {
            exponent = 1;
        }
        return res;
    }

    bool is_fraction() const { return variables.empty(); }

    constexpr explicit operator Fraction() const {
        assert(is_fraction());
        return coefficient;
    }
};

constexpr algebra::Variable operator*(const algebra::Fraction& value, const algebra::Variable& variable) { return variable * value; }

constexpr algebra::Variable operator/(const algebra::Fraction& value, const algebra::Variable& variable) { return variable / value; }

namespace std {
    inline algebra::Variable abs(algebra::Variable variable) {
        variable.coefficient = abs(variable.coefficient);
        return variable;
    }

    inline string to_string(const algebra::Variable& variable) {
        auto convert = [](const algebra::Variable& var) -> string {
            std::string res;

            for (const auto& [name, exponent] : var.variables) {
                if (exponent != 1) {
                    res.push_back('(');
                }
                res.append(name);

                if (exponent != 1) {
                    res.push_back('^');

                    if (exponent.denominator != 1) {
                        res.push_back('(');
                        res.append(std::to_string(exponent)).push_back(')');
                    } else {
                        res.append(std::to_string(exponent));
                    }
                }
                if (exponent != 1) {
                    res.push_back(')');
                }
            }
            return res;
        };
        if (variable.variables.empty()) {
            return to_string(variable.coefficient);
        }
        if (variable.coefficient == 0) {
            return "0";
        }
        if (variable.coefficient == 1) {
            return convert(variable);
        }
        if (variable.coefficient == -1) {
            return '-' + convert(variable);
        }
        if (variable.coefficient.denominator != 1) {
            string res;

            if (abs(variable.coefficient.numerator) == 1) {
                if (variable.coefficient.numerator == -1) {
                    res.push_back('-');
                }
                res.append(convert(variable)).push_back('/');
                res.append(to_string(variable.coefficient.denominator));
            } else {
                res.append(to_string(variable.coefficient.numerator)).append(convert(variable)).push_back('/');
                res.append(to_string(variable.coefficient.denominator));
            }
            return res;
        }
        return to_string(variable.coefficient) + convert(variable);
    }
} // namespace std

inline std::ostream& algebra::operator<<(std::ostream& out, const Variable& variable) { return out << std::to_string(variable); }
