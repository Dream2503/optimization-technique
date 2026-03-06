#pragma once

class algebra::Polynomial {
public:
    std::vector<Variable> expression;

    constexpr Polynomial() = default;

    Polynomial(const Fraction& fraction) : Polynomial(Variable(fraction)) {}

    Polynomial(const Variable& variable) { expression.push_back(variable); }

    Polynomial(const std::vector<Variable>& terms) : expression(terms) { std::ranges::sort(expression); }

    Polynomial operator-() const {
        Polynomial res = *this;

        for (Variable& variable : res.expression) {
            variable *= -1;
        }
        return res;
    }

    Polynomial& operator+=(const Fraction& value) { return *this += Polynomial(value); }

    Polynomial operator+(const Fraction& value) const { return *this + Polynomial(value); }

    Polynomial& operator+=(const Variable& value) { return *this += Polynomial(value); }

    Polynomial operator+(const Variable& value) const { return *this + Polynomial(value); }

    Polynomial& operator+=(const Polynomial& value) {
        for (const Variable& variable : value.expression) {
            const auto itr = std::ranges::lower_bound(expression, variable, [](const Variable& lhs, const Variable& rhs) -> bool {
                const bool lhs_const = lhs.variables.empty(), rhs_const = rhs.variables.empty();

                if (lhs_const != rhs_const) {
                    return !lhs_const;
                }
                return lhs.variables < rhs.variables;
            });

            if (itr != expression.end() && itr->variables == variable.variables) {
                itr->coefficient += variable.coefficient;

                if (itr->coefficient == 0) {
                    expression.erase(itr);
                }
            } else if (variable.coefficient != 0) {
                expression.insert(itr, variable);
            }
        }
        return *this;
    }

    Polynomial operator+(const Polynomial& value) const {
        Polynomial res = *this;
        res += value;
        return res;
    }

    Polynomial& operator-=(const Fraction& value) { return *this += -value; }

    Polynomial operator-(const Fraction& value) const { return *this + -value; }

    Polynomial& operator-=(const Variable& value) { return *this += -value; }

    Polynomial operator-(const Variable& value) const { return *this + -value; }

    Polynomial& operator-=(const Polynomial& value) { return *this += -value; }

    Polynomial operator-(const Polynomial& value) const { return *this + -value; }

    Polynomial& operator*=(const Fraction& value) { return *this *= Polynomial(value); }

    Polynomial operator*(const Fraction& value) const { return *this * Polynomial(value); }

    Polynomial& operator*=(const Variable& value) { return *this *= Polynomial(value); }

    Polynomial operator*(const Variable& value) const { return *this * Polynomial(value); }

    Polynomial& operator*=(const Polynomial& value) { return *this = *this * value; }

    Polynomial operator*(const Polynomial& value) const {
        Polynomial res;

        for (const Variable& lhs : expression) {
            for (const Variable& rhs : value.expression) {
                res += lhs * rhs;
            }
        }
        return res;
    }

    Polynomial& operator/=(const Fraction& value) { return *this /= Variable(value); }

    Polynomial operator/(const Fraction& value) const { return *this / Variable(value); }

    Polynomial& operator/=(const Variable& value) {
        for (Variable& variable : expression) {
            variable /= value;
        }
        return *this;
    }

    Polynomial operator/(const Variable& value) const {
        Polynomial res = *this;
        res /= value;
        return res;
    }

    Polynomial substitute(const std::vector<std::pair<std::string, Fraction>>& values) const {
        Polynomial res;

        for (const Variable& variable : expression) {
            res += variable.substitute(values);
        }
        return res;
    }

    bool is_fraction() const { return expression.empty() || expression.size() == 1 && expression.front().is_fraction(); }

    bool is_variable() const { return expression.size() == 1; }

    explicit operator Fraction() const {
        assert(is_fraction());

        if (expression.empty()) {
            return Fraction(0);
        }
        return static_cast<Fraction>(static_cast<Variable>(*this));
    }

    explicit operator Variable() const {
        assert(is_variable());


        return expression[0];
    }
};

inline algebra::Polynomial operator+(const algebra::Variable& lhs, const algebra::Variable& rhs) { return algebra::Polynomial(lhs) + rhs; }

inline algebra::Polynomial operator+(const algebra::Variable& lhs, const algebra::Fraction& rhs) { return algebra::Polynomial(lhs) + rhs; }

inline algebra::Polynomial operator+(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return rhs + lhs; }

inline algebra::Polynomial operator-(const algebra::Variable& lhs, const algebra::Variable& rhs) { return -rhs + lhs; }

inline algebra::Polynomial operator-(const algebra::Variable& lhs, const algebra::Fraction& rhs) { return -rhs + lhs; }

inline algebra::Polynomial operator-(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return -rhs + lhs; }

inline algebra::Polynomial operator+(const algebra::Fraction& lhs, const algebra::Polynomial& rhs) { return rhs + lhs; }

inline algebra::Polynomial operator-(const algebra::Variable& lhs, const algebra::Polynomial& rhs) { return -rhs + lhs; }

inline algebra::Polynomial operator*(const algebra::Fraction& lhs, const algebra::Polynomial& rhs) { return rhs * lhs; }

inline algebra::Polynomial operator*(const algebra::Variable& lhs, const algebra::Polynomial& rhs) { return rhs * lhs; }

namespace std {
    inline string to_string(const algebra::Polynomial& polynomial) {
        string res;

        if (!polynomial.expression.empty()) {
            res.append(to_string(polynomial.expression.front()));

            for (const algebra::Variable& variable : polynomial.expression | std::views::drop(1)) {
                res.append(variable.coefficient < 0 ? " - " : " + ").append(to_string(abs(variable)));
            }
            return res;
        }
        return "0";
    }
} // namespace std

inline std::ostream& algebra::operator<<(std::ostream& out, const Polynomial& polynomial) { return out << std::to_string(polynomial); }
