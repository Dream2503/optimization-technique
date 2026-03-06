#pragma once

class algebra::Inequation {

public:
    RelationalOperator opr;
    Polynomial lhs, rhs;

    constexpr Inequation() = default;

    Inequation(const Polynomial& polynomial, const RelationalOperator opr, const Polynomial& rhs) : opr(opr), lhs(polynomial), rhs(rhs) {}

    Inequation& operator*=(const Fraction& value) { return *this *= Polynomial(value); }

    Inequation operator*(const Fraction& value) const { return *this * Polynomial(value); }

    Inequation& operator*=(const Variable& value) { return *this *= Polynomial(value); }

    Inequation operator*(const Variable& value) const { return *this * Polynomial(value); }

    Inequation& operator*=(const Polynomial& value) {
        lhs *= value;
        rhs *= value;
        return *this;
    }

    Inequation operator*(const Polynomial& value) const {
        Inequation inequation = *this;
        inequation *= value;
        return inequation;
    }

    Inequation& operator/=(const Fraction& value) { return *this /= Variable(value); }

    Inequation operator/(const Fraction& value) const { return *this / Variable(value); }

    Inequation& operator/=(const Variable& value) {
        lhs /= value;
        rhs /= value;
        return *this;
    }

    Inequation operator/(const Variable& value) const {
        Inequation inequation = *this;
        inequation /= value;
        return inequation;
    }

    Inequation invert() const {
        Inequation res = *this;
        res.lhs *= -1;
        res.opr = detail::invert_relational_operator(res.opr);
        res.rhs *= -1;
        return res;
    }

    Inequation swap() const {
        Inequation res = *this;
        std::swap(res.lhs, res.rhs);
        res.opr = detail::invert_relational_operator(res.opr);
        return res;
    }

    Inequation substitute(const std::vector<std::pair<std::string, Fraction>>& values) const {
        Inequation res = *this;
        res.lhs = res.lhs.substitute(values);
        res.rhs = res.rhs.substitute(values);
        return res;
    }

    Inequation solve_for(const Variable& variable) const {
        Inequation res;
        res.opr = opr;

        for (const Variable& var : lhs.expression) {
            if (var.basis() == variable) {
                res.lhs += var;
            } else {
                res.rhs += -var;
            }
        }
        for (const Variable& var : rhs.expression) {
            if (var.basis() == variable) {
                res.lhs += -var;
            } else {
                res.rhs += var;
            }
        }
        if (res.lhs.is_variable()) {
            res /= static_cast<Variable>(res.lhs).coefficient;
        }
        return res;
    }

    bool is_bool() const { return lhs.is_fraction() && rhs.is_fraction(); }

    explicit operator bool() const { return detail::evaluate_relational_operator(static_cast<Fraction>(lhs), opr, static_cast<Fraction>(rhs)); }
};

class algebra::Equation : public Inequation {
public:
    Equation(const Polynomial& polynomial, const Polynomial& rhs) : Inequation(polynomial, RelationalOperator::EQ, rhs) {}

    Equation swap() const {
        Equation res = *this;
        std::swap(res.lhs, res.rhs);
        return res;
    }
};

inline algebra::Inequation operator<(const algebra::Polynomial& lhs, const algebra::Polynomial& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::LT, rhs);
}

inline algebra::Inequation operator<=(const algebra::Polynomial& lhs, const algebra::Polynomial& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::LE, rhs);
}

inline algebra::Inequation operator>(const algebra::Polynomial& lhs, const algebra::Polynomial& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::GT, rhs);
}

inline algebra::Inequation operator>=(const algebra::Polynomial& lhs, const algebra::Polynomial& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::GE, rhs);
}

inline algebra::Equation operator==(const algebra::Polynomial& lhs, const algebra::Polynomial& rhs) { return algebra::Equation(lhs, rhs); }

inline algebra::Inequation operator<(const algebra::Variable& lhs, const algebra::Fraction& rhs) {
    return algebra::Polynomial(lhs) < algebra::Polynomial(rhs);
}

inline algebra::Inequation operator<=(const algebra::Variable& lhs, const algebra::Fraction& rhs) {
    return algebra::Polynomial(lhs) <= algebra::Polynomial(rhs);
}

inline algebra::Inequation operator>(const algebra::Variable& lhs, const algebra::Fraction& rhs) {
    return algebra::Polynomial(lhs) > algebra::Polynomial(rhs);
}

inline algebra::Inequation operator>=(const algebra::Variable& lhs, const algebra::Fraction& rhs) {
    return algebra::Polynomial(lhs) >= algebra::Polynomial(rhs);
}

inline algebra::Equation operator==(const algebra::Variable& lhs, const algebra::Fraction& rhs) {
    return algebra::Polynomial(lhs) == algebra::Polynomial(rhs);
}

inline algebra::Inequation operator<(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return (rhs > lhs).swap(); }

inline algebra::Inequation operator<=(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return (rhs >= lhs).swap(); }

inline algebra::Inequation operator>(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return (rhs < lhs).swap(); }

inline algebra::Inequation operator>=(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return (rhs <= lhs).swap(); }

inline algebra::Equation operator==(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return (rhs == lhs).swap(); }


inline algebra::Inequation operator<(const algebra::Variable& lhs, const algebra::Variable& rhs) {
    return algebra::Polynomial(lhs) < algebra::Polynomial(rhs);
}

inline algebra::Inequation operator<=(const algebra::Variable& lhs, const algebra::Variable& rhs) {
    return algebra::Polynomial(lhs) <= algebra::Polynomial(rhs);
}

inline algebra::Inequation operator>(const algebra::Variable& lhs, const algebra::Variable& rhs) {
    return algebra::Polynomial(lhs) > algebra::Polynomial(rhs);
}

inline algebra::Inequation operator>=(const algebra::Variable& lhs, const algebra::Variable& rhs) {
    return algebra::Polynomial(lhs) >= algebra::Polynomial(rhs);
}

inline algebra::Equation operator==(const algebra::Variable& lhs, const algebra::Variable& rhs) {
    return algebra::Polynomial(lhs) == algebra::Polynomial(rhs);
}

inline algebra::Inequation operator<(const algebra::Polynomial& lhs, const algebra::Fraction& rhs) { return lhs < algebra::Polynomial(rhs); }

inline algebra::Inequation operator<=(const algebra::Polynomial& lhs, const algebra::Fraction& rhs) { return lhs <= algebra::Polynomial(rhs); }

inline algebra::Inequation operator>(const algebra::Polynomial& lhs, const algebra::Fraction& rhs) { return lhs > algebra::Polynomial(rhs); }

inline algebra::Inequation operator>=(const algebra::Polynomial& lhs, const algebra::Fraction& rhs) { return lhs >= algebra::Polynomial(rhs); }

inline algebra::Equation operator==(const algebra::Polynomial& lhs, const algebra::Fraction& rhs) { return lhs == algebra::Polynomial(rhs); }

inline algebra::Inequation operator<(const algebra::Fraction& lhs, const algebra::Polynomial& rhs) { return (rhs > lhs).swap(); }

inline algebra::Inequation operator<=(const algebra::Fraction& lhs, const algebra::Polynomial& rhs) { return (rhs >= lhs).swap(); }

inline algebra::Inequation operator>(const algebra::Fraction& lhs, const algebra::Polynomial& rhs) { return (rhs < lhs).swap(); }

inline algebra::Inequation operator>=(const algebra::Fraction& lhs, const algebra::Polynomial& rhs) { return (rhs <= lhs).swap(); }

inline algebra::Equation operator==(const algebra::Fraction& lhs, const algebra::Polynomial& rhs) { return (rhs == lhs).swap(); }

inline algebra::Inequation operator<(const algebra::Polynomial& lhs, const algebra::Variable& rhs) { return lhs < algebra::Polynomial(rhs); }

inline algebra::Inequation operator<=(const algebra::Polynomial& lhs, const algebra::Variable& rhs) { return lhs <= algebra::Polynomial(rhs); }

inline algebra::Inequation operator>(const algebra::Polynomial& lhs, const algebra::Variable& rhs) { return lhs > algebra::Polynomial(rhs); }

inline algebra::Inequation operator>=(const algebra::Polynomial& lhs, const algebra::Variable& rhs) { return lhs >= algebra::Polynomial(rhs); }

inline algebra::Equation operator==(const algebra::Polynomial& lhs, const algebra::Variable& rhs) { return lhs == algebra::Polynomial(rhs); }

inline algebra::Inequation operator<(const algebra::Variable& lhs, const algebra::Polynomial& rhs) { return (rhs > lhs).swap(); }

inline algebra::Inequation operator<=(const algebra::Variable& lhs, const algebra::Polynomial& rhs) { return (rhs >= lhs).swap(); }

inline algebra::Inequation operator>(const algebra::Variable& lhs, const algebra::Polynomial& rhs) { return (rhs < lhs).swap(); }

inline algebra::Inequation operator>=(const algebra::Variable& lhs, const algebra::Polynomial& rhs) { return (rhs <= lhs).swap(); }

inline algebra::Equation operator==(const algebra::Variable& lhs, const algebra::Polynomial& rhs) { return (rhs == lhs).swap(); }

namespace std {
    inline string to_string(const algebra::Inequation& inequation) {
        string res = to_string(inequation.lhs);
        res.push_back(' ');
        res.append(to_string(inequation.opr)).push_back(' ');
        return res.append(to_string(inequation.rhs));
    }
} // namespace std

inline std::ostream& algebra::operator<<(std::ostream& out, const Inequation& inequation) { return out << std::to_string(inequation); }
