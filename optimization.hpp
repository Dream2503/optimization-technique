#pragma once
#include <iomanip>
#include <queue>
#include "tensor/tensor.hpp"

namespace optimization {
    inline algebra::detail::FormatSettings& GLOBAL_FORMATTING = algebra::GLOBAL_FORMATTING;

    enum class Optimization : bool { MINIMIZE, MAXIMIZE };
    enum class Solution : uint8_t { UNOPTIMIZED, OPTIMIZED, INFEASIBLE, UNBOUNDED, ALTERNATE };

    class LPP;
    std::ostream& operator<<(std::ostream&, const LPP&);

    class ComputationalTable;
    std::ostream& operator<<(std::ostream&, const ComputationalTable&);

    class IPP;
    class NLPP;
    class QPP;

    std::vector<std::map<algebra::Variable, algebra::Fraction>> basic_feasible_solutions(const std::vector<algebra::Equation>&);

    namespace detail {
        enum class SerialClass : uint8_t { LPP, COMPUTATIONAL_TABLE };
    }
} // namespace optimization

#include "src/lpp.hpp"
#include "src/computation_table.hpp"
#include "src/ipp.hpp"
#include "src/nlpp.hpp"
#include "src/qpp.hpp"