#pragma once
#include <filesystem>
#include <iomanip>
#include <map>
#include <queue>
#include "linear-algebra/linalg.hpp"

namespace optimization {
    inline algebra::detail::FormatSettings GLOBAL_FORMATTING;

    enum class Optimization : bool { MINIMIZE, MAXIMIZE };
    enum class Solution : uint8_t { UNOPTIMIZED, OPTIMIZED, INFEASIBLE, UNBOUNDED, ALTERNATE };
    class LPP;
    class ComputationalTable;
    class IPP;

    std::vector<std::map<algebra::Variable, algebra::Fraction>> basic_feasible_solutions(const std::vector<algebra::Equation>&);
} // namespace optimization

#include "src/lpp.hpp"
#include "src/computation_table.hpp"
#include "src/ipp.hpp"
