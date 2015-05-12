#ifndef __SIMPLEX_H__
#define __SIMPLEX_H__

#include <vector>
#include <string>

namespace simplex {

enum class Status {
    Solved = 0,
    Unbounded = 1,
    Infeasible = -1,
    FailedToStop = -2
};

inline std::string textStatus(Status s) {
    if (s == Status::Solved)
        return "Solved";
    if (s == Status::Unbounded)
        return "Unbounded";
    if (s == Status::Infeasible)
        return "Infeasible";
    if (s == Status::FailedToStop)
        return "FailedToStop";
    return "???";
}

class ExtraStatus {
    int val;
    ExtraStatus(int v) : val(v) { }
public:
    enum ExtraStatusBits {
        Normal = 0,
        IncompleteRank = 0x1,
        Inconsistent = 0x2
    };
    ExtraStatus() { }
    ExtraStatus(ExtraStatusBits v) : ExtraStatus(static_cast<int>(v)) { }
    ExtraStatus operator|(const ExtraStatus &b) const {
        return ExtraStatus{val | b.val};
    }
    ExtraStatus operator&(const ExtraStatus &b) const {
        return ExtraStatus{val & b.val};
    }
    ExtraStatus &operator|=(const ExtraStatus &b) {
        val |= b.val;
        return *this;
    }
    ExtraStatus &operator&=(const ExtraStatus &b) {
        val &= b.val;
        return *this;
    }
    bool operator==(const ExtraStatus &o) const {
        return this->val == o.val;
    }
    bool operator!=(const ExtraStatus &o) const {
        return this->val != o.val;
    }
    bool operator!() const {
        return !val;
    }
    explicit operator bool() const {
        return val;
    }
    ExtraStatus operator~() const {
        return ExtraStatus{~val};
    }
    std::string value() const {
        std::string ret = "";
        if (val == 0)
            return "Normal";
        if (val & IncompleteRank)
            ret = "IncompleteRank";
        if (ret != "")
            ret += " | Inconsistent";
        else
            ret += "Inconsistent";
        return ret;
    }
};

/* Solve standard problem of the form
 * min  c^T x
 * s.t. Ax = b
 *      x >= 0
 * Input: A in col-major order
 * Output: if p != nullptr then *p contains objective function in optimum
 * */
Status solveProblem(const std::vector<double> &A, const std::vector<double> &b, const std::vector<double> &c,
        std::vector<int> &basis, std::vector<double> &x, double *p = nullptr, ExtraStatus *es = nullptr, bool debug = false);

}

#endif
