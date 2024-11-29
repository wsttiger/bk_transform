#include <set>
#include <complex>
#include <iostream>
#include <unordered_set>

#include "spin_op.h"

std::unordered_set<std::size_t> occupation_set(std::size_t index) {
    std::unordered_set<std::size_t> indices;
    index += 1;
    indices.insert(index - 1);
    
    std::size_t parent = index & (index - 1);
    index -= 1;
    
    while (index != parent) {
        indices.insert(index - 1);
        index &= index - 1;
    }
    
    return indices;
}

std::unordered_set<std::size_t> parity_set(std::size_t index) {
    std::unordered_set<std::size_t> indices;
    
    while (index > 0) {
        indices.insert(index - 1);
        index &= index - 1;
    }
    
    return indices;
}

std::unordered_set<std::size_t> update_set(std::size_t index, std::size_t n_qubits) {
    std::unordered_set<std::size_t> indices;
    
    index += 1;
    index += index & -index;
    
    while (index <= n_qubits) {
        indices.insert(index - 1);
        index += index & -index;
    }
    
    return indices;
}

std::unordered_set<int> remainder_set(int index) {
    std::unordered_set<int> result;
    auto parity = parity_set(index);
    auto occupation = occupation_set(index);
    
    for (const auto& elem : parity) {
        if (occupation.find(elem) == occupation.end()) {
            result.insert(elem);
        }
    }
    return result;
}

template<typename T>
std::unordered_set<T> symmetric_difference(const std::unordered_set<T>& set1, const std::unordered_set<T>& set2) {
    std::unordered_set<T> result;
    for (const auto& elem : set1) {
        if (set2.find(elem) == set2.end()) result.insert(elem);
    }
    for (const auto& elem : set2) {
        if (set1.find(elem) == set1.end()) result.insert(elem);
    }
    return result;
}

std::unordered_set<std::size_t> remainder_set(std::size_t index) {
    return symmetric_difference(parity_set(index), occupation_set(index));
}

std::unordered_set<std::size_t> F_ij_set(std::size_t i, std::size_t j) {
    return symmetric_difference(occupation_set(i), occupation_set(j));
}

std::unordered_set<std::size_t> P0_ij_set(std::size_t i, std::size_t j) {
    return symmetric_difference(parity_set(i), parity_set(j));
}

std::unordered_set<std::size_t> P1_ij_set(std::size_t i, std::size_t j) {
    return symmetric_difference(parity_set(i), remainder_set(j));
}

std::unordered_set<std::size_t> P2_ij_set(std::size_t i, std::size_t j) {
    return symmetric_difference(remainder_set(i), parity_set(j));
}

std::unordered_set<std::size_t> P3_ij_set(std::size_t i, std::size_t j) {
    return symmetric_difference(remainder_set(i), remainder_set(j));
}

std::unordered_set<std::size_t> U_ij_set(std::size_t i, std::size_t j, std::size_t n_qubits) {
    return symmetric_difference(update_set(i, n_qubits), update_set(j, n_qubits));
}

std::unordered_set<std::size_t> alpha_set(std::size_t i, std::size_t j, std::size_t n_qubits) {
    std::unordered_set<std::size_t> result;
    auto set1 = update_set(i, n_qubits);
    auto set2 = parity_set(j);
    for (const auto& elem : set1) {
        if (set2.find(elem) != set2.end()) result.insert(elem);
    }
    return result;
}

std::unordered_set<std::size_t> U_diff_a_set(std::size_t i, std::size_t j, std::size_t n_qubits) {
    return symmetric_difference(U_ij_set(i, j, n_qubits), alpha_set(i, j, n_qubits));
}

std::unordered_set<std::size_t> P0_ij_diff_a_set(std::size_t i, std::size_t j, std::size_t n_qubits) {
    return symmetric_difference(P0_ij_set(i, j), alpha_set(i, j, n_qubits));
}

/**
 * Algebraic expressions for general products of the form a_i^d a_j term in
 * the Bravyi-Kitaev basis. These expressions vary in form depending on the
 * parity of the indices i and j, as well as on the overlaps between the parity
 * and update sets of the indices
 */
void seeley_richard_love(int i, int j, double coef, int n_qubits) {

    using double_complex = std::complex<double>;

    coef *= 0.25;

    cudaq::spin_op seeley_richard_love_result = 0.0 * cudaq::spin::i(0);
    // Case 0
    if (i == j) {  // Simplifies to the number operator
        cudaq::spin_op ops;
        // Create occupation set operator
        for (int index : occupation_set(i)) {
            ops *= cudaq::spin::z(index);
        }
        seeley_richard_love_result +=  coef * 2.0 * cudaq::spin_op();
        seeley_richard_love_result += -coef * 2.0 * ops;
    }

    // Case 1
    else if (i % 2 == 0 && j % 2 == 0) {
        // Create padding operators
        cudaq::spin_op x_pad;
        for (int index : U_diff_a_set(i, j, n_qubits)) {
            x_pad *= cudaq::spin::x(index);
        }
        cudaq::spin_op y_pad;
        for (int index : alpha_set(i, j, n_qubits)) {
            y_pad *= cudaq::spin::y(index);
        }
        cudaq::spin_op z_pad;
        for (int index : P0_ij_diff_a_set(i, j, n_qubits)) {
            z_pad *= cudaq::spin::z(index);
        }

        // Combine padding operators
        cudaq::spin_op left_pad = x_pad * y_pad * z_pad;

        // Create four operator combinations
        cudaq::spin_op op1 = left_pad * cudaq::spin::y(j) * cudaq::spin::x(i);
        cudaq::spin_op op2 = left_pad * cudaq::spin::x(j) * cudaq::spin::y(i);
        cudaq::spin_op op3 = left_pad * cudaq::spin::x(j) * cudaq::spin::x(i);
        cudaq::spin_op op4 = left_pad * cudaq::spin::y(j) * cudaq::spin::y(i);

        if (i < j) {
            seeley_richard_love_result += double_complex( coef, 0.0) * op1;
            seeley_richard_love_result += double_complex(-coef, 0.0) * op2;
            seeley_richard_love_result += double_complex(0.0, -coef) * op3;
            seeley_richard_love_result += double_complex(0.0, -coef) * op4;
        } else {  // whenever i < j introduce phase of -j
            seeley_richard_love_result += double_complex(0.0, -coef) * op1;
            seeley_richard_love_result += double_complex(0.0,  coef) * op2;
            seeley_richard_love_result += double_complex(-coef, 0.0) * op3;
            seeley_richard_love_result += double_complex(-coef, 0.0) * op4;
        }
    }

    // Case 2
    else if (i % 2 == 1 && j % 2 == 0 && _parity_set(j).contains(i)) {
        cudaq::spin_op x_pad;
        for (int index : U_diff_a_set(i, j, n_qubits)) {
            x_pad *= cudaq::spin::x(index);
        }
        cudaq::spin_op y_pad;
        for (int index : alpha_set(i, j, n_qubits)) {
            y_pad *= cudaq::spin::y(index);
        }

        // Combine padding operators
        cudaq::spin_op left_pad = x_pad * y_pad;

        cudaq::spin_op right_pad_1;
        std::unordered_set<std::size_t> P0_minus_alpha;
        std::set_difference(P0_ij_set(i, j).begin(), 
                            P0_ij_set(i, j).end(), 
                            alpha(i, j, n_qubits).begin(), 
                            alpha(i, j, n_qubits).end(), 
                            std::inserter(P0_minus_alpha, P0_minus_alpha.begin()));
        for (int index : P0_minus_alpha) {
            right_pad_1 *= cudaq::spin::z(index);
        }
        cudaq::spin_op right_pad_2;
        std::unordered_set<std::size_t> P2_minus_alpha;
        std::set_difference(P2_ij_set(i, j).begin(), 
                            P2_ij_set(i, j).end(), 
                            alpha(i, j, n_qubits).begin(), 
                            alpha(i, j, n_qubits).end(), 
                            std::inserter(P2_minus_alpha, P2_minus_alpha.begin()));
        for (int index : P2_minus_alpha) {
            right_pad_2 *= cudaq::spin::z(index);
        }

        cudaq:spin_op op1 = left_pad * cudaq::spin::y(j) * cudaq::spin::x(i) * right_pad1;
        cudaq:spin_op op2 = left_pad * cudaq::spin::x(j) * cudaq::spin::x(i) * right_pad1;
        cudaq:spin_op op3 = left_pad * cudaq::spin::x(j) * cudaq::spin::y(i) * right_pad2;
        cudaq:spin_op op4 = left_pad * cudaq::spin::y(j) * cudaq::spin::y(i) * right_pad2;

        if (i < j) {
        } else {
        }
    }
}

int main() {
    int n_qubits = 4;
    cudaq::spin_op op;
    cudaq::spin_op identity;
    op += cudaq::spin::x(4)*cudaq::spin::z(5);
    op -= identity;
    // seeley_richard_love(1, 1, 2.5, n_qubits);
    seeley_richard_love(2, 4, 2.5, n_qubits);
    return 0;
}

