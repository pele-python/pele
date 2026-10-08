#ifndef CONSTEXPR_FOR
#define CONSTEXPR_FOR

#include <cstddef>
#include <utility>

template <std::size_t N>
struct num {
  static const constexpr auto value = N;
};

template <class F, std::size_t... Is>
void for_(F func, std::index_sequence<Is...>) {
  (func(num<Is>{}), ...);
}

template <std::size_t N, typename F>
void for_(F func) {
  for_(func, std::make_index_sequence<N>());
}

#endif  // !CONSTEXPR_FOR