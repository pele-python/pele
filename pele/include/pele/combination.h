#ifndef _COMBINATION_GENERATOR_H_
#define _COMBINATION_GENERATOR_H_

/*code adapted from http://coliru.stacked-crooked.com/view?id=c11dc445d3f7d49a415e3aa0478d7ba2-542192d2d8aca3c820c7acc656fa0c68
 * and here https://stackoverflow.com/questions/25138049/c-variable-number-of-nested-loops
 * This returns the combinations by obtaining the sorted permutations of the array elements
 * */
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace pele{

template<class iterator_type>
class combination_generator {
    iterator_type first, last;
    std::vector<bool> use;
    size_t r;
    typedef typename std::iterator_traits<iterator_type>::value_type element_type;
public:
    combination_generator(iterator_type first_, iterator_type last_, size_t r_)
    : first(first_),
      last(last_) ,
      r(r_)
    {
        use.resize(std::distance(first, last), false);
        if (r > use.size())
            throw std::domain_error("pele::combination_generator: can't select more elements than exist for combination");
        std::fill(use.end()-r, use.end(), true);
    }

    template<class output_iterator>
    bool operator()(output_iterator result)
    {
        iterator_type c=first;
        for (size_t i = 0; i<use.size(); ++i,++c) {
            if (use[i]){
                *result++ = *c;
            }
        }
        return std::next_permutation(use.begin(), use.end());
    }
};

template<class iterator_type>
combination_generator<iterator_type> make_combination_generator(iterator_type first, iterator_type last, size_t r)
{return combination_generator<iterator_type>(first, last, r);}

}

#endif
