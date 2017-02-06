#ifndef PYGMIN_BASE_INTERACTION_H
#define PYGMIN_BASE_INTERACTION_H


namespace pele {

struct BaseInteraction {
    const Array<double> m_radii;

    BaseInteraction()
        : m_radii(0)
    { }

    BaseInteraction(pele::Array<double> const & radii)
        : m_radii(radii.copy())
    { }
};

}

#endif // PYGMIN_BASE_INTERACTION_H
