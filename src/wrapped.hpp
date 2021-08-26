#pragma once

template <typename T, 
    typename = std::enable_if_t<std::is_integral<T>::value>
>
struct Wrapped
{
public:
    Wrapped(T start, T max):
        m_data(start), m_max(max)
        { }
    Wrapped() = delete;
    Wrapped& operator++()
    {
        m_data = (m_data + 1) % m_max;
        return *this;
    }
    Wrapped& operator++(int)
    {
        Wrapped<T> tmp(m_data, m_max);
        ++(*this);
        return tmp;
    }
    Wrapped& operator=(const Wrapped& rhs)
    {
        m_data = rhs.m_data;
        return *this;
    }
    Wrapped& operator=(const T& rhs)
    {
        m_data = rhs;
        return *this;
    }
    T& GetValue()
    {
        return m_data;
    }
    const T& GetValue() const
    {
        return m_data;
    }
private:
    T m_data = 0.0;
    const T m_max;
};