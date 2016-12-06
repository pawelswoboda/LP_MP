#ifndef LP_MP_STATIC_IF_HXX
#define LP_MP_STATIC_IF_HXX

// static if implementation taken from boost mailing list: lists.boost.org/Archives/boost/2014/08/216607.php

namespace LP_MP {

namespace aux {
    struct identity
    {
        template<class T>
        T operator()(T&& x) const
        {
            return std::forward<T>(x);
        }
    };

    template< bool Condition >
    struct static_if_statement
    {
        template< typename F > void then ( F const& f ) { f(identity()); }
        template< typename F > void else_ ( F const& ) { }
    };

    template< >
    struct static_if_statement<false>
    {
        template< typename F > void then ( F const& ) { }
        template< typename F > void else_ ( F const& f ) { f(identity()); }
    };
}

template< bool Condition, typename F >
aux::static_if_statement<Condition> static_if ( F const& f )
{
   aux::static_if_statement<Condition> if_;
   if_.then(f);
   return if_;
} 

} // end namespace LP_MP

#endif // LP_MP_STATIC_IF_HXX
