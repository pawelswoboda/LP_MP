#ifndef LP_MP_SPINLOCK_HXX
#define LP_MP_SPINLOCK_HXX 

#include <atomic>

#if defined(_MSC_VER) && _MSC_VER >= 1310 && ( defined(_M_IX86) || defined(_M_X64) )

extern "C" void _mm_pause();

#define LP_MP_SMT_PAUSE _mm_pause();

#elif defined(__GNUC__) && ( defined(__i386__) || defined(__x86_64__) )

#define LP_MP_PAUSE __asm__ __volatile__( "rep; nop" : : : "memory" );

#endif


namespace LP_MP { 

  class spinlock {
    std::atomic_flag locked = ATOMIC_FLAG_INIT ;
    public:
    void lock() {
      while (locked.test_and_set(std::memory_order_acquire)) { 
        LP_MP_PAUSE
      }
  }
  void unlock() {
    locked.clear(std::memory_order_release);
  }
};

} // end namespace LP_MP

#endif // LP_MP_SPINLOCK_HXX
