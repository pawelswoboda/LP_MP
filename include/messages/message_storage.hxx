#ifndef LP_MP_MESSAGE_STORAGE_HXX
#define LP_MP_MESSAGE_STORAGE_HXX

#include "LP_MP.h"

// various ways to store messages

namespace LP_MP {

class StandardMessageStorage : public std::vector<REAL>
{
public:
   StandardMessageStorage(const INDEX msg_size) : std::vector<REAL>(msg_size,0.0) {}
}; 

template<INDEX N>
class FixedMessageStorage : public std::array<REAL,N>
{
public:
   FixedMessageStorage(const INDEX msg_size) { std::array<REAL,N>::fill(0.0); }
}; 

class MessageStorageSIMD : public Vc::Memory<REAL_SIMD>
{ 
public:
   using Memory = Vc::Memory<REAL_SIMD>;
   MessageStorageSIMD(const INDEX msg_size) : Vc::Memory<REAL_SIMD>(msg_size)
   {
      for(INDEX i=0; i<Memory::vectorsCount(); ++i) {
         Memory::vector(i) = 0.0; // Vc::Zero
      }
   }  
   const size_t size() const { return Memory::entriesCount(); }
   REAL& operator[](const INDEX i) { return Memory::scalar(i); }
   const REAL operator[](const INDEX i) const { return Memory::scalar(i); }
}; 


// if left and right factor hold the reparametrizations explicitly, then there is no need to record the messages explicitly
// do zrobienia: not supported yet
class EmptyMessageStorage
{
public:
   EmptyMessageStorage(const INDEX msg_size) {}
};

} // end namespace LP_MP

#endif // LP_MP_MESSAGE_STORAGE_HXX
