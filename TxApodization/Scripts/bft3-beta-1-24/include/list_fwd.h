#ifndef LIST_FWD_H
#define LIST_FWD_H
// std::list forward declaration
namespace std
{
     template<typename T> class allocator;
     template<typename T, typename A = allocator<T> > class list;
}
#endif
