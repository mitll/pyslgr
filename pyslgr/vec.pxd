#vec.pxd
cdef extern from "vec.hpp":
    cdef cppclass vec[T]:
        vec() except +
        vec(const vec &old_vec) except +
        vec(int l) except +    
        T *data
        int len
        
        
