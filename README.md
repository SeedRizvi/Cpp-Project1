## University Project
The following was an assignment completed as a part of university project for the course "Computing for Mechatronic Engineers". This task was my first major experience with C++ and taught me how to write programs based on extensive project specifications as well as learning a wide range of libraries provided in C++. 

The outline of the task can be seen below.
## Task

A `SmallMatrix` is a small-storage-optimised matrix whose elements are allocated on the stack if the number of elements is less than 144 allowing fast read/write speeds. If the number of elements is 144 or greater, then its contents are allocated on the heap. Accessing heap memory is slower than accessing stack memory due to the extra level of indirection of accessing heap-allocated data by first looking up its address on the stack.

When the size of the `SmallMatrix` increases above its small-size threshold, it is able to switch from being stack-allocated to heap-allocated. However, when the size of the `SmallMatrix` decreases below its small-size threshold, it **remains** heap-allocated as there are no real-life performance benefits to be gained from switching back-and-forth between memory spaces.

The task is to implement a `SmallMatrix` class as required by the given [specification](#specification). The interface and behaviour of `SmallMatrix` is the same regardless if it is stack-allocated or heap-allocated.

The `SmallMatrix` interface is provided in `SmallMatrix.hpp`. You are also given `SmallMatrix.cpp` for implementation and `test_small_matrix.cpp` for a simple testing suite. You may make modifications to each of these files as you see fit and according to the specification.

