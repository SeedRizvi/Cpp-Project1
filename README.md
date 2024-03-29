## University Project
The following was an assignment completed as a part of university project for the course "Computing for Mechatronic Engineers". This task was my first major experience with C++ and taught me how to write programs based on extensive project specifications as well as learning a wide range of libraries provided in C++. 

The outline of the task can be seen below.
## Task

A `SmallMatrix` is a small-storage-optimised matrix whose elements are allocated on the stack if the number of elements is less than 144 allowing fast read/write speeds. If the number of elements is 144 or greater, then its contents are allocated on the heap. Accessing heap memory is slower than accessing stack memory due to the extra level of indirection of accessing heap-allocated data by first looking up its address on the stack.

When the size of the `SmallMatrix` increases above its small-size threshold, it is able to switch from being stack-allocated to heap-allocated. However, when the size of the `SmallMatrix` decreases below its small-size threshold, it **remains** heap-allocated as there are no real-life performance benefits to be gained from switching back-and-forth between memory spaces.

The task is to implement a `SmallMatrix` class as required by the given [specification](#specification). The interface and behaviour of `SmallMatrix` is the same regardless if it is stack-allocated or heap-allocated.

The `SmallMatrix` interface is provided in `SmallMatrix.hpp`. You are also given `SmallMatrix.cpp` for implementation and `test_small_matrix.cpp` for a simple testing suite. You may make modifications to each of these files as you see fit and according to the specification.

## Specification
The specification for `SmallMatrix` is summarised below:

<table>
    <tr>
        <th>Method</th>
        <th>Description</th>
        <th>Usage</th>
        <th>Exceptions</th>
    </tr>
    <tr>
        <td><code>SmallMatrix()</code></td>
        <td>A constructor which initialises an empty matrix with no rows and no columns. </td>
        <td><pre><code>SmallMatrix m;</code></pre></td>
        <td>None</td>
    </tr>
    <tr>
        <td><code>SmallMatrix(int, int)</code></td>
        <td>A constructor which initialises a zero matrix with the dimensions given by <code>mNumRows</code> and <code>mNumCols</code>.</td>
        <td><pre><code>SmallMatrix m(7, 4);</code></pre></td>
        <td>None</td>
    </tr>
    <tr>
        <td><code>SmallMatrix(int, int, double)</code></td>
        <td>A constructor which intialises a matrix whose elements are all initialised with the given value, and has the dimensions given by <code>mNumRows</code> and <code>mNumCols</code>. </td>
        <td><pre><code>SmallMatrix m(7, 4, 42.2);</code></pre></td>
        <td>None</td>
    </tr>
    <tr>
        <td><code>SmallMatrix(std::initializer_list&lt;std::initializer_list&lt;double&gt;&gt; const&)</code></td>
        <td><s>A constructor which initialises a matrix with a given initialiser list of initialiser list of doubles i.e. a 2D initialiser list of doubles. Each inner initialiser list represents a single row where each element in the inner initialiser list represents a column.</s> <b>GIVEN</b></td>
        <td><pre><code>SmallMatrix m({
    {1.0, 2.0, 3.0, 4.0},
    {5.0, 6.0, 7.0, 8.0},
});</code></pre></td>
        <td>Throws <code>invalid_argument</code> if the initialiser list is not rectangular i.e. each row does not have the same number of columns.</td>
    </tr>
    <tr>
        <td><code>SmallMatrix(SmallMatrix const&)</code></td>
        <td>Copy constructor.</td>
        <td><pre><code>SmallMatrix m1;
SmallMatrix m2(m1);</code></pre></td>
        <td>None</td>
    </tr>
    <tr>
        <td><code>SmallMatrix(SmallMatrix&&)</code></td>
        <td>Move constructor. Specified object should be invalidated after move.</td>
        <td><pre><code>SmallMatrix m1;
SmallMatrix m2(std::move(m1));</code></pre></td>
        <td>None</td>
    </tr>
    <tr>
        <td><code>SmallMatrix& operator=(SmallMatrix const&)</code></td>
        <td>Copy assignment.</td>
        <td><pre><code>SmallMatrix m1;
SmallMatrix m2;
m2 = m1;</code></pre></td>
        <td>None</td>
    </tr>
    <tr>
        <td><code>SmallMatrix& operator=(SmallMatrix&&)</code></td>
        <td>Move assignment. Specified object should be invalidated after move.</td>
        <td><pre><code>SmallMatrix m1;
SmallMatrix m2;
m2 = std::move(m1);</code></pre></td>
        <td>None</td>
    </tr>
    <tr>
        <td><code>~SmallMatrix()</code></td>
        <td>Destructor.</td>
        <td></td>
        <td>None</td>
    </tr>
    <tr>
        <td><code>double& operator()(int, int)</code></td>
        <td>Returns the reference of the matrix element at the specified row and column index. The order of access is: <code>(row, col)</code></td>
        <td><pre><code>SmallMatrix m(1, 1);
m(0, 0) = 24.4;</pre></code></td>
        <td>Throws <code>out_of_range</code> if the specified row and column is outside the range <code>[0, max_row)</code> and <code>[0, max_col)</code> respectively.<br><br>
        Throws <code>out_of_range</code> if the matrix has no rows and no columns.</td>
    </tr>
    <tr>
        <td><code>const double& operator()(int, int) const</code></td>
        <td>Returns the constant reference of the matrix element at the specified row and column index. It is guaranteed that the returned element is not modified.</td>
        <td><pre><code>SmallMatrix m(1, 1);
m(0, 0);</pre></code></td>
        <td>Throws <code>out_of_range</code> if the specified row and column is outside the range <code>[0, max_row)</code> and <code>[0, max_col)</code> respectively.<br><br>
        Throws <code>out_of_range</code> if the matrix has no rows and no columns.</td>
    </tr>
    <tr>
        <td><code>std::vector&lt;double*&gt; row(int)</code></td>
        <td>Returns a vector of pointers to each of the elements of the row of the matrix at the specified row index.</td>
        <td><pre><code>SmallMatrix m(1, 1);
auto r = m.row(0);
r[0] = 2.2;</pre></code></td>
        <td>Throws <code>out_of_range</code> if the specified row index is outside the range <code>[0, max_row)</code>.</td>
    </tr>
    <tr>
        <td><code>std::vector&lt;double const*&gt; row(int) const</code></td>
        <td>Returns a vector of pointers to each of the elements of constant type of the row of the matrix at the specified row index.</td>
        <td><pre><code>SmallMatrix m(1, 1);
m.row(0);</pre></code></td>
        <td>Throws <code>out_of_range</code> if the specified row index is outside the range <code>[0, max_row)</code>.</td>
    </tr>
    <tr>
        <td><code>std::vector&lt;double*&gt; col(int)</code></td>
        <td>Returns a vector of pointers to each of the elements of the column of the matrix at the specified column index.</td>
        <td><pre><code>SmallMatrix m(1, 1);
auto c = m.col(0);
r[0] = 2.2;</pre></code></td>
        <td>Throws <code>out_of_range</code> if the specified column index is outside the range <code>[0, max_col)</code>.</td>
    </tr>
    <tr>
        <td><code>std::vector&lt;double const*&gt; col(int) const</code></td>
        <td>Returns a vector of pointers to each of the elements of constant type of the column of the matrix at the specified column index.</td>
        <td><pre><code>SmallMatrix m(1, 1);
m.col(0);</pre></code></td>
        <td>Throws <code>out_of_range</code> if the specified column index is outside the range <code>[0, max_col)</code>.</td>
    </tr>
    <tr>
        <td><code>std::pair&lt;int, int&gt; size() const</code></td>
        <td>Returns the size of the matrix where the first of the pair is the number of rows and the second of the pair is the number of columns.</td>
        <td><pre><code>SmallMatrix m(1, 1);
auto s = m.size();
s.first;
s.second;</pre></code></td>
        <td>None</td>
    </tr>
    <tr>
        <td><code>bool isSmall() const</code></td>
        <td>Returns true if the matrix is using a small-storage-optimised data structure.</td>
        <td><pre><code>SmallMatrix m(1, 1);
s.isSmall();</pre></code></td>
        <td>None</td>
    </tr>
    <tr>
        <td><code>void resize(int, int)</code></td>
        <td>Resizes the matrix to the new number of rows and new number of columns. If any matrix dimension is increased, then the newly created dimension is zero-initialised. If any matrix dimension is decreased, then its previously-allocated elements are truncated.</td>
        <td><pre><code>SmallMatrix m(1, 1);
s.resize(100, 100);</pre></code></td>
        <td>Throws <code>out_of_range</code> if the specified row or column index is negative.</td>
    </tr>
    <tr>
        <td><code>void insertRow(int, std::vector<double> const&)</code></td>
        <td>Inserts a row at the specified row index. If the number of columns in the matrix is zero, then the matrix is resized to match the size of the specified row vector.</td>
        <td><pre><code>SmallMatrix m(2, 4);
s.insertRow(0, {1, 2, 3, 4});</pre></code></td>
        <td>Throws <code>out_of_range</code> if the specified row index is outside the range <code>[0, max_row]</code>.<br><br>Throws <code>invalid_argument</code> if the size of the specified vector is not equal to the number of columns in the matrix.</td>
    </tr>
    <tr>
        <td><code>void insertCol(int, std::vector<double> const&)</code></td>
        <td>Inserts a column at the specified column index. If the number of rows in the matrix is zero, then the matrix is resized to match the size of the specified column vector.</td>
        <td><pre><code>SmallMatrix m(3, 2);
s.insertCol(2, {1, 2, 3});</pre></code></td>
        <td>Throws <code>out_of_range</code> if the specified column index is outside the range <code>[0, max_col]</code>.<br><br>Throws <code>invalid_argument</code> if the size of the specified vector is not equal to the number of rows in the matrix.</td>
    </tr>
    <tr>
        <td><code>void eraseRow(int)</code></td>
        <td>Erases the row at the specified row index.</td>
        <td><pre><code>SmallMatrix m(3, 2);
s.eraseRow(2);</pre></code></td>
        <td>Throws <code>out_of_range</code> if the specified row index is outside the range <code>[0, max_row)</code>.</td>
    </tr>
    <tr>
        <td><code>void eraseCol(int)</code></td>
        <td>Erases the column at the specified column index.</td>
        <td><pre><code>SmallMatrix m(3, 2);
s.eraseCol(1);</pre></code></td>
        <td>Throws <code>out_of_range</code> if the specified column index is outside the range <code>[0, max_col)</code>.</td>
    </tr>
    <tr>
        <td><code>friend bool operator==(SmallMatrix const&, SmallMatrix const&)</code></td>
        <td>Returns true if all of the elements in the left-hand side matrix are equal to its positionally-corresponding element in the right-hand side matrix. Otherwise, false.</td>
        <td><pre><code>SmallMatrix m1({{1, 2, 3}, {4, 5, 6}});
SmallMatrix m2({{1, 2, 3}, {4, 5, 6}});
m1 == m2;</pre></code></td>
        <td>None.</td>
    </tr>
    <tr>
        <td><code>friend bool operator!=(SmallMatrix const&, SmallMatrix const&)</code></td>
        <td>Returns true if any of the elements in the left-hand side matrix are not equal to its positionally-corresponding element in the right-hand side matrix. Otherwise, false.</td>
        <td><pre><code>SmallMatrix m1({{1, 2, 3}, {4, 5, 7}});
SmallMatrix m2({{1, 2, 3}, {4, 5, 6}});
m1 != m2;</pre></code></td>
        <td>None.</td>
    </tr>
    <tr>
        <td><code>friend SmallMatrix operator+(SmallMatrix const&, SmallMatrix const&)</code></td>
        <td>Returns the matrix result of the element-wise addition of the two specified matrices.</td>
        <td><pre><code>SmallMatrix m1({{1, 2}, {3, 4}, {5, 6}});
SmallMatrix m2({{1, 2}, {3, 4}, {5, 6}});
auto r = m1 + m2;</pre></code></td>
        <td>Throws <code>invalid_argument</code> if the number of rows and columns on the left-hand side is not equal to the number of rows and columns on the right-hand side respectively.</td>
    </tr>
    <tr>
        <td><code>friend SmallMatrix operator-(SmallMatrix const&, SmallMatrix const&)</code></td>
        <td>Returns the matrix result of the element-wise subtraction of the two specified matrices.</td>
        <td><pre><code>SmallMatrix m1({{1, 2}, {3, 4}, {5, 6}});
SmallMatrix m2({{1, 2}, {3, 4}, {5, 6}});
auto r = m1 - m2;</pre></code></td>
        <td>Throws <code>invalid_argument</code> if the number of rows and columns on the left-hand side is not equal to the number of rows and columns on the right-hand side respectively.</td>
    </tr>
    <tr>
        <td><code>friend SmallMatrix operator*(SmallMatrix const&, SmallMatrix const&)</code></td>
        <td>Returns the matrix result of the matrix multiplication of the two specified matrices.</td>
        <td><pre><code>SmallMatrix m1({{1, 2}, {3, 4}, {5, 6}});
SmallMatrix m2({{1, 2}, {3, 4}});
auto r = m1 * m2;</pre></code></td>
        <td>Throws <code>invalid_argument</code> if the number of columns on the left-hand side is not equal to the number of rows on the right-hand side.</td>
    </tr>
    <tr>
        <td><code>friend SmallMatrix operator*(double, SmallMatrix const&)</code></td>
        <td>Returns the matrix result of the scalar multiplication of the the specified scalar value and specified matrix.</td>
        <td><pre><code>SmallMatrix m({{1, 2}, {3, 4}, {5, 6}});
auto r = 42.2 * m;</pre></code></td>
        <td>None.</td>
    </tr>
    <tr>
        <td><code>friend SmallMatrix operator*(SmallMatrix const&, double)</code></td>
        <td>Returns the matrix result of the scalar multiplication of the the specified scalar value and specified matrix.</td>
        <td><pre><code>SmallMatrix m({{1, 2}, {3, 4}, {5, 6}});
auto r = m * 42.2;</pre></code></td>
        <td>None.</td>
    </tr>
    <tr>
        <td><code>SmallMatrix& operator+=(SmallMatrix const&)</code></td>
        <td>Returns *this after the element-wise addition of *this and the specified matrix. This operation is equivalent to <code>*this = *this + m</code>.</td>
        <td><pre><code>SmallMatrix m1({{1, 2}, {3, 4}, {5, 6}});
SmallMatrix m2({{1, 1}, {1, 1}, {, 1}});
auto m1 += m2;</pre></code></td>
        <td>Throws <code>invalid_argument</code> if the number of rows and columns of *this is not equal to the number of rows and columns of the specified matrix respectively.</td>
    </tr>
    <tr>
        <td><code>SmallMatrix& operator-=(SmallMatrix const&)</code></td>
        <td>Returns *this after the element-wise subtraction of *this and the specified matrix. This operation is equivalent to <code>*this = *this - m</code>.</td>
        <td><pre><code>SmallMatrix m1({{1, 2}, {3, 4}, {5, 6}});
SmallMatrix m2({{1, 1}, {1, 1}, {, 1}});
auto m1 -= m2;</pre></code></td>
        <td>Throws <code>invalid_argument</code> if the number of rows and columns of *this is not equal to the number of rows and columns of the specified matrix respectively.</td>
    </tr>
    <tr>
        <td><code>SmallMatrix& operator*=(SmallMatrix const&)</code></td>
        <td>Returns *this after the matrix multiplication of *this and the specified matrix. This operation is equivalent to <code>*this = *this * m</code>.</td>
        <td><pre><code>SmallMatrix m1({{1, 2}, {3, 4}, {5, 6}});
SmallMatrix m2({{1, 2}, {3, 4}});
auto m1 *= m2;</pre></code></td>
        <td>Throws <code>invalid_argument</code> if the number of columns of *this is not equal to the number of rows of the specified matrix.</td>
    </tr>
    <tr>
        <td><code>SmallMatrix& operator*=(double)</code></td>
        <td>Returns *this after the scalar multiplication of *this and the specified scalar value. This operation is equivalent to <code>*this = *this * s</code>.</td>
        <td><pre><code>SmallMatrix m({{1, 2}, {3, 4}, {5, 6}});
auto m1 *= 42.2;</pre></code></td>
        <td>None.</td>
    </tr>
    <tr>
        <td><code>friend SmallMatrix transpose(SmallMatrix const&)</code></td>
        <td>Returns the result of the tranpose on the specified matrix.</td>
        <td><pre><code>SmallMatrix m({{1, 2, 3}, {4, 5, 6}});
auto r = transpose(m);</pre></code></td>
        <td>None.</td>
    </tr>
    <tr>
        <td><code>friend std::ostream& operator&lt;&lt;(std::ostream&, SmallMatrix const&)</code></td>
        <td>Writes the contents of the matrix to the output stream.</td>
        <td><pre><code>SmallMatrix m({{1, 2, 3}, {4, 5, 6}});
std::cout &lt;&lt; m;</pre></code></td>
        <td>None.</td>
    </tr>
</table>
