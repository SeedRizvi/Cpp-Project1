    #include "SmallMatrix.hpp"

    namespace mtrn2500 {

    SmallMatrix::SmallMatrix() 
        : mNumRows{0},
        mNumCols{0},
        mIsLargeMatrix{mNumRows * mNumCols >= mSmallSize} {
        // Make it so user cannot access any 144x144 array rows or cols.
        
    }

    SmallMatrix::SmallMatrix(int numRows, int numCols) 
        : mNumRows(numRows),
        mNumCols(numCols),
        mIsLargeMatrix(mNumRows * mNumCols >= mSmallSize) {

        if (mIsLargeMatrix) {
            mHeapData.resize(mNumRows); 
        }
        fillMatrix();
    }

    SmallMatrix::SmallMatrix(int numRows, int numCols, double value)
        : mNumRows(numRows),
        mNumCols(numCols),
        mIsLargeMatrix(mNumRows * mNumCols >= mSmallSize) {

        if (mIsLargeMatrix) {
            mHeapData.resize(mNumRows); 
        }
        fillMatrix(value);
    }

    SmallMatrix::SmallMatrix(std::initializer_list<std::initializer_list<double>> const& il)
        : mNumRows(il.size()),
        mNumCols(il.begin() == il.end() ? 0 : il.begin()->size()),
        mIsLargeMatrix(mNumRows * mNumCols >= mSmallSize) {
        if (std::adjacent_find(il.begin(), il.end(), [](auto const& lhs, auto const& rhs) {
                return lhs.size() != rhs.size();
            }) != il.end()) {
            throw std::invalid_argument("Rows have different sizes.");
        }

        if (mIsLargeMatrix) {
            mHeapData.resize(mNumRows);
        }

        int row_index{0};
        for (auto const& row : il) {
            if (mIsLargeMatrix) {
                mHeapData.at(row_index).resize(mNumCols);
                std::copy(row.begin(), row.end(), mHeapData.at(row_index).begin());
            } else {
                std::transform(row.begin(), row.end(), mStackData.at(row_index).begin(),
                            [](auto const& e) { return e; });
            }
            row_index++;
        }
    }

    SmallMatrix::SmallMatrix(SmallMatrix const& sm) { // Copy Constructor
        mNumRows = sm.mNumRows;
        mNumCols = sm.mNumCols;
        mIsLargeMatrix = sm.mIsLargeMatrix;
        mStackData = sm.mStackData;
        mHeapData  =  sm.mHeapData;
    }

    SmallMatrix::SmallMatrix(SmallMatrix&& sm) { // Move Constructor
        if (sm.mIsLargeMatrix) {
            mHeapData = sm.mHeapData;
        } else {
            mStackData = sm.mStackData;
        }
        // Invalidating original object
        sm.fillMatrix();
        mNumRows = std::exchange(sm.mNumRows, 0);
        mNumCols = std::exchange(sm.mNumCols, 0);
        mIsLargeMatrix = std::exchange(sm.mIsLargeMatrix, false);
    } // P  

    SmallMatrix& SmallMatrix::operator=(SmallMatrix const& sm) { // Copy assignment
        if (this !=  &sm) {
            //assign rows, cols, mStackData, mHeapData
            mNumRows = sm.mNumRows;
            mNumCols = sm.mNumCols;
            mIsLargeMatrix = sm.mIsLargeMatrix;
            mStackData = sm.mStackData;
            mHeapData  =  sm.mHeapData;
        }   
        return *this; 
    } // P  

    SmallMatrix& SmallMatrix::operator=(SmallMatrix&& sm) { // Move assignment
        if (this != &sm) {
            if (sm.mIsLargeMatrix) {
                mHeapData = sm.mHeapData;
            } else {
                mStackData = sm.mStackData;
            }
            // Invalidating original object
            sm.fillMatrix();
            mNumRows = std::exchange(sm.mNumRows, 0);
            mNumCols = std::exchange(sm.mNumCols, 0);
            mIsLargeMatrix = std::exchange(sm.mIsLargeMatrix, false);
        }
        
        return *this; 
    } // P  

    SmallMatrix::~SmallMatrix() {} // Nothing specified due to no dynamic memory allocation.

    double& SmallMatrix::operator()(int numRow, int numCol) {
        if (mNumRows == 0 or mNumCols == 0) {
            throw std::out_of_range("This matrix has no rows or columns.\n");

        } else if ((numRow >= mNumRows or numRow < 0) or (numCol >= mNumCols or numCol < 0)) {
            throw std::out_of_range("Specified row or column is out of bounds.\n");  

        } else {
            if (mIsLargeMatrix) {
                return mHeapData.at(numRow).at(numCol);
            } else {
                return mStackData.at(numRow).at(numCol);
            }
        }   
    } // P  

    const double& SmallMatrix::operator()(int numRow, int numCol) const {
        if (mNumRows == 0 or mNumCols == 0) {
            throw std::out_of_range("This matrix has no rows or columns.\n");

        } else if ((numRow >= mNumRows or numRow < 0) or (numCol >= mNumCols or numCol < 0)) {
            throw std::out_of_range("Specified row or column is out of bounds.\n");  

        } else {
            if (mIsLargeMatrix) {
                return mHeapData.at(numRow).at(numCol);
            } else {
                return mStackData.at(numRow).at(numCol);
            }
        }
    } // P  

    std::vector<double*> SmallMatrix::row(int numRow) {
        if ((numRow >= mNumRows or numRow < 0)) {
            throw std::out_of_range("Specified column is out of bounds.\n"); 
        }
        std::vector<double*> result_row;
        for (int col_index{0}; col_index < mNumCols; col_index++) {
            if (mIsLargeMatrix) {
                result_row.push_back(&mHeapData.at(numRow).at(col_index));
            } else {
                result_row.push_back(&mStackData.at(numRow).at(col_index));
            }
            
        }
        return result_row;
    } // CR

    std::vector<double const*> SmallMatrix::row(int numRow) const {
        if ((numRow >= mNumRows or numRow < 0)) {
            throw std::out_of_range("Specified column is out of bounds.\n"); 
        }
        std::vector<double const*> result_row;
        for (int col_index{0}; col_index < mNumCols; col_index++) {
            if (mIsLargeMatrix) {
                result_row.push_back(&mHeapData.at(numRow).at(col_index));
            } else {
                result_row.push_back(&mStackData.at(numRow).at(col_index));
            }
        }
        return result_row;
    } // CR

    std::vector<double*> SmallMatrix::col(int numCol) {
        if ((numCol >= mNumCols or numCol < 0)) {
            throw std::out_of_range("Specified column is out of bounds.\n"); 
        }
        std::vector<double*> result_col;
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            if (mIsLargeMatrix) {
                result_col.push_back(&mHeapData.at(row_index).at(numCol));
            } else {
                result_col.push_back(&mStackData.at(row_index).at(numCol));
            }
        }
        return result_col;
    } // CR

    std::vector<double const*> SmallMatrix::col(int numCol) const {
        if ((numCol >= mNumCols or numCol < 0)) {
            throw std::out_of_range("Specified column is out of bounds.\n"); 
        }
        std::vector<double const*> result_col;
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            if (mIsLargeMatrix) {
                result_col.push_back(&mHeapData.at(row_index).at(numCol));
            } else {
                result_col.push_back(&mStackData.at(row_index).at(numCol));
            }
        }
        return result_col;
    } // CR

    std::pair<int, int> SmallMatrix::size() const { // P  
        std::pair<int, int> Size(mNumRows, mNumCols);
        return Size;
    } 

    bool SmallMatrix::isSmall() const { // P  
        return !mIsLargeMatrix;
    }

    void SmallMatrix::resize(int numRows, int numCols) {
        if ((numRows < 0) or (numCols < 0)) {
            throw std::out_of_range("Resize rows or columns cannot be negative.");
        } else if ((numRows == 0) and (numCols == 0)) {
            // auto m1 = SmallMatrix();
            // *this = std::move(m1);
            *this = std::move(SmallMatrix());
        }
        else if (mIsLargeMatrix) {
            if (numRows < mNumRows) {
                // mHeapData.erase(mHeapData.begin() + numRows, mHeapData.begin() + mNumRows);
                mHeapData.erase(mHeapData.begin() + numRows, mHeapData.end());
                mNumRows = numRows;
            } 
            if (numCols < mNumCols) {
                for (int row_index{0}; row_index < mNumRows; row_index++) {
                    mHeapData.at(row_index).erase(mHeapData.at(row_index).begin() 
                        + numCols, mHeapData.at(row_index).end());
                }
                mNumCols = numCols;
            }
            if (numRows > mNumRows) {
                mHeapData.resize(numRows);
                fillRow(mNumRows, numRows);
                mNumRows = numRows;
            }
            if (numCols > mNumCols) {
                // do something
                fillCol(mNumCols, numCols);
                mNumCols = numCols;
            }
        }
        else { // Make sure this only runs if all other if's are not true
            if (numRows * numCols >= mSmallSize) {
                //Convert from Stack to heapdata if new size is Large.
                std::cout << "Calling stack to heap conversion\n";
                stackToHeap(numRows, numCols);
                mNumRows = numRows;
                mNumCols = numCols;
                mIsLargeMatrix = true;
            } else {
                if (numRows < mNumRows) {
                    fillRow(numRows, mNumRows);
                    mNumRows = numRows;
                }
                if (numCols < mNumCols) {
                    fillCol(numCols, mNumCols);
                    mNumCols = numCols;
                }
                if (numRows > mNumRows) {
                    fillRow(mNumRows, numRows);
                    mNumRows = numRows;
                }
                if (numCols > mNumCols) {
                    fillCol(mNumCols, numCols);
                    mNumCols = numCols;
                }
            }
        }
    } // DN

    void SmallMatrix::insertRow(int numRow, std::vector<double> const& row) {} // HD

    void SmallMatrix::insertCol(int numCol, std::vector<double> const& col) {} // HD

    void SmallMatrix::eraseRow(int numRow) {} // HD

    void SmallMatrix::eraseCol(int numCol) {} // HD

    bool operator==(SmallMatrix const& lhs, SmallMatrix const& rhs) {
        if ((lhs.mNumRows != rhs.mNumRows) or (lhs.mNumCols != rhs.mNumCols)) {
            return false;
        } 
        for (int row_index{0}; row_index < lhs.mNumRows; row_index++) {
            if (lhs.mIsLargeMatrix and rhs.mIsLargeMatrix) {
                if (lhs.mHeapData.at(row_index) != rhs.mHeapData.at(row_index)) {
                    return false;
                }
            } else {
                if (lhs.mStackData.at(row_index) != rhs.mStackData.at(row_index)) {
                    return false;
                }
            }
        }
        return true;
    } // P  

    bool operator!=(SmallMatrix const& lhs, SmallMatrix const& rhs) {
        if ((lhs.mNumRows != rhs.mNumRows) or (lhs.mNumCols != rhs.mNumCols)) {
            return true;
        } 
        for (int row_index{0}; row_index < lhs.mNumRows; row_index++) {
            if (lhs.mIsLargeMatrix and rhs.mIsLargeMatrix) {
                if (lhs.mHeapData.at(row_index) != rhs.mHeapData.at(row_index)) {
                    return true;
                }
            } else {
                if (lhs.mStackData.at(row_index) != rhs.mStackData.at(row_index)) {
                    return true;
                }
            }
        }
        return false;
    } // P  

    SmallMatrix operator+(SmallMatrix const& lhs, SmallMatrix const& rhs) {
        if ((lhs.mNumRows != rhs.mNumRows) or (lhs.mNumCols != rhs.mNumCols)) {
            throw std::invalid_argument("Provided matrices are not of equal size");
        } 
        SmallMatrix result(lhs.mNumRows, lhs.mNumCols);
        result.matrixAddition(lhs, rhs);
        return result;
    } // CR

    SmallMatrix operator-(SmallMatrix const& lhs, SmallMatrix const& rhs) {
        if ((lhs.mNumRows != rhs.mNumRows) or (lhs.mNumCols != rhs.mNumCols)) {
            throw std::invalid_argument("Provided matrices are not of equal size");
        } 
        SmallMatrix result(lhs.mNumRows, lhs.mNumCols);
        result.matrixSubtraction(lhs, rhs);
        return result;
    } // CR

    SmallMatrix operator*(SmallMatrix const& lhs, SmallMatrix const& rhs) {
        if (lhs.mNumCols != rhs.mNumRows) {
            throw std::invalid_argument("LHS num columns and RHS num rows do not match.");
        }
        SmallMatrix result(lhs.mNumRows, rhs.mNumCols);
        result.matrixMultiply(lhs, rhs);
        return result;
    } // DN

    SmallMatrix operator*(double s, SmallMatrix const& sm) {
        SmallMatrix result(sm.mNumRows, sm.mNumCols);
        result.scalarMultiply(s, sm);
        return result;
    } // CR

    SmallMatrix operator*(SmallMatrix const& sm, double s) {
        SmallMatrix result(sm.mNumRows, sm.mNumCols);
        result.scalarMultiply(s, sm);
        return result;
    } // CR

    SmallMatrix& SmallMatrix::operator+=(SmallMatrix const& sm) {
        if ((this->mNumRows != sm.mNumRows) or (this->mNumCols != sm.mNumCols)) {
            throw std::invalid_argument("Provided matrices are not of equal size");
        } 
        this->matrixAddition(sm);
        return *this;
    } // CR

    SmallMatrix& SmallMatrix::operator-=(SmallMatrix const& sm) {
        if ((this->mNumRows != sm.mNumRows) or (this->mNumCols != sm.mNumCols)) {
            throw std::invalid_argument("Provided matrices are not of equal size");
        } 
        this->matrixSubtraction(sm);
        return *this;
    } // CR

    SmallMatrix& SmallMatrix::operator*=(SmallMatrix const& sm) {
        if (this->mNumCols != sm.mNumRows) {
            throw std::invalid_argument("LHS num columns and RHS num rows do not match.");
        }
        this->matrixMultiply(sm);
        return *this;
    } // DN

    SmallMatrix& SmallMatrix::operator*=(double s) {
        this->scalarMultiply(s);
        return *this;
    } // CR

    SmallMatrix transpose(SmallMatrix const& sm) { return {}; } // DN

    std::ostream& operator<<(std::ostream& os, SmallMatrix const& sm) {
        os << "[\n";
        for (int row_index{0}; row_index < sm.mNumRows; row_index++) {
            if (sm.mIsLargeMatrix) {
                os << "  [ ";
                std::copy(sm.mHeapData.at(row_index).begin(),
                    sm.mHeapData.at(row_index).end(), 
                    std::ostream_iterator<double>(os," "));
                os << "]\n";
            } else {
                auto first = sm.mStackData.at(row_index).begin();
                auto last = sm.mStackData.at(row_index).begin() + sm.mNumCols;
                os << "  [ ";
                std::copy(first, last, std::ostream_iterator<double>(os," "));
                os << "]\n";
            }
        }
        os << "]\n";
        return os;
    } // DN

    //----------------------------CLASS HELPER FUNCTIONS------------------------
    
    void SmallMatrix::printNumRowCol() {
        std::cout << "Rows: " << mNumRows << " Cols: " << mNumCols << '\n';
    }

    void SmallMatrix::printMatrix() {
        if (mIsLargeMatrix) {
            for (auto row : mHeapData) {
                printVector(row);
            }
        } else {
            for (int i{0}; i < mNumRows; i++) {
                printVector(mStackData.at(i), mNumCols);
            }
        }
    }

    void SmallMatrix::fillMatrix() {
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            if (mIsLargeMatrix) {
                mHeapData.at(row_index).resize(mNumCols);
                std::fill(mHeapData.at(row_index).begin(), mHeapData.at(row_index).end(), 0);
            } else {
                mStackData.at(row_index).fill(0);
            }
        }
    }

    void SmallMatrix::fillMatrix(double value) {
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            if (mIsLargeMatrix) {
                mHeapData.at(row_index).resize(mNumCols);
                std::fill(mHeapData.at(row_index).begin(), mHeapData.at(row_index).end(), value);
            } else {
                mStackData.at(row_index).fill(value);
            }
        }
    }

    void SmallMatrix::fillRow(int first, int last) {
        for (int row_index{first}; row_index < last; row_index++) {
            if (mIsLargeMatrix) {
                mHeapData.at(row_index).resize(mNumCols);
                std::fill(mHeapData.at(row_index).begin(), mHeapData.at(row_index).end(), 0);
            } else {
                std::fill(mStackData.at(row_index).begin(),
                    mStackData.at(row_index).begin() + last, 0);
            }
        }
    }

    void SmallMatrix::fillCol(int first, int last) {
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            if (mIsLargeMatrix) {
                mHeapData.at(row_index).resize(last);
                std::fill(mHeapData.at(row_index).begin() + first,
                    mHeapData.at(row_index).end(), 0);
            } else {
                std::fill(mStackData.at(row_index).begin() + first,
                    mStackData.at(row_index).begin() + last, 0);
            }
        }
    }

    void SmallMatrix::scalarMultiply(double num, SmallMatrix const& sm)  {
        for (int row_index{0}; row_index < sm.mNumRows; row_index++) {
            if (sm.mIsLargeMatrix) {
                auto first_col = sm.mHeapData.at(row_index).begin();
                auto last_col = sm.mHeapData.at(row_index).end();
                auto result_col = mHeapData.at(row_index).begin();
                std::transform(first_col, last_col, result_col, [num](auto& i) {return num*i;});
            } else {
                auto first_col = sm.mStackData.at(row_index).begin();
                auto last_col = sm.mStackData.at(row_index).end();
                auto result_col = mStackData.at(row_index).begin();
                std::transform(first_col, last_col, result_col, [num](auto& i) {return num*i;});
            }
        }
    }

    void SmallMatrix::scalarMultiply(double num)  {
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            if (mIsLargeMatrix) {
                auto first_col = mHeapData.at(row_index).begin();
                auto last_col = mHeapData.at(row_index).end();
                auto result_col = mHeapData.at(row_index).begin();
                std::transform(first_col, last_col, result_col, [num](auto& i) {return num*i;});
            } else {
                auto first_col = mStackData.at(row_index).begin();
                auto last_col = mStackData.at(row_index).end();
                auto result_col = mStackData.at(row_index).begin();
                std::transform(first_col, last_col, result_col, [num](auto& i) {return num*i;});
            }
        }
    }

    void SmallMatrix::matrixSubtraction(SmallMatrix const& lhs, SmallMatrix const& rhs) {
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            for (int col_index{0}; col_index < mNumCols; col_index++) {
                if (mIsLargeMatrix) {
                    mHeapData.at(row_index).at(col_index) = 
                        lhs.mHeapData.at(row_index).at(col_index) -
                        rhs.mHeapData.at(row_index).at(col_index);
                } else {
                    mStackData.at(row_index).at(col_index) = 
                        lhs.mStackData.at(row_index).at(col_index) -
                        rhs.mStackData.at(row_index).at(col_index);
                }
            }
        }
    }

    void SmallMatrix::matrixSubtraction(SmallMatrix const& sm) {
        //
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            for (int col_index{0}; col_index < mNumCols; col_index++) {
                if (mIsLargeMatrix) {
                    mHeapData.at(row_index).at(col_index) -= 
                        sm.mHeapData.at(row_index).at(col_index);
                } else {
                    mStackData.at(row_index).at(col_index) -= 
                        sm.mStackData.at(row_index).at(col_index);
                }
            }
        }
    }

    void SmallMatrix::matrixAddition(SmallMatrix const& lhs, SmallMatrix const& rhs) {
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            for (int col_index{0}; col_index < mNumCols; col_index++) {
                if (mIsLargeMatrix) {
                    mHeapData.at(row_index).at(col_index) = 
                        lhs.mHeapData.at(row_index).at(col_index) +
                        rhs.mHeapData.at(row_index).at(col_index);
                } else {
                    mStackData.at(row_index).at(col_index) = 
                        lhs.mStackData.at(row_index).at(col_index) +
                        rhs.mStackData.at(row_index).at(col_index);
                }
            }
        }
    }

    void SmallMatrix::matrixAddition(SmallMatrix const& sm) {
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            for (int col_index{0}; col_index < mNumCols; col_index++) {
                if (mIsLargeMatrix) {
                    mHeapData.at(row_index).at(col_index) += 
                        sm.mHeapData.at(row_index).at(col_index);
                } else {
                    mStackData.at(row_index).at(col_index) += 
                        sm.mStackData.at(row_index).at(col_index);
                }
            }
        }
    }

    void SmallMatrix::matrixMultiply(SmallMatrix const& lhs, SmallMatrix const& rhs) {
        auto numCols1 = lhs.size().second;
        // operation for each row and col in resulting matrix
        for (int i{0}; i < mNumRows; i++) {
            for (int j{0}; j < mNumCols; j++) {
                for (int k{0}; k < numCols1; k++) {
                    if (mIsLargeMatrix) {
                        mHeapData.at(i).at(j) += 
                            lhs.mHeapData.at(i).at(k) * rhs.mHeapData.at(k).at(j);
                    } else {
                        mStackData.at(i).at(j) += 
                            lhs.mStackData.at(i).at(k) * rhs.mStackData.at(k).at(j);
                    }
                }
            }
        }
    }

    void SmallMatrix::matrixMultiply(SmallMatrix const& sm) {
        auto lhs(*this);
        // operation for each row and col in resulting matrix
        for (int i{0}; i < this->mNumRows; i++) {
            for (int j{0}; j < this->mNumCols; j++) {
                for (int k{0}; k < mNumCols; k++) {
                    if (this->mIsLargeMatrix) {
                        mHeapData.at(i).at(j) += 
                            lhs.mHeapData.at(i).at(k) * sm.mHeapData.at(k).at(j);
                    } else {
                        mStackData.at(i).at(j) += 
                            lhs.mStackData.at(i).at(k) * sm.mStackData.at(k).at(j);
                    }
                }
                if (mIsLargeMatrix) {mHeapData.at(i).at(j) -= lhs.mHeapData.at(i).at(j);}
                else {mStackData.at(i).at(j) -= lhs.mStackData.at(i).at(j);}
            }
        }
    }

    void SmallMatrix::stackToHeap(int numRows, int numCols) {
        // Copying current stack data to heap data
        mHeapData.resize(mNumRows);
        std::cout << "Beginning stack to heap copy operation\n";
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            mHeapData.at(row_index).resize(mNumCols);
            std::copy(mStackData.at(row_index).begin(), mStackData.at(row_index).begin() + mNumCols, mHeapData.at(row_index).begin()); 
        }

        // Expanding heap data to required size and zero initialisation
        mHeapData.resize(numRows);
        fillRow(mNumRows, numRows);
        fillCol(mNumCols, numCols);
    }

    }  // namespace mtrn2500


    //----------------------------HELPER FUNCTIONS------------------------------
    void printVector(const std::vector<double> vec) {
        for (int unsigned i{0}; i < vec.size(); i++) {
            std::cout << vec.at(i) << " ";
        }
        std::cout << '\n';
    }

    void printVector(const std::array<double, 144> &arr, int size) {
        for (int i{0}; i < size; i++) {
            std::cout << arr.at(i) << " ";
        }
        std::cout << '\n';
    }