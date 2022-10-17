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

    SmallMatrix::~SmallMatrix() {} // P

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

    std::vector<double*> SmallMatrix::row(int numRow) { return {}; } // CR

    std::vector<double const*> SmallMatrix::row(int numRow) const { return {}; } // CR

    std::vector<double*> SmallMatrix::col(int numCol) { return {}; } // CR

    std::vector<double const*> SmallMatrix::col(int numCol) const { return {}; } // CR

    std::pair<int, int> SmallMatrix::size() const { // P  
        std::pair<int, int> Size(mNumRows, mNumCols);
        return Size;
    } 

    bool SmallMatrix::isSmall() const { // P  
        return !mIsLargeMatrix;
    }

    void SmallMatrix::resize(int numRows, int numCols) {} // DN

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

    SmallMatrix operator+(SmallMatrix const& lhs, SmallMatrix const& rhs) { return {}; } // CR

    SmallMatrix operator-(SmallMatrix const& lhs, SmallMatrix const& rhs) { return {}; } // CR

    SmallMatrix operator*(SmallMatrix const& lhs, SmallMatrix const& rhs) { return {}; } // DN

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

    SmallMatrix& SmallMatrix::operator+=(SmallMatrix const& sm) { return *this; } // CR

    SmallMatrix& SmallMatrix::operator-=(SmallMatrix const& sm) { return *this; } // CR

    SmallMatrix& SmallMatrix::operator*=(SmallMatrix const& sm) { return *this; } // DN

    SmallMatrix& SmallMatrix::operator*=(double s) { return *this; } // CR

    SmallMatrix transpose(SmallMatrix const& sm) { return {}; } // DN

    std::ostream& operator<<(std::ostream& os, SmallMatrix const& sm) { return os; } // DN

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
        for (int row_index{0}; row_index < mNumRows; row_index++) { // Changed cond from mNumRows to Cols
            if (mIsLargeMatrix) {
                mHeapData.at(row_index).resize(mNumCols);
                std::fill(mHeapData.at(row_index).begin(), mHeapData.at(row_index).end(), 0);
            } else {
                mStackData.at(row_index).fill(0);
            }
        }
    }

    void SmallMatrix::fillMatrix(double value) {
        for (int row_index{0}; row_index < mNumRows; row_index++) { // Changed cond from mNumRows to Cols
            if (mIsLargeMatrix) {
                mHeapData.at(row_index).resize(mNumCols);
                std::fill(mHeapData.at(row_index).begin(), mHeapData.at(row_index).end(), value);
            } else {
                mStackData.at(row_index).fill(value);
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