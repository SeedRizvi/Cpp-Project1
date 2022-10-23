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
                mIsLargeMatrix = true;
                stackToHeap(numRows, numCols);
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

    void SmallMatrix::insertRow(int numRow, std::vector<double> const& row) {
        auto resized = false;
        if (mNumCols == 0) {
            this->resize(mNumRows + 1, row.size());
            resized = true;
        }
        if ((numRow < 0) or (numRow > mNumRows)) {
            throw std::out_of_range("Specified row is not in [0,maxRows].");
        } else if (row.size() != static_cast<unsigned int>(mNumCols)) {
            throw std::invalid_argument("Row size does not match matrix row size.");
        } else if (resized == false) {
            this->resize(mNumRows + 1, mNumCols);
        }
        auto const smCopy = *this;
        if (mIsLargeMatrix) {
            // Insert row into mHeapData
            mHeapData.at(numRow) = row;
            std::copy(smCopy.mHeapData.begin() + numRow, smCopy.mHeapData.end(),
                mHeapData.begin() + numRow + 1);
        } else {
            // Insert row into mStackData array
            std::copy(row.begin(), row.end(), mStackData.at(numRow).begin());
            std::copy(smCopy.mStackData.begin() + numRow, smCopy.mStackData.end() - 1,
                mStackData.begin() + numRow + 1);
                // Plus one to write after new row
        }                                         
    } // HD

    void SmallMatrix::insertCol(int numCol, std::vector<double> const& col) {
        auto resized = false;
        if (mNumRows == 0) {
            this->resize(col.size(), mNumCols + 1);
            resized = true;
        }
        if ((numCol < 0) or (numCol > mNumCols)) {
            throw std::out_of_range("Specified col is not in [0,maxCols].");
        } else if (col.size() != static_cast<unsigned int>(mNumRows)) {
            throw std::invalid_argument("Col size does not match matrix col size.");
        } else if (resized == false) {
            this->resize(mNumRows, mNumCols + 1);
        }
        auto const smCopy = *this;
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            if (mIsLargeMatrix) {
                // Insert col at row_index into mHeapData
                mHeapData.at(row_index).at(numCol) = col.at(row_index);
                auto copyFirstCol = smCopy.mHeapData.at(row_index).begin();
                auto copyLastCol = smCopy.mHeapData.at(row_index).end();
                std::copy(copyFirstCol + numCol, copyLastCol, 
                    mHeapData.at(row_index).begin() + numCol + 1);
            } else {
                // Insert col at row_index into mStackData
                mStackData.at(row_index).at(numCol) = col.at(row_index);
                auto copyFirstCol = smCopy.mStackData.at(row_index).begin();
                auto copyLastCol = smCopy.mStackData.at(row_index).end();
                std::copy(copyFirstCol + numCol, copyLastCol - 1, 
                    mStackData.at(row_index).begin() + numCol + 1);
                    // Plus one to write after new col
            }
        }   
    } // HD

    void SmallMatrix::eraseRow(int numRow) {
        if  ((numRow < 0) or (numRow >= mNumRows)) {
            throw std::out_of_range("Specified row is out of range");
        }
        if (mIsLargeMatrix) {
            // Plus one to ensure only a single element is erased
            mHeapData.erase(mHeapData.begin() + numRow, mHeapData.begin() + numRow + 1);
            mNumRows -= 1;
        } else {
            auto newMatrix = SmallMatrix(this->mNumRows - 1, this->mNumCols);
            auto firstRow = mStackData.begin();
            auto lastRow = mStackData.end();
            // Copy first row till the row to be deleted
            std::copy(firstRow, firstRow + numRow, newMatrix.mStackData.begin());
            // Copy from 1 row after row to be deleted till end of StackData
            // Thus not copying numRow and hence deleting it after std::move
            std::copy(firstRow + numRow + 1, lastRow, 
                newMatrix.mStackData.begin() + numRow);
            *this = std::move(newMatrix);
        }
    } // HD

    void SmallMatrix::eraseCol(int numCol) {
        if  ((numCol < 0) or (numCol >= mNumCols)) {
            throw std::out_of_range("Specified row is out of range");
        }
        auto newMatrix = SmallMatrix(this->mNumRows, this->mNumCols - 1);  
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            if (mIsLargeMatrix) {
                auto first = mHeapData.at(row_index).begin();
                auto oneAfterNumCol = mHeapData.at(row_index).begin() + numCol + 1;
                mHeapData.at(row_index).erase(first + numCol, oneAfterNumCol);
            } else {
                auto firstCol = mStackData.at(row_index).begin();
                auto lastCol = mStackData.at(row_index).end();
                // For each row, copy first column till the column to be deleted
                std::copy(firstCol, firstCol + numCol, 
                newMatrix.mStackData.at(row_index).begin());
                // Copy from 1 column after column to be deleted for each row
                // Thus not copying numCol and deleting it after std::move
                std::copy(firstCol + numCol + 1, lastCol, 
                newMatrix.mStackData.at(row_index).begin() + numCol);
            }
        }
        if (mIsLargeMatrix) {
            mNumCols -= 1;
        } else {
            *this = std::move(newMatrix);
        }
    } // HD

    bool operator==(SmallMatrix const& lhs, SmallMatrix const& rhs) {
        if ((lhs.mNumRows != rhs.mNumRows) or (lhs.mNumCols != rhs.mNumCols)) {
            return false;
        }
        auto epsilon{0.0000001}; 
        for (int row_index{0}; row_index < lhs.mNumRows; row_index++) {
            for (int col_index{0}; col_index < lhs.mNumCols; col_index++) {
                // Vector-Array comparison
                if (lhs.mIsLargeMatrix and (not rhs.mIsLargeMatrix)) {
                    auto lhsValue = lhs.mHeapData.at(row_index).at(col_index);
                    auto rhsValue = rhs.mStackData.at(row_index).at(col_index);
                    if (abs(lhsValue - rhsValue) > epsilon) {return false;}
                } 
                // Array-Vector comparison
                else if ((not lhs.mIsLargeMatrix) and rhs.mIsLargeMatrix) {
                    auto lhsValue = lhs.mStackData.at(row_index).at(col_index);
                    auto rhsValue = rhs.mHeapData.at(row_index).at(col_index);
                    if (abs(lhsValue - rhsValue) > epsilon) {return false;}
                } 
                // Vector-Vector comparison
                else if (lhs.mIsLargeMatrix and rhs.mIsLargeMatrix) {
                    auto lhsValue = lhs.mHeapData.at(row_index).at(col_index);
                    auto rhsValue = rhs.mHeapData.at(row_index).at(col_index);
                    if (abs(lhsValue - rhsValue) > epsilon) {return false;}
                } 
                // Array-Array comparison
                else {
                    auto lhsValue = lhs.mStackData.at(row_index).at(col_index);
                    auto rhsValue = rhs.mStackData.at(row_index).at(col_index);
                    if (abs(lhsValue - rhsValue) > epsilon) {return false;}
                }
            }
        }
        return true;
    } // P  

    bool operator!=(SmallMatrix const& lhs, SmallMatrix const& rhs) {
        if ((lhs.mNumRows != rhs.mNumRows) or (lhs.mNumCols != rhs.mNumCols)) {
            return true;
        } 
        auto epsilon{0.0000001}; 
        for (int row_index{0}; row_index < lhs.mNumRows; row_index++) {
            for (int col_index{0}; col_index < lhs.mNumCols; col_index++) {
                // Vector-Array comparison
                if (lhs.mIsLargeMatrix and (not rhs.mIsLargeMatrix)) {
                    auto lhsValue = lhs.mHeapData.at(row_index).at(col_index);
                    auto rhsValue = rhs.mStackData.at(row_index).at(col_index);
                    if (abs(lhsValue - rhsValue) > epsilon) {return true;}
                } 
                // Array-Vector comparison
                else if ((not lhs.mIsLargeMatrix) and rhs.mIsLargeMatrix) {
                    //do something
                    auto lhsValue = lhs.mStackData.at(row_index).at(col_index);
                    auto rhsValue = rhs.mHeapData.at(row_index).at(col_index);
                    if (abs(lhsValue - rhsValue) > epsilon) {return true;}
                } 
                // Vector-Vector comparison
                else if (lhs.mIsLargeMatrix and rhs.mIsLargeMatrix) {
                    auto lhsValue = lhs.mHeapData.at(row_index).at(col_index);
                    auto rhsValue = rhs.mHeapData.at(row_index).at(col_index);
                    if (abs(lhsValue - rhsValue) > epsilon) {return true;}
                } 
                // Array-Array comparison
                else {
                    auto lhsValue = lhs.mStackData.at(row_index).at(col_index);
                    auto rhsValue = rhs.mStackData.at(row_index).at(col_index);
                    if (abs(lhsValue - rhsValue) > epsilon) {return true;}
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

    SmallMatrix transpose(SmallMatrix const& sm) {
        auto result = SmallMatrix(sm.mNumCols, sm.mNumRows);
        for (int row_index{0}; row_index < sm.mNumRows; row_index++) {
            for (int col_index{0}; col_index < sm.mNumCols; col_index++) {
                if (result.mIsLargeMatrix) {
                    result.mHeapData.at(col_index).at(row_index) = 
                        sm.mHeapData.at(row_index).at(col_index);
                } else {
                    result.mStackData.at(col_index).at(row_index) = 
                        sm.mStackData.at(row_index).at(col_index);
                }
            }
        }
        return result;
    } // DN

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
        auto result = SmallMatrix(lhs.mNumRows, sm.mNumCols);
        auto numCols1 = lhs.size().second;
        // operation for each row and col in resulting matrix
        for (int i{0}; i < result.mNumRows; i++) {
            for (int j{0}; j < result.mNumCols; j++) {
                for (int k{0}; k < numCols1; k++) {
                    if (result.mIsLargeMatrix) {
                        result.mHeapData.at(i).at(j) += 
                            lhs.mHeapData.at(i).at(k) * sm.mHeapData.at(k).at(j);
                    } else {
                        result.mStackData.at(i).at(j) += 
                            lhs.mStackData.at(i).at(k) * sm.mStackData.at(k).at(j);
                    }
                }
            }
        }
        *this = std::move(result);
    }

    void SmallMatrix::stackToHeap(int numRows, int numCols) {
        // Copying current stack data to heap data
        mHeapData.resize(mNumRows);
        for (int row_index{0}; row_index < mNumRows; row_index++) {
            mHeapData.at(row_index).resize(mNumCols);
            std::copy(mStackData.at(row_index).begin(), mStackData.at(row_index).begin() + mNumCols, mHeapData.at(row_index).begin()); 
        }
        // Expanding heap data to required size and zero initialisation
        mHeapData.resize(numRows);
        fillRow(mNumRows, numRows);
        mNumRows = numRows;
        fillCol(mNumCols, numCols);
        mNumCols = numCols;
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