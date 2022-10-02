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
        // If req matrix is large, change HeapData size to match
        // the required size.
        if (mIsLargeMatrix) {
            mHeapData.resize(mNumRows); 
        }
        // Throw out_of_range exception if trying to access row / col out of bounds.
        for (int row_index{0}; row_index < mNumRows; row_index++) { // Changed cond from mNumRows to Cols
            if (mIsLargeMatrix) {
                mHeapData.at(row_index).resize(mNumCols);
                std::fill(mHeapData.at(row_index).begin(), mHeapData.at(row_index).end(), 0);
            } else {
                mStackData.at(row_index).fill(0);
            }
        }
    }

    SmallMatrix::SmallMatrix(int numRows, int numCols, double value) {}

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
                mHeapData.at(row_index).reserve(mNumCols);
                std::copy(row.begin(), row.end(), mHeapData.at(row_index).begin());
            } else {
                std::transform(row.begin(), row.end(), mStackData.at(row_index).begin(),
                            [](auto const& e) { return e; });
            }
            row_index++;
        }
    }

    SmallMatrix::SmallMatrix(SmallMatrix const& sm) {}

    SmallMatrix::SmallMatrix(SmallMatrix&& sm) {}

    SmallMatrix& SmallMatrix::operator=(SmallMatrix const& sm) { return *this; }

    SmallMatrix& SmallMatrix::operator=(SmallMatrix&& sm) { return *this; }

    SmallMatrix::~SmallMatrix() {}

    double& SmallMatrix::operator()(int numRow, int numCol) { return mStackData[0][0]; }

    const double& SmallMatrix::operator()(int numRow, int numCol) const { return mStackData[0][0]; }

    std::vector<double*> SmallMatrix::row(int numRow) { return {}; }

    std::vector<double const*> SmallMatrix::row(int numRow) const { return {}; }

    std::vector<double*> SmallMatrix::col(int numCol) { return {}; }

    std::vector<double const*> SmallMatrix::col(int numCol) const { return {}; }

    std::pair<int, int> SmallMatrix::size() const {
        std::pair<int, int> Size(mNumRows, mNumCols);
        return Size;
    }

    bool SmallMatrix::isSmall() const {
        // if (mNumRows * mNumCols >= mSmallSize) {return false;} else {return true;}
        return !mIsLargeMatrix;
    }

    void SmallMatrix::resize(int numRows, int numCols) {}

    void SmallMatrix::insertRow(int numRow, std::vector<double> const& row) {}

    void SmallMatrix::insertCol(int numCol, std::vector<double> const& col) {}

    void SmallMatrix::eraseRow(int numRow) {}

    void SmallMatrix::eraseCol(int numCol) {}

    bool operator==(SmallMatrix const& lhs, SmallMatrix const& rhs) { return false; }

    bool operator!=(SmallMatrix const& lhs, SmallMatrix const& rhs) { return false; }

    SmallMatrix operator+(SmallMatrix const& lhs, SmallMatrix const& rhs) { return {}; }

    SmallMatrix operator-(SmallMatrix const& lhs, SmallMatrix const& rhs) { return {}; }

    SmallMatrix operator*(SmallMatrix const& lhs, SmallMatrix const& rhs) { return {}; }

    SmallMatrix operator*(double s, SmallMatrix const& sm) { return {}; }

    SmallMatrix operator*(SmallMatrix const& sm, double s) { return {}; }

    SmallMatrix& SmallMatrix::operator+=(SmallMatrix const& sm) { return *this; }

    SmallMatrix& SmallMatrix::operator-=(SmallMatrix const& sm) { return *this; }

    SmallMatrix& SmallMatrix::operator*=(SmallMatrix const& sm) { return *this; }

    SmallMatrix& SmallMatrix::operator*=(double s) { return *this; }

    SmallMatrix transpose(SmallMatrix const& sm) { return {}; }

    std::ostream& operator<<(std::ostream& os, SmallMatrix const& sm) { return os; }

    //----------------------------GETTER FUNCTIONS------------------------------
    
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
    }  // namespace mtrn2500


    //----------------------------HELPER FUNCTIONS------------------------------
    void printVector(const std::vector<double> vec) {
        for (int i{0}; i < vec.size(); i++) {
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
