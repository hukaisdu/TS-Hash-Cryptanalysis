#include<iostream>
#include<random>
#include<cmath>
#include<unordered_set>

using namespace std;

typedef int F2;
typedef int F3;



const int lambda = 16;
const int NN = 2 * lambda; // 64

const int HalfNN = NN / 2; // 32

const int TT = int( lambda / log2(3) ); // 32 / 1.5 = 24
const int RR = TT / 2; // 12 
                       //
                       //
/* utility functions */
template< typename T >
void printMatrix( T ** M, int NRow, int NCol )
{
    for ( int i = 0; i < NRow; i++ )
    {
        for ( int j = 0; j < NCol; j++ )
            cout << static_cast<int> ( M[i][j] );
        cout << endl;
    }
}

inline F3 mod3(F3 x) 
{ 
    return (x % 3 + 3) % 3; 
}

inline F3 inverse_mod3(F3 x) 
{ 
    return x == 1 ? 1 : 2; 
} 

void generateCombinations( vector<vector<F3>>& res, vector<F3>& current, int depth, int k) 
{
    if (depth == k) 
    {
        res.push_back(current);
        return;
    }
    for (F3 val = 0; val < 3; ++val) 
    {
        current[depth] = val;
        generateCombinations(res, current, depth + 1, k);
    }
}

bool EQ( int * A, int * B, int N )
{
    for ( int i = 0; i < N; i++ )
        if ( A[i]  != B[i] )
            return false;
    return true;
}

bool EQ( int * A, const vector<int> & B, int N )
{
    for ( int i = 0; i < N; i++ )
        if ( A[i]  != B[i] )
            return false;
    return true;
}

bool EQ( const vector<int> & A, int * B, int N )
{
    for ( int i = 0; i < N; i++ )
        if ( A[i]  != B[i] )
            return false;
    return true;
}

bool EQ( const vector<int> & A, const vector<int> & B, int N )
{
    for ( int i = 0; i < N; i++ )
        if ( A[i]  != B[i] )
            return false;
    return true;
}


/*
inline ostream & operator<< ( ostream & os, int t )  
{
    os << static_cast<int> ( t );
    return os;
}
*/

void printVec( int * t, int N )
{
    for ( int i = 0; i < N; i++ )
        cout << t[i];
    cout << endl;
}

void printVec( vector<int> & t, int N )
{
    for ( int i = 0; i < N; i++ )
        cout << t[i];
    cout << endl;
}


/* the wPRF implementation */

void MatrixMulVector2( F2 ** K, F2 * x, F2 * y, int N = NN )
{
    for ( int i = 0; i < N; i++ )
    {
        for ( int j = 0; j < N; j++ )
            y[i] = ( y[i] + K[i][j] * x[j] ) % 2;
    }
}

void MatrixMulVector3( F3 ** B, F3 * x, F3 * y, int T = TT, int N = NN )
{
    for ( int i = 0; i < T; i++ )
    {
        for ( int j = 0; j < N; j++ )
            y[i] = ( y[i] + B[i][j] * x[j] ) % 3;
    }
}

void wPRF( F2 ** K, F3 ** B, F2 * x, F3 * y )
{
    F2 * xx = new F2[NN];

    MatrixMulVector2( K, x, xx );
    MatrixMulVector3( B, xx, y );

    delete [] xx;
}

/* the wPRF implementation finished */


/* offline phase */

int rankF3( F3** B, int TT, int NN) 
{
    F3** mat = new F3*[TT];
    for (int i = 0; i < TT; ++i ) 
    {
        mat[i] = new F3[NN];
        for (int j = 0; j < NN; ++j) 
            mat[i][j] = B[i][j]; // 转换为非负模3值
    }

    int rank = 0;
    for (int col = 0; col < NN && rank < TT; ++col) 
    {
        // 寻找当前列主元行
        int pivot_row = -1;
        for (int i = rank; i < TT; ++i) 
        {
            if (mat[i][col] != 0) 
            {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) 
            continue; // 跳过全零列

        // 交换行（直接交换指针）
        //int* temp = mat[rank];
        //mat[rank] = mat[pivot_row];
        //mat[pivot_row] = temp;

        swap( mat[rank], mat[ pivot_row ] );

        // 主元归一化为1
        int inv = inverse_mod3( mat[rank][col] );
        for (int j = 0; j < NN; ++j) 
            mat[rank][j] = mod3(mat[rank][j] * inv);

        // 消去其他行
        for (int i = rank + 1; i < TT; ++i) 
        {
            if ( mat[i][col] == 0 ) 
                continue;
            int factor = mod3( -mat[i][col] );
            for (int j = col; j < NN; ++j) 
                mat[i][j] = mod3( mat[i][j] + ( factor * mat[rank][j] ) );
        }
        ++rank;
    }

    // 释放临时矩阵内存
    for (int i = 0; i < TT; ++i) {
        delete[] mat[i];
    }
    delete[] mat;

    return rank;
}

int rankF2(F2** KEY, int NN) 
{
    F2** mat = new F2*[NN];
    for (int i = 0; i < NN; ++i) 
    {
        mat[i] = new F2[NN];
        for (int j = 0; j < NN; ++j) 
            mat[i][j] = KEY[i][j];
    }

    int rank = 0;
    for (int col = 0; col < NN && rank < NN; ++col) 
    {
        int pivot_row = -1;
        for (int i = rank; i < NN; ++i) 
        {
            if (mat[i][col] != 0) 
            {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) 
            continue;

        swap( mat[rank], mat[ pivot_row ] );

        for (int i = 0; i < NN; ++i) 
        {
            if (i == rank || mat[i][col] == 0) 
                continue;
            for (int j = col; j < NN; ++j) 
                mat[i][j] = ( mat[i][j] + mat[rank][j] ) % 2;
        }
        ++rank;
    }

    for (int i = 0; i < NN; ++i) {
        delete[] mat[i];
    }
    delete[] mat;

    return rank;
}


// N is the size of the 
int Gauss_F2( F2 ** mat, int NRow, int NCol )   
{
    for ( int col = 0; col < NCol; col++ )
    {
        int currentRow = col;

        // current row, find the first pivot
        for ( int row = currentRow; row < NRow; row ++ )
        {
            if ( mat[row][col] == 0 )
                continue;
            else
                //swapRow( mat, NCol, currentRow, row );
                swap( mat[ currentRow ], mat[ row ] );
        }

        for ( int row = currentRow + 1; row < NRow; row ++ )
        {
            if ( mat[row][col] == 1 )
                for ( int j = 0; j < NCol; j++ )
                    mat[row][j] = ( mat[row][j] + mat[currentRow][j] ) % 2;
        }
    }

    // eliminate the above elements of the pivot
    for ( int col = 0; col < NCol; col++ )
    {
        if ( mat[col][col] == 1 )
        {
            for ( int i = 0; i < col; i++ ) // for all rows above pivot
            {
                if ( mat[i][col] == 1 )
                    for ( int j = 0; j < NCol; j++ )
                        mat[i][j] = ( mat[i][j] + mat[col][j] ) % 2; 
            }
        }
    }

    int rank = 0;
    for ( int i = 0; i < NCol; i++ )
        if ( mat[i][i] == 1 )
            rank += 1;

    return rank;
}

//void genGuesseDiffs ( F2 ** mat, int NRow, int NCol, F2 ** diff ) 
vector<vector<F2> > genGuesseDiffs ( F2 ** mat, int NRow, int NCol ) 
{
    //printMatrix( mat, NRow, NCol );
    //getchar();

    vector< vector<F2> > diff;

    F2 * X = new F2 [ NCol ]; 

    for ( int x = 0; x < ( 1 << NCol ); x++ )
    {
        for ( int i = 0; i < NCol; i++ )
            X[i] = x >> ( NCol - 1 - i ) & 0x1; 

        vector<F2> vec;

        for ( int i = 0; i < NRow; i++ )
        {
            F2 temp = 0;
            for ( int j = 0; j < NCol; j++ )
                if ( mat[i][j] == 1 )
                    temp = ( temp + X[j] ) % 2;   
            vec.push_back( temp );
        }

       // for ( int i = 0; i < HalfNN; i++ )
       //     cout << vec[i];
       // cout << endl;
       // getchar();

        diff.push_back( vec );
    }

    delete [] X;

    return diff;
}
/* offline phase finished */

/* solve a F3 linear equation systems  */

// T rows and N cols
vector<vector<F3>> solveF3(F3** A, F3* b, int T, int N) 
{
    vector<vector<F3>> aug(T, vector<F3>(N + 1));

    // init the augmented matrix
    for (int i = 0; i < T; ++i) 
    {
        for (int j = 0; j < N; ++j) 
            aug[i][j] = mod3(A[i][j]);
        aug[i][N] = mod3(b[i]);
    }

    vector<int> pivot_cols;
    int rank = 0;

    // 高斯消元
    for (int col = 0; col < N; ++col) 
    {
        int pivot = rank;
        while ( ( pivot < T ) && ( aug[pivot][col] == 0 ) ) 
            ++pivot;
        if (pivot == T) 
            continue; // free variable

        swap( aug[rank], aug[pivot] );

        int inv = inverse_mod3(aug[rank][col]);

        for (int j = 0; j <= N; ++j)
            aug[rank][j] = mod3( aug[rank][j] * inv );

        for (int i = 0; i < T; ++i) 
        {
            if ( ( i != rank ) && ( aug[i][col] != 0 ) ) 
            {
                int factor = mod3(-aug[i][col]);
                for (int j = col; j <= N; ++j)
                    aug[i][j] = mod3(aug[i][j] + factor * aug[rank][j]);
            }
        }
        pivot_cols.push_back(col);
        ++rank;
    }

    // 检查矛盾方程
    for (int i = rank; i < T; ++i)
        if (aug[i][T] != 0) 
            return {}; // no solution

    // 确定自由变量
    unordered_set<int> pivot_set(pivot_cols.begin(), pivot_cols.end());

    vector<int> free_cols;
    for (int col = 0; col < N; ++col)
        if (!pivot_set.count(col)) 
            free_cols.push_back(col);

    vector<vector<F3>> solutions;

    // 生成所有自由变量组合
    int k = free_cols.size();
    vector<vector<F3>> free_combos;
    if ( k > 0 ) 
    {
        vector<F3> tmp(k, 0);
        generateCombinations(free_combos, tmp, 0, k);

        // 构造解集
        for (auto& combo : free_combos) 
        {
            vector<F3> sol(N, 0);
            // 设置自由变量值
            for (int i = 0; i < k; ++i) 
                sol[free_cols[i]] = combo[i];

            // 回代求解主变量
            for (int i = rank-1; i >= 0; --i) 
            {
                int pc = pivot_cols[i];
                F3 sum = 0;
                for (int j = pc+1; j < N; ++j)
                    sum = mod3(sum + mod3(aug[i][j] * sol[j]));
                sol[pc] = mod3(aug[i][N] - sum);
            }
            solutions.push_back(sol);
        }
    }
    else
    {
        // 唯一解：直接回代所有主变量
        vector<F3> sol(N, 0);
        for (int i = rank - 1; i >= 0; --i) {
            int pc = pivot_cols[i];
            F3 sum = 0;
            for (int j = pc + 1; j < N; ++j)
                sum = mod3(sum + mod3(aug[i][j] * sol[j]));
            sol[pc] = mod3(aug[i][N] - sum);
        }
        solutions.push_back(sol);
    }


    // 唯一解的特殊处理
    //if ( ( k == 0 ) && ( !solutions.empty() ) )
    //    return { solutions[0] };

    return solutions;
}


// 辅助函数：递归生成自由变量的所有可能组合 (0/1)
void generateCombinationsF2( vector<vector<F2>>& res, vector<F2>& current, int depth, int k) 
{
    if (depth == k) 
    {
        res.push_back(current);
        return;
    }
    for (F2 val = 0; val < 2; ++val) 
    {
        current[depth] = val;
        generateCombinationsF2(res, current, depth + 1, k);
    }
}


vector<vector<F2>> solveF2(F2** M, F2* b, int T, int N) 
{
    vector<vector<F2>> aug(T, vector<F2>(N + 1));

    // 初始化增广矩阵 [M | b]
    for (int i = 0; i < T; ++i) {
        for (int j = 0; j < N; ++j)
            aug[i][j] = M[i][j] % 2;
        aug[i][N] = b[i] % 2;
    }

    vector<int> pivot_cols; // 主元列索引
    int rank = 0;

    // 高斯消元
    for (int col = 0; col < N; ++col) 
    {
        // 寻找当前列的主元行
        int pivot = rank;
        while ( ( pivot < T ) && ( aug[pivot][col] == 0 ) )
            ++pivot;
        if (pivot == T)
            continue; // 自由变量列，跳过处理

        // 交换主元行到当前秩行
        swap(aug[rank], aug[pivot]);

        // 消去其他行的当前列元素（模2加）
        for (int i = 0; i < T; ++i) 
        {
            if ( ( i != rank ) && ( aug[i][col] != 0 ) ) 
            {
                for (int j = col; j <= N; ++j)
                    aug[i][j] = (aug[i][j] + aug[rank][j]) % 2;
            }
        }
        pivot_cols.push_back(col);
        ++rank;
    }

    //for ( int i = 0; i < T; i++ )
   // {
    //    for ( int j = 0; j <= N; j++ )
    //        cout << aug[i][j];
    //    cout << endl;
    //}


    //getchar();


    // 检查矛盾方程：0x + 0y + ... = 1
    for (int i = rank; i < T; ++i) {
        if (aug[i][N] != 0)
            return {}; // 无解
    }

    // 确定自由变量列（非主元列）
    unordered_set<int> pivot_set(pivot_cols.begin(), pivot_cols.end());
    vector<int> free_cols;
    for (int col = 0; col < N; ++col) {
        if (!pivot_set.count(col))
            free_cols.push_back(col);
    }

    vector<vector<F2>> solutions;
    int k = free_cols.size(); // 自由变量数量

    if (k > 0) 
    {
        // 生成所有自由变量组合（共 2^k 种）
        vector<vector<F2>> free_combos;
        vector<F2> tmp(k, 0);
        generateCombinationsF2(free_combos, tmp, 0, k);

        // 构造解集
        for (auto& combo : free_combos) {
            vector<F2> sol(N, 0);
            // 设置自由变量值
            for (int i = 0; i < k; ++i)
                sol[free_cols[i]] = combo[i];

            // 回代求解主变量
            for (int i = rank - 1; i >= 0; --i) {
                int pc = pivot_cols[i]; // 主元列
                F2 sum = 0;
                for (int j = pc + 1; j < N; ++j)
                    sum = (sum + aug[i][j] * sol[j]) % 2;
                sol[pc] = (aug[i][N] - sum) % 2;
                sol[pc] = (sol[pc] + 2) % 2; // 确保非负
            }
            solutions.push_back(sol);
        }
    } else {
        // 唯一解：直接回代
        vector<F2> sol(N, 0);
        for (int i = rank - 1; i >= 0; --i) {
            int pc = pivot_cols[i];
            F2 sum = 0;
            for (int j = pc + 1; j < N; ++j)
                sum = (sum + aug[i][j] * sol[j]) % 2;
            sol[pc] = (aug[i][N] - sum) % 2;
            sol[pc] = (sol[pc] + 2) % 2; // 确保非负
        }
        solutions.push_back(sol);
    }

    return solutions;
}


/*
int main() {
    // 示例：x + 2y = 1, 2x + y = 2 (mod 3)
    vector<vector<F3>> A = {{1, 2}, {2, 2}};

    vector<F3> b = {1, 2};
    vector<F3> x = solveF3(A, b, 2);

    if (x.empty()) cout << "无解" << endl;
    else {
        cout << "解为：";
        for (int xi : x) cout << xi << " ";
        // 输出：解为：1 2 （即 x=1, y=2）
    }
    return 0;
}
*/

int main()
{
    random_device rd;
    mt19937 gen ( rd() );
    uniform_int_distribution<F2> dis(0, 1);
    uniform_int_distribution<F3> disF3(0, 2);

    cout << "-----------------------------------------------------------" << endl;
    cout << "Init the wPRF " << endl;
    cout << "PRF parameters: N = " << NN << " T = " << TT  << " lambda = " << lambda << endl;

    /* initialize the wPRF */ 
    /* generate the key */
    F2 ** KEY = new F2 * [NN];
    for ( int i = 0; i < NN; i++ ) 
        KEY[i] = new F2[NN];

    for ( int i = 0; i < NN; i++ )
        KEY[0][i]= dis( gen ); // init the key
                               //
    for ( int i = 1; i < NN; i++ )
        for ( int j = 0; j < NN; j++ ) 
            KEY[i][ j ] = KEY[0][ ( i + j ) % NN ];

    //printMatrix( KEY, NN, NN );

    // compute the rank of KEY
    int rank_key = rankF2( KEY, NN );   

    if ( rank_key != NN )
    {
        cout << "The key rank is not full, run again" << endl;
        return -1;
    }

    cout << "Init the key... Done! " << endl;
    for ( int i = 0; i < NN; i++ )
        cout << static_cast<int>( KEY[0][i] );
    cout << endl;

    F3 ** B = new F3 * [TT];
    for ( int i = 0; i < TT; i++ )
        B[i] = new F3[ NN ];

    for ( int i = 0; i < TT; i++ )
        for ( int j = 0; j < NN; j++ )
            B[i][j] = disF3( gen );

    int rankB = rankF3( B, TT, NN );

    if ( rankB < TT )
    {
        cout << "B does not reach TT rank " << endl;
        return -2;
    }

    cout << "Init the B matrix... Done! " << endl;
    printMatrix<F3> ( B, TT, NN );

    // attack parameter 
    cout << "-----------------------------------------------------------" << endl;
    cout << "Offline Phase " << endl; 

    cout << "Attack Rank R = " << RR << endl;

    F2 ** mat = new F2 * [ HalfNN ];
    for ( int i = 0; i < HalfNN; i++ )
        mat[i] = new F2[ HalfNN ];

    F2 ** matT = new F2 * [HalfNN];
    for ( int i = 0; i < HalfNN; i++ )
        matT[i] = new F2[ RR ];

    F2 ** V = new F2 * [HalfNN];
    for ( int i = 0; i < HalfNN; i++ )
        V[i] = new F2[ HalfNN ];

    cout << "Search for an input with difference rank = " << RR << endl;
    cout << "The expected complexity is (Proposition 1) 2^{R-1-N/2} = 2^" << ( RR - 1 - HalfNN ) << endl;

    F2 * Input1 = new F2[ NN ];
    F3 * Output1 = new F3[ TT ];

    F2 * Input2 = new F2[ NN ];
    F3 * Output2 = new F3[ TT ];

    long long query = 0;

    vector< vector<F2> > diff1; // diff space
    vector< vector<F2> > diff2; // diff space

    for ( query = 0; query < ( 1LL << 32 ); query++ )
    {
        // the input we are waiting for
        for ( int i = 0; i < NN; i++ )
            Input1[i] = dis( gen );

        for ( int i = 0; i < HalfNN; i++ )
            V[0][i] = ( Input1[i] + Input1[ i + HalfNN ] ) % 2;

        // the difference rotations
        for ( int row = 1; row < HalfNN; row++ )
            for ( int j = row; j < HalfNN + row; j++ )
                V[row][j % HalfNN] = V[0][ j - row ];

        for ( int i = 0; i < HalfNN; i++ )
            for ( int j = 0; j < HalfNN; j++ )
                mat[j][i] = V[i][j];

        int rank = Gauss_F2( mat, HalfNN, HalfNN );   

        if ( rank == RR )
        {
            for ( int i = 0; i < rank; i++ )
                for ( int j = 0; j < HalfNN; j++ )
                    matT[j][i] = mat[i][j]; 

            diff1 = genGuesseDiffs ( matT, HalfNN, rank ); 

            cout << "Find the first one! The practical query complexity is  " << query << " = 2^" << log2( query ) << endl;
            break;
        }
    }

    for ( query = 0; query < ( 1LL << 32 ); query++ )
    {
        // the input we are waiting for
        for ( int i = 0; i < NN; i++ )
            Input2[i] = dis( gen );

        for ( int i = 0; i < HalfNN; i++ )
            V[0][i] = ( Input2[i] + Input2[ i + HalfNN ] ) % 2;

        // the difference rotations
        for ( int row = 1; row < HalfNN; row++ )
            for ( int j = row; j < HalfNN + row; j++ )
                V[row][j % HalfNN] = V[0][ j - row ];

        for ( int i = 0; i < HalfNN; i++ )
            for ( int j = 0; j < HalfNN; j++ )
                mat[j][i] = V[i][j];

        int rank = Gauss_F2( mat, HalfNN, HalfNN );   

        if ( rank == RR )
        {
            for ( int i = 0; i < rank; i++ )
                for ( int j = 0; j < HalfNN; j++ )
                    matT[j][i] = mat[i][j]; 

            diff2 = genGuesseDiffs ( matT, HalfNN, rank ); 

            cout << "Find the second one! The practical query complexity is  " << query << " = 2^" << log2( query ) << endl;
            break;
        }
    }

    cout << "The first input is " << endl;
    cout << "I1 = \t "; 
    for ( int i = 0; i < NN; i++ )
        cout << static_cast<int> ( Input1[i] );
    cout << endl;

    cout << "The first output is " << endl;
    cout << "O1 = \t ";
    wPRF( KEY, B, Input1, Output1 );

    for ( int i = 0; i < TT; i++ )
        cout << static_cast<int> ( Output1[i] );
    cout << endl;

    cout << "Transform the output to the F3*-homomorphism output " << endl;
    for ( int i = 0; i < TT; i++ )
        for ( int j = 0; j < NN; j++ )
            Output1[i] = ( Output1[i] + B[i][j] ) % 3;

    cout << "The F3*-homomorphism of the first output is " << endl;
    cout << "O1' = \t";
    for ( int i = 0; i < TT; i++ )
        cout << static_cast<int> ( Output1[i] );
    cout << endl;

    cout << endl;

    cout << "The second input is " << endl;
    cout << "I2 = \t "; 
    for ( int i = 0; i < NN; i++ )
        cout << static_cast<int> ( Input2[i] );
    cout << endl;

    cout << "The second output is " << endl;
    cout << "O2 = \t ";
    wPRF( KEY, B, Input2, Output2 );

    for ( int i = 0; i < TT; i++ )
        cout << static_cast<int> ( Output2[i] );
    cout << endl;

    cout << "Transform the output to the F3*-homomorphism output " << endl;
    for ( int i = 0; i < TT; i++ )
        for ( int j = 0; j < NN; j++ )
            Output2[i] = ( Output2[i] + B[i][j] ) % 3;

    cout << "The F3*-homomorphism of the second output is " << endl;
    cout << "O2' = \t";
    for ( int i = 0; i < TT; i++ )
        cout << static_cast<int> ( Output2[i] );
    cout << endl;

    cout << "first mid " << endl;
    F2 * mid1 = new F2 [ NN ];
    MatrixMulVector2( KEY, Input1, mid1, NN );
    printVec( mid1, NN );

    cout << "second mid " << endl;
    F2 * mid2 = new F2 [ NN ];
    MatrixMulVector2( KEY, Input2, mid2, NN );
    printVec( mid2, NN );

    cout << "-----------------------------------------------------------" << endl;
    cout << "Online phase " << endl;

    // guess the values in the first half, and generate the second half
    // will guess the values of first ( HalfNN - TT ) variables
    //
    F3 * tempO = new F3[TT];
    F3 ** tempB =  new F3 * [TT];
    for ( int i = 0; i < TT; i++ )
        tempB[i] = new F3[TT];

    vector< vector<F2> > solutions;

    cout << "Online phase complexity: 2 x 2^" << ( RR + HalfNN - TT ) << endl;   

    //cout << diff1.size() << " " << diff2.size() << endl;
    
    vector< vector<F2>> solution1;
    vector< vector<F2>> solution2;

    //for ( int i = 0; i < 32; i++ )
    //{
    //    for ( int j = 0; j < HalfNN; j++ )
    //    {
    //        cout << diff1[i][j];
    //    }
    //    cout << endl;
    //    for ( int j = 0; j < HalfNN; j++ )
    //    {
    //        cout << diff2[i][j];
    //    }
    //    cout << endl;
    //}

    //getchar();

    for ( auto it : diff2 ) // for each difference
    //for ( int index = 0; index < ( 1 << RR ); index ++ )
    {
        // **********|vvvvvvvvvv|**********|vvvvvvvvvv
        for ( int x = 0; x < ( 1 << (HalfNN - TT) ); x++ ) // guess the values for the first HalfNN - TT 
        {
            // prepare the coefficients
            for ( int t = 0; t < TT; t++ )
                for ( int j = 0; j < TT; j++ )
                    tempB[t][j] = ( B[t][ HalfNN - TT + j ] + ( B[t][HalfNN + HalfNN - TT + j ] * ( it[ HalfNN - TT + j] + 1 ) ) ) % 3;

            // prepare the output1
            for ( int t = 0; t < TT; t++ )
            {
                F3 reminder = 0;

                for ( int i = 0; i < ( HalfNN - TT ); i++ )
                {
                    reminder = ( reminder + B[t][i] * ( ( x >> ( HalfNN - TT - 1 - i ) & 0x1 ) + 1 ) ) % 3; // + 1 to make it into F3*
                    reminder = ( reminder + B[t][i + HalfNN] * ( ( x >> ( HalfNN - TT - 1 - i ) & 0x1 ) + 1 ) * ( it[i] + 1 ) ) % 3; // F3* is multiplicative group
                }

                tempO[t] = mod3 ( Output1[t] - reminder ); 
            }

            auto res1 = solveF3( tempB, tempO, TT, TT ); 

            if ( res1.size() == 0 )
                continue;

            for ( auto res : res1 )
            {
                bool flag = true;
                for ( int i = 0; i < TT; i++ )
                {
                    if ( res[i] == 0 )
                    {
                        flag = false; 
                        break;
                    }
                }
                if ( flag == false )
                    continue;

                vector<F2> sol ( NN );
                for ( int i = 0; i < ( HalfNN - TT ); i++ )
                    sol[i] = x >> ( HalfNN - TT - 1 - i ) & 0x1;

                for ( int i = 0; i < TT; i++ )
                    sol[HalfNN - TT + i] = ( res[i] - 1 );

                for ( int i = 0; i < ( HalfNN - TT ); i++ )
                    sol[HalfNN + i] = ( ( x >> ( HalfNN - TT - 1 - i ) & 0x1 ) ^ it[i] );

                for ( int i = 0; i < TT; i++ )
                    sol[HalfNN + HalfNN - TT + i] =  static_cast<int> ( ( res[i] - 1 ) ^ ( it[HalfNN - TT + i]) );

                F2 ** tempmat = new F2 * [NN];
                for ( int i = 0; i < NN; i++ )
                    tempmat[i] = new F2[ NN ];

                for ( int i = 0; i < NN; i++ )
                    for ( int j = i; j < NN + i; j++ )
                        tempmat[i][ j % NN ] = sol[j - i];

                solution1.push_back( sol );
            }
        }

        // **********|vvvvvvvvvv|**********|vvvvvvvvvv
        for ( int x = 0; x < ( 1 << (HalfNN - TT) ); x++ ) // guess the values for the first HalfNN - TT 
        {
            // prepare the coefficients
            for ( int t = 0; t < TT; t++ )
                for ( int j = 0; j < TT; j++ )
                    tempB[t][j] = ( B[t][ HalfNN - TT + j ] + ( B[t][HalfNN + HalfNN - TT + j ] * ( it[ HalfNN - TT + j] + 1 ) ) ) % 3;

            // prepare the output1
            for ( int t = 0; t < TT; t++ )
            {
                F3 reminder = 0;

                for ( int i = 0; i < ( HalfNN - TT ); i++ )
                {
                    reminder = ( reminder + B[t][i] * ( ( x >> ( HalfNN - TT - 1 - i ) & 0x1 ) + 1 ) ) % 3; // + 1 to make it into F3*
                    reminder = ( reminder + B[t][i + HalfNN] * ( ( x >> ( HalfNN - TT - 1 - i ) & 0x1 ) + 1 ) * ( it[i] + 1 ) ) % 3; // F3* is multiplicative group
                }

                tempO[t] = mod3 ( Output2[t] - reminder ); 
            }

            auto res2 = solveF3( tempB, tempO, TT, TT ); 

            if ( res2.size() == 0 )
                continue;

            for ( auto res : res2 )
            {
                bool flag = true;
                for ( int i = 0; i < TT; i++ )
                {
                    if ( res[i] == 0 )
                    {
                        flag = false; 
                        break;
                    }
                }

                if ( flag == false )
                    continue;

                vector<F2> sol ( NN );
                for ( int i = 0; i < ( HalfNN - TT ); i++ )
                    sol[i] = x >> ( HalfNN - TT - 1 - i ) & 0x1;

                for ( int i = 0; i < TT; i++ )
                    sol[HalfNN - TT + i] = ( res[i] - 1 );

                for ( int i = 0; i < ( HalfNN - TT ); i++ )
                    sol[HalfNN + i] = ( ( x >> ( HalfNN - TT - 1 - i ) & 0x1 ) ^ it[i] );

                for ( int i = 0; i < TT; i++ )
                    sol[HalfNN + HalfNN - TT + i] =  static_cast<int> ( ( res[i] - 1 ) ^ ( it[HalfNN - TT + i]) );

                solution2.push_back( sol );
            }
        }
        // generate the linear equations
    }

    cout << solution1.size() << " " << solution2.size() << endl;

    vector<F2> realmid1 ( mid1, mid1 + NN );

    if ( find( solution1.begin(), solution1.end(), realmid1 ) != solution1.end() )
        cout << "Find it in solution1 " << endl;

    vector<F2> realmid2 ( mid2, mid2 + NN );
    if ( find( solution2.begin(), solution2.end(), realmid2 ) != solution2.end() )
        cout << "Find it in solution2 " << endl;

    //for ( auto it : solution2 )
    //{
    //    for ( int i = 0; i < NN; i++ )
    //        cout << it[i];
    //    cout << endl;
    //}

    //cout << "Real " << endl;

    //printVec( mid2, NN );


    F2 ** M = new F2 * [ 2 * NN ];
    for ( int i = 0; i <  2 * NN; i++ )
        M[i] = new F2[NN];

    for ( int row = 0; row < NN; row++ )
        for ( int col = row; col < NN + row; col++ )
            M[row][col % NN] = Input1[col - row];

    for ( int row = 0; row < NN; row++ )
        for ( int col = row; col < NN + row; col++ )
            M[NN + row][col % NN] = Input2[col - row ];

    F2 * mid = new F2[ 2 * NN ];

    for ( auto it1 : solution1 )
        for ( auto it2 : solution2 )
        {
            for ( int i = 0; i < NN; i++ )
                mid[i] = it1[i];
            for ( int i = 0; i < NN; i++ )
                mid[NN + i] = it2[i];

            auto res = solveF2( M, mid, 2 * NN, NN  );

            if ( res.size() > 0 )
            {
                cout << res.size() << endl;

                for ( auto it : res )
                    if ( EQ( it, KEY[0], NN ) ) 
                    {
                        cout << "Find it " << endl;

                        int ** KK = new F2 * [NN];
                        for ( int i = 0; i < NN; i++ )
                            KK[i] = new F2 [ NN ];

                        for ( int i = 0; i < NN; i++ )
                            for ( int j = 0; j < NN; j++ ) 
                                KK[i][ j ] = it[ ( i + j ) % NN ];

                        F3 * O1 = new F3 [ TT ];
                        F3 * O2 = new F3 [ TT ];

                        wPRF( KK, B, Input1, O1 );
                        printVec( O1, TT );

                        wPRF( KK, B, Input2, O2 );
                        printVec( O2, TT );
                    }
            }
        }
}
