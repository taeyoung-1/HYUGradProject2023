#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <unordered_map>
#include <utility>
using namespace std;
#define MAX 2147483646

int t_size;
int countdone;
vector<char> alphabet;
unordered_map<string, pair<vector<int>, vector<int>>> precomputed_values;

// 가장 기본적인 distance 계산 방식. 모든 셀 하나씩 계산
int compute_basic(const string& T, const string& P);
// k 영역 내의 모든 셀에 대해 하나씩 distance 계산.
int compute_k_diff(int _k, const string& T, const string& P);

void debugTable(const vector<vector<int>> table) {
  for (auto v : table) {
    for (auto i : v)
      cout << i << " ";
    cout << endl;
  }
}

void debugVector(const vector<int> vec) {
  for (auto i : vec)
    cout << i << " ";
}

pair<vector<int>, vector<int>> PrecomputeSingle(
    const string &row_string, // length t - 1 string
    const string &column_string, // length t - 1 string
    const vector<int> &row_offset_vector, // length t - 1 offset vector
    const vector<int> &column_offset_vector) { // length t - 1 offset vector
  // table initialization
  vector<vector<int>> table(t_size, vector<int>(t_size));
  table[0][0] = 0;
  for (int i = 1; i < t_size; i++) {
    table[0][i] = row_offset_vector[i - 1] + table[0][i - 1];
    table[i][0] = column_offset_vector[i - 1] + table[i - 1][0];
  }

  // table calculation
  for (int row = 1; row < t_size; row++) {
    for (int col = 1; col < t_size; col++) {
      int t = (row_string[col - 1] == column_string[row - 1]) ? 0 : 1;
      int diagonal = table[row - 1][col - 1] + t;
      int vertical = table[row - 1][col] + 1;
      int horizontal = table[row][col - 1] + 1;
      table[row][col] = min(diagonal, min(vertical, horizontal));
    }
  }

  // returning the result row and column offset vectors
  vector<int> row_offset_vector_output(t_size - 1);
  vector<int> column_offset_vector_output(t_size - 1);

  for (int i = 0; i < t_size - 1; i++) {
    row_offset_vector_output[i] = table[t_size - 1][i + 1] - table[t_size - 1][i];
    column_offset_vector_output[i] = table[i + 1][t_size - 1] - table[i][t_size - 1];
  }
  countdone++;

  string key = row_string + column_string;
  for (auto v: row_offset_vector) key.append(to_string(v));
  for (auto v: column_offset_vector) key.append(to_string(v));
  precomputed_values.insert(make_pair(key, make_pair(row_offset_vector_output, column_offset_vector_output)));

  // debug
  // cout << key << " ";
  // for (auto v: row_offset_vector_output) cout << to_string(v);
  // cout << " ";
  // for (auto v: column_offset_vector_output) cout << to_string(v);
  // cout << endl;

  return make_pair(row_offset_vector_output, column_offset_vector_output);
}

void PossibleOffsets(const string &rowstr, const string &colstr,
                    vector<int> row, vector<int> col,
                    int kr, int kc) {
  if (kr == 0 && kc == 0) {
    PrecomputeSingle(rowstr, colstr, row, col);
  } else if (kr == 0) {
    for (auto i : {-1, 0, 1}) {
      col.push_back(i);
      PossibleOffsets(rowstr, colstr, row, col, kr, kc - 1);
      col.pop_back();
    }
  } else {
    for (auto i : {-1, 0, 1}) {
      row.push_back(i);
      PossibleOffsets(rowstr, colstr, row, col, kr - 1, kc);
      row.pop_back();
    }
  }
}

void PossibleStringsOffsets(string rowstr, string colstr,
                            int kr, int kc) {
  if (kr == 0 && kc == 0) {
    PossibleOffsets(rowstr, colstr, {}, {}, t_size - 1, t_size - 1);
  } else if (kr == 0) {
    for (auto c : alphabet) PossibleStringsOffsets(rowstr, colstr + c, kr, kc - 1);
  } else {
    for (auto c : alphabet) PossibleStringsOffsets(rowstr + c, colstr, kr - 1, kc);
  }
}

pair<vector<int>, vector<int>> getPair(
    const string &row_string, // length t - 1 string
    const string &column_string, // length t - 1 string
    const vector<int> &row_offset_vector, // length t - 1 offset vector
    const vector<int> &column_offset_vector) { // length t - 1 offset vector
  string key = row_string + column_string;
  for (auto v : row_offset_vector) key.append(to_string(v));
  for (auto v : column_offset_vector) key.append(to_string(v));
  return precomputed_values.find(key)->second;
}

vector<int> flatten(const vector<vector<int>> &vec) {
  vector<int> flat_vec;
  flat_vec.reserve(vec.size()*vec[0].size());
  for(auto& v : vec)
    flat_vec.insert(flat_vec.end(),v.begin(),v.end());
  return flat_vec;
}

int compute_basic(const string& T, const string& P) {
  vector<int> table;
  table.reserve(T.length()+1);
  int dia = 0;

  // 초기 값 설정
  for(int col=0; col<T.length()+1; col++) 
    table[col] = col;
  
  // 이전 row의 값을 통해 현재 row의 값 계산
  for(int row=1; row<P.length()+1; row++) {
    dia = table[0]; 
    table[0] = row;
    
    for(int col=1; col<T.length()+1; col++) {
      int t = (T[col-1] == P[row-1]) ? 0 : 1;
      int diagonal = dia + t;
      int vertical = table[col] + 1;
      int horizontal = table[col-1] + 1;

      dia = table[col];
      table[col] = min(diagonal, min(vertical, horizontal));
    }
  }

  return table[T.length()];
}

int compute_k_diff(int _k, const string& T, const string& P) {
  vector<int> table;
  table.reserve(2*_k + 2);
  
  // 초기 값 설정
  table[_k] = 0; // diff = 0
  table[2*_k+1] = MAX;
  for(int col=1; col<_k+1; col++) {
    table[_k-col] = col;  // diff = -1 ~ -k
    table[_k+col] = col;  // diff = 1 ~ k
  }
  
  // 이전 row의 값을 통해 현재 row의 값 계산
  for(int row=1; row<P.length()+1; row++) {
    int col_0 = row - _k;
    if (col_0 > 0) { // table[0] 값 업데이트
      int t = (T[col_0-1] == P[row-1]) ? 0 : 1;
      table[0] = min(table[0]+t, table[1]+1);
    }
    for(int col=max(1,_k-row+1); col<2*_k+1; col++) {
      if(row+col-_k > T.length()) {
        break;
      }
      int t = (T[row+col-_k-1] == P[row-1]) ? 0 : 1;
      int diagonal = table[col] + t;
      int vertical = table[col+1] + 1;
      int horizontal = table[col-1] + 1;

      table[col] = min(diagonal, min(vertical, horizontal));
      
    }
  }

  return table[T.length()-P.length()+_k];
}
vector<int> compute_distance(vector<int> row_vec, const vector<int> &col_vec, const string &T, const string &P) {
  // table initialization
  vector<vector<int>> table(P.length()+1, vector<int>(T.length()+1));
  table[0][0] = 0;
  for (int i = 1; i<T.length()+1; i++)
    table[0][i] = row_vec[i - 1] + table[0][i - 1];
  for (int i = 1; i<P.length()+1; i++)
    table[i][0] = col_vec[i - 1] + table[i - 1][0];
  
  // table calculation
  for(int row = 1; row<P.length()+1; row++) {
    for(int col = 1; col<T.length()+1; col++) {
      int t = (T[col-1] == P[row-1]) ? 0 : 1;
      int diagonal = table[row-1][col-1] + t;
      int vertical = table[row-1][col] + 1;
      int horizontal = table[row][col-1] + 1;
      table[row][col] = min(diagonal, min(vertical, horizontal));
    }
  }
  
  // remain_vector
  for (int col = 0; col < T.length(); col++) 
    row_vec[col] = table[P.length()][col+1] - table[P.length()][col];
  
  return row_vec;
}

// Assumption: n = n'(k - 1), m = m'(k - 1) for positives n' m'.
void overAllTest(const string& T, const string& P) {
  int accumulation = 0;
  int m = P.length();
  int n = T.length();
  int row_blocks = m / (t_size - 1);
  int col_blocks = n / (t_size - 1);

  vector<vector<int>> prev_row_offset_vec(col_blocks, vector<int>(t_size-1,1)); // previous row
  vector<vector<int>> last_col_offset_vec(row_blocks, vector<int>(t_size-1,1)); // last column
  
  for (int row = 0; row < row_blocks; row++) {
    vector<int> col_vec = last_col_offset_vec[row];
    // col < wing+col_blocks - row 까지만 loop
    for (int col = 0; col < col_blocks; col++) {
      auto pair = getPair(T.substr(col * (t_size - 1), t_size - 1),
                          P.substr(row * (t_size - 1), t_size - 1),
                          prev_row_offset_vec[col],
                          col_vec);

      prev_row_offset_vec[col] = pair.first;
      col_vec = pair.second;
    }
    last_col_offset_vec[row] = col_vec;
  }
  accumulation += row_blocks * (t_size-1);
  
  // 나머지가 존재하는 경우
  int col_remain = m%(t_size-1);
  int row_remain = n%(t_size-1);

  vector<int> col_remain_vector(col_remain, 1);
  vector<int> row_remain_vector(row_remain, 1);
  
  // prev_row_offset_vec, last_col_offset_vec를 1차원 벡터로 변환
  vector<int> row_prev_vector = flatten(prev_row_offset_vec);
  vector<int> col_last_vector = flatten(last_col_offset_vec);

  string sub_p = P.substr(row_blocks*(t_size-1),col_remain);
  string sub_t = T.substr(col_blocks*(t_size-1),row_remain);
  string SUB_t = T.substr(0,col_blocks*(t_size-1));

  if (row_remain) {
    string SUB_p = P.substr(0,row_blocks*(t_size-1));
    
    row_remain_vector = compute_distance(row_remain_vector,col_last_vector, sub_t, SUB_p);
    
    SUB_t += sub_t;
    row_prev_vector.insert(row_prev_vector.end(),row_remain_vector.begin(),row_remain_vector.end());
  }

  if (col_remain) {
    col_remain_vector = compute_distance(col_remain_vector,row_prev_vector,sub_p,SUB_t);
    for (auto v : col_remain_vector) accumulation += v;
  }

  for (auto v : row_prev_vector) accumulation += v;

  cout << "edit distance: " << accumulation << endl;
}

int main(void) {
  int alphabet_size;
  cout << "alphabet set size: ";
  cin >> alphabet_size;

  cout << "alphabets: ";
  for (int i = 0; i < alphabet_size; ++i) {
    char c;
    cin >> c;
    alphabet.push_back(c);
  }

  cout << "t-block size: ";
  cin >> t_size;

  PossibleStringsOffsets("","",t_size-1,t_size-1);
  while(1) {
	  string T;
	  string P;
	  cout << "enter text: ";
	  cin >> T;
	  cout << "enter pattern: ";
	  cin >> P;
	  if(T == "q") return 0;
	  overAllTest(T,P);

  }
  // cout << countdone << endl;
  // cout << precomputed_values.size() << endl;
  // auto pair = (precomputed_values.find("GGGGGG111111")->second);
  return 0;
}