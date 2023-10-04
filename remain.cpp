#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <unordered_map>
#include <utility>
using namespace std;

int t_size;
vector<char> alphabet;
unordered_map<string, pair<vector<int>, vector<int>>> precomputed_values;

int countdone;

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
  for(auto& v : vec)
    flat_vec.insert(flat_vec.end(),v.begin(),v.end());
  return flat_vec;
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
  vector<int> row_offset_vector_output(t_size - 1);
  
  for (int col = 0; col < T.length(); col++) 
    row_vec[col] = table[P.length()][col+1] - table[P.length()][col];
  
  return row_vec;
}

// Assumption: n = n'(k - 1), m = m'(k - 1) for positives n' m'.
void overAllTest(int k, const string& T, const string& P) {
  int wing = (k - 1) / (t_size - 1) + 1;
  cout << "wing: " << wing << endl;
  int num_t_block_per_row = 2 * wing + 1;
  cout << "number of t blocks per row: " << num_t_block_per_row << endl;

  int accumulation = 0;
  int m = P.length();
  int n = T.length();
  int row_blocks = m / (t_size - 1);
  int col_blocks = n / (t_size - 1);

  // 저장해야 하는 col vector 개수
  // 각 row의 마지막 t-block의 (row,col)은 기본적으로 (row, 2*wing)
        // 2*wing == num_t_block_per_row - 1
  // 마지막 row에서는 (row_blocks-1, ((col_blocks-1) - ((row_blocks-1)-wing)))
  // 따라서 ((col_blocks-1) - ((row_blocks-1)-wing))와 (2*wing) 사이의 개수만큼 필요하다
  int num_prev_col_blocks = wing - col_blocks + row_blocks + 1;

  cout << "row_blocks: " << row_blocks << ", col_blocks: " << col_blocks << endl;

  vector<int> temp_col_offset_vec(t_size-1); // temporary column vector
  vector<vector<int>> prev_row_offset_vec(num_t_block_per_row, vector<int>(t_size-1,1)); // previous row
  vector<vector<int>> prev_col_offset_vec(num_prev_col_blocks); // last column
  fill(prev_row_offset_vec.begin(), prev_row_offset_vec.begin()+1, vector<int>(t_size-1, 0));
  
  for (int row = 0; row < row_blocks; row++) {
    fill(temp_col_offset_vec.begin(), temp_col_offset_vec.end(), 1);
    
    // col < wing+col_blocks - row 까지만 loop
    for (int col = max(0,wing-row); col < num_t_block_per_row; col++) {
      auto pair = getPair(T.substr((row + col - wing) * (t_size - 1), t_size - 1),
                          P.substr(row * (t_size - 1), t_size - 1),
                          (col == num_t_block_per_row - 1) ?
                            vector<int>(t_size - 1, 1) : prev_row_offset_vec[col + 1],
                          temp_col_offset_vec);

      prev_row_offset_vec[col] = pair.first;
      temp_col_offset_vec = pair.second;

      // row의 마지막 t-block이 text까지의 마지막 t-block 영역인 경우 column vector 저장
      if(row + col == wing + col_blocks - 1) {
        prev_col_offset_vec[row-(row_blocks-num_prev_col_blocks)] = pair.second;
        break;
      }
    }
    
    accumulation += t_size-1;
    for (auto v : prev_row_offset_vec[0]) accumulation += v;
  }

  // 마지막 t-block 마지막 셀의 edit distance
  for (int col=0; col<col_blocks-row_blocks+wing; col++) {
    for (auto v : prev_row_offset_vec[col+1]) accumulation += v;
  }

  // debug. 마지막 t-block 까지의 edit distance
  // cout << "distance: " << accumulation << endl;

  // 나머지가 존재하는 경우
  int col_remain = m%(t_size-1);
  int row_remain = n%(t_size-1);

  vector<int> col_remain_vector(col_remain, 1);
  vector<int> row_remain_vector(row_remain, 1);

  string sub_p = P.substr(row_blocks*(t_size-1),col_remain);
  string sub_t = T.substr(col_blocks*(t_size-1),row_remain);
  
  if (col_remain) { // pattern 나머지가 존재
    string SUB_t = T.substr((row_blocks-wing-1)*(t_size-1),(col_blocks-row_blocks+wing+1)*(t_size-1));
    vector<int> row_prev_vector = flatten(prev_row_offset_vec);
    
    col_remain_vector = compute_distance(col_remain_vector,row_prev_vector,sub_p,SUB_t);
  }
  
  if (row_remain) { // text 나머지가 존재
    string SUB_p = P.substr((row_blocks-num_prev_col_blocks)*(t_size-1),(num_prev_col_blocks*(t_size-1)));
    vector<int> col_prev_vector = flatten(prev_col_offset_vec);
    
    row_remain_vector = compute_distance(row_remain_vector,col_prev_vector, sub_t, SUB_p);
  }

  if(col_remain && row_remain) {
    row_remain_vector = compute_distance(row_remain_vector,col_remain_vector, sub_t, sub_p);   
    
    for (auto v : col_remain_vector) accumulation += v;
    for (auto v : row_remain_vector) accumulation += v;
  }
  else if(col_remain)
    for (auto v : col_remain_vector) accumulation += v;
  else if(row_remain)
    for (auto v : row_remain_vector) accumulation += v;

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
	  overAllTest(6,T,P);

  }
  // cout << countdone << endl;
  // cout << precomputed_values.size() << endl;
  // auto pair = (precomputed_values.find("GGGGGG111111")->second);
  return 0;
}
