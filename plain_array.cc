#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <fstream>
#include <chrono>

using namespace std;
using namespace chrono;
#define MAX 2147483646
#define NUM_CASE 50 // 계산하는 case의 개수
#define T_SIZE 4    // t-block의 크기
#define K_SIZE 7         // k 영역의 크기

vector<char> alphabet = {'A', 'T', 'G', 'C'};
unordered_map<string, pair<vector<int>, vector<int>>> precomputed_values;
pair<vector<int>, vector<int>> precomputed_alt[2985984];

// 가장 기본적인 distance 계산 방식. 모든 셀 하나씩 계산
int compute_basic(const string& T, const string& P);
// k 영역 내의 모든 셀에 대해 하나씩 distance 계산.
int compute_k_diff(const string& T, const string& P);
// t-block 단위로 distance 계산
int compute_russian(const string& T, const string& P);

int convertChar(char c) {
  switch (c)
  {
  case 'A':
    return 0;
    break;
  case 'T':
    return 1;
    break;
  case 'G':
    return 2;
    break;
  case 'C':
    return 3;
    break;
  default:
    cout << "Unknown character";
    exit(0);
  }
}

inline int computeIndex(
    const string &row_string, // length t - 1 string
    const string &column_string, // length t - 1 string
    const vector<int> &row_offset_vector, // length t - 1 offset vector
    const vector<int> &column_offset_vector) {
  int index = 0;
  index += 1 * (column_offset_vector[2] + 1) +
           3 * (column_offset_vector[1] + 1) +
           9 * (column_offset_vector[0] + 1);
  index += 27 * (row_offset_vector[2] + 1) +
           81 * (row_offset_vector[1] + 1) +
           243 * (row_offset_vector[0] + 1);
  index += 729 * (convertChar(column_string[2])) +
           4 * 729 * (convertChar(column_string[1])) +
           16 * 729 * (convertChar(column_string[0]));
  index += 64 * 729 * (convertChar(row_string[2])) +
           256 * 729 * (convertChar(row_string[1])) +
           1024 * 729 * (convertChar(row_string[0]));
  return index;
}

void PrecomputeSingle(
    const string &row_string, // length t - 1 string
    const string &column_string, // length t - 1 string
    const vector<int> &row_offset_vector, // length t - 1 offset vector
    const vector<int> &column_offset_vector) { // length t - 1 offset vector
  // table initialization
  vector<vector<int>> table(T_SIZE, vector<int>(T_SIZE));
  table[0][0] = 0;
  for (int i = 1; i < T_SIZE; i++) {
    table[0][i] = row_offset_vector[i - 1] + table[0][i - 1];
    table[i][0] = column_offset_vector[i - 1] + table[i - 1][0];
  }

  // table calculation
  for (int row = 1; row < T_SIZE; row++) {
    for (int col = 1; col < T_SIZE; col++) {
      int t = (row_string[col - 1] == column_string[row - 1]) ? 0 : 1;
      int diagonal = table[row - 1][col - 1] + t;
      int vertical = table[row - 1][col] + 1;
      int horizontal = table[row][col - 1] + 1;
      table[row][col] = min(diagonal, min(vertical, horizontal));
    }
  }

  // returning the result row and column offset vectors
  vector<int> row_offset_vector_output(T_SIZE - 1);
  vector<int> column_offset_vector_output(T_SIZE - 1);

  for (int i = 0; i < T_SIZE - 1; i++) {
    row_offset_vector_output[i] = table[T_SIZE - 1][i + 1] - table[T_SIZE - 1][i];
    column_offset_vector_output[i] = table[i + 1][T_SIZE - 1] - table[i][T_SIZE - 1];
  }

  precomputed_alt[computeIndex(row_string, column_string, row_offset_vector, column_offset_vector)]
    = make_pair(row_offset_vector_output, column_offset_vector_output);
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
    PossibleOffsets(rowstr, colstr, {}, {}, T_SIZE - 1, T_SIZE - 1);
  } else if (kr == 0) {
    for (auto c : alphabet) PossibleStringsOffsets(rowstr, colstr + c, kr, kc - 1);
  } else {
    for (auto c : alphabet) PossibleStringsOffsets(rowstr + c, colstr, kr - 1, kc);
  }
}

inline pair<vector<int>, vector<int>> getPair(
    const string &row_string, // length t - 1 string
    const string &column_string, // length t - 1 string
    const vector<int> &row_offset_vector, // length t - 1 offset vector
    const vector<int> &column_offset_vector) { // length t - 1 offset vector
  return precomputed_alt[computeIndex(row_string, column_string, row_offset_vector, column_offset_vector)];
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

int compute_k_diff(const string& T, const string& P) {
  vector<int> table;
  table.reserve(2*K_SIZE + 2);
  
  // 초기 값 설정
  table[K_SIZE] = 0; // diff = 0
  table[2*K_SIZE+1] = MAX;
  for(int col=1; col<K_SIZE+1; col++) {
    table[K_SIZE-col] = col;  // diff = -1 ~ -k
    table[K_SIZE+col] = col;  // diff = 1 ~ k
  }
  
  // 이전 row의 값을 통해 현재 row의 값 계산
  for(int row=1; row<P.length()+1; row++) {
    int col_0 = row - K_SIZE;
    if (col_0 > 0) { // table[0] 값 업데이트
      int t = (T[col_0-1] == P[row-1]) ? 0 : 1;
      table[0] = min(table[0]+t, table[1]+1);
    }
    for(int col=max(1,K_SIZE-row+1); col<2*K_SIZE+1; col++) {
      if(row+col-K_SIZE > T.length()) {
        break;
      }
      int t = (T[row+col-K_SIZE-1] == P[row-1]) ? 0 : 1;
      int diagonal = table[col] + t;
      int vertical = table[col+1] + 1;
      int horizontal = table[col-1] + 1;

      table[col] = min(diagonal, min(vertical, horizontal));
      
    }
  }

  return table[T.length()-P.length()+K_SIZE];
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
int compute_russian(const string& T, const string& P) {
  int accumulation = 0;
  int m = P.length()-1;
  int n = T.length()-1;
  int row_blocks = m / (T_SIZE - 1);
  int col_blocks = n / (T_SIZE - 1);

  vector<vector<int>> prev_row_offset_vec(col_blocks, vector<int>(T_SIZE-1,1)); // previous row
  vector<vector<int>> last_col_offset_vec(row_blocks, vector<int>(T_SIZE-1,1)); // last column
  for (int row = 0; row < row_blocks; row++) {
    vector<int> col_vec = last_col_offset_vec[row];
    // col < wing+col_blocks - row 까지만 loop
    for (int col = 0; col < col_blocks; col++) {
      auto pair = getPair(T.substr(col * (T_SIZE - 1), T_SIZE - 1),
                          P.substr(row * (T_SIZE - 1), T_SIZE - 1),
                          prev_row_offset_vec[col],
                          col_vec);

      prev_row_offset_vec[col] = pair.first;
      col_vec = pair.second;
    }
    last_col_offset_vec[row] = col_vec;
  }
  accumulation += row_blocks * (T_SIZE-1);
  
  // 나머지가 존재하는 경우
  int col_remain = m%(T_SIZE-1);
  int row_remain = n%(T_SIZE-1);

  vector<int> col_remain_vector(col_remain, 1);
  vector<int> row_remain_vector(row_remain, 1);
  
  // prev_row_offset_vec, last_col_offset_vec를 1차원 벡터로 변환
  vector<int> row_prev_vector = flatten(prev_row_offset_vec);
  vector<int> col_last_vector = flatten(last_col_offset_vec);

  string sub_p = P.substr(row_blocks*(T_SIZE-1),col_remain);
  string sub_t = T.substr(col_blocks*(T_SIZE-1),row_remain);
  string SUB_t = T.substr(0,col_blocks*(T_SIZE-1));

  if (row_remain) {
    string SUB_p = P.substr(0,row_blocks*(T_SIZE-1));
    
    row_remain_vector = compute_distance(row_remain_vector,col_last_vector, sub_t, SUB_p);
    
    SUB_t += sub_t;
    row_prev_vector.insert(row_prev_vector.end(),row_remain_vector.begin(),row_remain_vector.end());
  }

  if (col_remain) {
    col_remain_vector = compute_distance(col_remain_vector,row_prev_vector,sub_p,SUB_t);
    for (auto v : col_remain_vector) accumulation += v;
  }

  for (auto v : row_prev_vector) accumulation += v;

  return accumulation;
}

int compute_k_and_russian(const string& T, const string& P) {
  int wing = (K_SIZE - 1) / (T_SIZE - 1) + 1;
  int num_t_block_per_row = 2 * wing + 1;
  int accumulation = 0;
  int m = P.length()-1;
  int n = T.length()-1;
  int row_blocks = m / (T_SIZE - 1);
  int col_blocks = n / (T_SIZE - 1);

  if (n - m > K_SIZE) return -1;

  vector<int> temp_col_offset_vec(T_SIZE-1); // temporary column vector
  vector<vector<int>> prev_row_offset_vec(num_t_block_per_row, vector<int>(T_SIZE-1,1)); // previous row
  
  // text까지 t_block보다 k영역이 클 때, 마지막 t-block의 column 벡터 저장
    // 저장해야 하는 벡터 개수 == wing + 1 - col_blocks + row_blocks, (0 <= col_blocks - row_blocks <= wing);
  vector<int> prev_col_offset_vec;
  prev_col_offset_vec.reserve((wing+1)*(T_SIZE-1)); 
  fill(prev_row_offset_vec.begin(), prev_row_offset_vec.begin()+1, vector<int>(T_SIZE-1, 0));
  
  for (int row = 0; row < row_blocks; row++) {
    fill(temp_col_offset_vec.begin(), temp_col_offset_vec.end(), 1);
    
    // col < wing+col_blocks - row 까지만 loop
    for (int col = max(0,wing-row); col < num_t_block_per_row; col++) {
      auto pair = getPair(T.substr((row + col - wing) * (T_SIZE - 1), T_SIZE - 1),
                          P.substr(row * (T_SIZE - 1), T_SIZE - 1),
                          (col == num_t_block_per_row - 1) ?
                            vector<int>(T_SIZE - 1, 1) : prev_row_offset_vec[col + 1],
                          temp_col_offset_vec);

      prev_row_offset_vec[col] = pair.first;
      temp_col_offset_vec = pair.second;

      // row의 마지막 t-block이 text까지의 마지막 t-block 영역인 경우 column vector 저장
      if(row + col == wing + col_blocks - 1) {
        prev_col_offset_vec.insert(prev_col_offset_vec.end(),pair.second.begin(),pair.second.end());
        break;
      }
    }
    accumulation += T_SIZE-1;
    for (auto v : prev_row_offset_vec[0]) accumulation += v;
  }
  // debug. 마지막 t-block 까지의 edit distance
  // cout << "distance: " << accumulation << endl;

  // 나머지가 존재하는 경우
  int col_remain = m%(T_SIZE-1);
  int row_remain = n%(T_SIZE-1);

  vector<int> col_remain_vector(col_remain, 1);
  vector<int> row_remain_vector(row_remain, 1);
  
  // prev_row_offset_vec을 1차원 벡터로 변환
  vector<int> row_prev_vector;
  row_prev_vector.reserve((col_blocks-row_blocks+wing+2)*(T_SIZE-1));
  for(int i=1; i<col_blocks-row_blocks+wing+1;i++) {
    row_prev_vector.insert(row_prev_vector.end(),prev_row_offset_vec[i].begin(),prev_row_offset_vec[i].end());
  }

  string sub_p = P.substr(row_blocks*(T_SIZE-1),col_remain);
  string sub_t = T.substr(col_blocks*(T_SIZE-1),row_remain);
  string SUB_t = T.substr((row_blocks-wing)*(T_SIZE-1),(col_blocks-row_blocks+wing)*(T_SIZE-1));
             
  if (row_remain) {
    string SUB_p = P.substr(P.length()-prev_col_offset_vec.size()-col_remain,prev_col_offset_vec.size());
    
    row_remain_vector = compute_distance(row_remain_vector,prev_col_offset_vec, sub_t, SUB_p);
    
    SUB_t += sub_t;
    row_prev_vector.insert(row_prev_vector.end(),row_remain_vector.begin(),row_remain_vector.end());
  }

  if (col_remain) {
    col_remain_vector = compute_distance(col_remain_vector,row_prev_vector,sub_p,SUB_t);
    for (auto v : col_remain_vector) accumulation += v;
  }

  for (auto v : row_prev_vector) accumulation += v;

  return accumulation;
}

vector<string> read_file(const string& filename) {
  vector<string> output;
  output.reserve(NUM_CASE);
  string line;
  
  ifstream in_file(filename);
  if(!in_file.is_open()) {
    cerr << "Error opening file '" << filename << "'" << endl;
    exit(1);
  }

  while(getline(in_file, line))
    output.push_back(line);
  in_file.close();
  return output;
}

int main(void) {
  string text_filename("text.txt");
  string pattern_filename("pattern.txt");
  string out_filename("result.txt");

  // precomputing
  system_clock::time_point start_pre = system_clock::now();
  PossibleStringsOffsets("","",T_SIZE-1,T_SIZE-1);
  duration<double, ratio<1,1000>> time_precompute = system_clock::now() - start_pre;
  
  // testcase
  vector<string> Text = read_file(text_filename);
  vector<string> Pattern = read_file(pattern_filename);

  int accumulation = 0;
  string T;
	string P;
  
  // result 파일 쓰기
  ofstream out_file(out_filename);
  if(!out_file.is_open()) {
    cerr << "Could not open file '" << out_filename << "'" << endl;
    return EXIT_FAILURE;
  }

  // basic method
  system_clock::time_point start_basic = system_clock::now();
  out_file << "----------basic method----------" << endl;
  for(int i=1; i<Text.size(); i++) {
    accumulation = compute_basic(Text[i],Pattern[i]);
    out_file << accumulation << "\t";  
  }
  duration<double, milli> time_basic = system_clock::now() - start_basic;
  out_file << endl;
  out_file << "cost time: " << time_basic.count() << "msec" << endl;
  
  // k_difference method
  system_clock::time_point start_k = system_clock::now();
  out_file << "----------K difference method----------" << endl;
  for(int i=1; i<Text.size(); i++) {
    accumulation = compute_k_diff(Text[i],Pattern[i]);
    out_file << accumulation << "\t";  
  }
  duration<double, milli> time_k = system_clock::now() - start_k;
  out_file << endl;
  out_file << "cost time: " << time_k.count() << "msec" << endl;
  
  // four_russians method
  system_clock::time_point start_russian = system_clock::now();
  out_file << "----------Four Russians method----------" << endl;
  out_file << "precomputing time: " << time_precompute.count() << "msec" << endl;
  for(int i=1; i<Text.size(); i++) {
    accumulation = compute_russian(Text[i],Pattern[i]);
    out_file << accumulation << "\t";  
  }
  duration<double, milli> time_russian = system_clock::now() - start_russian;
  out_file << endl;
  out_file << "cost time: " << time_russian.count() << "msec" << endl;

  // four_russians and k_difference method
  system_clock::time_point start_k_russian = system_clock::now();
  out_file << "----------Four Russians and K difference method----------" << endl;
  out_file << "precomputing time: " << time_precompute.count() << "msec" << endl;
  for(int i=1; i<Text.size(); i++) {
    if(Text[i].length() < Pattern[i].length())
      accumulation = compute_k_and_russian(Pattern[i],Text[i]);
    else
      accumulation = compute_k_and_russian(Text[i],Pattern[i]);
    out_file << accumulation << "\t";  
  }
  duration<double, milli> time_k_russian = system_clock::now() - start_k_russian;
  out_file << endl;
  out_file << "cost time: " << time_k_russian.count() << "msec" << endl;

  
  out_file.close();

  return 0;
}
