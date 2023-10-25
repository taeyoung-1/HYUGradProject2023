#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <fstream>
#include <chrono>
#include <cmath>

using namespace std;
using namespace chrono;
#define MAX 2147483646
#define NUM_CASE 50 // 계산하는 case의 개수
#define T_SIZE 5    // t-block의 크기
#define K_SIZE 7         // k 영역의 크기
#define TEXT_LAST pow(4,T_SIZE-2)-1
#define OFFSET_ALL_1 pow(16,T_SIZE-1)-1
#define VECTOR_ALL_1 pow(4,T_SIZE-1)-1
#define VECTOR_ALL_0 (pow(4,T_SIZE-1)-1) / 3 * 2

unordered_map<long int, int> precomputed_values;

// // 두 문자열 T와 P 사이의 distance를 계산하는 함수들
// 가장 기본적인 distance 계산 방식. 모든 셀 하나씩 계산
int compute_basic(const string& T, const string& P);
// k 영역 내의 모든 셀에 대해 하나씩 distance 계산.
int compute_k_diff(const string& T, const string& P);
// t-block 단위로 distance 계산
int compute_russian(const string& T, const string& P);
// k 영역 내에서 t-block 단위로 계산할 때 나머지 부분을 계산하기 위한 함수 
vector<int> compute_distance(const vector<int> &row_vec, const vector<int> &col_vec, const string &T, const string &P);
// k 영역 내에서 t-block 단위로 distance 계산
int compute_k_and_russian(const string& T, const string& P);

// // 타입 변환 함수
// 두 개의 문자열에서 alphabet이 같은 부분의 위치를 long int 형태로 저장
long int stringToLong(const string &_text, const string &_pattern) {
  long int result = 0;
  for (int row=0; row<T_SIZE-1; ++row) {
    for (int col=0; col<T_SIZE-1; ++col) {
      if (_text[col] == _pattern[row]) result = (result << 1) + 1;
      else result = (result << 1);
    }
  }
  return result;
}

int vectorToInt(vector<int> offsetVector) {
  int result;
  for (int i = 0; i < offsetVector.size(); ++i) {
    result += offsetVector[i] + 2;
    result = result << 2;
  }
  return result;
}
vector<int> intToVector(int offset_vector) {
  vector<int> result_vector;
  result_vector.reserve(T_SIZE-1);
  for (int i=0; i<T_SIZE-1; ++i)
    result_vector.push_back(((offset_vector >> 2*(T_SIZE-2 - i)) & 3) - 2);
  return result_vector;
}

void PrecomputeSingle(
    long int index, // length t - 1 string
    int row_offset_vector, // length t - 1 offset vector
    int col_offset_vector) { // length t - 1 offset vector
  // table initialization
  vector<vector<int>> table(T_SIZE, vector<int>(T_SIZE));
  table[0][0] = 0;
  int row_offset, col_offset;
  
  for (int i = 1; i < T_SIZE; i++) {
    // case '1': -1, case '2': 0, case '3': +1 
    row_offset = ((row_offset_vector >> 2*(T_SIZE-1 - i)) & 3) - 2;
    col_offset = ((col_offset_vector >> 2*(T_SIZE-1 - i)) & 3) - 2;
    
    table[0][i] = row_offset + table[0][i - 1];
    table[i][0] = col_offset + table[i - 1][0];
  }

  // table calculation
  for (int row = 1; row < T_SIZE; row++) {
    for (int col = 1; col < T_SIZE; col++) {
      int t = (index >> ((T_SIZE-1)*(T_SIZE-1-row)+(T_SIZE-1-col))) & 1;
      int diagonal = table[row - 1][col - 1] + (1-t);
      int vertical = table[row - 1][col] + 1;
      int horizontal = table[row][col - 1] + 1;
      table[row][col] = min(diagonal, min(vertical, horizontal));
    }
  }

  // returning the result row and column offset vectors
  int row_output = 0;
  int col_output = 0;

  for (int i = 0; i < T_SIZE - 1; i++) {
    row_output = (row_output << 2);
    col_output = (col_output << 2);
    row_output += table[T_SIZE - 1][i + 1] - table[T_SIZE - 1][i] + 2;
    col_output += table[i + 1][T_SIZE - 1] - table[i][T_SIZE - 1] + 2;
  }

  long int key = (index << 4*(T_SIZE-1)) + (row_offset_vector << 2*(T_SIZE-1)) + col_offset_vector;
  int value = (row_output << 2*(T_SIZE-1)) + col_output;
  precomputed_values.insert(make_pair(key, value));
}

void PossibleOffsets(long int index,
                    int rowvec, int colvec,
                    int kr, int kc) {
  if (kr == 0 && kc == 0) {
    PrecomputeSingle(index, rowvec, colvec);
  } else if (kr == 0) {
    for (auto i : {1, 2, 3}) {
      int temp = colvec << 2;
      temp += i;
      PossibleOffsets(index, rowvec, temp, kr, kc - 1);
    }
  } else {
    for (auto i : {1, 2, 3}) {
      int temp = rowvec << 2;
      temp += i;
      PossibleOffsets(index, temp, colvec, kr - 1, kc);
    }
  }
}

long int findSameIndex(int _text, int _pattern) {
  long int index = 0;
  for(int row=0; row<T_SIZE-1; ++row) {
    for(int col=0; col<T_SIZE-1; ++col) {
      if ((_text >> 2*(T_SIZE-2-col) & 3) == (_pattern >> 2*(T_SIZE-2-row) & 3))
        index = (index << 1) + 1;
      else index = (index << 1);
    }
  }
  return index;
}

void PossibleStringsOffsets() {
  vector<long int> possibleString;
  possibleString.reserve(pow(2,4*(T_SIZE-1)-5));
  for (int pattern=0; pattern<VECTOR_ALL_1; ++pattern) {
    for (int text=0; text<TEXT_LAST; ++text) {
      long int sameIndex = findSameIndex(pattern,text);
      if (find(possibleString.begin(), possibleString.end(), sameIndex) == possibleString.end())
        possibleString.push_back(sameIndex);
    }
  }
  long int count = 0;
  for (auto v : possibleString)
    PossibleOffsets(v, 0, 0, T_SIZE-1, T_SIZE-1);
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
vector<int> compute_distance(const vector<int> &row_vec, const vector<int> &col_vec, const string &T, const string &P) {
  // table initialization
  vector<int> result;
  result.reserve(row_vec.size());
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
    result.push_back(table[P.length()][col+1] - table[P.length()][col]);
  
  return result;
}

// Assumption: n = n'(k - 1), m = m'(k - 1) for positives n' m'.
int compute_russian(const string& T, const string& P) {
  int accumulation = 0;
  int m = P.length();
  int n = T.length();
  if(isspace(P[m-1])) { // 라인 마지막에 공백이나 개행이 있는 경우 제외
    m -= 1;
    n -= 1;
  }
  int row_blocks = m / (T_SIZE - 1);
  int col_blocks = n / (T_SIZE - 1);
  
  vector<int> prev_row_offset_vec(col_blocks, VECTOR_ALL_1); // previous row
  vector<int> last_col_offset_vec(row_blocks, VECTOR_ALL_1); // last column
  for (int row = 0; row < row_blocks; row++) {
    int col_vec = VECTOR_ALL_1;
    // col < wing+col_blocks - row 까지만 loop
    for (int col = 0; col < col_blocks; col++) {
      auto key = stringToLong((T.substr(col * (T_SIZE - 1), T_SIZE - 1)),
                          (P.substr(row * (T_SIZE - 1), T_SIZE - 1)));
      key = (key << 4*(T_SIZE-1)) + (prev_row_offset_vec[col] << 2*(T_SIZE-1)) + col_vec;
      auto value = precomputed_values.find(key)->second;
      
      prev_row_offset_vec[col] = (value >> 2*(T_SIZE-1));
      col_vec = value & (int)VECTOR_ALL_1;
    }
    last_col_offset_vec[row] = col_vec;
  }
  accumulation += row_blocks * (T_SIZE-1);
  
  // 나머지가 존재하는 경우
  int col_remain = m%(T_SIZE-1);
  int row_remain = n%(T_SIZE-1);

  vector<int> col_remain_vector(col_remain, 1);
  vector<int> row_remain_vector(row_remain, 1);
  vector<int> row_prev_vector;
  row_prev_vector.reserve(col_blocks);
  vector<int> col_last_vector;
  col_last_vector.reserve(row_blocks);
  vector<int> temp;

  // prev_row_offset_vec, last_col_offset_vec를 1차원 벡터로 변환
  for (auto v : prev_row_offset_vec) {
    temp = intToVector(v);
    row_prev_vector.insert(row_prev_vector.end(), temp.begin(), temp.end());
  }
  for (auto v : last_col_offset_vec) {
    temp = intToVector(v);
    col_last_vector.insert(col_last_vector.end(), temp.begin(), temp.end());
  }
  
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
  int m = P.length();
  int n = T.length();
  if(isspace(P[m-1])) { // 라인 마지막에 공백이나 개행이 있는 경우 제외
    m -= 1;
    n -= 1;
  }
  int row_blocks = m / (T_SIZE - 1);
  int col_blocks = n / (T_SIZE - 1);

  if (n - m > K_SIZE) return -1;

  int temp_col_offset_vec; // temporary column vector
  vector<int> prev_row_offset_vec(num_t_block_per_row, VECTOR_ALL_1); // previous row
  vector<int> row_zero_vector;

  // text까지 t_block보다 k영역이 클 때, 마지막 t-block의 column 벡터 저장
    // 저장해야 하는 벡터 개수 == wing + 1 - col_blocks + row_blocks, (0 <= col_blocks - row_blocks <= wing);
  vector<int> prev_col_offset_vec;
  prev_col_offset_vec.reserve(wing+1);

  // prev_row_offset_vec[0]을 0으로 세팅
  prev_row_offset_vec[0] = VECTOR_ALL_0;

  for (int row = 0; row < row_blocks; row++) {
    temp_col_offset_vec = VECTOR_ALL_1;
    
    // col < wing+col_blocks - row 까지만 loop
    for (int col = max(0,wing-row); col < num_t_block_per_row; col++) {
      auto key = stringToLong((T.substr((row + col - wing) * (T_SIZE - 1), T_SIZE - 1)),
                          (P.substr(row * (T_SIZE - 1), T_SIZE - 1)));
      key = (key << 4*(T_SIZE-1)) + temp_col_offset_vec;
      key += (((col == num_t_block_per_row - 1) ? (int)VECTOR_ALL_1 : prev_row_offset_vec[col + 1]) << 2*(T_SIZE-1));
      auto value = precomputed_values.find(key)->second;
      
      prev_row_offset_vec[col] = (value >> 2*(T_SIZE-1));
      temp_col_offset_vec = value & (int)VECTOR_ALL_1;

      // row의 마지막 t-block이 text까지의 마지막 t-block 영역인 경우 column vector 저장
      if(row + col == wing + col_blocks - 1) {
        prev_col_offset_vec.push_back(temp_col_offset_vec);
        break;
      }
    }
    
    accumulation += T_SIZE-1;
    row_zero_vector = intToVector(prev_row_offset_vec[0]);
    
    for (auto v : row_zero_vector) accumulation += v;
  }
  
  // 나머지가 존재하는 경우
  int col_remain = m%(T_SIZE-1);
  int row_remain = n%(T_SIZE-1);
  
  vector<int> col_remain_vector(col_remain, 1);
  vector<int> row_remain_vector(row_remain, 1);
  
  string sub_p = P.substr(row_blocks*(T_SIZE-1),col_remain);
  string sub_t = T.substr(col_blocks*(T_SIZE-1),row_remain);
  string SUB_t;

  // prev_row_offset_vec을 1차원 벡터로 변환
  vector<int> temp, row_prev_vector, col_last_vector;
  row_prev_vector.reserve((col_blocks-row_blocks+wing+1)*(T_SIZE-1)+row_remain);
  if (wing > row_blocks) {
    for (int i= wing-row_blocks; i<col_blocks-row_blocks+wing+1;i++) {
      temp = intToVector(prev_row_offset_vec[i]);
      row_prev_vector.insert(row_prev_vector.end(), temp.begin(), temp.end());
    }
    SUB_t = T.substr(0,col_blocks*(T_SIZE-1));
  }
  else {
    for(int i=1; i<col_blocks-row_blocks+wing+1;i++) {
      temp = intToVector(prev_row_offset_vec[i]);
      row_prev_vector.insert(row_prev_vector.end(), temp.begin(), temp.end());
    }
    SUB_t = T.substr((row_blocks-wing)*(T_SIZE-1),(col_blocks-row_blocks+wing)*(T_SIZE-1));
  }

  col_last_vector.reserve(prev_col_offset_vec.size()*(T_SIZE-1));
  for(auto v : prev_col_offset_vec) {
    temp = intToVector(v);
    col_last_vector.insert(col_last_vector.end(), temp.begin(), temp.end());
  }


  if (row_remain) {
    string SUB_p = P.substr(P.length()-col_last_vector.size()-col_remain,col_last_vector.size());
    
    row_remain_vector = compute_distance(row_remain_vector,col_last_vector, sub_t, SUB_p);
    
    SUB_t += sub_t;
    row_prev_vector.insert(row_prev_vector.end(),row_remain_vector.begin(),row_remain_vector.end());
  }

  if (col_remain) {
    col_remain_vector = compute_distance(col_remain_vector,row_prev_vector,sub_p,SUB_t);
    for (auto v : col_remain_vector) accumulation += v;
  }
  
  for (auto v : row_prev_vector)  {
    accumulation += v;
  }
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
  string text_file[5];
  string pattern_file[5];
  for (int i=0; i<5; ++i) {
    text_file[i] = "text" + to_string(i) + ".txt";
    pattern_file[i] = "pattern" + to_string(i) + ".txt";
  }
  string out_filename("result.txt");

  // precomputing
  system_clock::time_point start_pre = system_clock::now();
  PossibleStringsOffsets();
  chrono::duration<double, ratio<1,1000>> time_precompute = system_clock::now() - start_pre;
  
  // result 파일 쓰기
  ofstream out_file(out_filename);
  if(!out_file.is_open()) {
    cerr << "Could not open file '" << out_filename << "'" << endl;
    return EXIT_FAILURE;
  }

  for (int filenum=0; filenum<5; ++filenum) {
    vector<string> Text = read_file(text_file[filenum]);
    vector<string> Pattern = read_file(pattern_file[filenum]);
    
    int accumulation = 0;
    string T;
	  string P;
    out_file << "==========" << filenum << "th file==========" << endl;
    // basic method
    system_clock::time_point start_basic = system_clock::now();
    out_file << "----------basic method----------" << endl;
    for(int i=0; i<Text.size(); i++) {
      accumulation = compute_basic(Text[i],Pattern[i]);
      out_file << accumulation << "\t";  
    }
    chrono::duration<double, milli> time_basic = system_clock::now() - start_basic;
    out_file << endl;
    out_file << "cost time: " << time_basic.count() << "msec" << endl;
  
    // k_difference method
    system_clock::time_point start_k = system_clock::now();
    out_file << "----------K difference method----------" << endl;
    for(int i=0; i<Text.size(); i++) {
      accumulation = compute_k_diff(Text[i],Pattern[i]);
      out_file << accumulation << "\t";  
    }
    chrono::duration<double, milli> time_k = system_clock::now() - start_k;
    out_file << endl;
    out_file << "cost time: " << time_k.count() << "msec" << endl;
    
    // four_russians method
    system_clock::time_point start_russian = system_clock::now();
    out_file << "----------Four Russians method----------" << endl;
    out_file << "precomputing time: " << time_precompute.count() << "msec" << endl;
    for(int i=0; i<Text.size(); i++) {
      accumulation = compute_russian(Text[i],Pattern[i]);
      out_file << accumulation << "\t";  
    }
    chrono::duration<double, milli> time_russian = system_clock::now() - start_russian;
    out_file << endl;
    out_file << "cost time: " << time_russian.count() << "msec" << endl;

    // four_russians and k_difference method
    system_clock::time_point start_k_russian = system_clock::now();
    out_file << "----------Four Russians and K difference method----------" << endl;
    out_file << "precomputing time: " << time_precompute.count() << "msec" << endl;
    for(int i=0; i<Text.size(); i++) {
      if(Text[i].length() < Pattern[i].length())
        accumulation = compute_k_and_russian(Pattern[i],Text[i]);
      else
        accumulation = compute_k_and_russian(Text[i],Pattern[i]);
      out_file << accumulation << "\t";  
    }
    chrono::duration<double, milli> time_k_russian = system_clock::now() - start_k_russian;
    out_file << endl;
    out_file << "cost time: " << time_k_russian.count() << "msec" << endl;
    out_file << endl;
  }
  // testcase 
  out_file.close();

  return 0;
}