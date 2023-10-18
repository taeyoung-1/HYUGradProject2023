#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <map>
#include <utility>
#include <fstream>
#include <chrono>
#include <cmath>

using namespace std;
using namespace chrono;
#define MAX 2147483646
#define NUM_CASE 50 // 계산하는 case의 개수
#define T_SIZE 4    // t-block의 크기
#define K_SIZE 7         // k 영역의 크기
#define VECTOR_ALL_1 pow(4,T_SIZE-1)-1
typedef pair<int,int> pInt;

int countdone;
vector<char> alphabet = {'A', 'T', 'G', 'C'};
map<pair<pInt,pInt>, pInt> precomputed_values;
// precomputed_values.reserve(10);

// 가장 기본적인 distance 계산 방식. 모든 셀 하나씩 계산
int compute_basic(const string& T, const string& P);
// k 영역 내의 모든 셀에 대해 하나씩 distance 계산.
int compute_k_diff(const string& T, const string& P);
// t-block 단위로 distance 계산
int compute_russian(const string& T, const string& P);


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

int stringToInt(const string &sequence) {
  int result = 0;
  for (int i=0; i < sequence.length(); i++) {
    result = result << 2;
    switch(sequence[i]) {
      case 'A': 
        result += 0;
        break;
      case 'C':
        result += 1;
        break;
      case 'G':
        result += 2;
        break;
      case 'T':
        result += 3;
        break;
      default: 
        break;
    }
  }

  return result;
}

int vectorToInt(vector<int> offsetVector) {
  int result;
  for (int i = 0; i < offsetVector.size(); ++i) {
    if (offsetVector[i] == 1) result += 3;
    else if (offsetVector[i] == 0) result += 2;
    else result += 1;
    result = result << 2;
  }
  return result;
}

void PrecomputeSingle(
    int row_string, // length t - 1 string
    int col_string, // length t - 1 string
    int row_offset_vector, // length t - 1 offset vector
    int col_offset_vector) { // length t - 1 offset vector
  // table initialization
  vector<vector<int>> table(T_SIZE, vector<int>(T_SIZE));
  table[0][0] = 0;
  int row_offset, col_offset;
  // cout << "first: " << row_offset_vector << "  " << col_offset_vector << endl;
  for (int i = 1; i < T_SIZE; i++) {
    // case '1': -1, case '2': 0, case '3': +1 
    row_offset = ((row_offset_vector >> 2*(T_SIZE-1 - i)) & 3) - 2;
    col_offset = ((col_offset_vector >> 2*(T_SIZE-1 - i)) & 3) - 2;
    
    table[0][i] = row_offset + table[0][i - 1];
    table[i][0] = col_offset + table[i - 1][0];
    // cout << "table[0][" << i << "]: " << table[0][i] << endl;
    // cout << "table[" << i << "][0]: " << table[i][0] << endl;
  }

  // table calculation
  for (int row = 1; row < T_SIZE; row++) {
    for (int col = 1; col < T_SIZE; col++) {
      int t = ((row_string & (3 << 2*(T_SIZE-1 - col))) == (col_string & (3 << 2*(T_SIZE-1 - row)))) ? 0 : 1;
      int diagonal = table[row - 1][col - 1] + t;
      int vertical = table[row - 1][col] + 1;
      int horizontal = table[row][col - 1] + 1;
      table[row][col] = min(diagonal, min(vertical, horizontal));
      // cout << "table[" << row << "][" << col << "]: " << table[row][col] << endl;
    }
  }

  // returning the result row and column offset vectors
  int row_output = 0;
  int col_output = 0;

  for (int i = 0; i < T_SIZE - 1; i++) {
    row_output = row_output << 2;
    col_output = col_output << 2;
    row_output += table[T_SIZE - 1][i + 1] - table[T_SIZE - 1][i] + 2;
    col_output += table[i + 1][T_SIZE - 1] - table[i][T_SIZE - 1] + 2;
  }
  countdone++;
  pInt key_str = make_pair(row_string,col_string);
  pInt key_vec = make_pair(row_offset_vector,col_offset_vector);
  precomputed_values.insert(make_pair(make_pair(key_str, key_vec), make_pair(row_output,col_output)));
}

void PossibleOffsets(int rowstr, int colstr,
                    int rowvec, int colvec,
                    int kr, int kc) {
  if (kr == 0 && kc == 0) {
    PrecomputeSingle(rowstr, colstr, rowvec, colvec);
  } else if (kr == 0) {
    for (auto i : {1, 2, 3}) {
      int temp = colvec << 2;
      temp += i;
      PossibleOffsets(rowstr, colstr, rowvec, temp, kr, kc - 1);
    }
  } else {
    for (auto i : {1, 2, 3}) {
      int temp = rowvec << 2;
      temp += i;
      PossibleOffsets(rowstr, colstr, temp, colvec, kr - 1, kc);
    }
  }
}

// T_SIZE <= 16 일 경우 정상작동
// i>j 인 경우, row와 col을 swap하면 해당 case 값을 구할 수 있으므로 계산이 1/2 정도 줄어듦
// BBBB와 ABCD의 값은 ABCD와 BBBB의 결과를 통해 알아낼 수 있음 
void PossibleStringsOffsets() {
  for (int i=0; i<(1<<2*(T_SIZE-1)); ++i) {
    for (int j=i; j<(1<<2*(T_SIZE-1)); ++j) {
      PossibleOffsets(i, j, 0, 0, T_SIZE-1, T_SIZE-1);
    }
  }
}

vector<int> intToVector(int offset_vector) {
  vector<int> result_vector;
  result_vector.reserve(T_SIZE-1);
  for (int i=0; i<T_SIZE-1; ++i) {
    result_vector.push_back(((offset_vector >> 2*(T_SIZE-2 - i)) & 3) - 2);
  }
  return result_vector;
}
pair<int, int> getPair(
    const string &row_string, // length t - 1 string
    const string &col_string, // length t - 1 string
    const int &row_offset_vector, // length t - 1 offset vector
    const int &col_offset_vector) { // length t - 1 offset vector
  pInt key_str;
  pInt key_vec;
  if (row_string.compare(col_string) > 0) {
    key_str = make_pair(stringToInt(col_string),stringToInt(row_string));
    key_vec = make_pair(col_offset_vector,row_offset_vector);
  }
  else {
    key_str = make_pair(stringToInt(row_string),stringToInt(col_string));
    key_vec = make_pair(row_offset_vector,col_offset_vector);
  }
  // cout << "key: " << row_string << "  " << col_string << "  " << row_offset_vector << "  " << col_offset_vector << endl; 
  pair<int, int> result = precomputed_values.find(make_pair(key_str,key_vec))->second;
  if (row_string.compare(col_string) > 0) 
    return make_pair(result.second,result.first);
  return make_pair(result.first,result.second);
}

vector<int> flatten(const int &vec) {
  // prev_row_offset_vec을 1차원 벡터로 변환

  vector<int> flat_vec = intToVector(vec);
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
  int m = P.length();
  int n = T.length();
  int row_blocks = m / (T_SIZE - 1);
  int col_blocks = n / (T_SIZE - 1);
  
  vector<int> prev_row_offset_vec(col_blocks, VECTOR_ALL_1); // previous row
  vector<int> last_col_offset_vec(row_blocks, VECTOR_ALL_1); // last column
  for (int row = 0; row < row_blocks; row++) {
    int col_vec = last_col_offset_vec[row];
    // col < wing+col_blocks - row 까지만 loop
    for (int col = 0; col < col_blocks; col++) {
      auto pair = getPair(T.substr(col * (T_SIZE - 1), T_SIZE - 1),
                          P.substr(row * (T_SIZE - 1), T_SIZE - 1),
                          prev_row_offset_vec[col],
                          col_vec);
      prev_row_offset_vec[col] = pair.first;
      col_vec = pair.second;
      // if (T.compare("ACTGTACG") == 0) {
      //   cout << pair.first << "\t" << pair.second << endl;
      // }

    }
    last_col_offset_vec[row] = col_vec;
  }
  accumulation += row_blocks * (T_SIZE-1);
  // cout << "acc: " << accumulation << endl;
  
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
  for (int i=0; i<col_blocks; ++i) {
    // cout << "int: " << prev_row_offset_vec[i] << endl;
    temp = intToVector(prev_row_offset_vec[i]);
    row_prev_vector.insert(row_prev_vector.end(), temp.begin(), temp.end());
    // cout << "vec: ";
    // for (auto v : temp) cout << v << " ";
    // cout << endl; 
  }
  for(int i=0; i<row_blocks; ++i) {
    temp = intToVector(last_col_offset_vec[i]);
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
  // cout << "result: " << accumulation << endl;
  return accumulation;
}

int compute_k_and_russian(const string& T, const string& P) {
  int wing = (K_SIZE - 1) / (T_SIZE - 1) + 1;
  int num_t_block_per_row = 2 * wing + 1;
  int accumulation = 0;
  int m = P.length();
  int n = T.length();
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
  prev_row_offset_vec[0] = VECTOR_ALL_1 * 2 / 3;
  
  for (int row = 0; row < row_blocks; row++) {
    temp_col_offset_vec = VECTOR_ALL_1;
    
    // col < wing+col_blocks - row 까지만 loop
    for (int col = max(0,wing-row); col < num_t_block_per_row; col++) {
      auto pair = getPair(T.substr((row + col - wing) * (T_SIZE - 1), T_SIZE - 1),
                          P.substr(row * (T_SIZE - 1), T_SIZE - 1),
                          (col == num_t_block_per_row - 1) ?
                            VECTOR_ALL_1 : prev_row_offset_vec[col + 1],
                          temp_col_offset_vec);

      prev_row_offset_vec[col] = pair.first;
      temp_col_offset_vec = pair.second;

      // row의 마지막 t-block이 text까지의 마지막 t-block 영역인 경우 column vector 저장
      if(row + col == wing + col_blocks - 1) {
        prev_col_offset_vec.push_back(pair.second);
        break;
      }
    }
    
    accumulation += T_SIZE-1;
    row_zero_vector = intToVector(prev_row_offset_vec[0]);

    for (auto v : row_zero_vector) accumulation += v;
  }
  // debug. 마지막 t-block 까지의 edit distance
  // cout << "distance: " << accumulation << endl;
  
  // 나머지가 존재하는 경우
  int col_remain = m%(T_SIZE-1);
  int row_remain = n%(T_SIZE-1);

  vector<int> col_remain_vector(col_remain, 1);
  vector<int> row_remain_vector(row_remain, 1);
  
  string sub_p = P.substr(row_blocks*(T_SIZE-1),col_remain);
  string sub_t = T.substr(col_blocks*(T_SIZE-1),row_remain);
  string SUB_t;

  // prev_row_offset_vec을 1차원 벡터로 변환
  vector<int> row_prev_vector, temp;
  row_prev_vector.reserve((col_blocks-row_blocks+wing+2)*(T_SIZE-1));
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