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

// Assumption: n = n'(k - 1), m = m'(k - 1) for positives n' m'.
void overAllTest(int k, const string& T, const string& P) {
  int wing = (k - 1) / (t_size - 1) + 1;
  cout << "wing: " << wing << endl;
  int num_t_block_per_row = 2 * wing + 1;
  cout << "number of t blocks per row: " << num_t_block_per_row << endl;

  int m = P.length();
  int n = T.length();
  int row_blocks = m / (t_size - 1);
  int col_blocks = n / (t_size - 1);
  cout << "row_blocks: " << row_blocks << ", col_blocks: " << col_blocks << endl;

  vector<int> temp_col_offset_vec(t_size - 1); // temporary column vector
  vector<vector<int>> prev_row_offset_vec(num_t_block_per_row); // previous row
  int accumulation = 0;
  fill(prev_row_offset_vec.begin(), prev_row_offset_vec.end(), vector<int>(t_size - 1, 1));

  for (int row = 0; row < row_blocks; row++) {
    fill(temp_col_offset_vec.begin(), temp_col_offset_vec.end(), 1);
    for (int col = 0; col < num_t_block_per_row; col++) {
      if (row + col < wing || row + col >= wing + col_blocks) continue;
      auto pair = getPair(T.substr((row + col - wing) * (t_size - 1), t_size - 1),
                          P.substr(row * (t_size - 1), t_size - 1),
                          (col == num_t_block_per_row - 1) ?
                            vector<int>(t_size - 1, 1) : prev_row_offset_vec[col + 1],
                          temp_col_offset_vec);

      // cout << "row: " << row << ", col: " << col
          //  << ", Tsub: " << T.substr((row + col - wing) * (t_size - 1), t_size - 1)
          //  << ", Psub: " << P.substr(row * (t_size - 1), t_size - 1)
          //  << ", row: ";
      // debugVector(get<0>(pair));
      // cout << ", col: ";
      // debugVector(get<1>(pair));
      // cout << endl;

      if (col == wing) {
        for (auto v : temp_col_offset_vec)
          accumulation += v;
        for (auto v : get<0>(pair))
          accumulation += v;
      }
      prev_row_offset_vec[col] = get<0>(pair);
      temp_col_offset_vec = get<1>(pair);
    }
  }

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

  PossibleStringsOffsets("", "", t_size - 1, t_size - 1);
  overAllTest(3, "CCACAGGATTGAATCGAACCATAGTAAATAGGGTGACCAACAGAGTAGACG", "ACAGGATAGGGGATCGAGATAGGAAAACCGATGGAGCGAAGGGTGCGAGCG");

  // cout << countdone << endl;
  // cout << precomputed_values.size() << endl;
  // auto pair = (precomputed_values.find("GGGGGG111111")->second);
  return 0;
}