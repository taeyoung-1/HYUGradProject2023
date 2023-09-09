#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
using namespace std;

int countDone = 0;

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
  cout << endl;
}

pair<vector<int>, vector<int>> PrecomputeSingle(
    const string &row_string, // length t - 1 string
    const string &column_string, // length t - 1 string
    const vector<int> &row_offset_vector, // length t offset vector
    const vector<int> &column_offset_vector, // length t offset vector
    int t) {
  // table initialization
  vector<vector<int>> table(t, vector<int>(t));
  table[0][0] = 0;
  for (int i = 1; i < t; i++) {
    table[0][i] = row_offset_vector[i] + table[0][i - 1];
    table[i][0] = column_offset_vector[i] + table[i - 1][0];
  }

  // table calculation
  for (int row = 1; row < t; row++) {
    for (int col = 1; col < t; col++) {
      int t = (row_string[row - 1] == column_string[col - 1]) ? 0 : 1;
      int diagonal = table[row - 1][col - 1] + t;
      int vertical = table[row - 1][col] + 1;
      int horizontal = table[row][col - 1] + 1;
      table[row][col] = min(diagonal, min(vertical, horizontal));
    }
  }

  // returning the result row and column offset vectors
  vector<int> row_offset_vector_output(t);
  vector<int> column_offset_vector_output(t);

  for (int i = 1; i < t; i++) {
    row_offset_vector_output[i] = table[t - 1][i] - table[t - 1][i - 1];
    column_offset_vector_output[i] = table[i][t - 1] - table[i - 1][t - 1];
  }

  countDone++;
  return make_pair(row_offset_vector_output, column_offset_vector_output);
}

void PossibleOffsets(const string &rowstr, const string &colstr,
                    vector<int> row, vector<int> col,
                    int kr, int kc, int t) {
  if (kr == 0 && kc == 0) {
    PrecomputeSingle(rowstr, colstr, row, col, t);
  } else if (kr == 0) {
    for (auto i : {-1, 0, 1}) {
      col.push_back(i);
      PossibleOffsets(rowstr, colstr, row, col, kr, kc - 1, t);
    }
  } else {
    for (auto i : {-1, 0, 1}) {
      row.push_back(i);
      PossibleOffsets(rowstr, colstr, row, col, kr - 1, kc, t);
    }
  }
}

void PossibleStringsOffsets(const vector<char> &alphabets,
                            string rowstr, string colstr,
                            int kr, int kc, int t) {
  if (kr == 0 && kc == 0) {
    PossibleOffsets(rowstr, colstr, {}, {}, t, t, t);
  } else if (kr == 0) {
    for (auto c : alphabets) {
      colstr.push_back(c);
      PossibleStringsOffsets(alphabets, rowstr, colstr, kr, kc - 1, t);
    }
  } else {
    for (auto c : alphabets) {
      rowstr.push_back(c);
      PossibleStringsOffsets(alphabets, rowstr, colstr, kr - 1, kc, t);
    }
  }
}

int main(void) {
  int alphabet_size;
  vector<char> alphabets;
  cout << "alphabet set size: ";
  cin >> alphabet_size;

  cout << "alphabets: ";
  for (int i = 0; i < alphabet_size; ++i) {
    char c;
    cin >> c;
    alphabets.push_back(c);
  }

  int t_size;
  cout << "t-block size: ";
  cin >> t_size;

  PossibleStringsOffsets(alphabets, "", "", t_size - 1, t_size -1, t_size);
  cout << countDone;
  return 0;
}