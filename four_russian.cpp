#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
using namespace std;

// precomputing


// padding
void padding(int t, string& P, string& T) {
    int n = P.length();
    int m = T.length();

    int front_pad, back_pad;
    int rem_n = n % (t-1);
    int rem_m = m % (t-1);

    char alphabet_pad = P[0];   // pattern 의 첫번째 문자로 padding

    front_pad = min(rem_n,rem_m);
    back_pad = abs(rem_n - rem_m);

    P.insert(0, front_pad, alphabet_pad);
    P.insert(P.length(), back_pad, alphabet_pad);
    T.insert(0, front_pad, alphabet_pad);
    T.insert(T.length(), back_pad, alphabet_pad);
}

// (1,1)에서 P[0], T[0] 문자가 대응, (1,0)에서 P[0], space 대응
int k_diff_and_four_russian(int k,int t, string P, string T) {
    int n = P.length(); // abcdefgh 0~7 length = 8
    int m = T.length();

    // (i,i) 와 (i,i-k-(t-1)) 사이의 t-block 개수
    int k_t = (k-1)/(t-1)+1;
    
    vector<int> row(t-1,+1); // t-block의 first row
    vector<vector<int> > list_row(2*k_t + 1, vector<int>(t-1,+1));

    for(int i=0; i<=n-(t-1); i+=t-1) { // padding 시, i<n 여도 계산 안하는 경우 존재(n==i+a)
        string sub_p = P.substr(i,t);

        vector<int> col(t-1,+1); // t-block의 first column
        for(int j=max(0,i-k_t*(t-1)); j<i+k+(t-1); j+=t-1) {
            
            string sub_t = T.substr(j,t);
            // row 값을 list의 어디에서 가져와야 하는지..
                // 처음 k_t번(for i 루프 횟수)은 list[0] 통해 list[0] 채움
                // 그 이후로는 첫번째 list 값은 버리고 list[1] 통해 list[0] 채움
                // list[-1] 에는 항상 +1 만으로 구성되어 있어야 함

            // offset 값도 list 벡터를 통해..?
                // 처음 k_t번(for i 루프 횟수)은 list[0] 위치에 (t-1)*count 값 들어감
                // 그 이후로는 list[0] 값 + precomputing 값 더한 결과를 list[1]에 채움
            // vector<int> list_offset(2*k_t + 1, 0);
            // for(int l=0; l < 2*k_t + 1; i++)
                // list_offset[l] = l * (t-1);
            
            // (n,m) 위치 값을 list_offset과 col, row 값을 통해 구함..
                // front_padding 값을 빼줘야 함.
                // 어떻게 값 구하는지..

            // precomputing 결과 가져오기
            // col 값은 업데이트하고, row 값은 list 벡터에 저장
            // get_offset(col,row,list_row,sub_p,sub_t);
            
        }
    }
}

int main(int argc, char** argv) {
    int k = atoi(argv[1]);
    int t = atoi(argv[2]);
    string pattern = argv[3];
    string text = argv[4];

    padding(t,pattern,text);

    k_diff_and_four_russian(k,t,pattern, text);


    for(int i=0; i < argc; i++)
        cout << "argv[" << i << "] : " << argv[i] << endl;
    cout << "pattern: " << pattern << endl;
    cout << "text: " << text << endl;
    
    return 0;
}