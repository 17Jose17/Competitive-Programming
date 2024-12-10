/*
  Problem: https://codeforces.com/problemset/problem/633/E
  Explication:
*/

#include <bits/stdc++.h>
using namespace std;

// Holi c:

#define ll long long int
#define ld long double
#define ii __int128
#define fi first
#define se second
#define pb push_back
#define all(v) v.begin(), v.end()

const int Inf = 1e9;
const ll mod =  998244353;
const ll INF = 1e18;
const int maxn = 1e6 + 8;

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	int n, k; cin>>n>>k;
	vector<int> v1(n), u1(n), v(n);
	for(int i = 0; i < n; i++){
	    int a; cin>>a; v1[i] = a * 100;
	}
	for(int i = 0; i < n; i++){
	    int a; cin>>a; u1[i] = a;
	}
	//Precalcular todos los p(li)
	multiset<int> side1, side2;
	int j = n - 1;
	for(int i = n - 1; i >= 0; i--){
	    side1.insert(v1[i]); side2.insert(u1[i]);
	    auto it = side1.end(); it--; int c = *it, d = *side2.begin();
	    while(j > i && min(v1[i], u1[i]) > min(c, d)){
	        side1.erase(side1.find(v1[j])); side2.erase(side2.find(u1[j])); j--;
	        it = side1.end(); it--; c = *it; d = *side2.begin();
	    }
	    auto at = side1.end(); at--;
	    v[i] = min(*at, *side2.begin());
	}
	sort(all(v));
	int a = n, b = n - k + 1; ld res = (ld) k / n;
	ld ans = 0;
	for(int i = 0; i <= n - k; i++){
	    ans += res * v[i]; a--; b--;
	    res = res * b / a;
	}
	cout<<setprecision(20)<<ans;
}
