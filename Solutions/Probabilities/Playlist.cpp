/*
  Problem: https://codeforces.com/problemset/problem/268/E
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
const ll mod = 998244353;
const ll INF = 1e18;
const int maxn = 3e2 + 5;
const ld eps = 1e-12;

bool cmp(pair<ld, ld> a, pair<ld, ld> b){
    ll a1 = a.fi, a2 = b.fi, b1 = a.se, b2 = b.se;
    return (a1 * 2 + a2) * (b1 * (100 - b2)) + (a1 + a2) * ((100 - b1) * b2) >
     (a2 * 2 + a1) * (b2 * (100 - b1)) + (a1 + a2) * ((100 - b2) * b1);
}

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	int n; cin>>n; 
	vector<pair<ld, ld>> v(n);
	for(int i = 0; i < n; i++){
	    int a, b; cin>>a>>b;
	    v[i] = {a, b};
	}
	sort(all(v), cmp);
	ld ans = 0, res = 0;
	for(int i = 0; i < n; i++){
	    ans += v[i].fi;
	    ans += res * ((ld)1 - (v[i].se / 100));
	    res += v[i].fi * v[i].se / 100;
	}
	cout<<setprecision(20)<<ans;
}
