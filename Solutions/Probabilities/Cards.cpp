/*
  Problem: https://codeforces.com/problemset/problem/1278/F
  Explication:
*/

#include <bits/stdc++.h>
using namespace std;

// Holi c:

#define ll long long int
#define ld long double
#define fi first
#define se second
#define pb push_back
#define all(v) v.begin(), v.end()

const int Inf = 1e9;
const ll mod = 998244353;
const ll INF = 1e18;
const int maxn = 5e3 + 5;

ll expBin(ll b, ll p){
    ll res = 1;
    while(p){
        if(p & 1) res = (res * b) % mod;
        b = (b * b) % mod;
        p >>= 1;
    }
    return res;
}

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	ll n, m, k; cin>>n>>m>>k;
	ll p = expBin(m, mod - 2);
	ll DP[maxn][maxn] = {}; DP[1][1] = 1;
	for(int i = 2; i <= 5000; i++){
	    for(int j = 1; j <= 5000; j++){
	        DP[i][j] = (DP[i - 1][j] * j + DP[i - 1][j - 1]) % mod;
	    }
	}
	ll pref[maxn] = {}; ll u = n, res = 1;
	for(int i = 1; i <= k; i++){
	    pref[i] = res = (res * u) % mod;
	    u--;
	}
	ll ans = 0;
	for(int i = 1; i <= k; i++){
	    ans += (((DP[k][i] * pref[i]) % mod) * expBin(p, i)) % mod;
	    if(ans >= mod) ans -= mod;
	}
	cout<<ans;
}
