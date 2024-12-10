/*
  Problem: https://codeforces.com/problemset/problem/1194/F
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
const ll mod = 1e9 + 7;
const ll INF = 1e18;
const int maxn = 2e5 + 5;

ll expBin(ll b, ll p){
    ll res = 1;
    while(p){
        if(p & 1) res = (res * b) % mod;
        b = (b * b) % mod;
        p >>= 1;
    }
    return res;
}

ll nk(ll a, ll b, vector<ll> & f, vector<ll> & i){
    if(b > a || b < 1) return 1;
    return (((f[a] * i[b]) % mod) * i[a - b]) % mod;
}

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	vector<ll> fact(maxn), inv(maxn), invFac(maxn);
	fact[0] = fact[1] = inv[0] = inv[1] = invFac[0] = invFac[1] = 1;
	for(int i = 2; i <= maxn - 3; i++){
        fact[i] = (fact[i - 1] * i) % mod;
        inv[i] = (inv[mod % i] * (mod - mod / i)) % mod;
        invFac[i] = (invFac[i - 1] * inv[i]) % mod;
	}
	int n; ll k; cin>>n>>k; vector<ll> v(n);
	for(auto & i : v) cin>>i;
	ll act = 0; int j = 0;
	for(int i = 0; i < n; i++){
	    if(act + v[i] > k) break;
	    act += v[i]; j = i;
	}
	if(j == n - 1 && act + n <= k){
	    cout<<n; return 0;
	}
	ll ans = 0; int l = 0;
	for(ll i = j; i >= 0; i--){
	    while(act + l <= k && l <= i + 1){
	        ll r = ((((i + 1) * expBin(2, n - i - 1)) % mod) * nk(i + 1, l, fact, invFac)) % mod;
	        ans += r; if(ans >= mod) ans -= mod;
	        l++;
	    }
	    
	    if(i != j && act + v[i + 1] + i + 2 > k){
	        ans += ((((i + 1) * expBin(2, n - i - 2)) % mod) * nk(i + 1, k - v[i + 1] - 1 - act + 1, fact, invFac)) % mod;
	        if(ans >= mod) ans -= mod;
	    }
	    act -= v[i];
	}
	cout<<(ans * expBin(expBin(2, n), mod - 2)) % mod;
}
