/*
  Problem: https://codeforces.com/problemset/problem/1418/E
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
const int maxn = 5e5 + 5;

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
	int n, q; cin>>n>>q; vector<ll> v(n);
	for(auto & i : v) cin>>i;
	vector<ll> pref(n); sort(all(v)); pref[0] = v[0];
	for(int i = 1; i < n; i++) pref[i] = pref[i - 1] + v[i]; 
	while(q--){
	    ll a, b, ans = 0; cin>>a>>b;
	    auto at = lower_bound(all(v), b) - v.begin();
	    if(a >= n || a + at - 1 >= n){ cout<<0<<'\n'; continue; }
	    ll e = pref[n - 1] - pref[at] + v[at], f = pref[at] - v[at];
	    if(at + a >= n){
	        e = 0;
	    }
	    if(v[0] > b){
	        f = 0; e = pref[n - 1];
	    }
	    ll res = (((nk(n - at - 1, a, fact, invFac) * fact[n - a - at]) % mod) * fact[a]) % mod;
	    res = (((res * fact[n]) % mod) * invFac[n - at]) % mod;
	    e %= mod; f %= mod;
	    ans += (res * e) % mod;
	    ans %= mod;
	    res = (((nk(n - at, a, fact, invFac) * fact[a]) % mod) * fact[n - a - at + 1]) % mod;
	    res = (((res * fact[n]) % mod) * invFac[n - at + 1]) % mod;
	    ans += (res * f) % mod;
	    ans %= mod;
	    cout<<(ans * invFac[n]) % mod<<'\n';
	}
}
