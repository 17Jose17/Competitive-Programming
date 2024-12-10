/*
  Problem: https://codeforces.com/problemset/problem/1096/E
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
const int maxn = 1e4 + 5;

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
    if(a == b) return 1;
    if(b > a || b < 0) return 0;
    return (((f[a] * i[b]) % mod) * i[a - b]) % mod;
}

ll principleInclusionExclusion(vector<ll> & fact, vector<ll> & invFac, int n, int s, int r){
    ll res = 0;
    for(int i = 0; i <= n; i++){
        if(i & 1) res -= (1LL * nk(n, i, fact, invFac) * nk(s + n - i * r - 1, n - 1, fact, invFac)) % mod;
        else res += (1LL *nk(n, i, fact, invFac) * nk(s + n - i * r - 1, n - 1, fact, invFac)) % mod;
        if(res >= mod) res -= mod; if(res < 0) res += mod;
    }
    return res;
}

ll solve(vector<ll> & fact, vector<ll> & invFac){
    int n, r, s; cin>>n>>s>>r;
    ll ans = 0;
    for(int i = r; i <= s; i++){
        for(int j = 1; j <= n; j++){
            ans += (((nk(n - 1, j - 1, fact, invFac) * expBin(j, mod - 2)) % mod) *
             principleInclusionExclusion(fact, invFac, n - j, s - i * j, i)) % mod;
            if(ans >= mod) ans -= mod;
        }
    }
    return (ans * expBin(nk(s + n - r - 1, n - 1, fact, invFac), mod - 2)) % mod;
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
	cout<<solve(fact, invFac);
} 
