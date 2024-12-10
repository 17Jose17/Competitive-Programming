/*
  Problem: Rare Coins
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
const int maxn = 3e5 + 5;
const int maxm = 1e6 + 5;

struct block{
  ll total = 0;
};
 
vector<block> tree1(maxn * 4), tree2(maxn * 4);
 
void build(vector<block> & tree, vector<block> & t, int v, int tl, int tr){
    if(tl == tr){
        tree[v] = t[tl];
    }else{
        int tm = (tl + tr) / 2;
        build(tree, t, v * 2, tl, tm);
        build(tree, t, v * 2 + 1, tm + 1, tr);
        tree[v].total = tree[v * 2].total + tree[v * 2 + 1].total;
    }
}
 
block query(vector<block> & tree, int v, int tl, int tr, int l, int r){
    block y;
	if(l > r) return y;
    if(l == tl && r == tr){
        return tree[v];
    }
    int tm = (tl + tr) / 2;
    block a = query(tree, v * 2, tl, tm, l, min(r, tm)), b = query(tree, v * 2 + 1, tm + 1, tr, max(l, tm + 1), r), c;
    c.total = a.total + b.total;
    return c;
}

ll nk(ll a, ll b, vector<ll> & f, vector<ll> & i){
    if(b > a) return 0;
    return (((f[a] * i[b]) % mod) * i[a - b]) % mod;
}

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
	vector<ll> fact(maxm), inv(maxm), invFac(maxm);
	fact[0] = fact[1] = inv[0] = inv[1] = invFac[0] = invFac[1] = 1;
	for(int i = 2; i <= maxm - 3; i++){
        fact[i] = (fact[i - 1] * i) % mod;
        inv[i] = (inv[mod % i] * (mod - mod / i)) % mod;
        invFac[i] = (invFac[i - 1] * inv[i]) % mod;
	}
	int n, q; cin>>n>>q;
	vector<ll> v(n), u(n);
	for(auto & i : v) cin>>i; for(auto & i : u) cin>>i;
	ll gol = 0, plt = 0;
	for(auto e : v) gol += e; for(auto e : u) plt += e;
	vector<block> x, y; for(auto e : v){ block a; a.total = e; x.pb(a); } for(auto e : u){ block a; a.total = e; y.pb(a); }
	build(tree1, x, 1, 0, n - 1); build(tree2, y, 1, 0, n - 1);
	vector<ll> pref(maxm);
	pref[0] = 1;
	for(int i = 1; i < maxm; i++){
	    pref[i] = (pref[i - 1] + nk(plt, i, fact, invFac)) % mod;
	}
	while(q--){
	    int l, r; cin>>l>>r;
	    auto g = query(tree1, 1, 0, n - 1, l - 1, r - 1), p = query(tree2, 1, 0, n - 1, l - 1, r - 1);
	    ll a = g.total, b = p.total, c = gol - g.total, d = plt - p.total;
	    if(a + b <= c){
	        cout<<0<<" "; continue;
	    }
	    if(a > c + d){
	        cout<<1<<" "; continue;
	    }
	    ll e = a - c - 1;
	    ll ans = pref[b + e];
	    ans = (ans * expBin(expBin(2, b + d), mod - 2)) % mod;
	    cout<<ans<<" ";
	}
}
