/*
  Problem: https://codeforces.com/problemset/problem/26/D
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

ld calc(int a, int b, int c, int d){
    vector<ld> v; ld t = d;
    for(ld i = c + d; i > c; i--){
        v.pb(i / max((ld)1.0, t)); t--;
    }
    vector<ld> u; t = b;
    for(ld i = a + b; i > a; i--){
        u.pb(i / max((ld)1.0, t)); t--;
    }
    a = c - a - 1; t = d - b - 1; b = d - b - 1;
    ld res = 1;
    for(int i = 0; i < min((int)v.size(), (int)u.size()); i++) res *= u[i] / v[i];
    if(v.size() > u.size()){
        for(int i = u.size(); i < v.size(); i++) res /= v[i];
    }else if(u.size() > v.size()){
        for(int i = v.size(); i < u.size(); i++) res *= u[i];
    }
    return res;
}

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	int n, m, k; cin>>n>>m>>k;
	if(m <= k){ cout<<1; return 0; }
	if(m > n + k){ cout<<0; return 0; }
	int a = 0; ld ans = 0;
	ans = calc(n + k + 1, m - k - 1, n, m);
	cout<<setprecision(20)<<1.0 - ans;
}
