/*
  Problem: https://codeforces.com/gym/103640/problem/A
  Explication (Spanish): https://hackmd.io/@JoOrI6xWTeCaxNiU83oinw/SJK6Of6rp
*/

#include <bits/stdc++.h>
using namespace std;

// Holi c:

#define ll long long int
#define ld long long int
#define fi first
#define se second
#define pb push_back
#define all(v) v.begin(), v.end()

const ld eps = 0, inf = numeric_limits<ld>::max();
bool geq(ld a, ld b){return a-b >= -eps;}     //a >= b
bool leq(ld a, ld b){return b-a >= -eps;}     //a <= b
bool ge(ld a, ld b){return a-b > eps;}        //a > b
bool le(ld a, ld b){return b-a > eps;}        //a < b
bool eq(ld a, ld b){return abs(a-b) <= eps;}  //a == b
bool neq(ld a, ld b){return abs(a-b) > eps;}  //a != b

int BIT[401] = {};

void update(int pos, int value){
    while(pos < 802){
        BIT[pos] += value;
        pos |= pos + 1;
    }
}

//range query, [0, r]
int queryPrefix(int r){
    int res = 0;
    while(r >= 0){
        res += BIT[r];
        r = (r & (r + 1)) - 1;
    }
    return res;
}

//range query, [l, r]
int query(int l, int r){
    return queryPrefix(r) - queryPrefix(l - 1);
}


struct point{
	ld x, y;
	point(): x(0), y(0){}
	point(ld x, ld y): x(x), y(y){}

	point operator+(const point & p) const{return point(x + p.x, y + p.y);}
	point operator-(const point & p) const{return point(x - p.x, y - p.y);}
	point operator*(const ld & k) const{return point(x * k, y * k);}
	point operator/(const ld & k) const{return point(x / k, y / k);}

	point operator+=(const point & p){*this = *this + p; return *this;}
	point operator-=(const point & p){*this = *this - p; return *this;}
	point operator*=(const ld & p){*this = *this * p; return *this;}
	point operator/=(const ld & p){*this = *this / p; return *this;}

	point rotate(const ld & a) const{return point(x*cos(a) - y*sin(a), x*sin(a) + y*cos(a));}
	point perp() const{return point(-y, x);}
	ld dot(const point & p) const{return x * p.x + y * p.y;}
	ld cross(const point & p) const{return x * p.y - y * p.x;}
	ld norm() const{return x * x + y * y;}
	ld length() const{return sqrt(x * x + y * y);}
	point unit() const{return (*this) / length();}

	bool operator==(const point & p) const{return eq(x, p.x) && eq(y, p.y);}
	bool operator!=(const point & p) const{return !(*this == p);}
	bool operator<(const point & p) const{return le(x, p.x) || (eq(x, p.x) && le(y, p.y));}
	bool operator>(const point & p) const{return ge(x, p.x) || (eq(x, p.x) && ge(y, p.y));}
	bool half(const point & p) const{return le(p.cross(*this), 0) || (eq(p.cross(*this), 0) && le(p.dot(*this), 0));}
};

istream &operator>>(istream &is, point & p){return is >> p.x >> p.y;}
ostream &operator<<(ostream &os, const point & p){return os << "(" << p.x << ", " << p.y << ")";}

void polarSort(vector<pair<point, ll>> & P, const point & o, const point & v){
	sort(P.begin(), P.end(), [&](const pair<point,ll> & a, const pair<point, ll> & b){
		return point((a.fi - o).half(v), 0) < point((b.fi - o).half(v), (a.fi - o).cross(b.fi - o));
	});
}

ld Area_t(point p0, point p1, point p2){
	return abs((p1.x - p0.x) * (p2.y - p1.y) - (p2.x - p1.x) * (p1.y - p0.y));
}

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	ll s, n; cin>>s>>n; 
	s = s * 2;
    point vp[500] = {};
	for(int i = 0; i < n; i++)
		cin>>vp[i];
	ll res = 0, res_n = 0, res_n1 = 0;
	for(int i = 0; i < n; i++){
	    for(int j = i + 1; j < n; j++){
	        vector<pair<point, ll>> vps, vps1;
	        ll vt[400] = {}; ll vts[400] = {};
	        int tam_vt = 0, tam_vts = 0;
	        for(int l = 0; l < n; l++){
	            ll p = (vp[j] - vp[i]).cross(vp[l] - vp[i]);
	            if(p < 0){
	                ll u = Area_t(vp[i], vp[j], vp[l]);
	                vt[tam_vt] = u;
	                vps.pb({vp[l], u});
	                tam_vt++;
	            }else if(p > 0){
	                ll u = Area_t(vp[i], vp[j], vp[l]);
	                vts[tam_vts] = u;
	                vps1.pb({vp[l], u});
	                tam_vts++;    }   }
	                
	        if(!tam_vt || !tam_vts) continue;
	        sort(vt, vt + tam_vt); sort(vts, vts + tam_vts);
	        if(vt[tam_vt - 1] + vts[tam_vts - 1] < s) continue;
	        
	        int t = tam_vts - 1;
            for(int l1 = 0; l1 < tam_vt; l1++){
                if(t == tam_vts) t--;
                while(vt[l1] + vts[t] >= s && t > -1)
	                t--;    
	            t++;
	            res += tam_vts - t;    }
	        
	        vector<ll> areaTs;
	        for(int l1 = 0; l1 < tam_vts; l1++){
	            areaTs.pb(vts[l1]);
	            }
	   
            for(int i = 0; i <= 400; i++) BIT[i] = 0;

	        polarSort(vps, vp[i], vp[j] - vp[i]); polarSort(vps1, vp[i], vp[j] - vp[i]);
	            
	        int l = vps1.size()- 1;
	        ld maxx = 0;

            for(int l1 = vps.size()-1; l1 > -1; l1--){
                
                if(s > vps[l1].se + vts[tam_vts-1]) continue;
                
                auto it = lower_bound(areaTs.begin(), areaTs.end(), s - vps[l1].se);
                int lx = (it - areaTs.begin());
                
                while(((vp[i] - vps[l1].fi).cross(vps1[l].fi - vps[l1].fi) > 0) && l > -1){
                    auto at = lower_bound(areaTs.begin(), areaTs.end(), vps1[l].se);
                    auto at1 =  at - areaTs.begin();
                    update(at1, 1);
                    maxx = max(maxx, vps1[l].se);
                	l--;
                }
                	
				if(maxx + vps[l1].se < s) continue;
                res_n1 += query(lx, n + 1);
            }
		    
            for(int i = 0; i <= 400; i++) BIT[i] = 0;
            
            polarSort(vps, vp[j], vp[i] - vp[j]); polarSort(vps1, vp[j], vp[i] - vp[j]);
            
            l = maxx = 0;
            for(int l1 = 0; l1 < vps.size(); l1++){
                
                if(s > vps[l1].se + vts[tam_vts-1]) continue;
                
                auto it = lower_bound(areaTs.begin(), areaTs.end(), s - vps[l1].se);
                auto lx = (it - areaTs.begin());
                
                while(((vp[j] - vps[l1].fi).cross(vps1[l].fi - vps[l1].fi) < 0) && l < vps1.size()){
                    auto at = lower_bound(areaTs.begin(), areaTs.end(), vps1[l].se);
                    auto at1 = at - areaTs.begin();
                    update(at1, 1);
                    maxx = max(maxx, vps1[l].se);
                	l++;
                }
                	
				if(maxx + vps[l1].se < s) continue;
                res_n1 += query(lx, n + 1);
            }
	    }
	}
    cout<<(res + res_n1) / 2;
}
