/*
  Problem: https://codeforces.com/contest/442/problem/E
  Explication:
*/

#include <bits/stdc++.h>

#define endl '\n'
#define fi first
#define se second
#define MOD(n,k) ( ( ((n) % abs(k)) + abs(k) ) % abs(k))
#define forn(i,n) for (int i = 0; i < int(n); i++)
#define forr(i,a,b) for (int i = int(a); i <= int(b); i++)
#define all(v) v.begin(), v.end()
#define pb push_back

using namespace std;

typedef long long ll;
typedef long double ld;
typedef pair<int, int> ii;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<ii> vii;

ld eps = 1e-9, inf = 1e9, inf_coord = 1e9;

bool geq(ld a, ld b){return a-b >= -eps;}     //a >= b
bool leq(ld a, ld b){return b-a >= -eps;}     //a <= b
bool ge(ld a, ld b){return a-b > eps;}        //a > b
bool le(ld a, ld b){return b-a > eps;}        //a < b
bool eq(ld a, ld b){return abs(a-b) <= eps;}  //a == b
bool neq(ld a, ld b){return abs(a-b) > eps;}  //a != b

struct point{
	ld x, y;
	int idx = -1;
	point() {}
	point(ld x, ld y): x(x), y(y) {}
	
	point operator+(const point & p) const{return point(x + p.x, y + p.y);}

	point operator-(const point & p) const{return point(x - p.x, y - p.y);}

	point operator*(ld k) const{return point(x * k, y * k);}

	point operator/(ld k) const{return point(x / k, y / k);}

	ld dot(const point & p) const{return x * p.x + y * p.y;}

	ld cross(const point & p) const{return x * p.y - y * p.x;}

	ld norm() const{return x * x + y * y;}

	ld length() const{return sqrtl(x*x + y*y);}
	
	point unit() const{return (*this) / length();}

	point perpendicular() const{return point(-y, x);}

	bool operator<(const point & p) const{
		if(eq(x, p.x)) return le(y, p.y);
		return le(x, p.x);
	}

	bool operator>(const point & p) const{
		if(eq(x, p.x)) return ge(y, p.y);
		return ge(x, p.x);
	}

	bool operator==(const point & p) const{return eq(x, p.x) && eq(y, p.y);}

	bool operator!=(const point & p) const{return !(*this == p);}
};

istream&operator>>(istream & s, point & p){
	return s >> p.x >> p.y;
}

ostream &operator<<(ostream &os, const point & p){
	return os << "(" << p.x << ", " << p.y << ")";
}

int sgn(ld x){
	if(ge(x, 0)) return 1;
	if(le(x, 0)) return -1;
	return 0;
}

vector<point> convexHull(vector<point> P){
	sort(P.begin(), P.end());
	P.erase(unique(P.begin(), P.end()), P.end());
	if(P.size() <= 2) return P;
	vector<point> L, U;
	for(int i = 0; i < P.size(); i++){
		while(L.size() >= 2 && leq((L[L.size() - 2] - P[i]).cross(L[L.size() - 1] - P[i]), 0)){
			L.pop_back();
		}
		L.push_back(P[i]);
	}
	for(int i = P.size() - 1; i >= 0; i--){
		while(U.size() >= 2 && leq((U[U.size() - 2] - P[i]).cross(U[U.size() - 1] - P[i]), 0)){
			U.pop_back();
		}
		U.push_back(P[i]);
	}
	L.pop_back();
	U.pop_back();
	L.insert(L.end(), U.begin(), U.end());
	return L;
}

vector<point> convexHullWithColinear(vector<point> P){
	sort(P.begin(), P.end());
	P.erase(unique(P.begin(), P.end()), P.end());
	if(P.size() <= 2) return P;
	vector<point> L, U;
	for(int i = 0; i < P.size(); i++){
		while(L.size() >= 2 && le((L[L.size() - 2] - P[i]).cross(L[L.size() - 1] - P[i]), 0)){
			L.pop_back();
		}
		L.push_back(P[i]);
	}
	for(int i = P.size() - 1; i >= 0; i--){
		while(U.size() >= 2 && le((U[U.size() - 2] - P[i]).cross(U[U.size() - 1] - P[i]), 0)){
			U.pop_back();
		}
		U.push_back(P[i]);
	}
	L.pop_back();
	U.pop_back();
	L.insert(L.end(), U.begin(), U.end());
	return L;
}

point intersectLines(const point & a1, const point & v1, const point & a2, const point & v2){
	ld det = v1.cross(v2);
	return a1 + v1 * ((a2 - a1).cross(v2) / det);
}

int intersectLineSegmentInfo(const point & a, const point & v, const point & c, const point & d){
	point v2 = d - c;
	ld det = v.cross(v2);
	if(eq(det, 0)){
		if(eq((c - a).cross(v), 0)){
			return -1; //infinity points
		}else{
			return 0; //no point
		}
	}else{
		return sgn(v.cross(c - a)) != sgn(v.cross(d - a)); //1: single point, 0: no point
	}
}

pair<point, ld> getCircle(const point & m, const point & n, const point & p){
	point c = intersectLines((n + m) / 2, (n - m).perpendicular(), (p + n) / 2, (p - n).perpendicular());
	ld r = (c - m).length();
	return {c, r};
}

const point inf_pt(inf, inf);

struct QuadEdge{
	point origin;
	QuadEdge* rot = nullptr;
	QuadEdge* onext = nullptr;
	bool used = false;
	QuadEdge* rev() const{return rot->rot;}
	QuadEdge* lnext() const{return rot->rev()->onext->rot;}
	QuadEdge* oprev() const{return rot->onext->rot;}
	point dest() const{return rev()->origin;}
};

QuadEdge* make_edge(const point & from, const point & to){
	QuadEdge* e1 = new QuadEdge;
	QuadEdge* e2 = new QuadEdge;
	QuadEdge* e3 = new QuadEdge;
	QuadEdge* e4 = new QuadEdge;
	e1->origin = from;
	e2->origin = to;
	e3->origin = e4->origin = inf_pt;
	e1->rot = e3;
	e2->rot = e4;
	e3->rot = e2;
	e4->rot = e1;
	e1->onext = e1;
	e2->onext = e2;
	e3->onext = e4;
	e4->onext = e3;
	return e1;
}

void splice(QuadEdge* a, QuadEdge* b){
	swap(a->onext->rot->onext, b->onext->rot->onext);
	swap(a->onext, b->onext);
}

void delete_edge(QuadEdge* e){
	splice(e, e->oprev());
	splice(e->rev(), e->rev()->oprev());
	delete e->rev()->rot;
	delete e->rev();
	delete e->rot;
	delete e;
}

QuadEdge* connect(QuadEdge* a, QuadEdge* b){
	QuadEdge* e = make_edge(a->dest(), b->origin);
	splice(e, a->lnext());
	splice(e->rev(), b);
	return e;
}

bool left_of(const point & p, QuadEdge* e){
	return ge((e->origin - p).cross(e->dest() - p), 0);
}

bool right_of(const point & p, QuadEdge* e){
	return le((e->origin - p).cross(e->dest() - p), 0);
}

__int128_t det3(__int128_t a1, __int128_t a2, __int128_t a3, __int128_t b1, __int128_t b2, __int128_t b3, __int128_t c1, __int128_t c2, __int128_t c3) {
	return a1 * (b2 * c3 - c2 * b3) - a2 * (b1 * c3 - c1 * b3) + a3 * (b1 * c2 - c1 * b2);
}

bool in_circle(const point & a, const point & b, const point & c, const point & d) {
	__int128_t det = -det3(b.x, b.y, b.norm(), c.x, c.y, c.norm(), d.x, d.y, d.norm());
	det += det3(a.x, a.y, a.norm(), c.x, c.y, c.norm(), d.x, d.y, d.norm());
	det -= det3(a.x, a.y, a.norm(), b.x, b.y, b.norm(), d.x, d.y, d.norm());
	det += det3(a.x, a.y, a.norm(), b.x, b.y, b.norm(), c.x, c.y, c.norm());
	return det > 0;
}

pair<QuadEdge*, QuadEdge*> build_tr(int l, int r, vector<point> & P){
	if(r - l + 1 == 2){
		QuadEdge* res = make_edge(P[l], P[r]);
		return make_pair(res, res->rev());
	}
	if(r - l + 1 == 3){
		QuadEdge *a = make_edge(P[l], P[l + 1]), *b = make_edge(P[l + 1], P[r]);
		splice(a->rev(), b);
		int sg = sgn((P[l + 1] - P[l]).cross(P[r] - P[l]));
		if(sg == 0)
			return make_pair(a, b->rev());
		QuadEdge* c = connect(b, a);
		if(sg == 1)
			return make_pair(a, b->rev());
		else
			return make_pair(c->rev(), c);
	}
	int mid = (l + r) / 2;
	QuadEdge *ldo, *ldi, *rdo, *rdi;
	tie(ldo, ldi) = build_tr(l, mid, P);
	tie(rdi, rdo) = build_tr(mid + 1, r, P);
	while(true){
		if(left_of(rdi->origin, ldi)){
			ldi = ldi->lnext();
			continue;
		}
		if(right_of(ldi->origin, rdi)){
			rdi = rdi->rev()->onext;
			continue;
		}
		break;
	}
	QuadEdge* basel = connect(rdi->rev(), ldi);
	auto valid = [&basel](QuadEdge* e){return right_of(e->dest(), basel);};
	if(ldi->origin == ldo->origin)
		ldo = basel->rev();
	if(rdi->origin == rdo->origin)
		rdo = basel;
	while(true){
		QuadEdge* lcand = basel->rev()->onext;
		if(valid(lcand)){
			while(in_circle(basel->dest(), basel->origin, lcand->dest(), lcand->onext->dest())){
				QuadEdge* t = lcand->onext;
				delete_edge(lcand);
				lcand = t;
			}
		}
		QuadEdge* rcand = basel->oprev();
		if(valid(rcand)){
			while(in_circle(basel->dest(), basel->origin, rcand->dest(), rcand->oprev()->dest())){
				QuadEdge* t = rcand->oprev();
				delete_edge(rcand);
				rcand = t;
			}
		}
		if(!valid(lcand) && !valid(rcand))
			break;
		if(!valid(lcand) || (valid(rcand) && in_circle(lcand->dest(), lcand->origin, rcand->origin, rcand->dest())))
			basel = connect(rcand, basel->rev());
		else
			basel = connect(basel->rev(), lcand->rev());
	}
	return make_pair(ldo, rdo);
}

vector<vector<point>> delaunay(vector<point> P){
	sort(P.begin(), P.end());
	auto res = build_tr(0, (int)P.size() - 1, P);
	QuadEdge* e = res.first;
	vector<QuadEdge*> edges = {e};
	while(le((e->dest() - e->onext->dest()).cross(e->origin - e->onext->dest()), 0))
		e = e->onext;
	auto add = [&P, &e, &edges](){
		QuadEdge* curr = e;
		do{
			curr->used = true;
			P.push_back(curr->origin);
			edges.push_back(curr->rev());
			curr = curr->lnext();
		}while(curr != e);
	};
	add();
	P.clear();
	int kek = 0;
	while(kek < (int)edges.size())
		if(!(e = edges[kek++])->used)
			add();
	vector<vector<point>> ans;
	for(int i = 0; i < (int)P.size(); i += 3){
		ans.push_back({P[i], P[i + 1], P[i + 2]});
	}
	return ans;
}

pair<point, point> divPlane(point u, point v){
    point res, res1; point p, q;
    point aux = (u - v).perpendicular();
    p = (u + v) / 2;
    q = p + aux.unit();
    res = p; res1 = q;
    if((p - u).cross(q - u) < 0) swap(res, res1);
    return {res, res1};
}

vector<point> cutPolygon(const vector<point> & P, const point & a, const point & v){
	//returns the part of the convex polygon P on the left side of line a+tv
	int n = P.size();
	vector<point> lhs;
	for(int i = 0; i < n; ++i){
		if(geq(v.cross(P[i] - a), 0)){
			lhs.push_back(P[i]);
		}
		if(intersectLineSegmentInfo(a, v, P[i], P[(i+1)%n]) == 1){
			point p = intersectLines(a, v, P[i], P[(i+1)%n] - P[i]);
			if(p != P[i] && p != P[(i+1)%n]){
				lhs.push_back(p);
			}
		}
	}
	return lhs;
}

vector<vector<point>> voronoi (vector<point> P) {
    
    vector<point> u;
    u.pb({-inf_coord, -inf_coord}); u.pb({inf_coord, -inf_coord}); u.pb({inf_coord, inf_coord});
    u.pb({-inf_coord, inf_coord});
    
    if(P.size() == 1){
        vector<vector<point>> a;
        a.pb(u);
        return a;
    }
    if(P.size() == 2){
        vector<point> v1 = u, v2 = u;
        auto at = divPlane(P[0], P[1]);
        v1 = cutPolygon(v1, at.fi, at.se - at.fi); v2 = cutPolygon(v2, at.fi, at.fi - at.se);
        vector<vector<point>> a;
        a.pb(v1); a.pb(v2);
        return a;
    }
	vector<vector<point>> cells(P.size());
	forn (i, P.size()) {
		P[i].idx = i;
	}
	
	auto dt = delaunay(P);
	
	auto ch = convexHull(P);
	if (ch.size() > 2) {
		
		for (auto &tri : dt) {
			point c = getCircle(tri[0], tri[1], tri[2]).fi;
			for (auto &p : tri) {
				cells[p.idx].pb(c);
			}
		}
		
		ch = convexHullWithColinear(P);
		forn (i, ch.size()) {
			point &a = ch[i];
			point &b = ch[(i + 1) % ch.size()];
			point mit = (a + b) / 2;
			point c = (mit + (a - mit).perpendicular().unit() * inf_coord);
			
			cells[a.idx].pb(c);
			cells[b.idx].pb(c);
		}
	} else if (ch.size() == 2) {
		sort(all(P));
		
		forn (i, (int)P.size() - 1) {
			point del = (P[i + 1] - P[i]) / 2;
			point a = P[i] + (del + (P[i] - del).perpendicular().unit() * inf_coord);
			point b = P[i] + (del + (P[i] - del).perpendicular().unit() * -inf_coord);
			
			forn (j, 2) {
				cells[P[i + j].idx].pb(a);
				cells[P[i + j].idx].pb(b);
			}
		}
		
		point dir = (P.back() - P[0]).unit() * inf_coord;
		forn (j, 2) {
		    cells[P[0].idx].pb(cells[P[0].idx][j] - dir);
		    cells[P.back().idx].pb(cells[P.back().idx][j] + dir);
		}
	} else {
		cells[0] = {
			{inf_coord, inf_coord},
			{-inf_coord, inf_coord},
			{-inf_coord, -inf_coord},
			{inf_coord, -inf_coord},
		};
	}
	
	for (auto &cell : cells) {
		cell = convexHull(cell);
	}
		
	return cells;
}

bool pointInLine(const point & a, const point & v, const point & p){
	//line a+tv, point p
	return eq((p - a).cross(v), 0);
}

bool pointInSegment(const point & a, const point & b, const point & p){
	//segment ab, point p
	return pointInLine(a, b - a, p) && leq((a - p).dot(b - p), 0);
}

int main(){
	ios_base::sync_with_stdio(false);  cin.tie(NULL); cout.tie(0);
	int h, w, n; cin>>w>>h>>n;
	vector<point> v; vector<int> freq(n, 0);
	point h0(0, 0), h1(w, 0), h2(w, h), h3(0, h);
	for(int i = 0; i < n; i++){
	    int a, b; cin>>a>>b; point c(a, b);
	    if(v.size()){
	        int fl = 1;
	        for(int j = 0; j < v.size(); j++){
	            if(v[j] == c) fl = 0, freq[j]++;
	        }
	        if(fl) v.pb(c);
	    }else v.pb(c);
	}
	ld ans = 0;
	auto u = voronoi(v);
	for(int i = 0; i < u.size(); i++){ u[i] = cutPolygon(u[i], h0, h1 - h0); u[i] = cutPolygon(u[i], h1, h2 - h1);
	    u[i] = cutPolygon(u[i], h2, h3 - h2); u[i] = cutPolygon(u[i], h3, h0 - h3);
	    for(int l = 0; l < u[i].size(); l++) ans = max(ans, (v[i] - u[i][l]).length());
	}
	for(int i = 0; i < u.size(); i++){
	    if(freq[i]) continue;
	    vector<point> friends;
	    for(int j = 0; j < u.size(); j++){
	        if(i == j) continue;
	        for(int l = 0; l < u[i].size(); l++){
	            int l1 = (l + 1) % (int)u[i].size();
	            if(eq((v[i] - (u[i][l] + u[i][l1]) / 2).length(),(v[j] - (u[i][l] + u[i][l1]) / 2).length())) friends.pb(v[j]);
	        }
	    }
	    sort(all(friends));
	    for(int l = 0; l < friends.size(); l++){
	        auto v1 = u[i];
	        for(int l1 = 0; l1 < friends.size(); l1++){
	            if(l == l1) continue;
	            auto at = divPlane(friends[l], friends[l1]);
	            v1 = cutPolygon(v1, at.fi, at.se - at.fi);
	        }
	        for(int l1 = 0; l1 < v1.size(); l1++){
	            ans = max(ans, (v1[l1] - friends[l]).length());
	        }
	    }
	}
	cout<<setprecision(20)<<ans;
}
