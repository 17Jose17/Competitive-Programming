ll lagrange(vector<ll> v, ll x){
        
        int n = v.size();
        
        ll fact[maxn] = {}, inv[maxn] = {}, inv_fac[maxn] = {};
        fact[0] = fact[1] = inv[0] = inv[1] = inv_fac[0] = inv_fac[1] = 1;
        
        for(int i = 2; i <= maxn - 4; i++){
                fact[i] = (fact[i - 1] * i) % mod;
                inv[i] = (inv[mod % i] * (mod - mod / i)) % mod;
                inv_fac[i] = (inv_fac[i - 1] * inv[i]) % mod;
        }
        
        ll preff[maxn] = {}, suff[maxn] = {};
        preff[0] = suff[0] = 1;
        
        for(int i = 1; i <= maxn - 4; i++){
                preff[i] = (preff[i - 1] * (x - i)) % mod;
                suff[i] = (suff[i - 1] * (x - n + i - 1)) % mod;
        }
        
        ll ans = 0;
        for(int i = 1; i <= n; i++){
                if((n + i - 1) & 1) ans += (((((((v[i - 1] * inv_fac[n - i]) % mod) * inv_fac[i - 1]) % mod) * preff[i - 1]) % mod) * suff[n - i]) % mod;
                else ans -= (((((((v[i - 1] * inv_fac[n - i]) % mod) * inv_fac[i - 1]) % mod) * preff[i - 1]) % mod) * suff[n - i]) % mod;
                if(ans < 0) ans += mod;
        }
        return (ans % mod);
}

struct Node {
	ll l, r, res, id, lazy, islazy;
	Node(ll l = 0, ll r = 0, ll res = 1e18, ll id = -1, ll lazy = 0, ll islazy = 0): l(l), r(r), res(res), id(id), lazy(lazy), islazy(islazy){}
 
	void apply(ll v){
 
		res += v;
		lazy += v; 
		islazy = 1;
 
	}
 
	void reset(){
		lazy = 0; 
		islazy = 0;
	}
 
};
 
Node merge(Node a, Node b){
	if(a.res <= b.res){
		return Node(a.l, b.r, a.res, a.id, 0, 0);
	}
	return Node(a.l, b.r, b.res, b.id, 0, 0);
}
 
template <typename T> struct SegmentTree{
	vector<Node> ST;
	int N;
 
	SegmentTree(int n, const vector<T> & values): N(n){
		ST.resize(5 * N);
		build(1, 1, N, values);
	}
 
	void init_leaf(int curr, T value, int idx){
		ST[curr] = Node(idx, idx, value, idx, 0, 0);
	}
 
	void build(int curr, int l, int r, const vector<T> & values){
		//ST[curr].l = l, ST[curr].r = r;
		if(l == r){
			init_leaf(curr, values[l - 1], l);
		}else{
			int mid = l + (r - l) / 2;
			build(2 * curr, l, mid, values);
			build(2 * curr + 1, mid + 1, r, values);
			ST[curr] = merge(ST[2 * curr], ST[2 * curr + 1]);
		}
	}
 
	void push(int curr){
		if(ST[curr].islazy){
			ST[curr * 2].apply(ST[curr].lazy);
			ST[curr * 2 + 1].apply(ST[curr].lazy);
			ST[curr].reset();
		}
	}
 
	Node get(int curr){
		return ST[curr];
	}
 
	void update(int curr, int l, int r, int ql, int qr, T value){
		if(l > qr || r < ql){
			return;
		}else if(ql <= l && r <= qr){
			ST[curr].apply(value);
		}else{
			push(curr);
			int mid = (l + r) / 2;
			update(curr * 2, l, mid, ql, qr, value);
			update(curr * 2 + 1, mid + 1, r, ql, qr, value);
			ST[curr] = merge(ST[curr * 2], ST[curr * 2 + 1]);
		}
	}
 
	Node query(int curr, int l, int r, int ql, int qr){
		if(l > qr || r < ql)
			return Node();
		else if(ql <= l && r <= qr){
			return get(curr);
		}else{
			push(curr);
			int mid = l + (r - l) / 2;
			return merge(query(2 * curr, l, mid, ql, qr), query(2 * curr + 1, mid + 1, r, ql, qr));
		}
	}
 
	int walkingl(int curr, int l, int r, int ql, int qr, T value){
		if(l > qr || r < ql || ST[curr].res >= value) return -1;
		else if(l == r) if(ST[curr].res < value) return r; else return -1;
		else{
			push(curr);
			int mid = (l + r) / 2;
			int re = walkingl(curr * 2, l, mid, ql, qr, value);
			if(re != -1) return re;
			return walkingl(curr * 2 + 1, mid + 1, r, ql, qr, value);
		}
	}
 
	int walkingr(int curr, int l, int r, int ql, int qr, T value){
		if(l > qr || r < ql || ST[curr].res >= value) return -1;
		else if(l == r) if(ST[curr].res < value) return r; else return -1;
		else{
			push(curr);
			int mid = (l + r) / 2;
			int re = walkingr(curr * 2 + 1, mid + 1, r, ql, qr, value);
			if(re != -1) return re;
			return walkingr(curr * 2, l, mid, ql, qr, value);
		}
	}
  
	void update(int ql, int qr, T value){
		update(1, 1, N, ql, qr, value);
	}
 
	Node query(int ql, int qr){
		return query(1, 1, N, ql, qr);
	}
 
	int walkingl(int ql, int qr, T value){
		return walkingl(1, 1, N, ql, qr, value);
	}
 
	int walkingr(int ql, int qr, T value){
		return walkingr(1, 1, N, ql, qr, value);
	}
 
};
