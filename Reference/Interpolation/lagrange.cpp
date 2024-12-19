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
