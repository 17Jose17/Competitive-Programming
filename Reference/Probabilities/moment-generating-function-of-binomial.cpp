momentGeneratingBinomial(ll n, ll m, ll k){
        ll p = expBin(m, mod - 2);
        ll DP[k + 5][k + 5] = {}; DP[1][1] = 1;
        for(int i = 2; i <= k + 2; i++){
                for(int j = 1; j <= k + 2; j++){
                        DP[i][j] = (DP[i - 1][j] * j + DP[i - 1][j - 1]) % mod;
                }
        }
        ll pref[k + 5] = {}; 
        ll u = n, res = 1;
        for(int i = 1; i <= k; i++){
                pref[i] = res = (res * u) % mod;
                u--;
        }
        ll ans = 0;
        for(int i = 1; i <= k; i++){
                ans += (((DP[k][i] * pref[i]) % mod) * expBin(p, i)) % mod;
                if(ans >= mod) ans -= mod;
        }
        return ans;
}
