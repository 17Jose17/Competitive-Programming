// Factorizar n numeros en O(nlogA)
vector<int> fastFact(int n){
        vector<int> is(n + 1, 0);
        is[0] = is[1] = 0;
        for(int i = 4; i <= n; i += 2) is[i] = 2;
                for(int i = 3; i <= n; i += 2){
                        if(is[i] == 0){
                        if((long long) i * i <= n)
                                for(int j = i * i; j <= n; j += 2 * i)
                                        if(!is[j]) is[j] = i;
                }
        }
        return is;
}
