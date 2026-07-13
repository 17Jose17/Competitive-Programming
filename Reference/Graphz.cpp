void dfs(int root, int pa, vector<vector<pair<int, int>>> & L, vector<int> & vis, vector<int> & tin, vector<int> & low, int & tim, set<tuple<int, int, int>> & b){
    vis[root] = 1;
    tin[root] = low[root] = tim++;
    bool p = false;
    for(auto i : L[root]){
        if(i.fi == pa && !p){
            p = true; continue;
        }
        if(vis[i.fi]){
            low[root] = min(low[root], tin[i.fi]);
        }else{
            dfs(i.fi, root, L, vis, tin, low, tim, b);
            low[root] = min(low[root], low[i.fi]);
            if(low[i.fi] > tin[root]){
                b.insert({min(root, i.fi), max(root, i.fi), i.se});
            }
        }
    }
}
