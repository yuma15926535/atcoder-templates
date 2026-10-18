#ifndef GG_LIS
#define GG_LIS
// ---- 最長増加部分列（LIS。文字列第二世代） ----
// 列は vector<整数>（int・long long。負・10^18 級もよい）でも string でもよい。strict = true で狭義（<）、false で広義（≤）
// 位置・添字は 0 始まり。「部分列」は飛び飛びに選んだもの（順番はそのまま）
namespace str_internal {
// a[0, len) は単調増加。strict なら a[k] ≥ x の最初の k、そうでなければ a[k] > x の最初の k
// （末尾より大きければすぐ返す。ばらばらな列では分岐ありの二分探索の方が速かった）
inline int lis_slot(const long long* a, int len, long long x, bool strict) {
    if (len == 0 || (strict ? a[len - 1] < x : a[len - 1] <= x)) return len;
    int lo = -1, hi = len - 1;
    if (strict) {
        while (hi - lo > 1) {
            const int mid = (lo + hi) >> 1;
            if (a[mid] < x) lo = mid;
            else hi = mid;
        }
    } else {
        while (hi - lo > 1) {
            const int mid = (lo + hi) >> 1;
            if (a[mid] <= x) lo = mid;
            else hi = mid;
        }
    }
    return hi;
}
// end[i] = i で終わる最長の長さ。prev を渡すと復元用の 1 つ前の添字も、last には最長の列の最後の添字
inline int lis_core(const std::vector<long long>& v, bool strict, std::vector<int>* end, std::vector<int>* prev, int* last = nullptr) {
    const int n = (int)v.size();
    std::vector<long long> tail(n);
    std::vector<int> tail_id(prev || last ? n : 0);
    if (end) end->resize(n);
    if (prev) prev->resize(n);
    int L = 0;
    for (int i = 0; i < n; ++i) {
        const int p = lis_slot(tail.data(), L, v[i], strict);
        tail[p] = v[i];
        if (end) (*end)[i] = p + 1;
        if (prev) (*prev)[i] = p > 0 ? tail_id[p - 1] : -1;
        if (prev || last) tail_id[p] = i;
        if (p == L) ++L;
    }
    if (last) *last = L > 0 ? tail_id[L - 1] : -1;
    return L;
}
// 大小を逆にする（~v = -v - 1 は溢れない）
inline std::vector<long long> flipped(std::vector<long long> v) {
    for (auto& x : v) x = ~x;
    return v;
}
}  // namespace str_internal

template <class Seq>
int lis_length(const Seq& a, bool strict = true) {  // 最長増加部分列の長さ  O(N log N)
    return str_internal::lis_core(str_internal::codes_of(a), strict, nullptr, nullptr);
}
template <class Seq>
std::vector<int> lis_restore(const Seq& a, bool strict = true) {  // 最長増加部分列を 1 つ（添字の列。昇順）  O(N log N)
    std::vector<int> prev;
    int cur = -1;
    const int L = str_internal::lis_core(str_internal::codes_of(a), strict, nullptr, &prev, &cur);
    std::vector<int> r(L);
    for (int k = L - 1; k >= 0; --k) r[k] = cur, cur = prev[cur];
    return r;
}
template <class Seq>
std::vector<int> lis_ending_at(const Seq& a, bool strict = true) {  // r[i] = a[i] で終わる最長の増加部分列の長さ  O(N log N)
    std::vector<int> end;
    str_internal::lis_core(str_internal::codes_of(a), strict, &end, nullptr);
    return end;
}
template <class Seq>
std::vector<int> lis_starting_at(const Seq& a, bool strict = true) {  // r[i] = a[i] から始まる最長の増加部分列の長さ  O(N log N)
    auto v = str_internal::flipped(str_internal::codes_of(a));
    std::reverse(v.begin(), v.end());  // 後ろから読むと「減少」なので大小を逆にして増加に
    std::vector<int> end;
    str_internal::lis_core(v, strict, &end, nullptr);
    std::reverse(end.begin(), end.end());
    return end;
}

struct LisMembers {
    int length = 0;                // LIS の長さ
    std::vector<bool> in_some;     // in_some[i] = a[i] を含む LIS がある（ABC354 F）
    std::vector<bool> in_all;      // in_all[i]  = すべての LIS が a[i] を含む
};
template <class Seq>
LisMembers lis_members(const Seq& a, bool strict = true) {  // 各位置が LIS に入りうるか・必ず入るか  O(N log N)
    const auto e = lis_ending_at(a, strict), s = lis_starting_at(a, strict);
    const int n = (int)e.size();
    LisMembers r;
    for (int i = 0; i < n; ++i) r.length = std::max(r.length, e[i]);
    r.in_some.assign(n, false), r.in_all.assign(n, false);
    std::vector<int> cnt(r.length + 1, 0);  // 長さの段 k（LIS の k 番目）に入りうる位置の数
    for (int i = 0; i < n; ++i)
        if (e[i] + s[i] - 1 == r.length) r.in_some[i] = true, ++cnt[e[i]];
    for (int i = 0; i < n; ++i) r.in_all[i] = r.in_some[i] && cnt[e[i]] == 1;
    return r;
}

template <class T, class Seq>
T count_lis(const Seq& a, bool strict = true) {  // LIS の数（添字の列として数える。T は ModInt998・long long（答えが溢れたら止める）・double など）  O(N log N)
    const CompressedValues cv = compress_values(a);
    const int n = (int)cv.ids.size(), K = cv.size();
    std::vector<int> len(K + 1, 0);  // フェニック木：値の番号の接頭辞での（最長の長さ, その数, 溢れたか）
    std::vector<T> cnt(K + 1, T(0));
    std::vector<char> over(K + 1, 0);
    auto merge = [&](int& l0, T& c0, bool& o0, int l1, const T& c1, bool o1) {
        if (l1 > l0) l0 = l1, c0 = c1, o0 = o1;
        else if (l1 == l0 && l1 > 0) str_internal::add_count(c0, o0, c1, o1);
    };
    for (int i = 0; i < n; ++i) {
        const int id = cv.ids[i];
        int bl = 0;
        T bc = T(0);
        bool bo = false;
        for (int k = strict ? id : id + 1; k > 0; k -= k & -k) merge(bl, bc, bo, len[k], cnt[k], over[k]);
        if (bl == 0) bc = T(1);
        ++bl;
        for (int k = id + 1; k <= K; k += k & -k) {
            bool o = over[k];
            merge(len[k], cnt[k], o, bl, bc, bo);
            over[k] = o;
        }
    }
    int bl = 0;
    T bc = T(0);
    bool bo = false;
    for (int k = K; k > 0; k -= k & -k) merge(bl, bc, bo, len[k], cnt[k], over[k]);
    if (bo) str_internal::count_overflow("count_lis");
    return n == 0 ? T(1) : bc;  // 空の列の LIS は空の 1 通り
}

template <class Seq>
int lds_length(const Seq& a, bool strict = true) {  // 最長減少部分列の長さ  O(N log N)
    return str_internal::lis_core(str_internal::flipped(str_internal::codes_of(a)), strict, nullptr, nullptr);
}
template <class Seq>
int min_increasing_cover(const Seq& a, bool strict = true) {  // 増加部分列（strict なら狭義）に分ける最小の数 = 最長の「広義 / 狭義」減少部分列（Dilworth。ABC134 E）
    return lds_length(a, !strict);
}
template <class Seq>
int min_decreasing_cover(const Seq& a, bool strict = true) {  // 減少部分列に分ける最小の数
    return lis_length(a, !strict);
}

template <class X, class Y>
std::vector<int> lis_2d(const std::vector<std::pair<X, Y>>& ps) {  // x も y も狭義に増える最長の列（添字。x の小さい順。ABC038 D）  O(N log N)
    const int n = (int)ps.size();
    std::vector<int> ord(n);
    for (int i = 0; i < n; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](int i, int j) {  // x の昇順、同じ x は y の降順（同じ x を 2 つ選ばないように）
        if (ps[i].first < ps[j].first) return true;
        if (ps[j].first < ps[i].first) return false;
        if (ps[j].second < ps[i].second) return true;
        if (ps[i].second < ps[j].second) return false;
        return i < j;
    });
    // y を < だけで順位に（double・unsigned long long などでも大小を保つ）
    std::vector<int> by(n), ry(n);
    for (int i = 0; i < n; ++i) by[i] = i;
    std::sort(by.begin(), by.end(), [&](int i, int j) { return ps[i].second < ps[j].second; });
    for (int k = 0, r = -1; k < n; ++k) {
        if (k == 0 || ps[by[k - 1]].second < ps[by[k]].second) ++r;
        ry[by[k]] = r;
    }
    std::vector<long long> ys(n);
    for (int k = 0; k < n; ++k) ys[k] = ry[ord[k]];
    std::vector<int> r = lis_restore(ys, true);
    for (int& k : r) k = ord[k];
    return r;
}
#endif  // GG_LIS
