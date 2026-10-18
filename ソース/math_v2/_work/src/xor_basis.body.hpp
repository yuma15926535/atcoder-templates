#ifndef GG_XOR_BASIS
#define GG_XOR_BASIS
// ---- XOR の基底（数学第二世代）  XorBasis（64 ビットの値の XOR で作れる値）・PrefixXorBasis（区間の要素）・Gf2System（mod 2 の連立方程式） ----
// 「作れる値」= いくつかの要素（0 個でもよい）の XOR。0 個なら 0
struct XorBasis {
    using u64 = unsigned long long;
    XorBasis() {
        for (int i = 0; i < 64; ++i) b_[i] = 0, comp_[i] = 0;
    }
    // x を加える。新しく作れる値が増えたら true（既に作れるなら false）  O(64)
    bool insert(u64 x) {
        const int id = inserted_++;
        u64 c = 0;  // 引いた基底がどの slot の XOR か
        // 簡約した基底は他の基底の先頭ビットを持たないので、x の「基底の先頭ビット」の位置の基底を全部 XOR すればよい
        for (u64 t = x & have_; t; t &= t - 1) {
            const int i = __builtin_ctzll(t);
            x ^= b_[i], c ^= comp_[i];
        }
        if (x == 0) return false;
        const int h = 63 - __builtin_clzll(x), s = rank_++;
        slot_id_[s] = id;
        c ^= 1ULL << s;
        for (u64 t = have_ >> h >> 1; t; t &= t - 1) {  // 上の基底から x の先頭ビットを消す（簡約を保つ）
            const int i = h + 1 + __builtin_ctzll(t);
            if (b_[i] >> h & 1) b_[i] ^= x, comp_[i] ^= c;
        }
        b_[h] = x, comp_[h] = c, have_ |= 1ULL << h;
        return true;
    }
    bool contains(u64 x) const { return reduce(x) == 0; }  // x が作れるか
    u64 reduce(u64 x) const {                              // min_{s は作れる値} (x XOR s)
        for (u64 t = x & have_; t; t &= t - 1) x ^= b_[__builtin_ctzll(t)];
        return x;
    }
    u64 min_xor(u64 x = 0) const { return reduce(x); }  // x XOR (作れる値) の最小（x = 0 なら 0）
    u64 max_xor(u64 x = 0) const {                      // x XOR (作れる値) の最大（先頭ビットが 0 のところの基底を全部 XOR）
        for (u64 t = ~x & have_; t; t &= t - 1) x ^= b_[__builtin_ctzll(t)];
        return x;
    }
    // 作れる値を小さい順に並べた k 番目（0 から。0 番目は 0。k < 2^rank）
    u64 kth(u64 k) const {
        if (rank_ < 64 && (k >> rank_) != 0) math_internal::fail("XorBasis::kth", "k < 2^rank で使う");
        u64 x = 0;
        int j = 0;
        for (u64 t = have_; t; t &= t - 1, ++j)
            if (k >> j & 1) x ^= b_[__builtin_ctzll(t)];
        return x;
    }
    Int128 count_le(u64 x) const {  // x 以下の作れる値の個数（kth が増えていくことを使い、k を上のビットから決める）  O(64)
        u64 k = 0, cur = 0;
        int j = rank_;
        for (int i = 63; i >= 0; --i) {
            if (!b_[i]) continue;
            --j;
            if ((cur ^ b_[i]) <= x) cur ^= b_[i], k |= 1ULL << j;  // k のビット j を立てて、下のビットは 0 の値
        }
        return (Int128)k + 1;
    }
    int rank() const { return rank_; }
    Int128 count() const { return (Int128)1 << rank_; }  // 作れる値の個数 2^rank
    int inserted() const { return inserted_; }          // insert を呼んだ回数
    // 0 個でない要素の XOR で 0 を作れるか（insert したのに増えなかった要素があるか）
    bool zero_by_nonempty() const { return inserted_ > rank_; }
    void merge(const XorBasis& o) {  // o の作れる値も作れるようにする（represent の番号は o の基底の分だけ増える）
        for (int i = 0; i < 64; ++i)
            if (o.b_[i]) insert(o.b_[i]);
    }
    std::vector<u64> basis() const {  // 基底（簡約済み。先頭ビットの大きい順）
        std::vector<u64> r;
        for (int i = 63; i >= 0; --i)
            if (b_[i]) r.push_back(b_[i]);
        return r;
    }
    // x を作る要素の番号（insert を呼んだ順に 0, 1, …）。作れなければ ok = false
    Maybe<std::vector<int>> represent(u64 x) const {
        u64 c = 0;
        for (int i = 63; i >= 0; --i)
            if (x >> i & 1) {
                if (!b_[i]) return {false, {}};
                x ^= b_[i], c ^= comp_[i];
            }
        std::vector<int> r;
        for (int s = 0; s < rank_; ++s)
            if (c >> s & 1) r.push_back(slot_id_[s]);
        std::sort(r.begin(), r.end());
        return {true, r};
    }

  private:
    u64 b_[64];        // b_[i]：先頭ビットが i の基底（なければ 0）
    u64 comp_[64];     // b_[i] がどの slot の XOR か
    int slot_id_[64];  // slot s の要素が何番目に insert されたか
    u64 have_ = 0;     // 基底がある先頭ビットの集合
    int rank_ = 0, inserted_ = 0;
};

// ---- 区間の要素の XOR：a_0, a_1, … を順に push_back して、今の末尾までの [l, size()) について答える ----
// クエリは右端の小さい順に並べ替えて答える（各ビットの基底に「いちばん右の番号」を残す）  push_back・各クエリ O(64)
struct PrefixXorBasis {
    using u64 = unsigned long long;
    PrefixXorBasis() {
        for (int i = 0; i < 64; ++i) b_[i] = 0, pos_[i] = -1;
    }
    void push_back(u64 x) {  // a_{size()} = x を足す
        int p = n_++;
        for (int i = 63; i >= 0 && x; --i) {
            if (!(x >> i & 1)) continue;
            if (!b_[i]) {
                b_[i] = x, pos_[i] = p;
                return;
            }
            if (pos_[i] < p) std::swap(b_[i], x), std::swap(pos_[i], p);  // 新しい（右の）方を基底に残す
            x ^= b_[i];
        }
    }
    int size() const { return n_; }
    bool contains(int l, u64 x) const {  // [l, size()) の要素（0 個でもよい）の XOR で x を作れるか
        check(l, "PrefixXorBasis::contains");
        for (int i = 63; i >= 0; --i)
            if ((x >> i & 1) && b_[i] && pos_[i] >= l) x ^= b_[i];
        return x == 0;
    }
    u64 max_xor(int l, u64 x = 0) const {  // x XOR ([l, size()) で作れる値) の最大
        check(l, "PrefixXorBasis::max_xor");
        for (int i = 63; i >= 0; --i)
            if (b_[i] && pos_[i] >= l && !(x >> i & 1)) x ^= b_[i];
        return x;
    }
    int rank(int l) const {  // [l, size()) の要素の基底の大きさ
        check(l, "PrefixXorBasis::rank");
        int r = 0;
        for (int i = 0; i < 64; ++i) r += b_[i] && pos_[i] >= l;
        return r;
    }
    XorBasis basis(int l) const {  // [l, size()) の要素の基底（XorBasis として。kth なども使える）
        check(l, "PrefixXorBasis::basis");
        XorBasis r;
        for (int i = 63; i >= 0; --i)
            if (b_[i] && pos_[i] >= l) r.insert(b_[i]);
        return r;
    }

  private:
    u64 b_[64];
    int pos_[64];
    int n_ = 0;
    void check(int l, const char* f) const {
        if (l < 0 || l > n_) math_internal::fail(f, "0 ≤ l ≤ size() で使う");
    }
};

// ---- mod 2 の連立一次方程式（変数 x_0〜x_{m−1} ∈ {0, 1}。ビットで持つ）----
// 式を add したあと rank / solve / kernel。消去は O(式の数 × rank × m / 64)
struct Gf2System {
    using u64 = unsigned long long;
    explicit Gf2System(int m) : m_(m), w_((m + 64) / 64) {
        if (m < 0) math_internal::fail("Gf2System", "m ≥ 0 で使う");
    }
    void add(const std::vector<int>& vars, int rhs) {  // vars の変数の XOR = rhs（同じ変数が 2 回あれば消える）
        std::vector<u64> row((size_t)w_, 0);
        for (int v : vars) {
            if (v < 0 || v >= m_) math_internal::fail("Gf2System::add", "変数の番号は 0〜m−1");
            row[(size_t)(v >> 6)] ^= 1ULL << (v & 63);
        }
        if (rhs & 1) row[(size_t)(m_ >> 6)] ^= 1ULL << (m_ & 63);
        rows_.push_back(std::move(row)), dirty_ = true;
    }
    void add_row(const std::vector<int>& coef, int rhs) {  // Σ coef_j x_j = rhs（coef は長さ m の 0/1）
        if ((int)coef.size() != m_) math_internal::fail("Gf2System::add_row", "係数の長さは m");
        std::vector<int> vars;
        for (int j = 0; j < m_; ++j)
            if (coef[(size_t)j] & 1) vars.push_back(j);
        add(vars, rhs);
    }
    int rank() { return eliminate(), (int)piv_.size(); }
    bool solvable() { return eliminate(), ok_; }
    int free_count() { return eliminate(), m_ - (int)piv_.size(); }  // 解があれば解の数は 2^free_count
    Maybe<std::vector<int>> solve() {                                // 解の 1 つ（自由な変数は 0）
        eliminate();
        if (!ok_) return {false, {}};
        std::vector<int> x((size_t)m_, 0);
        for (size_t i = 0; i < piv_.size(); ++i) x[(size_t)piv_[i]] = bit(rows_[i], m_);
        return {true, x};
    }
    std::vector<std::vector<int>> kernel() {  // 右辺を 0 にした方程式の解の基底（解は solve() + kernel の XOR で全部）
        eliminate();
        std::vector<bool> is_piv((size_t)m_, false);
        for (int c : piv_) is_piv[(size_t)c] = true;
        std::vector<std::vector<int>> res;
        for (int f = 0; f < m_; ++f) {
            if (is_piv[(size_t)f]) continue;
            std::vector<int> v((size_t)m_, 0);
            v[(size_t)f] = 1;
            for (size_t i = 0; i < piv_.size(); ++i) v[(size_t)piv_[i]] = bit(rows_[i], f);
            res.push_back(v);
        }
        return res;
    }
    int equations() const { return (int)rows_.size(); }

  private:
    int m_, w_;
    std::vector<std::vector<u64>> rows_;
    std::vector<int> piv_;
    bool dirty_ = true, ok_ = true;
    static int bit(const std::vector<u64>& r, int j) { return (int)(r[(size_t)(j >> 6)] >> (j & 63) & 1); }
    void eliminate() {  // 行の簡約階段形に
        if (!dirty_) return;
        dirty_ = false, piv_.clear();
        const int n = (int)rows_.size();
        int r = 0;
        for (int c = 0; c < m_ && r < n; ++c) {
            const int wc = c >> 6;
            const u64 bc = 1ULL << (c & 63);
            int p = r;
            while (p < n && !(rows_[(size_t)p][(size_t)wc] & bc)) ++p;
            if (p == n) continue;
            std::swap(rows_[(size_t)p], rows_[(size_t)r]);
            const u64* pr = rows_[(size_t)r].data();
            for (int i = 0; i < n; ++i) {
                if (i == r || !(rows_[(size_t)i][(size_t)wc] & bc)) continue;
                u64* ri = rows_[(size_t)i].data();
                for (int k = wc; k < w_; ++k) ri[k] ^= pr[k];
            }
            piv_.push_back(c), ++r;
        }
        ok_ = true;
        for (int i = r; i < n; ++i)
            if (bit(rows_[(size_t)i], m_)) ok_ = false;
    }
};
#endif  // GG_XOR_BASIS
