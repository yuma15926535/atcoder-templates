#ifndef GG_REVERSIBLE_DEQUE
#define GG_REVERSIBLE_DEQUE
// ---- 前後に足せて O(1) で反転できる列・回転と反転が O(1) の配列（文字列第二世代） ----
// T = char で文字列として使う（to_string・substr が string を返す）。位置は 0 始まり。大量に足す前に reserve(n) で置き場を取れる

// ---- ReversibleDeque<T>：前後への追加・削除と反転（環状バッファ ＋ 向きの印） ----
template <class T>
struct ReversibleDeque {
    using seq_type = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    ReversibleDeque() = default;
    template <class Seq, class = decltype(str_internal::seq(std::declval<const Seq&>()).size())>
    explicit ReversibleDeque(const Seq& s) { push_back(s); }

    int size() const { return n_; }
    bool empty() const { return n_ == 0; }
    void clear() { head_ = 0, n_ = 0, rev_ = false; }
    void reserve(int cap) {  // 置き場を先に cap 個分取る（大量に足す前に）
        if (cap > n_) grow(cap - n_);
    }
    void reverse() { rev_ = !rev_; }  // 全体を反転  O(1)
    bool reversed() const { return rev_; }

    void push_back(const T& x) {  // 末尾に 1 つ  O(1)（ならし）
        T v(x);  // x がこの列の要素でも、置き場を広げる前に写しておく（d.push_back(d.front()) など）
        grow(1);
        if (!rev_) buf_[(head_ + n_) & mask_] = std::move(v);  // 置き場の後ろ
        else head_ = (head_ - 1) & mask_, buf_[head_] = std::move(v);  // 反転中は置き場の前
        ++n_;
    }
    void push_front(const T& x) {  // 先頭に 1 つ  O(1)（ならし）
        T v(x);
        grow(1);
        if (rev_) buf_[(head_ + n_) & mask_] = std::move(v);
        else head_ = (head_ - 1) & mask_, buf_[head_] = std::move(v);
        ++n_;
    }
    template <class Seq, class = decltype(str_internal::seq(std::declval<const Seq&>()).size())>
    void push_back(const Seq& s0) {  // 末尾に列をそのままの順でつなげる（自分自身も可）  O(|s|)
        if constexpr (std::is_same<Seq, ReversibleDeque>::value)
            if (&s0 == this) return push_back(to_vector());
        const auto& s = str_internal::seq(s0);
        const std::size_t m = s.size();
        grow((int)m);
        for (std::size_t i = 0; i < m; ++i) push_back((T)s[i]);
    }
    template <class Seq, class = decltype(str_internal::seq(std::declval<const Seq&>()).size())>
    void push_front(const Seq& s0) {  // 先頭に列をそのままの順でつなげる（s + 今の列。自分自身も可）  O(|s|)
        if constexpr (std::is_same<Seq, ReversibleDeque>::value)
            if (&s0 == this) return push_front(to_vector());
        const auto& s = str_internal::seq(s0);
        const std::size_t m = s.size();
        grow((int)m);
        for (std::size_t i = m; i-- > 0;) push_front((T)s[i]);
    }
    T pop_back() {  // 末尾を取り出す  O(1)
        if (n_ == 0) str_internal::fail("ReversibleDeque::pop_back", "空");
        T x = back();
        if (!rev_) --n_;
        else head_ = (head_ + 1) & mask_, --n_;
        return x;
    }
    T pop_front() {  // 先頭を取り出す  O(1)
        if (n_ == 0) str_internal::fail("ReversibleDeque::pop_front", "空");
        T x = front();
        if (rev_) --n_;
        else head_ = (head_ + 1) & mask_, --n_;
        return x;
    }
    const T& front() const { return at(0, "ReversibleDeque::front"); }
    const T& back() const { return at(n_ - 1, "ReversibleDeque::back"); }
    const T& operator[](int i) const { return at(i, "ReversibleDeque::operator[]"); }  // 先頭から i 番目  O(1)
    T& operator[](int i) { return const_cast<T&>(at(i, "ReversibleDeque::operator[]")); }
    std::vector<T> to_vector() const {  // 先頭から順に  O(N)
        std::vector<T> r(n_);
        for (int i = 0; i < n_; ++i) r[i] = buf_[phys(i)];
        return r;
    }
    std::string to_string() const {  // 文字列で（T = char）  O(N)
        static_assert(std::is_same<T, char>::value, "to_string は T = char のときだけ（ほかは to_vector）");
        std::string r(n_, '\0');
        for (int i = 0; i < n_; ++i) r[i] = buf_[phys(i)];
        return r;
    }
    seq_type substr(int pos, int len) const {  // [pos, pos + len) を取り出す（T = char なら string）  O(len)
        str_internal::check_range("ReversibleDeque::substr", pos, (long long)pos + len, n_);
        seq_type r(len, T());
        for (int i = 0; i < len; ++i) r[i] = buf_[phys(pos + i)];
        return r;
    }

   private:
    int phys(int i) const { return (head_ + (rev_ ? n_ - 1 - i : i)) & mask_; }  // 先頭から i 番目の場所
    const T& at(int i, const char* func) const {
        if (i < 0 || i >= n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜" + std::to_string(n_ - 1) + " の外");
        return buf_[phys(i)];
    }
    void grow(int add) {
        if ((long long)n_ + add <= (long long)buf_.size()) return;
        std::size_t cap = buf_.empty() ? 16 : buf_.size();
        while (cap < (std::size_t)n_ + (std::size_t)add) cap <<= 1;
        std::vector<T> nb(cap);
        for (int i = 0; i < n_; ++i) nb[i] = buf_[(head_ + i) & mask_];  // 置き場の順のまま詰める（向きの印はそのまま）
        buf_.swap(nb);
        head_ = 0, mask_ = (int)cap - 1;
    }
    std::vector<T> buf_;
    int head_ = 0, n_ = 0, mask_ = -1;
    bool rev_ = false;
};

// ---- RotatableArray<T>：大きさが決まった配列の回転と反転を O(1)（ずれと向きを持つ） ----
template <class T>
struct RotatableArray {
    using seq_type = std::conditional_t<std::is_same<T, char>::value, std::string, std::vector<T>>;
    RotatableArray() = default;
    template <class Seq>
    explicit RotatableArray(const Seq& s0) {
        const auto& s = str_internal::seq(s0);
        a_.assign(s.begin(), s.end());
        n_ = (long long)a_.size();
    }
    int size() const { return (int)n_; }
    void rotate_left(long long k) {  // 左に k 回す（新しい [i] = 元の [i + k]。負なら右）  O(1)
        if (n_ == 0) return;
        k %= n_;
        off_ = ((off_ + dir_ * k) % n_ + n_) % n_;
    }
    void rotate_right(long long k) { rotate_left(-(k % (n_ == 0 ? 1 : n_))); }  // 右に k 回す（新しい [i] = 元の [i - k]）  O(1)
    void reverse() {  // 全体を反転  O(1)
        if (n_ == 0) return;
        off_ = ((off_ + dir_ * (n_ - 1)) % n_ + n_) % n_;
        dir_ = -dir_;
    }
    const T& operator[](int i) const { return a_[phys(i, "RotatableArray::operator[]")]; }  // 今の i 番目  O(1)
    T& operator[](int i) { return a_[phys(i, "RotatableArray::operator[]")]; }
    const T& get(int i) const { return a_[phys(i, "RotatableArray::get")]; }
    void set(int i, const T& v) { a_[phys(i, "RotatableArray::set")] = v; }
    std::vector<T> to_vector() const {  // 今の順に  O(N)
        std::vector<T> r(n_);
        for (int i = 0; i < (int)n_; ++i) r[i] = (*this)[i];
        return r;
    }
    std::string to_string() const {  // 文字列で（T = char）  O(N)
        static_assert(std::is_same<T, char>::value, "to_string は T = char のときだけ（ほかは to_vector）");
        std::string r((std::size_t)n_, '\0');
        for (int i = 0; i < (int)n_; ++i) r[i] = (*this)[i];
        return r;
    }

   private:
    std::size_t phys(int i, const char* func) const {
        if (i < 0 || i >= n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜" + std::to_string(n_ - 1) + " の外");
        long long p = off_ + dir_ * (long long)i;
        if (p >= n_) p -= n_;
        else if (p < 0) p += n_;
        return (std::size_t)p;
    }
    std::vector<T> a_;
    long long n_ = 0, off_ = 0, dir_ = 1;
};
#endif  // GG_REVERSIBLE_DEQUE
